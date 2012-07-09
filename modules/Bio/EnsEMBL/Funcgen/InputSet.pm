#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSet
#

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::InputSet - A module to represent InputSet object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::InputSet;

#Create an InputSet

my $inp_set = Bio::EnsEMBL::Funcgen::InputSet->new
                 (
	              -DBID            => $dbID,
                  -ADAPTOR         => $self,
                  -EXPERIMENT      => $exp,
                  -FEATURE_TYPE => $ftype,
                  -CELL_TYPE    => $ctype,
                  -FORMAT       => 'READ_FORMAT',
                  -VENDOR       => 'SOLEXA',
                  -NAME         => 'ExpSet1',
                  -REPLICATE    => 1,
                 );

# Add some InputSubsets

$inp_set->add_new_subsets($subset_name, $




=head1 DESCRIPTION

An InputSet object provides a generic container for any non-array based feature import, 
allowing tracking of file import via the status table and integration into Data and FeatureSets to
provide traceability to the source experiment from a given FeatureSet.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::InputSet;

use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Set;
use Bio::EnsEMBL::Analysis;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);


=head2 new



  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSet->new(
                                                                     -EXPERIMENT    => $exp,
                                                                     -FEATURE_TYPE  => $ftype,
                                                                     -CELL_TYPE     => $ctype,
                                                                     -FORMAT        => 'READ_FORMAT',
                                                                     -VENDOR        => 'SOLEXA',
                                                                     -NAME          => 'ExpSet1',
                                                                     -ANALYSIS      => $anal,
                                                                     -FEATURE_CLASS => 'annotated',
                                                                     );

  Do we want to define subsets likes this or are we more likely to add them one by one?

  Description: Constructor for InputSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::InputSet
  Exceptions : Throws if no Experiment defined
               Throws if CellType or FeatureType are not valid or stored
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  #Add set_type here to overwrite default ref parsing in Set::set_type
  #This need to stay like this until we patch the DB
  my $self = $class->SUPER::new(@_);	
 
  my ($exp, $format, $vendor, $rep)
    = rearrange(['EXPERIMENT', 'FORMAT', 'VENDOR', 'REPLICATE'], @_);
    
  if (! (ref $exp && $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') && $exp->dbID())){
	throw('Must specify a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }

  
  #These are set in Set, just validate here
  throw ('Must provide a FeatureType') if(! defined $self->feature_type);
  throw ('Must provide a CellType') if(! defined $self->cell_type);

  my $type = $self->feature_class;

  #Need to move these types to config

  if(! ($type && grep /^${type}$/, ('annotated', 'result', 'segmentation'))){
	throw("You must define a valid InputSet feature_class e.g. 'annotated' or 'result'");
  }

  if(($type eq 'result') &&
	 ($format ne 'SEQUENCING')){
	throw('InputSet does not yet support a result type InputSet which does not have the \'SEQUENCING\' format');
	
  }


  #if(! defined $self->analysis){
  ##default analysis hack for v47
  ##Set directly to avoid dbID boolean check
  #This is to support supporting_set cache in data_set?
  $self->{'analysis'} = Bio::EnsEMBL::Analysis->new
	(-logic_name => 'external',
	 -id       => 0,#??someone needs to rewrite analysis
	);
  
  #Change to direct setting for speed
  $self->{format}     = $format;
  $self->{vendor}     = $vendor;
  $self->{replicate}  = $rep;
  $self->{experiment} = $exp;
  $self->{subsets}    = {};
  
  return $self;
}


=head2 add_new_subset

  Arg [1]    : string - sub set name e.g. the file name (not path as we're restricted to 30 chars)
  Arg [2]    : Bio::EnsEMBL::Funcgen::InputSubset - optional
               If not defined will create a sparse InputSubset based on the name
  Example    : $expset->add_new_subset($ss_name, $exp_subset);
  Description: Adds input_subset
  Returntype : none
  Exceptions : Throws if set is already present
               Throws if InputSubset is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

#Do we still use the optional subset function?

sub add_new_subset {
  my ($self, $ss_name, $exp_sset) = @_;
	
  #Need to test $ss_name here
  if(! ($ss_name && ref(\$ss_name) eq 'SCALAR')){#ref($exp_sset) would be 'REF'
	throw('You must pass a InputSubset name');
  }

  if($self->get_subset_by_name($ss_name)){
	throw("Subset $ss_name is already present in this InputSet, maybe you need to alter the filename?");
  }

  if(defined $exp_sset){

	if(!(ref($exp_sset) && $exp_sset->isa('Bio::EnsEMBL::Funcgen::InputSubset') && $exp_sset->dbID())){
	  throw('InputSubsets must be valid and stored');
	}
  }
  else{
	
	$exp_sset = Bio::EnsEMBL::Funcgen::InputSubset->new(
														-name => $ss_name,
														-input_set => $self,
													   );
  }

  $self->{subsets}{$ss_name} = $exp_sset;

  return $self->{subsets}{$ss_name};
}


=head2 get_Experiment

  Example    : my $exp = $exp_set->get_Experiment();
  Description: Getter for the Experiment of this DataSet.
  Returntype : Bio::EnsEMBL::Fuuncgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Experiment{  return $_[0]->{experiment}; }


=head2 get_InputSubsets

  Example    : my @subsets = @{$exp_set->get_InputSubsets()};
  Description: Getter for the InputSubsets for this InputSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_InputSubsets{
  my ($self)  = shift;

  return [ values %{$self->{'subsets'}} ];
}




=head2 get_subset_by_name

  Example    : my $subsets = $exp_set->get_subset_by_name('subset1');
  Description: Getter for the subset of a given name for this InputSet.
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_subset_by_name{
  my ($self, $name) = @_;
  return (exists $self->{'subsets'}{$name}) ? $self->{'subsets'}{$name} : undef;
}


=head2 get_subset_names

  Example    : my @subset_names = @{$exp_set->get_subset_names()};
  Description: Getter for the subset names for this InputSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_subset_names{
  my ($self) = shift;
  return [ keys %{$self->{'subsets'}} ];
}




=head2 vendor

  Arg[1]     : String - vendor e.g. ILLUMINA
  Example    : my $iset_vendor = $iset->vendor;
  Description: Getter for the vendor attribute of this InputSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub vendor {  return $_[0]->{vendor}; }


=head2 format

  Arg[1]     : string - format i.e. product type/format
  Example    : my $iset_format = $iset->format;
  Description: Getter for the format attribute of this InputSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub format {  return $_[0]->{format}; }


=head2 replicate

  Arg[1]     : Integer - replicate 0 = merged or NA, >0 refers to individual replicate
  Example    : if($iset->replicate){ #Do something replicate specific in here }
  Description: Getter for the replicate attribute of this InputSet.
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub replicate {  return $_[0]->{replicate}; }



=head2 source_info

  Example    : my $source_info = $input_set->source_info;
  Description: Getter for the experiment source info i.e. [ $label, $url ]
  Returntype : Listref
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

#Currently handling redundant/absent InputSubset data

sub source_info{
  my $self = shift;

  if(! defined $self->{source_info}){
    #could have data_url as highest priority here
    #but we need to ensure removal when adding archive ids
    #so we link to the archive and not the old data url
    
    my $exp_group = $self->get_Experiment->experimental_group;
    my %source_info; #Handles redundant InputSubsets
    my ($proj_name, $proj_link, $source_label, $source_link);

    if($exp_group->is_project){
      $proj_name = $exp_group->name;
      $proj_link = $exp_group->url;
    }

    foreach my $isset(@{$self->get_InputSubsets}){

      if(defined $isset->archive_id ){
        $source_label = $isset->archive_id;
      
        if(! exists $source_info{$source_label}){
          $source_info{$source_label} = [$source_label, undef];
          #source_link can is undef here as archive_id overrides display url
          #undef links will automatically go to the SRA
        }
      }
      elsif(defined $proj_name){
        #$source_label = $self->experimental_group->name;
        $source_link  = $isset->display_url || $proj_link;
        
        if(! exists $source_info{$source_link}){
          $source_info{$source_link} = [$proj_name, $source_link];
        }
      }
    }
    
    $self->{source_info} = [values %source_info];
  }

  return $self->{source_info};
}



1;

