#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSet
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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
	                -DBID         => $dbID,
                  -ADAPTOR      => $self,
                  -EXPERIMENT   => $exp,
                  -FEATURE_TYPE => $ftype,
                  -CELL_TYPE    => $ctype,
                  -FORMAT       => 'READ_FORMAT',
                  -VENDOR       => 'SOLEXA',
                  -NAME         => 'SRR00000.fastq.gz',
                  -REPLICATE    => 1, # >0 for specific replicate or 0 for merged
                 );



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


my %valid_types = (
                   annotated      => undef,
                   result         => undef,
                   segmentation   => undef,
                   dna_methylation => undef,
                  );


#Need the underscore to separate the words
#for handling conversion to ucfirst feature class namespace


=head2 new



  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSet->new
                            (
                             -EXPERIMENT    => $exp,
                             -FEATURE_TYPE  => $ftype,
                             -CELL_TYPE     => $ctype,
                             -FORMAT        => 'READ_FORMAT',
                             -VENDOR        => 'SOLEXA',
                             -NAME          => 'SRR00000.fastq.gz',
                             -ANALYSIS      => $anal,
                             -FEATURE_CLASS => 'annotated',
                             -REPLICATE     => 1, # >0 for specific replicate or 0 for merged
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

  if(! (defined $type && exists $valid_types{$type})){
    throw("You must define a valid InputSet feature_class($type), one of: ".join("\t", keys %valid_types));
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


=head2 _add_new_subset

  Arg [1]    : Bio::EnsEMBL::Funcgen::InputSubset
  Example    : $input_set->_add_new_subset($input_subset);
  Description: Adds an InputSubset to this InputSet
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::InputSubset->new
  Status     : At Risk

=cut

sub _add_new_subset {
  my ($self, $inp_sset) = @_;

  #No validation here as this is all done in InputSubset->new

  $self->{subsets}{$inp_sset->name} = $inp_sset;
  return;
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


### DEPRECATED ###

sub add_new_subset {
  #throw as this will have already been done, and the validation is done implicitly in InputSubset->new
  throw('add_new_subset was deprecated in v69, _add_new_subset is now called directly from InputSubset->new');
}

1;

