#
# Ensembl module for Bio::EnsEMBL::Funcgen::InputSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::InputSet - A module to represent InputSet object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExpeimentalSet;

my $data_set = Bio::EnsEMBL::Funcgen::InputSet->new(
	                                                      -DBID            => $dbID,
							 					          -ADAPTOR         => $self,
                                                          -EXPERIMENT   => $exp,
                                                          -FEATURE_TYPE => $ftype,
                                                          -CELL_TYPE    => $ctype,
                                                          -FORMAT       => 'READ_FORMAT',
                                                          -VENDOR       => 'SOLEXA',
                                                          -NAME         => 'ExpSet1',
                                                          );



=head1 DESCRIPTION

An InputSet object provides a generic container for any non-array based feature import, 
allowing tracking of file import via the status table and integration into Data and FeatureSets to
provide traceability to the source experiment from a given FeatureSet.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::InputSet;

use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Set;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);#change to Set once we have implemented analysis properly


=head2 new



  Example    : my $eset = Bio::EnsEMBL::Funcgen::InputSet->new(
                                                                     -EXPERIMENT   => $exp,
                                                                     -FEATURE_TYPE => $ftype,
                                                                     -CELL_TYPE    => $ctype,
                                                                     -FORMAT       => 'READ_FORMAT',
                                                                     -VENDOR       => 'SOLEXA',
                                                                     -NAME         => 'ExpSet1',
                                                                     -ANALYSIS     => $anal,
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
  my $self = $class->SUPER::new(@_);#, '-SET_TYPE', 'experimental');
	
  #do we need to add $fg_ids to this?  Currently maintaining one feature_group focus.(combi exps?)
  my ($exp, $format, $vendor)
    = rearrange(['EXPERIMENT', 'FORMAT', 'VENDOR'], @_);
    
  if (! (ref $exp && $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') && $exp->dbID())){
	throw('Must specify a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }

  
  #These are set in Set, just validate here
  throw ('Must provide a FeatureType') if(! defined $self->feature_type);
  throw ('Must provide a CellType') if(! defined $self->cell_type);

  my $type = $self->type;

  #should type be feature type?
  #i.e. result/annotated
  #Or do we need something more specific to segregate sequencing and other non annotated 
  #imports
  #i.e. to tell us something about the amount of data, so 
  #we can tell wether to import the 0 window_size or not
  #Or should this just be left to import params
  #Where are we capturing the differences between result InputSets
  #Not that there are any at the moment

  if(! ($type && grep /^$type$/, ('sequencing', 'feature'))){
	#Should we genericise sequencing to input
	#and specify sequencing at input_set level
	#This would allow us to capture this for direct imports into a FeatureSet
	#We would then have some overlap between input type here and input_set table_name?
	
  	throw("You must define a valid InputSet type e.g. 'sequencing' or 'feature'");
  }


  #if(! defined $self->analysis){
  ##default analysis hack for v47
  ##Set directly to avoid dbID boolean check
  #This is to support supporting_set cache in data_set?
  $self->{'analysis'} = Bio::EnsEMBL::Analysis->new(-logic_name => 'external',
													-id       => 0,#??someone needs to rewrite analysis
												   );
  #}

  $self->format($format) if defined $format;
  $self->vendor($vendor) if defined $vendor;
  $self->{'experiment'} = $exp;
  $self->{'subsets'} = {};
  
  return $self;
}


=head2 add_new_subset

  Arg [1]    : string - sub set name e.g. the file name (not path as we're restricted to 30 chars)
  Example    : $expset->add_new_subset($ss_name, $exp_subset);
  Description: Adds input_subset
  Returntype : none
  Exceptions : Throws if set is already present
               Throws if InputSubset is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub add_new_subset {
  my ($self, $ss_name, $exp_sset) = @_;
	
  #Need to test $ss_name here
  if(! ($ss_name && ref(\$ss_name) eq 'SCALAR')){#ref to $exp_sset woulf be REF
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

  $self->{'subsets'}{$ss_name} = $exp_sset;

  return $self->{'subsets'}{$ss_name};
}


=head2 get_Experiment

  Example    : my $exp = $exp_set->get_Experiment();
  Description: Getter for the Experiment of this DataSet.
  Returntype : Bio::EnsEMBL::Fuuncgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Experiment{
  my $self = shift;

  return $self->{'experiment'};  
}


=head2 get_subsets

  Example    : my @subsets = @{$exp_set->get_subsets()};
  Description: Getter for the subsets for this InputSet.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_subsets{
  my ($self)  = shift;

  return [ values %{$self->{'subsets'}} ];
}

=head2 get_subset_by_name

  Example    : my $subsets = $exp_set->get_subset_by_name('subset1');
  Description: Getter for the subset of a given name for this InputSet.
  Returntype : Bio::EnsEMBL::Funcgen::ExpeirmentalSubset
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

  Arg[1]     : string - vendor 
  Example    : my $eset->vendor('SOLEXA');
  Description: Getter/Setter for the vendor attribute of this DataSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub vendor {
  my $self = shift;
     	
  $self->{'vendor'} = shift if @_;

  return $self->{'vendor'};
}


=head2 format

  Arg[1]     : string - format i.e. product type/format
  Example    : my $eset->format('DATASET1');
  Description: Getter/Setter for the format attribute of this InputSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub format {
  my $self = shift;
     	
  $self->{'format'} = shift if @_;
  
  return $self->{'format'};
}

1;

