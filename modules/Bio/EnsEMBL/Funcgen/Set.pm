#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::Set - A module to represent a base Set object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Set;

@INC = qw (Bio::EnsEMBL::Funcgen::Set)

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);

  
}

=head1 DESCRIPTION

A base Set object which provides access common methods available across all Funcgen Set classes.


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Set;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor for Set objects.
  Returntype : Bio::EnsEMBL::Funcgen::Set
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  #do we need to add $fg_ids to this?  Currently maintaining one feature_group focus.(combi exps?)
  my ($name, $anal, $ftype, $ctype)
    = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE'], @_);
  
  throw('Need to specify a name') if ! defined $name;

  if(defined $anal){#Move this to child Sets, and just set anal here
	$self->analysis($anal);
  }elsif($self->set_type ne 'experimental'){
	throw('Must pass a valid -analysis parameter for a '.ref($self));
  }

  $self->{'name'} = $name;
  $self->cell_type($ctype) if $ctype;
  $self->feature_type($ftype) if $ftype;

  #Set type here but don't validate, as this can be done selectively in the Child Sets, do this for analysis too?
  
  return $self;
}






=head2 name

  Example    : my $set->name('SET1');
  Description: Getter/Setter for the name of this Set.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub name {
  my $self = shift;

  return $self->{'name'};
}

=head2 cell_type

  Example    : my $dset_ctype_name = $dset->cell_type->name();
  Description: Getter for the cell_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : throws if arg not valid
  Caller     : General
  Status     : At Risk

=cut

sub cell_type {
  my ($self, $ctype) = @_;

  if(defined $ctype){

	if(! (ref($ctype) eq 'Bio::EnsEMBL::Funcgen::CellType'
		  && $ctype->dbID())){ 
	  throw('Must pass a valid stored Bio::EnsEMBL::Funcgen::CellType');
	}
	$self->{'cell_type'} = $ctype;
  }

  return $self->{'cell_type'};
}

=head2 feature_type

  Example    : my $dset_ftype_name = $dset->feature_type->name();
  Description: Getter for the feature_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if arg not valid
  Caller     : General
  Status     : At Risk

=cut

sub feature_type {
  my ($self, $ftype) = @_;
   
  if(defined $ftype){

	if(! (ref($ftype) eq 'Bio::EnsEMBL::Funcgen::FeatureType'
		  && $ftype->dbID())){ 
	  throw('Must pass a valid stored Bio::EnsEMBL::Funcgen::FeatureType');
	}
	$self->{'feature_type'} = $ftype;
  }

  		  
  return $self->{'feature_type'};
}

=head2 analysis

  Example    : my $anal_name = $set->analysis->logic_name();
  Description: Getter for the analysis attribute for this Set.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub analysis {
  my $self = shift;

  if(@_){
	throw('Must pass a valid stored Analysis') if (! (ref($_[0]) eq 'Bio::EnsEMBL::Analysis'
													  && $_[0]->dbID()));
	$self->{'analysis'} = shift;
  }
  
 
  return $self->{'analysis'};
}

=head2 display_label

  Example    : print $set->display_label();
  Description: Getter for the display_label attribute for this Set.
               This is more appropriate for teh predicted_features of the set.
               Use the individual display_labels for each raw result set.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;

  throw('Not yet implemented in the Base Set class');

  #Add display label in table?
  #Can we aborc ResultSet method into this?

  if(! $self->{'display_label'}){

	if($self->FeatureSet->feature_type->class() eq 'REGULATORY FEATURE'){
	  $self->{'display_label'} = 'Regulatory Features';
	}
	else{

	  $self->{'display_label'} = $self->feature_type->name()." -";
	  $self->{'display_label'} .= " ".($self->cell_type->display_label() || 
									   $self->cell_type->description()   ||
									   $self->cell_type()->name());
	  $self->{'display_label'} .= " Enriched Sites";
	}
  }
 
  return $self->{'display_label'};
}



=head2 set_type

  Example    : my $set_type = $set->set_type;
  Description: Getter for the Set type for this Set.
  Returntype : string e.g. result, feature, data or experimental
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub set_type {
  my $self = shift;
 
  my $type;
  my @namespace = split/\:\:/, ref($self);
  ($type = lc($namespace[$#namespace])) =~ s/set//;

  return $type;
}

1;

