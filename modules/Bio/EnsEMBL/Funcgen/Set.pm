#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::Set - A module to represent a base Set object.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Set;

my $data_set = Bio::EnsEMBL::Funcgen::Set->new(
	                                              -DBID            => $dbID,
							 					  -ADAPTOR         => $self,
                                                  -ANALYSIS        => $anal,
                                                  -FEATURE_TYPE    => $ftype,
                                                  -CELL_TYPE       => $ctype,
                                                  -DISPLAYABLE     => 1,
                                                  -NAME            => 'SET1',
                                               );



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
  my ($ftype, $ctype, $name, $analysis)
    = rearrange(['FEATURE_TYPE', 'CELL_TYPE', 'NAME', 'ANALYSIS'], @_);
  

  if (! ($ftype && $anal && $name)){
    throw('Need to specify a name, feature_type and analysis parameters');
  }
  


  #this clashes with Data::Set->product_feature_type
  #do we @INC Set in DataSet?
  #mandatory params would also be different
  #keep DataSet separate for now.
  
  $self->feature_type($ftype);
  $self->cell_type($ctype)     if $ctype;
  $self->analysis($anal)      if $anal;
  $self->name($name);	
  
  return $self;
}






=head2 name

  Example    : my $set->name('SET1');
  Description: Getter/Setter for the name of this DataSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub name {
  my $self = shift;
     	
  $self->{'name'} = shift if @_;

  return $self->{'name'};
}

=head2 cell_type

  Example    : my $dset_ctype_name = $dset->cell_type->name();
  Description: Getter for the cell_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub cell_type {
  my $self = shift;

  if(@_){
	throw('Must pass a valid stored CellType') if (! (ref($_[0]) eq 'Bio::EnsEMBL::Funcgen::CellType'
													  && $_[0]->dbID()));
	$self->{'cell_type'} = shift;
  }

  return $self->{'cell_type'};
}

=head2 feature_type

  Example    : my $dset_ftype_name = $dset->feature_type->name();
  Description: Getter for the feature_type for this DataSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_type {
  my $self = shift;
     		
  if(@_){
   throw('Must pass a valid stored FeatureType') if (! (ref($_[0]) eq 'Bio::EnsEMBL::Funcgen::FeatureType'
														&& $_[0]->dbID()));
   $self->{'feature_type'} = shift;
 }
  
  return $self->{'feature_type'};
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





1;

