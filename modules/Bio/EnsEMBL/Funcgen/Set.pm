#
# Ensembl module for Bio::EnsEMBL::Funcgen::Set
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
  
  #TYPE was never parsed here?
  #Only in inheritants that used it i.e. FeatureSet

  my ($name, $anal, $ftype, $ctype, $set_type, $fclass, $type)
    = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE', 'SET_TYPE', 'FEATURE_CLASS', 'TYPE'], @_);
  
  throw('Need to specify a name') if ! defined $name;

  $self->set_type($set_type);
  $self->feature_class($fclass);
  $self->feature_class($type) if $type;#Remove this when fully implemented
  $self->{'name'} = $name;
  $self->cell_type($ctype) if $ctype;
  $self->feature_type($ftype) if $ftype;

  if(defined $anal){
	$self->analysis($anal);
  }elsif($self->set_type ne 'input'){
	#Could move this to child Sets and just set analysis here
	#As with ftype
	throw('Must pass a valid -analysis parameter for a '.ref($self));
  }

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


=head2 feature_class

  Arg[0]     : string - feature class e.g. result, annotated, regulatory or external.
  Example    : my $fclass = $dset->feature_class;
  Description: Getter for the feature_type for this Set.
  Returntype : string 
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#Supercededs type method in FeatureSet

sub feature_class {
  my ($self, $fclass) = @_;
   
  if(defined $fclass){

	#Leave this an implement in inheritants
	#if(! grep /^${fclass}$/, ('annotated', 'result', 'external', 'regulatory')){
	#  throw("You have no supplied a valid feature class:\t$fclass");
	#}

	$self->{'feature_class'} = $fclass;
  }

  return $self->{'feature_class'};
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

 
  #Add display label in table?
  #Can we aborc ResultSet method into this?

  if(! $self->{'display_label'}){

	#if($self->product_FeatureSet->feature_type->class() eq 'Regulatory Feature'){
	#  $self->{'display_label'} = 'Regulatory Features';
	#}
	#else{

	#This only works for annotated/regulatory_feature sets and result sets
	#Move to other Set classes?

	$self->{'display_label'} = $self->feature_type->name()." -";
	$self->{'display_label'} .= " ".($self->cell_type->display_label() || 
									 $self->cell_type->description()   ||
									 $self->cell_type()->name());
	

	if($self->set_type eq 'result'){
	  $self->{'display_label'} .= " signal";
	}
	else{
	  $self->{'display_label'} .= " enriched sites";
	}
  }
 
  return $self->{'display_label'};
}



=head2 set_type

  Example    : my $set_type = $set->set_type;
  Description: Getter for the Set type for this Set.
  Returntype : string e.g. result, feature, data, input
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub set_type {
  my ($self, $set_type) = @_;
   
  if(defined $set_type){
	$self->{'_set_type'} = $set_type;
  }
  elsif(! defined $self->{'_set_type'}){
	my @namespace = split/\:\:/, ref($self);
	($self->{'_set_type'} = lc($namespace[$#namespace])) =~ s/set//;
	
  }

  return $self->{'_set_type'};
}

=head2 type

  Example    : my $type = $set->type;
  Description: Getter for the type for this Set.
               e.g. annotated, external, regulatory for FeatureSets
                    or 
                    array, sequencing for InputSets
               Currently not applicable to DataSets or ResultSets
  Exceptions : None
  Returntype : string 
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub type {
  my $self = shift;
   
  deprecate("Please use feature_class instead");
  
  return $self->feature_class(@_);

  #$self->{'feature_class'} = shift if @_;
  
  #return $self->{'feature_class'};
}



1;

