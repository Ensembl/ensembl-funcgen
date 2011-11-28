#
# Ensembl module for Bio::EnsEMBL::Funcgen::FeatureSet
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

Bio::EnsEMBL::FeatureSet - A module to represent FeatureSet.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::FeatureSet;

my $result_set = Bio::EnsEMBL::Funcgen::FeatureSet->new(

); 



=head1 DESCRIPTION

A FeatureSet object provides access to a set of feature predictions and their details, which may have been generated from a 
single or multiple Experiments with potentially differing analyses.  The FeatureSet itself will only have a single analysis 
which may be one or a combination of programs but will be represented by one analysis record.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::FeatureSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate);
use Bio::EnsEMBL::Funcgen::Set;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);

my %valid_classes = (
					 annotated    => undef,
					 regulatory   => undef,
					 external     => undef,
					 segmentation => undef,
					);

=head2 new

                           -name          => $name,
                                                                    -feature_type  => $ftype,
                                                                    -cell_type     => $ctype,
                                                                    -name          => $name,
              -description   => 'Release 3.1',
                                                           -display_label => 'Short name',
  -analysis      => $analysis,
  Arg [-EXPERIMENT_ID]     : Experiment dbID
  -dbid          => $dbid,
Arg [-ADAPTOR]

  Example    : my $feature = Bio::EnsEMBL::Funcgen::FeatureSet->new(
                                                                    -dbid          => $dbid,
                                                                    -analysis      => $analysis,
                                                                    -feature_type  => $ftype,
                                                                    -cell_type     => $ctype,
                                                                    -name          => $name,
                                                                    -feature_class => 'annotated',
                                                                    -description   => 'Release 3.1',
                                                                    -display_label => 'Short name',
                                                                    -experiment_id => $exp_id,
			                                                       ); 
  Description: Constructor for FeatureSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : Throws if FeatureType defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;	
  my $self = $class->SUPER::new(@_);
	
  my ($desc, $dlabel, $exp_id, $exp)
    = rearrange(['DESCRIPTION', 'DISPLAY_LABEL', 'EXPERIMENT_ID', 'EXPERIMENT'],@_);

  #Allow exp or exp_id to be passed to support storing and lazy loading

  #Mandatory params checks here (setting done in Set.pm)
  throw ('Must provide a FeatureType') if(! defined $self->feature_type);

  #explicit type check here to avoid invalid types being imported as NULL
  #subsequently throwing errors on retrieval
  my $type = $self->feature_class;

  if(! ($type && exists $valid_classes{$type}) ){
	throw('You must define a valid FeatureSet type e.g. '.
		  join(', ', keys %valid_classes));
  }

  #Direct assignment to prevent need for set arg test in method

  $self->{'description'}   = $desc   if defined $desc;
  $self->{'display_label'} = $dlabel if defined $dlabel;
  $self->{'experiment_id'} = $exp_id if defined $exp_id;

  if(defined $exp){
	#Exp obj is only passed during object storing
	#so let the adaptor do is_stored_and_valid
	$self->{'experiment'}    = $exp;
  }
  
  return $self;
}


=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  return bless ($_[1], $_[0]);
}


=head2 description

  Example    : print "Feature set description is:\t".$fset->description."\n";
  Description: Getter for the description of this FeatureSet. e.g. Release 3.1
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description {
  return $_[0]->{'description'};
}



=head2 display_label

  Example    : print $rset->display_label;
  Description: Getter for the display_label attribute for this FeatureSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;
  
  if(! $self->{'display_label'}){

	if($self->feature_type->class() eq 'Regulatory Feature'){
	  $self->{'display_label'} = $self->name;
	}
	else{
	  #This still fails here if we don't have a class or a cell_type set
	  
	  $self->{'display_label'} = $self->feature_type->name()." - ".$self->cell_type->name()." Enriched Sites";
	}
  }
	
  return $self->{'display_label'};
}



=head2 get_FeatureAdaptor

  Example    : 
  Description: Retrieves and caches FeatureAdaptor of feature_set type 
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::ucfirst($self->feature_class())FeatureAdaptor
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_FeatureAdaptor{
  my $self = shift;

  if(! exists $self->{'adaptor_refs'}){

	foreach my $valid_class(keys %valid_classes){
	  my $method = 'get_'.ucfirst($valid_class).'FeatureAdaptor';

	  $self->{'adaptor_refs'}{$valid_class} =  $self->adaptor->db->$method;
	}
  }
  
  return $self->{'adaptor_refs'}->{$self->feature_class()};

}



=head2 get_Features_by_Slice

  Example    : my @features = @{$FeatureSet->get_Features_by_Slice($slice)};
  Description: Retrieves all Features for this FeatureSet for a given Slice
  Returntype : ARRAYREF containing Features of the feature_set type i.e. Annotated, Regulatory or Supporting;
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Features_by_Slice{
  my ($self, $slice) = @_;

  return $self->get_FeatureAdaptor->fetch_all_by_Slice_FeatureSets($slice, [$self]);
}


=head2 get_Features_by_FeatureType

  Arg[0]     : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my @features = @{$FeatureSet->get_Features_by_FeatureType($ftype)};
  Description: Retrieves all Features for this FeatureSet for a given FeatureType
               or associated FeatureType. This is mainly used by external FeatureSets
               which can sometimes have more than one associated FeatureType.
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_Features_by_FeatureType{
  my ($self, $type) = @_;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureType_FeatureSets($type, [$self]);
}


=head2 get_all_Features

  Example    : my @features = @{$FeatureSet->get_all_Features};
  Description: Retrieves all Features for this FeatureSet
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_all_Features{
  my $self = shift;

  return $self->get_FeatureAdaptor->fetch_all_by_FeatureSets([$self]);
}




=head2 is_focus_set

  Args       : None
  Example    : if($fset->is_focus_set){ ... }
  Description: Returns true if FeatureSet is a focus set used in the RegulatoryBuild
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

sub is_focus_set{
  my $self = shift;

  if(! defined $self->{focus_set}){

	if(! defined $self->cell_type){
	  warn "FeatureSet without an associated CellType cannot be a focus set:\t".$self->name;
	  $self->{focus_set} = 0;
	}
	else{
	   $self->{focus_set} = $self->adaptor->fetch_focus_set_config_by_FeatureSet($self);
	 }
  }

  return $self->{focus_set};
}


=head2 is_attribute_set

  Args       : None
  Example    : if($fset->is_attribute_set){ ... }
  Description: Returns true if FeatureSet is a supporting/attribute(focus or not) set used in the RegulatoryBuild
  Returntype : Boolean
  Exceptions : Throws if meta entry not present
  Caller     : General
  Status     : At Risk

=cut

sub is_attribute_set{
  my $self = shift;

  if(! defined $self->{attribute_set}){

	if(! defined $self->cell_type){
	  warn "FeatureSet without an associated CellType cannot be a attribute set:\t".$self->name;
	  $self->{attribute_set} = 0;
	}
	else{
	   $self->{attribute_set} = $self->adaptor->fetch_attribute_set_config_by_FeatureSet($self);
	 }
  }

  return $self->{attribute_set};
}


=head2 get_Experiment

  Example    : my $exp = $FeatureSet->get_Experiment;
  Description: Retrieves the Experiment for this FeatureSet
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_Experiment{
  my $self = shift;

  if( (! defined $self->{'experiment'}) &&
	  (defined $self->{'experiment_id'}) ){
	$self->{'experiment'} = $self->adaptor->db->get_ExperimentAdaptor->fetch_by_dbID($self->{experiment_id});
  }


  return $self->{'experiment'};
}


=head2 source_label

  Example    : my $source_label = $fset->source_label;
  Description: Retrieves the source label this FeatureSet, used in zmenus
  Returntype : String
  Exceptions : None
  Caller     : Webcode
  Status     : At Risk - remove, to be done by webcode?

=cut

sub source_label{
  my $self = shift;


  if(! defined $self->{'source_label'}){
	my $exp          = $self->get_Experiment;
	my $source_label;

	if($exp){

	  $source_label = $exp->archive_id;
	  
	  my $exp_group = $exp->experimental_group;
	
	  if($exp_group && 
		 $exp_group->is_project){
		
		if($source_label){
		  $source_label .= ' '.$exp_group->name;
		}
		else{
		  $source_label = $exp_group->name;
		}
	  }
	}

	$self->{'source_label'} = $source_label || '';
  }

  return $self->{'source_label'};
}


1;
