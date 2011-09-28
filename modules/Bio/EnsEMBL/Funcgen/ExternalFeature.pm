#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExternalFeature
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

Bio::EnsEMBL::ExternalFeature - A module to represent an externally curated feature 
mapping from an external_db.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExternalFeature;

my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
    -DISPLAY_LABEL => $text,
    -FEATURE_SET   => $fset,
    -FEATURE_TYPE  => $ftype,
);



=head1 DESCRIPTION

An ExternalFeature object represents the genomic placement of an externally curated
feature from and DB external to Ensembl.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ExternalFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::SetFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::SetFeature);


=head2 new

 
  Arg [-FEATURE_SET]  : Bio::EnsEMBL::Funcgen::FeatureSet
  Arg [-FEATURE_TYPE] : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-START]        : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]: string - Display label for this feature
  Arg [-STRAND]       : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int - Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.
  Example             : my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
										                                  -SLICE         => $chr_1_slice,
									                                      -START         => 1_000_000,
                      				                                      -END           => 1_000_024,
									                                      -STRAND        => -1,
									                                      -DISPLAY_LABEL => $text,
									                                      -FEATURE_SET   => $fset,
                                                                          -FEATURE_TYPE  => $ftpe,
                                                                         );


  Description: Constructor for ExternalFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::ExternalFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  #Remove this method if we interdb_stable_id to SetFeature
  ($self->{'interdb_stable_id'}) = rearrange(['INTERDB_STABLE_ID'], @_);
 		
  return $self;
}

=head2 interdb_stable_id

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $idb_sid = $feature->interdb_stable_id();
  Description: Getter for the interdb_stable_id attribute for this feature.
               This is simply to avoid using internal db IDs for inter DB linking
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub interdb_stable_id {
  return $_[0]->{'interdb_stable_id'};
}




=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium risk

=cut

sub display_label {
  my $self = shift;
	    
  if(! $self->{'display_label'}  && $self->adaptor){
	
	$self->{'display_label'} = $self->feature_set->feature_type->name().' - ';
	$self->{'display_label'} .= $self->cell_type->name() if $self->cell_type();
	$self->{'display_label'} .= $self->feature_type->name() if(defined $self->{'feature_type'});
  }
	
  return $self->{'display_label'};
}



1;

