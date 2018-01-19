#
# Ensembl module for Bio::EnsEMBL::Funcgen::MirnaTargetFeature
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::MirnaTargetFeature - A module to represent an externally curated feature
mapping from an external_db.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MirnaTargetFeature;

my $feature = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
    -DISPLAY_LABEL => $text,
    -FEATURE_SET   => $fset,
    -FEATURE_TYPE  => $ftype,
);



=head1 DESCRIPTION

An MirnaTargetFeature object represents the genomic placement of an externally curated
feature from and DB external to Ensembl.

=cut

package Bio::EnsEMBL::Funcgen::MirnaTargetFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw(Bio::EnsEMBL::Funcgen::SetFeature);


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
  Example             : my $feature = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
                            -SLICE         => $chr_1_slice,
                            -START         => 1_000_000,
                            -END           => 1_000_024,
                            -STRAND        => -1,
                            -DISPLAY_LABEL => $text,
                            -FEATURE_SET   => $fset,
                            -FEATURE_TYPE  => $ftpe,

                                               );


  Description: Constructor for MirnaTargetFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::MirnaTargetFeature
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($accession, $evidence, $interdb_stable_id, $method, $supporting_information ) = 
    rearrange (['ACCESSION', 'EVIDENCE', 'INTERDB_STABLE_ID', 'METHOD', 'SUPPORTING_INFORMATION'], @_);

  if(! defined $accession){
    throw("Mandatory parameter -accession not defined");
  }
  $self->{accession}              = $accession;
  $self->{evidence}               = $evidence;
  $self->{method}                 = $method;
  $self->{supporting_information} = $supporting_information;

  #Remove this method if we interdb_stable_id to SetFeature
  $self->{'interdb_stable_id'} = $interdb_stable_id;

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
  Status     : At Risk - might move to SetFeature

=cut

sub interdb_stable_id {
  return $_[0]->{'interdb_stable_id'};
}

=head2 get_all_DBEntries

  Arg[1]     : string - External DB name e.g. ensembl_core_Gene
  Arg[2]     : string - External DB type
  Example    : my @dbentries = @{ $set_feature->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this SetFeature.
               This does _not_ include the corresponding translations
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the SetFeature (i.e. they have not already been added or
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general, get_all_DBLinks
  Status     : Stable - at risk move to storable

=cut


#We could add 3rd arg here which would be xref(info_)type e.g. Gene/Transcript etc.
#Move info_type to ox.linkage_type to sit along side linkage_annotated


sub get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  #Need to add tests for valid objects for xrefs

  # if not cached, retrieve all of the xrefs for this gene

  #This is not using the caching optimally
  #It seems for naive(ex_db_exp,ex_db_type) queries we create a naive cache
  #This means that further more specific queries will make another query and not use the cache


  if( (! defined $self->{$cache_name}) && $self->adaptor() ){

	my @tables = $self->adaptor->_tables;
	@tables = split/_/, $tables[0]->[0];#split annotated_feature
	my $object_type = join('', (map ucfirst($_), @tables));#change to AnnotatedFeature

    $self->{$cache_name} =
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
  }
  elsif( ! defined $self->{$cache_name} ){
	throw('You must have set and adaptor to be able to get_all_DBEntries');
  }


  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
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

	$self->{'display_label'}  = $self->feature_set->feature_type->name().' - ';
	$self->{'display_label'} .= $self->epigenome->name() if $self->epigenome();
	$self->{'display_label'} .= $self->feature_type->name() if(defined $self->{'feature_type'});
  }

  return $self->{'display_label'};
}

=head2 accession

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->accession();
  Description: Getter for the accession attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub accession {
  return $_[0]->{'accession'};
}

=head2 evidence

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->evidence();
  Description: Getter for the evidence attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub evidence {
  return $_[0]->{'evidence'};
}

=head2 method

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->method();
  Description: Getter for the method attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub method {
  return $_[0]->{'method'};
}

=head2 supporting_information

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->supporting_information();
  Description: Getter for the supporting_information attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub supporting_information {
  my ($self, $info) = @_;

  if(defined $info){
    $self->{supporting_information} = $info;
  
  }
  return $self->{'supporting_information'};
}


1;

