#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::MirnaTargetFeatureAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::DBSQL::Funcgen::MirnaTargetFeatureAdaptor - A database adaptor for fetching and
storing MirnaTargetFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_ExternalFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The MirnaTargetFeatureAdaptor is a database adaptor for storing and retrieving
MirnaTargetFeature objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::MirnaTargetFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::MirnaTargetFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor; #DBI sql_types import

use base qw( Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor );


=head2 fetch_all_by_gene_stable_id

    Arg [1]    : String $gene_stable_id - The stable id of the linked gene
    Example    : my $mirna_target_features =
    mtf_adaptor->fetch_all_by_gene_stable_id('ENSG00000000001');
    Description: Retrieves a list of mirna target features that are linked to the
                 given Gene.
    Returntype : arrayref of Bio::EnsEMBL::Funcgen::MirnaTargetFeature objects
    Exceptions : none
    Caller     : general
    Status     : Stable

=cut

sub fetch_all_by_gene_stable_id {
    my $self      = shift;
    my $gene_stable_id = shift;
    throw('Must specify a gene_stable_id') if !defined $gene_stable_id;
    my $constraint = "mitaf.gene_stable_id = ?";
    $self->bind_param_generic_fetch($gene_stable_id, SQL_VARCHAR);

    my $mirna_target_features = $self->generic_fetch($constraint);
    return $mirna_target_features;
}

=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return ([ 'mirna_target_feature', 'mitaf' ]);
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  return qw(
			mitaf.mirna_target_feature_id
      mitaf.seq_region_id
			mitaf.seq_region_start
      mitaf.seq_region_end
			mitaf.seq_region_strand
      mitaf.display_label
			mitaf.feature_type_id
      mitaf.analysis_id
      mitaf.gene_stable_id
      mitaf.accession
      mitaf.evidence
      mitaf.method
      mitaf.supporting_information
	   );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates MirnaTargetFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::MirnaTargetFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper) = @_;

  # Remap the feature coordinates to another coord system if a mapper was provided
  return [] if ($mapper); 
	# For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	# So if it not defined then we need to generate one derived from the species_name and
  # schema_build of the feature we're retrieving.


  my (@features, %anal_hash, %slice_hash, %sr_name_hash, %ftype_hash);
  my $slice_adaptor ||= $self->db->dnadb->get_SliceAdaptor;

    my $anal_adaptor  = $self->db->get_AnalysisAdaptor();
    my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor();
    my (
      $external_feature_id,
      $seq_region_id,
      $seq_region_start,
      $seq_region_end,
      $seq_region_strand,
      $anal_id,
      $gene_stable_id,
      $display_label,
      $ftype_id,
      $accession,
      $evidence,
      $method,
      $supporting_information,
    );

    $sth->bind_columns(
        \$external_feature_id,
        \$seq_region_id,
        \$seq_region_start,
        \$seq_region_end,
        \$seq_region_strand,
        \$display_label,
        \$ftype_id,
        \$anal_id,
        \$gene_stable_id,
        \$accession,
        \$evidence,
        \$method,
        \$supporting_information,
      );

  FEATURE: while ( $sth->fetch() ) {

	  #Get the FeatureSet/Type objects
	  $anal_hash{$anal_id}   = $anal_adaptor->fetch_by_dbID($anal_id)   if(! exists $anal_hash{$anal_id});
	  $ftype_hash{$ftype_id} = $ftype_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});
    $slice_hash{$seq_region_id} = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
	  # Get the slice object

	  push @features, Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new_fast
		({
		  'start'               => $seq_region_start,
		  'end'                 => $seq_region_end,
		  'strand'              => $seq_region_strand,
		  'slice'               => $slice_hash{$seq_region_id},
		  'analysis'            => $anal_hash{$anal_id},
      'gene_stable_id'      => $gene_stable_id,
		  'adaptor'             => $self,
		  'dbID'                => $external_feature_id,
		  'display_label'       => $display_label,
		  'feature_type'        => $ftype_hash{$ftype_id},
      'accession'           => $accession,
      'evidence'            => $evidence,
      'method'              => $method,
      'supporting_information' => $supporting_information,
		 });
	  }

	return \@features;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::MirnaTargetFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given MirnaTargetFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : Arrayref of stored ExternalFeatures
  Exceptions : Throws if a list of MirnaTargetFeature objects is not provided or if
               the Analysis, CellType and FeatureType objects are not attached or stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @mrna_fs) = @_;

	if (scalar(@mrna_fs) == 0) {
		throw('Must call store with a list of MirnaTargetFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO mirna_target_feature (
			seq_region_id,
      seq_region_start,
			seq_region_end,
      seq_region_strand,
      display_label,
      feature_type_id,
      analysis_id,
      gene_stable_id,
      accession,
      evidence,
      method,
      supporting_information
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();

  FEATURE: foreach my $mrna_f (@mrna_fs) {

	  if(! (ref($mrna_f) && $mrna_f->isa('Bio::EnsEMBL::Funcgen::MirnaTargetFeature') )) {
  		throw('Feature must be an MirnaTargetFeature object');
	  }

	  if ( $mrna_f->is_stored($db) ) {
  		#does not accomodate adding Feature to >1 feature_set
  		warning('MirnaTargetFeature [' . $mrna_f->dbID() . '] is already stored in the database');
  		next FEATURE;
	  }

	  my $seq_region_id;
	  ($mrna_f, $seq_region_id) = $self->_pre_store($mrna_f);

	  $sth->bind_param(1,  $seq_region_id,                  SQL_INTEGER);
	  $sth->bind_param(2,  $mrna_f->start,                  SQL_INTEGER);
	  $sth->bind_param(3,  $mrna_f->end,                    SQL_INTEGER);
	  $sth->bind_param(4,  $mrna_f->strand,                 SQL_TINYINT);
	  $sth->bind_param(5,  $mrna_f->display_label,          SQL_VARCHAR);
	  $sth->bind_param(6,  $mrna_f->get_FeatureType->dbID,     SQL_INTEGER);
    $sth->bind_param(7,  $mrna_f->analysis->dbID,         SQL_INTEGER);
	  $sth->bind_param(8,  $mrna_f->gene_stable_id,         SQL_VARCHAR);
    $sth->bind_param(9,  $mrna_f->accession,              SQL_VARCHAR);
    $sth->bind_param(10, $mrna_f->evidence,               SQL_VARCHAR);
    $sth->bind_param(11, $mrna_f->method,                 SQL_VARCHAR);
    $sth->bind_param(12, $mrna_f->supporting_information, SQL_VARCHAR);

	  $sth->execute();
	  $mrna_f->dbID( $self->last_insert_id );
	  $mrna_f->adaptor($self);
	}

  return \@mrna_fs;
}


sub fetch_all_by_accessions{
  my ($self, $accessions) = @_;
  my $params = {constraints => {accession=>$accessions}};
  return $self->generic_fetch($self->compose_constraint_query($params));

}

sub _constrain_accession {
  my ($self, $accessions) = @_;
  my $constraint = 'mitaf.accession IN(' . join (', ', map { qq!"$_"! } @$accessions).')';
  return $constraint;

}

1;
