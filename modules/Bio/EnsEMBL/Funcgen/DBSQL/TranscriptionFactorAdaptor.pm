#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor
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

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;    #sql_types barewords import
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::TranscriptionFactor;

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string, name of the transcription factor
  Example    : my $tf = tf_adaptor->fetch_by_name($name);
  Description: Fetches TranscriptionFactor object given it's name
  Returntype : Bio::EnsEMBL::Funcgen::TranscriptionFactor
  Exceptions : Throws if name is not specified
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name {
    my ( $self, $name ) = @_;

    throw('Must specify a name') if !defined $name;

    my $constraint = " tf.name = ?";

    $self->bind_param_generic_fetch( $name, SQL_VARCHAR );

    my $result = $self->generic_fetch($constraint);

    return $result->[0];
}

=head2 fetch_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : my $tf = tf_adaptor->fetch_by_FeatureType($ft);
  Description: Fetches TranscriptionFactor object
  Returntype : Bio::EnsEMBL::Funcgen::TranscriptionFactor
  Exceptions : Throws if no feature_type paramater is defined
               Throws if parameter passed is not a
                 Bio::EnsEMBL::Funcgen::FeatureType object
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_FeatureType {
    my ($self, $feature_type) = @_;

    throw('Must pass a FeatureType parameter') if !defined $feature_type;
    assert_ref($feature_type,
               'Bio::EnsEMBL::Funcgen::FeatureType',
               'FeatureType');

    my $constraint = " tf.feature_type_id = ?";
    $self->bind_param_generic_fetch($feature_type->dbID(), SQL_INTEGER);
    my $result = $self->generic_fetch($constraint);

    return $result->[0];
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
    return ( [ 'transcription_factor', 'tf' ] );
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
    return
        qw( tf.transcription_factor_id tf.name tf.feature_type_id tf.gene_stable_id
    );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates objects from an executed DBI statement handle.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
    my ( $self, $sth ) = @_;

    my ( @result, $transcription_factor_id, $name, $feature_type_id,
        $gene_stable_id );

    $sth->bind_columns( \$transcription_factor_id, \$name, \$feature_type_id,
        \$gene_stable_id );

    my $feature_type_adaptor = $self->db->get_adaptor('FeatureType');
    my %feature_type_cache;

    while ( $sth->fetch() ) {

        my $feature_type;
# TODO test if bugfix works as expected
        if ($feature_type_id) {
            if (!exists $feature_type_cache{$feature_type_id}) {
                $feature_type_cache{$feature_type_id}
                    = $feature_type_adaptor->fetch_by_dbID($feature_type_id);
            }
            $feature_type = $feature_type_cache{$feature_type_id};
        }

        my $transcription_factor
            = Bio::EnsEMBL::Funcgen::TranscriptionFactor->new(
            -dbID           => $transcription_factor_id,
            -NAME           => $name,
            -FEATURE_TYPE   => $feature_type,
            -GENE_STABLE_ID => $gene_stable_id,
            -ADAPTOR        => $self,
            );

        push @result, $transcription_factor;

    }

    return \@result;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
  Example    : $tf_adaptor->store($tf1, $tf2, $tmf3);
  Description: Stores given TranscriptionFactor objects in the database.
         Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @args = @_;

    my $sth = $self->prepare( "
      INSERT INTO transcription_factor
      (name, feature_type_id, gene_stable_id)
      VALUES (?, ?, ?)" );

    my $stored_TranscriptionFactor;

    foreach my $transcription_factor (@args) {
        assert_ref( $transcription_factor,
            'Bio::EnsEMBL::Funcgen::TranscriptionFactor',
            'TranscriptionFactor' );

        my $feature_type
            = $transcription_factor->get_FeatureType();
            my $feature_type_id;

        if ($feature_type) {
            $self->db->is_stored_and_valid(
                'Bio::EnsEMBL::Funcgen::FeatureType',
                $feature_type
            );
            $feature_type_id= $feature_type->dbID();
        }

        if (!(     $transcription_factor->dbID()
                && $transcription_factor->adaptor() == $self
            )
            )
        {
            #Check for previously stored TranscriptionFactors
            $stored_TranscriptionFactor
                = $self->fetch_by_name( $transcription_factor->name );

            if ( !$stored_TranscriptionFactor ) {
                $sth->bind_param( 1, $transcription_factor->name(),
                    SQL_VARCHAR );
                $sth->bind_param(
                    2,
                    $feature_type_id,
                    SQL_INTEGER
                );
                $sth->bind_param( 3, $transcription_factor->gene_stable_id(),
                    SQL_VARCHAR );

                $sth->execute();
                $transcription_factor->dbID( $self->last_insert_id );
                $transcription_factor->adaptor($self);
            }
            else {
                $transcription_factor = $stored_TranscriptionFactor;
                warn(     "Using previously stored TranscriptionFactor:\t"
                        . $transcription_factor->name()
                        . "\n" );
            }
        }
    }

    return \@args;
}


1;