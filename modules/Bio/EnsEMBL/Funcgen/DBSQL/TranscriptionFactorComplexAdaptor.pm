#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor
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

package Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;    #sql_types barewords import
use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor;
use Bio::EnsEMBL::Funcgen::TranscriptionFactor;
use Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex;

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_all_by_TranscriptionFactor

  Arg [1]    : Bio::EnsEMBL::Funcgen::TranscriptionFactor
  Example    : my @tfc = @{$tfc_adaptor->fetch_all_by_TranscriptionFactor($tf)};
  Description: Fetches all TrancriptionFactorComplex objects that contain the
               given TranscriptionFactor
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TrancriptionFactorComplex objects
  Exceptions : Throws if TranscriptionFactor is not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_TranscriptionFactor {
    my ( $self, $transcription_factor ) = @_;

    $self->db->is_stored_and_valid(
        'Bio::EnsEMBL::Funcgen::TranscriptionFactor',
        $transcription_factor );

    my $sth = $self->prepare(
        "SELECT transcription_factor_complex_id
         FROM transcription_factor_complex
            JOIN transcription_factor_complex_composition
                USING (transcription_factor_complex_id)
         WHERE transcription_factor_id=?"
    );

    $sth->execute( $transcription_factor->dbID );

    my @dbIDs;
    while ( my @row = $sth->fetchrow_array ) {
        push @dbIDs, $row[0];
    }
    my @transcription_factor_complexes = $self->fetch_all_by_dbID_list(\@dbIDs);

    return \@transcription_factor_complexes;
}

=head2 fetch_by_production_name

  Arg [1]    : string, production name of the transcription factor complex
  Example    : my $tfc = tfc_adaptor->fetch_by_production_name($name);
  Description: Fetches TranscriptionFactorComplex object given it's production name
  Returntype : Bio::EnsEMBL::TranscriptionFactorComplex
  Exceptions : Throws if name is not specified
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_production_name {
    my ( $self, $production_name ) = @_;

    throw('Must specify a production name') if !defined $production_name;

    my $constraint = " tfc.production_name = ?";

    $self->bind_param_generic_fetch( $production_name, SQL_VARCHAR );

    my $result = $self->generic_fetch($constraint);

    return $result->[0];
}

=head2 fetch_all_by_BindingMatrix

  Arg [1]    : Bio::EnsEMBL::Funcgen::BindingMatrix
  Example    : my $transcription_factor_complexes =
             :     tfc_adaptor->fetch_by_BindingMatrix($binding_matrix);
  Description: Fetches all TranscriptionFactorComplex objects associated with
             : a given BindingMatrix
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex objects
  Exceptions : Throws if BindingMatrix is not specified
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_BindingMatrix {
    my ( $self, $binding_matrix ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    my $sth = $self->prepare(
        "SELECT transcription_factor_complex_id FROM
         binding_matrix_transcription_factor_complex WHERE binding_matrix_id=?"
    );

    $sth->execute( $binding_matrix->dbID );

    my @dbIDs;
    while ( my @row = $sth->fetchrow_array ) {
        push @dbIDs, $row[0];
    }
    my $transcription_factor_complexes = $self->fetch_all_by_dbID_list(\@dbIDs);

    return $transcription_factor_complexes;
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
    return (
        [ 'transcription_factor_complex', 'tfc', ]
    );
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
    return qw( tfc.transcription_factor_complex_id
               tfc.production_name
               tfc.display_name
    );
}

=head2 _fetch_components

  Arg [1]    : Integer, transcription factor complex id
  Example    : None
  Description: Fetches a list of TranscriptionFactor objects (components)
               that belong to the transcription factor complex
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _fetch_components {
    my ( $self, $transcription_factor_complex_id ) = @_;

    my $sth = $self->prepare( "
      SELECT transcription_factor_id FROM
      transcription_factor_complex_composition
      WHERE transcription_factor_complex_id=?
      " );

    $sth->execute($transcription_factor_complex_id);

    my @transcription_factor_ids;

    while ( my @row = $sth->fetchrow_array ) {
        push @transcription_factor_ids, $row[0];
    }

    my $transcription_factor_adaptor
        = $self->db->get_adaptor('TranscriptionFactor');

    my @components;

    for my $id (@transcription_factor_ids) {
        push @components, $transcription_factor_adaptor->fetch_by_dbID($id);
    }

    return \@components;
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates objects from an executed DBI statement handle.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
    my ( $self, $sth ) = @_;

    my ( @result, $transcription_factor_complex_id,
        $production_name, $display_name );

    $sth->bind_columns( \$transcription_factor_complex_id,
        \$production_name, \$display_name );

    while ( $sth->fetch() ) {
        my $components
            = $self->_fetch_components($transcription_factor_complex_id);

        my $transcription_factor_complex
            = Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex->new(
            -dbID            => $transcription_factor_complex_id,
            -PRODUCTION_NAME => $production_name,
            -DISPLAY_NAME    => $display_name,
            -COMPONENTS      => $components,
            -ADAPTOR         => $self,
            );

        push @result, $transcription_factor_complex;
    }

    return \@result;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex objects
  Example    : $tf_adaptor->store($tf1, $tf2, $tmf3);
  Description: Stores given TranscriptionFactorComplex objects in the database.
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
      INSERT INTO transcription_factor_complex
      (production_name, display_name)
      VALUES (?,?)" );

    my $comp_sth = $self->prepare(
        "INSERT INTO transcription_factor_complex_composition
        (transcription_factor_complex_id,transcription_factor_id)
        VALUES(?,?)"
    );

    my $stored_TranscriptionFactorComplex;

    foreach my $transcription_factor_complex (@args) {
        assert_ref( $transcription_factor_complex,
            'Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex',
            'TranscriptionFactorComplex' );

        if (!(     $transcription_factor_complex->dbID()
                && $transcription_factor_complex->adaptor() == $self
            )
            )
        {
            $stored_TranscriptionFactorComplex
                = $self->fetch_by_production_name(
                $transcription_factor_complex->production_name );

            if ( !$stored_TranscriptionFactorComplex ) {

                $sth->bind_param( 1,
                    $transcription_factor_complex->production_name,
                    SQL_VARCHAR );
                $sth->bind_param( 2,
                    $transcription_factor_complex->display_name,
                    SQL_VARCHAR );

                $sth->execute();

                $transcription_factor_complex->dbID( $self->last_insert_id );
                $transcription_factor_complex->adaptor($self);

                for my $component (
                    @{ $transcription_factor_complex->components() } )
                {
                    $comp_sth->bind_param( 1,
                        $transcription_factor_complex->dbID, SQL_INTEGER );
                    $comp_sth->bind_param( 2, $component->dbID, SQL_INTEGER );
                    $comp_sth->execute();
                }
            }
            else {
                $transcription_factor_complex
                    = $stored_TranscriptionFactorComplex;
                warn( "Using previously stored TranscriptionFactorComplex:\t"
                        . $transcription_factor_complex->production_name()
                        . "\n" );
            }
        }
    }

    return \@args;
}

1;
