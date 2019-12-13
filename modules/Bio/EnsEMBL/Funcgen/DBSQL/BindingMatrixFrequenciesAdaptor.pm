#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor
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

package Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;    #sql_types barewords import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_all_by_BindingMatrix

  Arg [1]    : Bio::EnsEMBL::Funcgen::BindingMatrix
  Example    : my @bmf = @{$bmf_adaptor->fetch_all_by_BindingMatrix($bm)};
  Description: Fetches BindingMatrixFrequencies object given it's BindingMatrix
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies objects
  Exceptions : Throws if BindingMatrix is not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_BindingMatrix{
  my ($self, $binding_matrix) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::BindingMatrix', $binding_matrix);

  my $constraint = " bmf.binding_matrix_id = ?";

  $self->bind_param_generic_fetch($binding_matrix->dbID,    SQL_INTEGER);

  return $self->generic_fetch($constraint);
}

=head2 fetch_by_BindingMatrix_position_nucleotide

  Arg [1]    : Bio::EnsEMBL::Funcgen::BindingMatrix
  Arg [2]    : integer, position of the frequencies matrix
  Arg [3]    : string, nucleotide of the frequencies matrix
  Example    : my @bmf = @{$bmf_adaptor->fetch_by_BindingMatrix($bm, $position, $nucleotide)};
  Description: Fetches BindingMatrixFrequencies object given it's BindingMatrix, position and nucleotide
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies
  Exceptions : Throws if BindingMatrix is not valid
               Throws if position is not specified
               Throws if nucleotide is not specified
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_BindingMatrix_position_nucleotide{
  my ($self, $binding_matrix, $position, $nucleotide) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::BindingMatrix', $binding_matrix);
  throw('Must specify a position') if ! defined $position;
  throw('Must specify a nucleotide') if ! defined $nucleotide;

  my $constraint = " bmf.binding_matrix_id = ?";
  $constraint .= " AND bmf.position = ?";
  $constraint .= " AND bmf.nucleotide = ?";

  $self->bind_param_generic_fetch($binding_matrix->dbID,    SQL_INTEGER);
  $self->bind_param_generic_fetch($position,           SQL_INTEGER);
  $self->bind_param_generic_fetch($nucleotide,           SQL_VARCHAR);

  return $self->generic_fetch($constraint)->[0];
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
    return ( [ 'binding_matrix_frequencies', 'bmf' ] );
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
        qw( bmf.binding_matrix_frequencies_id bmf.binding_matrix_id bmf.position bmf.nucleotide bmf.frequency
    );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates objects from an executed DBI statement handle.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
    my ( $self, $sth ) = @_;

    my ( @result, $bm_frequencies_id, $binding_matrix_id, $position,
        $nucleotide, $frequency );

    $sth->bind_columns( \$bm_frequencies_id, \$binding_matrix_id, \$position,
        \$nucleotide, \$frequency );

    my $binding_matrix_adaptor = $self->db->get_adaptor('BindingMatrix');
    my %binding_matrix_cache;

    while ( $sth->fetch() ) {

        if ( !exists $binding_matrix_cache{$binding_matrix_id} ) {
            $binding_matrix_cache{$binding_matrix_id}
                = $binding_matrix_adaptor->fetch_by_dbID($binding_matrix_id);
        }

        my $bm_frequencies
            = Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
            -dbID           => $bm_frequencies_id,
            -BINDING_MATRIX => $binding_matrix_cache{$binding_matrix_id},
            -POSITION       => $position,
            -NUCLEOTIDE     => $nucleotide,
            -FREQUENCY      => $frequency,
            -ADAPTOR        => $self,
            );

        push @result, $bm_frequencies;

    }

    return \@result;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies objects
  Example    : $bm_frequencies_adaptor->store($bmf1, $bmf2, $bmf3);
  Description: Stores given BindingMatrixFrequencies objects in the database.
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
      INSERT INTO binding_matrix_frequencies
      (binding_matrix_id, position, nucleotide, frequency)
      VALUES (?, ?, ?, ?)" );

    my $stored_bm_frequencies;

    foreach my $bm_frequencies (@args) {
        assert_ref( $bm_frequencies,
            'Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies',
            'BindingMatrixFrequencies' );
        $self->db->is_stored_and_valid(
            'Bio::EnsEMBL::Funcgen::BindingMatrix',
            $bm_frequencies->get_BindingMatrix
        );

        if (!( $bm_frequencies->dbID() && $bm_frequencies->adaptor() == $self )
            )
        {
            #Check for previously stored BindingMatrixFrequencies
            $stored_bm_frequencies
                = $self->fetch_by_BindingMatrix_position_nucleotide(
                $bm_frequencies->get_BindingMatrix(),
                $bm_frequencies->position(),
                $bm_frequencies->nucleotide()
                );

            if ( !$stored_bm_frequencies ) {

                $sth->bind_param( 1,
                    $bm_frequencies->get_BindingMatrix()->dbID(), SQL_INTEGER );
                $sth->bind_param( 2, $bm_frequencies->position(),
                    SQL_INTEGER );
                $sth->bind_param( 3, $bm_frequencies->nucleotide(),
                    SQL_VARCHAR );
                $sth->bind_param( 4, $bm_frequencies->frequency(),
                    SQL_INTEGER );

                $sth->execute();
                $bm_frequencies->dbID( $self->last_insert_id );
                $bm_frequencies->adaptor($self);

            }
            else {
                $bm_frequencies = $stored_bm_frequencies;
                warn( "Using previously stored BindingMatrixFrequencies:\t"
                        . $bm_frequencies->binding_matrix_frequencies()
                        ->dbID()
                        . "\n" );

            }
        }
    }

    return \@args;
}


1;
