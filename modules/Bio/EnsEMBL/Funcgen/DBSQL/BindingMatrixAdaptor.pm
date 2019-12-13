#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
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

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor - A database adaptor for fetching and
storing Funcgen BindingMatrix objects.

=head1 SYNOPSIS

my $matrix_adaptor = $db->get_BindingMatrixAdaptor();
my @matrices = @{$matrix_adaptor->fetch_all_by_name("MA0122.1")};

=head1 DESCRIPTION

The BindingMatrixAdaptor is a database adaptor for storing and retrieving
Funcgen BindingMatrix objects.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::BindingMatrix

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );

use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrix::Constants qw ( :all );
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types barewords import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string - name of Matrix
  Example    : my $matrix = $matrix_adaptor->fetch_by_name('MA0122.1');
  Description: Fetches matrix objects given a name
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix object
  Exceptions : Throws if no name if defined
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name {
    my ( $self, $name ) = @_;
    throw('Must specify a BindingMatrix name') if !defined $name;

    my $constraint = ' bm.name = ? ';
    $self->bind_param_generic_fetch( $name, SQL_VARCHAR );

    my $result = $self->generic_fetch($constraint);

    return $result->[0];
}

=head2 fetch_by_stable_id

  Arg [1]    : string - stable ID of Matrix
  Example    : my $matrix = $matrix_adaptor->fetch_by_stable_id('ENSPFM001');
  Description: Fetches matrix objects given a stable_id
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix object
  Exceptions : Throws if no stable_id is defined
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_stable_id {
    my ( $self, $stable_id ) = @_;
    throw('Must specify a BindingMatrix stable_id') if !defined $stable_id;

    my $constraint = ' bm.stable_id = ? ';
    $self->bind_param_generic_fetch( $stable_id, SQL_VARCHAR );

    my $result = $self->generic_fetch($constraint);

    return $result->[0];
}

=head2 fetch_all_by_TranscriptionFactor

  Arg [1]    : Bio::EnsEMBL::Funcgen::TranscriptionFactor object
  Example    : my $matrices =
                 $matrix_adaptor->fetch_all_by_TranscriptionFactor($tf)
  Description: Fetches all Binding Matrices that are linked to the given
                 Transcription Factor
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Exceptions : Throws if no transcription_factor paramater is defined
               Throws if parameter passed is not a
                 Bio::EnsEMBL::Funcgen::TranscriptionFactor object
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_TranscriptionFactor {
    my ($self, $transcription_factor) = @_;

    throw('Must pass a TranscriptionFactor parameter') if !defined $transcription_factor;
    assert_ref($transcription_factor,
               'Bio::EnsEMBL::Funcgen::TranscriptionFactor',
               'TranscriptionFactor');

    my $constraint = ' tf.transcription_factor_id = ? ';
    $self->bind_param_generic_fetch( $transcription_factor->dbID(), SQL_INTEGER );

    my $result = $self->generic_fetch($constraint);

    return $result;
}

sub _left_join {
    my $self = shift;

    return (['binding_matrix_transcription_factor_complex', "bmtfc.binding_matrix_id = bm.binding_matrix_id"],
            ['transcription_factor_complex_composition', "tfcc.transcription_factor_complex_id=bmtfc.transcription_factor_complex_id"],
            ['transcription_factor', "tf.transcription_factor_id=tfcc.transcription_factor_id"]
    );
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
    return([ 'binding_matrix', 'bm' ],
           [ 'binding_matrix_transcription_factor_complex', 'bmtfc' ],
           [ 'transcription_factor_complex_composition', 'tfcc'],
           [ 'transcription_factor', 'tf'],
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
  my $self = shift;

    return qw(
      bm.binding_matrix_id
      bm.name
      bm.threshold
      bm.source
      bm.stable_id
    );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates objects from an executed DBI statement handle.
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;

	my (@result, $matrix_id, $name, $thresh, $source, $stable_id);
	$sth->bind_columns(\$matrix_id, \$name, \$thresh, \$source, \$stable_id);

  while ( $sth->fetch() ) {

    my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new
    (
     -dbID         => $matrix_id,
     -NAME         => $name,
     -THRESHOLD    => $thresh,
     -SOURCE       => $source,
     -STABLE_ID    => $stable_id,
     -ADAPTOR      => $self,
    );

	  push @result, $matrix;

	}

	return \@result;
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::BindingMatrix objects
  Example    : $matrix_adaptor->store($m1, $m2, $m3);
  Description: Stores given Matrix objects in the database.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @args = @_;

    my $sth
        = $self->prepare(
        "INSERT INTO binding_matrix (name, threshold, source, stable_id) VALUES (?, ?, ?, ?)"
        );


    my $stored_matrix;

    foreach my $matrix (@args) {
        assert_ref( $matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
            'BindingMatrix' );

        if ( !( $matrix->dbID() && $matrix->adaptor() == $self ) ) {

            #Check for previously stored BindingMatrix
            $stored_matrix = $self->fetch_by_name( $matrix->name() );

            if ( !$stored_matrix ) {

                $sth->bind_param( 1, $matrix->name(),      SQL_VARCHAR );
                $sth->bind_param( 2, $matrix->threshold(), SQL_DOUBLE );
                $sth->bind_param( 3, $matrix->source(),    SQL_VARCHAR );
                $sth->bind_param( 4, $matrix->stable_id(), SQL_VARCHAR );

                $sth->execute();

                $matrix->dbID( $self->last_insert_id );
                $matrix->adaptor($self);

                if ($matrix->{elements}){
                  $self->_store_frequencies($matrix);
                }

                if ($matrix->{associated_transcription_factor_complexes}){
                  $self->_store_binding_matrix_transcription_factor_complex($matrix);
                }

            }
            else {
                $matrix = $stored_matrix;
                warn(     "Using previously stored Matrix:\t"
                        . $matrix->name()
                        . "\n" );
            }
        }
    }

    return \@args;
}

sub _store_frequencies {
    my ( $self, $binding_matrix ) = @_;

    if($binding_matrix->unit ne FREQUENCIES){
        throw("Can not store a matrix with units other than FREQUENCIES");
    }

    my $binding_matrix_frequencies_adaptor
        = $self->db->get_adaptor('BindingMatrixFrequencies');

    for my $frequency ( @{ $binding_matrix->{elements} } ) {
        $frequency->{binding_matrix} = $binding_matrix;
        $binding_matrix_frequencies_adaptor->store($frequency);
    }
}

sub _store_binding_matrix_transcription_factor_complex {
    my ( $self, $matrix ) = @_;

    assert_ref( $matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    my $select_sth = $self->prepare(
        "SELECT * FROM binding_matrix_transcription_factor_complex
         WHERE binding_matrix_id=? AND transcription_factor_complex_id=?"
    );

    my $store_sth = $self->prepare(
        "INSERT INTO binding_matrix_transcription_factor_complex
        (binding_matrix_id, transcription_factor_complex_id) VALUES(?, ?)"
    );

    for
      my $complex ( @{ $matrix->{associated_transcription_factor_complexes} } )
    {

        assert_ref( $complex,
            'Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex',
            'TranscriptionFactorComplex' );

        $select_sth->execute( $matrix->dbID, $complex->dbID );

        # check if the binding_matrix has already been associated with
        # the transcription_factor_complex
        if ( !$select_sth->fetchrow_array ) {
            $store_sth->bind_param( 1, $matrix->dbID,  SQL_INTEGER );
            $store_sth->bind_param( 2, $complex->dbID, SQL_INTEGER );

            $store_sth->execute();
        }
    }
}

1;
