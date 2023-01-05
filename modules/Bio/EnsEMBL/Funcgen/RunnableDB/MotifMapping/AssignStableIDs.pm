=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2023] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::AssignStableIDs;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::Utils::SqlHelper;

sub run {
    my $self = shift;

    my $species = $self->param_required('species');
    my $stable_id_mapping_dir = $self->param_required('stable_id_mapping_dir');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $previous_funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor( $species . '_previous_version',
                                               'funcgen' );

    my $max_stableID = get_single_sql_value_from_db(
    $previous_funcgen_adaptor->dbc,
    'SELECT MAX(stable_id) FROM motif_feature');

    my ( $stableID_prefix, $stableID_counter ) = $max_stableID =~ /(\D+)(\d+)/;

    my $max_motif_feature_dbID = get_single_sql_value_from_db(
        $funcgen_adaptor->dbc,
        'SELECT MAX(motif_feature_id) FROM motif_feature'
    );

    my $intersection = $stable_id_mapping_dir . '/intersection.sorted.out';
    my $outfile = $stable_id_mapping_dir . '/new_stableIDs.out';

    open my $mappings_fh, '<', $intersection;
    open my $outfile_fh,  '>', $outfile;

    for (
        my $motif_feature_dbID = 1 ;
        $motif_feature_dbID <= $max_motif_feature_dbID ;
        $motif_feature_dbID++
    )
    {
        my $lines = get_next_chunk($motif_feature_dbID, $mappings_fh);

        my $mapping_exists = 0;
        my ( $stableID, $current_binding_matrix_id, $previous_binding_matrix_id,
            $existing_stableID );

        for my $line ( @{$lines} ) {
            chomp $line;
            my @fields = split /\t/, $line;

            $current_binding_matrix_id  = $fields[6];
            $previous_binding_matrix_id = $fields[13];
            $existing_stableID          = $fields[14];

            if (    $previous_binding_matrix_id ne '.'
                and $current_binding_matrix_id == $previous_binding_matrix_id )
            {
                $mapping_exists = 1;
                last;
            }
        }

        if ($mapping_exists) {
            $stableID = $existing_stableID;
        }
        else {
            $stableID_counter++;
            $stableID = $stableID_prefix . $stableID_counter;
        }
        $outfile_fh->print( $motif_feature_dbID . "\t" . $stableID . "\n" );
    }

    close $mappings_fh;
    close $outfile_fh;

}


sub get_next_chunk {
    my ($motif_feature_id, $mappings_fh) = @_;
    my @lines;

    my $pos;

    while ( readline $mappings_fh ) {
        my $line = $_;
        chomp $line;
        my @fields = split /\t/, $line;

        if ( $motif_feature_id == $fields[3] ) {
            push @lines, $line;
            $pos = tell $mappings_fh;
        }
        else {
            seek $mappings_fh, $pos, 0;
            last;
        }
    }

    return \@lines;
}


sub get_single_sql_value_from_db {
    my ($dbc, $sql_query) = @_;

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
        -DB_CONNECTION => $dbc
    );

    my $rows = $helper->execute_simple(
        -SQL => $sql_query
    );

    return($rows->[0]);
}


1;
