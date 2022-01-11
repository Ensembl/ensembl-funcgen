=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2022] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::CreateWebBed;

use strict;
use autodie;
use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Hash::Util qw( lock_keys );

sub run {
    my $self = shift;

    my $species = $self->param_required('species');
    my $web_bed_dir = $self->param_required('web_bed_dir');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $core_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor( $species , 'core' );

    my $slice_adaptor = $core_adaptor->get_SliceAdaptor;

    # Include non-reference seq regions, important for mouse.
    my $slices = $slice_adaptor->fetch_all('toplevel', undef, 1);

    my %seq_region_id_to_name_lookup;

    for my $slice (@$slices) {
        $seq_region_id_to_name_lookup{$slice->get_seq_region_id} = $slice->seq_region_name;
    }

    lock_keys(%seq_region_id_to_name_lookup);

    my $sql = "
        SELECT
            mf.seq_region_id,
            mf.seq_region_start - 1 chromStart,
            mf.seq_region_end chromEnd,
            mf.stable_id name,
            ROUND(mf.score, 2) score,
            IF(mf.seq_region_strand=1,'+','-') strand,
                  mf.seq_region_start - 1 thickStart,
            mf.seq_region_end thickEnd,
            '0,0,0' itemRgb,
            bm.stable_id binding_matrix,
            mf.motif_feature_id
        FROM
            motif_feature mf
                        JOIN
                binding_matrix bm USING (binding_matrix_id)
  ";

    # Prevent memory issues from buffering
    $funcgen_adaptor->dbc->db_handle->{mysql_use_result} = 1;

    my $helper = $funcgen_adaptor->dbc->sql_helper;
    open my $fh, '>', $web_bed_dir . '/web_track.bed' ;

    my $last_line_printed;

    $helper->execute_no_return(
      -SQL      => $sql,
      -CALLBACK => sub {
        my $row  = shift;

        $row->[0] = $seq_region_id_to_name_lookup{$row->[0]};

        my $line = join "\t", @$row;

        if ($line ne $last_line_printed) {
          $fh->print( $line );
          $fh->print( "\n"  );
          $last_line_printed = $line;
        }
      },
    );

}

1;
