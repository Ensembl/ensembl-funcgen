=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::CreateMotifFeaturePeakLoadingFiles

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2024] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::CreateMotifFeaturePeakLoadingFiles;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use autodie;

sub run {
    my $self = shift;

    my $species = $self->param('species');

    my $funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

    my $Peaks_intersection_dir =
        $self->param_required('Peaks_intersection_dir');
    my $filtered_bed_dir = $Peaks_intersection_dir . '/filtered/';

    open my $MF_table_file, '>',
         $Peaks_intersection_dir . '/motif_feature_table.tsv';
    open my $MFP_table_file, '>',
         $Peaks_intersection_dir . '/motif_feature_peak_table.tsv';

    opendir(my $bed_dh, $filtered_bed_dir);
    my @files = readdir $bed_dh;

    shift @files;
    shift @files;

    my $motif_feature_dbID = 1;

    my $binding_matrix_adaptor = $funcgen_adaptor->get_adaptor('BindingMatrix');
    my $slice_adaptor          = $core_adaptor->get_adaptor('Slice');

    my %bmName2bmIDcache;
    my %seqRegionName2seqRegionIDcache;
    my %MFdbIDcache;

    for my $file (@files) {
        my ($binding_matrix_stable_id) = split /___/, $file;
        $binding_matrix_stable_id =~ s/\.sorted\.bed//ig;
        if (!exists $bmName2bmIDcache{$binding_matrix_stable_id}) {
            my $binding_matrix =
                $binding_matrix_adaptor->fetch_by_stable_id($binding_matrix_stable_id);
            $bmName2bmIDcache{$binding_matrix_stable_id} = $binding_matrix->dbID;
        }

        open my $fh, '<', $Peaks_intersection_dir . '/filtered/' . $file;
        while (readline $fh) {
            my $line = $_;
            chomp $line;

            my $overlap         = parse_overlap_line($line);
            my $seq_region_name = $overlap->{seq_region_name};

            if (!exists $seqRegionName2seqRegionIDcache{$seq_region_name}) {
                my $slice =
                    $slice_adaptor->fetch_by_name($seq_region_name);
                $seqRegionName2seqRegionIDcache{$seq_region_name} =
                    $slice->get_seq_region_id;
            }

            my $unique_constraint =
                $bmName2bmIDcache{$binding_matrix_stable_id}
                    . $seqRegionName2seqRegionIDcache{$seq_region_name}
                    . $overlap->{mf_start}
                    . $overlap->{mf_strand};

            if (!exists $MFdbIDcache{$unique_constraint}) {
                $MFdbIDcache{$unique_constraint} = $motif_feature_dbID;
                $motif_feature_dbID++;
            }

            print_MF_table_line(
                $overlap,
                $MF_table_file,
                $MFdbIDcache{$unique_constraint},
                $bmName2bmIDcache{$binding_matrix_stable_id},
                $seqRegionName2seqRegionIDcache{$seq_region_name}
            );

            print_MFP_table_line($overlap, $MFP_table_file,
                                 $MFdbIDcache{$unique_constraint});

        }
    }
}

sub parse_overlap_line {
    my $line = shift;

    my @fields = split /\t/, $line;

    my ($motif_score, $seq_region_name) = split "_", $fields[3];

    return {
        'seq_region_name' => $seq_region_name,
        'mf_start'        => $fields[1],
        'mf_end'          => $fields[2],
        'mf_score'        => $motif_score,
        'mf_strand'       => $fields[5],
        'peak_dbID'       => $fields[9],
    }

}

sub print_MF_table_line {
    my ($overlap, $fh, $motif_feature_dbID, $binding_matrix_dbID,
        $seq_region_id)
        = @_;

    my $strand = $overlap->{mf_strand};
    if ($strand eq '+') {
        $strand = 1;
    }
    elsif ($strand eq '-') {
        $strand = -1;
    }

    print $fh $motif_feature_dbID . "\t"
        . $binding_matrix_dbID . "\t"
        . $seq_region_id . "\t"
        . $overlap->{mf_start} . "\t"
        . $overlap->{mf_end} . "\t"
        . $strand . "\t"
        . $overlap->{mf_score} . "\n";

}

sub print_MFP_table_line {
    my ($overlap, $fh, $motif_feature_dbID) = @_;

    print $fh $motif_feature_dbID . "\t" . $overlap->{peak_dbID} . "\n";
}

1;
