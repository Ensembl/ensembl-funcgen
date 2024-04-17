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

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::CreateMotifFeatureRegulatoryFeatureLoadingFiles;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use autodie;

use feature qw/say/;
use List::Util qw(max);


sub run {
    my $self = shift;

    my $species = $self->param('species');

    my $funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

    # create caches for minimizing db connections
    my %chr2seqRegionID; # cache for converting chromosome names to region IDs
    my %seqRegionID2chr; # cache for converting region IDs to chromosome names

    my $slice_adaptor = $core_adaptor->get_adaptor('Slice');
    my $slices = $slice_adaptor->fetch_all('toplevel');
    for my $slice (@{$slices}) {
        $chr2seqRegionID{ $slice->seq_region_name } = $slice->get_seq_region_id;
        $seqRegionID2chr{ $slice->get_seq_region_id } = $slice->seq_region_name;
    }

    my $RF_intersection_dir =
        $self->param_required('RF_intersection_dir');
    my $filtered_bed_dir = $RF_intersection_dir . '/filtered/';

    open my $MF_table_file, '>',
         $RF_intersection_dir . '/motif_feature_table.tsv';
    open my $MFRF_table_file, '>',
         $RF_intersection_dir . '/motif_feature_regulatory_feature_table.tsv';

    opendir(my $bed_dh, $filtered_bed_dir);
    my @files = readdir $bed_dh;

    shift @files;
    shift @files;

    my %bmName2bmIDcache;
    my %MFdbIDcache;

    my $motif_feature_adaptor = $funcgen_adaptor->get_adaptor('MotifFeature');
    my $mf_iterator = $motif_feature_adaptor->fetch_Iterator;

    while (my $mf = $mf_iterator->next) {
        my $unique_constraint =
            $mf->get_BindingMatrix->dbID . '_'
                . $mf->slice->get_seq_region_id . '_'
                . $mf->start . '_'
                . $mf->strand;

        $MFdbIDcache{$unique_constraint} = $mf->dbID;
    }

    my $motif_feature_dbID = max values %MFdbIDcache;
    $motif_feature_dbID++;

    my $binding_matrix_adaptor = $funcgen_adaptor->get_adaptor('BindingMatrix');

    for my $file (@files) {
        say $file;
        my ($binding_matrix_stable_id, $epigenome_name) = split /___/, $file;
        $epigenome_name =~ s/\.filtered//ig;
        $epigenome_name =~ s/RegulatoryFeatures_//ig;
        if (!exists $bmName2bmIDcache{$binding_matrix_stable_id}) {
            my $binding_matrix =
                $binding_matrix_adaptor
                    ->fetch_by_stable_id($binding_matrix_stable_id);
            $bmName2bmIDcache{$binding_matrix_stable_id} =
                $binding_matrix->dbID;
        }

        open my $fh, '<', $filtered_bed_dir . $file;
        while (readline $fh) {
            my $line = $_;
            chomp $line;

            my $overlap = parse_overlap_line($line, \%chr2seqRegionID);
            my $seq_region_name = $overlap->{seq_region_name};

            my $unique_constraint =
                $bmName2bmIDcache{$binding_matrix_stable_id} . '_'
                    . $overlap->{seq_region_id} . '_'
                    . $overlap->{mf_start} . '_'
                    . $overlap->{mf_strand};

            if (!exists $MFdbIDcache{$unique_constraint}) {
                $MFdbIDcache{$unique_constraint} = $motif_feature_dbID;
                $motif_feature_dbID++;
            }

            print_MF_table_line(
                $overlap, $MF_table_file,
                $MFdbIDcache{$unique_constraint},
                $bmName2bmIDcache{$binding_matrix_stable_id}
            );

            print_MFRF_table_line($overlap, $MFRF_table_file,
                                  $MFdbIDcache{$unique_constraint});

        }

    }

}

sub parse_overlap_line {
    my ($line, $chr2seqRegionID)  = @_;

    my @fields = split /\t/, $line;

    return {
        'seq_region_id' => $chr2seqRegionID->{ $fields[0] },
        'mf_start'      => $fields[1],
        'mf_end'        => $fields[2],
        'mf_score'      => $fields[3],
        'mf_strand'     => $fields[4],
        'rf_dbID'       => $fields[8],
        'epigenome_id'  => $fields[9]
    }

}

sub print_MF_table_line {
    my ( $overlap, $fh, $motif_feature_dbID, $binding_matrix_dbID ) = @_;

    print $fh $motif_feature_dbID . "\t"
        . $binding_matrix_dbID . "\t"
        . $overlap->{seq_region_id} . "\t"
        . $overlap->{mf_start} . "\t"
        . $overlap->{mf_end} . "\t"
        . $overlap->{mf_strand} . "\t"
        . $overlap->{mf_score} . "\n";

}

sub print_MFRF_table_line {
    my ( $overlap, $fh, $motif_feature_dbID ) = @_;

    my $string = $motif_feature_dbID . "\t" . $overlap->{rf_dbID};

    if ( $overlap->{epigenome_id} ) {
        $string .= "\t" . $overlap->{epigenome_id} . "\t" . '1';
    }
    else {
        $string .= "\t" . '\N' . "\t" . '0';
    }
    $string .= "\n";
    print $fh $string;
}


1;
