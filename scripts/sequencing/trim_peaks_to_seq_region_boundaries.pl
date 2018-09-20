#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 405ress or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 trim_peaks_to_seq_region_boundaries.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

perl scripts/sequencing/trim_peaks_to_seq_region_boundaries.pl \
  --peak_calling H1ESC_H4K20me1_ccat_histone_ENCODE \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --dry_run 1 \
  --species homo_sapiens

perl scripts/sequencing/trim_peaks_to_seq_region_boundaries.pl \
  --peak_calling IMR90_H3K79me2_ccat_histone_Roadmap_Epigenomics \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --dry_run 1 \
  --species homo_sapiens

perl scripts/sequencing/trim_peaks_to_seq_region_boundaries.pl \
  --peak_calling IMR90_H4K8ac_ccat_histone_Roadmap_Epigenomics \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --dry_run 1 \
  --species homo_sapiens

peak_callings=`r2 -N mnuhn_testdb2_homo_sapiens_funcgen_91_38 -e "select name from peak_calling"`

for peak_calling in $peak_callings
do
  echo Processing $peak_calling
  perl scripts/sequencing/trim_peaks_to_seq_region_boundaries.pl \
    --peak_calling $peak_calling \
    --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --dry_run 0 \
    --species homo_sapiens
done

perl scripts/sequencing/trim_peaks_to_seq_region_boundaries.pl \
  --peak_calling Placenta_H3K9me3_ccat_histone_Roadmap_Epigenomics \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --dry_run 1 \
  --species homo_sapiens

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $peak_calling_name;
my $registry;
my $dry_run;
my $species;
my $tempdir;

my %config_hash = (
  'peak_calling'   => \$peak_calling_name,
  'registry'       => \$registry,
  'dry_run'        => \$dry_run,
  'species'        => \$species,
  'tempdir'        => \$tempdir,
);

my $result = GetOptions(
  \%config_hash,
  'peak_calling=s',
  'registry=s',
  'dry_run=s',
  'species=s',
  'tempdir=s',
);

$Data::Dumper::Maxdepth = 1;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'PeakCalling' );

my $peak_calling = $peak_calling_adaptor->fetch_by_name($peak_calling_name);

if (! defined $peak_calling) {
  die("Couldn't find peak calling with name: $peak_calling_name");
}

my $peak_adaptor  = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Peak' );
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'Slice' );

my $peaks_file = "$tempdir/peaks.${peak_calling_name}.bed";

use File::Path qw(make_path remove_tree);
make_path($tempdir);

open my $peaks_fh, '>', $peaks_file or die("Can't open file $peaks_file!");

$peak_adaptor->_bulk_export_to_bed_by_PeakCalling($peak_calling, $peaks_fh);
$peaks_fh->close;

open my $peaks_fh, '<', $peaks_file;

my %seq_region_name_to_slice;

while (my $current_line = <$peaks_fh>) {

  chomp $current_line;
  my $hash = $peak_adaptor->_parse_bed_line($current_line);
  
  my $seq_region_name  = $hash->{seq_region_name};
  my $seq_region_start = $hash->{seq_region_start};
  my $seq_region_end   = $hash->{seq_region_end};
  my $peak_id          = $hash->{peak_id};
  
  if (! exists $seq_region_name_to_slice{$seq_region_name}) {
    my $slice = $slice_adaptor->fetch_by_region(undef, $seq_region_name);
    $seq_region_name_to_slice{$seq_region_name} = $slice;
  }
  
  my $current_slice = $seq_region_name_to_slice{$seq_region_name};
  
  my $needs_fixing = undef;
  
  if ($seq_region_start < $current_slice->seq_region_start) {
    print "not ok start, start = $seq_region_start for peak with id $peak_id\n";
    $needs_fixing = 1;
  }
  if ($seq_region_end > $current_slice->seq_region_end) {
    print "not ok end, end of peak is $seq_region_end, end of slice is " . $current_slice->seq_region_end .  "\n";
    $needs_fixing = 1;
  }
  
  if ($needs_fixing) {
    fix_peak_boundaries($peak_id, $current_slice);
  }
}

$peaks_fh->close;

sub fix_peak_boundaries {

  my $peak_id = shift;
  my $slice   = shift;
  
  my $peak = $peak_adaptor->fetch_by_dbID($peak_id);
  
  if ($peak->seq_region_start < $slice->start ) {
    $peak->seq_region_start($slice->start);
  }
  if ($peak->seq_region_end > $slice->end ) {
    $peak->seq_region_end($slice->end);
  }
  
  if (! $dry_run) {
    print "Updating peak\n";
    eval {
        $peak_adaptor->update($peak);
    };
    if ($@) {
        my $is_duplicate_entry = $@ =~ /DBD::mysql::st execute failed: Duplicate entry/;
        if (! $is_duplicate_entry) {
            confess($@);
        }
        $logger->info("Peak already exists with the trimmed values, so deleting this one.");
        $peak_adaptor->_delete($peak);
    }
  }
}

$logger->finish_log;
exit(0);
