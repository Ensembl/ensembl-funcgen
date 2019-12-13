#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

perl scripts/sequencing/remove_intermediary_bam_files.pl \
  --experiment GM19193_WCE_ChIP-Seq_control_TF_no1_ENCODE86 \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens



experiments=`r2 -N mnuhn_testdb2_homo_sapiens_funcgen_91_38 -e "select name from experiment"`

for experiment in $experiments
do
  echo Processing $experiment
  perl scripts/sequencing/remove_intermediary_bam_files.pl \
    --experiment $experiment \
    --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --species homo_sapiens \
    --data_root_dir /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/dbfiles
done

perl scripts/sequencing/remove_intermediary_bam_files.pl \
  --experiment IMR90_H3K27ac_ChIP-Seq_RoadmapEpigenomics85 \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --data_root_dir /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/dbfiles

perl scripts/sequencing/remove_intermediary_bam_files.pl \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species mus_musculus \
  --data_root_dir /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/dbfiles



=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $experiment_name;
my $registry;
my $species;
my $data_root_dir;

my %config_hash = (
  'experiment'    => \$experiment_name,
  'registry'      => \$registry,
  'species'       => \$species,
  'data_root_dir' => \$data_root_dir,
);

my $result = GetOptions(
  \%config_hash,
  'experiment=s',
  'registry=s',
  'species=s',
  'data_root_dir=s',
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Experiment' );
my $alignment_adaptor  = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Alignment' );

my @experiment_list;

if (defined $experiment_name) {
  my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
  push @experiment_list, $experiment;
} else {
  @experiment_list = @{$experiment_adaptor->generic_fetch};
  #die;
}

# print Dumper(map { $_->name } @experiment_list);
# die;
foreach my $experiment (@experiment_list) {
  remove_intermediary_bam_files_by_Experiment($experiment);
}

$logger->finish_log;

sub remove_intermediary_bam_files_by_Experiment {

  my $experiment = shift;

  my $alignments_with_duplicates = $alignment_adaptor->fetch_all_with_duplicates_by_Experiment($experiment);

  $logger->info("Found " . scalar @$alignments_with_duplicates . " alignments with duplicates.\n");

  ALIGNMENT:
  foreach my $alignment_with_duplicates (@$alignments_with_duplicates) {

    if (! $alignment_with_duplicates->has_bam_DataFile) {
      $logger->info($alignment_with_duplicates->name . " has no bam file.\n");
      next ALIGNMENT;
    }
    my $bam_file = $alignment_with_duplicates->fetch_bam_DataFile;
    
    my $full_file_name = join '/',
      $data_root_dir,
      $bam_file->relative_ftp_site_path
    ;
    if (! -e $full_file_name) {
      $logger->error("Couldn't find alignment $full_file_name!\n");
      $logger->finish_log;
      die;
    }
    $logger->info("Deleting $full_file_name\n");
    unlink($full_file_name);
    $alignment_with_duplicates->_delete_bam_file_from_db;
    $bam_file->_delete_from_db;
  }
}


