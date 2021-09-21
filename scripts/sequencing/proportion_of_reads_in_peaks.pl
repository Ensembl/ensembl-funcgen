#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

=head1 proportion_of_reads_in_peaks.pl

=head1 SYNOPSIS

=head1 DESCRIPTION
=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;
use File::Path qw( make_path );

my $peak_file;
my $temp_dir;
my $peak_calling;
my $registry;
my $bam_file;

my %config_hash = (
  "peak_file"    => \$peak_file,
  "temp_dir"     => \$temp_dir,
  "peak_calling" => \$peak_calling,
  "registry"     => \$registry,
  "bam_file"     => \$bam_file,
);

my $result = GetOptions(
  \%config_hash,
  'peak_file=s',
  'temp_dir=s',
  'bam_file=s',
  'peak_file=s',
  'peak_calling=s',
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

if (! -e $peak_file) {
  $logger->error("The peak file $peak_file does not exist!\n");
}
if (! -e $bam_file) {
  $logger->error("The signal bam file $bam_file does not exist!\n");
}

if ($temp_dir) {
  if (-d $temp_dir) {
    $logger->info("Temporary directory $temp_dir exists.\n");
  } else {
    $logger->info("Temporary directory $temp_dir does not exist, so will create it.\n");
    make_path($temp_dir);
  }
}

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd run_system_cmd );

my $READS_IN_PEAKS_BAM_FILE = "${bam_file}.peaks.bam";

my $current_cmd = "bedtools intersect -abam $bam_file -b $peak_file > $READS_IN_PEAKS_BAM_FILE";
$logger->info("Running $current_cmd\n", 0, 1);
run_system_cmd($current_cmd, 0, 1);

$current_cmd = "samtools view -c $READS_IN_PEAKS_BAM_FILE";
$logger->info("Running $current_cmd\n", 0, 1);
my $num_reads_in_peaks = run_backtick_cmd($current_cmd);

$current_cmd = "samtools view -c $bam_file";
$logger->info("Running $current_cmd\n", 0, 1);
my $num_reads_in_total = run_backtick_cmd($current_cmd);

unlink($READS_IN_PEAKS_BAM_FILE);

my $proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total;

$logger->info("Proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total = $proportion_of_reads_in_peaks\n");


die;

__END__

Bio::EnsEMBL::Registry->load_all($registry);

my $logic_name = 'Proportion of reads in peaks';
my @proportion_of_reads_in_peaks_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'proportion_of_reads_in_peaks.pl',
        -parameters      => '',
        -description     => 'Computation of the proportion of reads in peaks',
        -display_label   => 'Proportion of reads in peaks',
        -displayable     => undef,
);

my $analysis_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Analysis' );
my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

if (! $analysis && ! $dry_run) {
      $logger->info("No analysis with logic name $logic_name found. Creating one.\n");
      $analysis = Bio::EnsEMBL::Analysis->new(@proportion_of_reads_in_peaks_analysis_details);
      $analysis_adaptor->store($analysis);
}

my $analysis_id;
if (! $analysis && $dry_run) {
  $analysis_id = 'No analysis inserted in dry run mode';
} else {
  $analysis_id = $analysis->dbID;
}
$logger->info("The analysis_id is: $analysis_id.\n");





my $sql_processor;
if ($dry_run) {
  $sql_processor = sub {
    my $sql = shift;
    $logger->info("Dry run: " . $sql . "\n");
  };
} else {
  $sql_processor = sub {
    my $sql = shift;
    $logger->info("Running: " . $sql . "\n");
    $dbc->do($sql);
  };
}

create_table({ 
  sql_processor => $sql_processor,
});

create_insert_sql({
  analysis_id         => $analysis_id,
  feature_set_id      => $feature_set_id,
  prop_reads_in_peaks => $proportion_of_reads_in_peaks,
  total_reads         => $num_reads_in_total,
  sql_processor       => $sql_processor,
});

$logger->finish_log;


=head2 create_insert_sql
=cut
sub create_insert_sql {

  my $param = shift;
  
  my $analysis_id         = $param->{analysis_id};
  my $feature_set_id      = $param->{feature_set_id};
  my $prop_reads_in_peaks = $param->{proportion_of_reads_in_peaks};
  my $total_reads         = $param->{total_reads};
  my $sql_processor       = $param->{sql_processor};  

  my $sql = "INSERT INTO feature_set_qc_prop_reads_in_peaks ("
    . "analysis_id, "
    . "feature_set_id, "
    . "prop_reads_in_peaks, "
    . "total_reads"
    . ")  VALUES ("
  . (
    join ', ', (
      $analysis_id,
      $feature_set_id,
      $proportion_of_reads_in_peaks,
      $num_reads_in_total
      )
    )
  . ");";
  $sql_processor->($sql);
}

=head2 create_table
=cut
sub create_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

my $sql = <<SQL
CREATE TABLE if not exists `peak_calling_qc_prop_reads_in_peaks` (
  `peak_calling_qc_prop_reads_in_peaks_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id`        int(10) unsigned,
  `feature_set_id` int(10) unsigned NOT NULL,
  `prop_reads_in_peaks`       double default NULL,
  `total_reads`      int(10) default NULL
  PRIMARY KEY (`peak_calling_qc_prop_reads_in_peaks_id`)
);
SQL
;
  $sql_processor->($sql);
}

