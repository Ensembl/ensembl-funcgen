#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 load_argenrich_file.pl

=head1 SYNOPSIS

# Wiggletools needs a sorted bed file:
sort -k1,1 /lustre/scratch110/ensembl/funcgen/il4/ersaCTTV/CTTV020_epigenomes_of_cell_lines_PILOT/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt > homo_sapiens_.CCAT_chr_lengths.sorted.txt

/software/ensembl/funcgen/argenrichformregions.pl homo_sapiens_.CCAT_chr_lengths.sorted.txt

TEMPDIR=/lustre/scratch109/ensembl/funcgen/mn1/temp/argenrich
mkdir -p $TEMPDIR

cp /warehouse/ensembl10/funcgen/alignments/homo_sapiens/GRCh38/3526/UT7:hist:BR2_H3K4me3_3526_bwa_samse_1.bam \
  /warehouse/ensembl10/funcgen/alignments/homo_sapiens/GRCh38/3526/UT7:hist:BR2_WCE_3526_bwa_samse_1.bam \
  $TEMPDIR

samtools index $TEMPDIR/UT7:hist:BR2_H3K4me3_3526_bwa_samse_1.bam
samtools index $TEMPDIR/UT7:hist:BR2_WCE_3526_bwa_samse_1.bam

# Number of reads in signal
ipsz_parameter=$(samtools view -c $TEMPDIR/UT7:hist:BR2_H3K4me3_3526_bwa_samse_1.bam)

# Number of reads in control
inputsz_parameter=$(samtools view -c $TEMPDIR/UT7:hist:BR2_WCE_3526_bwa_samse_1.bam)

/nfs/users/nfs_m/mn1/work_dir_faang/argenrich.R --args plot=TRUE outdir=$TEMPDIR \
  ipsz=$ipsz_parameter inputsz=$inputsz_parameter \
  ip=$TEMPDIR/UT7:hist:BR2_H3K4me3_3526_bwa_samse_1.bam \
  input=$TEMPDIR/UT7:hist:BR2_WCE_3526_bwa_samse_1.bam \
  outfile=argenrich.stdout.txt

./scripts/sequencing/load_argenrich_qc_file.pl \
  --argenrich_file /lustre/scratch109/ensembl/funcgen/mn1/ersa/debug/F36P:hist:BR2_H3K27me3_3526/argenrich_outfile.txt \
  --control_result_set_id 1 \
  --signal_result_set_id 2 --user ensadmin --pass xxx --host ens-genomics1 --dbname mn1_faang_tracking_homo_sapiens_funcgen_81_38

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $argenrich_file;
my $result_set_id;
my $dry_run;
my $user;
my $pass;
my $host;
my $dbname;
my $control_result_set_id;
my $signal_result_set_id;

my %config_hash = (
  "argenrich_file"        => \$argenrich_file,
  "result_set_id"         => \$result_set_id,
  'control_result_set_id' => \$control_result_set_id,
  'signal_result_set_id'  => \$signal_result_set_id,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'host'            => \$host,
  'dbname'          => \$dbname,
);

my $result = GetOptions(
  \%config_hash,
  'result_set_id=s',
  'control_result_set_id=s',
  'signal_result_set_id=s',
  'argenrich_file=s',
  'dry_run',
  'user=s',
  'pass=s',
  'host=s',
  'dbname=s',
);

die unless(-e $argenrich_file);
# die unless($control_result_set_id);
die unless($signal_result_set_id);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my @tracking_db_connection_details = (
  -user     => $user,
  -pass     => $pass,
  -host     => $host,
  -dbname   => $dbname,
);
my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(@tracking_db_connection_details);
my $analysis_id;

if ($dry_run) {
  $analysis_id = 'Dryrun';
} else {
  $analysis_id = create_analysis_if_not_exists($dbc);
}

open IN, $argenrich_file;

my %other_column_name_for = (
   'Control enrichment stronger than ChIP at bin' => '',
   'Zero-enriched IP, maximum difference at bin' => '',
   'PCR amplification bias in Input, coverage of 1% of genome' => '',
);

my %key_value_pairs;
LINE: while (my $current_line = <IN>) {
  chomp $current_line;
  #print " - $current_line\n";
  my @f = split '=', $current_line;
  next LINE unless(@f == 2);
  
  $key_value_pairs{$f[0]} = $f[1];
}
use Hash::Util qw( lock_hash );
lock_hash(%key_value_pairs);

my $sql = qq(insert into result_set_qc_chance (
      signal_result_set_id, analysis_id, p, q, divergence, z_score, percent_genome_enriched, input_scaling_factor, differential_percentage_enrichment,
      control_enrichment_stronger_than_chip_at_bin,
      zero_enriched_ip_maximum_difference_at_bin,
      pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome
    ) values (
    $signal_result_set_id,
    $analysis_id,
    $key_value_pairs{'p'}, 
    $key_value_pairs{'q'}, 
    $key_value_pairs{'divergence'}, 
    $key_value_pairs{'z_score'}, 
    $key_value_pairs{'percent_genome_enriched'}, 
    $key_value_pairs{'input_scaling_factor'}, 
    $key_value_pairs{'differential_percentage_enrichment'},
    $key_value_pairs{'Control enrichment stronger than ChIP at bin'},
    $key_value_pairs{'Zero-enriched IP, maximum difference at bin'},
    $key_value_pairs{'PCR amplification bias in Input, coverage of 1% of genome'}
  )
);

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

#print Dumper(\%key_value_pairs);
create_table({
  sql_processor => $sql_processor
});
$sql_processor->($sql);

$logger->finish_log;
exit;

sub create_analysis_if_not_exists {

  my $dbc = shift;
  my $logic_name = 'QC Chance';
  my @phantom_peak_analysis_details = (
	  -logic_name      => $logic_name,
	  -program         => 'argenrich.R',
	  -parameters      => '',
	  -description     => '',
	  -display_label   => 'Chance',
	  -displayable     => undef,
  );
  
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -DBCONN => $dbc,  
  );
  my $analysis_adaptor = $dba->get_AnalysisAdaptor();
  my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  if (! $analysis && ! $dry_run) {
	$logger->info("No analysis with logic name $logic_name found. Creating one.");
	$analysis = Bio::EnsEMBL::Analysis->new(@phantom_peak_analysis_details);
	$analysis_adaptor->store($analysis);
  }
  my $analysis_id = $analysis->dbID;
  return $analysis_id;
}

=head2 create_table
=cut
sub create_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

my $sql = <<SQL
 CREATE TABLE if not exists `result_set_qc_chance` (
  `result_set_qc_chance_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `signal_result_set_id` int(10),
  `analysis_id`        int(10) unsigned,
  `p` double default NULL,
  `q` double default NULL,
  `divergence` double default NULL,
  `z_score` double default NULL,
  `percent_genome_enriched` double default NULL,
  `input_scaling_factor` double default NULL,
  `differential_percentage_enrichment` double default NULL,
  `control_enrichment_stronger_than_chip_at_bin`double default NULL,
  `zero_enriched_ip_maximum_difference_at_bin`double default NULL,
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome`double default NULL,
  `path` varchar(100) NOT NULL,
  PRIMARY KEY (`result_set_qc_chance_id`)
);
SQL
;
  $sql_processor->($sql);
}




