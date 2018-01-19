#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
my $port;
my $dbname;
my $signal_result_set_id;
my $work_dir;

my %config_hash = (
  "argenrich_file"        => \$argenrich_file,
  "result_set_id"         => \$result_set_id,
  'signal_result_set_id'  => \$signal_result_set_id,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'port'            => \$port,
  'host'            => \$host,
  'dbname'          => \$dbname,
  'work_dir'        => \$work_dir,
);

my $result = GetOptions(
  \%config_hash,
  'result_set_id=s',
  'signal_result_set_id=s',
  'argenrich_file=s',
  'dry_run',
  'user=s',
  'pass=s',
  'port=s',
  'host=s',
  'dbname=s',
  'work_dir=s',
);

die unless(-e $argenrich_file);
die unless($signal_result_set_id);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my @tracking_db_connection_details = (
  -user     => $user,
  -pass     => $pass,
  -host     => $host,
  -port     => $port,
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
  
  $f[1] = 'null' if ($f[1] eq 'NA');
  
  $key_value_pairs{$f[0]} = $f[1];
}
use Hash::Util qw( lock_hash );
lock_hash(%key_value_pairs);

my $sql = qq(insert into result_set_qc_chance (
      signal_result_set_id, analysis_id, p, q, divergence, z_score, percent_genome_enriched, input_scaling_factor, differential_percentage_enrichment,
      control_enrichment_stronger_than_chip_at_bin,
      first_nonzero_bin_at,
      pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome, path
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
    $key_value_pairs{'PCR amplification bias in Input, coverage of 1% of genome'},
    '$work_dir'
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
-- Not really that important
-- See slide 38 on 
-- http://www.ebi.ac.uk/seqdb/confluence/download/attachments/18483313/UCL_ChIPseq_Wilder.pptx?version=2&modificationDate=1442910347000&api=v2
-- dashed green line
--
  `p` double default NULL,
-- Not really that important
  `q` double default NULL,
--
-- This is the main statistic.
-- It is a scaled version of differential_percentage_enrichment. The reason 
-- is that the exact location is important and that is not reflected in 
-- differential_percentage_enrichment.
-- 
--
  `divergence` double default NULL,
--
-- Distance from the mean, if the distribution was standardised to a normal distribution
--
  `z_score` double default NULL,
--
-- Distance between dashed green line and 1
--
  `percent_genome_enriched` double default NULL,
--
-- A suggestion on how to scale the control to equal the background noise 
-- in the signal
--
  `input_scaling_factor` double default NULL,
--
-- It is the greates distance between the cumulative coverage lines of the 
-- control and the signal when plotted into a graph.
--
  `differential_percentage_enrichment` double default NULL,
--
-- Usually the two curves would meet at one. If there is an enrichment in 
-- the control, then this reports the bin number where this happens.
-- A bit visible in diagram d, slide 40 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053734/figure/F2/
--
  `control_enrichment_stronger_than_chip_at_bin`double default NULL,
-- 
-- After sorting the bins from the signal, this is the rank of the first non zero bin.
--
  `first_nonzero_bin_at`double default NULL,
--
-- Proportion of control reads in the highest 1 percent of the bins. The expected value would be 0.01, but only
-- greater deviations from that are reported.
--
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome`double default NULL,
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_chance_id`)
);
SQL
;
  $sql_processor->($sql);
}




