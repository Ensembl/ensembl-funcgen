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

=head1 load_phantom_peak_file.pl

=head1 SYNOPSIS

export R_LIBS=/software/ensembl/funcgen/R-modules

/software/R-3.2.2/bin/Rscript /software/ensembl/funcgen/spp_package/run_spp.R \
  -c=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/K562:hist:BR1_H3K4me3_3526_bwa_samse_1.bam \
  -savp -out=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams

./scripts/sequencing/load_phantom_peak_file.pl  \
    --result_set_id 20  \
    --result_file /lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/test \
    --dry_run \
    --user ensro --host ens-genomics2 --dbname mn1_faang_tracking_homo_sapiens_funcgen_81_38 \

=head1 DESCRIPTION

Downloads data for input sets registered in the data tracking database.

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $result_file;
my $result_set_id;
my $dry_run;
my $user;
my $pass;
my $host;
my $port;
my $dbname;
my $work_dir;
my $bam_file;

my %config_hash = (
  "result_file"   => \$result_file,
  "result_set_id" => \$result_set_id,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'port'            => \$port,
  'host'            => \$host,
  'dbname'          => \$dbname,
  'work_dir'        => \$work_dir,
  'bam_file'        => \$bam_file,
);

my $result = GetOptions(
  \%config_hash,
  'result_set_id=s',
  'result_file=s',
  'dry_run',
  'user=s',
  'pass=s',
  'port=s',
  'host=s',
  'dbname=s',
  'work_dir=s',
  'bam_file=s',
);

if (! $result_file) {
  die("The result_file parameter was not specified!");
}
if (! -e $result_file) {
  die("The result_file ($result_file) specified on the command line does not exist!");
}
if (! $result_set_id) {
  die("The result_set_id parameter was not specified!");
}

my @tracking_db_connection_details = (
    -user     => $user,
    -pass     => $pass,
    -port     => $port,
    -host     => $host,
    -dbname   => $dbname,
);

my $logic_name = 'phantom peak quality tools';
my @phantom_peak_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'run_spp.R',
        -parameters      => '',
        -description     => 'Computes enrichment and quality measures and fragment lengths for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data',
        -display_label   => 'phantom peak quality tools',
        -displayable     => undef,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(@tracking_db_connection_details);
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

my $sql_processor;
if ($dry_run) {
  $sql_processor = sub {
    my $sql = shift;
    $logger->info($sql . "\n");
  };
} else {
  $sql_processor = sub {
    my $sql = shift;
    $logger->info($sql . "\n");
    $dbc->do($sql);
  };
}

create_table({
  sql_processor => $sql_processor
});

open IN, $result_file;

while (my $current_line = <IN>) {
  chomp $current_line;
  #print " - $current_line\n";
  my @f = split "\t", $current_line;
#   print Dumper(\@f);
  (
    my $filename,
    my $numReads,
    my $estFragLenTriple,
    my $corr_estFragLenTriple,
    my $phantomPeak,
    my $corr_phantomPeak,
    my $argmin_corr,
    my $min_corr,
    my $NSC,
    my $RSC,
    my $QualityTag,
  ) = @f;
  
  (
    my $estFragLen,
    my $estFragLen2,
    my $estFragLen3,
  ) = split ',', $estFragLenTriple;
  
  if ($estFragLen2 eq '') {
    $estFragLen2 = 'null';
  }
  if ($estFragLen3 eq '') {
    $estFragLen3 = 'null';
  }
  
  (
    my $corr_estFragLen,
    my $corr_estFragLen2,
    my $corr_estFragLen3,
  ) = split ',', $corr_estFragLenTriple;
  
  if ($corr_estFragLen2 eq '') {
    $corr_estFragLen2 = 'null';
  }
  if ($corr_estFragLen3 eq '') {
    $corr_estFragLen3 = 'null';
  }

  sub quote {
    my $string = shift;
    return '"' . $string . '"'
  }
  
  # This can happen, if the data is very bad. In that case the quality tag 
  # will indicate that with
  #
  if ($RSC eq 'Inf') {
    $RSC = 'null';
    die unless ($QualityTag == 2);
  }
  
  my $sql = "INSERT ignore INTO result_set_qc_phantom_peak ("
#  . "result_set_qc_phantom_peak_id, "
  . "result_set_id, "
  . "analysis_id, "
  . "filename, "
  . "numReads, "
  . "estFragLen, "
  . "estFragLen2, "
  . "estFragLen3, "
  . "corr_estFragLen, "
  . "corr_estFragLen2, "
  . "corr_estFragLen3, "
  . "phantomPeak, "
  . "corr_phantomPeak, "
  . "argmin_corr, "
  . "min_corr, "
  . "NSC, "
  . "RSC, "
  . "QualityTag, "
  . "path "
  . ")  VALUES ("
  . (
    join ', ', (
        $result_set_id,
        $analysis_id,
        #quote($filename),
        quote($bam_file),
        $numReads,

        $estFragLen,
        $estFragLen2,
        $estFragLen3,

        $corr_estFragLen,
        $corr_estFragLen2,
        $corr_estFragLen3,

        $phantomPeak,
        $corr_phantomPeak,
        $argmin_corr,
        $min_corr,
        $NSC,
        $RSC,
        $QualityTag,
        quote($work_dir)
      )
    )
  . ");";
   $sql_processor->("$sql");
}

=head2 create_table
=cut
sub create_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

# Test for 
# - the specificity of the size selection step
# - enrichment. (Poor enrichment would lead to higher backgound leading to lower correlation.)
  
my $sql = <<SQL
CREATE TABLE if not exists `result_set_qc_phantom_peak` (
  `result_set_qc_phantom_peak_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id`        int(10) unsigned,
  `result_set_id` int(10) unsigned NOT NULL,
  `filename` varchar(512) NOT NULL,
  `numReads` int(10) unsigned NOT NULL,
  `estFragLen`       double default NULL,
  `estFragLen2`      double default NULL,
  `estFragLen3`      double default NULL,
--
-- Seems to always be the same as the above three. Should be removed.
-- Check: Should be actual coorelations.
--
  `corr_estFragLen`  double default NULL,
  `corr_estFragLen2` double default NULL,
  `corr_estFragLen3` double default NULL,
-- Is an estimate of the read length
  `phantomPeak` int(10) unsigned NOT NULL,
  `corr_phantomPeak` double default NULL,
--
-- Shiftlength that minimizes correlation, (not important for qc purposes)
-- See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/figure/F5/
-- for correlation graphs.
--
  `argmin_corr` int(10),
--
-- The value at argmin_corr
--
  `min_corr` double default NULL,
--
-- Normalised strand cross correlation coefficient
--
-- Correlation found at the estimated fragment length (what should be in corr_estFragLen) divided by the minimum found correlation. (min_corr)
--
  `NSC`      double default NULL,
--
-- Relative strand cross correlation coefficient
--
-- This is the main statistic.
--
-- Ratio between correlations found from the fragment length and the phantom peak. Both values are corrected by subtracting min_corr.
-- It shows how much greater the correlation from the fragment length is compared to background correlation.
--
  `RSC`      double default NULL,
--
-- Quality values derived from the RSC
--
  `QualityTag` int(10),
  `path` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_phantom_peak_id`),
--  UNIQUE KEY `filename_idx` (`filename`)
  KEY `filename_idx` (`filename`)
);
SQL
;
  $sql_processor->($sql);
}




