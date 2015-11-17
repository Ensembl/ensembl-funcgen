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

=head1 load_phantom_peak_file.pl

=head1 SYNOPSIS

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
my $dbname;

my %config_hash = (
  "result_file"   => \$result_file,
  "result_set_id" => \$result_set_id,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'host'            => \$host,
  'dbname'          => \$dbname,
);

my $result = GetOptions(
  \%config_hash,
  'result_set_id=s',
  'result_file=s',
  'dry_run',
  'user=s',
  'pass=s',
  'host=s',
  'dbname=s',
);

die unless(-e $result_file);
die unless($result_set_id);

my @tracking_db_connection_details = (
    -user     => $user,
    -pass     => $pass,
    -host     => $host,
    -dbname   => $dbname,
);
my $logic_name = 'phantom peak quality tools';
my @flagstats_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'run_spp.R',
        -parameters      => '',
        -description     => 'Computes enrichment and quality measures and fragment lengths for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data',
        -display_label   => 'phantom peak quality tools',
        -displayable     => undef,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(@tracking_db_connection_details);
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -DBCONN => $dbc,  
);
my $analysis_adaptor = $dba->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

if (! $analysis && ! $dry_run) {
      $logger->info("No analysis with logic name $logic_name found. Creating one.");
      $analysis = Bio::EnsEMBL::Analysis->new(@flagstats_analysis_details);
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
  
  (
    my $corr_estFragLen,
    my $corr_estFragLen2,
    my $corr_estFragLen3,
  ) = split ',', $estFragLenTriple;
  

  sub quote {
    my $string = shift;
    return '"' . $string . '"'
  }
  
  my $sql = "INSERT INTO result_set_qc_phantom_peak ("
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
  . "QualityTag "
  . ")  VALUES ("
  . (
    join ', ', (
        $result_set_id,
        $analysis_id,
        quote($filename),
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

my $sql = <<SQL
CREATE TABLE `result_set_qc_phantom_peak` (
  `result_set_qc_phantom_peak_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id`        int(10) unsigned,
  `result_set_id` int(10) unsigned NOT NULL,
  `filename` varchar(100) NOT NULL,
  `numReads` int(10) unsigned NOT NULL,
  `estFragLen`       double default NULL,
  `estFragLen2`      double default NULL,
  `estFragLen3`      double default NULL,
  `corr_estFragLen`  double default NULL,
  `corr_estFragLen2` double default NULL,
  `corr_estFragLen3` double default NULL,
  `phantomPeak` int(10) unsigned NOT NULL,
  `corr_phantomPeak` double default NULL,
  `argmin_corr` int(10),
  `min_corr` double default NULL,
  `NSC`      double default NULL,
  `RSC`      double default NULL,
  `QualityTag` int(10),
  PRIMARY KEY (`result_set_qc_phantom_peak_id`),
  UNIQUE KEY `filename_idx` (`filename`)
);
SQL
;
  $sql_processor->($sql);
}




