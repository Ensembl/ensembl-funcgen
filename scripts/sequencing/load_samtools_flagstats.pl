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

=head1 load_samtools_flagstats.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

  ./scripts/sequencing/load_samtools_flagstats.pl \
    --result_set_id 3 \
    --flagstats_file /lustre/scratch109/ensembl/funcgen/mn1/qc/BR1_H3K27me3_3526_bwa_samse_1_2_3.alignment.log \
    --dry_run \
    --user ensro --host ens-genomics2 --dbname mn1_faang_tracking_homo_sapiens_funcgen_81_38 \

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $flagstats_file;
my $result_set_id;
my $dry_run;
my $user;
my $pass;
my $port;
my $host;
my $dbname;
my $work_dir;
my $bam_file;

my %config_hash = (
  'flagstats_file'  => \$flagstats_file,
  'result_set_id'   => \$result_set_id,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'port'            => \$port,
  'host'            => \$host,
  'dbname'          => \$dbname,
  'work_dir'        => \$work_dir,
  'bam_file'        => \$bam_file,
);

# Loading command line paramters into variables and into a hash.
my $result = GetOptions(
  \%config_hash,
  'flagstats_file=s',
  'result_set_id=s',
  'dry_run',
  'user=s',
  'pass=s',
  'port=s',
  'host=s',
  'dbname=s',
  'work_dir=s',
  'bam_file=s',
);

die unless(-e $flagstats_file);
die unless($result_set_id);

my @tracking_db_connection_details = (
    -user     => $user,
    -pass     => $pass,
    -port     => $port,
    -host     => $host,
    -dbname   => $dbname,
);
my $logic_name = 'flagstats';
my @flagstats_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'samtools',
        -parameters      => '-q 1 -F 4',
        -description     => 'samtools flagstats run for qc purposes',
        -display_label   => 'samtools flagstats',
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

create_flagstats_table({
  sql_processor => $sql_processor
});

create_insert_sql({
  analysis_id   => $analysis_id,
  result_set_id => $result_set_id,
  sql_processor => $sql_processor,
});

exit;

=head2 create_insert_sql
=cut
sub create_insert_sql {

  my $param = shift;
  
  my $analysis_id   = $param->{analysis_id};
  my $result_set_id = $param->{result_set_id};
  my $sql_processor = $param->{sql_processor};  

  open IN, $flagstats_file;
  
  while (my $current_line = <IN>) {
    chomp $current_line;
    my $recognized = $current_line =~ /^(\d+) \+ (\d) (.+)$/;
    if ($recognized) {
      my $qc_passed_reads = $1;
      my $qc_failed_reads = $2;
      my $category        = $3;
      
      # We make an exception here, because flagstats looks like this:
      #
      #   13503488 + 0 in total (QC-passed reads + QC-failed reads)
      #   0 + 0 secondary
      #   0 + 0 supplementary
      #   0 + 0 duplicates
      #   13214414 + 0 mapped (97.86%:-nan%)
      #   0 + 0 paired in sequencing
      #   0 + 0 read1
      #   0 + 0 read2
      #   0 + 0 properly paired (-nan%:-nan%)
      #   0 + 0 with itself and mate mapped
      #   0 + 0 singletons (-nan%:-nan%)
      #   0 + 0 with mate mapped to a different chr
      #   0 + 0 with mate mapped to a different chr (mapQ>=5)
      #
      # If we strip the bits in brackets, then the last two lines will have 
      # identical categories in the database.
      #
      if ($category ne 'with mate mapped to a different chr (mapQ>=5)') {
	# Remove the bits samtools puts in brackets. This makes it 
	# easier to use when it is loaded into the database.
	$category =~ s/\(.+\)//;
	chomp($category);
      }
      
      my $sql = "INSERT INTO result_set_qc_flagstats "
      . "(result_set_id, analysis_id, category, qc_passed_reads, qc_failed_reads, path, bam_file) "
      . "VALUES "
      . "($result_set_id, $analysis_id, '$category', $qc_passed_reads, $qc_failed_reads, '$work_dir', '$bam_file');";
      $sql_processor->($sql);
    } else {
      $logger->debug("Can't parse: " . $current_line . "\n");
    }
  }
  close (IN);
}

=head2 create_flagstats_table
=cut
sub create_flagstats_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

my $sql = <<SQL
CREATE TABLE if not exists `result_set_qc_flagstats` (
  `result_set_qc_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `result_set_id`      int(10) unsigned,
  `analysis_id`        int(10) unsigned,
  `category` varchar(100) NOT NULL,
  `qc_passed_reads`    int(10) unsigned,
  `qc_failed_reads`    int(10) unsigned,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`result_set_qc_id`),
  UNIQUE KEY `name_exp_idx` (`result_set_qc_id`,`category`)
);
SQL
;
  $sql_processor->($sql);
}

