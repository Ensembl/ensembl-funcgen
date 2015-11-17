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

=head1 load_fastqc_summary_file.pl

=head1 SYNOPSIS

drop table `result_set_qc_phantom_peak`;

CREATE TABLE `result_set_qc_phantom_peak` (
  `result_set_qc_phantom_peak_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
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
)

INPUT_SUBSET_ID=3436

TEMP_DIR=/lustre/scratch109/ensembl/funcgen/mn1/qc/${INPUT_SUBSET_ID}
mkdir -p $TEMP_DIR
GZIPPED_FASTQ_FILE=$(mysql -N -B $DB_MYSQL_ARGS -e "select local_url from input_subset_tracking where input_subset_id = ${INPUT_SUBSET_ID}")

echo $GZIPPED_FASTQ_FILE

fastqc -o $TEMP_DIR $GZIPPED_FASTQ_FILE

FASTQC_DIR=$(find $TEMP_DIR -maxdepth 1 -mindepth 1 -type d)

SUMMARY_FILE=$FASTQC_DIR/summary.txt

echo $SUMMARY_FILE

./scripts/sequencing/load_phantom_peak_file.pl --result_set_id $INPUT_SUBSET_ID --result_file /lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/test

=head1 DESCRIPTION

Downloads data for input sets registered in the data tracking database.

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;

my $result_file;
my $result_set_id;

my %config_hash = (
  "result_file"   => \$result_file,
  "result_set_id" => \$result_set_id,
);

my $result = GetOptions(
  \%config_hash,
  'result_set_id=s',
  'result_file=s',
);

die unless(-e $result_file);
die unless($result_set_id);

my $dry_run = 1;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $sql_processor;
if ($dry_run) {
  $sql_processor = sub {
    my $sql = shift;
    $logger->info($sql . "\n");
  };
} else {
#   $sql_processor = sub {
#     my $sql = shift;
#     $dbc->do($sql);
#   };
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
)
SQL
;
  $sql_processor->($sql);
}




