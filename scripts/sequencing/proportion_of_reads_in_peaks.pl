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

=head1 proportion_of_reads_in_peaks.pl

=head1 SYNOPSIS

SWEmbl
======

Convert the signal bam file, if necessary:

time samtools view -S -b /lustre/scratch110/ensembl/funcgen/il4/ersaCTTV/CTTV020_epigenomes_of_cell_lines_PILOT/alignments/homo_sapiens/GRCh38/3526/KU812:TF:BR1_CTCF_3526_bwa_samse_1_10_2_3_4_5_6_7_8_9.bam -o /lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/temp/KU812:TF:BR1_CTCF_3526_bwa_samse_1_10_2_3_4_5_6_7_8_9.really.bam

SIGNAL_BAM_FILE=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/temp/KU812:TF:BR1_CTCF_3526_bwa_samse_1_10_2_3_4_5_6_7_8_9.really.bam

PEAK_FILE=/lustre/scratch110/ensembl/funcgen/il4/ersaCTTV/CTTV020_epigenomes_of_cell_lines_PILOT/output/il4_CTTV_tracking_homo_sapiens_funcgen_81_38/peaks/KU812:TF:BR1_CTCF_3526/SWEmbl_R0005_IDR/KU812:TF:BR1_CTCF_3526_bwa_samse.SWEmbl_R0005_IDR.txt

TEMP_DIR=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/temp/$$
rm -rf $TEMP_DIR

./scripts/sequencing/proportion_of_reads_in_peaks.pl \
  --peak_file $PEAK_FILE \
  --temp_dir $TEMP_DIR \
  --bam_file $SIGNAL_BAM_FILE \
  --peak_caller swembl \
  --feature_set_id 123 \
  --user ensro \
  --host ens-genomics2 \
  --dbname mn1_faang_tracking_homo_sapiens_funcgen_81_38


CCAT
====

PEAK_FILE=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/output/il4_CTTV_tracking_homo_sapiens_funcgen_81_38/peaks/K562:hist:BR1_H3K27me3_3526/ccat_histone/K562:hist:BR1_H3K27me3_3526_bwa_samse.ccat_histone.significant.region

TEMP_DIR=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/temp/$$
rm -rf $TEMP_DIR

SIGNAL_BAM_FILE=/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/testbams/CTTV020_epigenomes_of_cell_lines_PILOT/temp/KU812:TF:BR1_CTCF_3526_bwa_samse_1_10_2_3_4_5_6_7_8_9.really.bam

./scripts/sequencing/proportion_of_reads_in_peaks.pl \
  --peak_file $PEAK_FILE \
  --temp_dir $TEMP_DIR \
  --bam_file $SIGNAL_BAM_FILE \
  --peak_caller ccat \
  --feature_set_id 123 \
  --user ensro \
  --host ens-genomics2 \
  --dbname mn1_faang_tracking_homo_sapiens_funcgen_81_38

=head1 DESCRIPTION


=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;
use File::Path qw(make_path);

my $peak_file;
my $bam_file;
my $temp_dir;
my $bed_file;
my $peak_caller;
my $user;
my $pass;
my $port;
my $host;
my $dbname;
my $feature_set_id;
my $dry_run;
my $mock_computation;

my %config_hash = (
  "peak_file"    => \$peak_file,
  "bam_file"     => \$bam_file,
  "temp_dir"     => \$temp_dir,
  "bed_file"     => \$bed_file,
  "peak_caller"  => \$peak_caller,
  "feature_set_id"  => \$feature_set_id,
  "mock_computation"  => \$mock_computation,
  'dry_run'         => \$dry_run,
  'user'            => \$user,
  'pass'            => \$pass,
  'port'            => \$port,
  'host'            => \$host,
  'dbname'          => \$dbname,
);

my $result = GetOptions(
  \%config_hash,
  'peak_file=s',
  'bam_file=s',
  'temp_dir=s',
  'bed_file=s',
  'peak_caller=s',
  'feature_set_id=s',
  'mock_computation=s',
  'dry_run',
  'user=s',
  'pass=s',
  'port=s',
  'host=s',
  'dbname=s',
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

if (! -e $peak_file) {
  $logger->error("The peak file $peak_file does not exist!\n");
}
if (! -e $bam_file) {
  $logger->error("The signal bam file $bam_file does not exist!\n");
}

$peak_caller = lc($peak_caller);

if ($peak_caller ne 'ccat' && $peak_caller ne 'swembl') {
  $logger->error("The parameter peak_caller must either be set to 'ccat' or 'swembl'! Got: $peak_caller\n");
}
if (! defined $feature_set_id) {
  $logger->error("The parameter feature_set_id was not provided!\n");
}

if ($temp_dir) {
  if (-d $temp_dir) {
    $logger->info("Temporary directory $temp_dir exists.\n");
  } else {
    $logger->info("Temporary directory $temp_dir does not exist, so will create it.\n");
    make_path($temp_dir);
  }
}

my @tracking_db_connection_details = (
    -user     => $user,
    -pass     => $pass,
    -host     => $host,
    -port     => $port,
    -dbname   => $dbname,
);
my $logic_name = 'Proportion of reads in peaks';
my @flagstats_analysis_details = (
        -logic_name      => $logic_name,
        -program         => 'proportion_of_reads_in_peaks.pl',
        -parameters      => '',
        -description     => 'Computation of the proportion of reads in peaks',
        -display_label   => 'Proportion of reads in peaks',
        -displayable     => undef,
);

my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(@tracking_db_connection_details);
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -DBCONN => $dbc,  
);
my $analysis_adaptor = $dba->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

if (! $analysis && ! $dry_run) {
      $logger->info("No analysis with logic name $logic_name found. Creating one.\n");
      $analysis = Bio::EnsEMBL::Analysis->new(@flagstats_analysis_details);
      $analysis_adaptor->store($analysis);
}

my $analysis_id;
if (! $analysis && $dry_run) {
  $analysis_id = 'No analysis inserted in dry run mode';
} else {
  $analysis_id = $analysis->dbID;
}
$logger->info("The analysis_id is: $analysis_id.\n");

if ($bed_file) {
  $logger->info("Writing peaks in bed format as specified on the command line to $bed_file\n", 0, 1);
} else {
  $bed_file = $temp_dir . '/' . 'peaks.bed';
  $logger->info("No file specified on the command line for output, using default $bed_file\n", 0, 1);
}

if ($peak_caller eq 'swembl') {

  convert_swembl_peak_to_bed({
    peak_file => $peak_file,
    bed_file  => $bed_file,
  });

} else {

  convert_ccat_peak_to_bed({
    peak_file => $peak_file,
    bed_file  => $bed_file,
  });

}

(
  my $proportion_of_reads_in_peaks,
  my $num_reads_in_peaks,
  my $num_reads_in_total,
) = run_counts({
  bam_file         => $bam_file,
  bed_file         => $bed_file,
  mock_computation => $mock_computation,
});

$logger->info("Proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total = $proportion_of_reads_in_peaks\n");

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

sub convert_ccat_peak_to_bed {

  my $param = shift;
  
  my $peak_file = $param->{peak_file};
  my $bed_file  = $param->{bed_file};

  open my $IN, $peak_file or die "Could not open file $peak_file $!";

  my @data = parse_ccat_peak_file($IN);

  close ($IN);
  
  #write_bed(*STDOUT, \@data);  
  write_bed_if_not_exists($bed_file, \@data);
}

=head1 parse_ccat_peak_file
=cut
sub parse_ccat_peak_file {
  my $IN = shift;
  
  my @column_caption = qw(seqid summit start end signal_reads ctrl_reads fold fdr);

  # Parse the data
  #
  my @data;
  while (my $current_line = <$IN>) {
    chomp $current_line;  
    my @f = split "\t", $current_line;  
    
    use List::MoreUtils qw( zip );
    my %current_row_data = zip(@column_caption, @f);
    
    push @data, \%current_row_data;
  }
  return @data;
}


=head2 convert_swembl_peak_to_bed
=cut
sub convert_swembl_peak_to_bed {

  my $param = shift;
  
  my $peak_file = $param->{peak_file};
  my $bed_file  = $param->{bed_file};

  open my $IN, $peak_file or die "Could not open file $peak_file $!";

  (
    my $header, 
    my $data
  ) = parse_swembl_peak_file($IN);

  close ($IN);

  write_bed_if_not_exists($bed_file, $data);
}

=head2 write_bed_if_not_exists
=cut
sub write_bed_if_not_exists {

  my $bed_file = shift;
  my $data     = shift;

# Change of plans: Always write the bed file. If the script dies when writing 
# out the bed file, it will use a faulty bed file when it is rerun.
#

#   if (! -e $bed_file) {
    open my $out, ">$bed_file" or die "Could not open file $bed_file $!";
    write_bed($out, $data);
    close($out);
    $logger->info("Written peaks to bed file.\n", 0, 1);
#   } else {
#     $logger->info("Bed file already exists. Reusing.\n");
#   }
}

=head2 run_counts
=cut
sub run_counts {
  my $param = shift;
  
  my $bam_file         = $param->{bam_file};
  my $bed_file         = $param->{bed_file};
  my $mock_computation = $param->{mock_computation};
  
  my $num_reads_in_peaks;
  my $num_reads_in_total;

  if ($mock_computation) {

    $num_reads_in_peaks = 2506364;
    $num_reads_in_total = 165055881;

  } else {
    my $current_cmd;
    
    my $READS_IN_PEAKS_BAM_FILE = "${bam_file}.peaks.bam";
    
    use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd run_system_cmd );

    $current_cmd = "bedtools intersect -abam $bam_file -b $bed_file > $READS_IN_PEAKS_BAM_FILE";
    $logger->info("Running $current_cmd\n", 0, 1);
    run_system_cmd($current_cmd, 0, 1);

    $current_cmd = "samtools view -c $READS_IN_PEAKS_BAM_FILE";
    $logger->info("Running $current_cmd\n", 0, 1);
    $num_reads_in_peaks = run_backtick_cmd($current_cmd);

    $current_cmd = "samtools view -c $bam_file";
    $logger->info("Running $current_cmd\n", 0, 1);
    $num_reads_in_total = run_backtick_cmd($current_cmd);
    
    unlink($READS_IN_PEAKS_BAM_FILE);
  }

  my $proportion_of_reads_in_peaks = $num_reads_in_peaks / $num_reads_in_total;
  return ($proportion_of_reads_in_peaks, $num_reads_in_peaks, $num_reads_in_total);
}

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
    . "total_reads, "
    . "path,"
    . "bam_file"
    . ")  VALUES ("
  . (
    join ', ', (
	$analysis_id,
	$feature_set_id,
	$proportion_of_reads_in_peaks,
	$num_reads_in_total,
	"'$temp_dir'",
	"'$bam_file'"
      )
    )
  . ");";
  $sql_processor->($sql);
}

=head2 write_bed
=cut
sub write_bed {

  my $out  = shift;
  my $data = shift;
  
  foreach my $current_peak (@$data) {
  
    my $start = $current_peak->{'start'};
    
    if ($start < 0) {
      $start = 0;
    }

    $out->print(join "\t", 
      $current_peak->{'seqid'},
      $start,
      $current_peak->{'end'});
    $out->print("\n");
  }
}

=head2 parse_swembl_peak_file
=cut
sub parse_swembl_peak_file {
  my $IN = shift;
  
  my %header;
  my @column_caption;
  my $file_pointer_in_data = undef;

  # Parse the header area and the column captions
  #
  LINE: while (
    (!$file_pointer_in_data)
    && (my $current_line = <$IN>)
  ) {
    chomp $current_line;
    
    my $is_header_line = $current_line =~ /^#(.+$)/;
    
    if ($is_header_line) {
      my $header_line = $1;
      (
	my $key,
	my $value
      ) = split "\t", $header_line;

      $header{$key} = $value;
      next LINE;
    }
    
    my $is_column_caption_line = $current_line =~ /^(Region.+$)/;
    
    if ($is_column_caption_line) {
      my $column_caption_line = $1;
      @column_caption = split "\t", $column_caption_line;
      
      # Hack: This column is missing in the output, so adding it here 
      # manually for now.
      #
      push @column_caption, 'p_value';
      
      $file_pointer_in_data = 1;
      
      # The above would give us this:
      #
      # @column_caption = qw(Region	Start pos.	End pos.	Count	Length	Unique pos.	Score	Ref. count	Max. Coverage	Summit	p_value);
      #
      # CCAT however uses "start" and "end" to denote the start and ends and 
      # "seqid". To get a consistent result, I am overriding the column 
      # captions with this:
      #
      @column_caption = qw(seqid	start	end	Count	Length	Unique pos.	Score	Ref. count	Max. Coverage	Summit	p_value);
      #
      # That way we can use the same sub for writing the bed file later.
      next LINE;
    }
    
    die("Can't parse line $current_line");
  }
  
  use Hash::Util qw( lock_hash );
  lock_hash(%header);

  # Parse the data
  #
  my @data;
  while (my $current_line = <$IN>) {
    chomp $current_line;  
    my @f = split "\t", $current_line;  
    
    use List::MoreUtils qw( zip );
    my %current_row_data = zip(@column_caption, @f);
    lock_hash(%current_row_data);
    push @data, \%current_row_data;
  }
  return (\%header, \@data);
}

=head2 create_table
=cut
sub create_table {

  my $param = shift;
  my $sql_processor = $param->{sql_processor};

my $sql = <<SQL
CREATE TABLE if not exists `feature_set_qc_prop_reads_in_peaks` (
  `feature_set_qc_prop_reads_in_peaks_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id`        int(10) unsigned,
  `feature_set_id` int(10) unsigned NOT NULL,
  `prop_reads_in_peaks`       double default NULL,
  `total_reads`      int(10) default NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`feature_set_qc_prop_reads_in_peaks_id`)
);
SQL
;
  $sql_processor->($sql);
}

