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

=head1 NAME

load_fastqc_summary_file.pl

=head1 SYNOPSIS

./scripts/sequencing/load_fastqc_summary_file.pl --input_subset_id 3 --summary_file /lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/alignments/homo_sapiens/GRCh38/3526/histone_control_fastqc/summary.txt | mysql $DB_MYSQL_ARGS

=head1 DESCRIPTION

Downloads data for input sets registered in the data tracking database.

=cut

use strict;
use Data::Dumper;
use Getopt::Long;

my $summary_file;
my $input_subset_id;

# Mapping of command line paramters to variables
my %config_hash = (
  "summary_file"    => \$summary_file,
  "input_subset_id" => \$input_subset_id,
);

# Loading command line paramters into variables and into a hash.
my $result = GetOptions(
  \%config_hash,
  'input_subset_id=s',
  'summary_file=s',
);

die unless(-e $summary_file);
die unless($input_subset_id);

open IN, $summary_file;

while (my $current_line = <IN>) {
  chomp $current_line;
  #print " - $current_line\n";
  my @f = split "\t", $current_line;
  #print Dumper(\@f);
  my $sql = "INSERT INTO input_subset_qc (input_subset_id,status,title,file_name) VALUES (".$input_subset_id.", '".$f[0]."', '".$f[1]."', '".$f[2]."')";
  
  print "$sql;\n";
}






