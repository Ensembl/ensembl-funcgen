#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

unpack_segmentation.pl

=head1 SYNOPSIS

unpack_segmentation.pl $dir

=head1 DESCRIPTION

Splits a concatenated segmentation into a set of bed files.

The directory given as parameter must contain one Bed file with all
the segmentations and one text file with the emissions.

=cut

use strict;
use File::Copy;

my $directory = $ARGV[0];
mkdir "$directory/split";
my @files = glob "$directory/*.bed";
if (scalar @files != 1) {
  print STDERR "Found ". scalar @files." in the directory, exiting...\n";
  exit 1;
}
my $file = $files[0];

system("cat $file | cut -f1 | grep -v track | sed -e 's/_.*\$//' | uniq | sort | uniq > $directory/split/cells.txt") && die; 
open my $fh, "<", "$directory/split/cells.txt";
while (my $cell = <$fh>) {
  chomp $cell;
  system("cat $file | grep '^${cell}_' | sed -e 's/^${cell}_//' > $directory/split/$cell.bed") && die;
}
close $fh;

unlink "$directory/split/cells.txt";

my @emissions = glob "directory/*.bed";
if (scalar @emissions != 1) {
  print STDERR "Found ". scalar @emissions." emission files, exiting...\n";
  copy($emissions[0], "$directory/split/emissions.txt");
}
