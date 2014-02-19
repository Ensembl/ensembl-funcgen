#!/usr/bin/env perl

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
