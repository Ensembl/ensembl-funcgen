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

convert_big_files.pl

=head1 SYNOPSIS

perl convert_big_files.pl $directory $name_map $chrom_lengths

=head1 DESCRIPTION

Used to change the chrosomome names in all the BigBed and BigWig files
contained in a directory (recursively)

To run, you will need to provide two files:
* A name mapping file: tab delimited, each line contains the name of
a chromsome followed by the replacement name
* A chromosome length file: tab delimited: each line contains a replacement 
chromosome name and the length of that chromosome

You must also ensure that the following executables are on your path:
* bigBedToBed
* bedToBigBed
* wigToBigWig
* bigWigToWig

=cut

use strict;
use File::Temp qw/tempfile/;

main();

sub main {
  scalar @ARGV == 3 || die ("Three arguments needed");

  my ($dir, $name_map, $chrom_lengths) = @ARGV;
  print "Directory=$dir\n";
  print "Chromosome name conversion =$name_map\n";
  print "Chromosome lengths =$chrom_lengths\n";

  foreach my $arg (@ARGV) {
    -e $arg || die("$arg does not exist!\n");
  }
  if (! -d $dir) {
    die ("$dir is not a directory\n");
  }
  if (! -f $name_map) {
    die ("$name_map is not a file\n");
  }
  if (! -f $chrom_lengths) {
    die ("$chrom_lengths is not a file\n");
  }

  my $new_name = read_hash($name_map);
  my @files = split (/\n/,  `find $dir`);
  foreach my $file (@files) {
    if ((substr($file, -3) eq '.bb') || (substr($file, -3) eq '.bw')) {
      convert_big_file($file, $new_name, $chrom_lengths);
    }
  }
}

sub read_hash {
  my ($file_name) = @_;
  my %hash = ();
  open my $fh, "<", $file_name;
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    scalar(@items) == 2 || die("Malformed line in $file_name:\n$line\n");
    $hash{$items[0]} = $items[1];
  }
  close $fh;
  return \%hash;
}

sub convert_big_file {
  my ($file_name, $new_name, $chroms) = @_;
  print "Converting $file_name\n";

  # Dumps content of big file into $unpacked
  my ($unpacked_fh, $unpacked) = tempfile();
  unpack_big_file($file_name, $unpacked);

  my ($converted_fh, $converted) = tempfile();
  # Copies content of $unpacked into $converted, converting names 
  convert_chromosome_names($unpacked_fh, $new_name, $converted_fh);
  close $unpacked_fh;
  unlink $unpacked;

  # Overwrites big file with new one created from $converted
  repack_big_file($converted, $chroms, $file_name);
  close $converted_fh;
  unlink $converted;
}

sub unpack_big_file {
  my ($source, $destination) = @_;
  my $exec = "";
  if (substr($source, -3) eq '.bb') {
    $exec = "bigBedToBed";
  } elsif (substr($source, -3) eq '.bw') {
    $exec = "bigWigToWig";
  } else {
    die "What's the type of $source?\n";
  }
  system("$exec $source $destination") && die("Failed to unpack $source\n");
}

sub repack_big_file {
  my ($source, $chrom_length_file, $destination) = @_;
  my $exec = "";
  if (substr($destination, -3) eq '.bb') {
    $exec = "bedToBigBed";
  } elsif (substr($destination, -3) eq '.bw') {
    $exec = "wigToBigWig";
  } else {
    die "What's the type of $destination?\n";
  }
  system("$exec $source $chrom_length_file $destination") && die("Failed to repack $source into $destination\n");
}

sub convert_chromosome_names {
  my ($in, $new_name, $out) = @_;
  while (my $line = <$in>) {
    chomp $line;
    my @items = split /\t/, $line;
    if ($line =~ /chrom=\S+/) {
      my @pair = split /=/, $1; 
      exists $new_name->{$pair[1]} || die ("Chromosome name $pair[1] not found in mapping file!\n");
      $pair[1] = $new_name->{$pair[1]};
      my $new_pair = join("=", @pair);
      my $new_line = $line;
      $new_line =~ s/$1/$new_pair/;
      print $out $new_line;
    } elsif (scalar(@items) > 2) {
      exists $new_name->{$items[0]} || die ("Chromosome name $items[0] not found in mapping file!\n");
      $items[0] = $new_name->{$items[0]};
      print $out join("\t", @items);
    } else {
      print $out $line;
    }
    print $out "\n";
  }
}
