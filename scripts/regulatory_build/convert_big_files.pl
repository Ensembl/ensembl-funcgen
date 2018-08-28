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

perl convert_big_files.pl $directory --species $SPECIES --from $FROM --to $TO

=head1 DESCRIPTION

Used to change the chrosomome names in all the BigBed and BigWig files
contained in a directory (recursively)

You need to specify the species of interest, what naming scheme (typically ensembl
or UCSC) you wish to convert to and from. If these parameters are not provided
it is assumed that the conversion is from Ensembl to UCSC. 

You must also ensure that the following executables are on your path:
* bigBedToBed
* bedToBigBed
* wigToBigWig
* bigWigToWig

=cut

use strict;
use Getopt::Long;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;

main();

sub main {
  my $options = get_options();
  my ($name_map, $lengths) = get_chrom_data($options);
  my @files = split (/\n/,  `find $options->{dir}`);
  foreach my $file (@files) {
    if ((substr($file, -3) eq '.bb') || (substr($file, -3) eq '.bw')) {
      convert_big_file($file, $name_map, $lengths);
    }
  }
}

sub get_options {
  my %options = ();

  GetOptions(\%options,
    "species=s",
    "dir=s",
    "from=s",
    "to=s",
    "host|h=s",
    "port|P=s",
    "user|u=s",
    "pass|p=s",
  );

  if (!defined $options{dir}) {
    die("Target directory not defined!\n");
  }

  $options{species} ||= 'Human';
  $options{host} ||= 'ensembldb.ensembl.org';
  $options{user} ||= 'anonymous';
  $options{from} ||= 'ensembl';
  $options{to} ||= 'ucsc';
  $options{registry} = 'Bio::EnsEMBL::Registry';
  $options{registry}->load_registry_from_db(
    -host    => $options{host},
    -user    => $options{user},
    -pass    => $options{pass},
    -port    => $options{port},
  );
  
  return \%options;
}

sub get_chrom_data {
  my ($options) = @_;
  my %hash = ();
  my ($fh, $name) = tempfile();

  my $slice_adaptor = $options->{registry}->get_adaptor($options->{species}, 'Core', 'Slice' );
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    my $keys;
    if ($options->{from} eq 'ensembl') {
      $keys = [ $slice->seq_region_name ];
    } else {
      my @array = map { $_->name } @{$slice->get_all_synonyms($options->{from})};
      $keys = \@array;
    }

    if ((!defined $keys) || scalar (@$keys) == 0 ) {
      die("Could not find sysnonym of ".$slice->seq_region_name." in $options->{from}");
    }

    my $value;
    if ($options->{to} eq 'ensembl') {
      $value = $slice->seq_region_name;
    } else {
      $value = $slice->get_all_synonyms($options->{to})->[0]->name;
    }

    if (!defined $value ) {
      die("Could not find sysnonym of ".$slice->seq_region_name." in $options->{to}");
    }

    foreach my $key (@$keys) {
      $hash{$key} = $value;
      print "$key\t=>\t$value\t".$slice->length." bp\n";
      print $fh "$value\t".$slice->length."\n";
    }
  }

  return \%hash, $name;
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
