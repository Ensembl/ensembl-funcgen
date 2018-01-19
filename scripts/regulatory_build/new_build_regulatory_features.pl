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

new_build_regulatory_features.pl

=head1 SYNOPSIS

Automatic run:
new_build_regulatory_features.pl -o ./ -a hg38 --chrom_lengths GRCh38.sizes -D $FUNCGEN_DB -h $HOST -P $PORT -u $USER -p $PASS

Where:
  * -o: directory for output
  * -a: UCSC assembly name (for trackHub)
  * --chrom_lengths: tab delimited file, each line contains a chromosome name followed by its length.
  * -h: host server
  * -D: host database
  * -P: host port
  * -u: host login
  * -p: host password

Manual run:
new_build_regulatory_features.pl -o ./ -d dump.txt -a hg38 -l chrom_sizes.txt -t tss.bed -g exons.bed

Where:
  * -o: directory for output
  * -a: UCSC assembly name (for trackHub)
  * -t: Bed file with TSS. Can be Ensembl transcript TSS, CAGE tags... [Overrides database info]
  * --chrom_lengths: tab delimited file, each line contains a chromosome name followed by its length.
  * -g: Bed file with Exons. A BioMart dump would work. [Overrides database info]
  * -dump: dump file (described below) [Added to database datasets]

=head1 DESCRIPTION

Generates regulatory features based on overlaps between segmentations
on different cell types.

In particular you will need a dump file which describes all the available
experimental data. The dump file is tab delimited, it contains two types of
entries:
* ChIPseq peaks
peak  $assay  $cell $location_bed_or_bigBed
The assay type is generally a TF antibody name or DNAse

* Segmentations
segmentation  $name $type $location
  - name: is just a free string, which will be used a directory name (avoid special characters and spaces).
  - type: is ChromHMM or Segway
  - location: directory which contains a bunch of bed files, each bedfile named $celltype.bed. In addition
  the directory must contain an emissions file (emissions*.txt for ChromHMM, *.tab for Segway).

=cut

use strict;
use warnings;
use File::Path qw(mkpath);
use File::Basename;
use File::Temp;
${File::Temp::KEEP_ALL} = 1;
use Storable;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use List::Util qw(sum);

########################################################
## Global constants
########################################################

# A bunch of arbitrary cutoffs
# All of these cutoff are applied to enrichment ratios (i.e. 1 == neutral)
our $weak_cutoff = 2;
our $strong_cutoff = 5;
our $very_strong_cutoff = 10;
# The functional labels are the labels of the build built from segmentation data
our @functional_labels = ('tss', 'proximal', 'distal', 'ctcf');
# The empirical labels are the labels of the build built from ChIP-seq peaks
our @empirical_labels = ('tfbs', 'dnase');
# The labels used in the build:
our @labels = (@functional_labels, @empirical_labels);
# Typical histone marks used to detect repression. Note that the strings are normalised
# via the clean_name function below
our @repressed_marks = ('H3K27ME3','H3K9ME3');
our @open_chromatin_assays = ('DNASE', 'DNASE1');
# These labels are used in annotating and coloring both the build and the segmentations, some are not
# present in both, only in one.
our %COLORS = (
  tss => "255,0,0",
  tfbs => "209,157,0",
  dnase => "255,252,4",
  proximal => "255,105,105",
  distal => "250,202,0",
  ctcf => "10,190,254",
  dead => "225,225,225",
  weak => "141,255,68",
  gene => "0,176,80",
  poised => "192,0,190",
  repressed => "127,127,127",
  na => "255,255,255"
);

# # These states describe the features at the cell-type level
# our %feature_states = (
#   active => 0,
#   poised => 1,
#   repressed => 2,
#   inactive => 3,
#   na => 4
# );

our $start_time = time;

main();

=head2 main

  Description: Top level pipeline.

=cut

sub main {
  print_log("Entering New Regulatory Build pipeline\n");
  # Read command line
  my $options = get_options();
  # Read dump file with file locations and metadata
  get_metadata($options);
  # Compute summaries of TF binding peaks
  print_log("Computing summaries of TF binding peaks\n");
  compute_tf_probs($options);
  # Compute summaries of segmentations
  print_log("Computing summaries of segmentations\n");
  extract_segmentation_state_summaries($options);
  # Select the segmentation states that best match a label
  print_log("Selecting the segmentation states that best match a label\n");
  label_segmentation_states($options);
  # Color the input segmentations based on the label assigments above
  print_log("Coloring the input segmentations based on the label assigments above\n");
  make_segmentation_bedfiles($options);
  # Set cutoffs to optimise the quality of the build (using TF binding as an indicator
  print_log("Setting cutoffs to optimise the quality of the build (using TF binding as an indicator\n");
  set_cutoffs($options);
  print_cutoffs($options);
  # Compute MultiCell features
  print_log("Computing MultiCell features\n");
  compute_regulatory_features($options);
  # Determine which features are active in which cell type
  print_log("Determine which features are active in which epigenome\n");

  # Problem:
  #
  # $options->{segmentations} is not being set!
  #
  print Dumper($options->{segmentations});
  if (@{$options->{segmentations}} == 0) {
    use Carp;
    confess("No segmentations have been set!");
  }
#   die();

  compute_states($options);
  # Prepare trackhub header files
  print_log("Preparing trackhub header files\n");
  make_track_hub($options);
  print_log("Exiting New Regulatory Build pipeline\n");
}

=head2 get_options

  Description: Option reading. Parses the command line, checks the options,
    creates directories, opens databases if necessary
  Returntype: hash ref

=cut

sub get_options {
  my $options = read_command_line();
  check_options($options);

  mkpath $options->{working_dir};
  mkpath $options->{trackhub_dir};

  if (defined $options->{host}) {
    $options->{db_adaptor} = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
	-user   => $options->{user},
	-dbname => $options->{db},
	-host   => $options->{host},
	-pass   => $options->{pass},
	-port   => $options->{port},
	-dnadb_user   => $options->{dnadb_user},
	-dnadb_dbname => $options->{dnadb_name},
	-dnadb_host   => $options->{dnadb_host},
	-dnadb_pass   => $options->{dnadb_pass},
	-dnadb_port   => $options->{dnadb_port}
    );
    $options->{dnadb_adaptor} = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
	-user   => $options->{dnadb_user},
	-dbname => $options->{dnadb_name},
	-host   => $options->{dnadb_host},
	-pass   => $options->{dnadb_pass},
	-port   => $options->{dnadb_port}
    );
  }

  return $options;
}

=head2 read_command_line

  Description: Option reading. Parses the command line
  Returntype: hash ref

=cut

sub read_command_line {
  my %options = ();

  GetOptions(\%options,
    "help=s",
    "tmp|t=s",
    "out|o=s",
    "dump=s",
    "assembly|a=s",
    "chrom_lengths|l=s",
    "genome_length|g=s",
    "tss|t=s",
    "exons|g=s",
    "host|h=s",
    "port|P=s",
    "db|D=s",
    "user|u=s",
    "pass|p=s",
    "dnadb_host=s",
    "dnadb_port=s",
    "dnadb_name=s",
    "dnadb_user=s",
    "dnadb_pass=s",
    "mask=s"
  );

  $options{output_dir} = $options{out};
  $options{working_dir} = $options{tmp};

  if (! exists $options{dump}) {
    die('"dump" is a mandatory parameter and has not been set.');
  }
  if (! -e $options{dump}) {
    die("The file " . $options{dump} . " does not exist!");
  }

  return \%options;
}

=head2 check_options

  Description: Option reading. Checks the options, opens databases if needed.
  Arg1: hash ref
  Returntype: Null

=cut

sub check_options {
  my $options = shift;

  die('Output directory not defined') unless defined $options->{out};
  die('Assembly name not defined') unless defined $options->{assembly};
  die('Chromosome lengths not provided') unless defined $options->{chrom_lengths};

  if (! defined $options->{host}) {
    die('Ensembl dumps not described') unless defined $options->{dump};
    die('TSS file not provided') unless defined $options->{tss};
    die('Exon file not provided') unless defined $options->{exons};
    die('Chromosome length file not provided') unless defined $options->{exons};
  } else {
    die('Database not defined') unless defined $options->{db};
    die('User account not defined') unless defined $options->{user};
    if (!defined $options->{dnadb_name}) {
      $options->{dnadb_name} = $options->{db};
      $options->{dnadb_name} =~ s/funcgen/core/;
    }
    if (!defined $options->{dnadb_user}) {
      $options->{dnadb_user} = $options->{user};
    }
    if (!defined $options->{dnadb_host}) {
      $options->{dnadb_host} = $options->{host};
    }
    if (!defined $options->{dnadb_pass}) {
      $options->{dnadb_pass} = $options->{pass};
    }
    if (!defined $options->{dnadb_port}) {
      $options->{dnadb_port} = $options->{port};
    }
  }

  if (!defined $options->{working_dir}) {
    $options->{working_dir} = "$options->{output_dir}/tmp/";
  }

  $options->{trackhub_dir} = "$options->{output_dir}/$options->{assembly}";
}

=head2 must_compute

  Description: Decides whether a computation is needed
    * As long as all checked files exist, they are not recomputed
    * As soon as a checked file needs to be recomputed, everything
      after that is recomputed
  Returntype: boolean

=cut

sub must_compute {
  my ($options, $file) = @_;
  if (defined $options->{must_compute}) {
    return 1;
  } elsif (!-e $file) {
    $options->{must_compute} = 1;
    return 1;
  } else {
    return 0;
  }
}

=head2 convert_to_bigBed

  Description: BigBed creation
  Arg1: path to file
  Returntype: undef
  Side effects:
    - creates indexed file (BigWig or BigBed)
    - deletes flat file

=cut

sub convert_to_bigBed {
  my ($options, $file) = @_;
  my $new = $file;
  $new =~ s/\.bed$/.bb/;
  run("bedToBigBed $file $options->{chrom_lengths} $new") ;
  unlink $file;
}

=head2 convert_to_bigWig

  Description: BigWig creation
  Arg1: path to file
  Returntype: undef
  Side effects:
    - creates indexed file (BigWig or BigWig)
    - deletes flat file

=cut

sub convert_to_bigWig {
  my ($options, $file) = @_;
  trim_bed_to_chrom_lengths($options, $file);
  my $new = $file;
  $new =~ s/\.wig$/.bw/;
  run("wigToBigWig $file $options->{chrom_lengths} $new") ;
  unlink $file;
}

=head2 run

  Description: Wrapper function for system calls
  Arg1: Command line
  Returntype: undef
  Side effects: Runs command, prints out error in case of failure

=cut

sub run {
  my ($cmd) = @_;
  print_log("Running $cmd\n");
  my $exit_code = system($cmd);
  if ($exit_code != 0) {
    use Carp;
    confess("Failure when running command\n$cmd\n")
  }
}

=head2 print_log

  Description: Print conveninence
  Arg1: String
  Returntype: undef
  Side effects: Prints string, with time stamp in front

=cut

sub print_log {
  my ($str) = @_;
  my $runtime = time - $start_time;
  print "[$runtime] $str";
}

=head2 clean_name

  Description: Removing unwanted characters for safe filenames
  Arg1: String
  Returntype: String

=cut

sub clean_name {
  my $string = shift;
  $string =~ s/[\-\(\)]//g;
#   $string =~ s/_.*//g;
# $string =~ s/ /_/g;
# $string =~ s/\+/_/g;
  $string = uc($string);
  $string =~ s/:/x/g;
  return $string;
}

=head2 deserialise

  Description: Reads a list of strings from a text file
  Arg1: filename
  Returntype: list ref

=cut

sub deserialise {
  my ($filename) = shift;
  my @data = ();
  my $fh;
  open $fh, "<", $filename;
  my $line;
  while ($line = <$fh>) {
    chomp $line;
    push @data, $line;
  }
  return \@data;
}

=head2 deserialise_hash

  Description: Reads a hash from a text file
  Arg1: filename
  Returntype: hash ref

=cut

sub deserialise_hash {
  my ($filename) = shift;
  my %data = ();
  my $fh;
  open $fh, "<", $filename;
  my $line;
  while ($line = <$fh>) {
    chomp $line;
    my @items = split("\t", $line);
    $data{$items[0]} = $items[1];
  }
  return \%data;
}

=head2 trim_bed_to_chrom_lengths

  Description: Trimming to chromosome length
    Ensures that all features of a bed file fit within the
    exact coordinates of a chromosome
    This is necessary for BigFile indexing, as ChromHMM
    rounds up coordinates to the nearest 200bp mark.
  Arg1: options hash
  Arg2: filepath
  Returntype: undef
  Side effects: Overwrites file with corrected coordinates

=cut

sub trim_bed_to_chrom_lengths {
  my ($options, $file) = @_;
  if (!defined $options->{length_hash}) {
    $options->{length_hash} = read_chrom_lengths($options);
  }
  my $lengths = $options->{length_hash};
  my ($fh, $out, $line);
  open $fh, "<", $file;
  open $out, ">", "$file.tmp";
  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    if (defined $lengths->{$items[0]} && $items[1] < $lengths->{$items[0]} && $items[2] > 0) {
      if ($items[2] >= $lengths->{$items[0]}) {
        $items[2] = $lengths->{$items[0]};
      }
      if ($items[1] < 0) {
        $items[1] = 0;
      }
      print $out join("\t", @items)."\n";
    }
  }
  close $out;
  close $fh;

  unlink $file;
  rename "$file.tmp", $file;

  if (-z $file) {
    print_log("WARNING: Filtered $file is now empty... is this normal?");
  }
}

=head2 read_chrom_lengths

  Description: reads chromosome lengths from flat file
  Arg1: options hash ref
  Returntype: hashref

=cut

sub read_chrom_lengths {
  my ($options) = @_;
  my ($fh, $line);

  open $fh, "<", $options->{chrom_lengths};
  my %lengths = ();

  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    $lengths{$items[0]} += $items[1];
  }
  close $fh;

  return \%lengths;
}

=head2 get_metadata

  Description: Parsing Ensembl dumps
    The Ensembl dumpfile is assumed to be a tab delimited
    file with two types of lines:
    peaks  $TF  $celltype  $path
    segmentation  $type  $path
    The type of a segmentation is currently only ChromHMM
    or Segway
    but other segmentation software could crop up.
    The path of a segmentation is the path to the directory
    containing all the output files.
  Arg1: options hash
  Returntype: undef
  Side effects: Fills up options hash

=cut

sub get_metadata {
  my $options = shift;

  my $location = "$options->{working_dir}/metadata";

  if (must_compute($options, $location)) {
    $options->{cell_type_tfs} = {};
    $options->{cell_type_open} = {};
    $options->{peak_calls} = {};
    $options->{segmentations} = [];

    read_dump($options);

    if (! defined $options->{genome_length}) {
      $options->{genome_length} = compute_genome_length($options);
    }
    warn("The genome length is: " . $options->{genome_length} . "\n");
    
    if (defined $options->{db_adaptor}) {
      fetch_metadata($options);
    }

    if (!defined $options->{tss}) {
      create_tss($options);
    }

    if (!defined $options->{exons}) {
      create_exons($options);
    }

    if (!defined $options->{mask} && defined $options->{host}) {
      create_mask($options);
    }
    my %hash = ();
    $hash{cell_type_tfs} = $options->{cell_type_tfs};
    $hash{cell_type_open} = $options->{cell_type_open};
    $hash{peak_calls} = $options->{peak_calls};
    $hash{segmentations} = $options->{segmentations};
    $hash{chrom_lengths} = $options->{chrom_lengths};
    $hash{tss} = $options->{tss};
    $hash{exons} = $options->{exons};
    $hash{mask} = $options->{mask};
    $hash{genome_length} = $options->{genome_length};
    store \%hash, $location;
  } else {
    my $prev_options = retrieve($location);
    $options->{cell_type_tfs} = $prev_options->{cell_type_tfs};
    $options->{cell_type_open} = $prev_options->{cell_type_open};
    $options->{peak_calls} = $prev_options->{peak_calls};
    $options->{segmentations} = $prev_options->{segmentations};
    $options->{chrom_lengths} = $prev_options->{chrom_lengths};
    $options->{tss} = $prev_options->{tss};
    $options->{exons} = $prev_options->{exons};
    $options->{mask} = $prev_options->{mask};
    $options->{genome_length} = $prev_options->{genome_length};
  }
}

=head2 read_dump

  Description: reads content of metadata file
  Arg1: options hash ref
  Returntype: undef
  Side effects: stores location of segmentation directories and chip-seq
    peak files in options

=cut

sub read_dump {
  my ($options) = @_;
  print_log("Reading $options->{dump}\n");
# die;
  open my $fh, "<", $options->{dump};
  while (my $line = <$fh>) {
    chomp $line;
    my @elems = split /\t/, $line;

    if ($elems[0] eq 'peaks') {
      record_peak_file($options, $elems[1], $elems[2], $elems[3]);
    } elsif ($elems[0] eq 'segmentation') {
      my $segmentation = {};
      $segmentation->{name} = $elems[1];
      $segmentation->{type} = $elems[2];
      $segmentation->{location} = $elems[3];

      push @{$options->{segmentations}}, $segmentation;
    }
  }
#   use Data::Dumper;
#   print Dumper($options->{segmentations});
#   die;
}

=head2 fetch_metadata

  Description: pulls out useful data from database
  Arg1: options hash ref
  Returntype: undef
  Side effects: Stores data into options

=cut

sub fetch_metadata {
  my ($options) = @_;
  print_log("Fetching metadata from funcgen DB\n");
  my @slices = sort {$a->seq_region_name cmp $b->seq_region_name} @{$options->{dnadb_adaptor}->get_SliceAdaptor->fetch_all('toplevel', undef, undef, 0)};
  foreach my $featureSet (@{$options->{db_adaptor}->get_adaptor("FeatureSet")->fetch_all_by_feature_class("annotated")}) {
    my $tf = $featureSet->feature_type->name;
    if ($tf =~ /^H[2-4][ABKZ]/) {
      # It's a histone, ignore
      next;
    }
    my $cell = $featureSet->epigenome->production_name;

    my $cell_display_label = $featureSet->epigenome->display_label;

    my $epigenome_is_excluded =
         ($cell_display_label eq "CD38- naïve B cell (CB)")
      || ($cell_display_label eq "CD38- naive B cell (VB)")
      || ($cell_display_label eq "CD4+ ab T cell (CB)")
      || ($cell_display_label eq "CD8+ ab T cell (VB)")
      || ($cell_display_label eq "EM CD8+ ab T cell (VB)")
      || ($cell_display_label eq "Naïve B cell (To)")
    ;
    if ($epigenome_is_excluded) {
      print_log("Skipping $cell_display_label, because it has been excluded from the regulatory build.\n");
      next;
    }

    record_peak_file($options, $tf, $cell, dump_peaks($options, $featureSet, $tf, $cell, \@slices));
    # close $fh;
  }
}

=head2 dump_peaks

  Description: dumps feature set into flat file
  Arg1: Bio::EnsEMBL::FeatureSet Object
  Arg2: TF name (string)
  Arg3: CellType name (String, processed by clean_name)
  Arg4: Array ref of Bio::EnsEMBL::Slice objects
  Returntype: String, file location
  Side effects: creates and fills file

=cut

sub dump_peaks {
  my ($options, $featureSet, $tf, $cell, $slices) = @_;
  my $dir = "$options->{working_dir}/peaks/$tf/$cell/";
  mkpath $dir;

  my $fh = File::Temp->new(DIR => $dir, SUFFIX => '.bed');
  my $filename = $fh->filename;
  foreach my $slice (@$slices) {
    foreach my $feature (sort {$a->start <=> $b->start} @{$featureSet->get_Features_by_Slice($slice)}) {
      print $fh join("\t", ($slice->seq_region_name, $feature->start - 1, $feature->end))."\n";
    }
  }

  return $fh->filename;
}

=head2 record_peak_file

  Description: stores location of Chip-Seq peak file
  Arg1: options hash ref
  Arg2: TF or histone mark name (string)
  Arg3: CellType name
  Arg4: File location
  Returntype: undef
  Side effects: stores location into relevant hash refs

=cut

sub record_peak_file {
  my ($options, $tf, $cell, $file) = @_;
  my $ctf = clean_name($tf);
  my $ccell = clean_name($cell);

  if (grep($_ eq $ctf, @open_chromatin_assays)) {
    if (!defined $options->{cell_type_open}->{$ccell}) {
      $options->{cell_type_open}->{$ccell} = [];
    }
    push @{$options->{cell_type_open}->{$ccell}}, $file;
  } else {
    if (!defined $options->{cell_type_tfs}->{$ccell}) {
      $options->{cell_type_tfs}->{$ccell} = [];
    }
    push @{$options->{cell_type_tfs}->{$ccell}}, $file;

    if (!defined $options->{peak_calls}->{$ctf}) {
      $options->{peak_calls}->{$ctf} = {};
    }
    if (!defined $options->{peak_calls}->{$ctf}->{$ccell}) {
      $options->{peak_calls}->{$ctf}->{$ccell} = [];
    }
    push @{$options->{peak_calls}->{$ctf}->{$ccell}}, $file;
  }
}

#
# Removed, because for human this returns three entries for the Y chromosome.
#
# Two of these are short. When these are used later for trimming bed files, 
# feaatures are lost.
#

# =head2 create_chrom_lengths
# 
#   Descriptions: store chromosome lengths from DB into flat file for BigFile compression
#   Arg1: options hash ref
#   Returntype: undef
#   Side effects: Stores location of newly created file in options
# 
# =cut
# 
# sub create_chrom_lengths {
#   my ($options) = @_;
#   $options->{chrom_lengths} = "$options->{working_dir}/chrom_lengths.txt";
#   open my $fh, ">", $options->{chrom_lengths};
#   fetch_chrom_lengths($options, $fh);
#   close $fh
# }

# =head2 fetch_chrom_lengths
# 
#   Description: stores chromosome lengths from DB into filehandle
#   Arg1: options hash ref
#   Arg2: file handle
#   Returntype: undef
#   Side effects: Writes into file handle
# 
# =cut
# 
# sub fetch_chrom_lengths {
#   my ($options, $fh) = @_;
#   print_log("Fetching chromosome lengths from core DB\n");
#   my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
#   my @slices = @{ $slice_adaptor->fetch_all('toplevel', undef, undef, 0) };
# 
#   foreach my $slice (@slices) {
#     print $fh join("\t", ($slice->seq_region_name(), $slice->end() - $slice->start())) . "\n";
#     $options->{genome_length} += $slice->length();
#   }
#    print_log("$options->{assembly} length: $options->{genome_length} \n");
# }

=head2 compute_genome_length

  Description: computes genome length
  Arg1: options hash ref
  Returntype: undef

=cut
sub compute_genome_length {

  my $options = shift;
  
  my $core_database_available = exists $options->{dnadb_adaptor};
  my $genome_length;

  if ($core_database_available) {
    $genome_length = compute_genome_length_by_querying_the_core_database($options);
  } else {
    die(
      "You have not specified the genome length on the command line with "
      . "the -genome_length parameter and no core database to "
      . "fetch this from!\n\n"
      . "The genome length can be computed by summing the lengths from your "
      . "chromosome length file, if you comment out this line from the script."
    );
    warn(
      "The genome length will be computed from the chromosome length file, "
      . "because no core database was provided and no genome length was "
      . "specified on the command line.\n"
    );
    $genome_length = compute_genome_length_from_chromosome_length_file($options);
    warn("The genome length computed is: $genome_length\n");
  }
  return $genome_length;
}

sub compute_genome_length_by_querying_the_core_database {
  my $options = shift;
  
  my $dnadb_adaptor = $options->{dnadb_adaptor};
  
  my $genome_container = $dnadb_adaptor->get_adaptor('GenomeContainer');
  my $genome_length    = $genome_container->get_ref_length;

  return $genome_length;
}

=head2 compute_genome_length_from_chromosome_length_file

  Description: computes genome length by adding the chromosome lengths
  Arg1: options hash ref
  Returntype: undef

=cut
sub compute_genome_length_from_chromosome_length_file {
  my $options = shift;

  if (!defined $options->{length_hash}) {
    $options->{length_hash} = read_chrom_lengths($options);
  }
  my @genome_sequences_lengths = values %{$options->{length_hash}};
  
  my $genome_length = 0;
  foreach my $current_genome_sequence_length (@genome_sequences_lengths) {
    $genome_length += $current_genome_sequence_length;
  }

  return $genome_length;
}

=head2 create_tss

  Description: stores TSS features from DB into flat file
  Arg1: options hash ref
  Returntype: undef
  Side effects: Stores location of new file into options hash ref

=cut

sub create_tss {
  my ($options) = @_;
  
  my $transcription_start_sites_file = "$options->{working_dir}/tss.bed";
  $options->{tss} = $transcription_start_sites_file;
  open my $fh, ">", $transcription_start_sites_file;
  fetch_tss($options, $fh);
  close $fh;
  if (! -e $transcription_start_sites_file) {
    use Carp;
    confess("The transcription start sites file ($transcription_start_sites_file) was not created!");
  }
}

=head2 fetch_tss

  Description: writes TSS locations from DB into filehandle
  Arg1: options hash ref
  Arg2: file handle
  Returntype: None
  Side effects: writes into filehandle,

=cut

sub fetch_tss {
  my ($options, $fh) = @_;
  print_log("Fetching TSSs from core DB\n");
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  my @tss_coords = ();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, 0) }) {
    foreach my $gene (@{$slice->get_all_Genes()}) {
      foreach my $transcript (@{$gene->get_all_Transcripts()}) {
	if ($transcript->strand() > 0) {
	  push @tss_coords, [$slice->seq_region_name(), $transcript->start() - 1, $transcript->start()];
        } else {
	  push @tss_coords, [$slice->seq_region_name(), $transcript->end() - 1, $transcript->end()];
	}
      }
    }
  }
  write_into_bedfile(\@tss_coords, $fh);
}

=head2 create_exons

  Description: Stores exon locations from DB into flatfile
  Arg1: options hash ref
  Returntype: undef
  Side effects: stores file location into options

=cut

sub create_exons {
  my ($options) = @_;
  $options->{exons} = "$options->{working_dir}/exons.bed";
  open my $fh, ">", $options->{exons};
  fetch_exons($options, $fh);
  close $fh
}

=head2 fetch_exons

  Description: Writes exon locations from DB into filehandle
  Arg1: options hash ref
  Arg2: filehandle
  Returntype: undef
  Side effects: writes into filehandle

=cut

sub fetch_exons {
  my ($options, $fh) = @_;
  my @exon_coords = ();
  print_log("Fetching exons from core DB\n");
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, 0) }) {
    foreach my $gene (@{$slice->get_all_Genes()}) {
      foreach my $transcript (@{$gene->get_all_Transcripts()}) {
        foreach my $exon (@{$transcript->get_all_Exons()}) {
	  push @exon_coords, [$slice->seq_region_name(), $exon->start() - 1, $exon->end()];
	}
      }
    }
  }
  write_into_bedfile(\@exon_coords, $fh);
}

=head2

  Description: stores masked regions into flat file
  Arg1: options hash ref
  Returntype: undef
  Side effects: stores file location into options

=cut

sub create_mask {
  my ($options) = @_;
  $options->{mask} = "$options->{working_dir}/mask.bed";
  open my $fh, ">", $options->{mask};
  fetch_mask($options, $fh);
  close $fh
}

=head2 fetch_mask

  Descriptions: writes masked regions from DB into file handle
  Arg1: options hash ref
  Arg2: filehandle
  Returntype: undef
  Side effects: writes into filehandle

=cut

sub fetch_mask {
  my ($options, $fh) = @_;
  my @mask_coords = ();
  print_log("Fetching ENCODE excluded regions from core DB\n");
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, 0) }) {
    foreach my $mask (@{$slice->get_all_MiscFeatures('encode_excluded')}) {
      push @mask_coords, [$slice->seq_region_name(), $mask->start() - 1, $mask->end()];
    }
  }
  write_into_bedfile(\@mask_coords, $fh);
}

=head2 write_into_bedfile

  Description: writes regions into sorted bed file
  Arg1: Regions: list ref of list refs which each contain [chromosoms, start, end, ...]
  Arg2: Filehandle
  Returntype: undef
  Side effects: write into filehandle

=cut

sub write_into_bedfile {
  my ($regions, $fh) = @_;
  my @sorted_regions = sort {comp_coords($a, $b)} @$regions;

  foreach my $region (@sorted_regions) {
    if ($region->[1] < $region->[2]) {
      print $fh join("\t", @{$region})."\n";
    }
  }
}

=head2 comp_coords

  Description: Simple comparator function for sorting regions
    Works on list refs, containing [chromosoms, start, end, ...]
  Arg1: $a: first element of comparison
  Arg2: $b: second element of comparison
  Returntype: -1 if $a before $b, 0 if overlap, 1 if $a after $b

=cut

sub comp_coords {
  my ($a, $b) = @_;
  if ($a->[0] ne $b->[0]) {
    return $a->[0] cmp $b->[0];
  } elsif ($a->[1] != $b->[1]) {
    return $a->[1] <=> $b->[1];
  } else {
    return $a->[2] <=> $b->[2];
  }
}

=head2 compute_tf_probs

  Description: Computing TF binding probs
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Returntype: undef
  Side effects: Writes into:
   * Celltype specific summaries in:
     $options->{working_dir}/celltype_tf/$celltype.bed
   * TF specific summaries in:
     $options->{trackhub_dir}/tfbs/$tf.wig
   * An overall TF summary, computed as a disjunction of the
   TF specific probabilities, treated as independent variables.
     $options->{trackhub_dir}/overview/all_tfbs.bw

=cut

sub compute_tf_probs {
  my $options = shift;

  if (!must_compute($options,"$options->{trackhub_dir}/overview/all_tfbs.bw")) {
    print_log("TF binding probs already there, skipping calculations...\n");
    return;
  }

  print_log("Computing TF binding tracks\n");
  compute_celltype_tf_sites($options);
  my $tf_probs = compute_antibody_specific_probs($options);
  
  if (@$tf_probs == 0) {
    use Carp;
    confess('No transcription factor probabilities found!');
  }
  
  compute_global_tf_prob($options, $tf_probs);
}

=head2 compute_celltype_tf_sites

  Description: Computing TF binding sites for each cell type
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Returntype: undef
  Side effects: Writes into:
   * Celltype specific summaries in:
     $options->{working_dir}/celltype_tf/$celltype.bed

=cut

sub compute_celltype_tf_sites {
  my ($options) = @_;
  my $celltype;

  mkdir "$options->{working_dir}/celltype_tf/";
  mkdir "$options->{working_dir}/celltype_dnase/";

  foreach $celltype (keys %{$options->{cell_type_tfs}}) {
    compute_celltype_tf_sites_2($options, $celltype);
  }
}

=head2 compute_celltype_tf_sites_2

  Description: Computing TF binding sites for a given cell type
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Arg2: CellType name (String, processed by clean_name)
  Returntype: undef
  Side effects: Writes into:
   * Celltype specific summaries in:
     $options->{working_dir}/celltype_tf/$celltype.bed

=cut

sub compute_celltype_tf_sites_2 {
  my ($options, $celltype) = @_;

  my $output2 = "$options->{working_dir}/celltype_tf/$celltype.bed";
  my $tf_exps = $options->{cell_type_tfs}->{$celltype};
  my $open_exps = $options->{cell_type_open}->{$celltype};

  if (defined $open_exps) {
    my $output1 = "$options->{working_dir}/celltype_dnase/$celltype.bed"; 
#     run("wiggletools write_bg $output1 unit sum " . join(" ", @{$open_exps}));
    run_unless_file_exists("wiggletools write_bg $output1 unit sum " . join(" ", @{$open_exps}), $output1);
    $options->{celltype_dnase}->{$celltype} = $output1;
    trim_bed_to_chrom_lengths($options, $output1);
#     run("wiggletools write_bg $output2 gt 1 sum ".join(" ", @{$tf_exps}). " $output1 ");
    run_unless_file_exists("wiggletools write_bg $output2 gt 1 sum ".join(" ", @{$tf_exps}). " $output1 ", $output2);
    
  } else {
#     run("wiggletools write_bg $output2 unit sum ".join(" ", @{$tf_exps}));
    run_unless_file_exists("wiggletools write_bg $output2 unit sum ".join(" ", @{$tf_exps}), $output2);
  }
  trim_bed_to_chrom_lengths($options, $output2);
  $options->{celltype_tfbs}->{$celltype} = $output2;
}

sub run_unless_file_exists {
  my $cmd  = shift;
  my $file = shift;
  
  if (! -e $file) {
    run($cmd);
  }
}

=head2 compute_antibody_specific_probs

  Description: Computing TF binding summaries for each TF
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Returntype: array ref of file locations
  Side effects: Writes into:
   * TF specific summaries in:
     $options->{trackhub_dir}/tfbs/$tf.wig

=cut

sub compute_antibody_specific_probs {
  my $options = shift;
  my @tfs = keys %{$options->{peak_calls}};
  my @tf_probs = ();
  my $tf;

  mkdir "$options->{trackhub_dir}/tfbs/";

  foreach $tf (@tfs) {
    my $res = compute_antibody_specific_prob($options, $tf);
    if (defined $res && !($tf eq "CTCF")) {
      push @tf_probs, $res;
    }
  }

  return \@tf_probs;
}

=head2 compute_antibody_specific_prob

  Description: Computing TF binding summaries for a given TF
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Arg2: TF name (string process by clean_name)
  Returntype: array ref of file locations
  Side effects: Writes into:
   * TF specific summaries in:
     $options->{trackhub_dir}/tfbs/$tf.wig

=cut

sub compute_antibody_specific_prob {
  my ($options, $tf) = @_;
  my $output = "$options->{trackhub_dir}/tfbs/$tf.wig";

  my @exps = @{antibody_regions($options, $tf)};
  my $exp_count = scalar(@exps);
  if ($exp_count > 0) {
    run("wiggletools write_bg $output scale " . (1 / $exp_count) . " sum " . join(" ", @exps));
    convert_to_bigWig($options, $output);
    $output =~ s/.wig/.bw/;
    return $output;
  } else {
    return undef;
  }
}

=head2 antibody_regions

  Description: create summary string of antibody regions
  Arg1: hash ref with values:
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Arg2: TF name (string process by clean_name)
  Returntype: wiggletools expression, string

=cut

sub antibody_regions {
  my ($options, $tf) = @_;
  my @cells = keys(%{$options->{peak_calls}->{$tf}});
  my @exps = ();

  foreach my $cell (@cells) {
     if (defined $options->{cell_type_open}->{$cell}) {
       my $tf_peaks = $options->{peak_calls}->{$tf}->{$cell};
       my $open_peaks = $options->{cell_type_open}->{$cell};
       push @exps, "gt 1 sum " . join(" ", @$tf_peaks) . " unit sum " . join(" ", @$open_peaks) . " : :";
     }
  }

  return \@exps;
}

=head2 compute_global_tf_prob

  Description: Computing global TF binding prob
  Arg1: hash ref with values:
   - $options->{trackhub_dir} (path to trackhub dir)
   - $options->{working_dir} (path to working tmp dir)
   - $options->{cell_type_open} (hash which assigns a list of peak calls to each cell type)
   - $options->{cell_type_tfs} (hash which assigns a list of peak calls to each cell type)
   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
  Arg2: array ref of file locations containing TF specific binding probs in BigWig format
  Returntype: undef
  Side effects: Writes an overall TF summary, computed as a disjunction of the
   TF specific probabilities, treated as independent variables.
     $options->{trackhub_dir}/overview/all_tfbs.bw

=cut

sub compute_global_tf_prob {
  my ($options, $probs) = @_;
  mkdir "$options->{trackhub_dir}/overview/";
  my $output = "$options->{trackhub_dir}/overview/all_tfbs.wig";

  # Probability of seeing at least one transcription factor
  #
#   run("wiggletools write_bg $output offset 1 scale -1 mult map offset 1 map scale -1 " . join(" ", @$probs)) ;

  if (@$probs == 0) {
    use Carp;
    confess("Probabilities are empty!");
  }

  run_unless_file_exists("wiggletools write_bg $output offset 1 scale -1 mult map offset 1 map scale -1 " . join(" ", @$probs), $output) ;

  convert_to_bigWig($options, $output);
}

=head2 extract_segmentation_state_summaries

  Description: Extracting segmentation summaries
    Computes for each segmentation and each state the
    number of celltypes with a given state at each position
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Returntype: undef
  Side effects: It fills out:
    * $segmentation->{states} List of strings
    * $segmentation->{celltypes} List of strings
    Writes into:
    * List of states in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt
    * List of cell-types in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt
    * Sorted copy of each cell type's segmentation:
      $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed
    * A summary statistic of each state:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw

=cut

sub extract_segmentation_state_summaries {
  my $options = shift;
  my $segmentation;

  mkdir "$options->{trackhub_dir}/segmentation_summaries/";
  mkdir "$options->{working_dir}/segmentation_summaries/";

  foreach $segmentation (@{$options->{segmentations}}) {
    extract_segmentation_state_summaries_2($options, $segmentation);
  }
}

=head2 extract_segmentation_state_summaries_2

  Description: Extracting segmentation summaries
    Computes for a given segmentation and each state the
    number of celltypes with a given state at each position
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Arg2: Segmentation name (string)
  Returntype: undef
  Side effects: It fills out:
    * $segmentation->{states} List of strings
    * $segmentation->{celltypes} List of strings
    Writes into:
    * List of states in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt
    * List of cell-types in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt
    * Sorted copy of each cell type's segmentation:
      $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed
    * A summary statistic of each state:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw

=cut

sub extract_segmentation_state_summaries_2 {
  my ($options, $segmentation) = @_;

  mkdir "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/";
  mkdir "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/";

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway' || $segmentation->{type} eq 'GMTK') {
    extract_ChromHMM_state_summaries($options, $segmentation);
  } else {
    die("Unknown segmentation format $segmentation->{type}\n");
  }
}

=head2 extract_ChromHMM_state_summaries

  Description: Extracting segmentation summaries
    Computes for a given ChromHMM, Segway or GMTK segmentation and each state the
    number of celltypes with a given state at each position
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Arg2: Segmentation name (string)
  Returntype: undef
  Side effects: It fills out:
    * $segmentation->{states} List of strings
    * $segmentation->{celltypes} List of strings
    Writes into:
    * List of states in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt
    * List of cell-types in segmentation:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt
    * Sorted copy of each cell type's segmentation:
      $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed
    * A summary statistic of each state:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw

=cut

sub extract_ChromHMM_state_summaries {
  my ($options, $segmentation) = @_;
  print_log("Going through output of segmentation $segmentation->{name}\n");
  my @bedfiles;
  if ($segmentation->{type} eq 'ChromHMM') {
    @bedfiles = glob "$segmentation->{location}/*_segments.bed";
  } else {
    @bedfiles = glob "$segmentation->{location}/*.bed";
  }
  $segmentation->{states} = extract_ChromHMM_states($options, $segmentation, \@bedfiles);
  $segmentation->{celltypes} = extract_ChromHMM_cells($options, $segmentation, \@bedfiles);

  foreach my $cell (keys %{$segmentation->{celltypes}}) {
    my $summary = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed";
    if (must_compute($options,$summary)) {
      run("sort -k1,1 -k2,2n $segmentation->{celltypes}->{$cell} > $summary");
      trim_bed_to_chrom_lengths($options, $summary);
    }
  }

  foreach my $state (@{$segmentation->{states}}) {
    extract_ChromHMM_state_summary($options, $segmentation, $state);
  }
}

=head2 extract_ChromHMM_states

  Description: Extracting list of states in a ChromHMM segmentation
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Arg2: Segmentation name (string)
  Arg3: Array ref of file locations containing cell type specific segmentation data
  Returntype: Array ref

=cut

sub extract_ChromHMM_states {
  my ($options, $segmentation, $files) = @_;
  print_log("Extracting state names\n");
  my $output = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt";
  if (must_compute($options,$output)) {
    run("cat " . join(" ", @$files) . " | grep -v '^track\>' | cut -f4 | uniq | sort | uniq > $output");
  }
  return deserialise($output);
}

=head2 extract_ChromHMM_states

  Description: Extracting list of cell types in a ChromHMM segmentation
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Arg2: Segmentation name (string)
  Arg3: Array ref of file locations containing cell type specific segmentation data
  Returntype: Array ref

=cut

sub extract_ChromHMM_cells {
  my ($options, $segmentation, $files) = @_;
  print_log("Extracting cell names\n");
  my $output = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt";
  if (must_compute($options,$output)) {
    open my $fh, ">", $output;
    foreach my $file (@$files) {
      my $cell = basename $file;
      $cell =~ s/.bed$//;
      $cell = clean_name($cell);
      print $fh "$cell\t$file\n";
    }
    close $fh;
  }
  return deserialise_hash($output);
}

=head2 extract_ChromHMM_state_summary

  Description: Extracting segmentation summaries
    Computes for a given ChromHMM, Segway or GMTK segmentation and a given state the
    number of celltypes with that state at each position
  Arg1: An $options hash ref with values:
    * $options->{trackhub_dir}
    * $options->{working_dir}
    * $options->{segmentations} List of hashes containing:
      . $segmentation->{name} String
      . $segmentation->{location} Path to directory
  Arg2: Segmentation name (string)
  Arg3: State name (string)
  Returntype: undef
  Side effects: Writes into:
    * Sorted copy of each cell type's segmentation:
      $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed
    * A summary statistic of each state:
      $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw

=cut

sub extract_ChromHMM_state_summary {
  my ($options, $segmentation, $state) = @_;
  my @celltypes = keys %{$segmentation->{celltypes}};
  my $cell;

  my @summaries = ();

  if (!must_compute($options,"$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw")) {
    print_log("Summary of state $state in segmentation $segmentation->{name} already computed, skipping...\n");
    return;
  } else {
    print_log("Computing summary of state $state in segmentation $segmentation->{name}.\n");
  }

  mkdir "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/";

  foreach $cell (@celltypes) {
    my $source = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed";
    my $summary = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$cell.bed";
    run("awk \'\$4 == \"$state\" \' $source > $summary") ;
    push @summaries, $summary;
  }

  my $summary = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.wig";
  run("wiggletools write_bg $summary sum " . join(" ", @summaries));
  convert_to_bigWig($options, $summary);
}

=head2 label_segmentation_states

  Description: Characterising segmentation states
    This first label assigns a label to each state of
    each segmentation, partly to color the segmentations
    partly to inform the next step.
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Side effects: It fills out:
    * $options->{assignments}->{$segmentation->{name}} Hash which assigns a label to each state
    It creates:
    * $options->{working_dir}/assignments.txt (Tab delimited file)

=cut

sub label_segmentation_states {
  my $options = shift;
  my ($segmentation, $state);

  if (!must_compute($options,"$options->{working_dir}/assignments.txt")) {
    print_log("Loading pre-existing assignments $options->{working_dir}/assignments.txt\n");
    $options->{assignments} = load_assignments( "$options->{working_dir}/assignments.txt");
  } else {
    print_log("Assigning states to functional labels\n");
    $options->{assignments} = {};
    foreach $segmentation (@{$options->{segmentations}}) {
      select_segmentation_states($options, $segmentation);
    }
    dump_assignments($options, "$options->{working_dir}/assignments.txt");
  }

  foreach $segmentation (@{$options->{segmentations}}) {
    my $assignments = $options->{assignments}->{$segmentation->{name}};
    foreach $state (keys %{$assignments}) {
      print_log("Assignment\t$segmentation->{name}\t$state\t$assignments->{$state}\n");
    }
  }
}

=head2 dump_assignments

  Description: Serialises state labelling
  Arg1: hash ref containing:
    - segmentations: array ref of segmentation names
    - assignments: hash ref mapping segmentation names to hash refs, which
        in turn map state names to labels
  Arg2: file location
  Returntype: undef
  Side effect: writes into file

=cut

sub dump_assignments {
  my ($options, $file) = @_;
  open my $fh, ">", $file;
  foreach my $segmentation (@{$options->{segmentations}}) {
    my $assignments = $options->{assignments}->{$segmentation->{name}};
    foreach my $state (@{$segmentation->{states}}) {
      print $fh "$segmentation->{name}\t$state\t$assignments->{$state}\n";
    }
  }
  close $fh;
}

=head2 load_assignments

  Description: Deserialises state labelling
  Arg1: file location
  Returntype: hash ref mapping segmentation names to hash refs, which
        in turn map state names to labels

=cut

sub load_assignments {
  my ($file) = @_;
  open my $fh, "<", $file;
  my $assignments = {};
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    if (!defined  $assignments->{$items[0]}) {
      $assignments->{$items[0]} = {};
    }
    $assignments->{$items[0]}->{$items[1]} = $items[2];
  }
  close $fh;
  return $assignments;
}

=head2 label_segmentation_states

  Description: Characterising segmentation states
    This first label assigns a label to each state of
    a given segmentation, partly to color the segmentation
    partly to inform the next step.
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation name
  Side effects: It fills out:
    * $options->{assignments}->{$segmentation->{name}} Hash which assigns a label to each state

=cut

sub select_segmentation_states {
  my ($options, $segmentation) = @_;
  my $state;

  $options->{assignments}->{$segmentation->{name}} = {};

  print_log("Entering segmentation $segmentation->{name}\n");

  compute_overlaps($options, $segmentation);

  foreach $state (@{$segmentation->{states}}) {
    label_segmentation_state($options, $segmentation, $state);
  }
}

=head2 label_segmentation_state

  Description: assigns a label to a state of
    a given segmentation, partly to color the segmentation
    partly to inform the next step.
  Arg1: A $options hash ref which contains:
    - assignments => hash ref:
      - segmentation name => empty hash ref
  Arg2: segmentation: hashref:
    - name => string
    - repressed_cutoff => numerical scalar
    - overlaps: hash ref:
      - ctcf => correlation (-1 <= x <= 1)
      - repressed => score
      - tfbs => enrichment score (> 0)
      - tss => enrichment score (> 0)
      - gene => enrichment score (> 0)
  Arg3: state name
  Side effects: It fills out:
    * $options->{assignments}->{$segmentation->{name}}->{state} with one of the different labels (See our @labels)

=cut

sub label_segmentation_state {
  my ($options, $segmentation, $state) = @_;
  my $overlaps = $segmentation->{overlaps};
  my $assignments = $options->{assignments}->{$segmentation->{name}};

  if ($overlaps->{ctcf}->{$state} > .25) {
    $assignments->{$state} = 'ctcf';
  } elsif ($overlaps->{repressed}->{$state} > $segmentation->{repressed_cutoff}) {
    if ($overlaps->{tfbs}->{$state} < $weak_cutoff) {
      $assignments->{$state} = 'repressed';
    } else {
      $assignments->{$state} = 'poised';
    }
  } elsif ($overlaps->{tfbs}->{$state} < 1 && $overlaps->{gene}->{$state} < 1) {
    $assignments->{$state} = 'dead';
  } elsif ($overlaps->{tfbs}->{$state} < $weak_cutoff && $overlaps->{gene}->{$state} < $weak_cutoff) {
    $assignments->{$state} = 'weak';
  } elsif ($overlaps->{tss}->{$state} > $very_strong_cutoff) {
    $assignments->{$state} = 'tss';
  } elsif ($overlaps->{gene}->{$state} > $overlaps->{tfbs}->{$state}) {
    $assignments->{$state} = 'gene';
  } elsif ($overlaps->{tss}->{$state} > $weak_cutoff) {
    $assignments->{$state} = 'proximal';
  } else {
    $assignments->{$state} = 'distal';
  }
}

=head2 label_segmentation_states

  Description: Characterising segmentation states
    with a number of scores
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
    - type => string
    - states => hash ref of strings
  Side effects: It fills out:
    * options->overlaps: hashref

=cut

sub compute_overlaps {
  my ($options, $segmentation) = @_;

  mkpath "$options->{working_dir}/overlaps/$segmentation->{name}";

  foreach my $test ('ctcf', 'tss', 'gene', 'tfbs') {
    compute_overlap_scores($options, $segmentation, $test);
  }
  if ($segmentation->{type} eq 'ChromHMM') {
    compute_ChromHMM_repressed_scores($options, $segmentation);
  } elsif ($segmentation->{type} eq 'Segway') {
    compute_Segway_repressed_scores($options, $segmentation);
  } elsif ($segmentation->{type} eq 'GMTK') {
    compute_GMTK_repressed_scores($options, $segmentation);
  } else {
    print STDERR "Could not recognize segmentation type $segmentation->{type}\n";
    exit 1;
  }

  open my $out, ">", "$options->{working_dir}/overlaps/$segmentation->{name}/summary.txt";
  print $out "Segmentation\tstate\tctcf\ttss\tgene\ttfbs\trepressed\n";
  foreach my $state (@{$segmentation->{states}}) {
    print $out "$segmentation->{name}\t$state";
    foreach my $test ('ctcf', 'tss', 'gene', 'tfbs','repressed') {
      chomp $segmentation->{overlaps}->{$test}->{$state};
      print $out "\t$segmentation->{overlaps}->{$test}->{$state}";
    }
    print $out "\n";
  }
  close $out;
}

=head2 compute_overlap_scores

  Description: Characterising segmentation states
    with a score
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
    - type => string
    - states => hash ref of strings
  Arg3: name of test (one of ctcf, tss, tfbs, gene or repressed)
  Side effects: It fills out:
    * options->overlaps: hashref

=cut

sub compute_overlap_scores {
  my ($options, $segmentation, $test) = @_;
  my $state;
  my $output = "$options->{working_dir}/overlaps/$segmentation->{name}/$test.txt";
  if (!must_compute($options,$output)) {
    print_log("Retrieving pre-computed overlaps for test $test\n");
    $segmentation->{overlaps}->{$test} = retrieve($output);
  } else {
    print_log("Computing overlaps for test $test\n");

    $segmentation->{overlaps}->{$test} = {};
    foreach $state (@{$segmentation->{states}}) {
      compute_overlap_score($options, $segmentation, $test, $state);
    }
    store $segmentation->{overlaps}->{$test}, $output;
  }
}

=head2 compute_overlap_scores

  Description: Characterising a segmentation state
    with a score
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
  Arg3: name of test (one of ctcf, tss, tfbs, gene or repressed)
  Arg3: name of state
  Side effects: It fills out:
    * segmentation->{overlaps}: hashref

=cut

sub compute_overlap_score {
  my ($options, $segmentation, $test, $state) = @_;

  my $reference;
  my $file = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";

  if ($test eq 'ctcf') {
    $reference = "$options->{trackhub_dir}/tfbs/CTCF.bw";
  } elsif ($test eq 'gene') {
    $reference = $options->{exons};
  } elsif ($test eq 'tfbs') {
    $reference = "$options->{trackhub_dir}/overview/all_tfbs.bw";
  } elsif ($test eq 'tss') {
    $reference = $options->{tss};
  } else {
    die "Unknown test\n";
  }

  if ($test eq 'ctcf') {
    print_log("Running wiggletools pearson $reference $file\n");
    $segmentation->{overlaps}->{$test}->{$state} = `wiggletools pearson $reference $file`;
  } else {
    my $genome_length = $options->{genome_length};
    $segmentation->{overlaps}->{$test}->{$state} = compute_enrichment_between_files($reference, $file, $genome_length);
  }

  chomp $segmentation->{overlaps}->{$test}->{$state};
}

=head2 compute_enrichment_between_files

  Description: compute enrichment between two bed files
  Arg1: reference file location
  Arg2: query file location
  Returntype: scalar

=cut

sub compute_enrichment_between_files {
  my ($reference, $file, $genome_length) = @_;

  if ($genome_length == 0) {
    use Carp;
    confess("genome_length parameter is zero! ($genome_length)");
  }
  
  my $auc = `wiggletools AUC mult $reference unit $file`;
  my $breadth = `wiggletools AUC unit $file`;
  my $ref_auc = `wiggletools AUC $reference`;

  if ($breadth == 0) {
    return 0;
  } else {
    return ($auc / $breadth) / ($ref_auc / $genome_length );
  }
}

=head2 compute_ChromHMM_repressed_scores

  Description: Characterising ChromHMM segmentation states
    with a repression score
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
  Side effects: It fills out:
    * segmentation->{overlaps}->{repressed}: hashref

=cut

sub compute_ChromHMM_repressed_scores {
  my ($options, $segmentation) = @_;
  my @files = glob "$segmentation->{location}/emissions*.txt";
  if (scalar @files != 1) {
    print STDERR "! Problem finding file $segmentation->{location}/*.txt";
    exit 1;
  }
  my $file = pop @files;
  my ($fh, $line);
  open $fh, "<", $file;
  $line = <$fh>;
  chomp $line;
  my @items = split /\t/, $line;
  my @columns = ();
  for (my $index = 1; $index < scalar @items; $index++) {
    if (grep($_ eq clean_name($items[$index]), @repressed_marks)) {
      push @columns, $index;
    }
  }

  my $max = 0;

  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    my $sum = 0;
    foreach my $column (@columns) {
       $sum += $items[$column];
    }
    $segmentation->{overlaps}->{repressed}->{"E$items[0]"} = $sum;
    print_log("Emission\t$segmentation->{name}\trepressed\t$items[0]\t$sum\n");
    if ($sum > $max) {
      $max = $sum;
    }
  }
  close $fh;

  $segmentation->{repressed_cutoff} = $max / 3;
}

=head2 compute_GMTK_repressed_scores

  Description: Characterising GMTK segmentation states
    with a repression score
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * GMTK segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
  Side effects: It fills out:
    * segmentation->{overlaps}->{repressed}: hashref
      - state name => scalar

=cut

sub compute_GMTK_repressed_scores {
  my ($options, $segmentation) = @_;
  my @files = glob "$segmentation->{location}/*.tab";
  if (scalar @files != 1) {
    print STDERR "! Problem finding file $segmentation->{location}/*.tab";
    exit 1;
  }
  my $file = pop @files;
  my ($fh, $line);
  open $fh, "<", $file;
  $line = <$fh>;
  chomp $line;
  my @columns = split /\t/, $line;

  my $max = 0;

  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    if (!grep($_ eq clean_name($items[0]), @repressed_marks)) {
      next;
    }
    for (my $index = 0; $index < scalar(@columns); $index++) {
      $segmentation->{overlaps}->{repressed}->{$columns[$index]} += $items[$index + 1];
      if ($segmentation->{overlaps}->{repressed}->{$columns[$index]} > $max) {
        $max = $segmentation->{overlaps}->{repressed}->{$columns[$index]};
      }
    }
  }
  close $fh;

  foreach my $state (@{$segmentation->{states}}) {
    print_log("Emission\t$segmentation->{name}\trepressed\t$state\t$segmentation->{overlaps}->{repressed}->{$state}\n");
  }
  $segmentation->{repressed_cutoff} = $max / 3;
}

=head2 compute_Segway_repressed_scores

  Description: Characterising Segway segmentation states
    with a repression score
  Arg1: A $options hash ref which contains:
    - $options->{working_dir}
    - $options->{trackhub_dir}
    - $options->{exons} Path to bed file with exonic regions
    - $options->{tss} Path to bed file with known TSSs
    - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
    - $options->{segmentations} where each value is a hash containing:
      . $segmentation->{name}
      . $segmentation->{states}
    It expects:
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
    * Segway segmentations: $segmentation->{location}/*.tab (1 file) Emission file
    * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
  Arg2: segmentation hashref
    - name => string
  Side effects: It fills out:
    * segmentation->{overlaps}->{repressed}: hashref

=cut

sub compute_Segway_repressed_scores {
  my ($options, $segmentation) = @_;
  my @files = glob "$segmentation->{location}/*.tab";
  if (scalar @files != 1) {
    print STDERR "! Problem finding file $segmentation->{location}/*.tab";
    exit 1;
  }
  my $file = pop @files;
  my ($fh, $line);
  open $fh, "<", $file;
  $line = <$fh>;

  my $max = 0;

  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    if (!grep($_ eq clean_name($items[1]), @repressed_marks)) {
      next;
    }
    if (!defined $segmentation->{overlaps}->{repressed}->{$items[0]}) {
      $segmentation->{overlaps}->{repressed}->{$items[0]} = 0;
    }
    $segmentation->{overlaps}->{repressed}->{$items[0]} += $items[2];
    if ($segmentation->{overlaps}->{repressed}->{$items[0]} > $max) {
      $max = $segmentation->{overlaps}->{repressed}->{$items[0]};
    }
  }
  close $fh;

  foreach my $state (@{$segmentation->{states}}) {
    print_log("Emission\t$segmentation->{name}\trepressed\t$state\t$segmentation->{overlaps}->{repressed}->{$state}\n");
  }

  $segmentation->{repressed_cutoff} = $max / 3;
}

=head2 make_segmentation_bedfiles

  Description: Preparing segmentation bedfiles
    This function adds in color information to the bedfiles,
    based on the assignment of states done above
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations} list which each element is a hash containing:
      . $segmentation->{name} String
      . $segmentation->{type} (ChromHMM or Segway)
      . $segmentation->{states} List of strings
      . $segmentation->{celltypes} List of strings
    It expects files:
    * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
  Returntype: undef
  Side effects: it creates:
  * $options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb

=cut

sub make_segmentation_bedfiles {
  my $options = shift;
  my $segmentation;
  mkdir "$options->{trackhub_dir}/segmentations/";
  foreach $segmentation (@{$options->{segmentations}}) {
    make_segmentation_bedfiles_2($options, $segmentation);
  }
}

=head2 make_segmentation_bedfiles_2

  Description: Preparing segmentation bedfiles for a given segmentation
    This function adds in color information to the bedfiles,
    based on the assignment of states done above
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations} list which each element is a hash containing:
      . $segmentation->{name} String
      . $segmentation->{type} (ChromHMM or Segway)
      . $segmentation->{states} List of strings
      . $segmentation->{celltypes} List of strings
    It expects files:
    * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
  Arg2: segmentation hash ref:
    name => string
    celltypes => array ref of strings
  Returntype: undef
  Side effects: it creates:
  * $options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb

=cut

sub make_segmentation_bedfiles_2 {
  my ($options, $segmentation) = @_;
  my $celltype;

  if (must_compute($options,"$options->{trackhub_dir}/segmentations/$segmentation->{name}")) {
    mkdir "$options->{trackhub_dir}/segmentations/$segmentation->{name}";
  }

  foreach $celltype (keys %{$segmentation->{celltypes}}) {
    make_segmentation_bedfile($options, $segmentation, $celltype);
  }
}

=head2 make_segmentation_bedfile

  Description: Preparing segmentation bedfiles for a given segmentation
    and a given cell type
    This function adds in color information to the bedfiles,
    based on the assignment of states done above
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations} list which each element is a hash containing:
      . $segmentation->{name} String
      . $segmentation->{type} (ChromHMM or Segway)
      . $segmentation->{states} List of strings
      . $segmentation->{celltypes} List of strings
    It expects files:
    * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
  Arg2: segmentation hash ref:
    name => string
    celltypes => array ref of strings
  Arg3: CellType name (string processed with clean_name)
  Returntype: undef
  Side effects: it creates:
  * $options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb

=cut

sub make_segmentation_bedfile {
  my ($options, $segmentation, $celltype) = @_;

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway' || $segmentation->{type} eq 'GMTK') {
    make_ChromHMM_bedfile($options, $segmentation, $celltype);
  } else {
    die("Unknown segmentation type $segmentation->{type}\n");
  }
}

=head2 make_ChromHMM_bedfile

  Description: Preparing ChromHMM, Segway or GMTK bedfiles for a given segmentation
    and a given cell type
    This function adds in color information to the bedfiles,
    based on the assignment of states done above
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations} list which each element is a hash containing:
      . $segmentation->{name} String
      . $segmentation->{type} (ChromHMM or Segway)
      . $segmentation->{states} List of strings
      . $segmentation->{celltypes} List of strings
    It expects files:
    * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
  Arg2: segmentation hash ref:
    name => string
    celltypes => array ref of strings
  Arg3: CellType name (string processed with clean_name)
  Returntype: undef
  Side effects: it creates:
  * $options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb

=cut

sub make_ChromHMM_bedfile {
  my ($options, $segmentation, $celltype) = @_;
  my $output = "$options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb";

  if (must_compute($options,$output)) {
    my @files = ();
    foreach my $state (@{$segmentation->{states}}) {
      push @files, make_ChromHMM_state_bedfile($options, $segmentation, $celltype, $state);
    }

    $output =~ s/.bb$/.bed/;
    run("sort -m ".join(" ", @files)." -k1,1 -k2,2n > $output");
    convert_to_bigBed($options, $output);
  } else {
    print_log("Colorised segmentation $output already exists, skipping...\n");
  }
}

=head2 make_ChromHMM_state_bedfile

  Description: Preparing ChromHMM, Segway or GMTK bedfiles for a given segmentation,
    a given cell type and a given state
    This function adds in color information to the bedfiles,
    based on the assignment of states done above
  Arg1: A $options hash which contains:
    - $options->{working_dir}
    - options->{assignments}: hashref:
      segmentation name => hashref:
        state name => label
    It expects files:
    * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
  Arg2: segmentation hash ref:
    name => string
  Arg3: CellType name (string processed with clean_name)
  Returntype: file location
  Side effects: it creates a file and writes into it

=cut

sub make_ChromHMM_state_bedfile {
  my ($options, $segmentation, $celltype, $state) = @_;
  my $input = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed";
  die ("Could not find $input\n") unless -e $input;
  my $output = "$options->{working_dir}/segmentations/$segmentation->{name}/$state/$celltype.bed";

  mkpath "$options->{working_dir}/segmentations/$segmentation->{name}/$state/";
  my $assigned_label = $options->{assignments}->{$segmentation->{name}}->{$state};
  my $color = $COLORS{$assigned_label};
  my $awk_string = "BEGIN {OFS=\"\\t\"} {\$4 = \"${state}_${assigned_label}_\"NR; \$5=1000; \$6=\".\"; \$7=\$2; \$8=\$3; \$9=\"$color\"; print}";
  run("awk \'$awk_string\' $input > $output") ;

  return $output;
}

=head2

  Description: Setting cutoffs
    This functions first selects which states are directly useful
    to infer TF binding, then determines an optimal
    cutoff so as to best fit the overall transcription factor binding
    probability:
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations} list which each element is a hash containing:
      . $segmentation->{name} String
      . $segmentation->{type} (ChromHMM or Segway)
      . $segmentation->{states} List of strings
      . $segmentation->{celltypes} List of strings
    It expects:
    * $options->{trackhub_dir}/overview/all_tfbs.bw
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
  Returntype: undef
  Side effects: It computes:
  * $options->{selected_states}->{$segmentation->{name}}->{$label}
    List of state names that given an appropriate cutoff can have a 5x enrichment
    in the reference marker
  * $options->{cutoffs}->{$segmentation->{name}}->{$label} Integer
    Cutoff on the sum of values of selected states that maximises
    the F-score of detection of TF binding
  It creates:
  * $options->{working_dir}/selected serialised hash (see above)
  * $options->{working_dir}/cutoffs serialised hash (see above)

=cut

sub set_cutoffs {
  my $options = shift;
  my $cutoffs_file = "$options->{working_dir}/cutoffs";
  my $selected_file = "$options->{working_dir}/selected";

  if (must_compute($options,$cutoffs_file) || must_compute($options,$selected_file)) {
    $options->{cutoffs} = {};
    $options->{selected_states} = {};
    foreach my $segmentation (@{$options->{segmentations}}) {
      select_segmentation_cutoffs($options, $segmentation);
    }
    store $options->{cutoffs}, $cutoffs_file;
    store $options->{selected_states}, $selected_file;
  } else {
    $options->{cutoffs} = retrieve $cutoffs_file;
    $options->{selected_states} = retrieve $selected_file;
  }
}

=head print_cutoffs

  Desription: print cutoffs into logs
  Returntype: undef
  Side effects: prints into stdout

=cut

sub print_cutoffs {
  my $options = shift;
  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach my $label (keys %{$options->{selected_states}->{$segmentation->{name}}}) {
      print_log("Selected\t$segmentation->{name}\t$label\t". join(" ", keys %{$options->{selected_states}->{$segmentation->{name}}->{$label}})."\n");
    }
  }

  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach my $label (keys %{$options->{cutoffs}->{$segmentation->{name}}}) {
      print_log("Cutoff\t$segmentation->{name}\t$label\t$options->{cutoffs}->{$segmentation->{name}}->{$label}\n");
    }
  }
}

=head2 select_segmentation_cutoffs

  Description: Setting cutoffs
    This functions first selects which states of a given segmentation are directly useful
    to infer TF binding, then determines an optimal
    cutoff so as to best fit the overall transcription factor binding
    probability:
  Arg1: A $options hash which contains:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    It expects:
    * $options->{trackhub_dir}/overview/all_tfbs.bw
    * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
  Arg2: Segmentation, hash ref:
    . $segmentation->{name} String
    . $segmentation->{type} (ChromHMM or Segway)
    . $segmentation->{states} List of strings
    . $segmentation->{celltypes} List of strings
  Returntype: undef
  Side effects: It computes:
  * $options->{selected_states}->{$segmentation->{name}}->{$label}
    List of state names that given an appropriate cutoff can have a 5x enrichment
    in the reference marker
  * $options->{cutoffs}->{$segmentation->{name}}->{$label} Integer
    Cutoff on the sum of values of selected states that maximises
    the F-score of detection of TF binding
  It creates:
  * $options->{working_dir}/selected serialised hash (see above)
  * $options->{working_dir}/cutoffs serialised hash (see above)

=cut

sub select_segmentation_cutoffs {
  my ($options, $segmentation) = @_;
  $options->{selected_states}->{$segmentation->{name}} = {};
  $options->{cutoffs}->{$segmentation->{name}} = {};

  foreach my $label (@functional_labels) {
    $options->{selected_states}->{$segmentation->{name}}->{$label} = select_relevant_states($options, $segmentation, $label);
    $options->{cutoffs}->{$segmentation->{name}}->{$label} = select_segmentation_cutoff($options, $segmentation, $label);
  }
}

=head2 select_relevant_states

  Description: for a given segmentations, returns weighted list of relevant states with a given label
  Arg1: options: hashref
  Arg2: Segmentation, hash ref:
    . $segmentation->{name} String
    . $segmentation->{type} (ChromHMM or Segway)
    . $segmentation->{states} List of strings
    . $segmentation->{celltypes} List of strings
  Arg3: label, string contained in @labels
  Returntype: hashref:
    - state => weight scalar

=cut

sub select_relevant_states {
  my ($options, $segmentation, $label) = @_;
  my %states = ();
  my $assignments = $options->{assignments}->{$segmentation->{name}};
  foreach my $state (@{$segmentation->{states}}) {
    if ($assignments->{$state} eq $label) {
      print_log("Selecting $state $label\n");
      my $weight = test_relevance($options, $segmentation, $label, $state);
      if ($weight > 0) {
        $states{$state} = $weight;
      }
    }
  }
  print_log("Selected\t$segmentation->{name}\t$label\t". join(" ", keys %states)."\n");
  return \%states;
}

=head2 select_relevant_states

  Description: for a given segmentation, a given state and a given label, returns confidence in that state's labelling
  Arg1: options: hashref
  Arg2: Segmentation, hash ref:
    . $segmentation->{name} String
    . $segmentation->{type} (ChromHMM or Segway)
    . $segmentation->{states} List of strings
    . $segmentation->{celltypes} List of strings
  Arg3: label, string contained in @labels
  Arg4: state, string
  Returntype: weight scalar

=cut

sub test_relevance {
  my ($options, $segmentation, $label, $state) = @_;
  my $reference = "$options->{trackhub_dir}/overview/all_tfbs.bw";
  my $file = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";

  my $genome_length = $options->{genome_length};
  
  for (my $i = 1; $i < scalar(keys %{$segmentation->{celltypes}}); $i++) {
    my $enrichment = compute_enrichment_between_files($reference, "gt $i $file", $genome_length);
    print_log("Enrichment\t$state\t$i:\t$enrichment\n");
    if ($enrichment == 0) {
      return 0;
    } elsif (compute_enrichment_between_files($reference, "gt $i $file", $genome_length) > $weak_cutoff) {
      return $i;
    }
  }

  return 0;
}

=head2 select_segmentation_cutoff

  Description: for a given segmentation, and a given label, determine cutoff to call a new feature
  Arg1: options: hashref
  Arg2: Segmentation, hash ref:
    . $segmentation->{name} String
    . $segmentation->{type} (ChromHMM or Segway)
    . $segmentation->{states} List of strings
    . $segmentation->{celltypes} List of strings
  Arg3: label, string contained in @labels
  Returntype: cutoff scalar

=cut

sub select_segmentation_cutoff {
  my ($options, $segmentation, $label) = @_;
  print_log("Setting cutoff for $label...\n");

  my $tfbs = "$options->{trackhub_dir}/overview/all_tfbs.bw";
  my $tfbs_auc = `wiggletools AUC $tfbs`;

  my $celltype_count = scalar(keys %{$segmentation->{celltypes}});
  my $max_weight = undef;
  my $min_step = 1;
  my @files = ();
  my @weights = ();
  my %hash = %{$options->{selected_states}->{$segmentation->{name}}->{$label}};
  foreach my $state (keys %hash) {
    my $cutoff = $hash{$state};
    push @files, "scale ".(1/$cutoff)." $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";
    $min_step /= $cutoff;
    if ((! defined $max_weight) || $max_weight < 1 / $cutoff) {
      $max_weight = 1 / $cutoff;
    }
  }

  if (scalar @files == 0) {
    return 0;
  }

  my $files_string = join(" ", @files);

  my $last_fscore = 0;
  my $last_cutoff = 0;
  # Searching for cutoff which produces max fscore
  for (my $i = 0; $i < $celltype_count * $max_weight; $i += $min_step) {
    my $score = fscore($tfbs, $tfbs_auc, $files_string, $i);
    print_log("CutoffTest\t$label\t.\t$score\n");

    # If current value smaller than previous, just got past peak (assuming convexity of fscore(i) )
    if ($score < $last_fscore) {
      print_log("cutoff set at $last_cutoff\n");
      return $last_cutoff;
    } else {
      $last_fscore = $score;
      $last_cutoff = $i;
    }
  }

  print_log("cutoff set at MAX\n");
  return $celltype_count * $max_weight - $min_step;
}

=head2 fscore

  Description: compute fscore of a composite wiggletools function with cutoff wrt general TF binding
  Arg1: location of TFBS summary
  Arg2: wiggletools expression
  Arg3: coverage cutoff
  Returntype: scalar

=cut

sub fscore {
  my ($tfbs, $tfbs_auc, $files_string, $i) = @_;
  # Compute fscore for that cutoff
  my $auc = `wiggletools AUC mult $tfbs gt $i sum $files_string`;
  my $sens = $auc / $tfbs_auc;

  my $overlap = `wiggletools AUC gt 0 mult $tfbs gt $i sum $files_string`;
  my $breadth = `wiggletools AUC gt $i sum $files_string`;
  my $spec;
  if ($breadth == 0) {
    $spec = 0;
  } else {
    $spec = $overlap / $breadth;
  }

  if ($sens + $spec == 0) {
    return 0;
  } else {
    return 2 * ($spec * $sens) / ($spec + $sens);
  }
}

=head2 compute_regulatory_features

  Description: This is the secret sauce:
  Foreach function ('tss', 'proximal', 'distal', 'ctcf'):
    Foreach segmentation:
      - Sum the summaries of all the states associated to
      that function
      - Select the regions where the sum is greater than
      the cutoff (determined above)
    Compute the union of these regions

  We then apply heuristic rules:
  * TSSs which have no experimental validation (CAGE)
  are downgraded to proximal
  * Proximal enhancers which overlap with a TSS are added
  to the TSS's whiskers (and do not exist separately)
  * TFBS+DNase regions which do not overlap any chromatin
  region are labelled 'tfbs'
  * DNAse regions which do not overlap any of the above
  are labelled 'DNAse'

  Arg1: $options: hash ref
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentations}: List of hashes with:
      . $segmentation->{name}
    - $options->{selected_states}->{$segmentation->{name}}->{$label}: List of strings
    - $options->{cutoffs}->{$segmentation->{name}}->{$label}: Scalar
    - $options->{tss} Path to BEd file with experimentally validated TSS
  Expected files:
  * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/overview/RegBuild.bb

=cut

sub compute_regulatory_features {
  my $options = shift;
  my $bed_output = "$options->{trackhub_dir}/overview/RegBuild.bed";
  my $bb_output = "$options->{trackhub_dir}/overview/RegBuild.bb";

  mkdir "$options->{working_dir}/build/";

  if (!must_compute($options,$bb_output)) {
    print_log("$bb_output already exists, skipping...\n");
    return;
  } else {
    print_log("Computing Regulatory Build...\n");
  }

  #############################################
  ## Precompute everything independently
  #############################################

  my $tss_tmp = compute_initial_regions($options, "tss");
  my $proximal_tmp = compute_initial_regions($options, "proximal");
  my $distal_tmp = compute_initial_regions($options, "distal");
  my $ctcf_tmp = compute_initial_regions($options, "ctcf");

  # Compute TF binding
  my $tfbs_signal = "$options->{trackhub_dir}/overview/all_tfbs.bw";
  my $tfbs_tmp = "$options->{working_dir}/build/tfbs.tmp.bed";
  my $awk_tfbs = make_awk_command("tfbs");
  run("wiggletools write_bg - unit $tfbs_signal | $awk_tfbs > $tfbs_tmp");

  #Compute open dnase
  my $dnase_tmp = "$options->{working_dir}/build/dnase.tmp.bed";
  my @dnase_files = glob "$options->{working_dir}/celltype_dnase/*.bed";
  my $awk_dnase = make_awk_command("dnase");
  run("wiggletools write_bg - unit sum ".join(" ", @dnase_files)." | $awk_dnase > $dnase_tmp");

  #############################################
  ## Percolate overlapping features up each level
  #############################################

  # DNAse sites are merged to overlapping TFBS sites
  my $tfbs_tmp2 = undef;
  my $remove_tfbs = "";
  if (defined $tfbs_tmp) {
    $tfbs_tmp2 = "$options->{working_dir}/build/tfbs.tmp2.bed";
    expand_boundaries([$dnase_tmp], $tfbs_tmp, $tfbs_tmp2);
    $remove_tfbs = " | bedtools intersect -wa -v -a stdin -b $tfbs_tmp2";
  }

  # DNAse and TFBS sites are merged to overlapping distal sites
  my $distal_tmp2 = undef;
  my $remove_distal = "";
  if (defined $distal_tmp) {
    $distal_tmp2 = "$options->{working_dir}/build/distal.tmp2.bed";
    expand_boundaries([$dnase_tmp, $tfbs_tmp], $distal_tmp, $distal_tmp2);
    $remove_distal = " | bedtools intersect -wa -v -a stdin -b $distal_tmp2";
  }

  # DNAse, TFBS and distal sites are merged to overlapping proximal sites
  my $proximal_tmp2 = undef;
  my $remove_proximal = "";
  if (defined $proximal_tmp) {
    $proximal_tmp2 = "$options->{working_dir}/build/proximal.tmp2.bed";
    expand_boundaries([$dnase_tmp, $tfbs_tmp, $distal_tmp], $proximal_tmp, $proximal_tmp2);
    $remove_proximal = " | bedtools intersect -wa -v -a stdin -b $proximal_tmp2";
  }

  # DNAse, TFBS, distal and proximal sites are merged to overlapping TSS sites
  my $tss_tmp2 = undef;
  my $remove_tss = "";
  if (defined $tss_tmp) {
    $tss_tmp2 = "$options->{working_dir}/build/tss.tmp2.bed";
    expand_boundaries([$dnase_tmp, $tfbs_tmp, $distal_tmp, $proximal_tmp], $tss_tmp, $tss_tmp2);
    $remove_tss = " | bedtools intersect -wa -v -a stdin -b $tss_tmp2";
  }

  #############################################
  ## Apply mask
  #############################################

  my $awk_mask = "";
  if (defined $options->{mask} && -s $options->{mask} > 0) {
    $awk_mask = "| bedtools intersect -wa -v -a stdin -b $options->{mask}";
  }
  my $awk_contract = "awk 'BEGIN {OFS=\"\t\"} \$3 > \$2 + 1 {\$2 += 1; \$3 -= 1;} \$7 < \$2 {\$7 = \$2} \$8 > \$3 {\$8 = \$3} {print}'";
  my $final_filter = "$awk_mask | $awk_contract";

  #############################################
  ## Find non overlapping features for each level
  #############################################

  # All features that overlap known TSS are retained
  my $tss = undef;
  my $demoted = undef;
  if (defined $tss_tmp2) {
    $tss = "$options->{working_dir}/build/tss.bed";
    run("bedtools intersect -u -wa -a $tss_tmp2 -b $options->{tss} $final_filter | sort -k1,1 -k2,2n > $tss");

    # All features that do not go into a demoted file
    $demoted = "$options->{working_dir}/build/demoted_tss.bed";
    run("bedtools intersect -v -wa -a $tss_tmp2 -b $options->{tss} $final_filter | sort -k1,1 -k2,2n > $demoted");

    $remove_tss = " | bedtools intersect -wa -v -a stdin -b $tss";
  }

  # Unaligned proximal sites are retained
  my $proximal = undef;
  if (defined $proximal_tmp2 or defined $demoted) {
    $proximal = "$options->{working_dir}/build/proximal.bed";
    my $awk_proximal = make_awk_command("proximal", 0);
    my $files = join(" ", grep defined, ($proximal_tmp2, $demoted));
    run("wiggletools write_bg - unit sum $files | $awk_proximal $remove_tss $final_filter > $proximal");
  }

  # Unaligned distal sites are retained
  my $distal = undef;
  if (defined $distal_tmp2) {
    $distal = "$options->{working_dir}/build/distal.bed";
    run("cat $distal_tmp2 $remove_tss $remove_proximal $final_filter > $distal");
  }

  # Unaligned TFBS sites are retained
  my $tfbs = "$options->{working_dir}/build/tfbs.bed";
  run("cat $tfbs_tmp2 $remove_tss $remove_proximal $remove_distal $final_filter > $tfbs");

  # Unaligned DNAse sites are retained
  my $dnase = "$options->{working_dir}/build/dnase.bed";
  run("bedtools intersect -wa -v -a $dnase_tmp -b $tfbs_tmp2 $remove_tss $remove_proximal $remove_distal $final_filter > $dnase");

  #############################################
  ## CTCF computed independently
  #############################################

  my $ctcf = undef;
  if (defined $ctcf_tmp) {
    $ctcf = "$options->{working_dir}/build/ctcf.bed";
    run("cat $ctcf_tmp $final_filter > $ctcf");
  }

  #############################################
  ## Merge
  #############################################
  my $files = join(" ", grep defined, ($tss, $proximal, $distal, $ctcf, $dnase, $tfbs));
  run("sort -m $files -k1,1 -k2,2n > $bed_output");
  convert_to_bigBed($options, $bed_output);
}

=header2 expand_boundaries

  Description: expands the boundaries of a bed file
    to encompass any overlapping regions contained in an array of bed files
    Writes the results into a final output file
  Arg1: source_files: array ref of file locations containing extensions
  Arg2: target_file: bed file location containing seed regions
  Arg3: output: location of destination file
  Returntype: undef
  Side effects: writes into output

=cut

sub expand_boundaries {
  my ($source_files, $target_file, $output) = @_;
  my $awk_move_boundaries = "awk 'BEGIN {OFS=\"\\t\"} \$4 != name {if (name) {print chr, start, end, name, 1000, \".\", thickStart, thickEnd, rgb; } chr=\$1; start=\$2; end=\$3; name=\$4; thickStart=\$7; thickEnd=\$8; rgb=\$9} \$10 == chr && \$11+1 < start {start=\$11+1} \$10 == chr && \$12-1 > end {end=\$12-1} END {if (chr) {print chr, start, end, name, 1000, \".\", thickStart, thickEnd, rgb}}'";

  my @defined_source_files = grep defined,  @{$source_files};
  if (scalar @defined_source_files) {
    run("cat ".join(" ", @{$source_files}). " | bedtools intersect -wa -wb -a $target_file -b stdin | $awk_move_boundaries > $output");
    run("cat ".join(" ", @{$source_files}). " | bedtools intersect -v -wa -a $target_file -b stdin >> $output");
    run("sort -k1,1 -k2,2n $output > $output.sorted");
    run("mv $output.sorted $output");
  } else {
    # No extension, just copy
    run("cp $target_file $output");
  }
}

=header2

  Description: computes regions based on regions where a weighted
    sum of segmentation states are above the required threshold
  Arg1: options hash ref
    - working_dir => string
    - selected_states => hash ref:
      - segmentation name => hash ref:
        - label string => hash ref
          - state string => scalar
    - cutoffs => hash ref:
      - segmentation name => hash ref:
        - label string => scalar
  Arg2: label string
  Returntype: file location

=cut

sub compute_initial_regions {
  my ($options, $label) = @_;
  my $output = "$options->{working_dir}/build/$label.tmp.bed";
  my $wiggletools_cmd = "wiggletools write_bg - unit sum ";
  my $wiggletools_params = "";
  foreach my $segmentation (@{$options->{segmentations}}) {
    $wiggletools_params .= weighted_summary_definition($options, $label, $segmentation);
  }
  if (length $wiggletools_params == 0) {
    return undef;
  } else {
    my $awk_cmd = make_awk_command($label);
    run("$wiggletools_cmd $wiggletools_params | $awk_cmd > $output");
    return $output;
  }
}

=header2

  Description: Creates a wiggletools command which computes the weighted sum of segmentation
    states which share a given label
  Arg1: options hash ref:
    - selected_states => hash ref:
      - segmentation name => hash ref:
        - label string => hash ref
          - state string => scalar
    - cutoffs => hash ref:
      - segmentation name => hash ref:
        - label string => scalar
  Arg2: label string
  Arg3: segmentation hash ref:
    - name => string
  Returntype: string

=cut

sub weighted_summary_definition {
  my ($options, $label, $segmentation) = @_;
  my $hash = $options->{selected_states}->{$segmentation->{name}}->{$label};
  my $cutoff = $options->{cutoffs}->{$segmentation->{name}}->{$label};
  my @summaries = ();
  foreach my $state (keys %$hash) {
    push @summaries, "scale ".(1/$hash->{$state})." $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";
  }
  if (scalar @summaries > 0) {
    return " gt $cutoff sum " . join(" ", @summaries) . " :";
  } else {
    return "";
  }
}

=header2 make_awk_command

  Description: Creates an awk command that takes in a BedGraph and converts into a Bed9
  Arg1: label assigned to the regions of the bedgraph
  Arg2: boolean to decide whether to expand the regions by 1 bp in each direction to compute overlaps
  Returntype: string

=cut

sub make_awk_command {
  my ($label, $expand) = @_;
  my $expansion = "if (\$2 > 0) \$2 -= 1; \$3 += 1;";
  if (defined $expand && $expand == 0) {
    $expansion = "";
  }
  # Note that it IDs the lines with the line number
  return "awk \'BEGIN {OFS=\"\\t\"} {$expansion \$4=\"${label}_\"NR; \$5=1000; \$6=\".\"; \$7=\$2; \$8=\$3; \$9=\"$COLORS{$label}\"; print }\'";
}

=header2 compute_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know which
    Feature is active in which cell type
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{trackhub_dir}/projected_segmentations/$celltype.bed

=cut

sub compute_states {
  my $options = shift;
  my $segmentation;

  mkdir "$options->{trackhub_dir}/projected_segmentations/";
  mkdir "$options->{working_dir}/projected_segmentations/";

  foreach $segmentation (@{$options->{segmentations}}) {
    compute_segmentation_states($options, $segmentation);
  }
}

=header2 compute_segmentation_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know which
    Feature is active in which cell type in a given segmentation
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Arg2: segmentation
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{trackhub_dir}/projected_segmentations/$celltype.bed

=cut

sub compute_segmentation_states {
  my ($options, $segmentation) = @_;
  my $celltype;

  mkdir "$options->{working_dir}/projected_segmentations/$segmentation->{name}/";

  foreach $celltype (keys %{$segmentation->{celltypes}}) {
    if (must_compute($options,"$options->{trackhub_dir}/projected_segmentations/$celltype.bb")) {
      compute_celltype_state($options, $segmentation, $celltype);
    } else {
      print_log("Projected segmentation of $celltype already computed, skipping...\n");
    }
  }
}

=header2 compute_celltype_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know which
    Feature is active in a given cell type
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Arg2: segmentation hash ref
  Arg3: Bio::EnsEMBL::FuncgenCellType name
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{trackhub_dir}/projected_segmentations/$celltype.bed

=cut

sub compute_celltype_state {
  my ($options, $segmentation, $celltype) = @_;

  mkdir "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype";

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway' || $segmentation->{type} eq 'GMTK') {
    compute_ChromHMM_celltype_state($options, $segmentation, $celltype);
  } else {
    die ("Unknown segmentation type $segmentation->{type}");
  }
}

=header2 compute_ChromHMM_celltype_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know which
    Feature is active in a given cell type in a ChromHMM segmentation
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Arg2: segmentation hash ref
  Arg3: Bio::EnsEMBL::FuncgenCellType name
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{trackhub_dir}/projected_segmentations/$celltype.bed

=cut

sub compute_ChromHMM_celltype_state {
  my ($options, $segmentation, $celltype) = @_;
  my $output = "$options->{trackhub_dir}/projected_segmentations/$celltype.bed";

  foreach my $prelabel ((@functional_labels, 'repressed', 'poised')) {
    precompute_ChromHMM_label_state($options, $segmentation, $celltype, $prelabel);
  }

  my @bedfiles = ();
  foreach my $label (@labels) {
    push @bedfiles, compute_ChromHMM_label_state($options, $segmentation, $celltype, $label);
  }

  run("sort -m " . join(" ", @bedfiles) . " -k1,1 -k2,2n > $output") ;
  convert_to_bigBed($options, $output);
}

=header2 precompute_ChromHMM_label_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know where evidence
    of a given label can be found in a given cell
    type in a ChromHMM segmentation
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Arg2: segmentation hash ref
  Arg3: Bio::EnsEMBL::FuncgenCellType name
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.bed$

=cut

sub precompute_ChromHMM_label_state {
  my ($options, $segmentation, $celltype, $label) = @_;
  my $output = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.bed";
  my @files = ();
  foreach my $state (@{$segmentation->{states}}) {
    if ($options->{assignments}->{$segmentation->{name}}->{$state} eq $label) {
      push @files, "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed";
    }
  }

  if (scalar @files > 0) {
    run("sort -m " . join(" ", @files). " -k1,1 -k2,2n > $output");
  } elsif ($label eq "repressed" || $label eq "poised") {
    run("echo > $output");
  }
}

=header2 compute_ChromHMM_celltype_states

  Description: Computing cell-type specific state
    Given the Regulatory Build, you want to know which
    Feature of a given label is active in a given cell type in a ChromHMM segmentation
  Arg1: $options: hash containing:
    - $options->{trackhub_dir}
    - $options->{working_dir}
    - $options->{segmentation}
  Arg2: segmentation hash ref
  Arg3: Bio::EnsEMBL::FuncgenCellType name
  Arg4: label name
  Expected files:
  * $options->{working_dir}/build/$label.bed
  * $options->{working_dir}/celltype_tf/$celltype.bed
  * $options->{working_dir}/celltype_dnase/$celltype.bed
  Side effects: Created files:
  * $options->{trackhub_dir}/projected_segmentations/$celltype.bed
  * $options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.bed$

=cut

sub compute_ChromHMM_label_state {
  my ($options, $segmentation, $celltype, $label) = @_;
  my $output = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.final.bed";
  my $reference = "$options->{working_dir}/build/$label.bed";

  ## Hack to resolve issue https://www.ebi.ac.uk/panda/jira/browse/ENSREGULATION-417
  (my $epigenome_name = $celltype) =~ s/_25_SEGMENTS//;
  my $clean_epigenome_name = clean_name($epigenome_name);

  my $temp;
  if ($label eq 'tfbs') {
    #$temp = "$options->{working_dir}/celltype_tf/$celltype.bed";
    $temp = "$options->{working_dir}/celltype_tf/$clean_epigenome_name.bed";
  } elsif ($label eq 'dnase') {
    #$temp = "$options->{working_dir}/celltype_dnase/$celltype.bed";
    $temp = "$options->{working_dir}/celltype_dnase/$clean_epigenome_name.bed";
  } else {
    $temp = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.bed";
  }

  if (-e $reference && -s $reference && defined $temp && -e $temp && -s $temp) {
    my $repressed = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/repressed.bed";
    my $poised = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/poised.bed";

    # POISED => POISED
    my $temp2 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.poised.bed";
    my $poised_color = $COLORS{poised};
    run("bedtools intersect -wa -u -a $reference -b $poised | awk \'{\$9=\"$poised_color\"; print}\' > $temp2");

    # !POISED && ACTIVE && REPRESSED => POISED
    my $temp3 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.active_repressed.bed";
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -u -a stdin -b $temp | bedtools intersect -wa -u -a stdin -b $repressed | awk \'{\$9=\"$poised_color\"; print}\' > $temp3");

    # !POISED && !ACTIVE && REPRESSED => REPRESSED
    my $temp4 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.repressed.bed";
    my $repressed_color = $COLORS{repressed};
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -v -a stdin -b $temp | bedtools intersect -wa -u -a stdin -b $repressed | awk \'{\$9=\"$repressed_color\"; print}\' > $temp4");

    # !POISED && ACTIVE && !REPRESSED => ACTIVE
    my $temp5 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.active.bed";
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -u -a stdin -b $temp | bedtools intersect -wa -v -a stdin -b $repressed > $temp5");

    # !POISED && !ACTIVE && !REPRESSED => DEAD
    my $temp6 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$label.inactive.bed";
    my $dead_color = $COLORS{dead};
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -v -a stdin -b $temp | bedtools intersect -wa -v -a stdin -b $repressed | awk \'{\$9=\"$dead_color\"; print}\' > $temp6");

    # Merge all this into one bed file
    run("sort -m $temp2 $temp3 $temp4 $temp5 $temp6 -k1,1 -k2,2n > $output");
    return $output;
  } else {
    # Could not find evidence for or against, do not report any regions
    my $na_color = $COLORS{na};
    run("cat $reference | awk \'{\$9=\"$na_color\"; print}\' > $output");
    return $output;
  }
}

=head2

  Description: Track hub creation
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Returntype: undef
  Side effects: Writes into:
  * hub.txt
  * genomes.txt
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub {
  my $options = shift;
  print_log("Making track hub header files\n");
  make_track_hub_headers($options);
  make_track_hub_assembly($options);
}

=head2

  Description: Track hub header creation
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Returntype: undef
  Side effects: Writes into:
  * hub.txt
  * genomes.txt

=cut

sub make_track_hub_headers {
  my $options = shift;
  my $fh;

  open($fh, ">", "$options->{output_dir}/genomes.txt");
  print $fh "genome $options->{assembly}\n";
  print $fh "trackDb trackDb_$options->{assembly}.txt\n";
  close($fh);

  open($fh, ">", "$options->{output_dir}/hub.txt");
  print $fh "hub EnsemblRegulatoryBuild\n";
  print $fh "shortLabel Beta Regulatory Build\n";
  print $fh "longLabel Evidence summaries and provisional results for the new Ensembl Regulatory Build\n";
  print $fh "genomesFile genomes.txt\n";
  print $fh "email helpdesk\@ensembl.org\n";
  close($fh);
}

=head2 make_track_hub_assembly

  Description: Track hub creation
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_assembly {
  my $options = shift;

  my $output = "$options->{output_dir}/trackDb_$options->{assembly}.txt";
  open(my $file, ">", $output) || die("Could not open $output\n");

  make_track_hub_overview($options, $file);
  make_track_hub_segmentations($options, $file);
  make_track_hub_segmentation_summaries($options, $file);
  make_track_hub_projected_segmentations($options, $file);
  make_track_hub_tfbs($options, $file);

  close $file;
}

=head2 make_track_hub_overview

  Description: Track hub Overview composite track
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Arg2: filehandle
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_overview {
  my ($options, $file) = @_;

  print $file "track BuildOverview\n";
  print $file "compositeTrack on\n";
  print $file "visibility squish\n";
  print $file "shortLabel Build Overview\n";
  print $file "longLabel Ensembl Regulatory Build overviews\n";
  print $file "subGroup1 source Source dimensions dimY=source Regulatory_Build=EnsRegBuild Transcription_Factor_Binding_Site_Peaks_Summary=TFBS\n";
  print $file "sortOrder source=+ \n";
  print $file "dragAndDrop subTracks\n";
  print $file "priority 1\n";
  print $file "type bigBed 9\n";
  print $file "noInherit on\n";
  print $file "html  overview/index\n";

  print $file "\n";

  print $file "\ttrack RegBuildOverview\n";
  print $file "\tparent BuildOverview on\n";
  print $file "\tbigDataUrl $options->{assembly}/overview/RegBuild.bb\n";
  print $file "\tshortLabel Regulatory Build\n";
  print $file "\tlongLabel Ensembl Regulatory label of regional function\n";
  print $file "\ttype bigBed 9\n";
  print $file "\titemRgb on\n";
  print $file "\tsubGroups source=EnsRegBuild\n";
  print $file "\tpriority 1\n";
  print $file "\tvisibility squish\n";

  print $file "\n";

  print $file "\ttrack TFBSSummary\n";
  print $file "\tparent BuildOverview on\n";
  print $file "\tbigDataUrl $options->{assembly}/overview/all_tfbs.bw\n";
  print $file "\tshortLabel TFBS Summary\n";
  print $file "\tlongLabel Summary of Transcription Factor Binding Site peaks from all cell types\n";
  print $file "\ttype bigWig\n";
  print $file "\tautoscale off\n";
  print $file "\tmaxHeightPixels 128:64:16\n";
  print $file "\tviewLimits 0:1\n";
  print $file "\tcolor 209,157,0\n";
  print $file "\tsubGroups source=TFBS\n";
  print $file "\tpriority 3\n";
  print $file "\tvisibility full\n";
}

=head2 make_track_hub_segmentations

  Description: Track hub Segmentations composite track
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Arg2: filehandle
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_segmentations {
  my ($options, $file) = @_;
  my $segmentation;
  foreach $segmentation (@{$options->{segmentations}}) {
    make_track_hub_segmentations_2($options, $file, $segmentation);
  }
}

=head2 make_track_hub_segmentations_2

  Description: Track hub Segmentations composite track for a
    given segmentation
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Arg2: filehandle
  Arg3: segmentation hashref
    - name => string
    - celltypes => arrayref of strings
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_segmentations_2 {
  my ($options, $file, $segmentation) = @_;

  my $name = $segmentation->{name};
  my @states = $segmentation->{states};
  my $state_count = scalar(@states);
  my @celltypes = keys %{$segmentation->{celltypes}};
  my $celltype_count = scalar(@celltypes);
  my $state_string = "";

  my $celltype_list = "";
  my $celltype;
  foreach $celltype (@celltypes) {
    my $short = $celltype;
    $short =~ s/[_-]//g;
    $celltype_list .= " $short=$celltype";
  }

  print $file "\n";

  print $file "track ${name}CellSegments\n";
  print $file "compositeTrack on\n";
  print $file "visibility dense\n";
  print $file "shortLabel ${name} Segmentations\n";
  print $file "longLabel Cell-type specific $name genome segmentations\n";
  print $file "subGroup1 cellType Cell_Type $celltype_list\n";
  print $file "subGroup2 source Source ERB=Ensembl_Regulatory_Build\n";
  print $file "dimensions dimY=cellType\n";
  print $file "sortOrder cellType=+\n";
  print $file "dragAndDrop subTracks\n";
  print $file "priority 5\n";
  print $file "type bigBed 9\n";
  print $file "noInherit on\n";

  foreach $celltype (@celltypes) {
    my $short = $celltype;
    $short =~ s/[_-]//g;

    print $file "\n";

    print $file "\ttrack ${name}_${celltype}_segments\n";
    print $file "\tbigDataUrl $options->{assembly}/segmentations/$segmentation->{name}/$celltype.bb\n";
    print $file "\tparent ${name}CellSegments off\n";
    print $file "\tshortLabel Seg. ${celltype}\n";
    print $file "\tlongLabel $segmentation->{type} segmentation of the $celltype genome after training on $celltype_count cell types\n";
    print $file "\ttype bigBed 9\n";
    print $file "\titemRgb on\n";
    print $file "\tsubGroups cellType=$short source=ERB\n";
  }
}

=head2 make_track_hub_segmentation_summaries

  Description: Track hub Segmentation_Summaries composite track
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Arg2: filehandle
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_segmentation_summaries {
  my ($options, $file) = @_;
  my $segmentation;

  foreach $segmentation (@{$options->{segmentations}}) {
    make_track_hub_segmentation_summaries_2($options, $file, $segmentation);
  }
}

=head2 make_track_hub_segmentation_summaries_2

  Description: Track hub Segmentation_Summaries composite track
    for a given segmentation
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
    Expects files:
    * $options->{trackhub_dir}/tfbs/$tf.bw
  Arg2: filehandle
  Arg3: segmentation hashref
    - name => string
    - celltypes => arrayref of strings
    - states => arrayref of strings
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_segmentation_summaries_2 {
  my ($options, $file, $segmentation) = @_;

  my $name = ${segmentation}->{name};
  my $states = $segmentation->{states};
  my $state_count = scalar(@$states);
  my @celltypes = keys %{$segmentation->{celltypes}};
  my $celltype_count = scalar(@celltypes);
  my $state_string = "";
  my $state;

  foreach $state (@$states) {
    $state_string .= " $state=$state";
  }

  print $file "\n";

  print $file "track ${name}SegmentationSummaries\n";
  print $file "compositeTrack on\n";
  print $file "visibility dense\n";
  print $file "shortLabel $name Segmentation Summaries\n";
  print $file "longLabel Overlap summaries of each state across $celltype_count segmentations\n";
  print $file "subGroup1 state State $state_string\n";
  print $file "subGroup2 source Source ERB=Ensembl_Regulatory_Build\n";
  print $file "dimensions dimY=state\n";
  print $file "dragAndDrop subTracks\n";
  print $file "priority 7\n";
  print $file "type bigWig\n";
  print $file "noInherit on\n";


  foreach $state (@$states) {
    my $color = $COLORS{$options->{assignments}->{$segmentation->{name}}->{$state}};

    print $file "\n";

    print $file "\ttrack ${name}_$state\n";
    print $file "\tbigDataUrl $options->{assembly}/segmentation_summaries/$segmentation->{name}/$state.bw\n";
    print $file "\tparent ${name}SegmentationSummaries off\n";
    print $file "\tshortLabel $state summary\n";
    print $file "\tlongLabel Overlap summary of $state state across $celltype_count segmentations\n";
    print $file "\ttype bigWig\n";
    print $file "\tautoscale off\n";
    print $file "\tmaxHeightPixels 64:32:16\n";
    print $file "\tviewLimits 0:$celltype_count\n";
    print $file "\tcolor $color\n";
    print $file "\tsubGroups state=$state source=ERB\n";
  }
}

=head2 make_track_hub_projected_segmentations

  Description: Track hub Segmentation_Summaries composite track
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
  Arg2: filehandle
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_projected_segmentations {
  my ($options, $file) = @_;
  my $celltype;

  my $celltype_list = "";
  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach $celltype (keys %{$segmentation->{celltypes}}) {
      my $short = $celltype;
      $short =~ s/[_-]//g;
      $celltype_list .= " $short=$celltype";
    }
  }

  print $file "\n";

  print $file "track ProjectedSegments\n";
  print $file "compositeTrack on\n";
  print $file "visibility squish\n";
  print $file "shortLabel Cell Activity\n";
  print $file "longLabel Ensembl Regulatory Build with cell type specific activity\n";
  print $file "subGroup1 cellType Cell_Type $celltype_list\n";
  print $file "subGroup2 source Source ERB=Ensembl_Regulatory_Build\n";
  print $file "dimensions dimY=cellType\n";
  print $file "sortOrder cellType=+\n";
  print $file "dragAndDrop subTracks\n";
  print $file "priority 6\n";
  print $file "type bigBed 9\n";
  print $file "noInherit on\n";

  foreach my $segmentation (@{$options->{segmentations}}) {
    my $name = $segmentation->{name};
    foreach $celltype (keys %{$segmentation->{celltypes}}) {
      my $short = $celltype;
      $short =~ s/[_-]//g;
      $celltype_list .= " $short=$celltype";

      print $file "\n";

      print $file "\ttrack ${name}_${celltype}_projected\n";
      print $file "\tbigDataUrl $options->{assembly}/projected_segmentations/$celltype.bb\n";
      print $file "\tparent ProjectedSegments off\n";
      print $file "\tshortLabel $celltype\n";
      print $file "\tlongLabel Projection of the Ensembl Regulatory build onto $celltype\n";
      print $file "\ttype bigBed 9\n";
      print $file "\titemRgb on\n";
      print $file "\tsubGroups cellType=$short source=ERB\n";
    }
  }
}

=head2 make_track_hub_tfbs

  Description: Track hub TFBS_Summaries composite track
  Arg1: $options: hash containing:
    - $options->{assembly} UCSC name of the assembly
    - $options->{trackhub_dir}
    - $options->{segmentations} List of hashes containing:
      . $segmentation->{name}
      . $segmentation->{states}
      . $segmentation->{celltypes}
  Arg2: filehandle
  Returntype: undef
  Side effects: Writes into:
  * $options->{trackhub_dir}/trackDb.txt

=cut

sub make_track_hub_tfbs {
  my ($options, $file) = @_;
  my $tf;
  my @tfs = ();

  foreach my $file (glob "$options->{trackhub_dir}/tfbs/*.bw") {
    $file =~ s/.bw//;
    push @tfs, basename $file;
  }

  my $tf_list = "";
  foreach $tf (@tfs) {
    my $short = $tf;
    $short =~ s/[_-]//g;
    $tf_list .= " $tf=$short";
  }

  print $file "\n";

  print $file "track TfbsPeaks\n";
  print $file "compositeTrack on\n";
  print $file "visibility full\n";
  print $file "shortLabel TFBS Peaks Summaries\n";
  print $file "longLabel Overlap summary of ChIPSeq binding peaks across available datasets\n";
  print $file "subGroup1 antibody Antibody $tf_list\n";
  print $file "subGroup2 source Source ERB=Ensembl_Regulatory_Build\n";
  print $file "dimensions dimY=antibody\n";
  print $file "sortOrder antibody=+\n";
  print $file "dragAndDrop subTracks\n";
  print $file "priority 9\n";
  print $file "type bigWig\n";
  print $file "noInherit on\n";

  foreach $tf (@tfs) {
    my $short = $tf;
    $short =~ s/-//g;

    print $file "\n";

    print $file "\ttrack $tf\n";
    print $file "\tbigDataUrl $options->{assembly}/tfbs/$tf.bw\n";
    print $file "\tparent TfbsPeaks off\n";
    print $file "\tshortLabel $tf \n";
    print $file "\tlongLabel Overlap summary of $tf ChipSeq binding peaks across available datasets\n";
    print $file "\ttype bigWig\n";
    print $file "\tautoscale off\n";
    print $file "\tmaxHeightPixels 64:32:16\n";
    print $file "\tviewLimits 0:1\n";
    my $color;
    if ($tf eq 'CTCF') {
      $color = $COLORS{ctcf};
    } else {
      $color = $COLORS{tfbs};
    }
    print $file "\tcolor $color\n";
    print $file "\tsubGroups antibody=$short source=ERB\n";
  }
}

1;
