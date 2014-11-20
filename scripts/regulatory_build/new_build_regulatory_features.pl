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

new_build_regulatory_features.pl

=head1 SYNOPSIS

Automatic run:
new_build_regulatory_features.pl -o ./ -a hg38 -u http://ngs.sanger.ac.uk/production/ensembl/regulation/GRCh38 -h $HOST -P $PORT -u $USER -p $PASS

Where:
	* -o: directory for output
	* -a: UCSC assembly name (for trackHub)
	* -u: Base url for trackHub location.
	* -h: host server
	* -D: host database
	* -P: host port
	* -u: host login 
	* -p: host password 

Manual run:
new_build_regulatory_features.pl -o ./ -d dump.txt -a hg38 -l chrom_sizes.txt -t tss.bed -g exons.bed -u http://ngs.sanger.ac.uk/production/ensembl/regulation/GRCh38

Where:
	* -o: directory for output
	* -a: UCSC assembly name (for trackHub)
	* -u: Base url for trackHub location.
	* -l: tab delimited file, each line contains a chromosome name followed
		by its length.  [Overrides database info]
	* -t: Bed file with TSS. Can be Ensembl transcript TSS, CAGE tags... [Overrides database info]
	* -g: Bed file with Exons. A BioMart dump would work. [Overrides database info]
	* -d: dump file (described below) [Added to database datasets]

=head1 DESCRIPTION

Generates regulatory features based on overlaps between segmentations
on different cell types.

In particular you will need a dump file which describes all the available 
experimental data. The dump file is tab delimited, it contains two types of
entries:
* ChIPseq peaks
peak	$assay	$cell	$location_bed_or_bigBed
The assay type is generally a TF antibody name or DNAse

* Segmentations
segmentation	$name	$type	$location
	- name: is just a free string, which will be used a directory name (avoid special characters and spaces).
	- type: is ChromHMM or Segway
	- location: directory which contains a bunch of bed files, each bedfile named $celltype.bed. In addition 
	the directory must contain an emissions file (emissions*.txt for ChromHMM, *.tab for Segway).

=cut

# TODO Allow 2 TFS to create feature in the absence of DNase

use strict;
use warnings;
use File::Path qw(mkpath);
use File::Basename;
use File::Temp;
use Storable;
use Data::Dumper;

########################################################
## Global constants 
########################################################

# A bunch of arbitrary cutoffs
# All of these cutoff are applied to enrichment ratios (i.e. 1 == neutral)
our $weak_cutoff = 2;
our $strong_cutoff = 5;
our $very_strong_cutoff = 10;
# The functions are the annotation labels of the build built from segmentation data
# Note that the order determines the order of the build. They must be in that order
our @functions = ('tss', 'proximal', 'distal', 'ctcf');
# The annotations are the labels used in the build. TFBS and DNAse are special in that they is not 
# built from segmentation data, rather from Chip-Seq peaks
our @annotations = (@functions, 'tfbs', 'dnase');
# Typical histone marks used to detect repression. Note that the strings are normalised 
# via the clean_name function below
our @repressed_marks = ('H3K27ME3');
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
  poised => "192,0,190"
  repressed => "127,127,127"
);

# These states describe the features at the cell-type level
our %states = (
  active => 0,
  poised => 1,
  repressed => 2,
  inactive => 3
);

########################################################
## Toplevel pipeline 
########################################################

main();

sub main {
  print "Entering New Regulatory Build pipeline\n";
  # Read command line
  my $options = get_options();
  # Read dump file with file locations and metadata
  get_metadata($options);
  # Compute summaries of TB binding peaks
  compute_tf_probs($options);
  # Compute summaries of segmentations
  extract_segmentation_state_summaries($options);
  # Select the segmentation states that best match a label
  label_segmentation_states($options);
  # Color the input segmentations based on the label assigments above
  make_segmentation_bedfiles($options);
  # Set cutoffs to optimise the quality of the build (using TF binding as an indicator_
  set_cutoffs($options);
  # Compute MultiCell features
  compute_regulatory_features($options);
  # Determine which features are active in which cell type
  compute_states($options);
  # Prepare trackhub header files
  make_track_hub($options);
  print "Exiting New Regulatory Build pipeline\n";
}

########################################################
## Option reading
## Parses the command line, checks the options
########################################################

sub get_options {
  my $options = read_command_line();
  check_options($options);
  return $options;
}

use Getopt::Long;

sub read_command_line {
  my %options = ();

  GetOptions(\%options, "help|h=s", "tmp|T=s", "out|o=s", "dump|d=s", "assembly|a=s", "chrom_lengths|l=s", "tss|t=s", "exons|g=s", "url|u=s", "host|h=s", "port|P=s", "db|D=s", "user|u=s", "pass|p=s", "mask|m=s");

  $options{output_dir} = $options{out};
  $options{working_dir} = $options{tmp};

  return \%options;
}

sub check_options {
  my $options = shift;

  die('Output directory not defined') unless defined $options->{out};
  die('Assembly name not defined') unless defined $options->{assembly};
  die('Base url not defined') unless defined $options->{url};

  if (! defined $options->{host}) {
    die('Ensembl dumps not described') unless defined $options->{dump};
    die('Chromosome lengths not provided') unless defined $options->{chrom_lengths};
    die('TSS file not provided') unless defined $options->{tss};
    die('Exon file not provided') unless defined $options->{exons};
    die('Chromosome length file not provided') unless defined $options->{exons};
  } else {
    die('Database not defined') unless defined $options->{db};
    die('User account not defined') unless defined $options->{user};
    $options->{db_adaptor} = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
	-user   => $options->{user},
	-dbname => $options->{db},
	-host   => $options->{host},
	-pass   => $options->{pass},
	-port   => $options->{port}
    );
    $options->{dnadb_adaptor} = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
	-user   => $options->{db_adaptor}->dnadb_user(),
	-dbname => $options->{db_adaptor}->dnadb_name(),
	-host   => $options->{db_adaptor}->dnadb_host(),
	-pass   => $options->{db_adaptor}->dnadb_pass(),
	-port   => $options->{db_adaptor}->dnadb_port()
    );
  }

  if (!defined $options->{working_dir}) {
    $options->{working_dir} = "$options->{output_dir}/tmp/";
  }
  mkpath $options->{working_dir};

  $options->{trackhub_dir} = "$options->{output_dir}/$options->{assembly}";
  mkpath $options->{trackhub_dir};
}

########################################################
## Compute forcing:
## * As long as all checked files exist, they are not recomputed
## * As soon as a checked file needs to be recomputed, everything
## downstream is recomputed
########################################################

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

########################################################
## BigBed creation
## Wrapper functions for BigFile indexing: 
## Params:
## - path to file
## Actions:
## - creates indexed file (BigWig or BigBed)
## - deletes flat file
########################################################

sub convert_to_bigBed {
  my ($options, $file) = @_;
  my $new = $file;
  $new =~ s/\.bed$/.bb/;
  run("bedToBigBed $file $options->{chrom_lengths} $new") ;
  unlink $file;
}

sub convert_to_bigWig {
  my ($options, $file) = @_;
  trim_bed_to_chrom_lengths($options, $file);
  my $new = $file;
  $new =~ s/\.wig$/.bw/;
  run("wigToBigWig $file $options->{chrom_lengths} $new") ;
  unlink $file;
}

########################################################
## System calls 
## Wrapper function for system calls
## Params:
## - Command line
## Actions:
## - Runs command, prints out error in case of failure
########################################################

sub run {
  my ($cmd) = @_;
  print "Running $cmd\n";
  system($cmd) && die("Failure when running command\n$cmd\n");
}

########################################################
## Removing unwanted characters 
## Quick string normalisation function tor remove weird 
## characters froms file names and remove variants
########################################################

sub clean_name {
  my $string = shift;
  $string =~ s/[\-\(\)]//g;
  $string =~ s/_.*//g;
  return uc($string);
}

########################################################
## Serialising functions
## Reads a list of strings from a text file
########################################################

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

########################################################
## Trimming to chromosome length 
## Ensures that all features of a bed file fit within the 
## exact coordinates of a chromosome
## This is necessary for BigFile indexing, as ChromHMM 
## rounds up coordinates to the nearest 200bp mark. 
## Params:
## - options hash
## - filepath
## Actions:
## - Overwrites file with corrected coordinates
########################################################

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
    if ($items[1] < $lengths->{$items[0]} && $items[2] > 0) {
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
}

sub read_chrom_lengths {
  my ($options) = @_;
  my ($fh, $line);

  open $fh, "<", $options->{chrom_lengths};
  my %lengths = ();

  while ($line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    $lengths{$items[0]} = $items[1]; 
  }
  close $fh;

  return \%lengths;
}

########################################################
## Parsing Ensembl dumps 
## The Ensembl dumpfile is assumed to be a tab delimited 
## file with two types of lines:
## peaks  $TF  $celltype  $path
## segmentation  $type  $path
## The type of a segmentation is currently only ChromHMM
## or Segway
## but other segmentation software could crop up.
## The path of a segmentation is the path to the directory
## containing all the output files.
## Params:
## - options hash
## Actions:
## - Fills up options hash
########################################################

sub get_metadata {
  my $options = shift;

  $options->{cell_type_tfs} = {};
  $options->{cell_type_open} = {};
  $options->{peak_calls} = {};
  $options->{segmentations} = [];

  if (defined $options->{dump}) {
    read_dump($options);
  }

  if (defined $options->{db_adaptor}) {
    fetch_metadata($options);
  }

  if (!defined $options->{chrom_lengths}) {
    create_chrom_lengths($options);
  }

  if (!defined $options->{tss}) {
    create_tss($options);
  }

  if (!defined $options->{exons}) {
    create_exons($options);
  }

  if (!defined $options->{mask}) {
    create_mask($options);
  }
}

sub read_dump {
  my ($options) = @_;
  print "Reading $options->{dump}\n";

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
}

sub fetch_metadata {
  my ($options) = @_;
  my @slices = sort {$a->name cmp $b->name} @{$options->{dnadb_adaptor}->get_SliceAdaptor->fetch_all('toplevel', undef, undef, False)};
  foreach my $featureSet (@{$options->{db_adaptor}->get_adaptor("FeatureSet")->fetch_all_by_feature_class("annotated")}) {
    my $tf = $featureSet->feature_type->name;
    my $cell = $featureSet->cell->name;
    my $dir = "$options->{working_dir}/peaks/$tf/$cell/";
    mkpath $dir;

    my $fh = File::Temp->new(DIR => $dir, SUFFIX => '.bed');
    foreach my $slice (@slices) {
      foreach my $feature (sort {$a->start <=> $b->start} @{$featureSet->get_Features_by_Slice($slice)});
	print $fh join("\t", ($slice->name, $feature->start - 1, $feature->end));
      }
    }
    record_peak_dir($options, $tf, $cell, $fh->filename);
    close $fh
  }
}

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

sub create_chrom_lengths {
  my ($options) = @_;
  $options->{chrom_lengths} = "$options->{working_dir}/chrom_lengths.txt";
  open my $fh, ">", $options->{chrom_lengths};
  fetch_chrom_lengths($options, $fh);
  close $fh
}

sub fetch_chrom_lengths {
  my ($options, $fh) = @_;
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  my @slices = @{ $slice_adaptor->fetch_all('toplevel', undef, undef, False) };

  foreach my $slice (@slices) {
    print $fh join("\t", ($slice->seq_region_name() - 1, $slice->end() - $slice->start())) . "\n";
  }
}

sub create_tss {
  my ($options) = @_;
  $options->{tss} = "$options->{working_dir}/tss.txt";
  open my $fh, ">", $options->{tss};
  fetch_tss($options, $fh);
  close $fh
}

sub fetch_tss {
  my ($options, $fh) = @_;
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  my @tss_coords = ();
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, False) }) {
    foreach my $gene (@{$slice->get_all_Genes()}) {
      foreach my $transcript (@{$gene->get_all_Transcripts()}) {
	if ($transcript->strand() > 0) {
	  push @tss_coords, [$slice->seq_regions_name(), $transcript->start() - 1, $transcript->start()];
        } else {
	  push @tss_coords, [$slice->seq_regions_name(), $transcript->end() - 1, $transcript->end()];
	}
      }
    }
  }

  my @sorted_tss_coords = sort {comp_coords($a, $b)} @tss_coords;

  foreach my $exons_coord (@sorted_tss_coords) {
    print $fh, join("\t", @{$tss_coord})."\n";
  }
}

sub create_exons {
  my ($options) = @_;
  $options->{exons} = "$options->{working_dir}/exons.txt";
  open my $fh, ">", $options->{exons};
  fetch_exons($options, $fh);
  close $fh
}

sub fetch_exons {
  my ($options, $fh) = @_;
  my @exon_coords = ();
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, False) }) {
    foreach my $gene (@{$slice->get_all_Genes()}) {
      foreach my $transcript (@{$gene->get_all_Transcripts()}) {
        foreach my $exon (@{$transcript->get_all_Exons()}) {
	  push @exon_coords, [$slice->seq_regions_name(), $exon$exon->start() - 1, $exon->end()];
	}
      }
    }
  }

  my @sorted_exon_coords = sort {comp_coords($a, $b)} @exon_coords;

  foreach my $exons_coord (@sorted_exon_coords) {
    print $fh, join("\t", @{$exon_coord})."\n";
  }
}

sub create_mask {
  my ($options) = @_;
  $options->{mask} = "$options->{working_dir}/mask.txt";
  open my $fh, ">", $options->{mask};
  fetch_mask($options, $fh);
  close $fh
}

sub fetch_mask {
  my ($options, $fh) = @_;
  my @mask_coords = ();
  my $slice_adaptor = $options->{dnadb_adaptor}->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, undef, False) }) {
    foreach my $mask (@{$slice->get_all_MiscFeatures('encode_excluded')}) {
      push @mask_coords, [$slice->seq_regions_name(), $mask$mask->start(), $mask->end() - 1];
    }
  }

  my @sorted_mask_coords = sort {comp_coords($a, $b)} @mask_coords;

  foreach my $mask_coords (@sorted_mask_coords) {
    print $fh, join("\t", @{$mask_coord})."\n";
  }
}


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

########################################################
## Computing TF binding probs 
##
## Params:
## * options hash with values:
##   - $options->{trackhub_dir} (path to trackhub dir)
##   - $options->{working_dir} (path to working tmp dir)
##   - $options->{celltype_open} (hash which assigns a list of peak calls to each cell type)
##   - $options->{celltype_tfs} (hash which assigns a list of peak calls to each cell type)
##   - $options->{peak_calls} (hash which assigns a list of peak calls to each )
## 
## Writes into:
## * Celltype specific summaries in:
##   $options->{working_dir}/celltype_tf/$celltype.bed
## * TF specific summaries in:
##   $options->{trackhub_dir}/tfbs/$tf.wig
## * An overall TF summary, computed as a disjunction of the 
## TF specific probabilities, treated as independent variables.
##   $options->{trackhub_dir}/overview/all_tfbs.bw
########################################################

sub compute_tf_probs {
  my $options = shift;

  if (!must_compute($options,"$options->{trackhub_dir}/overview/all_tfbs.bw")) {
    print "TF binding probs already there, skipping calculations...\n";
    return;
  }

  print "Computing TF binding tracks\n";
  compute_celltype_tf_sites($options);
  my $tf_probs = compute_antibody_specific_probs($options);
  compute_global_tf_prob($options, $tf_probs);
}

sub compute_celltype_tf_sites {
  my ($options) = @_;
  my $celltype;

  mkdir "$options->{working_dir}/celltype_tf/";
  mkdir "$options->{working_dir}/celltype_dnase/";

  foreach $celltype (keys %{$options->{cell_type_tfs}}) {
    compute_celltype_tf_sites_2($options, $celltype);
  }
}

sub compute_celltype_tf_sites_2 {
  my ($options, $celltype) = @_;

  my $output2 = "$options->{working_dir}/celltype_tf/$celltype.bed";  
  my $tf_exps = $options->{cell_type_tfs}->{$celltype};
  my $open_exps = $options->{cell_type_open}->{$celltype};

  if (defined $open_exps) {
    my $output1 = "$options->{working_dir}/celltype_dnase/$celltype.bed"; 
    run("wiggletools write_bg $output1 unit sum " . join(" ", @{$open_exps}));
    $options->{celltype_dnase}->{$celltype} = $output1;
    trim_bed_to_chrom_lengths($options, $output1);
    run("wiggletools write_bg $output2 unit mult sum ".join(" ", @{$tf_exps}). " : $output1 ");
  } else {
    run("wiggletools write_bg $output2 unit sum ".join(" ", @{$tf_exps}));
  }
  trim_bed_to_chrom_lengths($options, $output2);
  $options->{celltype_tfbs}->{$celltype} = $output2;
}

sub compute_antibody_specific_probs {
  my $options = shift;
  my @tfs = keys %{$options->{peak_calls}};
  my @tf_probs = ();
  my $tf;

  mkdir "$options->{trackhub_dir}/tfbs/";

  foreach $tf (@tfs) {
    my $res = compute_antibody_specific_prob($options, $tf);
    if (defined $res) {
      push @tf_probs, $res;
    }
  }

  return \@tf_probs;
}

sub compute_antibody_specific_prob {
  my ($options, $tf) = @_;
  my @cells = keys(%{$options->{peak_calls}->{$tf}});
  my $output = "$options->{trackhub_dir}/tfbs/$tf.wig";
  my @exps = ();

  foreach my $cell (@cells) {
     if (defined $options->{cell_type_open}->{$cell}) {
  my $tf_peaks = $options->{peak_calls}->{$tf}->{$cell};
  my $open_peaks = $options->{cell_type_open}->{$cell};
        push @exps, "unit mult sum " . join(" ", @$tf_peaks) . " : sum " . join(" ", @$open_peaks) . " : :";
     }
  }

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

sub compute_global_tf_prob {
  my ($options, $probs) = @_;
  mkdir "$options->{trackhub_dir}/overview/";
  my $output = "$options->{trackhub_dir}/overview/all_tfbs.wig";

  run("wiggletools write_bg $output offset 1 scale -1 mult map offset 1 map scale -1 " . join(" ", @$probs)) ;
  convert_to_bigWig($options, $output);
}

########################################################
## Extracting segmentation summaries 
## Computes for each segmentation and each state the 
## number of celltypes with a given state at each position
##
## Params:
## - An $options hash with values:
##   * $options->{trackhub_dir}
##   * $options->{working_dir}
##   * $options->{segmentations} List of hashes containing:
##     . $segmentation->{name} String
##     . $segmentation->{location} Path to directory
##
## It fills out:
##   * $segmentation->{states} List of strings 
##   * $segmentation->{celltypes} List of strings
##
## Writes into:
##   * List of states in segmentation:
##     $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt
##   * List of cell-types in segmentation:
##     $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt
##   * Sorted copy of each cell type's segmentation:
##     $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed
##   * A summary statistic of each state:
##     $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
########################################################

sub extract_segmentation_state_summaries {
  my $options = shift;
  my $segmentation;

  mkdir "$options->{trackhub_dir}/segmentation_summaries/";
  mkdir "$options->{working_dir}/segmentation_summaries/";

  foreach $segmentation (@{$options->{segmentations}}) {
    extract_segmentation_state_summaries_2($options, $segmentation);
  }
}

sub extract_segmentation_state_summaries_2 {
  my ($options, $segmentation) = @_;

  mkdir "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/";
  mkdir "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/";

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway') {
    extract_ChromHMM_state_summaries($options, $segmentation);
  } else {
    die("Unknown segmentation format $segmentation->{type}\n");
  }
}

sub extract_ChromHMM_state_summaries {
  my ($options, $segmentation) = @_;
  print "Going through output of segmentation $segmentation->{name}\n";
  my @bedfiles = glob "$segmentation->{location}/*.bed";
  $segmentation->{states} = extract_ChromHMM_states($options, $segmentation, \@bedfiles);
  $segmentation->{celltypes} = extract_ChromHMM_cells($options, $segmentation, \@bedfiles);

  foreach my $cell (@{$segmentation->{celltypes}}) {
    my $summary = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed";
    if (must_compute($options,$summary)) {
      run("sort -k1,1 -k2,2n $segmentation->{location}/$cell.bed > $summary") ;
    }
  }

  foreach my $state (@{$segmentation->{states}}) {
    extract_ChromHMM_state_summary($options, $segmentation, $state);
  }
}

sub extract_ChromHMM_states {
  my ($options, $segmentation, $files) = @_;
  print "Extracting state names\n";
  my $output = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/states.txt";
  if (must_compute($options,$output)) {
    run("cat " . join(" ", @$files) . " | grep -v '^track\>' | cut -f4 | uniq | sort | uniq > $output");
  }
  return deserialise($output);
}

sub extract_ChromHMM_cells {
  my ($options, $segmentation, $files) = @_;
  print "Extracting cell names\n";
  my $output = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/cells.txt";
  if (must_compute($options,$output)) {
    open my $fh, ">", $output;
    foreach my $file (@$files) {
      my $cell = basename $file;
      $cell =~ s/.bed$//;
      $cell = clean_name($cell);
      print $fh "$cell\n";
    }
    close $fh;
  }
  return deserialise($output);
}

sub extract_ChromHMM_state_summary {
  my ($options, $segmentation, $state) = @_;
  my $celltypes = $segmentation->{celltypes};
  my $cell;

  my @summaries = ();

  if (!must_compute($options,"$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw")) {
    print "Summary of state $state in segmentation $segmentation->{name} already computed, skipping...\n";
    return;
  } else {
    print "Computing summary of state $state in segmentation $segmentation->{name}.\n";
  }

  mkdir "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/";

  foreach $cell (@{$segmentation->{celltypes}}) {
    my $source = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$cell.bed";
    my $summary = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$cell.bed";
    run("awk \'\$4 == \"$state\" \' $source > $summary") ;
    push @summaries, $summary;
  }

  my $summary = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.wig";
  run("wiggletools write_bg $summary sum " . join(" ", @summaries));
  convert_to_bigWig($options, $summary);
}

########################################################
## Characterising segmentation states 
## This first annotation assigns a label to each state of 
## each segmentation, partly to color the segmentations 
## partly to inform the next step.
##
## Params:
## * A $options hash which contains:
##   - $options->{working_dir}
##   - $options->{trackhub_dir}
##   - $options->{exons} Path to bed file with exonic regions
##   - $options->{tss} Path to bed file with known TSSs
##   - $options->{peak_calls}->{CTCF}->{$cell} A list of filepaths to CTCF peak calls of celltype $cell
##   - $options->{segmentations} where each value is a hash containing:
##     . $segmentation->{name} 
##     . $segmentation->{states} 
##
## It expects:
## * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
## * ChromHMM segmentations: $segmentation->{location}/*.tab (1 file) Emission file
## * Segway segmentations: $segmentation->{location}/*.txt (1 file) Emission file
##
## It fills out:
## * $options->{assignments}->{$segmentation->{name}} Hash which assigns a label to each state
##
## It creates:
## * $options->{working_dir}/assignments.txt (Tab delimited file)
########################################################

sub label_segmentation_states {
  my $options = shift;
  my ($segmentation, $state);
  
  if (!must_compute($options,"$options->{working_dir}/assignments.txt")) {
    print "Loading pre-existing assignments $options->{working_dir}/assignments.txt\n";
    $options->{assignments} = load_assignments( "$options->{working_dir}/assignments.txt");
  } else {
    print "Assigning states to function\n";
    $options->{assignments} = {};
    foreach $segmentation (@{$options->{segmentations}}) {
      select_segmentation_states_2($options, $segmentation);
    }
    dump_assignments($options, "$options->{working_dir}/assignments.txt");
  }

  foreach $segmentation (@{$options->{segmentations}}) {
    my $assignments = $options->{assignments}->{$segmentation->{name}};
    foreach $state (keys %{$assignments}) {
      print "Assignment\t$segmentation->{name}\t$state\t$assignments->{$state}\n";
    }
  }
}

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

sub select_segmentation_states_2 {
  my ($options, $segmentation) = @_;
  my $state;

  $options->{assignments}->{$segmentation->{name}} = {};

  print "Entering segmentation $segmentation->{name}\n";

  compute_overlaps($options, $segmentation);

  foreach $state (@{$segmentation->{states}}) {
    label_segmentation_state($options, $segmentation, $state);
  }
}

sub label_segmentation_state {
  my ($options, $segmentation, $state) = @_;
  my $overlaps = $segmentation->{overlaps};
  my $assignments = $options->{assignments}->{$segmentation->{name}};

  if ($overlaps->{ctcf}->{$state} > .25) {
    $assignments->{$state} = 'ctcf';
  } elsif ($overlaps->{repressed}->{$state} > $segmentation->{repressed_cutoff}) {
    if ($overlaps->{tfbs}->{$state} < 1) {
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

sub compute_overlaps {
  my ($options, $segmentation) = @_;

  mkpath "$options->{working_dir}/overlaps/$segmentation->{name}";

  foreach my $test ('ctcf', 'tss', 'gene', 'tfbs') {
    compute_overlap_scores($options, $segmentation, $test);
  }
  if ($segmentation->{type} eq 'ChromHMM') {
    compute_ChromHMM_repressed_scores($options, $segmentation);
  } elsif ($segmentation->{type} eq 'Segway') {
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

sub compute_overlap_scores {
  my ($options, $segmentation, $test) = @_;
  my $state;
  my $output = "$options->{working_dir}/overlaps/$segmentation->{name}/$test.txt";
  if (!must_compute($options,$output)) {
    print "Retrieving pre-computed overlaps for test $test\n";
    $segmentation->{overlaps}->{$test} = retrieve($output);
  } else {
    print "Computing overlaps for test $test\n";

    $segmentation->{overlaps}->{$test} = {};
    foreach $state (@{$segmentation->{states}}) {
      compute_overlap_score($options, $segmentation, $test, $state);
    }
    store $segmentation->{overlaps}->{$test}, $output;
  }
}

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
    print "Running wiggletools pearson - $reference $file\n";
    $segmentation->{overlaps}->{$test}->{$state} = `wiggletools pearson - $reference $file`;
  } else {
    $segmentation->{overlaps}->{$test}->{$state} = compute_enrichment_between_files($reference, $file);
  }
  
  chomp $segmentation->{overlaps}->{$test}->{$state};
}

sub compute_enrichment_between_files {
  my ($reference, $file) = @_;

  my $auc = `wiggletools AUC - mult $reference unit $file`;
  my $breadth = `wiggletools AUC - unit $file`;
  my $ref_auc = `wiggletools AUC - $reference`;

  if ($breadth == 0) {
    return 0;
  } else {
    return ($auc / $breadth) / ($ref_auc / 3.09569e9);
  }
}

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
    print "Emission\t$segmentation->{name}\trepressed\t$items[0]\t$sum\n";
    if ($sum > $max) {
      $max = $sum;
    }
  }
  close $fh;

  $segmentation->{repressed_cutoff} = $max / 4;
}

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
    print "Emission\t$segmentation->{name}\trepressed\t$state\t$segmentation->{overlaps}->{repressed}->{$state}\n";
  }
  $segmentation->{repressed_cutoff} = $max / 4;
}

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
    print "Emission\t$segmentation->{name}\trepressed\t$state\t$segmentation->{overlaps}->{repressed}->{$state}\n";
  }

  $segmentation->{repressed_cutoff} = $max / 4;
}

########################################################
## Preparing segmentation bedfiles 
## This function adds in color information to the bedfiles,
## based on the assignment of states done above
##
## Params:
## * A $options hash which contains:
##   - $options->{trackhub_dir}
##   - $options->{working_dir}
##   - $options->{segmentations} list which each element is a hash containing:
##     . $segmentation->{name} String
##     . $segmentation->{type} (ChromHMM or Segway)
##     . $segmentation->{states} List of strings
##     . $segmentation->{celltypes} List of strings
##
## It expects files:
## * $options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed
##
## It creates:
## * $options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb
########################################################

sub make_segmentation_bedfiles {
  my $options = shift;
  my $segmentation;
  mkdir "$options->{trackhub_dir}/segmentations/";
  foreach $segmentation (@{$options->{segmentations}}) {
    make_segmentation_bedfiles_2($options, $segmentation);
  }
}

sub make_segmentation_bedfiles_2 {
  my ($options, $segmentation) = @_;
  my $celltype;

  if (must_compute($options,"$options->{trackhub_dir}/segmentations/$segmentation->{name}")) {
    mkdir "$options->{trackhub_dir}/segmentations/$segmentation->{name}";
  }

  foreach $celltype (@{$segmentation->{celltypes}}) {
    make_segmentation_bedfile($options, $segmentation, $celltype);
  }
}

sub make_segmentation_bedfile {
  my ($options, $segmentation, $celltype) = @_;

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway') {
    make_ChromHMM_bedfile($options, $segmentation, $celltype);
  } else {
    die("Unknown segmentation type $segmentation->{type}\n");
  }
}

sub make_ChromHMM_bedfile {
  my ($options, $segmentation, $celltype) = @_;
  my $output = "$options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bed";

  if (!must_compute($options,"$options->{trackhub_dir}/segmentations/$segmentation->{name}/$celltype.bb")) {
    print "Colorised segmentation $output already exists, skipping...\n";
    return;
  }

  my @files = ();
  foreach my $state (@{$segmentation->{states}}) {
    push @files, make_ChromHMM_state_bedfile($options, $segmentation, $celltype, $state);
  }

  run("sort -m ".join(" ", @files)." -k1,1 -k2,2n > $output");
  convert_to_bigBed($options, $output);
}

sub make_ChromHMM_state_bedfile {
  my ($options, $segmentation, $celltype, $state) = @_;
  my $input = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed";
  die ("Could not find $input\n") unless -e $input;
  my $output = "$options->{working_dir}/segmentations/$segmentation->{name}/$state/$celltype.bed";

  mkpath "$options->{working_dir}/segmentations/$segmentation->{name}/$state/";
  my $assigned_annotation = $options->{assignments}->{$segmentation->{name}}->{$state};
  my $color = $COLORS{$assigned_annotation};
  my $awk_string = "BEGIN {OFS=\"\\t\"} {\$4 = \"${state}_${assigned_annotation}_\"NR; \$5=1000; \$6=\".\"; \$7=\$2; \$8=\$3; \$9=\"$color\"; print}";
  run("awk \'$awk_string\' $input > $output") ;

  return $output;
}

########################################################
## Setting cutoffs 
## This functions first selects which states are directly useful
## to infer TF binding, then determines an optimal
## cutoff so as to best fit the overall transcription factor binding
## probability:
##
## Params:
## * A $options hash which contains:
##   - $options->{trackhub_dir}
##   - $options->{working_dir}
##   - $options->{segmentations} list which each element is a hash containing:
##     . $segmentation->{name} String
##     . $segmentation->{type} (ChromHMM or Segway)
##     . $segmentation->{states} List of strings
##     . $segmentation->{celltypes} List of strings
##
## It expects:
## * $options->{trackhub_dir}/overview/all_tfbs.bw
## * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw
##
## It computes:
## * $options->{selected_states}->{$segmentation->{name}}->{$function} 
##   List of state names that given an appropriate cutoff can have a 5x enrichment 
##   in the reference marker
## * $options->{cutoffs}->{$segmentation->{name}}->{$function} Integer
##   Cutoff on the sum of values of selected states that maximises 
##   the F-score of detection of TF binding 
##
## It creates:
## * $options->{working_dir}/selected serialised hash (see above)
## * $options->{working_dir}/cutoffs serialised hash (see above)
########################################################

sub set_cutoffs {
  my $options = shift;
  my $cutoffs_file = "$options->{working_dir}/cutoffs";
  my $selected_file = "$options->{working_dir}/selected";

  if (!must_compute($options,$cutoffs_file) && !must_compute($options,$selected_file)) {
    $options->{cutoffs} = retrieve $cutoffs_file;
    $options->{selected_states} = retrieve $selected_file;
  } else {
    $options->{cutoffs} = {};
    $options->{selected_states} = {};
    foreach my $segmentation (@{$options->{segmentations}}) {
      select_segmentation_cutoffs($options, $segmentation);
    }
    store $options->{cutoffs}, $cutoffs_file;
    store $options->{selected_states}, $selected_file;
  }

  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach my $function (keys %{$options->{selected_states}->{$segmentation->{name}}}) {
      print "Selected\t$segmentation->{name}\t$function\t". join(" ", @{$options->{selected_states}->{$segmentation->{name}}->{$function}})."\n";
    }
  }

  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach my $function (keys %{$options->{cutoffs}->{$segmentation->{name}}}) {
      print "Cutoff\t$segmentation->{name}\t$function\t$options->{cutoffs}->{$segmentation->{name}}->{$function}\n";
    }
  }
}

sub select_segmentation_cutoffs {
  my ($options, $segmentation) = @_;
  $options->{selected_states}->{$segmentation->{name}} = {};
  $options->{cutoffs}->{$segmentation->{name}} = {};

  foreach my $function (@functions) {
    $options->{selected_states}->{$segmentation->{name}}->{$function} = select_relevant_states($options, $segmentation, $function);
    $options->{cutoffs}->{$segmentation->{name}}->{$function} = select_segmentation_cutoff($options, $segmentation, $function);
  }
}

sub select_relevant_states {
  my ($options, $segmentation, $function) = @_;
  my %states = ();
  my $assignments = $options->{assignments}->{$segmentation->{name}};
  foreach my $state (@{$segmentation->{states}}) {
    if ($assignments->{$state} eq $function) {
      print "Selecting $state $function\n";
    }
    if ($assignments->{$state} eq $function) {
      my $weight = test_relevance($options, $segmentation, $function, $state));
      if ($weight > 0) {
        $states{$state} = $weight;
      }
    }
  }
  print "Selected\t$segmentation->{name}\t$function\t". join(" ", keys @states)."\n";
  return \%states;
}

sub test_relevance {
  my ($options, $segmentation, $function, $state) = @_;
  my $reference = "$options->{trackhub_dir}/overview/all_tfbs.bw"; 
  my $file = "$options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";

  for (my $i = 1; $i < scalar(@{$segmentation->{celltypes}}); $i++) {
    my $enrichment = compute_enrichment_between_files($reference, "gt $i $file");
    print "Enrichment\t$state\t$i:\t$enrichment\n";
    if ($enrichment == 0) {
      return 0;
    } elsif (compute_enrichment_between_files($reference, "gt $i $file") > $weak_cutoff) {
      return $i;
    }
  }

  return 0;
}

sub select_segmentation_cutoff {
  my ($options, $segmentation, $function) = @_;
  print "Setting cutoff for $function...\n";

  my $tfbs = "$options->{trackhub_dir}/overview/all_tfbs.bw";
  my $tfbs_auc = `wiggletools AUC - $tfbs`;

  my $celltype_count = scalar(@{$segmentation->{celltypes}});
  my $max_weight = undef;
  my $min_step = 1;
  my @files = ();
  my @weights = ();
  my %hash = %{$options->{selected_states}->{$segmentation->{name}}->{$function}};
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
  for (my $i = 0; $i < $celltype_count * $max_weight; $i += $min_step) {
    my $auc = `wiggletools AUC - mult $tfbs gt $i sum $files_string`;
    my $sens = $auc / $tfbs_auc; 

    my $overlap = `wiggletools AUC - gt 0 mult $tfbs gt $i sum $files_string`; 
    my $breadth = `wiggletools AUC - gt $i sum $files_string`;
    my $spec;
    if ($breadth == 0) {
      $spec = 0;
    } else {
      $spec = $overlap / $breadth; 
    }

    my $fscore;
    if ($sens + $spec == 0) {
      $fscore = 0;
    } else {
      $fscore = 2 * ($spec * $sens) / ($spec + $sens);
    }

    print "CutoffTest\t$function\t.\t$fscore\n";
    if ($fscore < $last_fscore) {
      print "cutoff set at $last_cutoff\n";
      return $last_cutoff;
    }
    $last_fscore = $fscore;
    $last_cutoff = $i;
  }

  print "cutoff set at MAX\n";
  return $celltype_count * $max_weight - $min_step;
}

########################################################
## Computing regions 
## This is the secret sauce:
## Foreach function ('tss', 'proximal', 'distal', 'ctcf'):
##   Foreach segmentation:
##     - Sum the summaries of all the states associated to
##     that function
##     - Select the regions where the sum is greater than
##     the cutoff (determined above)
##   Compute the union of these regions
##
## We then apply heuristic rules:
## * TSSs which have no experimental validation (CAGE)
## are downgraded to proximal
## * Proximal enhancers which overlap with a TSS are added
## to the TSS's whiskers (and do not exist separately)
## * TFBS+DNase regions which do not overlap any chromatin
## region are labelled 'tfbs'
## * DNAse regions which do not overlap any of the above 
## are labelled 'DNAse'
## 
## Params:
## * $options: hash ref
##   - $options->{trackhub_dir}
##   - $options->{working_dir}
##   - $options->{segmentations}: List of hashes with:
##     . $segmentation->{name}
##   - $options->{selected_states}->{$segmentation->{name}}->{$function}: List of strings
##   - $options->{cutoffs}->{$segmentation->{name}}->{$function}: Scalar
##   - $options->{tss} Path to BEd file with experimentally validated TSS
## Expected files:
## * $options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw 
## Writes into:
## * $options->{trackhub_dir}/overview/RegBuild.bb
########################################################

sub compute_regulatory_features {
  my $options = shift;
  my $bed_output = "$options->{trackhub_dir}/overview/RegBuild.bed";
  my $bb_output = "$options->{trackhub_dir}/overview/RegBuild.bb";

  mkdir "$options->{working_dir}/build/";

  if (!must_compute($options,$bb_output)) {
    print "$bb_output already exists, skipping...\n";
    return;
  } else {
    print "Computing Regulatory Build...\n";
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
  my $awk_dnase = make_awk_command("tfbs");
  run("wiggletools write_bg - unit $tfbs_signal | $awk_tfbs > $tfbs_tmp");

  #Compute open dnase
  my $dnase_tmp = "$options->{working_dir}/build/dnase.tmp.bed";
  my @files = glob "$options->{working_dir}/celltype_dnase/*.bed";
  my $awk_dnase = make_awk_command("dnase");
  run("wiggletools write_bg - unit sum ".join(" ", @files)." | $awk_dnase > $dnase_tmp");

  #############################################
  ## Percolate overlapping features up each level 
  #############################################

  # DNAse sites are merged to overlapping TFBS sites
  my $tfbs_tmp2 = "$options->{working_dir}/build/tfbs.tmp2.bed";
  expand_boundaries([$dnase_tmp], $tfbs_tmp, $tfbs_tmp2);

  # DNAse and TFBS sites are merged to overlapping distal sites
  my $distal_tmp2 = "$options->{working_dir}/build/distal.tmp2.bed";
  expand_boundaries([$dnase_tmp, $tfbs_tmp], $distal_tmp, $distal_tmp2);

  # DNAse, TFBS and distal sites are merged to overlapping proximal sites
  my $proximal_tmp2 = "$options->{working_dir}/build/proximal.tmp2.bed";
  expand_boundaries([$dnase_tmp, $tfbs_tmp, $distal_tmp], $proximal_tmp, $proximal_tmp2);

  # DNAse, TFBS, distal and proximal sites are merged to overlapping TSS sites
  my $tss_tmp2 = "$options->{working_dir}/build/tss.tmp2.bed";
  expand_boundaries([$dnase_tmp, $tfbs_tmp, $distal_tmp, $proximal_tmp], $tss_tmp, $tss_tmp2);

  #############################################
  ## Find non overlapping features for each level 
  #############################################

  # Convenience short hand
  my $awk_id = "awk '{\$4 = \$4\"_\"NR'; print;}";

  # All features that overlap known TSS are retained
  my $tss = "$output->{working_dir}/build/tss.bed";
  run("bedtools intersect -u -wa -a $tss_tmp2 -b $options->{tss} | $awk_id > $tss");

  # All features that do not go into a demoted file
  my $demoted = "$options->{working_dir}/build/demoted_tss.bed";
  run("bedtools intersect -v -wa -a $tss_tmp2 -b $options->{tss} > $demoted");

  # Unaligned proximal sites are retained
  my $proximal = "$options->{working_dir}/build/proximal.bed";
  run("bedtools unionbedg $proximal_tmp2 $demoted | bedtools intersect -wa -v -a stdin -b $tss_tmp2 | $awk_contract | $awk_id > $proximal");

  # Unaligned distal sites are retained
  my $distal = "$options->{working_dir}/build/distal.bed";
  run("bedtools intersect -wa -v -a $distal_tmp2 -b $tss_tmp2 | bedtools intersect -wa -v -a stdin -b $proximal_tmp2 | $awk_contract | $awk_id > $distal");

  # Unaligned TFBS sites are retained
  my $tfbs = "$options->{working_dir}/build/tfbs.bed";
  run("bedtools intersect -wa -v -a $tfbs_tmp2 -b $tss_tmp2 | bedtools intersect -wa -v -a stdin -b $proximal_tmp2 | bedtools intersect -wa -v -a stdin -b $distal_tmp2 | $awk_contract | $awk_id > $tfbs");

  # Unaligned DNAse sites are retained
  my $dnase = "$options->{working_dir}/build/dnase.bed";
  run("bedtools intersect -wa -v -a $dnase_tmp -b $tss_tmp2 | bedtools intersect -wa -v -a stdin -b $proximal_tmp2 | bedtools intersect -wa -v -a stdin -b $distal_tmp2 | bedtools intersect -wa -v -a stdin -b $tfbs | $awk_contract | $awk_id > $dnase");

  #############################################
  ## CTCF computed independently 
  #############################################

  my $awk_ctcf = make_awk_command("ctcf");
  my $ctcf = "$output->{working_dir}/build/ctcf.bed";
  run("$awk_ctcf $ctcf_tmp > $ctcf");

  #############################################
  ## Apply mask
  #############################################

  my $awk_mask = "";
  if (defined $options->{mask}) {
    $awk_mask = "| bedtools intersect -wa -v -a stdin -b $options->{mask}");
  }

  #############################################
  ## Merge
  #############################################
  run("sort -m $tss $proximal $distal $ctcf $dnase $tfbs -k1,1 -k2,2n $awk_mask > $bed_output");
  convert_to_bigBed($options, $bed_output);
}

sub expand_boundaries {
  my ($source_files, $target_file, $output);
  # Convenience short hand
  my $awk_expand = "awk 'BEGIN {OFS=\"\\t\"} {if (\$2 > 0) \$2 -= 1; \$3 += 1; print \$1,\$2,\$3,\$4;}'";
  my $awk_move_boundaries = "awk 'BEGIN {OFS=\"\\t\"} \$4 != name {if (name) {print chr, start, end, name, 1000, \".\", thickStart, thickEnd, rgb; } chr=\$1; start=\$2; end=\$3; name=\$4; thickStart=\$7; thickEnd=\$8; rgb=\$9} \$10 == chr && \$11+1 < start {start=\$11+1} \$10 == chr && \$12-1 > end {end=\$12-1} END {print chr, start, end, name, 1000, \".\", thickStart, thickEnd, rgb}'";

  run("cat ".join(" ", @{$source_files}). " | $awk_expand | bedtools intersect -loj -wa -wb -a $target_file -b stdin | $awk_move_boundaries > $output");
}

sub compute_initial_regions {
  my ($options, $function) = @_;
  my $output = "$output->{working_dir}/build/$function.tmp.bed";
  my $wiggletools_cmd = "wiggletools write_bg - unit sum ";
  foreach my $segmentation (@{$options->{segmentations}}) {
    $wiggletools_cmd .= function_definition($options, $function, $segmentation); 
  }
  my $awk_cmd = make_awk_command($function);
  run("$wiggletools_cmd | $awk_cmd > $output") ;
  return $output;
}

sub function_definition {
  my ($options, $function, $segmentation) = @_;
  my $hash = $options->{selected_states}->{$segmentation->{name}}->{$function};
  my $cutoff = $options->{cutoffs}->{$segmentation->{name}}->{$function};
  my @summaries = ();
  foreach my $state (keys %$hash) {
    push @summaries, "scale ".(1/$hash->{$state})." options->{trackhub_dir}/segmentation_summaries/$segmentation->{name}/$state.bw";
  }
  if (scalar @summaries > 0) {
    return " gt $cutoff sum " . join(" ", @summaries);
  } else {
    return "";
  }
}

sub make_awk_command {
  my ($function) = @_;
  # This awk one-liner takes in a BedGraph and converts into a Bed9
  return "awk \'BEGIN {OFS=\"\\t\"} {\$4=\"${function}_\"NR; \$5=1000; \$6=\".\"; \$7=\$2; \$8=\$3; \$9=\"$COLORS{$function}\"; print }\'";
}

########################################################
## Computing cell-type specific state
## Given the Regulatory Build, you want to know which 
## Feature is active in which cell type
## Params:
## * $options: hash containing:
##   - $options->{trackhub_dir}
##   - $options->{working_dir}
##   - $options->{segmentation}
## Expected files:
## * $options->{working_dir}/build/$annotation.bed
## * $options->{working_dir}/celltype_tf/$celltype.bed
## * $options->{working_dir}/celltype_dnase/$celltype.bed
## Created files:
## * $options->{trackhub_dir}/projected_segmentations/$celltype.bed
########################################################

sub compute_states {
  my $options = shift;
  my $segmentation;

  mkdir "$options->{trackhub_dir}/projected_segmentations/";
  mkdir "$options->{working_dir}/projected_segmentations/";

  foreach $segmentation (@{$options->{segmentations}}) {
    compute_segmentation_states($options, $segmentation);
  }
}

sub compute_segmentation_states {
  my ($options, $segmentation) = @_;
  my $celltype;

  mkdir "$options->{working_dir}/projected_segmentations/$segmentation->{name}/";

  foreach $celltype (@{$segmentation->{celltypes}}) {
    if (must_compute($options,"$options->{trackhub_dir}/projected_segmentations/$celltype.bb")) {
      compute_celltype_state($options, $segmentation, $celltype);
    } else {
      print "Projected segmentation of $celltype already computed, skipping...\n";
    }
  }
}

sub compute_celltype_state {
  my ($options, $segmentation, $celltype) = @_;

  mkdir "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype";

  if ($segmentation->{type} eq 'ChromHMM' || $segmentation->{type} eq 'Segway') {
    compute_ChromHMM_celltype_state($options, $segmentation, $celltype);
  } else {
    die ("Unknown segmentation type $segmentation->{type}");
  }
}

sub compute_ChromHMM_celltype_state {
  my ($options, $segmentation, $celltype) = @_;
  my $annotation;

  my @bedfiles = ();
  
  foreach $annotation (@annotations) {
    push @bedfiles, compute_ChromHMM_annotation_state($options, $segmentation, $celltype, $annotation);
  }

  my $output = "$options->{trackhub_dir}/projected_segmentations/$celltype.bed";
  run("sort -m " . join(" ", @bedfiles) . " -k1,1 -k2,2n > $output") ;
  convert_to_bigBed($options, $output);
}

sub compute_ChromHMM_annotation_state {
  my ($options, $segmentation, $celltype, $annotation) = @_;
  my $temp;
  my $process=1;

  if ($annotation eq 'tfbs') {
    $temp = "$options->{working_dir}/celltype_tf/$celltype.bed";
  } elsif ($annotation eq 'dnase') {
    $temp = "$options->{working_dir}/celltype_dnase/$celltype.bed";
  } else {
    my @files = ();
    foreach my $state (@{$segmentation->{states}}) {
      if ($options->{assignments}->{$segmentation->{name}}->{$state} eq $annotation) {
        push @files, "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$state/$celltype.bed";
      }
    }

    if (scalar @files == 0) {
      $process=0;
    } else {
      $temp = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.bed";
      run("sort -m " . join(" ", @files). " -k1,1 -k2,2n > $temp");
    }
  }

  if ($process && -e $temp) {
    my $reference = "$options->{working_dir}/build/$annotation.bed";
    my $repressed = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$celltype/repressed.bed";
    my $poised = "$options->{working_dir}/segmentation_summaries/$segmentation->{name}/$celltype/poised.bed";

    # POISED => POISED
    my $temp2 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.poised.bed";
    my $poised_color = $COLORS{poised};
    run("bedtools intersect -wa -u -a $reference -b $poised | awk \'{\$9=\"$poised_color\"; print}\' > $temp2");

    # !POISED && ACTIVE && REPRESSED => POISED
    my $temp3 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.active_repressed.bed";
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -u -a stdin -b $temp | bedtools intersect -wa -u -a stdin -b $repressed | awk \'{\$9=\"$poised_color\"; print}\' > $temp3");

    # !POISED && !ACTIVE && REPRESSED => REPRESSED
    my $temp4 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.repressed.bed";
    my $repressed_color = $COLORS{repressed};
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -v -a stdin -b $temp | bedtools intersect -wa -u -a stdin -b $repressed | awk \'{\$9=\"$repressed_color\"; print}\' > $temp4");

    # !POISED && ACTIVE && !REPRESSED => ACTIVE
    my $temp5 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.active.bed";
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -u -a stdin -b $temp | bedtools intersect -wa -v -a stdin -b $repressed > $temp5");

    # !POISED && !ACTIVE && !REPRESSED => DEAD
    my $temp6 = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.inactive.bed";
    my $dead_color = $COLORS{dead};
    run("bedtools intersect -wa -v -a $reference -b $poised | bedtools intersect -wa -v -a stdin -b $temp | bedtools intersect -wa -v -a stdin -b $repressed | awk \'{\$9=\"$dead_color\"; print}\' > $temp4");

    # Merge all this into one bed file
    my $output = "$options->{working_dir}/projected_segmentations/$segmentation->{name}/$celltype/$annotation.final.bed";
    run("sort -m $temp2 $temp3 $temp4 $temp5 $temp6 -k1,1 -k2,2n > $output");
    return $output;
  } else {
    # Could not find evidence for or against, do not report any regions
    return "";
  }
}

########################################################
## Track hub creation
##
## Params:
## * $options: hash containing:
##   - $options->{url} Base URL of track hub
##   - $options->{assembly} UCSC name of the assembly
##   - $options->{trackhub_dir}
##   - $options->{segmentations} List of hashes containing:
##     . $segmentation->{name}
##     . $segmentation->{states}
##     . $segmentation->{celltypes}
##
## Expects files:
## * $options->{trackhub_dir}/tfbs/$tf.bw
##
## Writes into:
## * hub.txt
## * genomes.txt
## * $options->{trackhub_dir}/trackDb.txt    
########################################################

sub make_track_hub {
  my $options = shift;
  print "Making track hub header files\n";
  make_track_hub_headers($options);
  make_track_hub_assembly($options);
}

sub make_track_hub_headers {
  my $options = shift;
  my $fh;

  open($fh, ">", "$options->{output_dir}/genomes.txt");
  print $fh "genome $options->{assembly}\n";
  print $fh "trackDb $options->{assembly}/trackDb.txt\n";
  close($fh);

  open($fh, ">", "$options->{output_dir}/hub.txt");
  print $fh "hub EnsemblRegulatoryBuild\n";
  print $fh "shortLabel Beta Regulatory Build\n";
  print $fh "longLabel Evidence summaries and provisional results for the new Ensembl Regulatory Build\n";
  print $fh "genomesFile genomes.txt\n";
  print $fh "email helpdesk\@ensembl.org\n";
  close($fh);
}

sub make_track_hub_assembly {
  my $options = shift;

  my $output = "$options->{trackhub_dir}/trackDb.txt";
  open(my $file, ">", $output) || die("Could not open $output\n");

  make_track_hub_overview($options, $file);
  make_track_hub_segmentations($options, $file);
  make_track_hub_segmentation_summaries($options, $file);
  make_track_hub_projected_segmentations($options, $file);
  make_track_hub_tfbs($options, $file);

  close $file;
}

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
  print $file "\tbigDataUrl $options->{url}/$options->{assembly}/overview/RegBuild.bb\n";
  print $file "\tshortLabel Regulatory Build\n";
  print $file "\tlongLabel Ensembl Regulatory annotation of regional function\n";
  print $file "\ttype bigBed 9\n";
  print $file "\titemRgb on\n";
  print $file "\tsubGroups source=EnsRegBuild\n";
  print $file "\tpriority 1\n";
  print $file "\tvisibility squish\n";

  print $file "\n";

  print $file "\ttrack TFBSSummary\n";
  print $file "\tparent BuildOverview on\n";
  print $file "\tbigDataUrl $options->{url}/$options->{assembly}/overview/all_tfbs.bw\n";
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

sub make_track_hub_segmentations {
  my ($options, $file) = @_;
  my $segmentation;
  foreach $segmentation (@{$options->{segmentations}}) {
    make_track_hub_segmentations_2($options, $file, $segmentation);
  }
}

sub make_track_hub_segmentations_2 {
  my ($options, $file, $segmentation) = @_;

  my $name = $segmentation->{name};
  my @states = $segmentation->{states};
  my $state_count = scalar(@states);
  my $celltypes = $segmentation->{celltypes};
  my $celltype_count = scalar(@$celltypes);
  my $state_string = "";

  my $celltype_list = "";
  my $celltype;
  foreach $celltype (@{$segmentation->{celltypes}}) {
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

  foreach $celltype (@{$segmentation->{celltypes}}) {
    my $short = $celltype;
    $short =~ s/[_-]//g;

    print $file "\n";

    print $file "\ttrack ${name}_${celltype}_segments\n";
    print $file "\tbigDataUrl $options->{url}/$options->{assembly}/segmentations/$segmentation->{name}/$celltype.bb\n";
    print $file "\tparent ${name}CellSegments off\n";
    print $file "\tshortLabel Seg. ${celltype}\n";
    print $file "\tlongLabel $segmentation->{type} segmentation of the $celltype genome after training on $celltype_count cell types\n";
    print $file "\ttype bigBed 9\n";
    print $file "\titemRgb on\n";
    print $file "\tsubGroups cellType=$short source=ERB\n";
  }
}

sub make_track_hub_segmentation_summaries {
  my ($options, $file) = @_;
  my $segmentation;

  foreach $segmentation (@{$options->{segmentations}}) {
    make_track_hub_segmentation_summaries_2($options, $file, $segmentation);
  }
}

sub make_track_hub_segmentation_summaries_2 {
  my ($options, $file, $segmentation) = @_;

  my $name = ${segmentation}->{name};
  my $states = $segmentation->{states};
  my $state_count = scalar(@$states);
  my $celltypes = $segmentation->{celltypes};
  my $celltype_count = scalar(@$celltypes);
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
    print $file "\tbigDataUrl $options->{url}/$options->{assembly}/segmentation_summaries/$segmentation->{name}/$state.bw\n";
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

sub make_track_hub_projected_segmentations {
  my ($options, $file) = @_;
  my $celltype;

  my $celltype_list = "";
  foreach my $segmentation (@{$options->{segmentations}}) {
    foreach $celltype (@{$segmentation->{celltypes}}) {
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
    foreach $celltype (@{$segmentation->{celltypes}}) {
      my $short = $celltype;
      $short =~ s/[_-]//g;
      $celltype_list .= " $short=$celltype";

      print $file "\n";

      print $file "\ttrack ${name}_${celltype}_projected\n";
      print $file "\tbigDataUrl $options->{url}/$options->{assembly}/projected_segmentations/$celltype.bb\n";
      print $file "\tparent ProjectedSegments off\n";
      print $file "\tshortLabel $celltype\n";
      print $file "\tlongLabel Projection of the Ensembl Regulatory build onto $celltype\n";
      print $file "\ttype bigBed 9\n";
      print $file "\titemRgb on\n";
      print $file "\tsubGroups cellType=$short source=ERB\n";
    }
  }
}

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
    print $file "\tbigDataUrl $options->{url}/$options->{assembly}/tfbs/$tf.bw\n";
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
