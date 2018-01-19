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

load_build.pl

=head1 SYNOPSIS

perl load_build.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38 --base_dir hg38/

=head1 DESCRIPTION

Loads a segmentation BigBed file annotated by the new regulatory build into the database.
Params:
  * base_dir: directory with assembly name (e.g. ./hg38) create by build
  * your usual Ensembl Funcgen MySQL params: host, user etc...

In short, the file $base_dir/segmentations/$segmentation/$cell_type_name.bb must exist

Also bigBedToBed must be on the commandline.

This script is pretty database intensive (think inserting 4Mo entries) so keep an eye
when running in parallel. In the good times, each job takes ~1h to run.

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Data::Dumper qw (Dumper);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Analysis;
use File::Temp qw/ tempfile tempdir /;
use Hash::Util qw( lock_hash );

${File::Temp::KEEP_ALL} = 1;

# our %rgb_state = (
#   '225,225,225' => 0, # dead
#   "192,0,190" => 2, # poised
#   "127,127,127" => 3, # repressed
#   "255,255,255" => 4, # NA
#
#   "255,0,0" => 1, # TSS
#   "209,157,0" => 1, # TFBS
#   "255,252,4" => 1, # DNase
#   "255,105,105" => 1, # Proximal
#   "250,202,0" => 1, # Distal
#   "10,190,254" => 1, # CTCF
# );
my %rgb_state = (
  '225,225,225' => 'INACTIVE',

  '255,0,0'     => 'ACTIVE', # TSS
  '209,157,0'   => 'ACTIVE', # TFBS
  '255,252,4'   => 'ACTIVE', # DNase
  '255,105,105' => 'ACTIVE', # Proximal
  '250,202,0'   => 'ACTIVE', # Distal
  '10,190,254'  => 'ACTIVE', # CTCF

  '192,0,190'   => 'POISED',
  '127,127,127' => 'REPRESSED',
  '255,255,255' => 'NA',
);
lock_hash(%rgb_state);

my %label_description= (
  'ctcf'     => 'CTCF Binding Site',
  'distal'   => 'Predicted enhancer',
  'proximal' => 'Predicted promoter flanking region',
  'tss'      => 'Predicted promoter',
  'tfbs'     => 'Transcription factor binding site',
  'dnase'    => 'Open chromatin region'
);
lock_hash(%label_description);

my $start_time = time;

main();

=head2 main

  Description: Overall process
  Returntype: undef

=cut

sub main {
  print_log("Getting options\n");
  my $options = get_options();
  print_log("Connecting to database\n");
  my $funcgen_db_adaptor = connect_db($options);

  use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor;
  my $regulatory_build_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor->new($funcgen_db_adaptor);
  my $current_regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

  if (! defined $current_regulatory_build) {
#     die("Couldn't find regulatory build in the database!");
  }
  if (defined $current_regulatory_build) {
    print "Found regulatory build: "
    . $current_regulatory_build->name
    . " " . $current_regulatory_build->version
    . " from  "
    . $current_regulatory_build->initial_release_date
    . " in the database.\n";
  }

#   print_log("Getting analysis\n");
#   my $analysis = get_analysis($funcgen_db_adaptor);
  print_log("Getting cell types\n");
  my $ctypes = get_cell_type_names($options->{base_dir}, $funcgen_db_adaptor);
  print_log("Getting stable ids\n");
  my $stable_id = get_stable_id($options, $funcgen_db_adaptor);
  print_log("Getting slices\n");
  my $slice = get_slices($funcgen_db_adaptor);
  print_log("Getting feature types\n");
  my $feature_type = get_feature_types($funcgen_db_adaptor);
  print_log("Counting active features\n");
  my $count_hash = compute_counts($options->{base_dir});

  print_log("Creating regulatory build object\n");

  my $is_small_update = defined $options->{small_update};

  my $new_regulatory_build = create_regulatory_build_object(
    $current_regulatory_build,
    $is_small_update,
    $funcgen_db_adaptor
  );

  # This sets the dbID of the regulatory build object. The regulatory
  # features are linked to that.
  #
  $regulatory_build_adaptor->store($new_regulatory_build);

  print_log("Creating regulatory_feature table\n");
  compute_regulatory_features($options, $ctypes, $feature_type, $stable_id, $count_hash, $slice, $funcgen_db_adaptor, $new_regulatory_build);

  print_log("Creating regulatory_annotation table\n");
  compute_regulatory_annotations($options);
  print_log("Updating meta table\n");

  if (defined $current_regulatory_build) {
    $current_regulatory_build->is_current(0);
    $regulatory_build_adaptor->update($current_regulatory_build);
  }

  $new_regulatory_build->is_current(1);
  $regulatory_build_adaptor->update($new_regulatory_build);
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

=head2 get_options

  Description: Command line options
  Returntype: hashreof

=cut

sub get_options {
  my %options = ();
  use Pod::Usage;
  GetOptions (
    \%options,
    "base_dir|b=s",
    "pass|p=s",
    "port=s",
    "host|h=s",
    "user|u=s",
    "dbname|d=s",
    "small_update",
  ) or pod2usage( -exitval => 1);
  defined $options{base_dir} || die ("You must define the base directory!\t--base_dir XXXX\n");
  defined $options{host} || die ("You must define the destination host!\t--host XXXX\n");
  defined $options{user} || die ("You must define the user login!\t--user XXXX\n");
  defined $options{dbname} || die ("You must define the database name!\t--dbname XXXX\n");
  return \%options;
}

=head2 connect_db

  Description: Connecting to the DB
  Arg1: options hash ref
  Returntype: Bio::EnsEMBL::Funcgen::DBAdaptor object

=cut

sub connect_db {
  my ($options) = @_;
  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $options->{host},
     -user   => $options->{user},
     -dbname => $options->{dbname},
     -pass   => $options->{pass},
     -port   => $options->{port},
     -dnadb_host => $options->{dnadb_host},
     -dnadb_port => $options->{dnadb_port},
     -dnadb_user => $options->{dnadb_user},
     -dnadb_pass => $options->{dnadb_pass},
     -dnadb_name => $options->{dnadb_name},
    );

  #Test connections
  $db->dbc->db_handle;
  $db->dnadb->dbc->db_handle;
  if(! defined $db->species){
    die("Could not get a valid species from $options->{dbname}, please check the meta species.production_name");
  }

  return $db;
}

=head2 get_analysis

  Description: Create/Get analysis for Build
  Arg1: Bio::EnsEMBL::Funcgen::DBAdaptor
  Returntype: Bio::EnsEMBL::Funcgen::Analysis object

=cut

sub get_analysis {
  my ($db) = @_;
  my $aa   = $db->get_AnalysisAdaptor();
  my $ana = $aa->fetch_by_logic_name('Regulatory_Build');

  if ( not defined $ana ) {
    my $analysis = Bio::EnsEMBL::Analysis->new
	(
	 -logic_name      => 'Regulatory_Build',
	 -db              => undef,
	 -db_version      => undef,
	 -db_file         => undef,
	 -program         => undef,
	 -program_version => undef,
	 -program_file    => undef, #Could use $program_name here, but this is only part of the build
	 -gff_source      => undef,
	 -gff_feature     => undef,
	 -module          => undef,
	 -module_version  => undef,
	 -parameters      => undef,
	 -created         => undef,
	 -description     => q({'reg_feats' => 'Features from <a href="/info/genome/funcgen/index.html" class="cp-external">Ensembl Regulatory Build</a>.','core' => 'Sites enriched for marks of open chromatin (e.g. DNase1) or transcription factor binding sites.','non_core' => 'Sites enriched for histone modifications or polymerase binding.'}),
	 -display_label   => 'Regulatory Build',
	 -displayable     => 1,
	 -web_data        => q({'type' => 'fg_regulatory_features', 'name' => 'Reg. Feats', 'display' =>'off', 'depth' => 10, 'default' => {'contigviewbottom' => 'normal', 'generegview' => 'normal'} }),
	);
    $aa->store($analysis);
    return $analysis;
  } else {
    return $ana;
  }
}

=head2 clean_name

  Description: Removing unwanted characters
    Quick string normalisation function tor remove weird
    characters froms file names and remove variants
  Arg1: String
  Returntype: String

=cut

sub clean_name {
  my $string = shift;
  $string =~ s/[\-\(\)]//g;
#   $string =~ s/_.*//g;
  $string = uc($string);
  $string =~ s/:/x/g;
  return $string;
}

=head2 get_cell_type_names

  Description: Get list of cell types used in the regulatory build.
  Arg1: The base_directory name
  Arg2: Bio::EnsEMBL::Funcgen::DBAdaptor object
  Returntype: array ref

=cut

sub get_cell_type_names {
  my ($base_dir, $db) = @_;

  my $epigenome_adaptor = $db->get_EpigenomeAdaptor();
  my %cell_type_from_clean = ();
  for my $epigenome (@{$epigenome_adaptor->fetch_all()}) {
    $cell_type_from_clean{clean_name($epigenome->production_name)} = $epigenome;
    defined $epigenome || die("Unrecognized cell type name $epigenome\n");
  }

#   my $debug_max = 1;

  my @cell_types = ();
  foreach my $file (glob "$base_dir/projected_segmentations/*.bb") {
    my $cell_type_name = basename $file;
#     $cell_type_name =~ s/\.bb//;
     $cell_type_name =~ s/_\d+_SEGMENTS\.bb//;
    if (!exists $cell_type_from_clean{$cell_type_name}) {

      my @known_cell_types = sort keys %cell_type_from_clean;
      my $known_cell_types_as_string = join "\n", map { ' - ' . $_ } @known_cell_types;

      die(
	"Celltype $cell_type_name has a projected segmentation in $file, but is unknown! Known cell types are:\n$known_cell_types_as_string"

      );
    }
    push @cell_types, $cell_type_from_clean{$cell_type_name};
#     last if (@cell_types == $debug_max);
  }
  return \@cell_types;
}

=head2 run

  Description: Convenience wrapper to run the commandline safely
  Arg1: command line command string
  Returntype: undef
  Side effects: runs command

=cut

sub run {
  my ($cmd) = @_;
  print_log($cmd . "\n");
  use Carp;
  system($cmd) && confess("Failed when running command:\n$cmd\n");
}

=head2

  Description: Assign stable ids to features
  Arg1: options hash ref
  Arg2: Bio::EnsEMBL::Funcgen::DBAdaptor
  Returntype: hashref

=cut

sub get_stable_id {
  my ($options, $db) = @_;

  my ($ofh, $old) = tempfile();
  sub cmp_features {
    if ($a->seq_region_name ne $b->seq_region_name) {
      return $a->seq_region_name cmp $b->seq_region_name;
    } else {
      return $a->seq_region_start <=> $b->seq_region_start;
    }
  }
  my $feature_set = $db->get_adaptor('FeatureSet')->fetch_by_name('RegulatoryFeatures:MultiCell_v'.$options->{old_version});
  my ($overlaps, $max_id);
  my ($fh, $new) = tempfile();
  close $fh;
  run("bigBedToBed $options->{base_dir}/overview/RegBuild.bb $new");

  if (defined $feature_set) {
    foreach my $feature (sort cmp_features @{$feature_set->get_all_Features()}) {
      print $ofh join("\t", ($feature->seq_region_name, $feature->bound_start, $feature->bound_end, $feature->feature_type->name, substr($feature->stable_id, 4))) . "\n";
    }
    close $ofh;

    ($overlaps, $max_id) = get_overlaps_between_files($old, $new);
  } else {
    ($overlaps, $max_id) = ([], 0);
  }

# Go through overlaps in order of increasing overlap length. This means that you should always
# overwrite an overlap with a later one.
  sub cmp_overlaps {
    return $a->[14] <=> $b->[14];
  }
  my %old_pref = ();
  my %new_pref = ();
  foreach my $entry (sort cmp_overlaps @{$overlaps}) {
    $old_pref{$entry->[4]} = $entry->[8];
    $new_pref{$entry->[8]} = $entry->[4];
  }

# Output new_id/old_id pairs if mutual best hits, else create new ids starting from
# the maximum of pre-existing IDs + 1
  my %stable_id_hash = ();
  my $next_free_id = $max_id + 1;
  open my $in, "<", $new;
  while (my $line = <$in>) {
    chomp $line;
    my @items = split /\t/, $line;
    my $new_id = $items[3];

    my $stable_id = undef;

    if (exists $new_pref{$new_id} && $old_pref{$new_pref{$new_id}} eq $new_id) {
      $stable_id = $new_pref{$new_id};
    } else {
      $stable_id = $next_free_id;
      $next_free_id += 1;
    }

    # Creating stable id string, composed of prefix + 11 digit integer, front padded with 0s
    my $species = $db->get_MetaContainer->get_production_name;
    if ($species eq 'homo_sapiens') {
      $stable_id_hash{$new_id} = "ENSR" . sprintf("%011d", $stable_id);
    } elsif ($species eq 'mus_musculus') {
      $stable_id_hash{$new_id} = "ENSMUSR" . sprintf("%011d", $stable_id);
    } else {
      # The general strategy (excluding human and mouse) is first letter of genus followed
      # by first two letters of species
      my @components = split('_', $species);
      my $triletter_code = substr($components[0], 0, 1) . substr($components[1], 0, 2);
      $stable_id_hash{$new_id} = "ENS" . uc($triletter_code) . 'R' . sprintf("%011d", $stable_id);
    }
  }
  close $in;
  unlink $new;
  unlink $old;

  return \%stable_id_hash;
}

=head2 get_overlaps_between_files

  Description: Computes overlaps between two builds contained in bed files
  Arg1: old build file location
  Arg2: new build file location
  Returntype: list ref containing: array ref of overlaps, maximum of the old build's stable ids
'
=cut

sub get_overlaps_between_files {
  my ($old, $new) = @_;

# Compute overlaps between regions defined in both file
  my ($fh, $filename) = tempfile();
  run("bedtools intersect -a $old -b $new -wo > $filename");

# Parse output
  my @overlaps = ();
  my $max_id = undef;
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    # Check for feature type incompatibility
    my @comps = split("_", $items[8]);
    my $label = $comps[0];
    if ($items[3] ne $label_description{$label}) {
      next;
    }
    push @overlaps, \@items;

    # Look for maximum ID number among the old features
    my $id = $items[4];
    if (!defined $max_id || $max_id < $id) {
      $max_id = $id;
    }
  }
  unlink $filename;

  return (\@overlaps, $max_id);
}

=head2 get_slices

  Description: Get slice for each chromosome
  Arg1: Bio::EnsEMBL::Funcgen::DBAdaptor object
  Returntype: hashref:
    - seq_region_name => Bio::EnsEMBL::Slice object

=cut

sub get_slices {
  my ($db) = @_;
  my $slices = $db->get_adaptor("Slice")->fetch_all('toplevel',undef,0,1);
  my %hash = ();
  foreach my $slice (@{$slices}) {
    $hash{$slice->seq_region_name} = $slice;
  }
  return \%hash;
}

=head2 compute_counts

  Description: Count the number of active features for each temporary
  ID across all cell types
  Arg1: Base directory
  Returntype: temp id => scalar

=cut

sub compute_counts {
  my ($base_dir) = @_;
  my $count_hash = {};
  foreach my $file (glob "$base_dir/projected_segmentations/*.bb") {
    count_active($file , $count_hash);
  }
  return $count_hash;
}

=head2 count_active

  Description: Count the number of active features for each temporary
  ID across one cell types
  Arg1: BigBed file location
  Arg2: Hashref
  Returntype: undef
  Side effects: updates hash ref

=cut

sub count_active {
  my ($filename, $count_hash) = @_;
  print_log("\tCounting in file $filename\n");
  my ($fh, $tmp_name) = tempfile();
  run("bigBedToBed $filename $tmp_name");
  while (my $line  = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;
    if (!defined $count_hash->{$name}) {
      $count_hash->{$name} = 0;
    }
    if ($rgb_state{$rgb} eq 'ACTIVE') {
      $count_hash->{$name} += 1;
    }
  }
  close $fh;
  unlink $tmp_name;
}

=head2 get_feature_types

  Description: Get a hashref from feature type name to FeatureType
   object
  Arg1: Bio::EnsEMBL::Funcgen::DBAdaptor object
  Returntype: Bio::EnsEMBL::Funcgen::FeatureType object

=cut

sub get_feature_types {
  my ($db) = @_;
  my $fta = $db->get_adaptor("FeatureType");
  my %long_name= (
    'ctcf'=>'CTCF Binding Site',
    'distal'=>'Enhancer',
    'proximal'=>'Promoter Flanking Region',
    'tss'=>'Promoter',
    'tfbs'=>'TF binding site',
    'dnase'=>'Open chromatin'
  );
  my %so_accession = (
    'ctcf'=>'SO:0001974',
    'distal'=>'SO:0000165',
    'proximal'=>'SO:0001952',
    'tss'=>'SO:0000167',
    'tfbs'=>'SO:0000235',
    'dnase'=>'SO:0001747'
  );
  my %so_name = (
    'ctcf'=>'CTCF_binding_site',
    'distal'=>'enhancer',
    'proximal'=>'promoter_flanking_region',
    'tss'=>'promoter',
    'tfbs'=>'TF_binding_site',
    'dnase'=>'open_chromatin_region'
  );

  my $feature_type = {};
  for my $label (('ctcf','distal','proximal','tss','tfbs','dnase')) {
    $feature_type->{$label} = $fta->fetch_by_name($long_name{$label});

    if (! defined $feature_type->{$label}) {
      my $ft = Bio::EnsEMBL::Funcgen::FeatureType->new(
	-name => $long_name{$label},
	-class => 'Regulatory Feature',
	-description => $label_description{$label},
	-analysis => undef,
	-so_name => $so_name{$label},
	-so_accession => $so_accession{$label}
      );
      $fta->store($ft);
      $feature_type->{$label} = $ft;
    }
  }
  return $feature_type;
}

=head2 compute_regulatory_features

  Description: Creates the actual RegulatoryFeature objects
  Arg1: Options hashref
  Arg2: Array ref: cell types
  Arg3: Hash ref: label -> FeatureType
  Arg3: Hash ref: temporary id -> stable id
  Arg4: Hash ref: temporary id -> count
  Arg5: hash ref: chromosome name -> seq_region_id
  Arg6: Bio::EnsEMBL::Funcgen::DBAdaptor object
  Arg7: Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Returntype: undef
  Side effects: writes into regulatory_feature table

=cut

sub compute_regulatory_features {
  my ($options, $cell_type, $feature_type, $stable_id, $count_hash, $slice, $db, $new_regulatory_build) = @_;

  my $regulatory_feature_adaptor = $db->get_adaptor("RegulatoryFeature");
   
  my $regulatory_features = load_regulatory_build($options->{base_dir}, $stable_id, $count_hash, $slice, $feature_type, $regulatory_feature_adaptor, $new_regulatory_build);
 
  print "Going through ".(scalar @$cell_type)." cell types\n";
  
  foreach my $current_cell_type (@$cell_type) {
    load_celltype_activity($options->{base_dir}, $current_cell_type, $regulatory_features);
  }

  $regulatory_feature_adaptor->store(values %$regulatory_features);
}

=head2 load_regulatory_build

  Description: loads the data from the build's BigBed files into the database
  Arg1: Location of base directory
  Arg2: hashref: bedfile id => stable id
  Arg3: hashref: bedfile id => integer
  Arg4: hashref: slice name => Bio::EnsEMBL::Slice
  Arg5: hashref: label => Bio::EnsEMBL::Funcgen::FeatureType
  Arg6: Bio::EnsEMBL::Funcgen::RegulatoryFeatureAdaptor
  Arg7: Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Returntype: hashref bedfile id => Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Side effects: writes into regulatory_feature table
'
=cut

sub load_regulatory_build {

  my ($base_dir, $stable_id, $count_hash, $slice, $feature_type, $regulatory_feature_adaptor, $new_regulatory_build) = @_;
  
  print_log("\tLoading Regulatory Features\n");

  my ($tmp, $tmp_name) = tempfile();
  my $bigbed = "$base_dir/overview/RegBuild.bb";
  run("bigBedToBed $bigbed $tmp_name");

  my $regulatory_features = process_regulatory_build_file($tmp, $stable_id, $count_hash, $slice, $feature_type, $regulatory_feature_adaptor, $new_regulatory_build);

  close $tmp;
  unlink $tmp_name;
  return $regulatory_features;
}


=head2 process_regulatory_build_file

  Description: loads the data from a Bed file into the database
  Arg1: filehandle into input file
  Arg3: hashref: bedfile id => stable id
  Arg4: hashref: bedfile id => integer
  Arg5: hashref: slice name => Bio::EnsEMBL::Slice
  Arg6: hashref: label => Bio::EnsEMBL::Funcgen::FeatureType
  Arg7: Bio::EnsEMBL::Funcgen::RegulatoryFeatureAdaptor
  Arg8: Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Returntype: hashref bedfile id => Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Side effects: writes into regulatory_feature table

=cut

sub process_regulatory_build_file {
  my ($fh, $stable_id, $count_hash, $slice, $feature_type, $rfa, $new_regulatory_build) = @_;

  my $regulatory_features = {};

  while (my $line = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split "\t", $line;
    my ($feature_type_str, $number) = split '_', $name;

    exists $feature_type->{$feature_type_str} || die("Could not find feature type for $feature_type_str\n".join("\t", keys %{$feature_type})."\n");
    exists $slice->{$chrom} || die("Could not find slice type for $chrom\n".join("\t", keys %{$slice})."\n");
    exists $stable_id->{$name} || die("Could not find stable ID for feature # $name\n");

    exists $count_hash->{$name} || die("Could not find count for feature # $name\n");

    my $regulatory_feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast({
      slice               => $slice->{$chrom},
      start               => $thickStart + 1,
      end                 => $thickEnd,
      strand              => 0,
      feature_type        => $feature_type->{$feature_type_str},
      _bound_lengths      => [$thickStart - $start, $end - $thickEnd],
      epigenome_count     => $count_hash->{$name},
      stable_id           => $stable_id->{$name},
#       analysis            => $analysis,
      regulatory_build_id => $new_regulatory_build->dbID,
    });

    $regulatory_features->{$name} = $regulatory_feature;
  }

  return $regulatory_features;
}

=head2 load_celltype_activity

  Description: loads the data from the build's BigBed files into the database
  Arg1: Location of base directory
  Arg2: Bio::EnsEMBL::Funcgen::Celltype
  Arg3: hashref: bedfile id => Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Returntype: undef
  Side effects: writes into regulatory_feature table
'
=cut

sub load_celltype_activity {
  my ($base_dir, $cell_type, $regulatory_features) = @_;

  print_log("\tProcessing data from cell type " . $cell_type->display_label . " (". $cell_type->name .")" . "\n");

  my $cell_type_name = clean_name($cell_type->production_name);
#   my $bigbed = "$base_dir/projected_segmentations/$cell_type_name.bb";
  my $bigbed = "$base_dir/projected_segmentations/${cell_type_name}_25_SEGMENTS.bb";
  my ($tmp, $tmp_name) = tempfile();
  run("bigBedToBed $bigbed $tmp_name");

  process_celltype_file($tmp, $cell_type, $regulatory_features);

  close $tmp;
  unlink $tmp_name;
}

=head2 process_celltype_file

  Description: loads the data from a Bed file into the database
  Arg1: filehandle into input file
  Arg2: Bio::EnsEMBL::Funcgen::Celltype
  Arg3: hashref: bedfile id => Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Returntype: undef
  Side effects: writes into regulatory_feature table

=cut

sub process_celltype_file {
  my ($fh, $cell_type, $regulatory_feature) = @_;
  while (my $line = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split "\t", $line;

    use Bio::EnsEMBL::Funcgen::RegulatoryActivity;
    my $regulatory_activity = Bio::EnsEMBL::Funcgen::RegulatoryActivity->new;
    $regulatory_activity->activity($rgb_state{$rgb});

    $regulatory_activity->_epigenome_id($cell_type->dbID);
    
    exists $regulatory_feature->{$name} || die("Could not find Regulatory Feature for $name\n");
    $regulatory_feature->{$name}->add_regulatory_activity($regulatory_activity);
  }
}

=head2 compute_regulatory_annotations

  Description: Assign motifs and annotations to regulatory features
  Arg1: options hash ref
  Returntype: undef
  Side effects: writes into regulatory_evidences table

=cut

sub compute_regulatory_annotations {
  my $options = shift;

  my ($tmp_fh, $cell_type_regulatory_features) = tempfile();
  my ($tmp_fh2, $annotations) = tempfile();
  my ($tmp_fh3, $motifs) = tempfile();
  my ($tmp_fh4, $out) = tempfile();
  close $tmp_fh;
  close $tmp_fh2;
  close $tmp_fh3;
  close $tmp_fh4;

  # Extract regulatory, annotated, and motif features into flat tab-delimited files.
  # The common format of these files is:
  # chrom	start	end	ID
  
  # "convert" in the "order by" clause makes MySql sort lexicographically. 
  # This is necessary, or bedtools will complain later.
  #
  
  run(
    qq(
	mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
	  select
	    rf.seq_region_id,
	    rf.seq_region_start - rf.bound_start_length as start,
	    rf.seq_region_end + rf.bound_end_length as end,
	    regulatory_activity_id,
	    epigenome_id
	  from
	    regulatory_feature rf
	    join regulatory_activity using (regulatory_feature_id)
	  order by
	    convert(rf.seq_region_id, char),
	    start;
	' \\
 	> $cell_type_regulatory_features
    )
  );

  run(
    qq(
	mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
	  select
	    af.seq_region_id,
	    af.seq_region_start,
	    af.seq_region_end,
	    af.annotated_feature_id,
	    fs.epigenome_id
	  FROM
	    annotated_feature af
	    LEFT JOIN feature_set fs USING(feature_set_id)
	  ORDER BY
	    convert(af.seq_region_id,    char),
	    af.seq_region_start;
	' \\
	> $annotations
    )
  );

  # run(
  #   qq(
	# mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e '
	#   select
	#     mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, fs.epigenome_id
	#   from
	#     motif_feature mf
	#     JOIN associated_motif_feature amf USING(motif_feature_id)
	#     JOIN annotated_feature af USING(annotated_feature_id)
	#     JOIN feature_set fs USING(feature_set_id)
	#   GROUP BY
	#     motif_feature_id,
	#     epigenome_id
	#   ORDER BY
	#     seq_region_id,
	#     seq_region_start;
	#  ' > $motifs
  #   )
  # );

  run(
    qq(
    mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
      select
        mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, "PHONY" 
      from
        motif_feature mf
      ORDER BY
        convert(seq_region_id,    char),
        seq_region_start;
  ' > $motifs
    )
  );

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -sorted -wa -wb -a $cell_type_regulatory_features -b $annotations | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"annotated\"} ' > $out");
  run("bedtools intersect -sorted -wa -wb -a $cell_type_regulatory_features -b $motifs      | awk 'BEGIN {OFS = \"\\t\"} {print \$4,\$9,\"motif\"} ' > $out");

  # Load into database
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -P $options->{port} -e 'TRUNCATE TABLE regulatory_evidence;'");
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -P $options->{port} -e 'LOAD DATA LOCAL INFILE \"$out\" INTO TABLE regulatory_evidence;'");

#   print "cell_type_regulatory_features = $cell_type_regulatory_features\n";
#   print "annotations = $annotations\n";
#   print "motifs = $motifs\n";
#   print "out = $out\n";

    # Unlink should be unnecessary, if the temporary files should go, unset KEEP_ALL at the beginning of the script.
#   unlink $cell_type_regulatory_features;
#   unlink $annotations;
#   unlink $motifs;
#   unlink $out;
}

=head2 create_regulatory_build_object

  Description: Updates data in metatable
  Arg1: Bio::EnsEMBL::Funcgen::RegulatoryBuild - The object representing the current regulatory build
  Arg2: boolean - Flag indicating whether this is to be considered a small update. This affects how the version string is incremented.
  Returntype: Bio::EnsEMBL::Funcgen::RegulatoryBuild
  Side effects: None

=cut

sub create_regulatory_build_object {
  my ($current_regulatory_build, $is_small_update, $funcgen_db_adaptor) = @_;

  use Bio::EnsEMBL::Funcgen::RegulatoryBuild;

  my $new_regulatory_build;

  my $current_build_version;
  my $current_build_initial_release_date;

  if (defined $current_regulatory_build) {
    $new_regulatory_build = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
        -name            => 'The new ' . $current_regulatory_build->name,
        -feature_type_id => $current_regulatory_build->feature_type_id,
        -analysis_id     => $current_regulatory_build->analysis_id,
        -is_current      => 0,
    );
    $current_build_version = $current_regulatory_build->version;
    $current_build_initial_release_date = $current_regulatory_build->initial_release_date;
  } else {
    
      use Bio::EnsEMBL::Analysis;
      use Bio::EnsEMBL::Funcgen::FeatureType;
      
      my $regulatory_build_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor->new($funcgen_db_adaptor);
      my $analysis_adaptor         = Bio::EnsEMBL::Analysis->new($funcgen_db_adaptor);
      my $feature_type_adaptor     = Bio::EnsEMBL::Funcgen::FeatureType->new($funcgen_db_adaptor);

      my $regulatory_build_logic_name        = 'Regulatory_Build';
      my $regulatory_build_feature_type_name = 'RegulatoryFeature';

      my $regulatory_build_analysis = $analysis_adaptor->fetch_by_logic_name($regulatory_build_logic_name);

      if (! defined $regulatory_build_analysis) {
        die();
      }

      my $regulatory_build_feature_type = $feature_type_adaptor->fetch_by_name($regulatory_build_feature_type_name);

      if (! defined $regulatory_build_feature_type) {
        die();
      }

      my $new_regulatory_build = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
              -name            => 'The Ensembl Regulatory Build',
              -feature_type    => $regulatory_build_feature_type,
              -analysis        => $regulatory_build_analysis,
              -is_current      => 0,
          );

#       print Dumper($regulatory_build);
      $regulatory_build_adaptor->store($new_regulatory_build);

      $current_build_version = '0.0';
      $current_build_initial_release_date = '';
  }

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
  # Seriously localtime, you're useless
  $year += 1900;
  $mon += 1;
  my ($main, $update);
  my $version = $current_build_version;
  if (defined $version) {
    ($main, $update) = split('.', $version);
    if ($is_small_update) {

      my $initial_release_date = $current_regulatory_build->initial_release_date;
      $new_regulatory_build->initial_release_date($initial_release_date);

      $update += 1;
    } else {

      $new_regulatory_build->initial_release_date("$year-$mon");

      $main += 1;
      $update = 0;
    }
  } else {

    $new_regulatory_build->initial_release_date("$year-$mon");
    $main = 1;
    $update = 0;
  }
  $new_regulatory_build->version(join('.', ($main, $update)));
  $new_regulatory_build->last_annotation_update("$year-$mon");
  return $new_regulatory_build;
}

1;
