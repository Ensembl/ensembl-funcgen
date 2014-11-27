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

load_build.pl

=head1 SYNOPSIS

perl load_new_assembly.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38 --dnadb_name homo_sapiens_core_76_38 --dnadb_user $user2 --dnadb_host $host2 --base_dir hg38/

=head1 DESCRIPTION

Loads a segmentation BigBed file annotated by the new regulatory build into the database.
Params:
	* base_dir: directory with assembly name (e.g. ./hg38) create by build
	* segmentation_name: name of segmentation as named in build
	* cell_type_name: name of cell type as named in build.

In short, the file $base_dir/segmentations/$segmentation/$cell_type_name.bb must exist

Also bigBedToBed must be on the commandline.

This script is pretty database intensive (think inserting 4Mo entries) so keep an eye 
when running in parallel. In the good times, each job takes ~1h to run. 

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Analysis;
use File::Temp qw/ tempfile tempdir /;

my $dead_rgb = '225,225,225';
my $poised = "192,0,190";
my $repressed = "127,127,127";

main();

####################################################
## Overview
####################################################
sub main {
  print "Getting options\n";
  my $options = get_options();
  print "Connecting to database\n";
  my $db = connect_db($options);
  print "Getting analysis\n";
  my $analysis = get_analysis($db);
  print "Getting cell types\n";
  my $ctypes = get_cell_type_names($options->{base_dir});
  print "Getting stable ids\n";
  my $stable_id = get_stable_id($options, $db);
  print "Getting slices\n";
  my $slice = get_slices($db);
  print "Getting feature sets\n";
  my $feature_set = get_regulatory_FeatureSets($analysis, $ctypes, $db);
  print "Getting feature types\n";
  my $feature_type = get_feature_types($db);
  print "Counting active features\n";
  my $count_hash = compute_counts($options->{base_dir});
  print "Creating regulatory_feature table\n";
  compute_regulatory_features($options, $feature_set, $feature_type, $stable_id, $count_hash, $slice);
  print "Creating regulatory_annotation table\n";
  compute_regulatory_annotations($options);
}

####################################################
## Command line options
####################################################
sub get_options {
  my %options = ();
  GetOptions (
    \%options,
    "base_dir|b=s",
    "pass|p=s",
    "port=s",
    "host|h=s",
    "user|u=s",
    "dbname|d=s",
  ) or pod2usage( -exitval => 1);
  defined $options{base_dir} || die ("You must define the base directory!\t--base_dir XXXX\n");
  defined $options{host} || die ("You must define the destination host!\t--base_dir XXXX\n");
  defined $options{user} || die ("You must define the user login!\t--base_dir XXXX\n");
  defined $options{dbname} || die ("You must define the database name!\t--base_dir XXXX\n");
  return \%options;
}

####################################################
## Connecting to the DB
####################################################
sub connect_db {
  my ($options) = @_;
  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $options->{host},
     -user   => $options->{user},
     -dbname => $options->{dbname},
     -pass   => $options->{pass},
     -port   => $options->{port},
     #-dnadb  => $options->{cdb},
     -group  => 'funcgen',#Should be set as default in adaptor new method
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

#####################################################
# Create/Get analysis for Build
#####################################################
# Params: A DBAdaptor
#####################################################

sub get_analysis {
### Check whether analysis is already stored
#TO DO Update the description text here? Use flat file import?
#my $program_name = ($0 =~ s'.*/''g);
  my ($db) = @_;

  my $aa  = $db->get_AnalysisAdaptor();
  my $ana = $aa->fetch_by_logic_name('Regulatory_Build');

  if ( not defined $ana ) {
    my $analysis = Bio::EnsEMBL::Analysis->new
      (
       -logic_name      => 'Regulatory_Build',
       -db              => 'NULL',
       -db_version      => 'NULL',
       -db_file         => 'NULL',
       -program         => 'NULL',
       -program_version => 'NULL',
       -program_file    => 'NULL', #Could use $program_name here, but this is only part of the build
       -gff_source      => 'NULL',
       -gff_feature     => 'NULL',
       -module          => 'NULL',
       -module_version  => 'NULL',
       -parameters      => 'NULL',
       -created         => 'NULL',
       -description     => q({'reg_feats' => 'Features from <a href="http://www.ensembl.org/info/docs/funcgen/index.html" class="cp-external">Ensembl Regulatory Build</a>.', 'core' => 'Sites enriched for marks of open chromatin (e.g. Dnase1) or transcription factor binding sites.  Used to define the Regulatory Feature core regions in the <a href="http://www.ensembl.org/info/docs/funcgen/index.html" class="cp-external">Ensembl Regulatory Build</a>.', 'non_core' => 'Sites enriched for histone modifications or polymerase binding.  Used to define Regulatory Feature bound regions in the <a href="http://www.ensembl.org/info/docs/funcgen/index.html" class="cp-external">Ensembl Regulatory Build</a>.'}),
       -display_label   => 'Regulatory Build',
       -displayable     => 1,
       -web_data        => q({'type' => 'fg_regulatory_features', 'name' => 'Reg. Feats',  'display' =>'off', 'depth' => 10, 'default' => {'contigviewbottom' => 'normal', 'generegview' => 'normal'} }),
      );

    $aa->store($analysis);
  }
  return $aa;
}

#####################################################
# Get list of cell types
#####################################################
# Params: The base_directory
#####################################################

sub get_cell_type_names{
  my ($base_dir) = @_;
  my @cell_types = ();
  foreach my $file (glob "$base_dir/projected_segmentations/*.bb") {
    my $cell_type = basename $file;
    $cell_type =~ s/\.bb//;
    print "\tCELL TYPE $cell_type\n";
    push @cell_types, $cell_type;
  }
  return \@cell_types;
}

#####################################################
# Convenience wrapper to run the commandline safely 
#####################################################
# Params: command line command
#####################################################

sub run {
  my ($cmd) = @_;
  print cmd;
  system($cmd) && die("Failed when running command:\n$cmd\n");
}


#####################################################
# Alternate assign stable ids to features
#####################################################
# Params: Bed file with previous Build
#####################################################

sub get_stable_id {
  my ($options, $db) = @_;
  my ($fh, $new) = tempfile();
  close $fh;

  run("bigBedToBed $options->{base_dir}/overview/RegBuild.bb $new");
  my ($overlaps, $max_id) = get_overlaps_between_files($old, $new);

# Go through overlaos in order of increasing overlap length. This means that you should always 
# overwrite an overlap with a later one. 
  sub cmp_overlaps {
    return $a->[14] <=> $b->[14];
  }
  my %old_pref = ();
  my %new_pref = ();
  foreach my $entry (sort cmp_overlaps @{$overlaps}) {
    $old_pref{$entry->[3]} = $entry->[8];
    $new_pref{$entry->[8]} = $entry->[3];
  }

# Output new_id/old_id pairs if mutual best hits, else create new ids starting from 
# the maximum of pre-existing IDs + 1
  my %stable_id_hash = ();
  my $next_free_id = $max_id + 1;
  open my $in, "<", $new;
  while (my $line = <$in>) {
    chomp $line;
    my @items = split /\t/, $line;
    my $new_id = @items[3];

    my $stable_id = undef;

    if (exists $new_pref{$new_id} && $old_pref{$new_pref{$new_id}} eq $new_id) {
      $stable_id = $new_pref{$new_id};
    } else {
      $stable_id = $next_free_id;
      $next_free_id += 1;
    }

    $stable_id_hash{$new_id} = $stable_id;
  }
  close $in;
  unlink $new;

  return \%stable_id_hash;
}

sub get_overlaps_between_files {
  my ($old, $new) = @_;
  
# Compute overlaps between regions defined in both file
  my ($fh, $filename) = tempfile();
  my $cmd = "bedtools intersect -a $old -b $new -wo > $filename";
  system($cmd) && die("Failed bedtools command\n");

# Parse output
  my @overlaps = ();
  my $max_id = undef;
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    # Check no feature type incompatibility
    if ($items[4] ne 'NONE') {
      # TODO map strings to feature type
      if (0) {
        next;  
      }
    }
    push @overlaps, \@items;

    # Look for maximum ID number among the old features
    my $id = substr($items[3], 4);
    if (!defined $max_id || $max_id < $id) {
      $max_id = $id;
    }
  }
  unlink $filename;
  
  return (\@overlaps, $max_id);
}

#####################################################
# Get slice for each chromosome
#####################################################
# Params: - DBAdaptor
#####################################################

sub get_slices {
  my ($db) = @_;
  my $slices = $db->get_adaptor("Slice")->fetch_all();
  my %hash = ();
  foreach my $slice (@{$slices}) {
    $hash{$slice->name} = $slice;
  }
  return \%hash;
}

#####################################################
# Create/Get FeatureSet for each cell type
#####################################################
# Params: - Analysis object
#         - Array ref of cell types
#         - DBAdaptor
#####################################################

sub get_regulatory_FeatureSets{
  my ($analysis, $ctypes, $db) = @_;
  my %rf_sets;
  my $fta = $db->get_FeatureTypeAdaptor();
  my $ftype = $fta->fetch_by_name('RegulatoryFeature');
  
  if (! $ftype) {
    $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
      (
       -name        => 'RegulatoryFeature',
       -description => 'Ensembl Regulatory Feature',
       -class       => 'Regulatory Feature',
      );
    
    ($ftype) = @{$fta->store($ftype)};
  }
  
  my (@dsets, @fsets);
  my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper();
  my $dsa = $db->get_DataSetAdaptor();
  my $fsa = $db->get_FeatureSetAdaptor();
  my $cta = $db->get_CellTypeAdaptor();
  my $ctype_ssets = get_cell_type_supporting_sets($ctypes, $cta, $fsa);

# make sure that attribute sets also contain focus sets 
# and ctype_fsets contain core fsets

  foreach my $ctype (@{$ctypes},('MultiCell')) {  
    print "\tCreating feature set for cell type $ctype\n";
    my ($desc, $dlabel, $fset_name);
    
    $fset_name = "RegulatoryFeatures:$ctype";
    $dlabel    = "Reg.Feats $ctype";

    if($ctype eq 'MultiCell') {
      $desc = 'Generic RegulatoryFeature focus regions';
    }
    else {
      $desc = "$ctype specific RegulatoryFeatures";
    }

    $helper->log("Defining FeatureSet:\t$fset_name");

    my $description;
    if ($ctype eq 'MultiCell') {
      $description = "MultiCell consensus RegulatoryFeatures";
    } else {
      $description = "$ctype specific RegulatoryFeatures";
    }

    my $dset = $helper->define_DataSet
      (
      -NAME                 => "RegulatoryFeatures:$ctype", 
      -FEATURE_CLASS        => 'regulatory', 
      -SUPPORTING_SETS      => $ctype_ssets->{$ctype},
      -DBADAPTOR            => $db,
      -FEATURE_SET_ANALYSIS => $analysis,#i.e. RegulatoryBuild
      -CELL_TYPE            => $cta->fetch_by_name($ctype),
      -FEATURE_TYPE         => $ftype,
      -ROLLBACK             => 'feature_set',
      -display_label        => "Reg.Feats $ctype",
      -description          => $description,
      ); 

    #Always overwrite in case we have redefined the sets
    &store_regbuild_meta_strings($dset, 1);
    $rf_sets{$ctype} = $dset->product_FeatureSet;

    #Set states
    #Move to Utils/RegBuilder.pm?
     
    push @dsets, $dset;

    foreach my $sset(@{$ctype_ssets->{$ctype}}){
      my $ss_dset = $dsa->fetch_by_product_FeatureSet($sset);

      if(! $ss_dset){
      die("Could not find DataSet for FeatureSet:\t".$sset->name);
      }
      
      push @dsets, $ss_dset;
    }

    push @dsets, ($dset->product_FeatureSet, @{$ctype_ssets->{$ctype}});  
  }

  #Set states
  foreach my $dset(@dsets){
    # TODO Remove hardcoded values?
    foreach my $ds_state(('DISPLAYABLE')){
      $dsa->store_status($ds_state, $dset);
    }
  }
    
  foreach my $fset(@fsets){
    # TODO Remove hardcoded values?
    foreach my $fs_state(('DISPLAYABLE','IMPORTED','IMPORTED_GRCh38','MART_DISPLAYABLE')) {
      $fsa->store_status($fs_state, $fset);
    }
  }

  $helper->log("Got RegulatoryFeature sets for CellTypes:\t".join(', ', keys %rf_sets));
  return \%rf_sets;
}

#This could move to DataSetAdaptor and be called from define_and_validate_sets
#If focus set info was stored in data set

#Thsi needs to be made available to the Helper for HC/Updating
#in case of archive failure

sub store_regbuild_meta_strings{
  my ($dset, $overwrite) = @_;

  my $ds_adaptor = $dset->adaptor;
  $ds_adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::DataSet', $dset);
  my ($sql, $meta_value, $reg_string, $cmd, $msg);
  my $fset = $dset->product_FeatureSet;

  if (! defined $fset ||
      $fset->feature_class ne 'regulatory') {
    die('You must provide a DataSet with an associated \'regulatory\' product FeatureSet');
  }

  my @ssets = @{$dset->get_supporting_sets};

  if (! @ssets) {
    ('You must provide a DataSet with associated supporting sets');
  }

  my $ctype = (defined $fset->cell_type) ?  $fset->cell_type->name : 'core';

  ### build and import regbuild strings by feature_set_id and feature_type_id

  #($meta_value) = $ds_adaptor->db->dbc->db_handle->selectrow_array("select meta_value from meta where meta_key='regbuild.${ctype}.feature_set_ids'");
  #$reg_string = join(',', map {$_->dbID} sort {$a->name cmp $b->name} @ssets);

  my %reg_strings = 
    (
     "regbuild.${ctype}.feature_set_ids" => join(',', map {
       $_->dbID} sort {$a->name cmp $b->name
                     } @ssets),
     "regbuild.${ctype}.feature_type_ids" => join(',', map {
       $_->feature_type->dbID} sort {$a->name cmp $b->name
                                   } @ssets),
    );
  
  my @ffset_ids;

  #this should be sorted to avoid string mismatches with the same contents.
  $reg_strings{"regbuild.${ctype}.focus_feature_set_ids"} = join(', ', @ffset_ids);

  foreach my $meta_key (keys %reg_strings) {
    ($meta_value) = $ds_adaptor->db->dbc->db_handle->selectrow_array("select string from regbuild_string where name='${meta_key}'");

    if (! defined $meta_value){
      eval { $ds_adaptor->db->dbc->do("insert into regbuild_string (name, string) values ('${meta_key}', '$reg_strings{${meta_key}}')") };
      die("Couldn't store $meta_key in regbuild_string table.\n$@") if $@;
    } 
    elsif ($meta_value ne $reg_strings{$meta_key}) {

      if ($overwrite) {
        warn "Overwriting old $meta_key regbuild_string entry:\t$meta_value\nwith:\t\t\t\t\t\t\t\t\t\t".$reg_strings{$meta_key}."\n";
        eval { $ds_adaptor->db->dbc->do("update regbuild_string set string='$reg_strings{${meta_key}}' where name ='${meta_key}'") };
        die("Couldn't store $meta_key in regbuild_string table.\n$@") if $@;
      }
      else{
        die "$meta_key already exists in regbuild_string table and does not match\nOld\t$meta_value\nNew\t$reg_strings{${meta_key}}\nPlease archive previous RegulatoryBuild.\n";
      }
    }
  }
  
  return \%reg_strings;
}

sub get_cell_type_supporting_sets {
  my ($ctypes, $cta, $fsa) = @_;
  my %ctype_ssets;
  $ctype_ssets{MultiCell} = [];
  foreach my $ctype (@{$ctypes}) {  
    print "\tSearching for supporting sets on cell type $ctype\n";
    my $CellType = $cta->fetch_by_name($ctype);
    my @ssets = ();
    print $ctype."\n";
    print $CellType."\n";
    foreach my $fs (@{$fsa->fetch_all_by_CellType($CellType)}) {
      if ($fs->feature_class eq 'annotated') {
        push @ssets, $fs;
      }
    }
    $ctype_ssets{$ctype} = \@ssets;
    push @{$ctype_ssets{MultiCell}}, @ssets;
  }
  return \%ctype_ssets;
}

#####################################################
# Count the number of active features for each temporary
# ID across all cell types
#####################################################
# Params: - Base directory
#####################################################

sub compute_counts {
  my ($base_dir) = @_;
  my $count_hash = {};
  foreach my $file (glob "$base_dir/projected_segmentations/*.bb") {
    count_active($file , $count_hash);
  }
  return $count_hash;
}

sub count_active {
  my ($filename, $count_hash) = @_;
  print "\tCounting in file $filename\n";
  my ($fh, $tmp_name) = tempfile();
  run("bigBedToBed $filename $tmp_name");
  while (my $line  = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;
    if (!defined $count_hash->{$name}) {
      $count_hash->{$name} = 0;
    }
    if ($rgb ne $dead_rgb && $rgb ne $poised && $rgb ne $repressed) {
      $count_hash->{$name} += 1;
    }
  }
  close $fh;
  unlink $tmp_name;
}

#####################################################
# Get a hashref from feature type name to FeatureType 
# object
#####################################################
# Params: - DBAdaptor object
#####################################################
 
sub get_feature_types {
  my ($db) = @_;
  my $fta = $db->get_adaptor("FeatureType");
  my %long_name= (
    'ctcf'=>'CTCF Binding Site',
    'distal'=>'Enhancer',
    'proximal'=>'Promoter Flanking Region',
    'tss'=>'Promoter',
    'tfbs'=>'TF binding site',
    'open'=>'Open chromatin'
  );
  my %description= (
    'ctcf'=>'CTCF Binding Site',
    'distal'=>'Predicted enhancer',
    'proximal'=>'Predicted promoter flanking region',
    'tss'=>'Predicted promoter',
    'tfbs'=>'Transcription factor binding site',
    'open'=>'Open chromatin region'
  );
  my %so_accession = (
    'ctcf'=>'SO:0001974',
    'distal'=>'SO:0000165',
    'proximal'=>'SO:0001952',
    'tss'=>'SO:0000167',
    'tfbs'=>'SO:0000235',
    'open'=>'SO:0001747'
  );
  my %so_name = (
    'ctcf'=>'CTCF_binding_site',
    'distal'=>'enhancer',
    'proximal'=>'promoter_flanking_region',
    'tss'=>'promoter',
    'tfbs'=>'TF_binding_site',
    'open'=>'open_chromatin_region'
  );

  my $feature_type = {};
  for my $label (('ctcf','distal','proximal','tss','tfbs','open')) {
    $feature_type->{$label} = $fta->fetch_by_name($long_name{$label});

    if (! defined $feature_type->{$label}) {
      my $ft = Bio::EnsEMBL::Funcgen::FeatureType->new(
	-name => $long_name{$label},
	-class => 'Regulatory Feature',
	-description => $description{$label},
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

#####################################################
# Creates the actual RegulatoryFeature objects
#####################################################
# Params: - Base directory
#         - feature_set: Hash ref: cell type name -> FeatureSet
#         - stable_id: Hash ref: temporary id -> new id
#         - count_hash: Hash ref: temporary id -> count
#         - seq_region_ids: Hash ref: chromosome name -> seq_region_id
#         - host
#         - user
#         - pass
#         - dbname
#####################################################

sub compute_regulatory_features {
  my ($options, $feature_set, $feature_type, $stable_id, $count_hash, $slice) = @_;
  foreach my $cell_type (keys %{$feature_set}) {
    load_celltype_build($options->{base_dir}, $feature_set->{$cell_type}, $stable_id, $count_hash, $slice, $cell_type, $feature_type);
  }
  # TODO Remove previous???
}

sub load_celltype_build {
  my ($base_dir, $feature_set, $stable_id, $count_hash, $slice, $cell_type, $feature_type) = @_;
  my ($tmp, $tmp_name) = tempfile();
  print "\tProcessing data from cell type $cell_type\n";
  my $bigbed;
  if ($cell_type eq 'MultiCell') {
    $bigbed = "$base_dir/overview/RegBuild.bb";
  } else {
    $bigbed = "$base_dir/projected_segmentations/$cell_type.bb";
  }
  run("bigBedToBed $bigbed $tmp_name");
  process_file($tmp, $feature_set, $stable_id, $count_hash, $slice, $feature_type);
  close $tmp;
  unlink $tmp_name;
}

sub process_file {
  my ($fh, $feature_set, $stable_id, $count_hash, $slice, $feature_type, $rfa) = @_;
  my @features = ();

  while (my $line = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;;
    my ($feature_type_str, $number) = split /_/, $name;
    my $has_evidence = 0;
    if ($rgb ne $dead_rgb && $rgb ne $poised && $rgb ne $repressed) {
      $has_evidence = 1; # TODO 4 way
    }

    exists $feature_type->{$feature_type_str} || die("Could not find feature type for $feature_type\n".join("\t", keys %{$feature_type})."\n");
    exists $slice->{$chrom} || die("Could not find slice type for $chrom\n".join("\t", keys %{$slice})."\n");
    exists $stable_id->{$number} || die("Could not find stable ID for feature # $number\n");
    exists $count_hash->{$number} || die("Could not find count for feature # $number\n");

    push @features, Bio::EnsEMBL::Funcgen::RegulatoryFeature->new(
      -SLICE => $slice->{$chrom};
      -START         => $thickStart + 1,
      -END           => $thickEnd,
      -STRAND        => 0,
      -DISPLAY_LABEL => '\\N',
      -FEATURE_SET   => $feature_set,
      -FEATURE_TYPE  => $feature_type->{$feature_type_str},
      -ATTRIBUTE_CACHE => \%attr_cache, # TODO 
                     => $thickStart - $start, # TODO
                     => $end - $thickEnd, # TODO
                     => $has_evidence, # TODO
                     => $count->{$number}, # TODO
                     => $stable_id->{$number} # TODO
    );

    if (scalar @features > 10000) {
      $rfa->store(@features);
      @features = ();
    }
  }

  if (scalar @features > 0) {
    $rfa->store();
  }
}

#####################################################
# Assign motifs and annotations to regulatory features
#####################################################
# Params: - host
#         - user
#         - pass
#         - dbname
#####################################################

sub compute_regulatory_annotations {
  my ($options) = @_;

  my ($tmp_fh, $regulatory_features) = tempfile();
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
  run("mysql -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf LEFT JOIN feature_set fs ON rf.feature_set_id=fs.feature_set_id ORDER by rf.seq_region_id, rf.seq_region_starti - rf.bound_start_length' > $regulatory_features");
  run("mysql -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select af.seq_region_id, af.seq_region_start, af.seq_region_end, af.annotated_feature_id, fs.cell_type_id FROM annotated_feature af LEFT JOIN feature_set fs ON af.feature_set_id=fs.feature_set_id ORDER BY af.seq_region_id, af.seq_region_start;' > $annotations");
  run("mysql -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id FROM motif_feature mf ORDER BY mf.seq_region_id, mf.seq_region_start;' > $motifs");

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $annotations | awk '\$5==\$10 {print \$4\"\t\"\$9\"\tannotated\"} ' > $out");
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $motifs | awk ' {print \$4\"\t\"\$9\"\tmotif\"} ' >> $out");

  # Load into database
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'LOAD DATA LOCAL INFILE \"$out\"} INTO TABLE regulatory_attribute;'");

  unlink $regulatory_features;
  unlink $annotations;
  unlink $motifs;
  unlink $out;
}

1;
