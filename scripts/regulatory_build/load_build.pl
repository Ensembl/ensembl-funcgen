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

perl load_new_assembly.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38 --base_dir hg38/

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
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Analysis;
use File::Temp qw/ tempfile tempdir /;
${File::Temp::KEEP_ALL} = 1;

my $dead_rgb = '225,225,225';
my $poised_rgb = "192,0,190";
my $repressed_rgb = "127,127,127";
my $na_rgb = "255,255,255";

our %label_description= (
  'ctcf'=>'CTCF Binding Site',
  'distal'=>'Predicted enhancer',
  'proximal'=>'Predicted promoter flanking region',
  'tss'=>'Predicted promoter',
  'tfbs'=>'Transcription factor binding site',
  'dnase'=>'Open chromatin region'
);
our $start_time = time;

main();

####################################################
## Overview
####################################################
sub main {
  print_log("Getting options\n");
  my $options = get_options();
  print_log("Ensuring previous build is archived\n");
  archive_previous_build($options);
  print_log("Connecting to database\n");
  my $db = connect_db($options);
  print_log("Getting analysis\n");
  my $analysis = get_analysis($db);
  print_log("Getting cell types\n");
  my $ctypes = get_cell_type_names($options->{base_dir}, $db);
  print_log("Getting stable ids\n");
  my $stable_id = get_stable_id($options, $db);
  print_log("Getting slices\n");
  my $slice = get_slices($db);
  print_log("Getting feature sets\n");
  my $feature_set = get_regulatory_FeatureSets($analysis, $ctypes, $db);
  print_log("Getting feature types\n");
  my $feature_type = get_feature_types($db);
  print_log("Counting active features\n");
  my $count_hash = compute_counts($options->{base_dir});
  print_log("Creating regulatory_feature table\n");
  compute_regulatory_features($options, $feature_set, $feature_type, $stable_id, $count_hash, $slice, $db);
  print_log("Creating regulatory_annotation table\n");
  compute_regulatory_annotations($options);
  print_log("Updating meta table\n");
  update_meta_table($options, $db);
}

########################################################
## Print conveninence
## Params:
## - String
## Actions:
## - Prints string, with time stamp in front
########################################################

sub print_log {
  my ($str) = @_;
  my $runtime = time - $start_time;
  print "[$runtime] $str";
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
    "small_update=s",
  ) or pod2usage( -exitval => 1);
  defined $options{base_dir} || die ("You must define the base directory!\t--base_dir XXXX\n");
  defined $options{host} || die ("You must define the destination host!\t--host XXXX\n");
  defined $options{user} || die ("You must define the user login!\t--user XXXX\n");
  defined $options{dbname} || die ("You must define the database name!\t--dbname XXXX\n");
  return \%options;
}

####################################################
## Archiving old build
####################################################
sub archive_previous_build {
  my ($options) = @_;
  my $connection = "mysql -u $options->{user} -h $options->{host} -D $options->{dbname}";
  if (defined $options->{port}) {
    $connection .= " -P $options->{port}";
  }
  if (defined $options->{pass}) {
    $connection .= " -p$options->{pass}";
  }
  run("$connection -e 'UPDATE data_set SET name = CONCAT(name, \".ARCHIVED\") WHERE name LIKE \"RegulatoryFeatures:%\" AND name NOT LIKE \"%.ARCHIVED\"'");
  run("$connection -e 'UPDATE feature_set SET name = CONCAT(name, \".ARCHIVED\") WHERE name LIKE \"RegulatoryFeatures:%\" AND name NOT LIKE \"%.ARCHIVED\"'");
  run("$connection -e 'UPDATE meta SET meta_key = CONCAT(meta_key, \".ARCHIVED\") WHERE meta_key LIKE \"regbuild.%\" AND meta_key NOT LIKE \"%.ARCHIVED\"'");
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
# Params: Bio::EnsEMBL::Funcgen::DBAdaptor
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
       -description     => q({'reg_feats' => 'Features from <a href="/info/genome/funcgen/index.html" class="cp-external">Ensembl Regulatory Build</a>.', 'core' => 'Sites enriched for marks of open chromatin (e.g. DNase1) or transcription factor binding sites.', 'non_core' =>  'Sites enriched for histone modifications or polymerase binding.'}),
       -display_label   => 'Regulatory Build',
       -displayable     => 1,
       -web_data        => q({'type' => 'fg_regulatory_features', 'name' => 'Reg. Feats',  'display' =>'off', 'depth' => 10, 'default' => {'contigviewbottom' => 'normal', 'generegview' => 'normal'} }),
      );

    $aa->store($analysis);
    return $analysis;
  } else {
    return $ana;
  }
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

#####################################################
# Get list of cell types
#####################################################
# Params: The base_directory name
#####################################################

sub get_cell_type_names {
  my ($base_dir, $db) = @_;
  my $cta = $db->get_CellTypeAdaptor();  
  my %cell_type_from_clean = ();
  for my $cell_type (@{$cta->fetch_all()}) {
    $cell_type_from_clean{clean_name($cell_type->name)} = $cell_type;
    defined $cell_type || die("Unrecognized cell type name $cell_type\n");
  }

  my @cell_types = ();
  foreach my $file (glob "$base_dir/projected_segmentations/*.bb") {
    my $cell_type_name = basename $file;
    $cell_type_name =~ s/\.bb//;
    if (!exists $cell_type_from_clean{$cell_type_name}) {
      die("Celltype $cell_type_name unknown!");
    }
    push @cell_types, $cell_type_from_clean{$cell_type_name}->name;
  }
  return \@cell_types;
}

#####################################################
# Convenience wrapper to run the commandline safely 
#####################################################
# Params: command line command string
#####################################################

sub run {
  my ($cmd) = @_;
  print_log($cmd . "\n");
  system($cmd) && die("Failed when running command:\n$cmd\n");
}


#####################################################
# Alternate assign stable ids to features
#####################################################
# Params: Bio::EnsEMBL::Funcgen::DBAdaptor
#####################################################

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
  foreach my $feature (sort cmp_features @{$db->get_adaptor('FeatureSet')->fetch_by_name('RegulatoryFeatures:MultiCell.ARCHIVED')->get_all_Features()}) {
    print $ofh join("\t", ($feature->seq_region_name, $feature->bound_start, $feature->bound_end, $feature->feature_type->name, substr($feature->stable_id, 4))) . "\n";
  }
  close $ofh;

  my ($fh, $new) = tempfile();
  close $fh;

  run("bigBedToBed $options->{base_dir}/overview/RegBuild.bb $new");
  my ($overlaps, $max_id) = get_overlaps_between_files($old, $new);

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
  unlink $old;

  return \%stable_id_hash;
}

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

#####################################################
# Get slice for each chromosome
#####################################################
# Params: - DBAdaptor
#####################################################

sub get_slices {
  my ($db) = @_;
  my $slices = $db->get_adaptor("Slice")->fetch_all('toplevel',undef,0,1);
  my %hash = ();
  foreach my $slice (@{$slices}) {
    $hash{$slice->seq_region_name} = $slice;
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

sub get_regulatory_FeatureSets {
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
    print_log("\tCreating feature set for cell type $ctype\n");
    my ($desc, $dlabel, $fset_name);
    
    $fset_name = "RegulatoryFeatures:$ctype";
    $dlabel    = "Reg.Feats $ctype";

    if($ctype eq 'MultiCell') {
      $desc = 'Consensus RegulatoryFeature regions';
    }
    else {
      $desc = "$ctype specific RegulatoryFeatures";
    }

    print_log("Defining FeatureSet:\t$fset_name");

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
    $dsa->store_status('DISPLAYABLE', $dset);
  }
    
  foreach my $fset(@fsets){
    foreach my $fs_state(('DISPLAYABLE','IMPORTED','IMPORTED_GRCh38','MART_DISPLAYABLE')) {
      $fsa->store_status($fs_state, $fset);
    }
  }

  print_log("Got RegulatoryFeature sets for CellTypes:\t".join(', ', keys %rf_sets));
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
    print_log("\tSearching for supporting sets on cell type $ctype\n");
    my $CellType = $cta->fetch_by_name($ctype);
    my @ssets = ();
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
  print_log("\tCounting in file $filename\n");
  my ($fh, $tmp_name) = tempfile();
  run("bigBedToBed $filename $tmp_name");
  while (my $line  = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;
    if (!defined $count_hash->{$name}) {
      $count_hash->{$name} = 0;
    }
    if ($rgb ne $dead_rgb && $rgb ne $poised_rgb && $rgb ne $repressed_rgb && $rgb ne $na_rgb) {
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
  my ($options, $feature_set, $feature_type, $stable_id, $count_hash, $slice, $db) = @_;
  my $rfa = $db->get_adaptor("RegulatoryFeature");
  foreach my $cell_type (keys %{$feature_set}) {
    load_celltype_build($options->{base_dir}, $feature_set->{$cell_type}, $stable_id, $count_hash, $slice, $cell_type, $feature_type, $rfa);
  }
}

sub load_celltype_build {
  my ($base_dir, $feature_set, $stable_id, $count_hash, $slice, $cell_type, $feature_type, $rfa) = @_;
  my ($tmp, $tmp_name) = tempfile();
  print_log("\tProcessing data from cell type $cell_type\n");
  my $bigbed;
  if ($cell_type eq 'MultiCell') {
    $bigbed = "$base_dir/overview/RegBuild.bb";
  } else {
    my $cell_type_name = clean_name($cell_type);
    $bigbed = "$base_dir/projected_segmentations/$cell_type_name.bb";
  }
  run("bigBedToBed $bigbed $tmp_name");
  process_file($tmp, $feature_set, $stable_id, $count_hash, $slice, $feature_type, $rfa);
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
    if ($rgb ne $dead_rgb && $rgb ne $poised_rgb && $rgb ne $repressed_rgb && $rgb ne $na_rgb) {
      $has_evidence = 1; # TODO 4 way
    }

    exists $feature_type->{$feature_type_str} || die("Could not find feature type for $feature_type_str\n".join("\t", keys %{$feature_type})."\n");
    exists $slice->{$chrom} || die("Could not find slice type for $chrom\n".join("\t", keys %{$slice})."\n");
    exists $stable_id->{$name} || die("Could not find stable ID for feature # $name\n");
    exists $count_hash->{$name} || die("Could not find count for feature # $name\n");

    push @features, Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast({
      slice         => $slice->{$chrom},
      start         => $thickStart + 1,
      end           => $thickEnd,
      strand        => 0,
      display_label => '\\N',
      set           => $feature_set,
      feature_type  => $feature_type->{$feature_type_str},
      _bound_lengths => [$thickStart - $start, $end - $thickEnd],
      has_evidence  => $has_evidence,
      cell_type_count => $count_hash->{$name},
      stable_id     => $stable_id->{$name},
      analysis      => $feature_set->analysis,
      #     adaptor       => $rfa,
      projected     => 0,
    });

    if (scalar @features > 10000) {
      $rfa->store(@features);
      @features = ();
    }
  }

  if (scalar @features > 0) {
    $rfa->store(@features);
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
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf LEFT JOIN feature_set fs USING(feature_set_id) LEFT JOIN cell_type USING(cell_type_id) WHERE cell_type.name LIKE \"MultiCell\" ORDER by rf.seq_region_id, rf.seq_region_start - rf.bound_start_length'| sort -k1,1 -k2,2n > $multicell_regulatory_features");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf LEFT JOIN feature_set fs USING(feature_set_id) LEFT JOIN cell_type USING(cell_type_id) WHERE cell_type.name NOT LIKE \"MultiCell\" ORDER by rf.seq_region_id, rf.seq_region_start - rf.bound_start_length'| sort -k1,1 -k2,2n > $cell_type_regulatory_features");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select af.seq_region_id, af.seq_region_start, af.seq_region_end, af.annotated_feature_id, fs.cell_type_id FROM annotated_feature af LEFT JOIN feature_set fs USING(feature_set_id) ORDER BY af.seq_region_id, af.seq_region_start;' | awk '$2 <= $3' > $annotations");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, fs.cell_type_id from motif_feature mf LEFT JOIN associated_motif_feature amf USING(motif_feature_id) LEFT JOIN annotated_feature af USING(annotated_feature_id) LEFT JOIN feature_set fs USING(feature_set_id) GROUP BY motif_feature_id, cell_type_id ORDER BY seq_region_id, seq_region_start;' > $motifs");

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -sorted -wa -wb -a $cell_type_regulatory_features -b $annotations | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"annotated\"} ' > $out");
  run("bedtools intersect -sorted -wa -wb -a $multicell_regulatory_features -b $motifs | awk 'BEGIN {OFS = \"\\t\"} {print \$4,\$9,\"motif\"} ' >> $out");
  run("bedtools intersect -sorted -wa -wb -a $cell_type_regulatory_features -b $motifs | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"motif\"} ' >> $out");

  # Load into database
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'TRUNCATE TABLE regulatory_attribute;'");
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'LOAD DATA LOCAL INFILE \"$out\" INTO TABLE regulatory_attribute;'");

  unlink $regulatory_features;
  unlink $annotations;
  unlink $motifs;
  unlink $out;
}

#####################################################
# Updates data in metatable 
#####################################################
# Params: - Bio::EnsEMBL::Funcgen::DBAdaptor
#####################################################

sub update_meta_table {
  my ($options, $db) = @_;
  my $mc = $db->get_MetaContainer();
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
  $year += 1900;
  $mon += 1;
  my $version = $mc->single_value_by_key('regbuild.version.ARCHIVED');
  if (defined $version) {
    my ($main, $update) = split(".", $version);
    if (defined $options->{small_update}) {
      $update += 1;
    } else {
      $main += 1;
      $update = 0;
      $mc->update_key_value('regbuild.initial_release_date', "$year-$mon");
    }
    $mc->update_key_value('regbuild.version', join(".", ($main, $update)));
    $mc->update_key_value('regbuild.last_annotation_update', "$year-$mon");
  } else {
    $mc->set_key_value('regbuild.initial_release_date', "$year-$mon");
    $mc->set_key_value('regbuild.version', "1.0");
    $mc->set_key_value('regbuild.last_annotation_update', "$year-$mon");
  }
}

1;
