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

perl load_new_assembly.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38 --dnadb_name homo_sapiens_core_76_38 --dnadb_user $user2 --dnadb_host $host2 --base_dir hg38/ --stable stable_ids.txt

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

main();

sub main {
  my ($base_dir, $stable_id_filename, $pass,$port,$host,$user,$dbname,
	  $dnadb_pass,$dnadb_port,$dnadb_host,$dnadb_user,$dnadb_name,
      $assembly_schema);

  GetOptions 
    (
     "base_dir|b=s"       => \$base_dir,
     "stable|s=s"       => \$stable_id_filename,
     "assembly_schema|a=s"       => \$assembly_schema,
     "pass|p=s"       => \$pass,
     "port=s"         => \$port,
     "host|h=s"       => \$host,
     "user|u=s"       => \$user,
     "dbname|d=s"     => \$dbname,
     "dnadb_pass|p=s"       => \$dnadb_pass,
     "dnadb_port=s"   => \$dnadb_port,
     "dnadb_host|h=s" => \$dnadb_host,
     "dnadb_user|u=s" => \$dnadb_user,
     "dnadb_name|d=s" => \$dnadb_name,
     'help|?'         => sub { pos2usage(-exitval => 0, ); }
    ) or pod2usage( -exitval => 1);

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $host,
     -user   => $user,
     -dbname => $dbname,
     -pass   => $pass,
     -port   => $port,
     #-dnadb  => $cdb,
     -group  => 'funcgen',#Should be set as default in adaptor new method
     -dnadb_host => $dnadb_host,
     -dnadb_port => $dnadb_port,
     -dnadb_user => $dnadb_user,
     -dnadb_pass => $dnadb_pass,
     -dnadb_name => $dnadb_name,
    );


#Test connections
  $db->dbc->db_handle;
  $db->dnadb->dbc->db_handle;
  if(! defined $db->species){
    die("Could not get a valid species from $dbname, please check the meta species.production_name");
  }

  defined $base_dir || die ("You must define the base directory!\t--base_dir XXXX\n");
  defined $stable_id_filename || die ("You must define the stable ID mapping file\t--stable XXX\n");

  print "Getting analysis\n";
  my $analysis = get_analysis($db);
  print "Getting cell types\n";
  my $ctypes = get_cell_type_names($base_dir);
  print "Getting stable ids\n";
  my $stable_id = get_stable_id($stable_id_filename);
  print "Getting seq region ids\n";
  my $seq_region_id = get_seq_region_ids($db, $assembly_schema);
  print "Getting feature set ids\n";
  my $feature_set = get_regulatory_FeatureSets($analysis, $ctypes, $db);
  print "Counting active features\n";
  my $count_hash = compute_counts($base_dir);
  print "Creating regulatory_feature table\n";
  compute_regulatory_features($base_dir, $feature_set, $stable_id, $count_hash, $seq_region_id, $host, $user, $pass, $dbname);
  print "Creating regulatory_annotation table\n";
  compute_regulatory_annotations($host, $user, $pass, $dbname);
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

  if ( defined $ana ) {
    print "\tPrexisting\n";
    return $ana;
  }

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
# Assign stable ids to features
#####################################################
# Params: PreProcessed file
#####################################################

sub get_stable_id {
  my ($filename) = @_;
  my %stable_id = ();
  open my $fh, "<", $filename;
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split /\t/, $line;
    $stable_id{$items[0]} = $items[1];
  }
  close $fh;
  return \%stable_id;
}

#####################################################
# Alternate assign stable ids to features
#####################################################
# Params: Bed file with previous Build
#####################################################

sub get_stable_id_2 {
  my ($old, $base_dir) = @_;
  my ($fh, $new) = tempfile();
  close $fh;

  run("bigBedToBed $base_dir/overview/RegBuild.bb $new");
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
# Get seq_region_ids for each chromosome
#####################################################
# Params: - DBAdaptor
#         - Assembly schema string
#####################################################

sub get_seq_region_ids {
  my ($db, $schema) = @_;
  my %hash = ();
  my $sql = "SELECT name, seq_region_id FROM seq_region WHERE schema_build='$schema'";
  my $rows = $db->dbc->db_handle->selectall_hashref($sql, 'name');

  foreach my $name (keys %{$rows}) {
    print "\tCHROM\t$name\t$rows->{$name}->{seq_region_id}\n";
    $hash{$name} = $rows->{$name}->{seq_region_id};
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
    if ($rgb ne $dead_rgb) {
      $count_hash->{$name} += 1;
    }
  }
  close $fh;
  unlink $tmp_name;
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
  my ($base_dir, $feature_set, $stable_id, $count_hash, $seq_region_ids, $host, $user, $pass, $dbname) = @_;
  my ($outfile, $out_filename) = tempfile();
  my $feature_type_id = get_feature_type_ids();
  foreach my $cell_type (keys %{$feature_set}) {
    load_celltype_build($base_dir, $feature_set->{$cell_type}->dbID, $stable_id, $count_hash, $seq_region_ids, $outfile, $cell_type, $feature_type_id);
  }
  close $outfile;
  # TODO Remove previous???
  # TODO Load with API
  run("mysql -h $host -u $user -p$pass -D $dbname -e 'LOAD DATA LOCAL INFILE \"$outfile\" INTO TABLE regulatory_feature;'");
  unlink $out_filename;
}

sub get_feature_type_ids {
  # TODO get IDs from DB
  my %feature_type_id = (
    'ctcf'=>'179034',
    'distal'=>'179037',
    'proximal'=>'179038',
    'tss'=>'179040',
    'tfbs'=>'179194',
    'open'=>'179195',
  );
  return \%feature_type_id;
}

sub load_celltype_build {
  my ($base_dir, $feature_set_id, $stable_id, $count_hash, $seq_region_ids, $outfile, $cell_type, $feature_type_id) = @_;
  my ($tmp, $tmp_name) = tempfile();
  print "\tProcessing data from cell type $cell_type\n";
  my $bigbed;
  if ($cell_type eq 'MultiCell') {
    $bigbed = "$base_dir/overview/RegBuild.bb";
  } else {
    $bigbed = "$base_dir/projected_segmentations/$cell_type.bb";
  }
  run("bigBedToBed $bigbed $tmp_name");
  process_file($tmp_name, $feature_set_id, $stable_id, $count_hash, $outfile, $seq_region_ids, $feature_type_id);
  close $tmp;
  unlink $tmp_name;
}

sub process_file {
  my ($filename, $feature_set_id, $stable_id, $count_hash, $outfile, $seq_region_ids, $feature_type_id) = @_;
  open my $fh, "<", $filename;
  while (my $line = <$fh>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;;
    # TODO Smarter filter for valid chromosomes
    if ($chrom eq "MT") {
      next;
    }
    my ($feature_type, $number) = split /_/, $name;
    my $uniqueID = $stable_id->{$name};
    my $bounded_start_length = $thickStart - $start;;
    my $bounded_end_length = $end - $thickEnd;
    my $has_evidence = 0;
    if ($rgb ne $dead_rgb) {    
      $has_evidence = 1;
    }
    my $count = $count_hash->{$name};
    exists $feature_type_id->{$feature_type} || die("Could not find feature type for $feature_type\n".join("\t", keys %{$feature_type_id})."\n");
    print $outfile join("\t", (0, $seq_region_ids->{$chrom}, $thickStart+1, $thickEnd, 0, '\\N', $feature_type_id->{$feature_type}, $feature_set_id, $uniqueID, '\\N', 0, $bounded_start_length, $bounded_end_length, $has_evidence, $count));
    print $outfile "\n";
  }
  close $fh;
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
  my ($host, $user, $pass, $dbname) = @_;

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
  run("mysql -NB -h $host -u $user -p$pass -D $dbname -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf LEFT JOIN feature_set fs ON rf.feature_set_id=fs.feature_set_id ORDER by rf.seq_region_id, rf.seq_region_starti - rf.bound_start_length' > $regulatory_features");
  run("mysql -NB -h $host -u $user -p$pass -D $dbname -e 'select af.seq_region_id, af.seq_region_start, af.seq_region_end, af.annotated_feature_id, fs.cell_type_id FROM annotated_feature af LEFT JOIN feature_set fs ON af.feature_set_id=fs.feature_set_id ORDER BY af.seq_region_id, af.seq_region_start;' > $annotations");
  run("mysql -NB -h $host -u $user -p$pass -D $dbname -e 'select mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id FROM motif_feature mf ORDER BY mf.seq_region_id, mf.seq_region_start;' > $motifs");

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $annotations | awk '\$5==\$10 {print \$4\"\t\"\$9\"\tannotated\"} ' > $out");
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $motifs | awk ' {print \$4\"\t\"\$9\"\tmotif\"} ' >> $out");

  # Load into database
  run("mysql -h $host -u $user -p$pass -D $dbname -e 'LOAD DATA LOCAL INFILE \"$out\" INTO TABLE regulatory_attribute;'");

  unlink $regulatory_features;
  unlink $annotations;
  unlink $motifs;
  unlink $out;
}

1;
