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

load_segmentation.pl

=head1 SYNOPSIS

load_segmentations.pl --host $host --user $user --pass $pass --base_dir hg38 --dbname homo_sapiens_funcgen_76_38 --segmentation $segmentation_name --cell_type $cell_type_name

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
use Bio::EnsEMBL::Utils::Exception         qw( throw stack_trace );

main();

sub main {
  my ($base_dir, $segmentation, $cell_type, $pass,$port,$host,$user,$dbname);

  GetOptions 
    (
     "base_dir|b=s"       => \$base_dir,
     "segmentation=s"       => \$segmentation,
     "cell_type=s"       => \$cell_type,
     "pass|p=s"       => \$pass,
     "port=s"         => \$port,
     "host|h=s"       => \$host,
     "user|u=s"       => \$user,
     "dbname|d=s"     => \$dbname,
     'help|?'         => sub { pos2usage(-exitval => 0, ); }
    ) or pod2usage( -exitval => 1);

  my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $host,
     -user   => $user,
     -dbname => $dbname,
     -pass   => $pass,
     -port   => $port,
     -group  => 'funcgen',#Should be set as default in adaptor new method
    );

#Test connections
  $db->dbc->db_handle;

#Test species
  my $species = $db->species;

  if(! defined $species){
    die("Could not get a valid species from $dbname, please check the meta species.production_name");
  }

  defined $base_dir || die ("You must define the base directory!\t--base_dir XXXX\n");
  defined $segmentation || die ("You must define the segmentation name!\t--segmentation XXXX\n");
  defined $cell_type || die ("You must define the cell type name!\t--cell_type XXXX\n");

  my $clean_cell_type = clean_name($cell_type);

  load_segmentation_features_from_file($db, $segmentation, $cell_type, "$base_dir/segmentations/$segmentation/$clean_cell_type.bb");
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
# Convenience wrapper to run the commandline safely 
#####################################################
# Params: command line command string
#####################################################

sub run {
  my ($cmd) = @_;
  print cmd;
  system($cmd) && die("Failed when running command:\n$cmd\n");
}

sub load_segmentation_features_from_file {
  my ($db, $segmentation, $cell_type, $file) = @_;
  print "Processing file $file\n";
  my $analysis = create_analysis($db, $segmentation, $cell_type);
  my $feature_set = create_feature_set($db, $segmentation, $cell_type, $analysis);
  my $feature_type = create_feature_type($db, $analysis);
  my $sa = $db->get_SliceAdaptor();
  my $sfa = $db->get_SegmentationFeatureAdaptor();
  my ($fh, $tmp) = tempfile();
  my $prev_chrom = undef;
  my $slice = undef;
  run("bigBedToBed $file $tmp");
  my @features = ();
  my $previous_label = undef;
  my $previous_seg_feat = undef;
  while (my $line = <$fh>) {
    chomp $line;
    my ($chr, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $rgb) = split /\t/, $line;
    my ($state, $label, $number) = split /_/, $name;
    defined $feature_type->{$label} || die("$label has no associated feature type!");

    if ($label eq $previous_label && $chr eq $prev_chrom) {
      $previous_seg_feat->end($end);
    } else {
      if ($chr ne $prev_chrom) {
        $slice = $sa->fetch_by_region('toplevel', $chr);
        $prev_chrom = $chr;
      }

      if (scalar(@features) == 10000) {
	print "FLUSH\n";
        $sfa->store(@features);
        @features = ();
      }

      if (!defined $feature_type->{$label}) {
        die("Unknown label $label\n");
      }

      $previous_seg_feat = Bio::EnsEMBL::Funcgen::SegmentationFeature->new
       (
	-SLICE         => $slice,
	-START         => $start + 1,
	-END           => $end,
	-STRAND        => 0,
	-FEATURE_SET   => $feature_set,
	-FEATURE_TYPE  => $feature_type->{$label},
      );
      $previous_label = $label;
      push @features, $previous_seg_feat;
    }
  }

  if (scalar(@features) > 0) {
    $sfa->store(@features);
  }
  close $fh;
  unlink $tmp;

}

sub create_analysis {
### Check whether analysis is already stored
  my ($db, $segmentation, $cell_type) = @_;

  my $logic_name = "Segmentation.$cell_type.$segmentation";
  my $aa  = $db->get_AnalysisAdaptor();
  my $ana = $aa->fetch_by_logic_name($logic_name);

  if ( defined $ana ) {
    return $ana;
  }

  $ana = Bio::EnsEMBL::Analysis->new
    (
     -logic_name      => $logic_name,
     -db              => 'NULL',
     -db_version      => 'NULL',
     -db_file         => 'NULL',
     -program         => "chromHMM/Segway",
     -program_version => 'NULL',
     -program_file    => 'NULL', #Could use $program_name here, but this is only part of the build
     -gff_source      => 'NULL',
     -gff_feature     => 'NULL',
     -module          => 'NULL',
     -module_version  => 'NULL',
     -parameters      => 'NULL',
     -created         => 'NULL',
     -display_label   => "$segmentation Segmentation for cell type $cell_type",
     -displayable     => 1,
    );

  $aa->store($ana);
  return $ana;
}

sub rollback_segmentation_FeatureSet {
  my ($db, $fset) = @_;
  $db->is_stored_and_valid( 'Bio::EnsEMBL::Funcgen::FeatureSet', $fset );
  my $table = 'segmentation_feature';
  warn("Rolling back $fset->feature_class FeatureSet:\t".$fset->name."\n");

  #Check whether this is a supporting set for another data_set
  my @dsets =  @{ $db->get_DataSetAdaptor->fetch_all_by_supporting_set($fset) };

  if (@dsets) {
    my $txt = $fset->name." is a supporting set of the following DataSets:\t" .
              join( ', ', ( map { $_->name } @dsets ) );
    throw( $txt ."\nPlease resolve by deleting dependant Feature/DataSets" ); 
    #This will currently not allow rollback if we have already associated it 
    #with a regulatory build
  }
  
  # Now do the rollback 
  $fset->adaptor->revoke_states($fset);
      
  #Remove object_xref records (but not xref which may be used by soemthing else)
  my $sql = "DELETE ox from object_xref ox, $table f where ox.ensembl_object_type='"
    .ucfirst( $fset->feature_class )."Feature' and ox.ensembl_id=f.${table}_id and ".
    "f.feature_set_id=". $fset->dbID;
  $db->rollback_table( $sql, 'object_xref', 'object_xref_id' );
    
  #Remove associated_feature_type records
  $sql = "DELETE aft from associated_feature_type aft, $table f where ".
    "f.feature_set_id=".$fset->dbID." and f.${table}_id=aft.table_id and ".
    "aft.table_name='".$fset->feature_class . "_feature'";
  $db->rollback_table( $sql, 'associated_feature_type' );
    
  #Remove features
  $sql = "DELETE f from $table f where f.feature_set_id=".$fset->dbID;
    
  warn $sql;  
  $db->rollback_table( $sql, $table, "${table}_id" );
}

sub create_feature_set {
  my ($db, $segmentation, $cell_type, $analysis) = @_;
  my $fsa = $db->get_FeatureSetAdaptor();
  my $fta = $db->get_FeatureTypeAdaptor();
  my $cta = $db->get_CellTypeAdaptor();
  my $name = "$cell_type $segmentation segmentation";
  my $previous = $fsa->fetch_by_name($name);
  my $ft = $fta->fetch_by_name('SegmentationState');
  if (defined $previous) {
    rollback_segmentation_FeatureSet($db, $previous);
    warn "Done rolling back";
    return $previous;
  }
  if (!defined $cell_type) {
    die "Undefined cell_type";
  }
  if (!defined $name) {
    die "Undefined name";
  }
  if (!defined $ft) {
    $ft = Bio::EnsEMBL::Funcgen::FeatureType->new(
      -name  => 'SegmentationState',
      -class => "Segmentation State",
      -description => 'Regulatory genome segment'
    );
    $fta->store($ft);
  }
  if (!defined $ft) {
    die "Undefined feature type with name 'SegmentationState'";
  }
  my $feature_set = Bio::EnsEMBL::Funcgen::FeatureSet->new(
    -name          => $name,
    -analysis      => $analysis,
    -feature_type  => $ft,
    -feature_class => 'segmentation',
    -cell_type     => $cta->fetch_by_name($cell_type),
    -description   => "$cell_type $segmentation segmentation",
    -display_label => "$cell_type segmentation",
    -states        => ['DISPLAYABLE','IMPORTED','IMPORTED_GRCh38','MART_DISPLAYABLE'],
  );
  $fsa->store($feature_set);
  return $feature_set;
}

sub create_feature_type {
  my ($db, $analysis) = @_;

  my %name = (
    tss => "Predicted Promoter with TSS",
    proximal => "Predicted Promoter Flank",
    distal => "Predicted Enhancer",
    ctcf => "CTCF enriched",
    dead => "Predicted heterochomatin",
    weak => "Predicted low activity",
    gene => "Predicted Transcribed Region",
    repressed => "Predicted Repressed",
    poised => "Predicted Poised"
  );

  my %description = (
    tss => "Predicted promoter region including transcription start site",
    proximal => "Predicted promoter flanking region",
    distal => "Predicted enhancer",
    ctcf => "CTCF enriched element",
    dead => "Predicted heterochromatin",
    weak => "Predicted low activity proximal to active states",
    gene => "Predicted transcribed region",
    repressed => "Predicted repressed region",
    poised => "Predicted poised region"
  );

  my $fta = $db->get_FeatureTypeAdaptor();
  my %existing_feature_type = ();
  foreach my $ft (@{$fta->fetch_all_by_Analysis($analysis)}) {
    $existing_feature_type{$ft->name} = $ft;
  }

  my %feature_type = ();
  foreach my $key (keys %name) {
    if (defined $existing_feature_type{$name{$key}}) {
      $feature_type{$key} = $existing_feature_type{$name{$key}};
    } else {
      my $ft = Bio::EnsEMBL::Funcgen::FeatureType->new
      (
        -name  => $name{$key},
        -class => "Segmentation State",
        -description => $description{$key},
        -analysis => $analysis,
      );
      $fta->store($ft);
      $feature_type{$key} = $ft;
    }

    if (!defined $feature_type{$key}) {
      die("Failed to create feature type for label $key\n");
    }
    if (!defined $feature_type{$key}->dbID) {
      die("Failed to store feature type for label $key\n");
    }
  }
  return \%feature_type;
}

1;
