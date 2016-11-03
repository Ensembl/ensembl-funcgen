use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;

my $species;
my $url;
my $segmentation_file_data;

GetOptions (
   'url=s'                    => \$url,
   'species=s'                => \$species,
   'segmentation_file_data=s' => \$segmentation_file_data,
);

Bio::EnsEMBL::Registry->load_registry_from_url($url, 1);

my $segmentation_logic_name   = 'Segmentation';
my $segmentation_display_name = 'Segmentation';
my $segmentation_description  = 'Segmentation of an epigenome';

my $segmentation_file_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'SegmentationFile');
my $analysis_adaptor          = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Analysis');
my $epigenome_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Epigenome');
my $regulatory_build_adaptor  = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryBuild');

my $analysis         = $analysis_adaptor->fetch_by_logic_name($segmentation_logic_name);
my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

open my $in, '<', $segmentation_file_data;

my @segmentation_file_list;
while (my $current_line = <$in>) {

  chomp $current_line;

  (my $epigenome_name, my $segmentation_file_name) = split "\t", $current_line;
  
  push @segmentation_file_list, {
    epigenome_name => $epigenome_name, 
    segmentation_file_name => $segmentation_file_name,
  };
}

$in->close;

if (! defined $analysis) {

  print "Not in database, so creating\n";

  use Bio::EnsEMBL::Analysis;
  $analysis = Bio::EnsEMBL::Analysis->new(
    -logic_name    => $segmentation_logic_name,
    -description   => $segmentation_description,
    -display_label => $segmentation_display_name,
    -displayable   => 1,
  );
  $analysis_adaptor->store($analysis);

} else {
  print "Was in database, so using existing object.\n";
}

use Bio::EnsEMBL::Funcgen::SegmentationFile;

# Check all epigenomes are in the database
#
foreach my $current_item (@segmentation_file_list) {

  my $epigenome_name         = $current_item->{epigenome_name};
  my $segmentation_file_name = $current_item->{segmentation_file_name};
  
  my $epigenome = $epigenome_adaptor->fetch_by_name($epigenome_name);
  
  if (! defined $epigenome) {
    die("Can't find epigenome with name $epigenome_name!");
  }
}

foreach my $current_item (@segmentation_file_list) {

  my $epigenome_name         = $current_item->{epigenome_name};
  my $segmentation_file_name = $current_item->{segmentation_file_name};
  
  my $epigenome = $epigenome_adaptor->fetch_by_name($epigenome_name);
  
  die if (! defined $epigenome);
  
  my $segmentation_file = Bio::EnsEMBL::Funcgen::SegmentationFile->new(
    -name             => 'Segmentation of ' . $epigenome->display_label,
    -analysis         => $analysis,
    -epigenome        => $epigenome,
    -regulatory_build => $regulatory_build,
    -file             => $segmentation_file_name,
    -file_type        => 'BIGBED',
  );

  print "$epigenome_name\t$segmentation_file_name\t" . $epigenome->description;
  print "\n";
  
  $segmentation_file_adaptor->store($segmentation_file);
}
