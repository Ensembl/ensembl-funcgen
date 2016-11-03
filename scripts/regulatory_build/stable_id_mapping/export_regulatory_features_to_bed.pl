use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;

my %options;

GetOptions (
  \%options,
  "url=s",
  "species=s",
  "chromosome=s",
  "stable_id_prefix=s",
  "outfile=s",
);

my $url              = $options{url};
my $species          = $options{species};
my $chromosome       = $options{chromosome};
my $stable_id_prefix = $options{stable_id_prefix};
my $outfile          = $options{outfile};

Bio::EnsEMBL::Registry->load_registry_from_url($url, 1);

my $regulatory_build_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryBuild');
my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryFeature');

my $current_regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

open my $out_fh, '>' . $outfile;

my $regulatory_feature_iterator = $regulatory_feature_adaptor->fetch_Iterator;

while (my $current_regulatory_feature = $regulatory_feature_iterator->next) {

  my $bed_string = join("\t", (
    $current_regulatory_feature->seq_region_name,
    $current_regulatory_feature->bound_start,
    $current_regulatory_feature->bound_end,
    $current_regulatory_feature->feature_type->name,
    remove_stable_id_prefix($stable_id_prefix, $current_regulatory_feature->stable_id),
    $current_regulatory_feature->dbID
  ));
  $out_fh->print($bed_string . "\n");
}
$out_fh->close;

sub remove_stable_id_prefix {

  my $stable_id_prefix = shift;
  my $stable_id = shift;
  
  $stable_id =~ s/^$stable_id_prefix//;
  return $stable_id;
}
