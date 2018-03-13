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
);
 
my $url          = $options{url};
my $species      = $options{species};
 
Bio::EnsEMBL::Registry->load_registry_from_url($url, 1);
 
my $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $dbc = $db_adaptor->dbc;
 
my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $dbc
);
 
$dbc-> do('update regulatory_feature set stable_id=NULL');

