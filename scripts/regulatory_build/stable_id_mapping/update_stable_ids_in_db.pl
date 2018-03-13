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
  "mapping_file=s",
);

my $url          = $options{url};
my $species      = $options{species};
my $mapping_file = $options{mapping_file};

Bio::EnsEMBL::Registry->load_registry_from_url($url, 1);

my $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $dbc = $db_adaptor->dbc;

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $dbc
);

open my $IN, $mapping_file;

$helper->batch(
  -SQL      => 'update regulatory_feature set stable_id=? where regulatory_feature_id=?',
  -CALLBACk => sub {
    my ( $sth, $dbc ) = @_;
    while (my $current_line = <$IN>) {
      chomp ($current_line);
      my @fields = split "\t", $current_line;
      
      my $regulatory_feature_id = $fields[0];
      my $stable_id             = $fields[1];
      
      $sth->execute( $stable_id, $regulatory_feature_id );
    }
  }
);
