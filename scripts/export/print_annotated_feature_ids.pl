#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

=head1

print_annotated_feature_ids.pl \
  --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --epigenome_production_name H1ESC \
  --feature_type_name H3K27ac

=cut

my $registry;
my $species;
my $epigenome_production_name;
my $feature_type_name;

GetOptions (
   'registry=s'                   => \$registry,
   'species=s'                    => \$species,
   'epigenome_production_name=s'  => \$epigenome_production_name,
   'feature_type_name=s'          => \$feature_type_name,
);

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'funcgen' );

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

$helper->execute_no_return(
  -SQL      => "select annotated_feature.annotated_feature_id from annotated_feature join feature_set using (feature_set_id) join epigenome using (epigenome_id) join feature_type using (feature_type_id) where epigenome.production_name=? and feature_type.name=?",
  -PARAMS => [ $epigenome_production_name, $feature_type_name ],
  -CALLBACK => sub {
    my $row = shift;
    my $annotated_feature_id = $row->[0];
    print $annotated_feature_id;
    print "\n";
    return;
  },
);
