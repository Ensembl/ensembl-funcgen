#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use v5.10;
use Data::Dumper;
use Carp;
$Data::Dumper::Maxdepth = 3;

=head1 

store_stable_id_mapping_statistics.pl \
    --species             mus_musculus \
    --registry            /homes/mnuhn/work_dir_dev_break/lib/ensembl-funcgen/registry.with_previous_version.pm \
    --mapping_report_file /nfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/temp_dir/regulatory_build/mus_musculus/stable_id_mapping/mapping_report.pl

=cut 

use strict;
use Getopt::Long;

my $species;
my $registry;
my $mapping_report_file;

GetOptions (
   'species=s'             => \$species,
   'registry=s'            => \$registry,
   'mapping_report_file=s' => \$mapping_report_file,
);

use Bio::EnsEMBL::Utils::Logger;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry               = " . $registry              . "\n");
$logger->info("species                = " . $species               . "\n");
$logger->info("mapping_report_file    = " . $mapping_report_file   . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $dbc = $funcgen_dba->dbc;

$dbc->do("delete from regulatory_build_statistic where statistic like '%mapping%';");

my $regulatory_build_adaptor = $funcgen_dba->get_RegulatoryBuildAdaptor;
my $current_regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

my $regulatory_build_statistics_adaptor = $funcgen_dba->get_RegulatoryBuildStatisticAdaptor;

my $parsed_emission = load_data_from_file($mapping_report_file);

use Hash::Util qw( lock_keys );

lock_keys(%$parsed_emission);

my @stable_id_mapping_report_keys = qw(
    total_number_regulatory_features
    new_stable_ids
    mapped_stable_ids
);

$regulatory_build_statistics_adaptor->store(
    Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
        -statistic           => 'stable_id_mapping_number_regulatory_features',
        -value               => $parsed_emission->{total_number_regulatory_features},
        -regulatory_build_id => $current_regulatory_build->dbID,
    )
);
$regulatory_build_statistics_adaptor->store(
    Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
        -statistic           => 'stable_id_mapping_new_stable_ids',
        -value               => $parsed_emission->{new_stable_ids},
        -regulatory_build_id => $current_regulatory_build->dbID,
    )
);
$regulatory_build_statistics_adaptor->store(
    Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
        -statistic           => 'stable_id_mapping_mapped_stable_ids',
        -value               => $parsed_emission->{mapped_stable_ids},
        -regulatory_build_id => $current_regulatory_build->dbID,
    )
);

print Dumper($parsed_emission);

sub load_data_from_file {

  my $file = shift;
  
  open my $fh, '<', $file or confess("Can't open file for reading $file");
  local $/;
  
  my $serialised = <$fh>;
  $fh->close;
  no strict;
  my $data = eval $serialised;
  use strict;
  
  return $data;
}
