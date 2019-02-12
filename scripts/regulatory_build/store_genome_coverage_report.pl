#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use v5.10;
use Data::Dumper;
use Carp;
$Data::Dumper::Maxdepth = 3;

=head1 

store_genome_coverage_report.pl \
    --species homo_sapiens \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb12.pm \
    --genome_coverage_report genome_coverage_report_2.pl

=cut 

use strict;
use Getopt::Long;

my $species;
my $registry;
my $genome_coverage_report_file;

GetOptions (
   'species=s'                => \$species,
   'registry=s'               => \$registry,
   'genome_coverage_report=s' => \$genome_coverage_report_file,
);

use Bio::EnsEMBL::Utils::Logger;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry               = " . $registry              . "\n");
$logger->info("species                = " . $species               . "\n");
$logger->info("genome_coverage_report_file    = " . $genome_coverage_report_file   . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $dbc = $funcgen_dba->dbc;

my $regulatory_build_adaptor = $funcgen_dba->get_RegulatoryBuildAdaptor;
my $current_regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

my $regulatory_build_statistics_adaptor = $funcgen_dba->get_RegulatoryBuildStatisticAdaptor;

my $parsed_statistics = load_data_from_file($genome_coverage_report_file);

use Hash::Util qw( lock_keys );

lock_keys(%$parsed_statistics);

print Dumper($parsed_statistics);

my $statistics = 
    [
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'open_chromatin_overlap_percent',
            -value               => $parsed_statistics->{'Open chromatin'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'open_chromatin_overlap_bp',
            -value               => $parsed_statistics->{'Open chromatin'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'promoter_flanking_overlap_percent',
            -value               => $parsed_statistics->{'Promoter Flanking Region'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'promoter_flanking_overlap_bp',
            -value               => $parsed_statistics->{'Promoter Flanking Region'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'promoter_overlap_percent',
            -value               => $parsed_statistics->{'Promoter'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'promoter_overlap_bp',
            -value               => $parsed_statistics->{'Promoter'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'tf_binding_overlap_percent',
            -value               => $parsed_statistics->{'TF binding site'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'tf_binding_overlap_bp',
            -value               => $parsed_statistics->{'TF binding site'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'enhancer_overlap_percent',
            -value               => $parsed_statistics->{'Enhancer'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'enhancer_overlap_bp',
            -value               => $parsed_statistics->{'Enhancer'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'ctcf_overlap_percent',
            -value               => $parsed_statistics->{'CTCF Binding Site'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'ctcf_overlap_bp',
            -value               => $parsed_statistics->{'CTCF Binding Site'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'regulatory_build_overlap_percent',
            -value               => $parsed_statistics->{'regulatory_build'}->{percent_overlap},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
        Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
            -statistic           => 'regulatory_build_overlap_bp',
            -value               => $parsed_statistics->{'regulatory_build'}->{total_overlap_size},
            -regulatory_build_id => $current_regulatory_build->dbID,
        ),
    ];

print Dumper($statistics);

foreach my $statistic (@$statistics) {

    my $statistic_name = $statistic->statistic;
    
    print "Storing $statistic_name\n";
    
    $dbc->do("delete from regulatory_build_statistic where statistic = '$statistic_name';");
    $regulatory_build_statistics_adaptor->store($statistic);
}

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
