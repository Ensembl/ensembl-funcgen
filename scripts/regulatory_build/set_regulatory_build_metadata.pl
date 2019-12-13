#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

  A script to compute the number of bases covered by regulatory features both in total and by feature type

set_regulatory_build_metadata.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb13.pm \
    --species homo_sapiens \
    --name "The new name" \
    --description "Superdescription" \
    --release_version 095

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "name=s",
    "description=s",
    "release_version=s",
    "sample_regulatory_feature_id=s",
 );

my $species                      = $options{'species'};
my $registry                     = $options{'registry'};
my $name                         = $options{'name'};
my $description                  = $options{'description'};
my $release_version              = $options{'release_version'};
my $sample_regulatory_feature_id = $options{'sample_regulatory_feature_id'};

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $funcgen_adaptor          = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $regulatory_build_adaptor = $funcgen_adaptor->get_RegulatoryBuildAdaptor;
my $regulatory_build         = $regulatory_build_adaptor->fetch_current_regulatory_build;

$logger->info("Regulatory build before updating:\n");
$logger->info(Dumper($regulatory_build));

if ($name) {

  $logger->info("New name: $name\n");
  $regulatory_build->name($name);

}

if ($description) {

  $logger->info("New description: $description\n");
  $regulatory_build->description($description);

}

if ($release_version) {

  $logger->info("New release_version: $release_version\n");
  $regulatory_build->release_version($release_version);

}

if ($sample_regulatory_feature_id) {

  $logger->info("New sample_regulatory_feature_id: $sample_regulatory_feature_id\n");
  $regulatory_build->sample_regulatory_feature_id($sample_regulatory_feature_id);

}

$logger->info("Regulatory build after updating:\n");
$logger->info(Dumper($regulatory_build));

$regulatory_build_adaptor->update($regulatory_build);

$logger->finish_log;
