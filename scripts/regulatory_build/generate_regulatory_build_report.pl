#!/usr/bin/env perl

=head1 

  generate_regulatory_build_report.pl \
    --species homo_sapiens \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
    --output_file ./test/regulatory_build.html

=cut

use strict;
use Getopt::Long;

my $species;
my $registry;
my $output_file;

GetOptions (
   'species=s'     => \$species,
   'registry=s'    => \$registry,
   'output_file=s' => \$output_file,
);

use Bio::EnsEMBL::Funcgen::Report::RegulatoryBuild;
my $regulatory_build_report = Bio::EnsEMBL::Funcgen::Report::RegulatoryBuild->new(
  -species      => $species,
  -registry     => $registry,
  -output_file  => $output_file,
);

$regulatory_build_report->generate_report;
