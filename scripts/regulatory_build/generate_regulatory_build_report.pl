#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

=head1 

  generate_regulatory_build_report.pl \
    --species homo_sapiens \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
    --output_file /homes/mnuhn/public_html/regulatory_build_stats/rb_grch38_testdb16_reporttest/homo_sapiens/regulatory_build.html

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

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry    = " . $registry    . "\n");
$logger->info("species     = " . $species     . "\n");
$logger->info("output_file = " . $output_file . "\n");

my $output_directory = dirname($output_file);

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $previous_version_funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species . '_previous_version', 'funcgen');

my $regulatory_build_adaptor
  = $funcgen_dba
    ->get_RegulatoryBuildAdaptor;

my $regulatory_build_statistics_adaptor 
  = $funcgen_dba
    ->get_RegulatoryBuildStatisticAdaptor;

my $regulatory_build_statistics_adaptor_previous_version 
  = $previous_version_funcgen_dba
    ->get_RegulatoryBuildStatisticAdaptor;

my $file = __FILE__;
use File::Basename qw( dirname basename );

my $template_dir = dirname($file) . '/../../templates/regulatory_build';
my $description_template = $template_dir . '/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

my $genome = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );
my $ref_length = $genome->get_ref_length;

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new(
  ABSOLUTE     => 1,
  RELATIVE     => 1,
  INCLUDE_PATH => $template_dir,
);

my $output;

use Number::Format qw( format_number );

my $de = new Number::Format(
    -thousands_sep   => ',',
    -decimal_point   => '.',
);

$tt->process(
    $description_template, 
    {
        regulatory_build_adaptor 
          => $regulatory_build_adaptor,
        
        regulatory_build_statistics_adaptor 
          => $regulatory_build_statistics_adaptor,
        
        regulatory_build_statistics_adaptor_previous_version 
          => $regulatory_build_statistics_adaptor_previous_version,
        
        dbc        => $funcgen_dba->dbc,
        species    => $species,
        ref_length => $ref_length,
        
        round_num => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
        },
      
        time => sub {
          return "" . localtime
        },
      
      format_number => sub {
        my $number = shift;
        if (! defined $number) {
          return '-'
        }
        if ($number eq '') {
          return '-'
        }
        #return 'foo';
        return $de->format_number($number);
      },


    },
    $output_file
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;


