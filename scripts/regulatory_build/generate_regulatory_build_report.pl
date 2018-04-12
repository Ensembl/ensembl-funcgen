#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

=head1 

time mysql $(r2-w details mysql) mnuhn_testdb5_mus_musculus_funcgen_91_38 -e "drop table peak_calling_statistic"

time mysql $(r2-w details mysql) mnuhn_testdb5_mus_musculus_funcgen_91_38 -e "
    create table 
        peak_calling_statistic as
    select 
        peak_calling_id, 
        sum(peak.seq_region_end - peak.seq_region_start + 1) as total_length,
        count(peak.peak_id) as num_peaks,
        avg(peak.seq_region_end - peak.seq_region_start + 1) as average_length
    from 
        peak_calling 
        left join peak using (peak_calling_id) 
    group by 
        peak_calling_id
    ;
"

perl scripts/regulatory_build/generate_regulatory_build_report.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_dev_break/lib/ensembl-funcgen/registry.with_previous_version.pm \
    --output_directory ./reports/

=cut

use strict;
use Getopt::Long;

my $species;
my $registry;
my $output_directory;

GetOptions (
   'species=s'          => \$species,
   'registry=s'         => \$registry,
   'output_directory=s' => \$output_directory,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry          = " . $registry         . "\n");
$logger->info("species           = " . $species          . "\n");
$logger->info("output_directory  = " . $output_directory . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $mouse_funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $regulatory_build_statistics_adaptor = $mouse_funcgen_dba->get_RegulatoryBuildStatisticAdaptor;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $description_template = dirname($file) . '/../../templates/regulatory_build/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

my $output_file = "$output_directory/regulatory_build_report.html";

my $genome = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );
my $ref_length = $genome->get_ref_length;

# die($ref_length);

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1, RELATIVE => 1 );

my $output;

use Number::Format qw( format_number );

my $de = new Number::Format(
    -thousands_sep   => ',',
    -decimal_point   => '.',
    -int_curr_symbol => 'DEM'
);

$tt->process(
    $description_template, 
    {
        regulatory_build_statistics_adaptor => $regulatory_build_statistics_adaptor,
        dbc        => $mouse_funcgen_dba->dbc,
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


