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

perl scripts/sequencing/generate_read_file_report.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./read_file_report





peak_calling_statistics.pl \
    --species          mouse_with_regbuild \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./reports/

perl scripts/segmentation//load_state_emissions_and_assignments.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --emissions_file   /hps/nobackup/production/ensembl/mnuhn/regbuild_small_test_set/segmentation/mus_musculus/learn_model/emissions_25.txt \
    --assignments_file /hps/nobackup/production/ensembl/mnuhn/regbuild_small_test_set/regulatory_build/mus_musculus/tmp/assignments.txt

perl scripts/segmentation/segmentation_statistics.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./reports/

perl scripts/sequencing/generate_read_file_report.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
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

my $experiment_adaptor = $mouse_funcgen_dba->get_ExperimentAdaptor;
my @signal_experiments  = $experiment_adaptor->_fetch_all_signal_experiments;
my @control_experiments = $experiment_adaptor->_fetch_all_control_experiments;

my @all_experiments = $experiment_adaptor->_fetch_all_experiments_with_read_files;

my $total_number_of_reads = 0;
my $total_file_size       = 0;

map {

  $total_number_of_reads += $_->sum_number_of_reads;
  $total_file_size       += $_->sum_read_file_sizes;

} @all_experiments;

# print " total_number_of_reads = $total_number_of_reads \n";
# print " total_file_size = $total_file_size \n";
# print "\n";
# exit;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $description_template = dirname($file) . '/../../templates/read_files/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1, RELATIVE => 1 );

my $output;

my $output_file = "$output_directory/read_file_report.html";

my $dbc = $mouse_funcgen_dba->dbc;

my $sql_1 = qq(
    select 
        logic_name,
        count(experiment_id) as count_experiments
    from 
        experiment 
        join read_file_experimental_configuration using (experiment_id) 
        join read_file using (read_file_id) 
        join analysis using (analysis_id) 
    group by 
        logic_name
);

my $read_file_1_counts = $dbc->db_handle->selectall_hashref($sql_1, 'logic_name');

# my $sql_2 = qq(
#     select 
#         experiment.name                            as experiment_name,
#         idr.max_peaks                              as idr_max_peaks,
#         idr.type                                   as idr_type,
#         idr.failed_idr_pairs                       as idr_failed_idr_pairs,
#         group_concat(distinct analysis.logic_name) as alignment_methods,
#         feature_type.name                          as feature_type_name, 
#         sum(read_file.file_size)                   as num_reads, 
#         avg(read_file.read_length)                 as avg_read_length, 
#         stddev(read_file.read_length)              as stddev_read_length, 
#         max(read_file.read_length)                 as max_read_length, 
#         min(read_file.read_length)                 as min_read_length, 
#         max(read_file.read_length) - min(read_file.read_length) as span
#     from 
#         feature_type 
#         join experiment using (feature_type_id) 
#         join read_file_experimental_configuration using (experiment_id) 
#         join read_file using (read_file_id) 
#         left join idr using (experiment_id)
#         left join alignment_read_file using (read_file_id)
#         left join alignment on (alignment.alignment_id = alignment_read_file.alignment_id and alignment.has_duplicates=true)
#         left join analysis on (alignment.analysis_id = analysis.analysis_id)
#     group by 
#         experiment.name, 
#         feature_type.name
# );

# my $read_file_2_counts = $dbc->db_handle->selectall_hashref($sql_2, 'experiment_name');

# die( Dumper($read_file_1_counts) );
use Number::Format qw( format_number );

my $de = new Number::Format(
    -thousands_sep   => ',',
    -decimal_point   => '.',
    -int_curr_symbol => 'DEM'
);

$tt->process(
    $description_template, 
    {
        signal_experiments  => \@signal_experiments,
        control_experiments => \@control_experiments,

        total_number_of_reads => $total_number_of_reads,
        total_file_size       => $total_file_size,

        fetch_idr => sub  {
          my $experiment = shift;
          return $experiment->_fetch_Idr
        },
        
        time => sub {
          return "" . localtime
        },
        bytes_to_gb => sub {
          my $size_in_bytes = shift;
          return ( 0 + $size_in_bytes ) / (1024 * 1024 * 1024)
        },
        dbc => $dbc,
        read_file_analysis_experiment_counts => $read_file_1_counts,
        
        round_num => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
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


