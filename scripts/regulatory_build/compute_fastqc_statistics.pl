#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

  compute_fastqc_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
    --species homo_sapiens \
    --output_directory /homes/mnuhn/public_html/regulatory_build_stats/rb_grch38_testdb16/homo_sapiens/

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
    "output_directory|o=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species          = $options{'species'};
my $registry         = $options{'registry'};
my $output_directory = $options{'output_directory'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

#my $bins = &generate_bins;

my $read_file_adaptor = $funcgen_adaptor->get_ReadFileAdaptor;
my $all_read_files = $read_file_adaptor->fetch_all;

my $datasets;

my $save_file = 'compute_fastqc_statistics';

if (0) {

  $datasets = compute_datasets($all_read_files);

  open my $out, '>', $save_file;
  $out->print(Dumper($datasets));
  $out->close;

} else {

  local $/;
  open my $in, $save_file;
  
  my $dumped = <$in>;
  
  no strict;
  $datasets = eval $dumped;
  use strict;
  
  $in->close;
}

my $file = __FILE__;
use File::Basename qw( dirname basename );

my $template_dir = dirname($file) . '/../../templates/';
#my $description_template = $template_dir . '/quality_checks/report_fastqc.html';
my $description_template = $template_dir . '/quality_checks/report_fastqc_barcharts.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

#my $output_file = "$output_directory/report_fastqc.html";
my $output_file = "$output_directory/report_fastqc_barcharts.html";
#my $output_file = "/homes/mnuhn/public_html/regulatory_build_stats/rb_grch38_testdb16/homo_sapiens/report_fastqc_barcharts.html";


use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new(
  ABSOLUTE     => 1,
  RELATIVE     => 1,
  INCLUDE_PATH => $template_dir,
);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

$tt->process(
    $description_template, 
    {
        title    => 'FastQC',
#        labels   => $bins,
        
#        distribution_y_axis_max => 350,
        
        datasets => $datasets,

        dbc        => $funcgen_dba->dbc,
        species    => $species,
        time => sub {
          return "" . localtime
        },


    },
    $output_file
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;
exit(0);

sub create_read_file_dataset {

  my $read_files                = shift;
  my $initial_dataset           = shift;
  my $read_file_filter_callback = shift;

  my $filtered_peak_callings = [ grep { $read_file_filter_callback->($_) } @$read_files ];
  
  my $final_dataset = {
      %$initial_dataset,
      fastqc_outcome => fetch_fastqc_counts($filtered_peak_callings),
  };
  
  return $final_dataset;
}

sub fetch_fastqc_counts {

  my $read_files = shift;
  
  my $pass_fail_warn_hash_factory = sub {
    return {
      'PASS' => 0,
      'FAIL' => 0,
      'WARN' => 0,
    };
  };
  
  my %fastqc_outcome_count = (

    basic_statistics             => $pass_fail_warn_hash_factory->(),
    per_base_sequence_quality    => $pass_fail_warn_hash_factory->(),
    per_tile_sequence_quality    => $pass_fail_warn_hash_factory->(),
    per_sequence_quality_scores  => $pass_fail_warn_hash_factory->(),
    per_base_sequence_content    => $pass_fail_warn_hash_factory->(),
    per_sequence_gc_content      => $pass_fail_warn_hash_factory->(),
    per_base_n_content           => $pass_fail_warn_hash_factory->(),
    sequence_length_distribution => $pass_fail_warn_hash_factory->(),
    sequence_duplication_levels  => $pass_fail_warn_hash_factory->(),
    overrepresented_sequences    => $pass_fail_warn_hash_factory->(),
    adapter_content              => $pass_fail_warn_hash_factory->(),
    kmer_content                 => $pass_fail_warn_hash_factory->(),

    run_failed                   => 0,

  );

  process_read_file_qc_values_from_peak_callings(
  
    $read_files,
    
    sub {
    
      my $fastqc = shift;
      
      my $basic_statistics              = $fastqc->basic_statistics;
      my $per_base_sequence_quality     = $fastqc->per_base_sequence_quality;
      my $per_tile_sequence_quality     = $fastqc->per_tile_sequence_quality;
      my $per_sequence_quality_scores   = $fastqc->per_sequence_quality_scores;
      my $per_base_sequence_content     = $fastqc->per_base_sequence_content;
      my $per_sequence_gc_content       = $fastqc->per_sequence_gc_content;
      my $per_base_n_content            = $fastqc->per_base_n_content;
      my $sequence_length_distribution  = $fastqc->sequence_length_distribution;
      my $sequence_duplication_levels   = $fastqc->sequence_duplication_levels;
      my $overrepresented_sequences     = $fastqc->overrepresented_sequences;
      my $adapter_content               = $fastqc->adapter_content;
      my $kmer_content                  = $fastqc->kmer_content;
      my $run_failed                    = $fastqc->run_failed;

      #$fastqc_outcome_count{$basic_statistics}++;

      $fastqc_outcome_count{basic_statistics}              {$basic_statistics}++;
      $fastqc_outcome_count{per_base_sequence_quality}     {$per_base_sequence_quality}++;
      $fastqc_outcome_count{per_tile_sequence_quality}     {$per_tile_sequence_quality}++;
      $fastqc_outcome_count{per_sequence_quality_scores}   {$per_sequence_quality_scores}++;
      $fastqc_outcome_count{per_base_sequence_content}     {$per_base_sequence_content}++;
      $fastqc_outcome_count{per_sequence_gc_content}       {$per_sequence_gc_content}++;
      $fastqc_outcome_count{per_base_n_content}            {$per_base_n_content}++;  
      $fastqc_outcome_count{sequence_length_distribution}  {$sequence_length_distribution}++;
      $fastqc_outcome_count{sequence_duplication_levels}   {$sequence_duplication_levels}++;
      $fastqc_outcome_count{overrepresented_sequences}     {$overrepresented_sequences}++;
      $fastqc_outcome_count{adapter_content}               {$adapter_content}++;
      $fastqc_outcome_count{kmer_content}                  {$kmer_content}++;
      
      if ($run_failed) {
        $fastqc_outcome_count{run_failed}++;
      }
    },
    
    sub {
      $fastqc_outcome_count{run_failed}++;
    }
  );
  return \%fastqc_outcome_count
}

sub process_read_file_qc_values_from_peak_callings {

  my $read_files                    = shift;
  my $read_file_callback            = shift;
  my $read_file_run_failed_callback = shift;

  PEAK_CALLING:
  foreach my $read_file (@$read_files) {

    my $fastqc = $read_file->fetch_FastQC;
    
    if ($fastqc) {
      $read_file_callback->($fastqc);
      next PEAK_CALLING;
    }
    $read_file_run_failed_callback->($fastqc);
  }
  return;
}

sub compute_datasets {

  my $all_read_files = shift;
  my @datasets;

  push 
    @datasets, 
    {
      title => 'All',
      all => create_read_file_dataset(

            $all_read_files,
            
            {
              title   => 'All',
              html_id => 'overall-all', 
              colour  => 'window.chartColors.gray',
            },

            sub { return 1; }
        ),
      narrow => create_read_file_dataset(

            $all_read_files,
            
            {
              title   => 'All narrow peaks',
              html_id => 'overall-narrow', 
              colour  => 'window.chartColors.gray',
            },

            sub { return (! shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks) }
        ),
      broad => create_read_file_dataset(

            $all_read_files,
            
            {
              title   => 'All broad peaks',
              html_id => 'overall-broad', 
              colour  => 'window.chartColors.gray',
            },

            sub { return shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks }
        ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Blueprint',
      all => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Blueprint',
            html_id => 'blueprint-all',
            colour  => 'window.chartColors.red',
          },

        sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      ),
      narrow => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Blueprint narrow peaks',
            html_id => 'blueprint-narrow',
            colour  => 'window.chartColors.red',
          },

        sub { 
          my $read_file = shift;

                  $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
            && (! $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks)
        }
      ),
      broad => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Blueprint broad peaks',
            html_id => 'blueprint-broad',
            colour  => 'window.chartColors.red',
          },

        sub { 
          my $read_file = shift;

              $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
            && $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks
        }
      ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'ENCODE',
      all => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'ENCODE',
            html_id => 'encode-all',
            colour  => 'window.chartColors.blue',
          },

          sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      ),
      narrow => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'ENCODE narrow peaks',
            html_id => 'encode-narrow',
            colour  => 'window.chartColors.blue',
          },

          sub { 
            my $read_file = shift;
            
                  $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
            && (! $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks)
          }
      ),
      broad => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'ENCODE broad peaks',
            html_id => 'encode-broad',
            colour  => 'window.chartColors.blue',
          },

          sub { 
            my $read_file = shift;
            
              $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
            && $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks
          }
      ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Roadmap Epigenomics',
      all => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Roadmap Epigenomics',
            html_id => 'roadmap_epigenomics',
            colour  => 'window.chartColors.green',
          },

          sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      ),
      narrow => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Roadmap Epigenomics narrow peaks',
            html_id => 'roadmap_epigenomics-narrow',
            colour  => 'window.chartColors.green',
          },

          sub { 
            my $read_file = shift;
            
                  $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
            && (! $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks)
          }
      ),
      broad => create_read_file_dataset(

          $all_read_files,
          
          {
            title   => 'Roadmap Epigenomics broad peaks',
            html_id => 'roadmap_epigenomics-broad',
            colour  => 'window.chartColors.green',
          },

          sub { 
            my $read_file = shift;
            
              $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
            && $read_file->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_FeatureType->creates_broad_peaks
          }
      ),
    }
  ;
  return \@datasets;
}

