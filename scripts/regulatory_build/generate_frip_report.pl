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

  generate_frip_report.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
    --species homo_sapiens \
    --output_file ./test/frip.html

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
    "output_file|o=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species          = $options{'species'};
my $registry         = $options{'registry'};
my $output_file      = $options{'output_file'};

my $output_directory = dirname($output_file);

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $bins = &generate_bins;

my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
my $all_peak_callings = $peak_calling_adaptor->fetch_all;

my @datasets;

push 
  @datasets, 
  {
    title => 'All',
    all => create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'All',
            colour  => 'window.chartColors.gray',
          },

          sub { return 1; }
      ),
    narrow => create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'All narrow peaks',
            colour  => 'window.chartColors.gray',
          },

          sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) }
      ),
    broad => create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'All broad peaks',
            colour  => 'window.chartColors.gray',
          },

          sub { return shift->fetch_Experiment->get_FeatureType->creates_broad_peaks }
      ),
  }
;

push 
  @datasets, 
  {
    title => 'Blueprint',
    all => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Blueprint',
          colour  => 'window.chartColors.red',
        },

      sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
    ),
    narrow => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Blueprint narrow peaks',
          colour  => 'window.chartColors.red',
        },

      sub { 
        my $peak_calling = shift;

                $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
          && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
      }
    ),
    broad => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Blueprint broad peaks',
          colour  => 'window.chartColors.red',
        },

      sub { 
        my $peak_calling = shift;

             $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
          && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
      }
    ),
  }
;

push 
  @datasets, 
  {
    title => 'ENCODE',
    all => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'ENCODE',
          colour  => 'window.chartColors.blue',
        },

        sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
    ),
    narrow => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'ENCODE narrow peaks',
          colour  => 'window.chartColors.blue',
        },

        sub { 
          my $peak_calling = shift;
          
                $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
          && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
        }
    ),
    broad => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'ENCODE broad peaks',
          colour  => 'window.chartColors.blue',
        },

        sub { 
          my $peak_calling = shift;
          
             $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
          && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
        }
    ),
  }
;

push 
  @datasets, 
  {
    title => 'Roadmap Epigenomics',
    all => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Roadmap Epigenomics',
          colour  => 'window.chartColors.green',
        },

        sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
    ),
    narrow => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Roadmap Epigenomics narrow peaks',
          colour  => 'window.chartColors.green',
        },

        sub { 
          my $peak_calling = shift;
          
                $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
          && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
        }
    ),
    broad => create_phantom_peak_dataset(

        $all_peak_callings,
        
        {
          title   => 'Roadmap Epigenomics broad peaks',
          colour  => 'window.chartColors.green',
        },

        sub { 
          my $peak_calling = shift;
          
             $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
          && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
        }
    ),
  }
;



my $file = __FILE__;
use File::Basename qw( dirname basename );

my $template_dir = dirname($file) . '/../../templates/';
my $description_template = $template_dir . '/quality_checks/frip.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

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
        title    => 'Fraction of reads in peaks',
        labels   => $bins,
        
        distribution_y_axis_max => 250,
        
        datasets => \@datasets,

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

sub create_phantom_peak_dataset {

  my $peak_callings                = shift;
  my $initial_dataset              = shift;
  my $peak_calling_filter_callback = shift;

  my $filtered_peak_callings = [ grep { $peak_calling_filter_callback->($_) } @$peak_callings ];
  
  my $final_dataset = {
      %$initial_dataset,
      counts             => count_frip_values_in_bins($filtered_peak_callings, $bins),
  };
  
  return $final_dataset;
}

sub generate_bins {

  my $min = 0;
  my $max = 1;
  my $num_bins = 50;

  my @bins
    =
      map { $_ / $num_bins * ($max - $min) + $min }
      0..$num_bins
    ;
  return \@bins
}

sub count_frip_values_in_bins {

  my $peak_callings = shift;
  my $bins          = shift;

  my $values = fetch_frip_values($peak_callings);
  my $counts = compute_bin_counts($bins, $values);
  
  return $counts;
}

sub compute_bin_counts {

  my $bins   = shift;
  my $values = shift;

  use Statistics::Descriptive;
  my $stat = Statistics::Descriptive::Full->new();

  $stat->add_data( @$values );
  my $f = $stat->frequency_distribution_ref( $bins );

  my @counts;
  for (sort {$a <=> $b} keys %$f) {
    push @counts, $f->{$_};
  }
  return \@counts;
}

sub fetch_frip_values {

  my $peak_callings = shift;
  
  my @frip_values;
  process_frip_qc_values_from_peak_callings(
  
    $peak_callings,
    
    sub {
      my $frip = shift;
      my $frip_value = $frip->frip;
      
      push @frip_values, $frip_value;
    },
    
    sub {}
  );
  return \@frip_values;
}

sub process_frip_qc_values_from_peak_callings {

  my $peak_callings            = shift;
  my $frip_callback            = shift;
  my $frip_run_failed_callback = shift;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $frip = $peak_calling->fetch_Frip;
    
    if ($frip) {
      $frip_callback->($frip);
      next PEAK_CALLING;
    }
    $frip_run_failed_callback->($frip);
  }
  return;
}
