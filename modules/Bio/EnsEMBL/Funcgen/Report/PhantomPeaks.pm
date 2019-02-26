package Bio::EnsEMBL::Funcgen::Report::PhantomPeaks;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::Report::binnedValueReport';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/report_phantom_peak_bar_chart.html';
  return $template;
}

sub _static_content {
  return {
    title => 'Phantom peak quality checks',
    distribution_y_axis_max => 350,
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
  my $all_peak_callings    = $peak_calling_adaptor->fetch_all;

  my $phantom_peak_datasets = $self->_compute_datasets_with_bin_counts($all_peak_callings);
  
  my $bins = $self->bins;
  
  return {
    labels   => $bins,
    datasets => $phantom_peak_datasets,
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,
  };
}

sub bins {
  my $self = shift;
  my $bins = $self->_generate_bins({
      min      =>  1,
      max      =>  2,
      num_bins => 50,
  });
  return $bins
}

sub _compute_datasets_with_bin_counts {

  my $self = shift;
  
  my $all_peak_callings = shift;

  my @datasets;

  push 
    @datasets, 
    {
      title => 'All consortia combined',
      all => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'Narrow and Broad data combined',
              html_id => 'overall-all', 
              colour  => 'window.chartColors.gray',
            },

            object_filter => sub { return 1; }
        }),
      narrow => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'All narrow peaks',
              html_id => 'overall-narrow', 
              colour  => 'window.chartColors.gray',
            },

            object_filter => sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) }
        }),
      broad => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'All broad peaks',
              html_id => 'overall-broad', 
              colour  => 'window.chartColors.gray',
            },

            object_filter => sub { return shift->fetch_Experiment->get_FeatureType->creates_broad_peaks }
        }),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Blueprint',
      all => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Blueprint',
            html_id => 'blueprint-all',
            colour  => 'window.chartColors.red',
          },

        object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Blueprint narrow peaks',
            html_id => 'blueprint-narrow',
            colour  => 'window.chartColors.red',
          },

        object_filter => sub { 
          my $peak_calling = shift;

                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
        }
      }),
      broad => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Blueprint broad peaks',
            html_id => 'blueprint-broad',
            colour  => 'window.chartColors.red',
          },

        object_filter => sub { 
          my $peak_calling = shift;

              $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
            && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
        }
      }),
    }
  ;

  push 
    @datasets, 
    {
      title => 'ENCODE',
      all => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'ENCODE',
            html_id => 'encode-all',
            colour  => 'window.chartColors.blue',
          },

          object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'ENCODE narrow peaks',
            html_id => 'encode-narrow',
            colour  => 'window.chartColors.blue',
          },

          object_filter => sub { 
            my $peak_calling = shift;
            
                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
          }
      }),
      broad => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'ENCODE broad peaks',
            html_id => 'encode-broad',
            colour  => 'window.chartColors.blue',
          },

          object_filter => sub { 
            my $peak_calling = shift;
            
              $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
            && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
          }
      }),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Roadmap Epigenomics',
      all => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Roadmap Epigenomics',
            html_id => 'roadmap_epigenomics',
            colour  => 'window.chartColors.green',
          },

          object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Roadmap Epigenomics narrow peaks',
            html_id => 'roadmap_epigenomics-narrow',
            colour  => 'window.chartColors.green',
          },

          object_filter => sub { 
            my $peak_calling = shift;
            
                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
          }
      }),
      broad => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Roadmap Epigenomics broad peaks',
            html_id => 'roadmap_epigenomics-broad',
            colour  => 'window.chartColors.green',
          },

          object_filter => sub { 
            my $peak_calling = shift;
            
              $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
            && $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks
          }
      }),
    }
  ;
  return \@datasets;
}

sub _compute_values_from_object_list {

  my $self          = shift;
  my $peak_callings = shift;

  my $quality_tag_values
    = $self->fetch_quality_tag_values($peak_callings);

  my @nsc_values;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $signal_alignment = $peak_calling->fetch_signal_Alignment;
    my $phantom_peak     = $signal_alignment->fetch_PhantomPeak;
    
    my $phantom_peak_run_failed = $phantom_peak->run_failed;
    
    if ($phantom_peak_run_failed) {
      next PEAK_CALLING;
    }
    
    my $nsc = $phantom_peak->nsc;
    push @nsc_values, $nsc;
  }

  my $bins = $self->bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@nsc_values);

  return {
    counts => $bin_counts,
    quality_tag_values => $quality_tag_values,
  };
}

sub fetch_quality_tag_values {

  my $self = shift;

  my $peak_callings = shift;
  
  my %quality_tags = (
    -2 => 0,
    -1 => 0,
     0 => 0,
     1 => 0,
     2 => 0,
     'run_failed' => 0,
  );
  
  $self->process_phantom_peak_qc_values_from_peak_callings(
  
    $peak_callings,
    
    sub {
      my $phantom_peak = shift;
      my $quality_tag = $phantom_peak->quality_tag;
      
      if (! exists $quality_tags{$quality_tag}) {
        $quality_tags{$quality_tag} = 1;
      } else {
        $quality_tags{$quality_tag}++;
      }
    },
    
    sub {
      $quality_tags{'run_failed'}++;
    }
  );
  
  my $quality_tag_list = [
    $quality_tags{'run_failed'},
    $quality_tags{-2},
    $quality_tags{-1},
    $quality_tags{0},
    $quality_tags{1},
    $quality_tags{2},
  ];
  
  return $quality_tag_list;
}

sub process_phantom_peak_qc_values_from_peak_callings {

  my $self = shift;

  my $peak_callings                    = shift;
  my $phantom_peak_callback            = shift;
  my $phantom_peak_run_failed_callback = shift;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $signal_alignment = $peak_calling->fetch_signal_Alignment;
    my $phantom_peak     = $signal_alignment->fetch_PhantomPeak;
    
    if (! defined $phantom_peak) {
      die('No phantom peak result for ' . $signal_alignment->name . " found!");
    }
    
    my $phantom_peak_run_failed = $phantom_peak->run_failed;
    
    if ($phantom_peak_run_failed) {
      $phantom_peak_run_failed_callback->($phantom_peak);
      next PEAK_CALLING;
    }
    $phantom_peak_callback->($phantom_peak);
  }
  return;
}

1;
