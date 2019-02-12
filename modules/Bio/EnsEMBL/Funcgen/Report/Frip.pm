package Bio::EnsEMBL::Funcgen::Report::Frip;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::Report::binnedValueReport';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/frip.html';
  return $template;
}

sub _static_content {
  return {
    title => 'Fraction of reads in peaks',
    distribution_y_axis_max => 250,
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
  my $all_peak_callings    = $peak_calling_adaptor->fetch_all;

  my $frip_datasets = $self->_compute_datasets_with_bin_counts($all_peak_callings);
  
  my $bins = $self->bins;
  
  return {
    labels   => $bins,
    datasets => $frip_datasets,
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,
  };
}

sub bins {
  my $self = shift;
  my $bins = $self->_generate_bins({
      min      =>  0,
      max      =>  1,
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
      title => 'All',
      all => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'All',
              colour  => 'window.chartColors.gray',
            },

            object_filter => sub { return 1; }
        }),
      narrow => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'All narrow peaks',
              colour  => 'window.chartColors.gray',
            },

            object_filter => sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) }
        }),
      broad => $self->_compute_dataset_with_bin_counts({

            object_list => $all_peak_callings,
            
            static_values => {
              title   => 'All broad peaks',
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
            colour  => 'window.chartColors.red',
          },

        object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Blueprint narrow peaks',
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
            colour  => 'window.chartColors.blue',
          },

          object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'ENCODE narrow peaks',
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
            colour  => 'window.chartColors.green',
          },

          object_filter => sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      }),
      narrow => $self->_compute_dataset_with_bin_counts({

          object_list => $all_peak_callings,
          
          static_values => {
            title   => 'Roadmap Epigenomics narrow peaks',
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

  my $self  = shift;
  my $peak_callings = shift;
  
  my @frip_values;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $frip = $peak_calling->fetch_Frip;
    
    if (! $frip) {
      next PEAK_CALLING;
    }
    
    my $frip_value = $frip->frip;
    push @frip_values, $frip_value;
  }
  
  my $bins = $self->bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@frip_values);

  return {
    counts => $bin_counts
  };
}

1;
