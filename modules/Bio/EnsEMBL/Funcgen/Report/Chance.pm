package Bio::EnsEMBL::Funcgen::Report::Chance;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::Report::binnedValueReport';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/chance.html';
  return $template;
}

sub _static_content {
  return {
    title => 'Chance',
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
  my $all_peak_callings    = $peak_calling_adaptor->fetch_all;

  my $differential_percentage_enrichment_datasets
    = $self->_differential_percentage_enrichment_datasets($all_peak_callings);

  my $datasets = $self->_divergence_datasets($all_peak_callings);
  
  my $percent_genome_enriched_datasets
    = $self->_percent_enrichment_datasets($all_peak_callings);
    
  return {
    
    datasets => $datasets,
    
    datasets_graph => {
      labels       => $self->divergence_bins,
      title        => 'Divergence',
      y_axis_label => 'Number of peak callings',
      y_axis_max   => 600,
      x_axis_label => '',
    },
    
    percent_genome_enriched_datasets => $percent_genome_enriched_datasets,
    
    percent_genome_enriched_graph => {
      labels       => $self->percent_genome_enriched_bins,
      title        => 'Percent genome enriched',
      y_axis_label => 'Number of peak callings',
      y_axis_max   => 120,
      x_axis_label => 'Percent enrichment',
    },

    differential_percentage_enrichment_datasets
      => $differential_percentage_enrichment_datasets,
    
    differential_percentage_enrichment_graph => {
      labels       => $self->percent_genome_enriched_bins,
      title        => 'Differential percentage enrichment',
      y_axis_label => 'Number of peak callings',
      y_axis_max   => 250,
      x_axis_label => 'Differential percentage enrichment',
    },

    dbc      => $funcgen_adaptor->dbc,
    species  => $species,
  };
}

sub divergence_bins {
  my $self = shift;
  my $bins = $self->_generate_bins({
      min      =>  0,
      max      =>  0.2,
      num_bins => 50,
  });
  return $bins
}

sub percent_genome_enriched_bins {
  my $self = shift;
  my $bins = $self->_generate_bins({
      min      =>   0,
      max      => 100,
      num_bins =>  50,
  });
  return $bins
}

sub differential_percentage_enrichment_bins {
  my $self = shift;
  return $self->percent_genome_enriched_bins;
}

sub _divergence_datasets {

  my $self = shift;
  my $all_peak_callings = shift;

  my $datasets = [
    {
      title => 'Chance',
      all => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'divergence',
              colour  => 'window.chartColors.red',
            },

            object_list => $all_peak_callings,

            compute_value_method => '_divergence_values'
        }),
    },
    {
      title => 'Chance',
      all => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'sqrt divergence',
              colour  => 'window.chartColors.green',
            },

            object_list => $all_peak_callings,

            compute_value_method => '_sqrt_divergence_values'
        }),
    }
  ];
  return $datasets;
}

sub _percent_enrichment_datasets {

  my $self = shift;
  my $all_peak_callings = shift;

  my $datasets = [
    {
      title => 'Percent genome enriched',
      all => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'All experiments',
              colour  => 'window.chartColors.gray',
            },

            object_list => $all_peak_callings,

            compute_value_method => '_percent_enrichment_values'
        }),
      narrow => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'Narrow peak experiments',
              colour  => 'window.chartColors.red',
            },

            object_list => $all_peak_callings,
            
            object_filter => sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) },

            compute_value_method => '_percent_enrichment_values'
        }),
      broad => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'Broad peak experiments',
              colour  => 'window.chartColors.green',
            },

            object_list => $all_peak_callings,
            
            object_filter => sub { return (shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) },

            compute_value_method => '_percent_enrichment_values'
        }),
    },
  ];
  return $datasets;
}

sub _differential_percentage_enrichment_datasets {

  my $self = shift;
  my $all_peak_callings = shift;

  my $datasets = [
    {
      title => 'Differential percentage enrichment',
      all => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'All experiments',
              colour  => 'window.chartColors.gray',
            },

            object_list => $all_peak_callings,

            compute_value_method => '_differential_percentage_enrichment_values'
        }),
      narrow => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'Narrow peak experiments',
              colour  => 'window.chartColors.red',
            },

            object_list => $all_peak_callings,
            
            object_filter => sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) },

            compute_value_method => '_differential_percentage_enrichment_values'
        }),
      broad => $self->_compute_dataset_with_bin_counts({

            static_values => {
              title   => 'Broad peak experiments',
              colour  => 'window.chartColors.green',
            },

            object_list => $all_peak_callings,
            
            object_filter => sub { return (shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) },

            compute_value_method => '_differential_percentage_enrichment_values'
        }),
    },
  ];
  return $datasets;
}

sub _differential_percentage_enrichment_values {

  my $self  = shift;
  my $peak_callings = shift;
  
  my @differential_percentage_enrichment;
  
  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $chance = $peak_calling->fetch_Chance;
    
    if (! $chance) {
      next PEAK_CALLING;
    }
    if ($chance->run_failed) {
      next PEAK_CALLING;
    }
    
    my $differential_percentage_enrichment = $chance->differential_percentage_enrichment;
    push @differential_percentage_enrichment, $differential_percentage_enrichment;
  }
  my $bins = $self->differential_percentage_enrichment_bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@differential_percentage_enrichment);

  return {
    counts => $bin_counts,
  };
}

sub _percent_enrichment_values {

  my $self  = shift;
  my $peak_callings = shift;
  
  my @chance_percent_genome_enriched;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $chance = $peak_calling->fetch_Chance;
    
    if (! $chance) {
      next PEAK_CALLING;
    }
    if ($chance->run_failed) {
      next PEAK_CALLING;
    }
    
    my $percent_genome_enriched = $chance->percent_genome_enriched;
    push @chance_percent_genome_enriched, $percent_genome_enriched;
  }
  
  my $bins = $self->percent_genome_enriched_bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@chance_percent_genome_enriched);

  return {
    counts => $bin_counts,
  };
}

sub _divergence_values {

  my $self  = shift;
  my $peak_callings = shift;
  
  my @chance_divergence;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $chance = $peak_calling->fetch_Chance;
    
    if (! $chance) {
      next PEAK_CALLING;
    }
    if ($chance->run_failed) {
      next PEAK_CALLING;
    }
    
    my $divergence = $chance->divergence;
    push @chance_divergence, $divergence;
  }
  
  my $bins = $self->divergence_bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@chance_divergence);

  return {
    counts => $bin_counts,
  };
}

sub _sqrt_divergence_values {

  my $self  = shift;
  my $peak_callings = shift;
  
  my @chance_divergence;
  my @chance_sqrt_divergence;

  PEAK_CALLING:
  foreach my $peak_calling (@$peak_callings) {

    my $chance = $peak_calling->fetch_Chance;
    
    if (! $chance) {
      next PEAK_CALLING;
    }
    if ($chance->run_failed) {
      next PEAK_CALLING;
    }
    
    my $divergence = $chance->divergence;
    push @chance_divergence, $divergence;
    
    use Math::Complex;
    my $sqrt_divergence = sqrt($divergence);
    push @chance_sqrt_divergence, $sqrt_divergence;
  }
  
  my $bins = $self->divergence_bins;

  my $bin_counts
    = $self->_count_values_per_bin($bins, \@chance_divergence);

  my $bin_counts_sqrt
    = $self->_count_values_per_bin($bins, \@chance_sqrt_divergence);
    
  return {
    counts => $bin_counts_sqrt,
  };
}

1;
