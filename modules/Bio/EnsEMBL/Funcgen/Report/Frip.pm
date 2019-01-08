package Bio::EnsEMBL::Funcgen::Report::Frip;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/frip.html';
  return $template;
}

sub _static_content {

  my $self = shift;
  
  return {
    title    => 'Fraction of reads in peaks',
    distribution_y_axis_max => 250,
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
  my $all_peak_callings    = $peak_calling_adaptor->fetch_all;

  my $frip_datasets = $self->create_frip_datasets($all_peak_callings);
  
  my $bins = $self->_generate_bins;
  
  return {
    labels   => $bins,
    datasets => $frip_datasets,
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,

  };
}

sub _generate_bins {

  my $self = shift;

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

sub create_frip_datasets {

  my $self = shift;
  my $all_peak_callings = shift;

  my @datasets;

  push 
    @datasets, 
    {
      title => 'All',
      all => $self->create_frip_dataset(

            $all_peak_callings,
            
            {
              title   => 'All',
              colour  => 'window.chartColors.gray',
            },

            sub { return 1; }
        ),
      narrow => $self->create_frip_dataset(

            $all_peak_callings,
            
            {
              title   => 'All narrow peaks',
              colour  => 'window.chartColors.gray',
            },

            sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) }
        ),
      broad => $self->create_frip_dataset(

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
      all => $self->create_frip_dataset(

          $all_peak_callings,
          
          {
            title   => 'Blueprint',
            colour  => 'window.chartColors.red',
          },

        sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      ),
      narrow => $self->create_frip_dataset(

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
      broad => $self->create_frip_dataset(

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
      all => $self->create_frip_dataset(

          $all_peak_callings,
          
          {
            title   => 'ENCODE',
            colour  => 'window.chartColors.blue',
          },

          sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      ),
      narrow => $self->create_frip_dataset(

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
      broad => $self->create_frip_dataset(

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
      all => $self->create_frip_dataset(

          $all_peak_callings,
          
          {
            title   => 'Roadmap Epigenomics',
            colour  => 'window.chartColors.green',
          },

          sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      ),
      narrow => $self->create_frip_dataset(

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
      broad => $self->create_frip_dataset(

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

  return \@datasets;
}

sub create_frip_dataset {
  
  my $self = shift;
  
  my $peak_callings                = shift;
  my $initial_dataset              = shift;
  my $peak_calling_filter_callback = shift;

  my $filtered_peak_callings = [ grep { $peak_calling_filter_callback->($_) } @$peak_callings ];
  
  my $bins = $self->_generate_bins;
  
  my $final_dataset = {
      %$initial_dataset,
      counts => $self->count_frip_values_in_bins($filtered_peak_callings, $bins),
  };
  
  return $final_dataset;
}

sub count_frip_values_in_bins {

  my $self          = shift;
  my $peak_callings = shift;
  my $bins          = shift;

  my $values = $self->fetch_frip_values($peak_callings);
  my $counts = $self->compute_bin_counts($bins, $values);
  
  return $counts;
}

sub compute_bin_counts {

  my $self   = shift;
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
  
  my $self = shift;
  my $peak_callings = shift;
  
  my @frip_values;
  $self->process_frip_qc_values_from_peak_callings(
  
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

  my $self = shift;
  
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

1;
