package Bio::EnsEMBL::Funcgen::Report::PhantomPeaks;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/report_phantom_peak_bar_chart.html';
  return $template;
}

sub _static_content {

  my $self = shift;
  
  return {
    title    => 'Phantom peak quality checks',
    distribution_y_axis_max => 350,
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
  my $all_peak_callings    = $peak_calling_adaptor->fetch_all;

  my $phantom_peak_datasets = $self->create_phantom_peak_datasets($all_peak_callings);
  
  my $bins = $self->_generate_bins;
  
  return {
    labels   => $bins,
    datasets => $phantom_peak_datasets,
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,

  };
}

sub create_phantom_peak_datasets {

  my $self = shift;
  
  my $all_peak_callings = shift;

  my @datasets;

  push 
    @datasets, 
    {
      title => 'All consortia combined',
      all => $self->create_phantom_peak_dataset(

            $all_peak_callings,
            
            {
              title   => 'Narrow and Broad data combined',
              html_id => 'overall-all', 
              colour  => 'window.chartColors.gray',
            },

            sub { return 1; }
        ),
      narrow => $self->create_phantom_peak_dataset(

            $all_peak_callings,
            
            {
              title   => 'All narrow peaks',
              html_id => 'overall-narrow', 
              colour  => 'window.chartColors.gray',
            },

            sub { return (! shift->fetch_Experiment->get_FeatureType->creates_broad_peaks) }
        ),
      broad => $self->create_phantom_peak_dataset(

            $all_peak_callings,
            
            {
              title   => 'All broad peaks',
              html_id => 'overall-broad', 
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
      all => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Blueprint',
            html_id => 'blueprint-all',
            colour  => 'window.chartColors.red',
          },

        sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      ),
      narrow => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Blueprint narrow peaks',
            html_id => 'blueprint-narrow',
            colour  => 'window.chartColors.red',
          },

        sub { 
          my $peak_calling = shift;

                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
        }
      ),
      broad => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Blueprint broad peaks',
            html_id => 'blueprint-broad',
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
      all => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'ENCODE',
            html_id => 'encode-all',
            colour  => 'window.chartColors.blue',
          },

          sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      ),
      narrow => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'ENCODE narrow peaks',
            html_id => 'encode-narrow',
            colour  => 'window.chartColors.blue',
          },

          sub { 
            my $peak_calling = shift;
            
                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'ENCODE' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
          }
      ),
      broad => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'ENCODE broad peaks',
            html_id => 'encode-broad',
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
      all => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Roadmap Epigenomics',
            html_id => 'roadmap_epigenomics',
            colour  => 'window.chartColors.green',
          },

          sub { shift->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      ),
      narrow => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Roadmap Epigenomics narrow peaks',
            html_id => 'roadmap_epigenomics-narrow',
            colour  => 'window.chartColors.green',
          },

          sub { 
            my $peak_calling = shift;
            
                  $peak_calling->fetch_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' 
            && (! $peak_calling->fetch_Experiment->get_FeatureType->creates_broad_peaks)
          }
      ),
      broad => $self->create_phantom_peak_dataset(

          $all_peak_callings,
          
          {
            title   => 'Roadmap Epigenomics broad peaks',
            html_id => 'roadmap_epigenomics-broad',
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

sub create_phantom_peak_dataset {

  my $self = shift;

  my $peak_callings                = shift;
  my $initial_dataset              = shift;
  my $peak_calling_filter_callback = shift;

  my $filtered_peak_callings = [ grep { $peak_calling_filter_callback->($_) } @$peak_callings ];
  
  my $bins = $self->_generate_bins;
  
  my $final_dataset = {
      %$initial_dataset,
      counts             => $self->count_nsc_values_in_bins($filtered_peak_callings, $bins),
      quality_tag_values => $self->fetch_quality_tag_values($filtered_peak_callings),
  };
  
  return $final_dataset;
}

sub _generate_bins {

  my $self = shift;

  my $min = 1;
  my $max = 2;
  my $num_bins = 50;

  my @bins
    =
      map { $_ / $num_bins * ($max - $min) + $min }
      0..$num_bins
    ;
  return \@bins
}

sub count_nsc_values_in_bins {

  my $self = shift;

  my $peak_callings = shift;
  my $bins          = shift;

  my $values = $self->fetch_nsc_values($peak_callings);
  my $counts = $self->compute_bin_counts($bins, $values);
  
  return $counts;
}

sub compute_bin_counts {

  my $self = shift;

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

sub fetch_nsc_values {

  my $self = shift;
  
  my $peak_callings = shift;
  
  my @nsc_values;
  $self->process_phantom_peak_qc_values_from_peak_callings(
  
    $peak_callings,
    
    sub {
      my $phantom_peak = shift;
      my $nsc = $phantom_peak->nsc;
      
      push @nsc_values, $nsc;
    },
    
    sub {}
  );
  return \@nsc_values;
}

sub fetch_phantom_peak_success_rate_stats {

  my $self = shift;

  my $peak_callings = shift;
  
  my $run_successful = 0;
  my $run_failed     = 0;
  
  $self->process_phantom_peak_qc_values_from_peak_callings(
  
    $peak_callings,
    
    sub { $run_successful++ },
    sub { $run_failed++     }
  );
  
  my $phantom_peak_success_rate_stats = {
  
    run_successful => $run_successful,
    run_failed     => $run_failed,
  
  };
  
  return $phantom_peak_success_rate_stats;
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
