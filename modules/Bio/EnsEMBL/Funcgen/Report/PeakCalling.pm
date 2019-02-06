package Bio::EnsEMBL::Funcgen::Report::PeakCalling;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub _constructor_parameters {
  my $self = shift;
  return {
    %{$self->SUPER::_constructor_parameters},
    peak_calling => 'peak_calling',
  };
}

sub peak_calling { return shift->_generic_get_or_set('peak_calling',     @_); }

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/description.txt';
  return $template;
}

sub template_dir {
  my $self = shift;
  my $template_dir = $self->template_base_dir . '/peak_calling_description';
  return $template_dir;
}

sub _dynamic_content {

  my $self = shift;

  return {
    peak_calling  => $self->peak_calling,
  };
}

sub _in_template_functions {
  my $self = shift;
  return {
  
    %{$self->SUPER::_in_template_functions},
    
    canonpath => sub {
      my $path = shift;
      return File::Spec->canonpath($path)
    },
    
    bool_to_yes_no => sub {
      my $boolean = shift;
      if ($boolean) {
        return 'yes'
      }
      return 'no'
    },
    
    round_percent => sub {
      my $number = shift;
      return sprintf("%.2f", $number) . '%';
    },
    
    default_round => sub {
      my $number = shift;
      return sprintf("%.2f", $number);
    },
    
    scientific_notation => sub {
      my $number = shift;
      return sprintf("%.2e", $number);
    },
    
    fetch_deduplicated_partial_alignments => sub {
      my $alignment = shift;
      return fetch_deduplicated_partial_alignments($alignment);
    },
    
    summarise_read_file_experimental_configurations_in_alignment => sub {
      my $alignment = shift;
      return summarise_read_file_experimental_configurations_in_alignment($alignment);
    },
    
  }
}

sub fetch_deduplicated_partial_alignments {
  my $alignment = shift;
  
  my $read_file_experimental_configurations = [ 
    sort by_read_file_experimental_configuration @{$alignment->fetch_all_ReadFileExperimentalConfigurations}
  ];
  my @alignments;
  
  RFEC:
  foreach my $read_file_experimental_configuration (@$read_file_experimental_configurations) {

    my $alignments = $alignment
      ->adaptor
      ->fetch_all_by_ReadFileExperimentalConfiguration(
        $read_file_experimental_configuration
      );
    
    my $partial_alignments = 
      [
        grep {
             ! $_->is_complete
          && ! $_->has_duplicates
        } @$alignments
      ];
    if (@$partial_alignments == 0) {
      next RFEC;
    }
    if (@$partial_alignments != 1) {
      die(
        "Unexpected number of alignments returned!\n"
        . Dumper($partial_alignments)
      )
    }
    push @alignments, $partial_alignments->[0];
  }
  return \@alignments;
}

sub by_read_file_experimental_configuration {
  return $a->biological_replicate <=> $b->biological_replicate
      || $a->technical_replicate  <=> $b->technical_replicate;
}

sub summarise_read_file_experimental_configurations_in_alignment {
  my $alignment = shift;
  return summarise_read_file_experimental_configurations($alignment->fetch_all_ReadFileExperimentalConfigurations);
}

sub summarise_read_file_experimental_configurations {
  my $read_file_experimental_configurations = shift;
  
    my $summary = join ' ',
      map {
        summarise_read_file_experimental_configuration($_)
      } sort
        by_read_file_experimental_configuration
        @$read_file_experimental_configurations;
  return $summary;
}

sub summarise_read_file_experimental_configuration {
  my $read_file_experimental_configuration = shift;
  
  return 
    '(BR' . $read_file_experimental_configuration->biological_replicate 
    . ', ' 
    . 'TR' . $read_file_experimental_configuration->technical_replicate
    . ')'
}

1;
