package Bio::EnsEMBL::Funcgen::Template::PeakCallingDescription;

use strict;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

our @EXPORT_OK = qw(
  PEAK_CALLING_TXT_TEMPLATE
  apply
);

sub apply {

  my $peak_calling = shift;

  use Template;
  my $tt = Template->new(
    ABSOLUTE     => 1,
    RELATIVE     => 1,
  );
  
  my $output;
  
  my $file = __FILE__;
  use File::Basename qw( dirname basename );

  my $template_dir = dirname($file) . '/../../../../../templates/peak_calling_description';
  my $description_template = $template_dir . '/description.txt';
  
  use Number::Format qw( format_number );
  
  $tt->process(
    $description_template, 
    {
      peak_calling  => $peak_calling,
      
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
      format_number => sub {
        my $number = shift;
        if (! defined $number) {
          return '-'
        }
        if ($number eq '') {
          return '-'
        }
        #return 'foo';
        return format_number($number);
      },
      fetch_deduplicated_partial_alignments => sub {
        my $alignment = shift;
        return fetch_deduplicated_partial_alignments($alignment);
      },
      summarise_read_file_experimental_configurations_in_alignment => sub {
        my $alignment = shift;
        return summarise_read_file_experimental_configurations_in_alignment($alignment);
      },
    },
    \$output
  )
      || die $tt->error;
  return $output;
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
      ->db
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
