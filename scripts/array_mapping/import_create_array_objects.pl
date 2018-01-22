#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Carp;

=head1

perl /nfs/users/nfs_m/mn1/various_scripts/probemapping/test_create_array_objects.pl \
  --array_name AGILENT \
  --parsed_probe_data /lustre/scratch109/ensembl/funcgen/array_mapping/parsed_AGILENT_probes.pl
  --output_file AGILENT

=cut

my $array_name;
my $parsed_probe_data;
my $output_file;

GetOptions (
   'array_name=s'        => \$array_name,
   'parsed_probe_data=s' => \$parsed_probe_data,
   'output_file=s'       => \$output_file,
);

use Bio::EnsEMBL::Funcgen::Config::ArrayInstance;
my $probemapping_config = Bio::EnsEMBL::Funcgen::Config::ArrayInstance->new();
$probemapping_config->read_and_check_config($array_name);

open my $output_fh, '>', $output_file;

my $process_probe_data = sub {

  my $probe_data = shift;
  
  ### set description to null if not informed
  my $probe_desc = $probe_data->{'-description'};
  if (! defined $probe_desc || $probe_desc eq '') {
	  $probe_data->{'-description'} = undef;
  }
  
  my $array_name      = $probe_data->{'-array'};
  my $array_chip_name = $probe_data->{'-array_chip'};
  
#   if ($array_name eq '') {
#     confess(
#       "No array name for\n" 
#       . Dumper($probe_data)
#     );
#   }
  
  if (! defined $array_name || $array_name eq '') {
    confess("Array name was not defined for probe: " . Dumper($probe_data));
  }

  my $array_data = $probemapping_config->get_ARRAY_PARAMS_by_array_name($array_name);
  
#   print Dumper($array_name);
#   print Dumper($array_data);
  
  use Bio::EnsEMBL::Funcgen::Array;
  my $array = Bio::EnsEMBL::Funcgen::Array->new(%$array_data);
  
  use Bio::EnsEMBL::Funcgen::ArrayChip;
  my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
    -name      => $array->name,
    -design_id => $array->name,
  );
  if (exists $probe_data->{'-probe_set'}) {
  
    use Bio::EnsEMBL::Funcgen::ProbeSet;
    my $probe_set = Bio::EnsEMBL::Funcgen::ProbeSet->new(
      -name => $probe_data->{'-probe_set'}
    );
    $probe_data->{'-probe_set'} = $probe_set;
  }
  
  $probe_data->{'-array'}      = $array;
  $probe_data->{'-array_chip'} = $array_chip;
  
  use Bio::EnsEMBL::Funcgen::Probe;
  my $probe = Bio::EnsEMBL::Funcgen::Probe->new(%$probe_data);
  
#   print(Dumper($probe));
  $output_fh->print(Dumper($probe));
};

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$parser->parse({
  data_dumper_file => $parsed_probe_data,
  call_back        => $process_probe_data,
});

$output_fh->close;
