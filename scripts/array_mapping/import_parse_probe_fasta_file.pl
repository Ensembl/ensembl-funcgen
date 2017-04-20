#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1

perl /nfs/users/nfs_m/mn1/various_scripts/probemapping/test_parse_probe_fasta_file.pl \
  --array_name AGILENT \
  --species homo_sapiens \
  --probe_directory /lustre/scratch109/ensembl/funcgen/array_mapping/ \
  --parsed_output /lustre/scratch109/ensembl/funcgen/array_mapping/parsed_AGILENT_probes.pl

=cut

my $array_name;
my $probe_file;
my $parsed_output;

GetOptions (
   'array_name=s'      => \$array_name,
   'probe_file=s'      => \$probe_file,
   'parsed_output=s'   => \$parsed_output,
);

if (! -e $probe_file) {
  die("Can't find probe file ${probe_file}!");
}

use Bio::EnsEMBL::Funcgen::Config::ArrayInstance;
my $probemapping_config = Bio::EnsEMBL::Funcgen::Config::ArrayInstance->new();
$probemapping_config->read_and_check_config($array_name);

use Bio::EnsEMBL::Funcgen::Parsers::ArrayFasta;
my $probemapping_parser = Bio::EnsEMBL::Funcgen::Parsers::ArrayFasta->new();

open my $outfile, '>', $parsed_output;

my $parser_parameters = {
  IIDREGEXP   => $probemapping_config->IIDREGEXP,
  IFIELDORDER => $probemapping_config->IFIELDORDER,
  probe_file  => $probe_file,
  call_back   => sub {
    my $probe_data = shift;
    
    assert_probe_data_ok($probe_data);
    
    if ($probe_data->{-sequence} eq '') {
      warn "Skipping probe, because it has no sequence:\n" . Dumper($probe_data);
      return;
    }
    $outfile->print(Dumper($probe_data));
  }
};

$probemapping_parser->validate_configuration($parser_parameters);
$probemapping_parser->run($parser_parameters);

$outfile->close;

sub assert_probe_data_ok {
  my $probe_data = shift;
  if (! defined $probe_data->{-sequence}) {
    die;
  }
}