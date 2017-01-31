#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Analysis;

=head1 create_probe_mapping_analyses.pl

create_probe_mapping_analyses.pl \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species mus_musculus
=cut

my $registry;
my $species;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
);

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $analysis_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');

my $analysis;

$analysis = Bio::EnsEMBL::Analysis->new(
  -logic_name      => 'ProbeAlign_transcript',
  -db              => undef,
  -db_version      => undef,
  -db_file         => undef,
  -program         => undef,
  -program_version => undef,
  -program_file    => undef,
  -gff_source      => undef,
  -gff_feature     => undef,
  -module          => undef,
  -module_version  => undef,
  -parameters      => undef,
  -created         => undef,
  -description     => 'Probe alignment',
  -display_label   => 'Probe alignment',
  -displayable     => 1,
  -web_data        => {'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }
);

$analysis_adaptor->store($analysis);

$analysis = Bio::EnsEMBL::Analysis->new(
  -logic_name      => 'ProbeAlign_genomic',
  -db              => undef,
  -db_version      => undef,
  -db_file         => undef,
  -program         => undef,
  -program_version => undef,
  -program_file    => undef,
  -gff_source      => undef,
  -gff_feature     => undef,
  -module          => undef,
  -module_version  => undef,
  -parameters      => undef,
  -created         => undef,
  -description     => 'Probe alignment',
  -display_label   => 'Probe alignment',
  -displayable     => 1,
  -web_data        => {'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }
);

$analysis_adaptor->store($analysis);

$analysis = Bio::EnsEMBL::Analysis->new(
  -logic_name      => 'probe2transcript',
  -db              => undef,
  -db_version      => undef,
  -db_file         => undef,
  -program         => undef,
  -program_version => undef,
  -program_file    => undef,
  -gff_source      => undef,
  -gff_feature     => undef,
  -module          => undef,
  -module_version  => undef,
  -parameters      => undef,
  -created         => undef,
  -description     => 'Microarray probes from manufacturers are aligned to the genome by Ensembl, if the probe sequences are provided. The mapping is a two-step procedure outlined <a href="/info/genome/microarray_probe_set_mapping.html">here</a>.',
  -display_label   => 'Probe2Transcript Annotation',
  -displayable     => undef,
  -web_data        => undef
);

$analysis_adaptor->store($analysis);

