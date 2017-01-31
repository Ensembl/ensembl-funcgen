#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1 store_unmapped_objects.pl

load_probeset_to_transcript_rejections.pl \
  --probeset_rejections_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp_integration_test/mus_musculus/probe2transcript/rejected_probesets.pl \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --analysis_logic_name ProbeAlign_transcript \
  --species mus_musculus

select * from unmapped_object join unmapped_reason using (unmapped_reason_id) where ensembl_object_type="ProbeSet" and summary_description!="Insufficient hits" and summary_description != "No transcript mappings" and summary_description != "Promiscuous ProbeSet"

select * from unmapped_object join unmapped_reason using (unmapped_reason_id) where ensembl_object_type="ProbeSet" and summary_description="Insufficient hits"


=cut

my $probeset_rejections_file;
my $registry;
my $species;
my $analysis_logic_name;
my $target_type;

GetOptions (
   'probeset_rejections_file=s' => \$probeset_rejections_file,
   'registry=s'                 => \$registry,
   'species=s'                  => \$species,
   'analysis_logic_name=s'      => \$analysis_logic_name,
);

if (! -e $probeset_rejections_file) {
  die("Can't find file ${probeset_rejections_file}!");
}

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $analysis_adaptor        = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
my $unmapped_object_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'unmappedobject');

use Bio::EnsEMBL::Analysis;
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $analysis_logic_name);

my $process_data = sub {

  my $probeset_rejection = shift;
  
  use Bio::EnsEMBL::UnmappedObject;
  my $unmapped_object = Bio::EnsEMBL::UnmappedObject->new (
    -type                => $probeset_rejection->{type},
    -analysis            => $analysis,
    -ensembl_id          => $probeset_rejection->{dbID},
    -ensembl_object_type => $probeset_rejection->{object_type},
    -external_db_id      => undef,
    -identifier          => $probeset_rejection->{stable_id},
    -summary             => $probeset_rejection->{summary},
    -full_desc           => $probeset_rejection->{full_description},
  );
  $unmapped_object_adaptor->store($unmapped_object);
};

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$parser->parse({
  data_dumper_file => $probeset_rejections_file,
  call_back        => $process_data,
});
