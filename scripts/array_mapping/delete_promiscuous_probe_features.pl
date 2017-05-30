#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use DBI qw(:sql_types);

=head1

  delete_promiscuous_probe_features.pl \
    --species  mus_musculus \
    --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.testruns.pm \
    --analysis_logic_name ProbeAlign_genomic \
    --max_allowed_hits_per_probe 100

=cut

my $registry;
my $species;
my $only_test;
my $analysis_logic_name;
my $max_allowed_hits_per_probe = 100;

GetOptions (
   'registry=s'                   => \$registry,
   'species=s'                    => \$species,
   'analysis_logic_name=s'        => \$analysis_logic_name,
   'max_allowed_hits_per_probe=s' => \$max_allowed_hits_per_probe,
   'only_test=s'                  => \$only_test,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

# Regarding inclusion of hit_id: 
#
# We can't remove probe features that were generated from different sequences
# , even if all else is the same. The reason is that probe2transcript queries 
# whether there were any matches of a particular probe feature to a specific 
# transcript. And this is ultimately used to assign probes to transcripts.
#
my $sql = "
  select 
      probe_id, 
      count(probe_feature_id) total_probe_features,
      sum(
        case
            when(source='transcript') then 1
            else 0
        end
      ) as num_transcript,
      sum(
        case
            when(source='genomic') then 1
            else 0
        end
      ) as num_genomic
  from 
      probe_feature 
  group by 
      probe_id 
  having 
      total_probe_features > $max_allowed_hits_per_probe
  order by 
      total_probe_features desc 
";

# Prevent memory issues from buffering (not so important here)
$funcgen_db_adaptor->dbc->db_handle->{mysql_use_result} = 1;

my $helper = $funcgen_db_adaptor->dbc->sql_helper;
  
$logger->info("Running $sql\n");

my @promiscuous_probe_data;

$helper->execute_no_return(
  -SQL      => $sql,
  -CALLBACK => sub {
    my $row  = shift;
    my $promiscuous_probe_data_description = {
      probe_id                            => $row->[0],
      num_probe_features                  => $row->[1],
      num_probe_features_from_transcripts => $row->[2],
      num_probe_features_from_genome      => $row->[3],
    };
    push @promiscuous_probe_data, $promiscuous_probe_data_description;
  },
);

my $number_of_promiscuous_probes_found = scalar @promiscuous_probe_data;

$logger->info("Found " . $number_of_promiscuous_probes_found . " promiscuous probes.\n");

my $delete_sql = 'delete from probe_feature where probe_id = ?;';
my $delete_by_probe_id = $funcgen_db_adaptor->dbc->prepare($delete_sql);

my $progressbar_id = $logger->init_progress($number_of_promiscuous_probes_found);
$logger->info("Deleting redundant probe features\n");

my $number_of_features_deleted = 0;
my $number_of_promiscuous_probes_processed = 0;

my $unmapped_object_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'unmappedobject');
my $analysis_adaptor        = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');

use Bio::EnsEMBL::Analysis;
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $analysis_logic_name);

foreach my $current_promiscuous_probe_data (@promiscuous_probe_data) {

  my $probe_id                            = $current_promiscuous_probe_data->{probe_id};
  my $num_probe_features                  = $current_promiscuous_probe_data->{num_probe_features};
  my $num_probe_features_from_transcripts = $current_promiscuous_probe_data->{num_probe_features_from_transcripts};
  my $num_probe_features_from_genome      = $current_promiscuous_probe_data->{num_probe_features_from_genome};

  use Bio::EnsEMBL::UnmappedObject;
  my $unmapped_object = Bio::EnsEMBL::UnmappedObject->new (
    -type                => 'array_mapping',
    -analysis            => $analysis,
    -ensembl_id          => $probe_id,
    -ensembl_object_type => 'Probe',
    -external_db_id      => undef,
    -identifier          => $probe_id,
    -summary             => 'Promiscuous probe',
    -full_desc           => "Probe exceeded maximum allowed number of mappings when combining genomic ($num_probe_features_from_genome matches) and transcript matches ($num_probe_features_from_transcripts matches). (Total: $num_probe_features matches)"
  );
  $unmapped_object_adaptor->store($unmapped_object);

  $delete_by_probe_id->bind_param(1, $probe_id, SQL_INTEGER);
  $number_of_features_deleted += $delete_by_probe_id->execute;

  $number_of_promiscuous_probes_processed++;

  $logger->log_progressbar($progressbar_id, $number_of_promiscuous_probes_processed);

}
$logger->info("\nDone.\n");
$logger->info("$number_of_features_deleted probe features from $number_of_promiscuous_probes_processed promiscuous probes were deleted.\n");
$logger->finish_log;
