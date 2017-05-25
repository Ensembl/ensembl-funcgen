#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use DBI qw(:sql_types);

=head1

delete_redundant_probe_features.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species  homo_sapiens \
  --only_test 1

=cut

my $registry;
my $species;
my $only_test;

GetOptions (
   'registry=s'  => \$registry,
   'species=s'   => \$species,
   'only_test=s' => \$only_test,
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
my $sql = '
  select
    group_concat(cast(probe_feature_id as char)) ids, 
    count(probe_feature_id) num_occurrences, 
    seq_region_id, 
    seq_region_start, 
    seq_region_end, 
    probe_id, 
    mismatches, 
    cigar_line,
    hit_id
  from
    probe_feature
  group by 
    seq_region_id, 
    seq_region_start, 
    seq_region_end, 
    probe_id, 
    mismatches, 
    cigar_line,
    hit_id
  having 
    num_occurrences > 1 
  order by 
    num_occurrences desc;
';

# Prevent memory issues from buffering
$funcgen_db_adaptor->dbc->db_handle->{mysql_use_result} = 1;

my $helper = $funcgen_db_adaptor->dbc->sql_helper;
  
$logger->info("Running $sql\n");

my $last_line_printed;

my @probe_feature_ids_to_delete;

$helper->execute_no_return(
  -SQL      => $sql,
  -CALLBACK => sub {
    my $row  = shift;
    
    my @redundant_probe_feature_ids = split ',', $row->[0];
    (
      my $id_of_probe_feature_to_keep,
      my @ids_of_probe_features_to_remove
    ) = @redundant_probe_feature_ids;
    push @probe_feature_ids_to_delete, @ids_of_probe_features_to_remove;
  },
);

my $number_of_redundant_probe_features = scalar @probe_feature_ids_to_delete;
$logger->info("Found " . $number_of_redundant_probe_features . " redundant probe features.\n");

if ($number_of_redundant_probe_features==0) {
  exit 0;
}
if ($only_test && $number_of_redundant_probe_features>0) {
  exit 1;
}

my $sql = 'delete from probe_feature where probe_feature_id = ?;';
my $delete_by_probe_feature_id = $funcgen_db_adaptor->dbc->prepare($sql);

my $progressbar_id = $logger->init_progress(scalar @probe_feature_ids_to_delete);
$logger->info("Deleting redundant probe features\n");

my $number_of_features_deleted = 0;

foreach my $id_of_probe_features_to_remove (@probe_feature_ids_to_delete) {

  $delete_by_probe_feature_id->bind_param(1, $id_of_probe_features_to_remove, SQL_INTEGER);
  $delete_by_probe_feature_id->execute;
  $number_of_features_deleted++;
  $logger->log_progressbar($progressbar_id, $number_of_features_deleted);

}
$logger->info("Done.\n");

$logger->finish_log;
