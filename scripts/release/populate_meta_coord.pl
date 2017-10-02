use strict;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

=head1 

    perl scripts/release/populate_meta_coord.pl --registry registry.pl --species homo_sapiens

=cut 

use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Hash::Util qw ( lock_hash );

my $registry;
my $species;
my $gff_file;
my $bed_file;
my $min_id;
my $max_id;
my $feature_set_id;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $core_dbc = $core_dba->dbc;

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $funcgen_dbc = $funcgen_dba->dbc;

my @feature_table = qw(
  peak
  external_feature
  mirna_target_feature
  motif_feature
  probe_feature
  regulatory_feature
);

# Global variable to save redundant lookups
#
my $seq_region_id_to_coord_system_id_lookup = fetch_core_seq_region_id_to_coord_system_id_lookup({
  core_dba => $core_dba 
});


my $sql = "truncate meta_coord;";
my $sth = $funcgen_dbc->prepare($sql);
$sth->execute;

$sql = "INSERT INTO meta_coord (table_name,coord_system_id,max_length) VALUES (?, ?, ?);";
my $populate_meta_coord_sth = $funcgen_dbc->prepare($sql);

foreach my $current_feature_table (@feature_table) {

  $logger->info("Processing: $current_feature_table\n");

  my $coord_system_id_to_longest_feature_on_it = create_hash_of_coord_system_id_to_longest_feature_on_it({
    core_dba      => $core_dba,
    funcgen_dba   => $funcgen_dba,
    feature_table => $current_feature_table,
  });
  
  foreach my $current_coord_system_id (keys %$coord_system_id_to_longest_feature_on_it) {
  
    my $length_of_longest_feature = $coord_system_id_to_longest_feature_on_it->{$current_coord_system_id};
    
    $logger->info("    Writing " . $current_feature_table . " " . $current_coord_system_id . " " . $length_of_longest_feature . "\n");
    
    $populate_meta_coord_sth->bind_param(1, $current_feature_table);
    $populate_meta_coord_sth->bind_param(2, $current_coord_system_id);
    $populate_meta_coord_sth->bind_param(3, $length_of_longest_feature);
    $populate_meta_coord_sth->execute;
  }
}

$logger->info("Done.\n");

=head2 fetch_length_of_longest_feature_for_each_seq_region
=cut
sub create_hash_of_coord_system_id_to_longest_feature_on_it {
  my $param = shift;

  my $core_dba      = $param->{core_dba};
  my $funcgen_dba   = $param->{funcgen_dba};
  my $feature_table = $param->{feature_table};
  
  my $length_of_longest_feature_for_each_seq_region = fetch_length_of_longest_feature_for_each_seq_region({
    funcgen_dba   => $funcgen_dba,
    feature_table => $feature_table,
  });

  my $seq_region_ids_with_features_on_it = fetch_all_seq_region_ids_with_features_on_it({
    funcgen_dba   => $funcgen_dba,
    feature_table => $feature_table,
  });

  my %seen_coord_system_ids;
  my %coord_system_id_to_longest_feature_on_it;

  foreach my $seq_region_id_with_features_on_it (@$seq_region_ids_with_features_on_it) {

    my $coord_system_id           = $seq_region_id_to_coord_system_id_lookup->{$seq_region_id_with_features_on_it};
    my $length_of_longest_feature = $length_of_longest_feature_for_each_seq_region->{$seq_region_id_with_features_on_it};
    
    $seen_coord_system_ids{$coord_system_id} = 1;
    
    if (! exists $coord_system_id_to_longest_feature_on_it{$coord_system_id}) {
      $coord_system_id_to_longest_feature_on_it{$coord_system_id} = $length_of_longest_feature;
    }
    if ($length_of_longest_feature > $coord_system_id_to_longest_feature_on_it{$coord_system_id}) {
      $coord_system_id_to_longest_feature_on_it{$coord_system_id} = $length_of_longest_feature;
    }
  }
  return \%coord_system_id_to_longest_feature_on_it;
}

=head2 fetch_length_of_longest_feature_for_each_seq_region
=cut
sub fetch_length_of_longest_feature_for_each_seq_region {

  my $param = shift;

  my $funcgen_dba   = $param->{funcgen_dba};
  my $feature_table = $param->{feature_table};

  my $sth = $funcgen_dbc->prepare("select seq_region_id, max(seq_region_end - seq_region_start + 1) as length_of_longest_feature from " . $feature_table . " group by seq_region_id");

  $sth->execute;
  
  my %result;
  while (my $current_hash = $sth->fetchrow_hashref) {
    
    my $seq_region_id = $current_hash->{seq_region_id};
    my $length_of_longest_feature = $current_hash->{length_of_longest_feature};
    
    $result{$seq_region_id} = $length_of_longest_feature;
  }
  return \%result;
}

=head2 fetch_all_seq_region_ids_with_features_on_it
=cut
sub fetch_all_seq_region_ids_with_features_on_it {

  my $param = shift;

  my $funcgen_dba   = $param->{funcgen_dba};
  my $feature_table = $param->{feature_table};

  my $sth = $funcgen_dbc->prepare("select distinct seq_region_id from " . $feature_table);

  $sth->execute;

  #   $VAR1 = [
  #             [
  #               '131196'
  #             ],
  #             [
  #               '131211'
  #             ],
  #   ...
  #   }
  #
  my $arrayref_raw = $sth->fetchall_arrayref;
  my @seq_region_ids = map { $_->[0] } @$arrayref_raw;

  return \@seq_region_ids;
}

=head2 fetch_core_seq_region_id_to_coord_system_id_lookup
=cut
sub fetch_core_seq_region_id_to_coord_system_id_lookup {

  my $param = shift;
  my $core_dba = $param->{core_dba};
  
  my $core_dbc = $core_dba->dbc;
  my $sth = $core_dbc->prepare('select distinct seq_region_id, coord_system_id from seq_region');
  $sth->execute;

  #   $VAR1 = {
  #             '4' => {
  #                     'coord_system_id' => '1',
  #                     'seq_region_id' => '4'
  #                   },
  #             '1' => {
  #                     'coord_system_id' => '1',
  #                     'seq_region_id' => '1'
  #                   },
  #             '3' => {
  #                     'coord_system_id' => '1',
  #                     'seq_region_id' => '3'
  #                   },
  #             '2' => {
  #                     'coord_system_id' => '1',
  #                     'seq_region_id' => '2'
  #                   },
  #             '5' => {
  #                     'coord_system_id' => '1',
  #                     'seq_region_id' => '5'
  #                   }
  #           ...
  #           };
  #
  my $hashref_raw = $sth->fetchall_hashref('seq_region_id');
  
  my %seq_region_id_to_coord_system_id_lookup;
  foreach my $seq_region_id (keys %$hashref_raw) {
  
    my $coord_system_id = $hashref_raw->{$seq_region_id}->{coord_system_id};
    $seq_region_id_to_coord_system_id_lookup{$seq_region_id} = $coord_system_id;
  
  }
  lock_hash(%seq_region_id_to_coord_system_id_lookup);
  return \%seq_region_id_to_coord_system_id_lookup;
}

