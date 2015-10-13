package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertProbeAlias;

use strict;
use Carp;
use base ('Bio::EnsEMBL::Hive::Process');

sub run { 
    my $self = shift;
    
    my $tracking_dba_hash = $self->param('tracking_dba_hash');
  
    my $dbc_tracking = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
    $dbc_tracking->connect();

    use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw( create_probe_adaptor );
    my $probe_dba = create_probe_adaptor($tracking_dba_hash);
    
    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
      -DB_CONNECTION => $dbc_tracking 
    );

    my $num_probe_seq_ids_in_probe_table =
      $helper->execute_single_result(
	-SQL => 'select count(distinct probe_seq_id) from probe',
    );

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->info("There are $num_probe_seq_ids_in_probe_table probe seq ids in the probe table.\n");


    my $sth = $dbc_tracking->prepare('select min(probe_id) as default_probe_id, count(probe_id) as num_synonyms, probe_seq_id from probe where probe_seq_id is not null group by probe_seq_id order by num_synonyms desc, default_probe_id asc');
    $sth->execute;

    my $sth_fetch_aliases = $dbc_tracking->prepare('select probe_id, name from probe where probe_seq_id=? and probe_id != ?');
    my $sth_store_alias   = $dbc_tracking->prepare('INSERT ignore INTO probe_alias (probe_id,alias,default_name) VALUES (?, ?, ?);');

    my $progressbar_id = $logger->init_progress($num_probe_seq_ids_in_probe_table, 100);
    $logger->info("Updating probe_alias table\n");
    
    my $i = 0;
    while (my $data = $sth->fetchrow_hashref) {

      # This is the data for the default probe
      #
      my $num_synonyms     = $data->{num_synonyms};
      my $probe_seq_id     = $data->{probe_seq_id};
      my $default_probe_id = $data->{default_probe_id};

      my $default_probe       = $probe_dba->fetch_by_dbID($default_probe_id);  
      my $default_probe_names = $default_probe->get_all_probenames;
      
      # Fetch all other probes sharing the same dna sequence
      #
      $sth_fetch_aliases->bind_param(1, $probe_seq_id);
      $sth_fetch_aliases->bind_param(2, $default_probe_id);
      $sth_fetch_aliases->execute;

      my $default_probe_name = shift @$default_probe_names;
      
      # Collect all known aliases for the default probe
      #
      my @alias;  
      foreach my $current_probe_name (@$default_probe_names) {
	push @alias, {
	  probe_id   => $default_probe_id,
	  probe_name => $current_probe_name
	};
      }
      
      while (my $aliases = $sth_fetch_aliases->fetchrow_hashref) {
	push @alias, {
	  probe_id   => $aliases->{probe_id},
	  probe_name => $aliases->{name}
	};
      }
      
      # Store the default probe
      #
      $sth_store_alias->bind_param(1, $default_probe_id);
      $sth_store_alias->bind_param(2, $default_probe_name);
      $sth_store_alias->bind_param(3, 1);
      $sth_store_alias->execute;
      
      # Store the aliases
      #
      foreach my $current_alias (@alias) {
      
	$sth_store_alias->bind_param(1, $current_alias->{probe_id});
	$sth_store_alias->bind_param(2, $current_alias->{probe_name});
	$sth_store_alias->bind_param(3, 0);
	$sth_store_alias->execute;
      }
      $i++;
      $logger->log_progressbar($progressbar_id, $i);
    }
    $logger->info("Done updating the probe_alias table.");
}

1;

