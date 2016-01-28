package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::CreateDB;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run { 
    my $self = shift;
    my $tracking_dba_hash = $self->param('tracking_dba_hash');
    
    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    
    my $dbc_name = $tracking_dba_hash->{-dbname};
    delete $tracking_dba_hash->{-dbname};
    
    use Bio::EnsEMBL::DBSQL::DBConnection;
    my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
    
    if ($self->db_exists($dbc, $dbc_name)) {    
      $logger->info("$dbc_name already exists, so not doing anything.\n");    
    } else {
      $logger->info("Creating $dbc_name and loading empty schema.\n");
      $self->create_empty_tracking_db($dbc, $dbc_name);
    }
    return;
}

sub db_exists {
    my $self    = shift;
    my $dbc     = shift;
    my $db_name = shift;

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
      -DB_CONNECTION => $dbc
    );
    my $count =
      $helper->execute_single_result(
	-SQL => qq(select count(*) from information_schema.schemata where schema_name = '$db_name'),
    );
    if ($count == 1) {
      return 1;
    }
    return;
}

sub create_empty_tracking_db {
    my $self     = shift;
    my $dbc      = shift;
    my $dbc_name = shift;
    
    $dbc->dbname($dbc_name);
    
    use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw(
      create_db_url_from_dbc
      create_db_server_url_from_dbc
    );
    my $url_db_server = create_db_server_url_from_dbc($dbc);
    my $url_db        = create_db_url_from_dbc($dbc);
    
    my $db_setup_commands = [
      qq% db_cmd.pl -url $url_db_server -sql "create database $dbc_name" %,
      qq% db_cmd.pl -url $url_db < sql/efg.sql %,
      qq% db_cmd.pl -url $url_db < sql/array2organism.sql %,
      qq% db_cmd.pl -url $url_db < sql/probe_seq.sql %,
      qq% db_cmd.pl -url $url_db < sql/probe_alias.sql %,
      qq% db_cmd.pl -url $url_db -sql "INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1, 'species.production_name','homo_sapiens')" %,
    ];

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
    
    foreach my $current_cmd (@$db_setup_commands) {
    
	$logger->info($current_cmd . "\n", undef, 1);
	run_system_cmd($current_cmd);
    }
    return;
}

1;
