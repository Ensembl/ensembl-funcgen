package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertExternalDb;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

sub run {   
    my $self = shift;
    
    my $tracking_dba_hash = $self->param('tracking_dba_hash');
    
    my $mapping_type = [
      'transcript',
      'genomic'
    ];

    my $species = 'homo_sapiens';

    my $externaldb = {
      'transcript' => {
	    db_name      => $species.'_core_Transcript',
	    display_name => 'EnsemblTranscript',
      },
      'genomic' => {
	    db_name      => $species.'_core_Genome',
	    display_name => 'EnsemblGenome',
      },
    };
    
    use Bio::EnsEMBL::DBSQL::DBConnection;
    my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
    
    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
      -DB_CONNECTION => $dbc
    );
    
    use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw( create_funcgen_adaptor create_dna_db_adaptor );
    my $funcgen_adaptor = create_funcgen_adaptor($tracking_dba_hash);    
    my $schema_build = $funcgen_adaptor->_get_schema_build($funcgen_adaptor->dnadb);
    
    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    
    foreach my $current_mapping_type (@$mapping_type) {
    
      my $db_name      = $externaldb->{$current_mapping_type}->{db_name};
      my $display_name = $externaldb->{$current_mapping_type}->{display_name};
      
      my $count = 
	$helper->execute_single_result(
	  -SQL    => "select count(external_db_id) from external_db where db_name=? and db_release=?",
	  -PARAMS => [ $db_name, $schema_build ],
      );

      if($count == 0) {      
	    $logger->info("Inserting external db entry for $db_name.\n");

	    my $insert_sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type)'.
	      " values(?, ?, 'KNOWNXREF', 1, 5, ?, 'MISC')";
	      use DBI qw(:sql_types);
	    $helper->batch(
	      -SQL => $insert_sql,	      
	      -data => [	      
		[
		  [ $db_name,      SQL_VARCHAR ],
		  [ $schema_build, SQL_VARCHAR ],
		  [ $display_name, SQL_VARCHAR ]
		]
	      ] 
	    );
      } else {
	  $logger->info("External db entry for $db_name found, no need to insert one.\n");
      }
    }
}
1;
