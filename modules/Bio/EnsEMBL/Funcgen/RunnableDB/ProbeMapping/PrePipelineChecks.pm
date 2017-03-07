package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::PrePipelineChecks;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();    
    my $error_msg;

#     #
#     # pipeline_repository_dir is needed to find other files
#     # 
#     my $pipeline_repository_dir = $self->param('pipeline_repository_dir');    
#     if (! $pipeline_repository_dir) {    
#       $error_msg .= "pipeline_repository_dir parameter has not been set.\n";
#     } elsif (! -d $pipeline_repository_dir) {    
#       $error_msg .= "The directory $pipeline_repository_dir does not exist.\n";
#     }
#     
#     if ($error_msg) {
#       die($error_msg);
#     }    
#     #
#     # Files we are expecting but might not be there
#     #    
#     my $expected_files = [
#       'sql/table.sql',
# #       'sql/array2organism.sql',
#       'sql/probe_seq.sql',
# #       'sql/probe_alias.sql'
#     ];
#     
#     foreach my $current_expected_file (@$expected_files) {
# 	my $file_name_with_full_path = $pipeline_repository_dir . '/' . $current_expected_file;
#     
# 	if (! -e $file_name_with_full_path) {
# 	  $error_msg .= "File $current_expected_file in $pipeline_repository_dir is missing.\n";
# 	}
#     }

    my @programs_expected_in_path = qw(
      import_parse_probe_fasta_file.pl
      import_create_array_objects.pl
      import_store_array_objects.pl
    );
    
    foreach my $current_program (@programs_expected_in_path) {
      system("which $current_program > /dev/null");
      if ($?) {
        $error_msg .= "Can't find $current_program in path. These scripts are part of the ensembl-funcgen repository.\n";
      }
    }

    my $mysql_bin = `which mysql`;
    if (! $mysql_bin) {
      $error_msg .= "Can't find mysql command in path. Mysql has to be in the path.\n";
    }
    if ($mysql_bin eq '/usr/bin/mysql') {
      $error_msg .= "If you are on the sanger farm, you are using the wrong mysql binary. This will not work with probe2transcript. Please use the one in /software/ensembl/central/bin by running\nexport PATH=/software/ensembl/central/bin/:\$PATH\n";
    }

    system("which db_cmd.pl > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find db_cmd.pl command in path. db_cmd.pl is a script from Ensembl ehive. Please make sure it is in the PATH.\n";
    }
    
    system("which sequence_dump.pl > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find db_cmd.pl command in path. db_cmd.pl is a script from ensembl-analysis. Please make sure it is in the PATH.\n";
    }
    
    system("which exonerate > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find exonerate in path. On the sanger farm there is an installation in /software/ensembl/compara/exonerate/. Please make sure exonerate is in the PATH.\n";
    }
    
    system("which dump_genes.pl > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find dump_genes.pl in path. dump_genes.pl is a script from scripts/export/ in the ensembl-funcgen repository. Please make sure it is in the PATH.\n";
    }
    
    system("which bedtools > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find bedtools in path. On the sanger farm there is an installation in /software/ensembl/funcgen/. Please make sure the bedtools binary is in the PATH.\n";
    }    

    #my $funcgen_dba_hash = $self->param('funcgen_dba_hash');
    
    # The database might not exist yet. Deleting this key makes the 
    # DBConnection object not connect to it.
    #
    #delete $funcgen_dba_hash->{-dbname};
    
    my $dbc;
    
#     $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$funcgen_dba_hash);
#     
#     eval {
#       $dbc->connect();
#      };
#      if ($@) {
#        $error_msg .= "Couldn't connect to the database server defined in funcgen_dba_hash. Please review the connection details.\n\n";
#        $error_msg .= "The error message was:\n\n$@";
#      }
#     my $tracking_dba_hash = $self->param('tracking_dba_hash');
#     
#     # The database might not exist yet. Deleting this key makes the 
#     # DBConnection object not connect to it.
#     #
#     delete $tracking_dba_hash->{-dbname};
#     
# #     $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
#     $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
#     
#     eval {
#       $dbc->connect();
#      };
#      if ($@) {
#        $error_msg .= "Couldn't connect to the database server defined in tracking_dba_hash. Please review the connection details.\n\n";
#        $error_msg .= "The error message was:\n\n$@";
#      }
#     
#     if ($error_msg) {
#       die($error_msg);
#     }
#     
#     use Bio::EnsEMBL::DBSQL::DBAdaptor;
#     eval {
#       my $dnadb =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
# 	  -dbname          => $tracking_dba_hash->{'-dnadb_name'},
# 	  -host            => $tracking_dba_hash->{'-dnadb_host'},
# 	  -port            => $tracking_dba_hash->{'-dnadb_port'},
# 	  -user            => $tracking_dba_hash->{'-dnadb_user'},
# 	  -species         => $tracking_dba_hash->{'-species'}, 
#       );
#     };
#     if ($@) {
#        $error_msg .= "Couldn't connect to the dnadb database defined in tracking_dba_hash. Please review the connection details.\n\n";
#        $error_msg .= "The error message was:\n\n$@";
#     }
    if ($error_msg) {
      die($error_msg);
    }
    $logger->info("Pre pipeline checks have completed successfully. You can now run the rest of the pipeline.\n");
}

1;
