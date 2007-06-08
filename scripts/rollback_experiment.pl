#!/usr/local/ensembl/bin/perl



### add pod here


use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($chips, $pass, $full_delete, @chips);
my ($exp_name, $host, $dbname, $help, $man, $log_msg);
my $port = 3306;
my $user = 'ensadmin';

GetOptions (
			"experiment|e=s"      => \$exp_name,
			"chip_ids|c=s"        => \@chips,
			"pass|p=s"            => \$pass,
			"port=s"              => \$port,
			"dbname|n=s"          => \$dbname,
			"host|h=s"            => \$host,
			"user|u=s"            => \$user,
			"full_delete|f"     => \$full_delete,
			#"data_version|d=s"   => \$data_version,
			"help|?"              => \$help,
			"man|m"               => \$man,
		   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

### Set up adaptors and FeatureSet 

throw('Must define a -dbname parameter') if ! $dbname;
throw('Must define a -dbhost parameter') if ! $host;
throw('Must define a -pass parameter') if ! $pass;
throw('Must define an -experiment_name parameter') if ! $exp_name;

$log_msg = "::\tRolling back experiment:\t$exp_name\n".
  "::\tOn:\t${host}:${port}:${dbname}\n";



if(@chips){
  @chips = split/,/, join('', @chips);
  $log_msg .= "::\tPerforming chip only roll back\n::\tChip IDs:\t".join("\t", @chips)."\n";
}else{
  $log_msg .= "::\tPerforming full chip roll back, all experiment data will be lost\n";
}

$log_msg .= "::\t";
$log_msg .= ($full_delete) ? "Removing all associated non-complex data/feature_sets\n" :
  "All associated data/feature_sets will persist\n";

print $log_msg;


#we need to warn if we have a data_set attached to any of the result_sets
#do not delete data sets, just provide info about orphaned feature/data_sets?
#we might not want to remove a feature set as it may be art of a combined experiment analysis

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $host,
													  -dbname => $dbname,
													  -user => $user,
													  -pass => $pass,
						   						      -port => $port,
													 );


my $exp_a = $db->get_ExperimentAdaptor();
my $rset_a = $db->get_ResultSetAdaptor();
my $dset_a = $db->get_DataSetAdaptor();
my $ec_a = $db->get_ExperimentalChipAdaptor();
my $fset_adaptor = $db->get_FeatureSetAdaptor();


my $exp = $exp_a->fetch_by_name($exp_name);
throw("Experiment $exp_name does not exist in the database") if ! defined $exp;

#do chips belong to experiment?
if(@chips){
  
  foreach my $chip(@chips){
	
	my $tmp_ec = $ec_a->fetch_by_unique_and_experiment_id($chip, $exp->dbID());

	throw("ExperimentalChip $chip is not part of the Experiment $exp_name") if ! defined $tmp_ec;
	#we could add a list of the ecs present here

  }
}


# we also need to remove any result sets which are entirely consituted by the achips
# also need to log which rsets have been removed

my %table_syns = (
				  channel => 'c',
				  experimental_chip => 'ec',
				 );

#could maybe get all rsets first, get all ids, so we can compare to linked data_sets
#skip delete data/feature_set if it is linked to another rset
my @rsets = @{$rset_a->fetch_all_by_Experiment($exp)};
my @rset_ids = map $_->dbID(), @rsets;
my ($sql, %shared_dsets);

foreach my $rset(@rsets){
  print "\n::\tChecking ResultSet:\t".$rset->name()."\n";
  my ($syn, %cc_ids);
  my $remove_rset = 1;

  foreach my $ec(@{$rset->get_ExperimentalChips()}){

	my $remove_cc = 1;
	

	if(@chips){#delete only @chips
	  my $uid = $ec->unique_id();

	  if(! grep/$uid/, @chips){#other chips present, don't remove rset/cc
		print "::\tResultSet contains other ExperimentalChip:\t".$ec->unique_id()."\n";
		$remove_cc = 0;
		$remove_rset = 0;
	  }
	}


	if($remove_cc){

	  if($rset->table_name eq 'experimental_chip'){
		$cc_ids{$ec->dbID()} = $rset->get_chip_channel_id($ec->dbID());
	  }
	  elsif($rset->table_name eq 'channel'){

		foreach my $chan(@{$ec->get_Channels()}){
		  $cc_ids{$chan->dbID()} = $rset->get_chip_channel_id($chan->dbID());
		}
	  }
	  else{
		throw('rollback_experiment.pl does not yet accomodate none-chip roll backs');
	  }
   	}
  }
	
  #clean result, chip_channel, experimental_chip and remove rset if required.
  #do in staged delete, otherwise pseudo sets will fail as there will be no r with corresponding cc_id?
  #do we have rsets on the pseuod level?  I think not?
  #channels will be deleted by association with rsets ecs
  #we may have an ec which is not part of an rset?
  #so delete separately?

  if(! keys %cc_ids){
	print "::\tResultSet does not contain specified ExperimentalChips\n";
  }else{#we have something to delete

	#remove everything in reverse order to avoid orphaning records
	#and making it impossible to link to them in the event of a failure

	my $table_name = $rset->table_name();
	my $syn = $table_syns{$table_name};


	print "::\tDeleting result chip_channel records for $table_name chip_channel_ids:\t".join(', ', values %cc_ids)."\n";
	$sql = "DELETE from result where chip_channel_id IN (".join(', ', values %cc_ids).")";	
	$db->dbc->do($sql);

	
	if(! $remove_rset){
	  print "::\tOther ExperimentalChips persist, skipping ResultSet delete for:\t".$rset->name()."\n";
	}else{
	  
	  #remove data and feature sets if not linked to other experiments
	  if($full_delete){
		
		foreach my $dset(@{$dset_a->fetch_all_by_ResultSet($rset)}){
		  
		  #has dset been seen with a lin kto another exp?
		  if(! exists $shared_dsets{$dset->dbID()}){
			my $delete = 1;

			foreach my $d_rset(@{$dset->get_ResultSets()}){
			  my $rset_id = $d_rset->dbID();

			  if(! grep/$rset_id/, @rset_ids){
				$delete = 0;
				
				print "::\tSkipping delete of shared DataSet with dbID ".$dset->dbID."\n";
				last;
			  }
			}

			if($delete){
			  
			  #delete feature_set first, so we don't ever have an orphaned feature_set
			  my $fset = $dset->feature_set();

			  if(defined $fset){
				print "::\tDeleting FeatureSet:\t".$fset->name()."\n";

				#delete status entries (should we do this first?)
				$sql = 'DELETE from status where table_name="feature_set" and table_id='.$fset->dbID();
				$db->dbc->do($sql) || throw("Failed to delete status entries for feature_set with dbID:\t".$fset->dbID());

				$sql = 'DELETE from feature_set where feature_set_id='.$fset->dbID();
				$db->dbc->do($sql) || throw("Failed to delete feature_set with dbID:\t".$fset->dbID());
			  }

			  print "::\tDeleting DataSet:\t".$dset->name()."\n";

			  #dset status entries
			  $sql = 'DELETE from status where table_name="data_set" and table_id='.$dset->dbID();
			  $db->dbc->do($sql) || throw("Failed to delete status entries for data_set with dbID:\t".$dset->dbID());
			  #now delete data_set
			  $sql = 'DELETE from data_set where data_set_id='.$dset->dbID();
			  $db->dbc->do($sql) || throw("Failed to delete data_set with dbID:\t".$dset->dbID());
			}
		  }
		}
	  }

	  #remove rset last, so we don't potentially orphan the data_set
	  
	  print "::\tRemoving ResultSet:\t".$rset->name()."\n";
	  #delete status entries (should we do this first?)
	  $sql = 'DELETE from status where table_name="result_set" and table_id='.$rset->dbID();
	  $db->dbc->do($sql);
	  
	  #don't join tables as we have empty result set and would fail
	  $sql = 'DELETE rs, cc from result_set rs, chip_channel cc where rs.result_set_id='.$rset->dbID().
		' and cc.result_set_id='.$rset->dbID();
	  $db->dbc->do($sql);
	}

	#do this very last so we don't orphan the ResultSet
	#roll back just result and chip_channel first
	print "::\tDeleting chip_channel records for $table_name chip_channel_ids:\t".join(', ', values %cc_ids)."\n";
	$sql = "DELETE from chip_channel where chip_channel_id IN (".join(', ', values %cc_ids).")";	
	$db->dbc->do($sql);

  }
}

#now do final clean up delete channels, ec, and experiemnt (inc mage_xml...and any other linked entries)

print "\n\n";

foreach my $ec(@{$exp->get_ExperimentalChips()}){
  my $delete = 1;

  if(@chips){
	if(! grep/$ec->unique_id()/, @chips){
	  $delete = 0;
	}
  }

  if($delete){

	#channels first
	foreach my $chan(@{$ec->get_Channels()}){
	  print "::\tDeleting channel records for ExperimentalChip:\t".$ec->unique_id().":".$chan->dye()."\n";

	  #and status entries
	  $sql = 'DELETE from status where table_name="channel" and table_id='.$chan->dbID();
	  $db->dbc->do($sql) || throw("Failed to delete status entries for channel with dbID:\t".$chan->dbID());
	  
	  $sql = 'DELETE from channel where channel_id='.$chan->dbID();
	  $db->dbc->do($sql) || throw("Failed to delete channel with dbID:\t".$chan->dbID());
	}

	print "::\tDeleting experimental_chip records for ExperimentalChip:\t".$ec->unique_id()."\n";
	#and status entries
	$sql = 'DELETE from status where table_name="experimental_chip" and table_id='.$ec->dbID();
	$db->dbc->do($sql) || throw("Failed to delete status entries for experimental_chip with dbID:\t".$ec->dbID());
	#now chip
	$sql = 'DELETE from experimental_chip where experimental_chip_id='.$ec->dbID();
	$db->dbc->do($sql) || throw("Failed to delete experimental_chip with dbID:\t".$ec->dbID());
  }
}


#delete exp if we aren't doing just a chip level rollback
if(! @chips){
  
  #xml first
  if($exp->mage_xml_id()){
	print "::\tDeleting mage_xml for Experiment:\t".$exp->name()."\n";
	$sql = 'DELETE from mage_xml where mage_xml_id='.$exp->mage_xml_id();
	$db->dbc->do($sql) || throw("Failed to delete mage_xml with dbID:\t".$exp->mage_xml_id());
  }
  
  #and other tables? experimental_design_type, experimental_variable
  
  print "::\tDeleting Experiment:\t".$exp->name()."\n";
  #then exp
  $sql = 'DELETE from experiment where experiment_id='.$exp->dbID();
  $db->dbc->do($sql) || throw("Failed to delete experimentl with dbID:\t".$exp->dbID);
}else{
  print "::\tWARNING:\tThe xml has not been changed to reflect the experimental chips you have deleted from this experiment\n::\tPlease update mage_xml manually\n";
}
