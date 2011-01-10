#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

 ensembl-efg rollback_experiment.pl

=head1 SYNOPSIS

 rollback_experiment.pl [options]

=head1 OPTIONS

 Mandatory
  -experiment|e    Name of array based experiment
  -feature_set     FeatureSet name
  -chip_ids|c      List of ExperimentalChip unique IDs (comma separated, no spaces)
  -pass|p          The MySQL password
  -dbname|n        Defines the eFG dbname if it is not standard
  -port            The port for the MySQL instance
  -host            The MySQL host
  -user|u          The MySQL user name.
  -full_delete|f   Performs a full delete, removing 'non-complex' feature and result sets associated with this experiment.
  -force_delete|d  Forces a full delete of all experimental information, even if it is part of a combinved data set.
  #-result_set      Name to give the raw/normalised result set.
  -help            Brief help message
  -man             Full documentation


=head1 DESCRIPTION

B<This program> removes all the imported data for a given experiment name.

=cut


#To do
# 1 Can we implement some of the Helper rollback methods here?
# 2 Handle Experiments with no ResultSets
# 3 Move this to Helper and make this handle Experiments, FeatureSets or ResultSets etc..
# 4 Needs updating to handle ExperimentalSets

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($chips, $pass, $full_delete, @chips, $fset_name);
my ($exp_name, $host, $dbname, $log_msg, $species, $force_delete);
my $port = $ENV{'EFG_PORT'};
my $user = $ENV{'EFG_WRITE_USER'};

my @tmp_args=@ARGV;

warn "@tmp_args";

GetOptions (
			"experiment|e=s"      => \$exp_name,
			"feature_set=s"       => \$fset_name,
			"chip_ids|c=s"        => \@chips,
			"pass|p=s"            => \$pass,
			"port=s"              => \$port,
			"dbname|n=s"          => \$dbname,
			"host=s"              => \$host,
			"species=s"           => \$species,
			"user|u=s"            => \$user,
			"full_delete|f"       => \$full_delete,
			"force_delete|d"      => \$force_delete,
			#"data_version|d=s"   => \$data_version,
			"help|?"              => sub { pos2usage(-exitval => 0, 
													 -message => "Params are:\t@tmp_args"); },
			"man|m"               => sub { pos2usage(-exitval => 0, 
													 -verbose => 2,
													 -message => "Params are:\t@tmp_args"); },
		   ) or pod2usage(
							 -exitval => 1,
							 -message => "Params are:\t@tmp_args"
							);

$port ||=3306;
throw('Must define a -user paramter') if ! $user;
throw('Must define a -dbname parameter') if ! $dbname;
throw('Must define a -host parameter') if ! $host;
throw('Must define a -pass parameter') if ! $pass;

if($fset_name && $exp_name){
  throw('Cannot sepcify -experiment and -feature_set, please choose one');
}
elsif(! ($fset_name || $exp_name)){
  throw('Must define an -experiment|feature_set parameter');
}

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

#Need to add dnadb_host as parameter
#my $dnadb =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#												 -host => 'ens-staging',
#												 -dbname => 'homo_sapiens_core_55_37',
#												 -user => 'ensro',
#												 #-pass => $pass,
#												 -port => $port,
#												 -species => $species,
#												);


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $host,
													  -dbname => $dbname,
													  -user => $user,
													  -pass => $pass,
						   						      -port => $port,
													  -species => $species,
													  #-dnadb => $dnadb,
													 );


my $exp_a = $db->get_ExperimentAdaptor();
my $rset_a = $db->get_ResultSetAdaptor();
my $dset_a = $db->get_DataSetAdaptor();
my $ec_a = $db->get_ExperimentalChipAdaptor();
my $fset_adaptor = $db->get_FeatureSetAdaptor();
my $helper = Bio::EnsEMBL::Funcgen::Utils::Helper->new(
													   no_log => 1,#default to STDOUT
													  );


#Need to move this to Helper.pm
if ($fset_name){
  
  my $fset = $fset_adaptor->fetch_by_name($fset_name);

  throw("Could not find FeatureSet $fset_name for rollback") if ! $fset;
  
  $helper->rollback_FeatureSet($fset, $force_delete, undef, $full_delete);
  
  #This currently just rollsback anntotated features, xrefs, states etc
  #Full delete would remove records from feature_set? (supporting sets and experiment)



}
else{#array experiment
  die("This script needs to be updated and is currently unsafe for array experiment roll backs");

  my $exp = $exp_a->fetch_by_name($exp_name);
  throw("Experiment $exp_name does not exist in the database") if ! defined $exp;
  
  #do chips belong to experiment?
  if(@chips){
	

	if($force_delete){
	  die "Cannot force delete when restricting to a chip list";
	}
	
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
  my ($sql, %cc_ids, %no_delete_cc_ids, @simple_rsets, @rollback_rsets, %rollback_dsets);
  


  #The only way we can roll back just individual result_sets is to protect the IMPORT sets
  #then we can reconfigure other result_sets if we want to.
  #Not yet implemented
  
  #we need to check whether any of these rsets are used in combined data_sets using other experiment data
  #We then need log their cc_ids, skip the result delete and if we encounter these cc_ids in another rset, then we only delete the cc records pertaining to that rset.
  
  
  foreach my $rset(@rsets){

	print "\n::\tChecking ResultSet:\t".$rset->name()."\n";
	my $simple_rset = 1;
	
	foreach my $dset(@{$dset_a->fetch_all_by_supporting_set($rset)}){
	  my @dsets;
	  
	  foreach my $d_rset(@{$dset->get_supporting_sets()}){
		#implicit that this will always be a ResultSet
		
		my $rset_id = $d_rset->dbID();
		
		if(! grep/$rset_id/, @rset_ids){
		  #Contains other rset i.e. combined data set
		  push @dsets, $dset;
		}
	  }
	  
	  
	  if (scalar(@dsets) == 0){
		#Not a complex dset, so can delete
		
		#But is this accounting for the specified chip_ids?
		#We need to filter these dsets based on the chip ids!!
		
		if(! exists $rollback_dsets{$dset->dbID}){
		  print "\n::\tIdentified DataSet for removal:\t".$dset->name()."\n";
		  $rollback_dsets{$dset->dbID} = $dset;
		}
	  }
	  else{
		#Found complex dset

		if($force_delete){
		  die "force delete not yet implemented";
		  #Is this sensible to remove underlying data from a combined data set?
		  $simple_rset = 1;
		}
		else{
		  $simple_rset = 0;
		  
		  print "::\tSkipping delete of ResultSet ".$rset->name." as it is used in the combined DataSets:\n\t"
			.join(', ', (map $_->name, @dsets))."\n";
		  
		  map {$no_delete_cc_ids{$_} = 1} @{$rset->chip_channel_ids};
		}
	  }
	}
	
	#Have to do this here as we may not have a DataSet for a given ResultSet
	push @simple_rsets, $rset if $simple_rset;
  }


  
  #Now we have an nr hash of data sets to remove, an nr list of result sets to remove and a nr hash of chip_channel_ids we don't want to remove from the result table
  #We haven't yet accunted for any user specified chip_ids

  #Now filter simple rsets for chip IDs
  my ($remove_rset, $remove_cc);

  foreach my $rset (@simple_rsets) {
	$remove_rset = 1;

	foreach my $ec (@{$rset->get_ExperimentalChips()}) {
	  $remove_cc = 1;
	
	  #Is this working?
	  if (@chips) {				#delete only @chips
		my $uid = $ec->unique_id();

		if (! grep/$uid/, @chips) {	#other chips present, don't remove rset/cc
		  print "::\tResultSet contains other ExperimentalChip:\t".$ec->unique_id()."\n";
		  $remove_cc = 0;
		  $remove_rset = 0;
		}
	  }

	  if ($remove_cc) {
		#Need to filter here on no_delete_cc_ids
	
		if ($rset->table_name eq 'experimental_chip') {
		  $cc_ids{$rset->get_chip_channel_id($ec->dbID())} = $ec;
		} elsif ($rset->table_name eq 'channel') {

		  foreach my $chan (@{$ec->get_Channels()}) {
			$cc_ids{$rset->get_chip_channel_id($chan->dbID())} = $chan;
		  }
		} else {
		  throw('rollback_experiment.pl does not yet accomodate non-chip roll backs');
		}
	  }
	}
	
	#clean result, chip_channel, experimental_chip and remove rset if required.
	#do in staged delete, otherwise pseudo sets will fail as there will be no r with corresponding cc_id?
	#do we have rsets on the pseuod level?  I think not?
	#channels will be deleted by association with rsets ecs
	#we may have an ec which is not part of an rset?
	#so delete separately?


	#Remove corresponding data_sets from hash if we are keeping this result_set
	if (! $remove_rset) {

	  foreach my $dset (@{$dset_a->fetch_all_by_supporting_set($rset)}) {
		delete $rollback_dsets{$dset->dbID} if exists  $rollback_dsets{$dset->dbID};
	  }
	}


	if (! keys %cc_ids) {
	  print "::\tResultSet does not contain specified ExperimentalChips\n";
	} else {					#we have something to delete

	  if (! $remove_rset) {
		print "::\tOther ExperimentalChips persist, skipping ResultSet delete for:\t".$rset->name()."\n";
	  } else {
		push @rollback_rsets, $rset;
	  }
	}
  }



  #Now remove filtered Features/DataSets
  if($full_delete){
	
	foreach my $dset(values %rollback_dsets){
	  
	  #delete feature_set first, so we don't ever have an orphaned feature_set
	  my $fset = $dset->product_FeatureSet();
	  
	  if(defined $fset){
	  print "::\tDeleting FeatureSet:\t".$fset->name()."\n";
	  
	  #delete status entries (should we do this first?)
	  $sql = 'DELETE from status where table_name="feature_set" and table_id='.$fset->dbID();
	  $db->dbc->do($sql) || throw("Failed to delete status entries for feature_set with dbID:\t".$fset->dbID());

	  	  
	  $sql = 'DELETE from '.$fset->type.'_feature where feature_set_id='.$fset->dbID();
	  $db->dbc->do($sql) || throw("Failed to delete ".$fset->type." features with feature_set dbID:\t".$fset->dbID());
	  
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



### Delete all results and chip_channel records in one go

#filter cc_ids first
foreach my $key(keys %no_delete_cc_ids){
  delete $cc_ids{$key} if exists $cc_ids{$key};
}

if(keys %cc_ids){
  
  print "::\tDeleting result, chip_channel and result_set records for:\t".join(', ', (map {$_->name} @rollback_rsets))."\n";
  
  #We could do with testing for the cc_ids here
  #So we can rollback an Rset if it have already been partially mangled
  $sql = "DELETE from result where chip_channel_id IN (".join(', ', keys %cc_ids).")";
  $db->dbc->do($sql);
  $sql = "DELETE from chip_channel where result_set_id IN (".join(', ', (map {$_->dbID} @rollback_rsets)).")";
  $db->dbc->do($sql);
  $sql = "DELETE from status where table_name='result_set' and table_id IN (".join(', ', (map {$_->dbID} @rollback_rsets)).")";
  $db->dbc->do($sql);
  $sql = "DELETE from result_set where result_set_id IN (".join(', ', (map {$_->dbID} @rollback_rsets)).")";	
  $db->dbc->do($sql);
}


#Where are we deleting the resultsets?

#now do final clean up delete channels, ec, and experiemnt (inc mage_xml...and any other linked entries)
#We should only delete those in the cc_ids hash as these have been filtered appropriately
#We will get some redundant activity here as were deleting the underlying channels explicitly
#Can we guarantee that the channels will deleted automatically? I think so.  Therefore skip channel elements

print "\n\n";

foreach my $cc(values %cc_ids){

  next if ! $cc->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip');

  #channels first
  foreach my $chan(@{$cc->get_Channels()}){
	print "::\tDeleting channel records for ExperimentalChip:\t".$cc->unique_id().":".$chan->dye()."\n";

	#and status entries
	$sql = 'DELETE from status where table_name="channel" and table_id='.$chan->dbID();
	$db->dbc->do($sql) || throw("Failed to delete status entries for channel with dbID:\t".$chan->dbID());
	  
	$sql = 'DELETE from channel where channel_id='.$chan->dbID();
	$db->dbc->do($sql) || throw("Failed to delete channel with dbID:\t".$chan->dbID());
	}

  print "::\tDeleting experimental_chip records for ExperimentalChip:\t".$cc->unique_id()."\n";
  #and status entries
  $sql = 'DELETE from status where table_name="experimental_chip" and table_id='.$cc->dbID();
  $db->dbc->do($sql) || throw("Failed to delete status entries for experimental_chip with dbID:\t".$cc->dbID());

  #now chip
  $sql = 'DELETE from experimental_chip where experimental_chip_id='.$cc->dbID();
  $db->dbc->do($sql) || throw("Failed to delete experimental_chip with dbID:\t".$cc->dbID());
}




#Delete experiment only if it is now empty.
$exp = $exp_a->fetch_by_name($exp->name);

if(! @{$exp->get_ExperimentalChips}){
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
  print "::\tWARNING:\tSkipping full experiment delete as some data still persists.  The xml has not been changed to reflect the experimental chips you have deleted from this experiment\n::\tPlease update mage_xml manually\n";
}

}
