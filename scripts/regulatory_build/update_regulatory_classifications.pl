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

update_regulatory_classifications.pl

=head1 SYNOPSIS



=head1 DESCRIPTION

This script update the regulatory_feature feature_type_id values based on the output of the 
regulatory classification pipeline. It also produces a simple report which counts the transitions
of feature_types between the previous and current builds.

=head1 OPTIONS

 Mandatory
  -species    Latin name as used in DB name or meta table e.g. homo_sapiens



 Optional



=head1 EXAMPLE


=head1 SEE ALSO



=cut



#To do
# 1 Add new/old DB support. To enable reports when old/archived feature sets have been overwritten/removed in the dev DB
# 2 Add -cell_types -skip_cell_types params

use strict;

use Pod::Usage;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
$| = 1; # auto flush stdout


#Helper params
#declared here to avoid single usage warning
$main::_log_file = undef;
$main::_tee      = 0;    #tee the output of this scripts


my ($host, $user, $port, $pass, $dbname, $species,
	$dnadb_host, $dnadb_user, $dnadb_port, $dnadb_pass, $dnadb_name,
	$dbprefix, $overwrite, $report_only, $skip_report);


my @tmp_args = @ARGV;

GetOptions(
		   'host=s'       => \$host,
           'user=s'       => \$user,
           'port=i'       => \$port,
           'pass=s'       => \$pass,
           'dbname=s'     => \$dbname,
		   'dnadb_host=s' => \$dnadb_host,
           'dnadb_user=s' => \$dnadb_user,
           'dnadb_port=i' => \$dnadb_port,
           'dnadb_pass=s' => \$dnadb_pass,
           'dnadb_name=s' => \$dnadb_name,
           'species=s'    => \$species,
		   'dbprefix=s'   => \$dbprefix,
		   'overwrite'    => \$overwrite,
		   'report_only'  => \$report_only,#In case we have already updated
		   'skip_report'  => \$skip_report,#In case of no stable IDs

		   #cell_types
		   #skip_cell_types

		   
		   #Helper params 
		   'tee'          => \$main::_tee,
		   'logfile=s'      => \$main::_log_file,
		   
 
           'help'         => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);



my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;

$Helper->log_header("Running on $0");
$Helper->log("Params are:\t@tmp_args");



# TEST MANDATORY PARAMS

if($skip_report && $report_only){
  die('You have specified mutually exclusive parameter -skip_report and -report_only. Please chose one or the other');
}


$dbprefix ||= 'annotation_'.$dbname;

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host    => $host,
													  -user    => $user,
													  -pass    => $pass,
													  -port    => $port,
													  -dbname  => $dbname,
													  -species => $species,
													  
													  #dnadb params
													  -dnadb_host => $dnadb_host,
													  -dnadb_user => $dnadb_user,
													  -dnadb_pass => $dnadb_pass,
													  -dnadb_port => $dnadb_port,
													  -dnadb_name => $dnadb_name,
													   

													  );

#Test connections
$db->dbc;
$db->dnadb->dbc;

my $fset_adaptor = $db->get_FeatureSetAdaptor;


#Get all current and archive regulatory feature sets

my $sql ="SELECT name from feature_set where name like 'RegulatoryFeatures:%' and name not rlike 'RegulatoryFeatures:.*v[0-9]+\$'";
my @fset_names = @{$db->dbc->db_handle->selectcol_arrayref($sql)};
my (%fsets, @rows);

foreach my $fset_name(@fset_names){

  my $fset = $fset_adaptor->fetch_by_name($fset_name);

  next if $fset->cell_type->name eq 'MultiCell';

  if(! $fset){
	die("Failed to fetch current FeatureSet:\t$fset_name");
  }


  $sql = "SELECT name from feature_set where name like '${fset_name}_v%'";
  @rows = @{$db->dbc->db_handle->selectcol_arrayref($sql)};

  my $arch_fset;

  if(scalar(@rows) == 0){
	warn "Could not find archived FeatureSet for $fset_name\n";
  }
  elsif(scalar(@rows) < 1){
	#fail report this instead of dying here?
	die("Could not find unique archived FeatureSet for $fset_name:\t@rows\nPlease remove all but the latest");
  }
  else{
	$arch_fset = $fset_adaptor->fetch_by_name($rows[0]);
  
	if(! $arch_fset){
	  die("Failed to fetch archive FeatureSet:\t".$rows[0]);
	}
  }

  $fsets{$fset->name} = {
						 new => $fset,
						 old => $arch_fset,
						};
}


my $annot_db;


$Helper->log_header('Updating RegulatoryFeature Classifications') if ! $report_only;

my @fail_report;


foreach my $fset_name(keys %fsets){
  my $nfset = $fsets{$fset_name}{new};
  my $ofset = $fsets{$fset_name}{old};
  
  ($annot_db = "${dbprefix}_".$nfset->cell_type->name) =~ s/-/_/g;
  $sql = "SHOW databases like '${annot_db}'";
  @rows = @{$db->dbc->db_handle->selectcol_arrayref($sql)};
  
  if(! @rows){
	die("Could not find $annot_db DB\nPlease copy regulatory classification DBs to $host or speficy correct -dbprefix");
  }


  if(! $report_only){
	$Helper->log("\n"); #Just to get some nice spacing in the log
	$Helper->log('Updating '.$nfset->name);

	#Check if we already have annotations

	$sql = 'SELECT count(*) from regulatory_feature where feature_set_id='.$nfset->dbID.
	  ' AND feature_type_id!='.$nfset->feature_type->dbID;
	my ($count) = $db->dbc->db_handle->selectrow_array($sql);

	if($count){

	  if(! $overwrite){
		$Helper->log("Found $count existing classifications");
		$Helper->log('Skipping '.$nfset->name.' update');

		push @fail_report, 'Failed to update '.$nfset->name.
		  " due to $count existing classifications\nPlease specify -overwrite";
		#delete $fsets{$fset_name}; Don't skip the report as we have logged to fail report
		next;
	  }
	  else{
		$Helper->log("Overwriting $count existing classifications");
	  }
	}

	#Test for regulatory_features_classified table
	$sql = "use $annot_db";
	$db->dbc->do($sql);
	$sql = "show tables like 'regulatory_features_classified'";
	my @tables = @{$db->dbc->db_handle->selectcol_arrayref($sql)};

	if(scalar(@tables) != 1){
	  my $message = "Failed to find $annot_db.regulatory_features_classified\nSkipping $fset_name update and report\n";
	  $Helper->log_error($message);
	  push @fail_report, $message;
	  delete $fsets{$fset_name};
	  next;
	}


	$sql = "update ${annot_db}.regulatory_features_classified crf, ${dbname}.regulatory_feature rf ".
	  "SET rf.feature_type_id=crf.feature_type_id where crf.regulatory_feature_id=rf.regulatory_feature_id";
	
	my $row_cnt = $db->dbc->do($sql);
	$Helper->log("Updated $row_cnt regulatory_feature classifications");
	#Will this return the update count
	
	$sql = "use $dbname";
	$db->dbc->do($sql);
  }

  #Now check for any which have been missed
  
  $sql = 'SELECT count(*) from regulatory_feature where feature_set_id='.$nfset->dbID.
	' AND feature_type_id='.$nfset->feature_type->dbID;
  my ($count) = $db->dbc->db_handle->selectrow_array($sql);

  if($count){

	$Helper->log_error("Found $count $fset_name entries with no classification");
	#build fail report
	delete $fsets{$fset_name};
	push @fail_report, 'Failed to report '.$nfset->name.
	  " due to $count entries with no classifications\n".
		"Please specify rectify manually\nOr specify -overwrite or remove -report_only to update the classifications";
  }
}


#Now do reports
if($skip_report){
  $Helper->log_header("Skipping classification report");
}
else{

  
  $Helper->log_header("Analyzing/Optimizing tables");
  
  #old vs new annotation can only be done if stable ID mapping has been performed.
  foreach my $table(qw(regulatory_feature feature_type feature_set)){
	$Helper->log("\t$table");
	
	$sql = "ANALYZE table $table";
	$db->dbc->do($sql);
	$sql = "OPTIMIZE table $table";
	$db->dbc->do($sql);
  }
  
  $Helper->log("\n");
  
  $sql = "DROP TABLE IF EXISTS tmp_new_vs_old_reg_feature_types";
  $db->dbc->do($sql);
  $sql = "CREATE table tmp_new_vs_old_reg_feature_types(
	`stable_id` mediumint(8) unsigned DEFAULT NULL,
	`old_feature_type_name`  varchar(100) DEFAULT NULL,
	`new_feature_type_name`  varchar(100) DEFAULT NULL,
	`feature_set_id`  int(10) unsigned NOT NULL,
	PRIMARY KEY (`stable_id`, `feature_set_id`),
    KEY (`feature_set_id`),
    KEY old_feature_type_name (`old_feature_type_name`),
  KEY new_feature_type_name (`new_feature_type_name`)
)  ENGINE=MyISAM DEFAULT CHARSET=latin1;";
  $db->dbc->do($sql);

  

  foreach my $fset_name(keys %fsets){
	my $nfset = $fsets{$fset_name}{new};
	my $ofset = $fsets{$fset_name}{old};
	
	if(! $ofset){
	  warn 'Archive set not found for '.$nfset->name.". Skipping report\n";
	  next;
	}
	
	#Test for NULL stable IDs!
	$sql = 'SELECT count(*) from regulatory_feature where feature_set_id='.
	$nfset->dbID.' AND stable_id is NULL';
	my ($count) = $db->dbc->db_handle->selectrow_array($sql);
	
	if($count){
	  my $msg = "Found $count $fset_name entries with no stable IDs";
	  $Helper->log_error($msg);
	  push @fail_report, $msg."\nSkipping report for ".$nfset->name;
	  next;
	}
	
	
	$Helper->log('Building report info for '.$nfset->name);
	
	#Store new feature_set_id with old annotattion
	$sql = 'INSERT into tmp_new_vs_old_reg_feature_types SELECT stable_id, ft.name, NULL, '.
	  $nfset->dbID.' from regulatory_feature rf, feature_type ft '.
		' WHERE rf.feature_type_id=ft.feature_type_id AND feature_set_id='.$ofset->dbID;
	$db->dbc->do($sql);
	
	#Update the old reg feat records
	$sql = 'UPDATE tmp_new_vs_old_reg_feature_types t, regulatory_feature rf, feature_type ft set t.new_feature_type_name=ft.name '.
	  ' WHERE rf.feature_type_id=ft.feature_type_id AND rf.stable_id=t.stable_id and rf.feature_set_id=t.feature_set_id '.
		' and rf.feature_set_id='.$nfset->dbID;
	$db->dbc->do($sql);
	
	#Add the new reg feat records
	$sql = 'INSERT IGNORE into tmp_new_vs_old_reg_feature_types SELECT stable_id, NULL, ft.name, '.
	  $nfset->dbID.' from regulatory_feature rf, feature_type ft where rf.feature_type_id=ft.feature_type_id AND feature_set_id='.$nfset->dbID;
	$db->dbc->do($sql);
	
	
	#Check for NULL feature_type names here?
	#No these may be valid new/deleted features
	
  }

  #Now do the report
  $sql = 'SELECT fs.name, t.old_feature_type_name, t.new_feature_type_name, count(*) from tmp_new_vs_old_reg_feature_types t, feature_set fs '.
  ' WHERE t.feature_set_id=fs.feature_set_id '.
	' GROUP by fs.name, t.old_feature_type_name, t.new_feature_type_name';
  #warn $sql;#May need to explain this if it is not performant
  @rows = @{$db->dbc->db_handle->selectall_arrayref($sql)};
  

  $sql = "DROP TABLE tmp_new_vs_old_reg_feature_types";
  $db->dbc->do($sql);
  
  
  my $header = ['FeatureSet', 'Old FeatureType', 'New FeatureType', 'Count'];
  unshift @rows, $header;


  $Helper->log_header('RegulatoryFeature Classification Report');
  
  foreach my $row(@rows){
	
	my $count = pop @$row;
	
	my $line = '';
	
	foreach my $i(0..$#{$row}){
	  
	  if(! defined $row->[$i]){
		
		if($i == 1){
		  $row->[$i] = 'New';
		}
		elsif($i == 2){
		  $row->[$i] = 'Deleted';
		}
		else{
		  die("Cannot have NULL value for field ".$header->[$i]."\n");
		}
	  }
	  
	  $line .= sprintf ('%-40s', $row->[$i]);
	}
	
	$Helper->log($line.sprintf ('%-8s', $count));
  }
}

if(@fail_report){
  $Helper->log_header('Failure Report');
  map $Helper->log_error($_), @fail_report;
  #This will require Helper check in, which will screw rollback stuff?
  exit;
}

1;
