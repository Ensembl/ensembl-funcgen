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


=head1 OPTIONS

 Mandatory
  -species    Latin name as used in DB name or meta table e.g. homo_sapiens



 Optional



=head1 EXAMPLE


=head1 SEE ALSO



=cut

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
	$dbprefix, $overwrite, $report_only);


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
		   'report_only'  => \$report_only,

		   #cell_types
		   #skip_cell_types

		   
		   #Helper params 
		   'tee'          => \$main::_tee,
		   'logfile'      => \$main::_log_file,
		   
 
           'help'         => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);



my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;

$Helper->log_header("Running on $0");
$Helper->log("Params are:\t@tmp_args");



# TEST MANDATORY PARAMS 




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

my $sql ="SELECT name from feature_set where name like 'RegulatoryFeatures:%' and name not rlike 'RegulatoryFeatures:.*v[0-9]+$'";
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

  if(scalar(@rows) != 1){
	die("Found more than one archived FeatureSet:\t@rows\nPlease remove all but the latest");
  }
  
  my $arch_fset = $fset_adaptor->fetch_by_name($rows[0]);

  if(! $arch_fset){
	die("Failed to fetch archive FeatureSet:\t".$rows[0]);
  }

  $fsets{$fset->name} = {
						 new => $fset,
						 old => $arch_fset,
						};
}


my $annot_db;


$Helper->log_header('Updating RegulatoryFeature Classifications') if $report_only;

my @fail_report;


foreach my $fset_name(keys %fsets){
  my $nfset = $fsets{$fset_name}{new};
  my $ofset = $fsets{$fset_name}{old};
  $Helper->log("\n"); #Just to get some nice spacing in the log

  $annot_db = "${dbprefix}_${dbname}_".$nfset->cell_type_name;
  $sql = "SHOW databases like '${annot_db}'";
  @rows = @{$db->dbc->db_handle->selectcol_arrayref($sql)};
  
  if(! @rows){
	die("Could not find $annot_db DB\nPlease copy regulatory classification DBs to $host");
  }


  if(! $report_only){
	$Helper->log('Updating '.$nfset->name);

	#Check if we already have annotations

	$sql = 'SELECT count(*) from regualtory_feature where feature_set_id='.$nfset->dbID.
	  ' AND feature_type_id!='.$nfset->feature_type_id;
	my $count = @{$db->dbc->db_handle->selectrow_array($sql)};

	if($count){

	  #build fail report unless overwrite

	  $Helper->log("Found $count existing classifications");
	  
	  if(! $overwrite){
		$Helper->log('Skipping update and report');

		push @fail_report, 'Failed to update/report '.$nfset->name.
		  " due to $count existing classifications\nPlease specify -overwrite";
		delete $fsets{$fset_name};
		next;
	  }
	}



	$sql = "update regulatory_features_classified.$annot_db crf, $dbname.regulatory_feature rf ".
	  "SET rf.feature_type_id=crf.feature_type_id where crf.regulatory_feature_id=rf.regulatory_feature_id";
	
	my $row_cnt = $db->dbc->do($sql);
	$Helper->log("Updated $row_cnt regulatory_feature classifications");
	#Will this return the update count

  }

  #Now check for any which have been missed
  
  $sql = 'SELECT count(*) from regualtory_feature where feature_set_id='.$nfset->dbID.
	' AND feature_type_id='.$nfset->feature_type_id;
  my $count = @{$db->dbc->db_handle->selectrow_array($sql)};

  if($count){

	$Helper->log("Found $count entries with no classification");
	#build fail report
	delete $fsets{$fset_name};
	push @fail_report, 'Failed to report '.$nfset->name.
	  " due to $count entries with no classifications\n".
		"Please specify rectify manually or specify -overwrite to try again";
  }
}


#Now do reports



#old vs new annotation can only be done if stable ID mapping has been performed.
foreach my $table(qw(regulatory_feature feature_type feature_set)){
  $sql = "ANALYZE table $table";
  $db->dbc->do($sql);
  $sql = "OPTIMIZE table $table";
  $db->dbc->do($sql);
}

$sql = "DROP TABLE IF EXISTS tmp_new_vs_old_reg_feature_types";
$db->dbc->do($sql);
$sql = "CREATE table tmp_new_vs_old_reg_feature_types(
	`stable_id` mediumint(8) unsigned DEFAULT NULL,
	`old_feature_type_name`  varchar(100) DEFAULT NULL,
	`new_feature_type_name`  varchar(100) DEFAULT NULL,
	`feature_set_id`  int(10) unsigned NOT NULL,
	PRIMARY KEY (`stable_id`),
    KEY (`feature_set_id`),
    KEY old_feature_type_name (`old_feature_type_name`),
  KEY new_feature_type_name (`new_feature_type_name`)
)  ENGINE=MyISAM DEFAULT CHARSET=latin1;";
$db->dbc->do($sql);

#Do we need to add keys for feature_types?
#We are going to group by feature_set.name, old_feature_type_name, new_feature_type_name



foreach my $fset_name(keys %fsets){
  my $nfset = $fsets{$fset_name}{new};
  my $ofset = $fsets{$fset_name}{old};

  $sql = 'INSERT into tmp_new_vs_old_reg_feature_types SELECT stable_id, ft.name, NULL, '.
	$nfset->dbID.' from regulatory_feature rf, feature_type ft '.
	  ' WHERE rf.feature_type_id=ft.feature_type_id AND feature_set_id='.$ofset->dbID;
  $db->dbc->do($sql);

  #Update the old reg feat records
  $sql = 'UPDATE tmp_new_vs_old_reg_feature_types t, regulatory_feature rf, feature_type ft set t.new_feature_type_name=ft.name '.
	' WHERE rf.feature_type_id=ft.feature_type_id AND rf.stable_id=t.stable_id and rf.feature_set_id='.$nfset->dbID;
  $db->dbc->do($sql);
  
  #Add the new reg feat records
  $sql = 'INSERT IGNORE into tmp_new_vs_old_reg_feature_types SELECT stable_id, NULL, ft.name, '.
	$nfset->dbID.' from regulatory_feature rf, feature_type ft where rf.feature_type_id=ft.feature_type_id AND feature_set_id='.$nfset->dbID;
  $db->dbc->do($sql);
}

#Now do the report
$sql = 'SELECT fs.name, t.old_feature_type_name, t.new_feature_type_name, t.count(*) from tmp_new_vs_old_reg_feature_types t, feature_set fs '.
  ' GROUP by fs.name, t.old_feature_type_name, t.new_feature_type_name';

@rows = @{$db->dbc->db_handle->selectall_arrayref($sql)};

unshift @rows, [qw(FeatureSet Old New Count)];


$Helper->log_header('Classification Report');

foreach my $row(@rows){
  
  my $count = pop @$row;

  my $line = '';

  foreach my $field(@$row){
	$line .= sprintf ('%-40s', $field);
  }

  $Helper->log($line.sprintf ('%-8s', $count));
}
