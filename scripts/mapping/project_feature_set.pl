#!/software/bin/perl

=head1 NAME

projection_feature_set.pl - Projects a feature set to a new genome assembly.

=head1 SYNOPSIS

project_feature_set.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  --is_component, --is-component      script is called from a wrapper script

  -i, --interactive                   run script interactively (default: true)
  -n, --dry_run, --dry                don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script uses Bio::EnsEMBL::Utils::AssemblyProjector to project features from an old 
assembly to a new assembly

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Nathan Johnson <njohnson@ebi.ac.uk>, based on the test_projection.pl script written by 
Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team.

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::AssemblyProjector;
use Bio::EnsEMBL::Utils::Exception qw( throw );

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => ""
);

$conf->parse_options
  (
   'host=s' => 1,
   'port=n' => 0,
   'user=s' => 1,
   'pass=s' => 0,
   'dbname=s' => 1,
   'old_cdbname=s' => 1,
   'cdbname=s' => 0,
   'cdbport=n' => 0,
   'cdbuser=s' => 0,
   'cdbpass=s' => 0,
   'cdbhost=s' => 0,
   'clobber' => 0,
   'slice=s' => 0,
   'species=s' => 0,
   'ignore_length' => 0,
   'force_store' => 1,
   'feature_set=s' => 1,
   'old_assembly=s' => 1,
   'new_assembly=s' => 1,
  );



#assign assemblies for regex embedding
my $new_assembly = $conf->param('new_assembly');
my $old_assembly = $conf->param('old_assembly');

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE    => $conf->param('logfile'),
  -LOGPATH    => $conf->param('logpath'),
  -LOGAPPEND  => $conf->param('logappend'),
  -VERBOSE    => $conf->param('verbose'),
);

# initialise log
$logger->init_log($conf->list_param_values);

# connect to database and get adaptors

#we need to make this optional so we can autogenerate the

my $old_cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host   => 'ensembldb',
   -user   => 'anonymous',
   -dbname => $conf->param('old_cdbname'),
   -group  => 'core',
  );

my ($cdb);

if(defined $conf->param('cdbname') || defined $conf->param('cdbhost') || defined $conf->param('cdbport') ){
 
  $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(
	 -host   => $conf->param('cdbhost') || $conf->param('host'),
	 -port   => $conf->param('cdbport') || $conf->param('port'),
	 -user   => $conf->param('cdbuser') || $conf->param('user'),
	 -pass   => $conf->param('cdbpass') || $conf->param('pass'),
	 -dbname => $conf->param('cdbname'),
	 -group  => 'core',
	);
}elsif(! defined  $conf->param('species')){

  ##????????????????????????
  throw('Must provide a species name if no core DB has been defined');
}

my $efg_db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  (
   -host   => $conf->param('host'),
   -port   => $conf->param('port'),
   -user   => $conf->param('user'),
   -pass   => $conf->param('pass'),
   -dbname => $conf->param('dbname'),
   -species =>  $conf->param('species'),
   #-group  => 'funcgen',#this is set by default no?
   -dnadb  => $old_cdb,
);



my $ap = Bio::EnsEMBL::Utils::AssemblyProjector->new
  (
   -OLD_ASSEMBLY    =>  $old_assembly,
   -NEW_ASSEMBLY    =>  $new_assembly,
   -CHECK_LENGTH    => 0,
   -MERGE_FRAGMENTS => 1,
  );


my $sa = $efg_db->get_SliceAdaptor;
my $new_sa = $cdb->get_SliceAdaptor;
my %setfeat_adaptor = (
					   'annotated'  => $efg_db->get_AnnotatedFeatureAdaptor(),
					   'regulatory' => $efg_db->get_RegulatoryFeatureAdaptor(),
					   'external'   => $efg_db->get_ExternalFeatureAdaptor(),
					  );

#can't dnadb set here as this will mess up the efg_db dnadb for the fetch on the old DB

my $fset_adaptor = $efg_db->get_FeatureSetAdaptor();
my $fset = $fset_adaptor->fetch_by_name($conf->param('feature_set'));

									   
if(! defined $fset){
  throw("Cannot findFeatureSet:\t".$conf->param('feature_set')."\n");
}



#Do this slice wise to avoid possibility of running out of memmory for large sets
#or we could pull everythign back based on the feature_set?
my @slices;

if($conf->param('slice')){
  @slices = ($sa->fetch_by_name($conf->param('slice')));

  warn('Need to test vs old assembly here');
}
else{
  #we can get toplevel for old assembly so we have to do for all?
  #Just do for chromosome now.
  warn 'Hardcoded for just chromosome slices on old assembly due to now toplevel';
  @slices = @{$sa->fetch_all('chromosome', $old_assembly)};


  if(! @slices){
	throw("There are no slices available for the old assembly $old_assembly. Maybe it is not present in the DB?");
  }
}

#should check if we can get coord_system on new assembly here, or will this be caught by mapper?


$logger->info("Projecting features on ".scalar(@slices)." slices from assembly $old_assembly to assembly $new_assembly\n");

my ($new_slice, $slice_name, $failed_cnt, $wrong_length_cnt, $old_slice);
my ($new_fslice, $old_fslice, $stored_cnt, @old_feats, @new_feats);
my $total_failed = 0;
my $total_wrong_length = 0;
my $total_stored = 0;

my $length_txt = ($conf->param('ignore_length')) ? 'features were skipped due to mismatched lengths in the new assembly' : 'features were stored despite mismatched lengths in the new assembly';



#This will not currently work as the seq_region_id are most likely going to be different between the current cdb
#and the db on which the features were loaded.
#unless we are already mapping for the correct seq_region_id for the assemmbly
#This is failing because the coord systems are not equal
#can we get around this by resetting the dnadb?
#this will most likely cause problems as the seq_region_ids won't be consistent
#we could simply use the old DB to retrieve them and then set the new DB when we are storing
#this should be fine so long as we are changing the slice?



foreach my $slice(@slices){
  $fset_adaptor->db->dnadb($old_cdb);
  throw('Not yet implemented projection segmant based clobber') if  $conf->param('clobber');
 
  #we need to try and pull back all features from projection segments first to see if there are already features projected
  #we could also print out the projection segments for each slice we are projecting to
  #could also have force_store, i.e. no clobber but load

  @new_feats = ();
  @old_feats = @{$fset->get_Features_by_Slice($slice)};
  $failed_cnt = 0;
  $wrong_length_cnt = 0;
  $stored_cnt = 0;

  ($slice_name = $slice->name) =~ s/$old_assembly/$new_assembly/;
  $new_slice = $new_sa->fetch_by_name($slice_name);#recreate old slice in new DB  
  $logger->info('Projecting '.scalar(@old_feats)." features from slice:\t".$slice->name()."\n");


  #do we have to switch the dnadb here?
  #yes we do as the internal caches only work for the current assembly/coordsystem
  $fset->adaptor->db->dnadb($cdb);
  @new_feats = @{$fset->get_Features_by_Slice($new_slice)};
  $logger->info('Found '.scalar(@new_feats)." features on new slice\t:".$new_slice->name()."\n");
  
  if(@new_feats){

	if($conf->param('clobber')){
	  my $sql = 'DELETE from '.$fset->type().'_feature WHERE  '.$fset->type().'_feature_id IN('.join(',', map $_->dbID(), @new_feats).')';
	  $logger->info("Clobbering features from FeatureSet:\t".$conf->param('feature_set')." on ".$new_slice->name."\n");
	  $efg_db->dbc_db_handle->do($sql);
	}
	elsif($conf->param('force_store')){
	  $logger->info('WARNING: Projecting '.scalar(@old_feats).' old features onto '.scalar(@new_feats)." new features\n");
	}else{
	  throw('Must define -clobber or -force_store');
	}
  }

  #reset dnadb and new_feats
  $fset->adaptor->db->dnadb($old_cdb);

  @new_feats = ();

  foreach my $old_feat(@{$fset->get_Features_by_Slice($slice)}){
	$old_fslice = $old_feat->feature_Slice();
	($slice_name = $old_fslice->name) =~ s/$old_assembly/$new_assembly/;
    $new_fslice = $new_sa->fetch_by_name($slice_name);#recreate old slice in new DB

	if(! $new_fslice){
	  $failed_cnt ++;
	  $logger->info('Failed to project feature with dbID '.$old_feat->dbID().' from '.$old_fslice->name()."\n");
	  next;
	}
	elsif($old_fslice->length() != $new_fslice->length()){
	  $logger->info('Found length mismatch for dbID '.$old_feat->dbID.":\t".$conf->param('old_assembly').' '.
					$old_fslice->length().' vs '.$conf->param('new_assembly').' '.$new_fslice->length()."\n");
	  next if ! $conf->param('ignore_length');
	}


	$old_feat->slice($new_slice);
	$old_feat->start($new_fslice->start);
	$old_feat->end($new_fslice->end);
	$old_feat->strand($new_fslice->strand);
	$old_feat->{'dbID'} = undef;
	$old_feat->{'adaptor'} = undef;

	#
	push @new_feats, $old_feat;
	$stored_cnt ++;
  }

  #need to reset dnadb here??????????????????????????
  $setfeat_adaptor{$fset->type}->db->dnadb($cdb);
  $setfeat_adaptor{$fset->type}->store(@new_feats);

  #test just to be sure
 
  my @tmp = @{$fset->get_Features_by_Slice($new_slice)};
		   

  $logger->info("Projected $stored_cnt features from ".$slice->name.' to new assembly '.
				$conf->param('new_assembly')."\n");

  $logger->info("$wrong_length_cnt $length_txt\n");
  $logger->info("$failed_cnt features failed to project to the new assembly\n");
  $logger->info('Retrieved '.scalar(@tmp).' features from new slice '.$new_slice->name()."\n");
  

  $total_failed += $failed_cnt;
  $total_wrong_length += $wrong_length_cnt;
  $total_stored += $stored_cnt;

}

$logger->info('Total projected feature to new assembly '.$conf->param('new_assembly').":\t$total_stored\n");
$logger->info("$total_wrong_length $length_txt\n");
$logger->info("$total_failed features failed to project to the new assembly\n");
# finish logfile
$logger->finish_log;

1;
