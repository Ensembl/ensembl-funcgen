#Nath

#Done
#Remove median & get_date and implement EFGUtils when migrating to eFG
#Implemented strand check and unmapped object for anti-sense mapping
#Changed transcript_probeset_count cache to use dbID instead of complete name
#Renamed check_names_match to validate_arrays and extended to perform probeset warning, external_db name guess and returns cache, user array param validation and association check
#Merged delete_unampped objects into delete_existing_xrefs due to dependency. Now throws if arrays param set until force_delete specified, deletes using IN list
#cache_arrays_per_probeset now uses arrays param to reduce size of cache
#removed unmapped_object dbID = 2000>
#Change all probe specific unmapped object to reeflect the individual probe rather than the probeset
#Updated logs
#Updated docs
#Added control of promiscuous probesets and unmapped objects
#Added array level stats
#added slice and transcript test modes
#Fix extension to only extend if there are no UTRs!!! Now much more flexible!
# Handle UTR/Extension overlap bug, wasn't catch feature across the UTR>extension border
# Remove Anti-sense UnmappedObjects?
# Test for external_dbs in healtcheck and force insert manually or prompt to change species name
#Implement Helper for logs etc.
# Handle non-probeset probes
# Add Probe level DBEntries to enable ProbeSet view? No just add ProbeFeature IDs to linkage annotation?
# Or just leave out, altho it will not be obvious exactly which ones will have been used in the mapping
# No we need to do both, so we can identify exactly which feature mapped
# But also know how many times a probe might have mapped to another transcript, to give a quality score!
# Pod::Usage for help
#Species specific check existing and delete based on external_db_name
#Change docs to Pod
#Can remove object_name and object_key from caches if we disable logs and just depend on xref/unmapped objects.
#Else we need to maintain them so we have names rather than dbIDs in the logs

#To do

# 1. Reimpliment validate arrays, see old script?
# 2. Add unannotated UTR clipping dependant on nearest neighbour
# 3. Extend UTRs to default length is they are less than defaults, so long as they don't overlap neighbour, 
#    then use annotated if present or clip to neighbour start/end if not, also accounting for default UTRs 
#    in the neighbour.
# 4. Separate UTR multipliers for 3' and 5'?
# 5. Implement incremental update from list of stable IDs. Consider unmapped probe changes etc. 
# 6. Parallelise by probeset chunks, can't do this by chromosome slices as we need to know genomewide 
#    counts for a given probeset. Calc UTRs then submit chunks jobs to farm
#    Chunk by retrieving all probesets and sorting an array of probeset names, then splice the array 
#    according to the number of chunks. We're still going to have retrieve all the transcripts and retrieve 
#    all probes for each, so we are really not gaining anything!! The only gain we can make is by chunking 
#    by slice, but then we need to know how many times something has mapped. Can we do some clean up afterwards? 
#    Let's add a clean up mode which simply deletes all probe sets which map too many times. We would need to 
#    ignore this threshold as we were mapping!!! So we don't delete and then mess up the counts for post run 
#    clean up.
# 7. There is no reason to have separate probe and xref DBs???
# 8. Validate array format against arrays specified? May want to just use an array format as a template???
# 9. Add mismatch filter for ProbeTranscriptAlign xrefs as match rules can differ between alignment and 
#    annotation
# 10.Handle ProbeAlign mismatch vs overlap mis match. Currently the overlap calculation is naive to the 
#    presence of alignment mis-matches.  Which means there is a possiblity of including probes with a total 
#    sequence mismatch of (align mismatch + overlap mismatch). This has always been the case.
# 11.Move ProbeAlign unmapped object storage to write_output, then this will not get written in test mode and 
#    we won't get duplication should the job fail halfway through
# 12.Enable probesets to have different sizes on different arrays, see notes in cache_arrays_per_object


#Ensembl Genomes stuff
# TEST Registry usage required as species will come from same DB
# In which case we need to take a species param for each of the transcript, array and xref DBs
# Also need to have -no_delete option, which will allow running of pipeline with previous xrefs stored?
# Or can we force delete to be species specific? We would need to do this anyway to support updating of species asynchronously
# We probably need to think about this for the 1st stage too, 
# but will be easy as we just need to dump the correct top level sequence
# Validate species against registry alias and use this to generate species_core_Gene DB rather than ensembl_core_Gene
# patch other efg DBs and alter External parsers accordingly.
# Can't rely on Registry as species aliases may not be present or loaded


# Issues
# Cannot account for running non-linked arrays which may use the same probe/set name.  This may cause failure if the probeset sizes are different. Xrefs and counts should be unaffected as we base these on the probe_set_ids not the names. This is not really an issue as unlinked arrays should not be run together
# Cannot currently handle probesets with different sizes between arrays, defaults to lowest probeset size to be permissive. See todo 12.


=head1 NAME

probe2transcript.pl

=head1 SYNOPSIS

This script performs probe(set) to transcript mapping based on a few simple parameters. Overlap analysis 
of ProbeFeatures is performed and annotations are stored as xrefs for individual ProbeFeatures, Probes or 
ProbeSets as a whole. Any probe(set)s which fail the mapping procedure are by default stored in the 
UnmappedObject tables and a logfile is also written.

e.g. perl probe2transcript.pl --species $SPECIES --transcript_dbname $DNADB_NAME --transcript_host $DNADB_HOST --transcript_port $DNADB_PORT --transcript_user $DNADB_USER --xref_host $DB_HOST --xref_dbname $DB_NAME --xref_user $DB_USER --xref_pass $DB_PASS --calculate_utrs --utr_multiplier 1 --arrays HT_MG-430A  -vendor AFFY -format AFFY_UTR


=head1 DESCRIPTION

More wordy description here

promiscuous probes
completed unmapped probes
linked arrays

This is generally executed by the eFG array mapping environment

=head1 OPTIONS

 Mandatory
  -species    Latin name as used in DB name or meta table e.g. homo_sapiens
  -vendor     Array vendor e.g. AFFY, ILLUMINA etc
  -format     Array format e.g. AFFY_UTR, AFFY_ST, ILLUMINA_WG. Sets default array configuration
  -arrays     List (space separated) of array names.

 Array configuration:
  -linked_arrays       Boolean(0|1) - For probe(set)s which exist on multiple array e.g. AFFY
  -probeset_arrays     Boolean(0|1) - For arrays which contain probesets rather than just single probes
  -sense_interrogation Boolean(0|1) - Sets interrogation strand, normally 1 for AFFY but 0 for AFFY_ST i.e. anti-sense

 Mapping rules:

  -mismatches                  Maximum number of mismatches allowed per probe
  -calculate_utrs              This calculates the default unannotated UTR extension defined by
                               the greater of either the mean or the median of all annotated UTRs
                               This is overridden by the following extend options
  -unannotated_5_utr           Default extension for transcripts with unannotated 5' UTRs
  -unannotated_3_utr           Default extension for transcripts with unannotated 3' UTRs
  -annotated_5_prime_extend    Default bp extension for all transcripts with annotated 5' UTRs
  -annotated_3_prime_extend    Default bp extension for all transcripts with annotated 3' UTRs
  -utr_multiplier              Defines UTR extension as multiple of annotated UTR. This is overridden 
                               if the above options are set.
  -max_transcripts             Maximum number of transcripts probe(set) can map to before we call it 
                               promiscuous, default is 100.
  -threshold                   This is the fraction(0-1) of probes within a probeset required to call it 
                               mapped, default is 0.5


 Running modes:

  -delete           Delete all pre-existing array xrefs generated by probe2transcript for given arrays
  -no_triage        Does not load UnmappedObjects (still writes to log file)
  -parallelise      Not yet implemented, will chunk and submit to farm
  -clean_up         Not yet implemented, will perform post parallelised run clean up
		 
 DB connection parameters, registry will override direct connection:

  -reg_verbose    Turns on verbose output when loading the registry

  -reg_host
  -reg_port
  -reg_user
  -reg_pass

  or

  -reg_file

  or

  -transcript_host          Mandatory
  -transcript_user          Mandatory
  -transcript_port    
  -transcript_pass          Mandatory
  -transcript_dbname        Mandatory
  -transcript_multi_species
  -transcript_species_id

  -xref_host          Mandatory
  -xref_user          Mandatory
  -xref_port
  -xref_pass          Mandatory
  -xref_dbname        Mandatory
  -xref_multi_species
  -xref_species_id

  -probe_host           Probe parameters
  -probe_user           default to xref
  -probe_port           DB paramters
  -probe_pass           if not
  -probe_dbname         specified
  -probe_multi_species
  -probe_species_id

 Testing:

  -test_transcripts   Number of transcripts to perform a test run on
  -slice              Name of test slice to perform a test run on 
  -transcript         Test transcript stable ID
  -no_store           Not yet implemented	
		 
 Other options:
  -tee                Tees output to STDOUT
  -filename           Sets name to be used in output and logfile, default is xref_dbname_probe2transcript.log|out
  -help               Prints this POD documentation and exits


=head1 EXAMPLE


=head1 SEE ALSO

ensembl-functgenomics/scripts/environments/arrays.env

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut



use strict;

use Pod::Usage;
use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (median mean get_date);


$| = 1; # auto flush stdout

my $reg = 'Bio::EnsEMBL::Registry';

#Helper params
#declared here to avoid single usage warning
$main::_log_file = undef;
$main::_tee      = 0;    #tee the output of this scripts


my ($transcript_host, $transcript_user, $transcript_pass, $transcript_dbname,
    $probe_host, $probe_user, $probe_pass, $probe_dbname, $load_from_db,
    $xref_host, $xref_user, $xref_pass, $xref_dbname, $calc_utrs,
    $test_transcripts, @arrays_names, $delete, $no_triage, $test_slice);
#$annotated_utrs);

my ($probe_db, $xref_db, $transcript_db, %promiscuous_objects, %transcripts_per_object, @unmapped_objects, $um_obj,
	%transcript_ids , %transcript_feature_info, %arrays_per_object, %probeset_sizes, @transcripts, %arrays,
   %array_xrefs, %transcript_xrefs, $test_transcript_sid, $clean_up, $parallelise, $filename, $sql, @array_names);

my $reg_verbose = 0;
my ($species, $reg_file, $reg_host, $reg_user, $reg_pass, $reg_port);

my ($transcript_multi_species, $xref_multi_species, $probe_multi_species) = (0,0,0);
my ($transcript_species_id, $xref_species_id, $probe_species_id) = (1,1,1);


# Default options
my $transcript_port = 3306; 
my $probe_port = 3306; 
my $xref_port = 3306;
my $max_mismatches = 1;
my $sense_interrogation;
#my $vendor = 'AFFY';
#my $format = 'AFFY_IVT';#Rename this UTR?
my ($vendor, $format, $multi_species);

my %array_config = (
					probeset_arrays      => undef,
					linked_arrays        => undef,
					sense_interrogation  => undef,
				   );

my %utr_extends = (
				   3 => undef,	  #three => 2000,
				   5 => undef,
				  );

#Set these to undef so we can allow no estimated UTR
#By setting these to 0 using the params
my %unannotated_utrs = (
						3 => undef,	  #three => 2000,
						5 => undef,
					   );


my ($utr_multiplier);

#my $unannotated_utr_length = 2000;
my $max_transcripts = 100;
my $mapping_threshold = 0.5;
#based on the current UTR length

my @tmp_args = @ARGV;

GetOptions(
		   'transcript_host=s'      => \$transcript_host,
           'transcript_user=s'      => \$transcript_user,
           'transcript_port=i'      => \$transcript_port,
           'transcript_pass=s'      => \$transcript_pass,
           'transcript_dbname=s'    => \$transcript_dbname,
           'transcript_species_id=i' => \$transcript_species_id,
           'transcript_multi_species' => \$transcript_multi_species,
		   'probe_host=s'           => \$probe_host,
           'probe_user=s'           => \$probe_user,
           'probe_port=i'           => \$probe_port,
           'probe_pass=s'           => \$probe_pass,
           'probe_dbname=s'         => \$probe_dbname,
           'probe_species_id=i' => \$probe_species_id,
           'probe_multi_species' => \$probe_multi_species,
		   'xref_host=s'            => \$xref_host,
           'xref_user=s'            => \$xref_user,
           'xref_port=i'            => \$xref_port,
           'xref_pass=s'            => \$xref_pass,
           'xref_dbname=s'          => \$xref_dbname,
           'xref_species_id=i' => \$xref_species_id,
           'xref_multi_species' => \$xref_multi_species,
		   'reg_file=s'             => \$reg_file,
		   'reg_host=s'             => \$reg_host,
		   'reg_user=s'             => \$reg_user,
		   'reg_pass=s'             => \$reg_pass,
		   'reg_port=i'             => \$reg_port,
		   'reg_verbose'            => \$reg_verbose,
		   'species=s'              => \$species,
		   'vendor=s'               => \$vendor,
		   'format=s'               => \$format,
		   'mismatches=i'           => \$max_mismatches,
           'annotated_3_prime_extend=s'       => \$utr_extends{3},
		   'annotated_5_prime_extend=s'       => \$utr_extends{5},
		   'calculate_utrs'         => \$calc_utrs,
		   'unannotated_5_utr=s'    => \$unannotated_utrs{5},
		   'unannotated_3_utr=s'    => \$unannotated_utrs{3},
		   'utr_multiplier=s'       => \$utr_multiplier,#Make this for 5 and 3?
		   'test_transcripts=i'     => \$test_transcripts,
		   'max_transcripts=i'      => \$max_transcripts,
		   'threshold=s'            => \$mapping_threshold,
		   'arrays=s{,}'            => \@array_names, # this should take 1 or more space separate array names WARNING experimental feature!
	
		   'delete'                 => \$delete,
		   #'force_delete'           => \$force_delete,
		   'no_triage'              => \$no_triage,
		   #'health_check'           => \$health_check,
		   'parallelise'            => \$parallelise,
		   'clean_up'               => \$clean_up,
		   'slice=s'                => \$test_slice,#Only for testing purposes!
		   'transcript=s'           => \$test_transcript_sid,

		   'linked_arrays=i'          => \$array_config{linked_arrays},
		   'probeset_arrays=i'        => \$array_config{probeset_arrays},
		   'sense_interrogation=i'    => \$array_config{sense_interrogation},

		   #Helper params
		   'tee'                    => \$main::_tee,
		   'filename'               => #\$main::_log_file,
		   
		   #add a reduced log to minimize memory usage?
           'help'                   => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);





#exit if unkown options specified? TEST!

#change this to just @ARGV?
#@arrays = split(/,/,join(',',@arrays));#?

# Have default setting flags handled here which will not override other param
# in default 'set'



#Make arrays mandatory?
#Not sensible to AFFY AFFY_ST and ILLUMINA at same time!

#Set log type so we are no over writing to the same files for different 
#format, or custom formats
my $log_type = $format || $$;
$filename ||= "${xref_dbname}_${log_type}_probe2transcript";
$main::_log_file ||=  "./${filename}.log";
my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;
my $hostname = `hostname`;
chomp($hostname);
$Helper->log_header('Running on probe2transcript.pl on: '.$hostname, 0, 'append_date');
$Helper->log("Params are:\t@tmp_args");


die("It is not wise to run all available arrays at the same time\nYou must supply a list of array names using -arrays, i.e. for all or a subset of a given array format(e.g. AFFY_UTR, AFFY_ST, ILLUMINA_WG)") if(! @array_names);


#OTHER MANDATORY PARAMS HERE?

#Now set some array config

my %array_format_config = (
						   AFFY_UTR => {
									   probeset_arrays         => 1,
									   linked_arrays      => 0,
									   sense_interrogation  => 0,
									   
									  },
						   AFFY_ST => {
									   probeset_arrays         => 1,
									   linked_arrays      => 1,
									   sense_interrogation  => 1,
									  },
						   
						   ILLUMINA_WG => {
										probeset_arrays        => 0,
										linked_arrays     => 1,
										sense_interrogation => 0,
										  },
						   
						   AGILENT => {
										   probeset_arrays        => 0,
										   linked_arrays     => 1,
										   sense_interrogation => 0,			   
										  },					


						     PHALANX => {
										 probeset_arrays        => 0,
										 linked_arrays     => 1,
										 sense_interrogation => 0,
										 },

						   CODELINK => {
										probeset_arrays        => 0,
										linked_arrays     => 1,
										sense_interrogation => 0,
									   },

						   LEIDEN => {
									  probeset_arrays        => 0,
									  linked_arrays     => 1,
									  sense_interrogation => 0,
									 },


						  );
die ('Must supply a -vendor parameter e.g. AFFY') if ! $vendor;

if(defined $format && ! exists $array_format_config{$format}){
  die("-format is not valid:\t$format\nMust specify valid format e.g. AFFY_UTR, AFFY_ST, ILLUMINA_WG.\nOr maybe you want to use -probeset_arrays, -linked_arrays and -sense_interrogation to define the format parameterss?\n");
}


if(! ($array_config{probeset_arrays} && 
	  $array_config{linked_arrays} && 
	  $array_config{sense_interrogation})
   && ! $format){
  die('You must specify a valid format parameter if you are not using -probeset_arrays, -linked_arrays and -sense_interrogation\n');
}


foreach my $key(keys %array_config){

  if(! defined $array_config{$key}){
  
	if(exists $array_format_config{$format}){
	  $array_config{$key} = $array_format_config{$format}{$key};
	}
	else{
	  #This should never happen
	  die("Cannot find default $key config for $format format");
	}
  }
}

my $xref_object = ($array_config{probeset_arrays}) ? 'ProbeSet' : 'Probe';


#we need to do a check here on utr_length and unannotated_utr_length

#Let's just have one species!?
#What if we want to use different DBs for the xref, and probe DB?
#Then we'll just have to use the old method and specify different dbnames

if(! $species){
  die("Must provide a -species");
}


if($reg_file || $reg_host){

  #Some sanity checks
  if($reg_host && $reg_file){
	die('You have specified confliciting parameters -reg_file -reg_host');
  }

  if($probe_dbname || $xref_dbname || $transcript_dbname){
	die('You have specified conflicting paramters -reg_dbname and -probe_dbname or -transcript_dbname or -xref_dbname');
  }
  
  #load the reg
  if($reg_file){
	
	$Helper->log("Loading registry from:\t".$reg_file);
	$reg->load_all($reg_file, $reg_verbose);
  }
  else{#reg_host

	if(! ($reg_user && $reg_pass)){
	  die('Must provide at least a -reg_user -reg_pass (optional -reg_port -reg_verbose) if loading from db');
	}
	
	$reg->load_registry_from_db(
								-host    => $reg_host,
								-port    => $reg_port || 3306,
								-user    => $reg_user,
								-pass    => $reg_pass,
								-verbose => $reg_verbose,
							   );
  }	

  #Now load the db adaptors

  $transcript_db = $reg->get_DBAdaptor($species, 'Core');
  $xref_db       = $reg->get_DBAdaptor($species, 'Funcgen');
  #$probe_db      = $reg->get_DBAdaptor($probe_species, 'Funcgen');
  $probe_db = $xref_db;
}
else{#load dbs from params

  if(!$transcript_user || !$transcript_dbname || !$transcript_host){
	warn "You must specify a -transcript_user -transcript_dbname -transcript_host\n";
	usage();
  }

  $transcript_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
		-host    => $transcript_host,
		-port    => $transcript_port,
		-user    => $transcript_user,
		-pass    => $transcript_pass,
		-dbname  => $transcript_dbname,
		-species => $species,
		-multispecies_db => $transcript_multi_species,
		-species_id => $transcript_species_id
	);


  #Change this to only check for dbname?
  #quite likely core and efg DB will be on different machines

  if(!$xref_user || !$xref_dbname || !$xref_host){
	warn "You must specify a -xref_user -xref_dbname and -xref_host\n";
	usage();
  }

  $xref_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
		-host => $xref_host,
		-port   => $xref_port,
		-user   => $xref_user,
		-pass   => $xref_pass,
		-dbname => $xref_dbname,
		-species => $species,
		-multispecies_db => $xref_multi_species,
		-species_id => $xref_species_id
	);

  if ($probe_host && $probe_dbname && $probe_user) {
	
	$probe_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
		-host   => $probe_host,
		-port   => $probe_port,
		-user   => $probe_user,
		-pass   => $probe_pass,
		-dbname => $probe_dbname,
		-species => $species,
		-multispecies_db => $probe_multi_species,
		-species_id => $probe_species_id
	);
	
  } else{
		$Helper->log("No probe DB params specified, defaulting to xref params");
		$probe_db = $xref_db;
  }
}





#Test the DBs here before starting
$transcript_db->dbc->db_handle;
#print $transcript_db->species."\n";
$xref_db->dbc->db_handle;
#print $xref_db->species."\n";
$probe_db->dbc->db_handle;
#print $probe_db->species."\n";
#Registry was doing the species increment trick.



#Grab species ID for healtcheck delete and check
my $species_id = 1;

if($xref_db->is_multispecies){
  #Can probably do this via the MetaContainer?
  my $sql     = "SELECT species_id from meta where meta_key='species.db_name' and meta_value='$species'";
  ($species_id) = $xref_db->dbc->db_handle->selectrow_array($sql);

  if(! $species_id){
		#This should never occur as we have generated the xref_db with this species
		die("Could not find species_id for $species");
  }
}


#Check for external_db records for species DBs

my $schema_build = $xref_db->_get_schema_build($transcript_db);
#Should we allow a param to over ride this?
#This should be used in all the DBEntry and UnmappedObject records
my ($edb_name, $transc_edb_name, $transc_edb_id, $transc_edb_display_name, $edb_display);

#for my $edb_type('Transcript', 'Species'){
#  my $found_edb_id = 0;

#  if($edb_type eq 'Transcript'){
$edb_name                = "${species}_core_Transcript";
$transc_edb_name         = $edb_name;
$transc_edb_display_name = "EnsemblTranscript";
$edb_display             = $transc_edb_display_name;
	
#  }else{
#	#This is used for storing completely Unmapped probes
#	$edb_name = 'ensembl_core_Species';
#	$edb_display = 'Ensembl Species';
#  }


$sql = "SELECT external_db_id, db_release from external_db where db_name='$edb_name'";
#warn "$sql";

my @versions = @{$xref_db->dbc->db_handle->selectall_arrayref($sql)};
$sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('.
  "'${edb_name}', '${schema_build}', 'KNOWNXREF', 1, 5, '$edb_display', 'MISC')";
my @tmp;
  
foreach my $row(@versions){
  my ($edb_id, $version) = @$row;
  push @tmp, $version;
  
  if($schema_build eq $version){ 
	$transc_edb_id  = $edb_id;
	last;
  }
  #$transc_edb_id  = $edb_id if($edb_type eq 'Transcript');
  #$species_edb_id = $edb_if if($edb_type eq 'Species');
}


if(! $transc_edb_id){
  $sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('.
	"'${edb_name}', '${schema_build}', 'KNOWNXREF', 1, 5, '$edb_display', 'MISC')";
  die("Could not find current external_db $edb_name $schema_build from available versions:\t @tmp\nMaybe you have mis-spelt the -trans-species or you may need to manually add the external_db to the table and master file:\n\n$sql\n\n");
}




#Can we make this account for ArrayChips
#So we can do a staged run and not have to delete all previous xrefs
#if we subsequently get some for a previously unavailable array_chip

#This validates arrays
#Why would we ever want to write the xrefs to a different DB?
#my %array_name_cache =  %{&validate_arrays($probe_db, $xref_db)};

#Validate array names
my $array_adaptor = $xref_db->get_ArrayAdaptor;


my $array_format;

foreach my $name(@array_names){

  my $array = $array_adaptor->fetch_by_name_vendor($name, $vendor);

  if(! $array){
	die("Could not find $vendor $name array in DB");
  }

  #Now test array is of same format
  $array_format ||= $array->format();

  if($array_format ne $array_format){
	die('You must not map arrays of different formats in the same process');
	#This is really to keep the run time/mem usage down for a given process
  }



  $arrays{$name} = $array;
}



#$delete = $force_delete if $force_delete;

#Merge these as they are related and we need to force unmapped check first

if($delete){
  delete_existing_xrefs($xref_db);
}
else{
  check_existing_and_exit($probe_db, $xref_db);
}


#if ($health_check){
#  $Helper->log('Have you migrated your probe_features from a different DB?\n'.
#			   "If yes, check meta_coord entries and of.analysis_id=a.analysis_id\n".
#			   "Current probe_feature analyses are:\nanalysis_id\tlogic_name");
  
  #This can happen if you have simply migrated the tables from another DB
  #also use Healtchecker for meta_coord update.
  #Need to re-write the Migrate function in arrays.env
  
  
$sql = 'SELECT distinct pf.analysis_id, a.logic_name from probe_feature pf left join analysis a on pf.analysis_id=a.analysis_id';
my @analysis_info = @{$probe_db->dbc->db_handle->selectall_arrayref($sql)};

foreach my $record(@analysis_info){
  my ($anal_id, $lname) = @$record;
  #$Helper->log("$anal_id\t\t$lname");
  
  if(! $lname){
	die("Found probe_feature analysis without a corresponding analysis entry:\t$lname");
  }
	
#  if(! ($lname =~ /_ProbeAlign/ || $lname =~ /_ProbeTranscriptAlign/)){
#	die("Found unexpected/mismatched analysis entry in probe_feature table:\t$lname");
  #}
}


my $transcript_adaptor = $transcript_db->get_TranscriptAdaptor();
my $slice_adaptor = $transcript_db->get_SliceAdaptor();
my $probe_feature_adaptor = $probe_db->get_ProbeFeatureAdaptor();
my $dbentry_adaptor = $xref_db->get_DBEntryAdaptor();
my $analysis_adaptor = $xref_db->get_AnalysisAdaptor();
my $unmapped_object_adaptor = $xref_db->get_UnmappedObjectAdaptor();

my $analysis = get_or_create_analysis($analysis_adaptor);

my $i = 0;
my $last_pc = -1;

my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new();

open (OUT, ">${filename}.out");

#no defined extends and no multiplier would simply be just the annotated utr, which we don't want
#no defined extends and mulplier would be multiplier for 5 and 3 prime
#one extend defined (can be 0) and multiplier would be hard extend(or no extend/annotated) for one and multiplier for other
#We currently can't set different multipliers for 3 and 5 prime.

if(defined $utr_extends{3} && defined $utr_extends{5} && $utr_multiplier){#Can't have all
  die('You cannot set both -3/5_prime_extend values and a -utr_multiplier');
}
elsif(!(defined  $utr_extends{3} || defined $utr_extends{5} || defined $utr_multiplier)){#Can't have none
  die("You must set some extension rules using  -3/5_prime_extend values and/or -utr_multiplier\n".
	  "Set -utr_multiplier to 0 and omit -3/5_prime_extend to run against UTR only\n");
}
else{
	$Helper->log("You have specified -utr_multiplier and a -3|5_prime_extend, -3|5_prime_extend will override where appropriate");
}

#Multiplier needs to be a pisitive real number
if(defined $utr_multiplier && $utr_multiplier !~ /-*[0-9]+[\.0-9]*/){
  die("-utr_multiplier must be a positive real number:\t$utr_multiplier");
}

#Validate extend params
for my $end('3', '5'){
  
  if(defined $utr_extends{$end} && $utr_extends{$end} =~ /\D+/){
	die("Invalid -${end}_prime_extend parameter(".$utr_extends{$end}.").  Must be a number(bp)");
  }
}

#Can't use both default unannotated UTRs and calculate
if($unannotated_utrs{5} && $unannotated_utrs{3} && $calc_utrs){
  die('You cannot set both -unannotated_5/3utr values and -calc_utrs');
}



#Validate unannotated defaults
for my $end('3', '5'){

  if(defined $unannotated_utrs{$end}){
	
	if($unannotated_utrs{$end} =~ /^\D+$/){
	  die("Invalid -unannotated_${end}_utr parameter(".$unannotated_utrs{$end}.").  Must be a number(bp)");
	
	}
	else{
	  $Helper->log("Setting $end unannotated UTR length to ".$unannotated_utrs{$end});
	}
  }
  else{
	if(! $calc_utrs){
	  $Helper->log("Defaulting unannotated $end UTR length to 0");
	  $unannotated_utrs{$end} = 0;
	}
  }
}


if($test_slice && $test_transcript_sid){
  die('Can only run in one test mode, please specify -slice or -transcript');
}
  

#Need to warn here if max_transcripts??

if($test_slice){
  $Helper->log("Running in test mode with slice:\t$test_slice\n".
			   "WARNING:\tPromiscuous probesets will not be caught! Calculated UTRs will be wrong!");


  #Need to add better text here when we have implemented parallel runs
  my $slice = $slice_adaptor->fetch_by_name($test_slice);
  die("Could not get slice from the DB:\t$slice") if ! defined $slice;
  @transcripts = @{$transcript_adaptor->fetch_all_by_Slice($slice)};
}
elsif($test_transcript_sid){
  $Helper->log("Running test mode with transcript:\t$test_transcript_sid\n".
			   "WARNING:\tPromiscuous probeset will not be caught!");

  @transcripts = ($transcript_adaptor->fetch_by_stable_id($test_transcript_sid));
}
else{
  @transcripts = @{$transcript_adaptor->fetch_all()};
}

my $total =  scalar(@transcripts);

die('Could not find any transcripts') if $total == 0;
  
$Helper->log("Identified ".scalar(@transcripts)." transcripts for probe mapping");
$Helper->log("Allowed mismatches = $max_mismatches");

if($calc_utrs){
  $Helper->log('Calculating default UTR lengths from greatest of max median|mean', 0, 'append_date');
  #We actually want to extend the potentially conservative ensembl UTRs
  #to the calculated default if they are shorter, but only if this does 
  #not cause overlap with a neighbouring gene.
  #Do not implement UTR extension until clipping is in place
  #This will require knowledge of genomic context
  #What is fastest solution here?
  #1. Run in a slice context to fetch all genes, then we know the previous transcript
  #and can easily access the next transcript
  #2. Simply generate an extended slice from the transcript and pull back genes
  #We would have to either do this for every trans or just for the longest
  #
  #1 is probably most efficient altho' and will actually reduce the memory usage by
  #chunking by chromosome
  my ($five_utr, $three_utr, @five_lengths, @three_lengths);
  my ($three_median, $five_median, $remainder, $three_mean, $five_mean);
  my $five_cnt  = 0;
  my $three_cnt = 0;
  my $five_zero_cnt = 0;
  my $three_zero_cnt = 0;
 
  foreach my $transcript(@transcripts){

	 
	if(! defined $unannotated_utrs{3}){
	  $three_utr = $transcript->three_prime_utr;

	  if(defined $three_utr){
		$three_cnt++;
		push @three_lengths, $three_utr->length;
	  }
	  else{
		$three_zero_cnt++;
	  }
	}


	if(! defined $unannotated_utrs{5}){
	  $five_utr  = $transcript->five_prime_utr;

	  if(defined $five_utr){
		$five_cnt++;
		push @five_lengths, $five_utr->length;
	  }
	  else{
		$five_zero_cnt++;	
	  }
	}
  }

 
  if(! defined $unannotated_utrs{5}){
	$Helper->log("Seen $five_cnt 5' UTRs, $five_zero_cnt have length 0");

	if($five_cnt){
	  $five_mean = mean(\@five_lengths);
	  ($five_mean, $remainder) = split/\./, $five_mean;
	  $five_mean++ if $remainder =~ /^[5-9]/;
	  @five_lengths = sort {$a <=> $b} @five_lengths;
	  $five_median = median(\@five_lengths);
	  $unannotated_utrs{5}  = ($five_mean  > $five_median)  ? $five_mean  : $five_median;
	  $Helper->log("Calculated default unannotated 5' UTR length:\t".$unannotated_utrs{5});
	}
	else{
	  die("Found no 5' UTRs, you must specify a -unannotated_5_utr");
	}
  }

   if(! defined $unannotated_utrs{3}){
	 $Helper->log("Seen $three_cnt 3' UTRs, $three_zero_cnt have length 0");

	 if($three_cnt){

	   $three_mean = mean(\@three_lengths);
	   ($three_mean, $remainder) = split/\./, $three_mean;
	   $three_mean++ if $remainder =~ /^[5-9]/;
	   @three_lengths = sort {$a <=> $b} @three_lengths;
	   $three_median = median(\@three_lengths);
	   $unannotated_utrs{3} = ($three_mean > $three_median) ? $three_mean : $three_median;
	   $Helper->log("Calculated default unannotated 3' UTR length:\t".$unannotated_utrs{3});
	 }
	 else{
	   die("Found no 3' UTRs, you must specify a -unannotated_3_utr");
	 }
   }

  $Helper->log('Finished calculating unannotated UTR lengths', 0, 'append_date');
}



if($parallelise){

  #Need to test for $test_slice
  #Don't count max mapping per transcript if test slice

  die('Not yet implemented, need to tidy up unmapped object etc');

  #we need to test for which probe2transcript_new.pl

  my @slices = $slice_adaptor->fetch_all('toplevel', undef, 1);
  
  foreach my $slice(@slices){
	#my $bsub_cmd = "bsub -J $xref_dbname.".$slice->seq_region_name.".probe2transcript -q normal -R "select[mem>8000] rusage[8000=mem]" -M 8000000 -o ${DB_HOME}/mapper.out -e ${DB_HOME}/mapper.err time perl ${SRC}/ensembl/misc-scripts/probe_mapping/probe2transcript.pl --transcript_dbname $DBNAME --transcript_host $DBHOST --transcript_port $DBPORT --transcript_user ensro --xref_host $ODBHOST --xref_dbname $ODBNAME --xref_user ensadmin --xref_pass ensembl --utr_length annotated -calculate_utrs";
	
  }
}



my %flanks = (
			  5 => 0,
			  3 => 0,
			 );

### Not getting any xrefs? Try checking the meta_coord table if you have migrated the probe_features
#Need to healtcheck this!!

#Do caching first to highlight any problem with arrays before running mapping
cache_arrays_per_object($probe_db);

my ($linkage_annotation, $pc, $transcript_slice, $slice);
my (@exons, $num_exons, $first_exon, $last_exon, $probe_features, $probe, $probe_id, $feature_id);
my ($probeset_id, $probeset_name, $transcript_key, $log_name, $transcript_sid);
my $failed_extend = 0;

$Helper->log("Performing overlap analysis. % Complete:");

foreach my $transcript (@transcripts) {

  $pc = int ((100 * $i) / $total);
  $i++;

  if ($pc > $last_pc) {
	$Helper->log("$pc ", 0, 0, 'no_return');
    $last_pc = $pc;
  }

  last if ($test_transcripts && $i >= $test_transcripts);
  $transcript_sid = $transcript->stable_id();

   
  #Handle UTR extensions
  #The UTRs themselves are already included in the transcript/exons!!
  foreach my $end('5', '3'){	
	#so we only need to consider if we have extend, multplier or a default unannotated
	#These could all be set to 0

	if($utr_extends{$end} || $utr_multiplier || $unannotated_utrs{$end}){
	  my $method = ($end == 5) ? 'five' : 'three';
	  $method .= '_prime_utr';
	  my $utr = $transcript->$method;

	  if(defined $utr){
		
		#Hard extend takes presidence over multiplier
		if(defined $utr_extends{$end}){
		  $flanks{$end} = $utr_extends{$end};
		}
		else{#mulitplier
		  $flanks{$end} = $utr->length * $utr_multiplier;
		}
	  }
	  elsif(defined  $unannotated_utrs{$end}){#Use unannotated default
		#which can either bet set(can be 0) or calculated
		$flanks{$end} = $unannotated_utrs{$end};
	  }
	}
  }
  
  $transcript_slice = $transcript->feature_Slice();
  $slice            = $transcript_slice->expand($flanks{5}, $flanks{3}) if $flanks{5} || $flanks{3};

  #TODO as $slice is undefined if the above statement does not run!
  #$slice           ||= $transcript_slice;

  if(! $slice){
	$failed_extend++;
	$slice = $transcript_slice;
  }


  #We currently don't account for extension of next transcript!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #We need to do this in a slice context, looking ahead at the next transcript to see
  #if there would be any over lap, then removing any extensions
  #And also take into account strandedness of neighbouring transcripts
  #So we we need to process the next transcript on each strand from a different gene
  #How will miRNA genes effect this?
  #foreach my $end('five', 'three'){
  #	if($flanks{$end}){
  
  #  }

  #Register exon and utr loci
  $rr->flush();
  @exons     = @{$transcript->get_all_Exons};
  $num_exons = $#exons;
  #These exons start/ends are always local
  #But we pass probe_feature seq_region vals to the rr?!

  foreach my $e (0..$#exons) {
	$rr->check_and_register('exonic', $exons[$e]->seq_region_start, $exons[$e]->seq_region_end);
	
	$first_exon = $exons[$e] if $e == 0;
	$last_exon  = $exons[$e] if $e == $num_exons;
  }

  #These should include the first and last exon also, otherwise we may omit features which lie on the
  #boundary between the EXON/UTR and the extension
 
 
  #Sub slice from the full slice so we don't have to worry about CoordSystem name
  #Just omit either three or five depending on flanks
  my %exonutrs;#Clean

  if ($transcript->strand == 1) {
	$exonutrs{3} = [$last_exon->seq_region_start, $slice->end]   if $flanks{3};
	$exonutrs{5}  = [$slice->start, $first_exon->seq_region_end] if $flanks{5};	
  } 
  else {#-1 reverse strand
	$exonutrs{3} = [$slice->start(), $first_exon->seq_region_end] if $flanks{3};
	$exonutrs{5}  = [$last_exon->seq_region_start, $slice->end]   if $flanks{5};
  }

  #Test whether this overlaps with a neighbouring transcript on the same strand
  #for each of the exonutr slices
  #And set a flag to so we can record this in the Unmapped object
  foreach my $end('3', '5'){
	$rr->check_and_register("${end}_exonutr", @{$exonutrs{$end}}) if $flanks{$end};
  }

  ### Now map each feature to the transcript
  #This works on the assumption that probesets are identical between arrays
  #i.e. if a probe set is present on different arrays, their probes are identical.
  $probe_features = $probe_feature_adaptor->fetch_all_by_Slice_Arrays($slice, [values(%arrays)]);
  my %probe_feature_xrefs;

  foreach my $feature (@$probe_features) {
 	$probe          = $feature->probe;
	$feature_id     = $feature->dbID;
	$probe_id       = $probe->dbID;

	#Set some probe/probeset vars
	if($array_config{probeset_arrays}){
	  $probeset_id    = $probe->probeset->dbID;
	  $probeset_name  = $probe->probeset->name;
	  #Due to resitrcting the collapse within probesets
	  #we can only ever have one probe set name here.
	  $log_name       = $transcript_sid."\t(${probeset_name})\t${probe_id}";
	  $transcript_key = $transcript_sid.":".$probeset_id;
	}
	else{
	  #clean just in case
	  #It would also be nice to print the probe name here
	  #as there is generally only one
	  #but we can't guarantee that this won't change.
	  #Leaving this out makes the logs next to useless
	  #For none probeset probes
	  #Let's just assume that this will never be AFFY
	  #It will fail if it is.

	  #We may have probe present on 2 arrays, so we must pass array name here?

	  #HACK
	  #Let's just get all probenames and make sure they are all the same
	  #Otherwise we need to either collapse separately
	  #Or account for probes with different names
	  #my %pnames;
	  #map $pnames{$_} = undef, @{$probe->get_all_probenames};

	#  my $pname;

#	  if(keys %pnames == 1){
#		#die("Found inconsistent probe names between arrays:\t".$probe->dbID.' has names '.join(', ', keys %pnames).'.');
#		($pname) = keys %pnames;
#	  }
	##  else{

	#	#This will not match arrays_per_object
	#	#And so will not get logged if orphan!
	#	$pname = 'MULTI_NAME'
	#  }

	  $probeset_id    = '';
	  $probeset_name  = '';
	  #Probe arrays are collapsed slightly differently as there is no probeset to restrict
	  #the collapse within, any probe across all the available arrays may be associated
	  #i.e. we may have more than one probe name for the same sequence
	  $log_name       = $transcript_sid."\t(".join(',', @{$arrays_per_object{$probe_id}{names}}).")\t${probe_id}";
	  $transcript_key = $transcript_sid.':'.$probe_id;	
	}


	#Do we need these as UnmappedObjects or in the log at all?
	if($array_config{sense_interrogation}){
	  
	  if($transcript->seq_region_strand == $feature->seq_region_strand){
		print OUT "Unmapped sense ".$log_name."\n";
		next;
	  }
	}
	elsif ($transcript->seq_region_strand != $feature->seq_region_strand){
	  print OUT "Unmapped anti-sense ".$log_name."\n";
	  next;
	}

	my $cigar_line = $feature->cigar_line;

	if($cigar_line =~ /D/){#ProbeTranscriptAlign
	  #Do we skip this and just get all in one go by the external_id
	  #e.g. the transcript ID?
	  
	  #my @xrefs = @{$dbentry_adaptor->fetch_all_by_Transcript($transcript)};
	  #This will not return the ProbeFeature objects?
	  #But we only need the probe ID to build the cache
	  #so

	  my @dbentries = @{$feature->get_all_Transcript_DBEntries};


	  if(@dbentries){
		my $txref;
		

		foreach my $dbentry(@dbentries){

		  if($dbentry->primary_id eq $transcript_sid){
			$txref = 1;
			
			#We need to implement mismatch checking here as we can define different rules for alignment
			#and annotation

			#Add the probe count
			if ($array_config{probeset_arrays}) {
			  
			  if (! $transcript_feature_info{$transcript_key}{$probe_id}) {
				$transcript_feature_info{$transcript_key}{$probe_id} = 1;
			  } else {
				$transcript_feature_info{$transcript_key}{$probe_id}++;
			  }
			} 
			else {
			  $transcript_feature_info{$transcript_key}{$probe_id} ||= ();
			  push @{$transcript_feature_info{$transcript_key}{$probe_id}}, 'exon-exon match';
			}
			
			#No need to add an xref as we already have one
			#last;#likely only ever to be one DBEntry/ProbeFeature
		  }
		}
		
		if(! $txref){
		  print OUT "Unmapped Gapped ProbeFeature ".$log_name."\n";
		  
		  if (!$no_triage) {		
			$um_obj = new Bio::EnsEMBL::UnmappedObject(
													   -type       => 'probe2transcript',
													   -analysis   => $analysis,
													   -identifier => $transcript_sid,
													   -external_db_id => $transc_edb_id,
													   -summary    => "Unmapped Gapped ProbeFeature",
													   -full_desc  => "Gapped ProbeFeature did not match transcript structure",
													   -ensembl_object_type => 'ProbeFeature',
													   -ensembl_id => $feature_id,
													  );
			&cache_and_load_unmapped_objects($um_obj);
			
		  }
		}
	  }
	}
	else{#ProbeAlign

	  #Handle flank mismatches
	  my $five_mismatch  = 0;
	  my $three_mismatch = 0;
	  my $feature_start  = $feature->seq_region_start;
	  my $feature_end    = $feature->seq_region_end;
	  
	  if($cigar_line =~ /(^[0-9]+m)/){
		($five_mismatch = $1) =~ s/m//;
		$feature_start += $five_mismatch;
	  }
	  
	  if($cigar_line =~ /([0-9]+m$)/){
		($three_mismatch = $1) =~ s/m//;
		$feature_end -= $three_mismatch;
	  }
	
	  #Now takes into account overhangs/end mis matches and internal mismatches??
	  #This is naive to the internal mismatch as it could be anywhere in the alignment
	  #Leave like this for now as it is very unlikely to happen and is currently on the conservative side
	  #Likely only to happen if the max mismatches are set high and there is a mis match close to but not at an end
	  my $min_overlap  = ($probe->length - $max_mismatches + ($feature->mismatchcount - $five_mismatch - $three_mismatch));
	  my $exon_overlap = $rr->overlap_size('exonic', $feature_start, $feature_end);

	  #This could map like this
	  #exonutrextension
	  #---probe--------
	  #Where is maps to last exon and to the exonutrextension
	  #If we have no extension and it was just 1 bp overlap then we have to take acount of the potential 1bp mismatch position
	  #As could cause it to be mapped or not
	  #So we are accounting for end mismatches already, but not internal matches
	

	  ### Pre compute some values for the flanks
	  #Can only be one of 4 values
	  #0          = Not flank
	  #1          = Flank used but no overlap
	  #three|five = Flank used and overlap present
	  my $flank_end     = 0;
	  my $flank_overlap = 0;

	  #This assumes we only have one flank overlap!!!!!!!!!!!!!!!!!
	  #Fine...unless we had a probe the full length of the transcript!!!!!

	  foreach my $end ('3', '5') {
		$flank_overlap = $rr->overlap_size("${end}_exonutr", $feature_start, $feature_end) if $flanks{$end};

		if ($flank_overlap) {
		  $flank_end = $end;
		  last;
		}
	  }
	

	  ### Test overlaps
	  if (($exon_overlap >= $min_overlap) ||
		  ($flank_overlap >= $min_overlap)) {

		#We don't need the feature dbIDs here(to add to the ProbeSet/Probe linkage_annotation
		#if we are storing the ProbeFeature xrefs
		#But we do want to pass the linkage annotation for single Probes
		#It is possible that we may get two matches of a single probe to a transcript
		#We should then change the linkage annotation to double/triple hit?
		#$transcript_feature_count{$transcript_key}{$dbID} ||= ();
		#push @{$transcript_feature_count{$transcript_key}{$dbID}}, $feature->dbID;
		#we inc here but we never use the count, just record the fact that is has mapped
		#Can we add this to the mapping annotation?
		#Do we need to store DBEntries for individual probes?
		#If we want to display them in ProbeSet panel yes.
		#As we currently only store Unmapped info for Probes
		#And ProbeSet for mapped info.

	 
		if ($exon_overlap && $flank_overlap) {
		  $linkage_annotation = "exon/${flank_end}' flank boundary";
		} elsif ($exon_overlap) {
		  $linkage_annotation = 'exon match';
		} else {				#only flank over lap
		  $linkage_annotation = "${flank_end}' flank";
		}

		#Count the probe & add the ProbeFeature xref
		if ($array_config{probeset_arrays}) {
		
		  if (! $transcript_feature_info{$transcript_key}{$probe_id}) {
			$transcript_feature_info{$transcript_key}{$probe_id} = 1;
		  } else {
			$transcript_feature_info{$transcript_key}{$probe_id}++;
		  }
		} 
		else {
		  $transcript_feature_info{$transcript_key}{$probe_id} ||= ();
		  push @{$transcript_feature_info{$transcript_key}{$probe_id}}, $linkage_annotation;
		}

		add_xref($transcript_sid, $feature_id, 'ProbeFeature', $linkage_annotation);
	  } 
	  else {					# must be intronic, intron-exon, intron-first|lastexon, five_prime_flank, three_prime_flank
	  	  
		#Ignore the first last exon classes and count these as intron-exon
		#We only know if is intron-exon if we have a negative result for the flank extension
		#So we could call this exon boundary instead?
		#We only lose this if we don't extend at all.

		#Deal with intronic first
		my ($summary, $region);
	  
		if (!( $exon_overlap || $flank_overlap)) {
		  $summary = 'intronic';
		  $region  = 'intronic region';
		} elsif ($exon_overlap) {

		  if (! $flank_overlap) {
			#This can be no overlap or no flank
			#e.g. exon-intron boundary or 5|3' exon boundary
			#No way of telling so just call both of these exon boundary
			$summary = 'exon boundary';
			$region  = $summary;
		  } else {
			#This has to be a flank boundary
			$summary = "${flank_end}' flank boundary";
			$region  = $summary;
		  }
		} else {				#only flank over lap
		  $summary = "${flank_end}' flank boundary";
		  $region  = $summary;
		}


		#No this could also over lap the 5' or 3' ends and/or intron depending on exon size?
		#What about exon/utr overlap, this would be a valid hit but would not be caught by this overlap method?!
		#Unless overlap method takes into account adjacent exons for UTRs
		#UTRs included in exons?! So not point in doing UTR overlap 
		#Unless we're extending, then we need to be mindful probes which may overlap short UTR and last/first coding exon
		#e.g.
		#last exon
		#      UTR
		#         extension
		#     probe

		#The probe in the above alignment would not be caught!
		#We need to treat the 5' and 3' + extensions separately and redundantly
		#So we know where this has failed?
		#Still won't know exactly where the fail is?

		#Is we get exon and a UTR result we know it must be the first or last exon
		#If we get just a UTR result we know it is just the UTR or the very 5' or 3' end
		#if we get just an exon we know it is an internal exon or an exon/intron
		#If none of the above we know it is fully intronic

		#This only counts if we do some extension
		#Otherwise the exonutr and exon values will be exactly the same!
		#So we need to test the flanks value before making this assumption!

		#Use probe dbID instead of name as this can vary across arrays
		print OUT "Unmapped $summary ".$log_name."\n";
	  
		if (!$no_triage) {		
		  $um_obj = new Bio::EnsEMBL::UnmappedObject(
													 -type       => 'probe2transcript',
													 -analysis   => $analysis,
													 -identifier => $transcript_sid,
													 -external_db_id => $transc_edb_id,
													 -summary    => "Unmapped $summary",
													 -full_desc  => "Probe mapped to $region of transcript",
													 -ensembl_object_type => 'ProbeFeature',
													 -ensembl_id => $feature_id,
													);
		  &cache_and_load_unmapped_objects($um_obj);
		}
	  }
	}
  }

  #Now get ProbeFeature IdentityXrefs for each transcript

  #my @id_xrefs = @{$dbentry_adaptor->fetch_all_by_Transcript($transcript)};


  #We could list the probe feature ID xref'd to this transcript and then pull them back
  #But what we want to do is pull back the full DBEntries directly
  #Leave this to Ian and just do thw work around?
  #Will this give us the full xref info?
  #We will have to write the get_all_DBEntries method in ProbeFeature
  #Are they implemented for all other xreffable classes? Just for SetFeatures.
  #Not for Probe, ProbeFeature or FeatureType. Put this in storable and restrict to appropriate classes?
  #ProbeFeature and all features have to inherit from Feature, so would need replicate this between Bio::EnsEMBL::Funcgen::Feature and funcgen Storable. Or can we just put the Funcgen storable first in the ISA array, so it overrides just for those functions

  #Do not store extra DBEntry here as we already have it in the IDXref?
  #Just count and deal with promiscuous probes
  #Can we get the linkage_annotation from above by match the feature ids in a cache?

  #No as we are storing IDXref on Probes not ProbeFeatures
  #Can we ditch the IDXref and just have the cigarline in probe_feature
  #We lose some information this way i.e. mismatches? True full length alignment
  #Or are we capturing this in probe_feature now?

  #Yes we need to change ProbeTranscriptAlign to use normal xrefs to ProbeFeatures rather than 
  #IDXrefs to Probes
  #This will maintain all the 'annotation' xrefs at the xref object level.
  #And all the contributing feature xrefs at the feature level.

  #Must make sure we change the rollback methods to account for this!
  
  #What is the most efficient way of retrieving the transcript link?
  #If there is an existing ProbeFeature DBEntry then count and skip overlap?
  #Skip overlap if prob feature cigar line has has D



  print OUT "\n";
}

$Helper->log("Failed to extend $failed_extend transcripts") if $failed_extend;

$Helper->log_header("Writing @array_names Xrefs", 0, 'append_date');
my $um_cnt = 0;

# now loop over all the mappings and add xrefs for those that have a suitable number of matches
#values can be a simple count or an array of annotations depending on probeset_arrays

foreach my $key (keys %transcript_feature_info) {

  my ($transcript_sid, $ensembl_id) = split (/:/, $key);
   
  #ensembl_id pname can be either probeset or probe name
  my $probeset_size = $probeset_sizes{$ensembl_id};
  #This should always be 1 for non probeset_arrays


  
  #This is the distinct number of probes, not features!
  #i.e. probe could hit twice, do we need to handle this?
  #For non-probeset arrays the key in  %{transcript_feature_info{xref_object_id}}
  #Has will be the same as the xref_object_id i.e. a probe_id
  my $hits = ($array_config{probeset_arrays}) ?  scalar(keys %{$transcript_feature_info{$key}}) 
	: scalar(@{$transcript_feature_info{$key}{$ensembl_id}});
  my $id_names = $ensembl_id.'('.join(',', @{$arrays_per_object{$ensembl_id}{names}}).')';

  if (($hits / $probeset_size) >= $mapping_threshold) {
	#This is inc'ing an undef?

	#We also need to report xref_name here for logs
	$transcripts_per_object{$ensembl_id}++;
	

	###This is why we can't run on chr slices
	#as we need to know how many transcripts have mapped to a given probeset
	#We need to chunk on probesets to parallelise
	
	
	if ($transcripts_per_object{$ensembl_id} <= $max_transcripts) {
	  
	  #So we're passing these tests for xrefs of the opposite strand
	  
	  #Can we annotate the type of Probe match? e.g. exon, exon-exon, exon-utr, utr-exon

	  #We can do this above!!
	  #But this won't account for promiscuous probes
	  #So we need to add something different to transcript_feature_count
	  #current push feature dbID, but we don't use this
	  #we already store the feature xrefs individually, but there is no way 
	  #of associating them with an probeset level xref.
	  #We could store all the probe feature dbID in the linkage annotation
	  #This is not ideal but would make it explicit which probes are contributing to a given xref
	  #This gives slight redundancy between the ProbeFeature and Probe/ProbeSet xrefs
	  #Would need to tweak the data model to handle this properly.


	  if($xref_object eq 'ProbeSet'){
		$linkage_annotation = "${hits}/${probeset_size} in ProbeSet";
	  }
	  else{

		#Hits here is number of distinct hits for a given probe dbIDs
		#Not features
		#Therefore for non-probeset arrays this will always be 1?

		if($hits > 1){
		  $linkage_annotation = "Probe matches $hits times";
		}
		else{
		  $linkage_annotation = 'Probe matches '.$transcript_feature_info{$key}{$ensembl_id}->[0];
		}
	  }

	  add_xref($transcript_sid, $ensembl_id, $xref_object, $linkage_annotation);
	  print OUT "$id_names\t$transcript_sid\tmapped\t${hits}/$probeset_size\n";
	  
	}
	else {
	  #Change this to print at end so we can add the reall total number of transcripts
	  print OUT "$id_names\t$transcript_sid\tpromiscuous\t${hits}/$probeset_size\tCurrentTranscripts".$transcripts_per_object{$ensembl_id}."\n";
	  push @{$promiscuous_objects{$ensembl_id}}, $transcript_sid;
	}
	
  } 
  else {
	#This will never happen for single probes


	print OUT "$id_names\t$transcript_sid\tinsufficient\t${hits}/${probeset_size} in ProbeSet\n";
	
	if (!$no_triage) {
	  
	  #Can/should we concentrate all unmapped info into one record
	  #Currently getting one for each probe and each probeset?
	  #Or is this a problem with array association i.e. Collapse not warked properly?
	  
	  $um_obj = new Bio::EnsEMBL::UnmappedObject(
												 -type       => 'probe2transcript',
												 -analysis   => $analysis,
												 -identifier => $transcript_sid,
												 -external_db_id => $transc_edb_id,
												 -summary    => "Insufficient hits",
												 -full_desc  => "Insufficient number of hits ${hits}/${probeset_size} in ProbeSet",
												 -ensembl_object_type => 'ProbeSet',
												 -ensembl_id => $probeset_id,
												);
	  
	  &cache_and_load_unmapped_objects($um_obj);
	}
  }
  # }
}


# Find probesets that don't match any transcripts at all, write to log file
#Why is this sub'd, we only call it once?
log_orphan_probes();


#Now update promiscuous probesets
$Helper->log("Updating ".scalar(keys %promiscuous_objects).' promiscuous probesets', 0 , 'append_date');

foreach my $object_id(keys %promiscuous_objects){

  #First delete mapped object_xrefs
  #As there is a chance that probes might be xreffed to a non-transcript entity
  #Deleting ox and x at the same time would orphan any non-transcript ox's
  $xref_db->dbc()->do("DELETE ox FROM object_xref ox, xref x, external_db edb WHERE edb.db_name='$transc_edb_name' and edb.db_release='$schema_build' and edb.external_db_id=x.external_db_id AND x.xref_id=ox.xref_id AND ox.ensembl_object_type='ProbeSet' and ox.ensembl_id='$object_id'");

  #Any other oxs?
  #if(! @{ $xref_db->dbc->db_handle->selectall_arrayref("SELECT ox.object_xref_id from object_xref ox, xref x WHERE x.dbprimary_acc='$probeset' AND x.xref_id=ox.xref_id")}){
	#Then delete xref
#	 $xref_db->dbc()->do("DELETE FROM xref WHERE dbprimary_acc='$probeset'");
#  }
  

  #Now load all unmapped objects
  #One for all arrays rather than one for each
  foreach my $transcript_sid(@{$promiscuous_objects{$object_id}}){

	$um_obj = new Bio::EnsEMBL::UnmappedObject
	  (
	   -type       => 'probe2transcript',
	   -analysis   => $analysis,
	   -identifier => $transcript_sid,
	   -summary    => "Promiscuous $xref_object",
	   -full_desc  => $xref_object.'maps to '.$transcripts_per_object{$object_id}.' transcripts (max 100)',
	   -ensembl_object_type => $xref_object,
	   -ensembl_id => $object_id,
	   -external_db_id => $transc_edb_id,
	  );
		
	&cache_and_load_unmapped_objects($um_obj);
  }

  warn 'Do we need to modify %transcript_xrefs here so we dont summary report mappings which have been deleted';

}



close (OUT);

# upload triage information if required
if (!$no_triage) {

  #$Helper->log("Uploading remaining unmapped objects to the xref DB", 0, 'append_date');
  $unmapped_object_adaptor->store(@unmapped_objects);
  $um_cnt += scalar(@unmapped_objects);
  $Helper->log("Loaded a total of $um_cnt UnmappedObjects to xref DB");
}

#Can we do this with some SQL to save memory here?

foreach my $aname(keys %array_xrefs){
  $Helper->log($aname." total xrefs mapped:\t".$array_xrefs{$aname});
}

$Helper->log('Mapped '. scalar(keys(%transcript_xrefs))."/$total transcripts ", 0, 'append_date');

#if($annotated_utrs){

#  my $method = ($calc_utrs) ? 'Calculated' : 'Hard';

#  for my $flank('5', '3'){
#	print $method." unannnotated  ${flank}_prime default UTR length used for ".$utr_defaults{$flank}." UTRs\n";
#  }
#}

$Helper->log_header("Top 5 most mapped transcripts:");

#sort keys with respect to values.
my @tids = sort { $transcript_xrefs{$b} <=>  $transcript_xrefs{$a} } keys %transcript_xrefs;
my @tcounts = sort { $b <=> $a } values %transcript_xrefs;

for my $i(0..4){
  $Helper->log("$tids[$i] mapped $tcounts[$i] times");
}

#Most (un)mapped probesets?

#$transcripts_per_object{$ensembl_id} <= $max_transcripts

$Helper->log_header("Top 5 most mapped ${xref_object}s(inc. promoscuous):");
my @xo_ids = sort { $transcripts_per_object{$b} <=>  $transcripts_per_object{$a} } keys %transcripts_per_object;
my @xo_counts = sort { $b <=> $a } values %transcripts_per_object;
my $num_ids = scalar(@xo_counts);

for my $i(0..4){
  $Helper->log("$xo_ids[$i] mapped $xo_counts[$i] times");
}


$Helper->log_header("Top 5 most mapped ${xref_object}s(no promiscuous):");
#Now we need to grep out those counts which are less than the max_transcripts
#As these are ordered we can then simply splice the relative elements from the id array
@xo_counts = grep($_ <= $max_transcripts, @xo_counts);
$num_ids = $num_ids - scalar(@xo_counts);
splice(@xo_ids, 0, $num_ids);

for my $i(0..4){
  $Helper->log("$xo_ids[$i] mapped $xo_counts[$i] times");
}

$Helper->log_header("Completed Transcript $xref_object annotation for @array_names", 0, 'append_date');



# ----------------------------------------------------------------------

# only loads unless cache hits size limit

sub cache_and_load_unmapped_objects{
  my @um_obj = @_;

  push @unmapped_objects, @um_obj;

  if(scalar(@unmapped_objects) >10000){
	#This is setting the dbID of the unmapped obj, not the unmapped reason
	#Is this really working?
	#This won't do anything as the dbID attr is ignored on store
	#What was this trying to solve?
	#Was Martin mysqlimporting and overwriting data?
	#Surely this would fail on import with duplicate keys?
	#Remove for now
	#if($first_cache){
	#  $unmapped_objects[0]->dbID('2000');
	#  $first_cache = 0;
	#}

	$um_cnt += scalar(@unmapped_objects);

	$unmapped_object_adaptor->store(@unmapped_objects);
	@unmapped_objects = ();
  }
}





# ----------------------------------------------------------------------

sub log_orphan_probes {

  if($test_slice || $test_transcript_sid){
	my $tmp = $test_slice || $test_transcript_sid;
	$Helper->log("Skipping log_orphan_probes as we are running on a test on:\t${tmp}");
	return;
  }

  $Helper->log("Logging probesets that don't map to any transcripts", 0, 'append_date');

  my ($object_id, $object_name);

  foreach my $ensembl_id(keys %arrays_per_object) {
	
    if (!$transcripts_per_object{$ensembl_id}) {

	  #($object_id, $object_name) = split/:/, $object_key;
	  
	  my $names = join(',', @{$arrays_per_object{$ensembl_id}{names}});

	  #Do we need to add dbID here in case of redundant unlined probe/set names?
      print OUT "$ensembl_id($names)\tNo transcript mappings\n";

      if (!$no_triage){
		
		#Shouldn't this be using the species external db?
		#Can we just not use an identifier and then link to the transcript DB?


		$um_obj = new Bio::EnsEMBL::UnmappedObject(-type           => 'probe2transcript',
												   -external_db_id => $transc_edb_id,
												   -analysis       => $analysis,
												   -identifier     => 'NO_TRANSCRIPT_MAPPINGS',
												   -ensembl_id     => $ensembl_id,
												   -ensembl_object_type => $xref_object,
												   -summary             => 'No transcript mappings',
												   -full_desc           => $xref_object.' did not map to any transcripts');
		
		&cache_and_load_unmapped_objects($um_obj);
      }
    }
  }
}

# ----------------------------------------------------------------------

#Change this to mappable xrefs?
#as this can be both probe or probeset?
#Will probes with same name on illumina arrays be the same?
#i.e. Are we collapsing none affy arrays?
#No!
#So we need a flag which turns this off.

sub cache_arrays_per_object {
  my $db = shift;




  $Helper->log("Caching arrays per $xref_object", 0, 'append_date');
  my $sql;#do not need distinct on count as we're grouping by array?

  if($array_config{probeset_arrays}){
	$sql = 'SELECT ps.probe_set_id, ps.name, a.name, count(p.probe_id) FROM probe p, probe_set ps, array a, array_chip ac WHERE a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and p.probe_set_id=ps.probe_set_id and a.name in ("'.join('", "', @array_names).'") GROUP BY p.probe_set_id, a.name';
	#We are not accounting for on plate replicate probes here
	#These should really be removed from the hit calculation
	#But we don't know which are replciates at that point
  }
  else{#Non probeset arrays
	#We don't really need the count at all
	#We need this for unmapped probes
	$sql = 'SELECT p.probe_id, GROUP_CONCAT(p.name SEPARATOR "#"), a.name, count(p.probe_id) FROM probe p, array a, array_chip ac WHERE a.array_id=ac.array_id and ac.array_chip_id=p.array_chip_id and a.name in ("'.join('", "', @array_names).'") GROUP BY p.probe_id, a.name';
	#GROUP_CONCAT handles nr probes with the same seq but different names

  }

  my $sth = $db->dbc()->prepare($sql);
  my ($object_id, $array, $probeset_size, $object_key, $object_name, @arrays, @names);
  $sth->execute();
  $sth->bind_columns(\$object_id, \$object_name, \$array, \$probeset_size);
  my $first_record = 1;
  my ($last_object_id, $last_probeset_size);


  while($sth->fetch()){

    if ($object_id eq $last_object_id || $first_record) {
	  
	  if(! $first_record && ($probeset_size !=  $last_probeset_size)){
		
		if(! $array_config{probeset_arrays}){
		  #If probeset size is different, this is due to identical seq
		  #Therefore we only have one probe record, but multiple names
		  $probeset_size = 1;
		  push @names, (split/#/, $object_name);
		}
		else{

		  #This is still causing problems
		  #As we can have different sizes between arrays for the same probeset
		  #e.g. RT arrays(~850 genes) have 20 probes for the same set which only have 16 on the comparable
		  #RG genomewide array
		  #This means we need to handle two different annotations!
		  #This would require being able to determine whether a give probe is available on 
		  #a given array, e.g. Probe->arrays
		  #This would be massivly redundant for most probes, so do we need a flag at the array/probeset 
		  #level to say this is a member of a variable set?
		  #Then need to generate array specific annotations

		  warn("Found probeset(dbID=$object_id) with differing size between arrays:\t@arrays($last_probeset_size) and $array($probeset_size)\n");
		  $probeset_size = $last_probeset_size if $last_probeset_size < $probeset_size;
		}
	  }
	  else{
		push @names, $object_name;
	  }
	
	  push @arrays, $array;
	  $first_record = 0;
	  $last_probeset_size = $probeset_size;
	  $last_object_id     = $object_id;
	} 
	else {
	  #Populate last entry
	  @{$arrays_per_object{$last_object_id}{arrays}} = @arrays;
	  @{$arrays_per_object{$last_object_id}{names}}  = @names;
	  $probeset_sizes{ $last_object_id }             = $last_probeset_size;
	  
	  #Start new entry
	  $last_probeset_size = $probeset_size;
	  $last_object_id     = $object_id;
	  @arrays = ($array);
	  @names  = ($object_name);
    }
  }

  #Now deal with last record
  @{$arrays_per_object{$last_object_id}{arrays}} = @arrays;
  @{$arrays_per_object{$last_object_id}{names}}  = @names;
  $probeset_sizes{ $last_object_id }             = $last_probeset_size;

  $sth->finish();
}


# ----------------------------------------------------------------------

sub pc {

  my ($a, $total) = @_;

  return "?" if (!$total);

  my $number = 100 * $a / $total;
  my $pad = "";
  $pad .= " " if ($number < 100);
  $pad .= " " if ($number < 10);

  return $pad . sprintf "%3.2f", $number;

}

# ----------------------------------------------------------------------

#We need to update this in line with new xref schema!

sub add_xref {

  my ($transcript_sid, $ensembl_id, $object_type, $linkage_annotation) = @_;
  #Maybe add linkage_type here?

 
  #Some counts
  #Add key on type of xref?
  #Will this warn as not defined?


  
  #Here the key is not the ensembl_id!
  #The key is probe/setdbID:name|MULTI_NAME

  #This defines the value as undef is the key doesn't exist
  #Were counting all xrefs here not just Probe/ProbeSet xrefs
  #But also ProbeFeature xrefs which is why we're getting probe feature ids in 
  #what was initially a probe/probeset cache

  #We don't need to count the ProbeFeature xrefs
  #We are capturing this info in UnmappedObjects

  if($object_type ne 'ProbeFeature'){
	
	foreach my $array(@{$arrays_per_object{$ensembl_id}{arrays}}){
	  $array_xrefs{$array}++;
	}
  }

  #$transcript_xrefs{$transcript_sid}{$object_type}++;
  #Only count ProbeSet xrefs
  $transcript_xrefs{$transcript_sid}++ if $object_type eq $xref_object;

  my $dbe = new Bio::EnsEMBL::DBEntry
	(
	 #xref data
	 -primary_id           => $transcript_sid,
	 -display_id           => $Helper->get_core_display_name_by_stable_id($transcript_db, $transcript_sid, 'transcript'),
	 #-version              => $schema_build,#xref.version This is implied through release of of the edb
	 #And should be actual version transcript sid?

	 #object_xref data
	 #-info_type            => "Transcript",
	 -info_type => 'MISC',
	 -info_text            => 'TRANSCRIPT',#?
	 -linkage_annotation   => $linkage_annotation,
	 -analysis             => $analysis,

	 #external_db data
	 -dbname               => $transc_edb_name,
	 -release              => $schema_build,

	 #edb stuff useless in a storage context?!
	 #-primary_id_linkable  => 1,
	 #-display_id_linkable  => 0,
	 #-priority             => 5,
	 #-db_display_name      => $transc_edb_display_name,
	 
	);


  #Can we add ox linkage annotation to this for specific transcript xref i.e. Score or how many probes hit?
  #is info_text generic for probeset xref entry or specific to a ox transcript?
  #No this is wrong!!!!!!! We are currently storing the first ox's number of probes in the xref for all ox's
  #We need to put this in the ox linkage annotation
  #DBEntryAdaptor does no handle storing or retrieving linkage annotation!!!???

  #This will enter blank entries if object_type does not enum!!

  $dbentry_adaptor->store($dbe, $ensembl_id, $object_type);
  
}

# ----------------------------------------------------------------------

# Delete existing xrefs & object xrefs & unmapped objects. Use user-specified arrays if
# defined, otherwise all arrays.
# Don't restrict to db_version as this would result in DBEntries/UnmappedObjects for old
# releases persisting.


sub delete_existing_xrefs {
  my $xref_db = shift;

  warn "Should we back up her before delete or leave this to the env?";

  my ($sql, $row_cnt);
  my $array_names = '"'.join('", "', @array_names).'"';
  my $text = "Deleting $species";
  $text .= ($array_names) ? "($array_names)" : 'ALL';
  $Helper->log("$text unmapped records and xrefs for probe2transcript...this may take a while");

  #can we pass an arrayref of ArrayChips here instead of doing it once foreach?
  
  foreach my $array(values(%arrays)){
	
		foreach my $ac(@{$array->get_ArrayChips}){
			$Helper->rollback_ArrayChip($ac, 'probe2transcript');
		}
  }

  return;
}


# ----------------------------------------------------------------------

# Check if there are already xrefs defined, and exit if there are.
# Use user-specified arrays if defined, otherwise all arrays.
# Assumes external_db.dbname == probe_array.name

sub check_existing_and_exit {
  $Helper->log_header('Checking existing Xrefs');
  warn "Need to add check for ProbeFeature xrefs or unmapped objects here to catch half done mapping";

  #We need to check for ProbeFeature xrefs and UnamppedObjects here
  #in case it has crashed half way through.

  #Or can we set status?
  #Or should we just check for ProbeFeature xrefs instead?

  my $probe_join = ($array_config{probeset_arrays}) ? 'p.probe_set_id' : 'p.probe_id';

  
  my $xref_sth = $xref_db->dbc()->prepare("SELECT COUNT(*) FROM xref x, object_xref ox, external_db e, probe p, array_chip ac, array a WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name ='${transc_edb_name}' and ox.ensembl_object_type='$xref_object' and ox.ensembl_id=${probe_join} and ox.linkage_annotation!='ProbeTranscriptAlign' and p.array_chip_id=ac.array_chip_id and ac.array_id=a.array_id and a.name=?");


  foreach my $array (@array_names) {

    $xref_sth->execute($array);
    my $cnt = $xref_sth->fetchrow_array();#Does not return array if only one value?!
	
    if ($cnt > 0) {
      print "Array $array already has $cnt xrefs, exiting.\nThere may be other arrays with xrefs. Use -delete to remove them if required.\n";
      exit(1);
    }

  }

  $xref_sth->finish();

}


sub get_or_create_analysis {

  my ($analysis_adaptor) = @_;

  my $analysis = $analysis_adaptor->fetch_by_logic_name("probe2transcript");

  if (!$analysis) {

    my $id = $analysis_adaptor->store(new Bio::EnsEMBL::Analysis(-logic_name    => 'probe2transcript',
																 -program       => 'probe2transcript.pl',
																 -description   => 'ProbeSet to Transcript mapping',
																 -displayable   => '0'));

    $analysis = $analysis_adaptor->fetch_by_logic_name("probe2transcript");

  }

  return $analysis;

}


sub usage {

  print << "EOF";

  Maps probe probes to transcripts.

  perl $0 {options}

  Options ([..] indicates optional):

  READING TRANSCRIPTS:

   --transcript_host            The database server to read transcripts from.
		
   [--transcript_port]          The port to use for reading transcripts. Defaults to 3306.
		
   --transcript_user            Database username for reading transcripts.
		
   --transcript_pass            Password for transcript_user, if required.
		
   --transcript_dbname          Database name to read transcripts from.

   [--transcript_multi_species] Indicates that the transcript database is multi-species

   [--transcript_species_id]    Species ID to use top access multi-species data

  READING ARRAYS:

   --probe_host            The database server to read probe features from.
		
   [--probe_port]          The port to use for reading probe features. Defaults to 3306.
		
   --probe_user            Database username for reading probe featuress.
		
   --probe_pass            Password for probe_user, if required.
		 
   --probe_dbname          Database name to read probe features from.

   [--probe_multi_species] Indicates that the transcript database is multi-species

   [--probe_species_id]    Species ID to use top access multi-species data

  WRITING XREFS:

   --xref_host            The database server to write xrefs to.
		
   [--xref_port]          The port to use for writing xrefs.. Defaults to 3306.
		
   --xref_user            Database username for xrefs. Must allow writing.
		
   --xref_pass            Password for xref_user, if required.
		
   --xref_dbname          Database name to write xrefs to.

   [--xref_multi_species] Indicates that the transcript database is multi-species

   [--xref_species_id]    Species ID to use top access multi-species data

  Note that if no probe_host, xref_host etc is specified, probe features will be read from,
  and xrefs written to, the database specified by the transcript_* parameters.

  Also, triage information will be written to the unmapped_object & unmapped_reason tables
  in the xref database, unless the -no_triage option is specified.

  GENERAL MAPPING OPTIONS:

  [--mismatches]      Allow up to this number of mismatches, inclusive.
                      Defaults to 1.
		
  [--utr_length]      Search this many bases downstream of the transcript
                      coding region as well. Defaults to 2000. Specify 'annotated' to use annotated lengths.
		
  [--max_probesets]   Don't store mappings to any 'promiscuous' probesets that map
                      to more than this number of transcripts. Defaults to 100.

  [--arrays]          Mandatory. Space separated list of arrays of same format e.g. AFFY or AFFY_ST etc.

  [--threshold]       Fraction of probes per probeset that have to map. Default 0.5 

  MISCELLANEOUS:

  [--delete]          Delete existing xrefs and object_xrefs, and entries in unmapped_object.
                      No deletion is done by default.

  [--force_delete]    Forces deletion of all unmapped object info even if using a subset of arrays.

  [--max_transcripts] Only use this many transcripts. Useful for debugging.

  [--no_triage]       Don't write to the unmapped_object/unmapped_reason tables.

  [--health_check]    Only do sanity checks, then stop. Useful for capthing errors before nohuping the process proper.

  [--help]            This text.


EOF

  exit(0);

}
