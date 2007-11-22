
=head1 NAME

Bio::EnsEMBL::Funcgen::Importer
  
=head1 SYNOPSIS

my $imp = Bio::EnsEMBL::Funcgen::Importer->new(%params);
$imp->register_experiment();


=head1 DESCRIPTION

B<This program> is the main class coordinating import of Arrays and experimental data.
It utilises several underlying definitions classes specific to array vendor, array class and
experimental group.  

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.ukAll sounds pretty good to me.  Maybe I'd just add that a wiggle track is available on contigview too?




=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Importer;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_date open_file run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::Parsers::ArrayDesign;
use Bio::EnsEMBL::Funcgen::Parsers::Sanger;
use Bio::EnsEMBL::Funcgen::Parsers::Nimblegen;
use Bio::EnsEMBL::Funcgen::Parsers::Solexa;
use Bio::EnsEMBL::Funcgen::Parsers::Bed;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::MAGE::XMLUtils;
use strict;
use vars qw(@ISA);

use Data::Dumper;


################################################################################

=head2 new

 Description : Constructor method

 Arg  [1]    : hash containing optional attributes:
                    -name     Name of Experiment(dir) 
                    -format   of array e.g. Tiled(default)
                    -vendor   name of array vendor
                    -description of the experiment
                    -pass DB password
		    -host DB host
		    -user  DB user
		    -port  DB port
                    -ssh  Flag to set connection over ssh via forwarded port to localhost (default = 0); remove?
                    -group    name of experimental/research group
                    -location of experimental/research group
                    -contact  e/mail address of primary contact for experimental group
                    -species 
                    -data_version  schema_build of the corresponding dnadb (change name to mirror meta_entry)
                    -recover Recovery flag (default = 0)
                    -data_dir  Root data directory (default = $ENV{'EFG_DATA'})
                    -output_dir review these dirs ???????
                    -input_dir  ?????????
                    -import_dir  ???????
                    -norm_dir    ??????
                    -fasta dump FASTA flag (default =0)
                    -array_set Flag to treat all chip designs as part of same array (default = 0)
                    -array_name Name for array set
                    -array_file Path of array file to import for sanger ENCODE array
                    -result_set_name  Name to give the raw and normalised result sets (default uses experiment and analysis name)
                    -norm_method  Normalisation method (Nimblegen default = VSN_GLOG or $ENV{'NORM_METHOD'})
                    -dbname Override for autogeneration of funcgen dbaname
                    -reg_config path to local registry config file (default = ~/ensembl.init || undef)
                    -design_type MGED term (default = binding_site_identification) get from meta/MAGE?
                    -verbose
 ReturnType  : Bio::EnsEMBL::Funcgen::Importer
 Example     : my $Exp = Bio::EnsEMBL::Importer->new(%params);
 Exceptions  : throws if mandatory params are not set or DB connect fails
 Caller      : General
 Status      : Medium - potential for %params names to change, remove %attrdata?

=cut

################################################################################

sub new{
  my ($caller) = shift;

  my $reg = "Bio::EnsEMBL::Registry";
  my $class = ref($caller) || $caller;

  my ($name, $format, $vendor, $group, $location, $contact, $species,
	  $array_name, $array_set, $array_file, $data_dir, $result_files,
	  $ftype_name, $ctype_name, $exp_date, $desc, $user, $host, $port, 
	  $pass, $dbname, $db, $data_version, $design_type, $output_dir, $input_dir,
	  $farm, $ssh, $fasta, $recover, $reg_config, $write_mage, $update_xml, 
	  $no_mage, $eset_name, $norm_method, $old_dvd_format, $feature_analysis)
	= rearrange(['NAME', 'FORMAT', 'VENDOR', 'GROUP', 'LOCATION', 'CONTACT', 'SPECIES', 
				 'ARRAY_NAME', 'ARRAY_SET', 'ARRAY_FILE', 'DATA_DIR', 'RESULT_FILES',
				 'FEATURE_TYPE_NAME', 'CELL_TYPE_NAME', 'EXPERIMENT_DATE', 'DESCRIPTION',
				 'USER', 'HOST', 'PORT', 'PASS', 'DBNAME', 'DB', 'DATA_VERSION', 'DESIGN_TYPE',
				 'OUTPUT_DIR', 'INPUT_DIR',	#to allow override of defaults
				 'FARM', 'SSH', 'DUMP_FASTA', 'RECOVER', 'REG_CONFIG', 'WRITE_MAGE', 
				 'UPDATE_XML', 'NO_MAGE', 'EXPERIMENTAL_SET_NAME', 'NORM_METHOD', 'OLD_DVD_FORMAT',
				 'FEATURE_ANALYSIS'], @_);

  
  #### Define parent parser class based on vendor
  throw("Mandatory argument -vendor not defined") if ! defined $vendor;
  my $parser_type = ucfirst(lc($vendor));
  unshift @ISA, 'Bio::EnsEMBL::Funcgen::Parsers::'.$parser_type;
  #change this to be called explicitly from the load script?

  #### Create object from parent class

  my $self = $class->SUPER::new(@_);
    
  #### Set vars and test minimum mandatory params for any import type

  $self->{'name'} = $name if $name;
  $self->vendor(uc($vendor));	#already tested
  $self->{'format'} = uc($format) || 'TILED'; #remove default?
  $self->group($group) if $group;
  $self->location($location) if $location;
  $self->contact($contact) if $contact;
  $self->{'species'} = $species || throw('Mandatory param -species not met');
  $self->array_name($array_name) if $array_name;
  $self->array_set($array_set) if $array_set;
  $self->array_file($array_file) if $array_file;
  $self->{'data_dir'} = $data_dir || $ENV{'EFG_DATA'};
  $self->result_files($result_files)if $result_files; #Sanger specific ???
  $self->{'feature_type_name'} = $ftype_name if $ftype_name;#make mandatory?
  $self->{'cell_type_name'} = $ctype_name if $ctype_name;#make madatory?
  $self->experiment_date($exp_date) if $exp_date;
  $self->description($desc) if $desc;
  $self->{'user'} = $user || throw('Mandatory param -user not met');
  $self->{'host'} = $host || 'localhost';
  $self->{'port'} = $port || '3306';
  $self->{'pass'} = $pass || throw('Mandatory param -pass not met');
  $self->dbname($dbname) if $dbname; #overrides autogeneration of dbname
  $self->db($db) if $db;		#predefined efg db
  $self->{'data_version'} = $data_version || throw('Mandatory param -data_version not met');
  $self->{'design_type'} = $design_type || 'binding_site_identification'; #remove default?
  $self->{'output_dir'} = $output_dir if $output_dir; #config default override
  $self->input_dir($input_dir) if $input_dir; #config default override
  $self->farm($farm) if $farm;
  $self->{'ssh'} = $ssh || 0;
  $self->{'_dump_fasta'} = $fasta || 0;
  $self->{'recover'} = $recover || 0;
  #check for ~/.ensembl_init to mirror general EnsEMBL behaviour
  $self->{'reg_config'} = $reg_config || ((-f "$ENV{'HOME'}/.ensembl_init") ? "$ENV{'HOME'}/.ensembl_init" : undef);
  $self->{'update_xml'} = $update_xml || 0;
  $self->{'write_mage'} = $write_mage || 0;
  $self->{'no_mage'} = $no_mage || 0;
  $self->{'experimental_set_name'} = $eset_name if $eset_name;
  $self->{'old_dvd_format'} = $old_dvd_format || 0;



  #Move all type and analysis validation here
  #use same attr?

  #Will a general norm method be applicable fo all imports?
  $self->{'norm_method'} = $norm_method || $ENV{'NORM_METHOD'};
 
 
 
  if ($self->vendor ne 'NIMBLEGEN'){
	$self->{'no_mage'} = 1;
	warn "Hardcoding no_mage for non-NIMBLEGEN imports";
  }

  #Set vendor specific attr dependent vars
  $self->set_config();

  my $host_ip = '127.0.0.1';#is this valid for all localhosts?

  ### LOAD AND RE-CONFIG REGISTRY ###
  if(defined $db){

	#need to load and set db in registry?


  }
  elsif (! defined $self->{'_reg_config'} && ! %Bio::EnsEMBL::Registry::registry_register) {
	
	#current ensembl DBs
	$reg->load_registry_from_db(
								-host => "ensembldb.ensembl.org",
								-user => "anonymous",
								-verbose => $self->verbose(),
							   );
	
	
	#why this all a bit backwards? doc please
	#should we define the eFG DB first based on the params
	#Then check the dnadb and reset if required?
	
	
	#Get standard FGDB
	$self->db($reg->get_DBAdaptor($self->species(), 'funcgen'));
      
	#reset species to standard alias to allow dbname generation
	$self->species($reg->get_alias($self->species()));
      
	#configure dnadb

	#this should be in DBAdaptor?
	#set_dnadb_by_data_version
     
	if (! $self->db() || ($self->data_version() ne $self->db->_get_schema_build($self->db()))) {
	
	  if ($self->{'ssh'}) {
	  

		$host = `host localhost`; #mac specific? nslookup localhost wont work on server/non-PC 
		#will this always be the same?

		if (! (exists $ENV{'EFG_HOST_IP'})) {
		  warn "Environment varaible EFG_HOST_IP not set for ssh mode, defaulting to $host_ip for $host";
		} else {
		  $host_ip = $ENV{'EFG_HOST_IP'};
		}
		 
		if ($self->host() ne 'localhost') {
		  warn "Overriding host ".$self->host()." for ssh connection via localhost($host_ip)";
		}
	  }

	  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-host => 'ensembldb.ensembl.org',
												-user => 'anonymous',
												-dbname => $self->species()."_core_".$self->data_version(),
												-species => $self->species(),
											   );
	} else {
	  $db = $self->db->dnadb();
	}
      
      
	$self->{'dbname'} ||= $self->species()."_funcgen_".$self->data_version();
      
	#generate and register DB with local connection settings
	$db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													   -user => $self->user(),
													   -host => ($self->{'ssh'}) ? $host_ip : $self->host(),
													   -port => $self->port(),
													   -pass => $self->pass(),
													   #we need to pass dbname else we can use non-standard dbs
													   -dbname => $self->dbname(),
													   -dnadb  => $db,
													   -species => $self->species(),
													  );
      
      
	#Redefine Fungen DB in registry
	#dnadb already added to reg via SUPER::dnadb method		
	$reg->add_DBAdaptor($self->species(), 'funcgen', $db);
	$self->db($reg->get_DBAdaptor($self->species(), 'funcgen'));
      
	throw("Unable to connect to local Funcgen DB\nPlease check the DB connect parameters and make sure the db is appropriately named") if( ! $self->db());
      
  } else {						#from config
	$reg->load_all($self->{'_reg_config'}, 1);
	$self->db($reg->get_DBAdaptor($self->species(), 'funcgen'));
	#we also need to override the registry if the dbname doesn't match that in the registry
	#we still need to reset dnadb here
	#do we need to then reset the efg DB in the reg?
  }


  #check analyses/feature_type/cell_type
  
  if($feature_analysis){
	my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($feature_analysis);
 	throw("The Feature Analysis $feature_analysis does not exist in the database") if(!$fanal);
	$self->feature_analysis($fanal);
  }




  $self->debug(2, "Importer class instance created.");
  $self->debug_hash(3, \$self);
    
  return ($self);
}


#init method kept separate from new due to differing madatory check and set up

=head2 init_array_import

  Example    : $self->init_import();
  Description: Initialises import by creating working directories 
               and by storing the Experiemnt
  Returntype : none
  Exceptions : warns and throws depending on recover and Experiment status 
  Caller     : general
  Status     : at risk - merge with register_array_design

=cut

sub init_array_import{
  # we need to define which paramters we'll be storing
  #use the logic names of the analyses as the field headers

  #need to test for vendor here

  #Sanger, NIMBLEGEN(no design_id issue, could get from the ndf, but we want it in the DesignNotes.txt)
  #Then we can change the Array/Chip generation to soley use the DesignNotes.txt rather than SampleKey 
  #which is experiment specific
  #or eFG format.


}


=head2 init_experiment_import

  Example    : $self->init_import();
  Description: Initialises import by creating working directories 
               and by storing the Experiemnt
  Returntype : none
  Exceptions : warns and throws depending on recover and Experiment status 
  Caller     : general
  Status     : at risk - merge with register exeriment

=cut

sub init_experiment_import{
  my ($self) = shift;

  foreach my $tmp ("name", "group", "data_dir") {
    throw("Mandatory arg $tmp not been defined") if (! defined $self->{$tmp});
  }
  #Should we separate path on group here too, so we can have a dev/test group?
  
  #Set and validate input dir
  $self->{'input_dir'} ||= $self->get_dir("data").'/input/'.$self->vendor().'/'.$self->name();
  throw('input_dir is not defined or does not exist ('.$self->get_dir('input').')') 
      if(! -d $self->get_dir('input')); #Helper would fail first on log/debug files

  $self->create_output_dirs('raw', 'norm', 'cache');

  throw("No result_files defined.") if (! defined $self->result_files());
  if (@{$self->result_files()}) {
    $self->log("Found result files arguments:\n\t".join("\n\t", @{$self->result_files()}));
  }


  if ($self->{'feature_type_name'}) {
    my $ftype = $self->db->get_FeatureTypeAdaptor->fetch_by_name($self->{'feature_type_name'});
    
    if (! $ftype) {
      throw("FeatureType '".$self->{'feature_type_name'}."' is not valid or is not present in the DB\n".
			"Please import using the import_type.pl script");
    }
    
    $self->feature_type($ftype);
    
  }

  if ($self->{'cell_type_name'}) {
    my $ctype = $self->db->get_CellTypeAdaptor->fetch_by_name($self->{'cell_type_name'});

	if(! $ctype){
	  throw("CellType '".$self->{'cell_type_name'}."' is not valid or is not present in the DB\n".
			"Please import using the import_type.pl script");
    }

    $self->cell_type($ctype);
  }

  
 
  #check for cell||feature and warn if no met file supplied?


  my $norm_anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($self->norm_method);

  #should we list the valid analyses?
  throw($self->norm_method.' is not a valid analysis') if ! $norm_anal;
  $self->norm_analysis($norm_anal);

  
  #check for ENV vars?
  #R_LIBS
  #R_PATH if ! farm
  #R_FARM_PATH 




  
  #Need to import to egroup here if not present and name, location & contact specified
  $self->validate_group();


  warn "Need to check env vars here?";
	
  #fetch experiment
  #if recovery and ! experiment throw 
  #else if ! experiment new and store
  
  #rename instance to name? Do we need this composite fetch?
  my $exp_adaptor = $self->db->get_ExperimentAdaptor();
  
  #print "XXXXXXXXXX featching by name ".$self->name()."\n";
  
  my $exp = $exp_adaptor->fetch_by_name($self->name());	#, $self->group());
  #should we not just do store here, as this will return the experiment if it has already been stored?




  


  #have this in a separate method?
  #we want it to generate if not present or flag specified
  #i.e. only skip generation if flag not defined 

  #separate script?  Hard to decouple meta import from Tab2MAGE generation
  #can we pass Tab2MAGE mode?
  #if false we write the template and quit, if true we read xml and populate ExperimentalChips accordingly
  #this should only write template when xml present if we specify write flag.
  #this needs to be deciphered in init_import

  


  #so we want three modes?
  #generate Tab2MAGE and exit, default if no file present of flag?
  #read Tab2MAGE, import xml and update replicates and continue. default if no xml record present or overwrite flag specified
  #skip all of above as we have already done it or are not bothered

  #change get_import_ResultSet to sets based on replicate info
  #validate this against xml?

  #then exit 
  #unless mage_xml specified or present unless generate mage_tab specified in which case start from scratch?

  #Do we really want to read this all the time?  Only if specified?
  #How can we know, test whether xml is stored in DB, overwrite but warn if already present?
  #what to do during recovery?

  #what about retrospectively changing replicates in the DB?
  #what impact would his have on result_sets and or norm results
  #none if we have the standard chip by chip norm on the 'all chip' import set
  #only affects exps which have mixed targets or design types
  #for all others we just need to update the result sets accordingly and reimport the mage-xml

  
  #~/src/bin/tab2mage.pl -e E-TABM-TEST.txt -k -t ./ -c -d PairData/


  #this should only write if not present and xml record not present
  #then if we define write mage we rewrite template only if not present

  #so first thing is to try and retrieve xml from DB
  #if true and write_mage is false skip everything
  #this would only happen on recovery
  #how are we going to handle changing xml midway through recovery

  #we should probably add an option to load a custom xml file

  #Need to validate result sets

  #can we read XML directly from DB via MAGE, NO :(
  
  #how is recovery going to be affected

  my $xml = $exp_adaptor->fetch_mage_xml_by_experiment_name($self->name());# if $self->{'write_xml'};

  #if(! $self->{'write_mage'}){#else must be a recovery? unless it's specified in error
	
	#surely we will only get xml if we are recovering!
	#this depends on when we write the xml? has to be after exp import?
	#do we have to import exp to write magetab?
	#nope, this is too messy to implement a dry run, as we'redepending on DB queries during the array/exp chip import

	#just run with recovery on ?  which kinda makes it a bit useless, not for first import when we write the mage tab
	#this means we'll have to update the replicate fields
	
	#we need some way of turning 'recovery' on if we know we've imported, maybe don't have to do this
	#if all recovery based methods after this are non critical
	#The main problem is having to roll back ExpChips which haven't actually been imported
	#We just need to create a 'LOADING' status
	#then we can omit some of the recovery methods? and not have to define recovery unless we are truly recovering
	#what about recovery methods during meta import, define another flag, or just turn recovery on?
	#turn recovery on before generating(in 2nd real import) exp, if we have xml file?
	#what if someone defines a custom file and we haven't actually imported yet?
	#This will enable someone to append/overwrite exp/result data if exp names are duplicated

	#should have full import flag to force people to be aware of possible overwrite issue?
	#this will just get turned on by default all the time rendering it useless?
	#This is a hard one to bullet proof.  Just write it and worry about protection later

	#do we have xml_id(would have to update experiment table) or experiment_id(does not follow primary key naming convention).

	#we could have this as a no fail, which would enable a custom xml file


#	if( ! $self->run_system_cmd('mysql '.$self->db->connect_string()." -e '$sql' > ".$self->get_config('mage_xml_file'), 1)){
	  #only write template if write_mage defined
#	  $got_xml = 1;
#	}
#  }

  #DO NOT CHANGE THIS LOGIC!
  #write mage if we specify or we don't have a the final xml or the template
  #recovery is turned on to stop exiting when previously stored chips
  #are found from the 'write_mage' run.
  #this does mean than if you import without running the write_mage step 
  #(This is now not possible as we test the DB for xml, not for the presence of a file?)
  #you could potentially be overwriting someone elses experiment info
  #To get around the problem of rolling back every chip, we need to add teh LOADING status
  #which should be removed once imported.

  if( ! $self->{'no_mage'}){
  
	if($self->{'write_mage'} || !( -f $self->get_config('tab2mage_file') || $xml)){
	  $self->{'write_mage'} = 1;
	  $self->backup_file($self->get_config('tab2mage_file'));
	}
	elsif($xml && (! $self->{'update_xml'})){
	  $self->{'recover'} = 1;
	  $self->{'skip_validate'} = 1;
	}
	elsif( -f $self->get_config('tab2mage_file')){#logic dictates this has to be true
	  #run tab2mage and import xml
	  #update replicate info
	  #do we need to validate xml vs meta info in read methods?
	  #turn recovery on?
	  #back up xml if present? Or just recreate?
	  #can we allow custom xml files here
	  $self->backup_file($self->get_config('mage_xml_file'));
	  
	  #do experiment check script first?
	  
	  my $cmd = 'tab2mage.pl -e '.$self->get_config('tab2mage_file').' -k -t '.$self->get_dir('output').' -c -d '.$self->get_dir('results');
	  
	  $self->log('Reading tab2mage file');
	  
	  my $t2m_exit_code = run_system_cmd($cmd, 1);#no exit flag due to non-zero exit codes
	  
	  warn "tab2mage exit code is  $t2m_exit_code"; 
	  
	  if(! ($t2m_exit_code > -1) && ($t2m_exit_code <255)){
		throw("tab2mage failed.  Please check and correct:\t".$self->get_config('tab2mage_file')."\n...and try again");
	  }
	  
	  #rename file
	  #$cmd = 'mv '.$self->get_dir('output').'/\{UNASSIGNED\}.xml '.$self->get_dir('output').'/'.$self->name().'.xml';
	  #$self->run_system_command($cmd);
	  $self->{'recover'} = 1;
	}#else{
	#this is true if we delete the tab2mage file and specify write_mage
	#we should write the mage file and exit as normal
	#	throw('Grrr, this should never be true, check the mage logic');
	#  }
  	
	
	#now we just need to use read_meta/write_mage in read methods to figure out if we need to validate or just write the template
	#we should read again if not write mage
	#should validate mage
  }
	
  
  if ($self->recovery() && ($exp)) {
    $self->log("Using previously stored Experiment:\t".$exp->name);
  } elsif ((! $self->recovery()) && $exp) {
    throw("Your experiment name is already registered in the database, please choose a different \"name\", this will require renaming you input directory, or specify -recover if you are working with a failed/partial import.");
    #can we skip this and store, and then check in register experiment if it is already stored then throw if not recovery
  } else {						# (recover && exp) || (recover  && ! exp) 
 
    
    $exp = Bio::EnsEMBL::Funcgen::Experiment->new(
												  -GROUP => $self->group(),
												  -NAME  => $self->name(),
												  -DATE  => $self->experiment_date(),
												  -PRIMARY_DESIGN_TYPE => $self->design_type(),
												  -DESCRIPTION => $self->description(),
												  -ADAPTOR => $self->db->get_ExperimentAdaptor(),
												 );
    
    ($exp) =  @{$exp_adaptor->store($exp)};
  }
  
  
  $self->experiment($exp);
  
  #remove and add specific report, this is catchig some Root stuff
  #$self->log("Initiated efg import with following parameters:\n".Data::Dumper::Dumper(\$self));
  
							
 #if we write Exp info in init then a reimport will fail unless recover is set?




  return;
}

=heead init_tab2mage_export

  Example    : $self->init_tab2mage_export;
  Description: Writes the standard experiment section of the tab2mage file
  Returntype : FileHandle
  Exceptions : ???
  Caller     : general
  Status     : at risk

=cut

sub init_tab2mage_export{
  my $self = shift;

  $self->backup_file($self->get_config('tab2mage_file')) if(-f $self->get_config('tab2mage_file'));

  my $t2m_file = open_file($self->get_config('tab2mage_file'), '>');

  #reformat this
  my $exp_section = "experiment section\ndomain\t".(split/@/, $self->contact())[1]."\naccession\t\n".
	"quality_control\tbiological_replicate\nexperiment_design_type\tbinding_site_identification\n".
	  "name\t".$self->name()."\nrelease_date\t\nsubmission_date\t\nsubmitter\t???\n".
		"submitter_email\t???\ninvestigator\t???\ninvestigator_email\t???\norganization\t???\naddress\t".
		  "???\npublication_title\t\nauthors\t\njournal\t\nvolume\t\nissue\t\npages\t\nyear\t\npubmed_id\t\n";

  my $protocol_section = "Protocol section\naccession\tname\ttext\tparameters\n";

  foreach my $protocol(sort (keys %{$self->get_config('protocols')})){
	$protocol_section .= $self->get_config('protocols')->{$protocol}->{'accession'}.
	  "\t".$self->get_config('protocols')->{$protocol}->{'name'}.
		"\t".$self->get_config('protocols')->{$protocol}->{'text'}."\t";

	$protocol_section .= (defined $self->get_config('protocols')->{$protocol}->{'parameters'}) ?
	  $self->get_config('protocols')->{$protocol}->{'parameters'}."\t\n" : "\t\n";
  }

  #File[raw]	Array[accession]	Array[serial]	Protocol[grow]	Protocol[treatment]	Protocol[extraction]	Protocol[labeling]	Protocol[hybridization]	Protocol[scanning]	BioSource	Sample	Extract	LabeledExtract	Immunoprecipitate	Hybridization	BioSourceMaterial	SampleMaterial	ExtractMaterial	LabeledExtractMaterial	Dye	BioMaterialCharacteristics[Organism]	BioMaterialCharacteristics[BioSourceType]	BioMaterialCharacteristics[StrainOrLine]	BioMaterialCharacteristics[CellType]	BioMaterialCharacteristics[Sex]	FactorValue[StrainOrLine]	FactorValue[Immunoprecipitate]


  #Need to do this bit better?
  #have array of fields.  We can then populate a hash in the read method based on field names, then use the array to print in order

  my $hyb_header = "\nHybridization section\n".join("\t", @{$self->hybridisation_fields()});

  print $t2m_file $exp_section."\n".$protocol_section."\n".$hyb_header."\n";

  return $t2m_file;
}


#Move to MAGE package?

sub hybridisation_fields{
  my $self = shift;

  return ['File[raw]', 'Array[accession]', 'Array[serial]', 
		  (map 'Protocol['.$_.']', (sort (keys %{$self->get_config('protocols')}))),
		  'BioSource', 'Sample', 'Extract', 'LabeledExtract', 'Immunoprecipitate', 'Hybridization', 
		  'BioSourceMaterial', 'SampleMaterial', 'ExtractMaterial', 'LabeledExtractMaterial',
		  'Dye', 'BioMaterialCharacteristics[Organism]', 'BioMaterialCharacteristics[BioSourceType]',	
		  'BioMaterialCharacteristics[StrainOrLine]', 'BioMaterialCharacteristics[CellType]', 
		  'BioMaterialCharacteristics[Sex]', 'FactorValue[StrainOrLine]', 'FactorValue[Immunoprecipitate]'];
}


=head2 validate_group

  Example    : $self->validate_group();
  Description: Validates groups details
  Returntype : none
  Exceptions : throws if insufficient info defined to store new Group and is not already present
  Caller     : general
  Status     : Medium - check location and contact i.e. group name clash?

=cut

sub validate_group{
  my ($self) = shift;

  my $group_ref = $self->db->fetch_group_details($self->group());

  if (! $group_ref) {
    if ($self->location() && $self->contact()) {
      $self->db->import_group($self->group(), $self->location, $self->contact());
    } else {
      throw("Group ".$self->group()." does not exist, please specify a location and contact to register the group");
    }
  }
  
  return;
}

=head2 create_output_dirs

  Example    : $self->create_output_dirs();
  Description: Does what it says on the tin, creates dirs in 
               the root output dir foreach @dirnames, also set paths in self
  Arg [1]    : mandatory - list of dir names
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Medium - add throw?

=cut

sub create_output_dirs{
  my ($self, @dirnames) = @_;
	

  #unshift @dirnames, "";#create base output dir
  #now done in control script due to log being generated first


  foreach my $name (@dirnames) {

	if($name eq 'cache'){
	  $self->{"${name}_dir"} = $ENV{'EFG_DATA'}.'/caches/'.$self->dbname() if(! defined $self->{"${name}_dir"});
	}
	else{
	  $self->{"${name}_dir"} = $self->get_dir("output")."/${name}" if(! defined $self->{"${name}_dir"});
	}

	if(! (-d $self->get_dir($name) || (-l $self->get_dir($name)))){
	  $self->log("Creating directory:\t".$self->get_dir($name));
	  mkdir $self->get_dir($name) || throw('Failed to create directory:    '. $self->get_dir($name));
	  chmod 0744, $self->get_dir($name);
	}
  }
  
  return;
}



#move this to SolexaDefs/ExperimentalSetDefs?

=head2 experimental_set_name
  
  Example    : my $esset_name = $imp->experimental_set_name();
  Description: Getter/Setter for experimental_set_name
  Arg [1]    : optional - ExperimentalSet name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub experimental_set_name{
  my $self = shift;

  $self->{'experimental_set_name'} = shift if @_;

  return $self->{'experimental_set_name'};
}



=head2 vendor
  
  Example    : $imp->vendor("NimbleGen");
  Description: Getter/Setter for array vendor
  Arg [1]    : optional - vendor name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub vendor{
  my ($self) = shift;

  if(@_){
	$self->{'vendor'} = shift;
	$self->{'vendor'} = uc($self->{'vendor'});
  }

  return $self->{'vendor'};
}


=head2 feature_type
  
  Example    : $imp->feature_type($ftype);
  Description: Getter/Setter for Experiment FeatureType
  Arg [1]    : optional - Bio::EnsEMBL::Funcgen::FeatureType
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if arg is not valid or stored
  Caller     : general
  Status     : at risk

=cut

sub feature_type{
  my ($self) = shift;

  if (@_) {
    my $ftype = shift;
    
    #do we need this as we're checking in new?
    if (! ($ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType') && $ftype->dbID())) {
      throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::FeatureType");
    }
  
    $self->{'feature_type'} = $ftype;
  }

  return $self->{'feature_type'};
}

=head2 feature_analysis
  
  Example    : $imp->feature_analysis($fanal);
  Description: Getter/Setter for Analysis used for creating the imported Features
  Arg [1]    : optional - Bio::EnsEMBL::Analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : Throws if arg is not valid or stored
  Caller     : general
  Status     : at risk

=cut

sub feature_analysis{
  my ($self) = shift;

  if (@_) {
    my $fanal = shift;
    
    #do we need this as we're checking in new?
    if (! ($fanal->isa('Bio::EnsEMBL::Analysis') && $fanal->dbID())) {
      throw("Must pass a valid stored Bio::EnsEMBL::Analysis");
    }
  
    $self->{'feature_analysis'} = $fanal;
  }

  return $self->{'feature_analysis'};
}

=head2 norm_analysis
  
  Example    : $imp->norm_analysis($anal);
  Description: Getter/Setter for the normalisation analysis
  Arg [1]    : optional - Bio::EnsEMBL::Analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : Throws if arg is not valid or stored
  Caller     : general
  Status     : at risk

=cut

sub norm_analysis{
  my ($self) = shift;

  if (@_) {
    my $anal = shift;
    
    #do we need this as we're checking in new?
    if (! (ref($anal) && $anal->isa('Bio::EnsEMBL::Analysis') && $anal->dbID())) {
      throw("Must pass a valid stored Bio::EnsEMBL::Analysis");
    }
  
    $self->{'norm_analysis'} = $anal;
  }

  return $self->{'norm_analysis'};
}



=head2 cell_type
  
  Example    : $imp->cell_type($ctype);
  Description: Getter/Setter for Experiment CellType
  Arg [1]    : optional - Bio::EnsEMBL::Funcgen::CellType
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if arg is not valid or stored
  Caller     : general
  Status     : at risk

=cut

sub cell_type{
  my ($self) = shift;

  if (@_) {
    my $ctype = shift;
    
    #do we need this as we're checking in new?
    if (! ($ctype->isa('Bio::EnsEMBL::Funcgen::CellType') && $ctype->dbID())) {
      throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::CellType");
    }
  
    $self->{'cell_type'} = $ctype;
  }

  return $self->{'cell_type'};
}

=head2 array_file
  
  Example    : my $array_file = $imp->array_file();
  Description: Getter/Setter for sanger/design array file
  Arg [1]    : optional - path to adf or gff array definition/mapping file
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub array_file{
  my ($self) = shift;
  $self->{'array_file'} = shift if(@_);
  return $self->{'array_file'};
}

=head2 array_name
  
  Example    : my $array_name = $imp->array_name();
  Description: Getter/Setter for array name
  Arg [1]    : optional string - name of array
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub array_name{
  my ($self) = shift;
  $self->{'array_name'} = shift if(@_);
  return $self->{'array_name'};
}


=head2 array_set
  
  Example    : $imp->array_set(1);
  Description: Getter/Setter for array set flag
  Arg [1]    : optional boolean - treat all array chips as the same array
  Returntype : boolean
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub array_set{
  my ($self) = shift;
  $self->{'array_set'} = shift if(@_);
  return $self->{'array_set'};
}


=head2 add_Array

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array
  Example    : $self->add_Array($array);
  Description: Setter for array elements
  Returntype : none
  Exceptions : throws if passed non Array or if more than one Array set
  Caller     : Importer
  Status     : At risk - Implement multiple arrays? Move to Experiment?

=cut

sub add_Array{
  my $self = shift;

  #do we need to check if stored?
  if (! $_[0]->isa('Bio::EnsEMBL::Funcgen::Array')) {
    throw("Must supply a Bio::EnsEMBL::Funcgen::Array");
  } elsif (@_) {
    push @{$self->{'arrays'}}, @_;
  }
   
  throw("Does not yet support multiple array imports") if(scalar (@{$self->{'arrays'}}) > 1);
  #need to alter read_probe data at the very least
  
  return;
}



=head2 arrays
  
  Example    : foreach my $array(@{$imp->arrays}){ ...do an array of things ...};
  Description: Getter for the arrays attribute
  Returntype : ARRAYREF
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub arrays{
  my $self = shift;

  if(! defined $self->{'arrays'}){
	$self->{'arrays'} = $self->db->get_ArrayAdaptor->fetch_all_by_Experiment($self->experiment());
  }

  return $self->{'arrays'};
}


=head2 data_dir
  
  Example    : $imp->data_dir($ENV{'EFG_DATA'});
  Description: Getter/Setter for root data directory
  Arg [1]    : optional - default $ENV{'EFG_DATA'}
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub data_dir{
  my ($self) = shift;
  $self->{'data_dir'} = shift if(@_);
  return $self->{'date_dir'} || $ENV{'EFG_DATA'};
}

=head2 input_dir
  
  Example    : $imp->input_dir($dir);
  Description: Getter/Setter for input directory for an experiment
  Arg [1]    : optional - input directory path
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut


sub input_dir{
  my ($self) = shift;
  $self->{'input_dir'} = shift if(@_);


  #not implemented, need to convert all get_dir calls to check if method exists or use VendorDefs
  warn "Not implmented input dir..need to cahnge all get_dir calls when we've implemented VendorDefs";

  return $self->{'input_dir'};
}

=head2 output_dir
  
  Example    : $imp->output_dir($dir);
  Description: Getter/Setter for input directory for an experiment
  Arg [1]    : optional - output directory path
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : deprecated

=cut


sub output_dir{
  my ($self) = shift;
  $self->{'output_dir'} = shift if(@_);

  throw("Deprecated, use get_dir('output') instead");

  #not implemented, need to convert all get_dir calls to check if method exists or use VendorDefs
  warn "Not implmented output dir..need to cahnge all get_dir calls when we've implemented VendorDefs";

  return $self->{'output_dir'};
}



=head2 location
  
  Example    : $imp->vendor("Hinxton");
  Description: Getter/Setter for group location
  Arg [1]    : optional - location
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub location{
  my ($self) = shift;
  $self->{'location'} = shift if(@_);
  return $self->{'location'};
}


=head2 contact
  
  Example    : my $contact = $imp->contact();
  Description: Getter/Setter for the group contact
  Arg [1]    : optional - contact name/email/address
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub contact{
  my ($self) = shift;
  $self->{'contact'} = shift if(@_);
  return $self->{'contact'};
}

=head2 name
  
  Example    : $imp->name('Experiment1');
  Description: Getter/Setter for the experiment name
  Arg [1]    : optional - experiment name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub name{
  my ($self) = shift;	
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

=head2 result_files
  
  Example    : $imp->result_files(\@files);
  Description: Getter/Setter for the result file paths
  Arg [1]    : Listref of file paths
  Returntype : Listref
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut


sub result_files{
  my ($self) = shift;	
  $self->{'result_files'} = shift if(@_);
  return $self->{'result_files'};
}




=head2 verbose
  
  Example    : $imp->verbose(1);
  Description: Getter/Setter for the verbose flag
  Arg [1]    : optional - 0 or 1
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub verbose{
  my ($self) = shift;	
  $self->{'verbose'} = shift if(@_);
  return $self->{'verbose'};
}

=head2 data_version
  
  Example    : my $schema_build = $imp->data_version();
  Description: Getter/Setter for the data version
  Arg [1]    : optional - schema and build version e.g. 41_36c
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk - rename to mirror MetaConatiner method implement reset_dbnadb?

=cut



sub data_version{
  my ($self) = shift;	

  if (@_) {
    $self->{'data_version'} = shift;
    #have reset_dnadb here?
    #Can only do this if we set data_version directly in new
    #rather than calling this method
    #as reset_dnadb assumes db is set
  }

  return $self->{'data_version'};
}


=head2 experiment_date
  
  Example    : $imp->experiment_date('2006-11-02');
  Description: Getter/Setter for the experiment date
  Arg [1]    : optional - date string in yyyy-mm-dd
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk 

=cut



sub experiment_date{
  my ($self) = shift;	

  if (@_) {
    my $date = shift;

    if ($date !~ /[0-9]{4}-[0-9]{2}[0-9]{2}/) {
      throw('Parameter -experiment_date needs to fe in the format: YYYY-MM-DD');
    }

    $self->{'experiment_date'} = $date;
  } elsif ($self->vendor() eq "nimblegen" && ! defined $self->{'experiment_date'}) {
    $self->{'experiment_date'} = &get_date("date", $self->get_config("chip_file")),
  }

  return $self->{'experiment_date'};
}



=head2 group
  
  Example    : my $exp_group = $imp->group();
  Description: Getter/Setter for the group name
  Arg [1]    : optional - group name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub group{
  my ($self) = shift;	
  $self->{'group'} = shift if(@_);
  return $self->{'group'};
}

=head2 dbname
  
  Example    : my $exp_group = $imp->group();
  Description: Getter/Setter for the group name
  Arg [1]    : optional - group name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbname{
  my ($self) = shift;	
  
  $self->{'dbname'} = shift if @_;


  if(! defined  $self->{'dbname'}){
	warn 'Need to guess dbname here?';
  }
    
  return $self->{'dbname'};
}

=head2 recovery
  
  Example    : if($imp->recovery()){ ....do recovery code...}
  Description: Getter/Setter for the recovery flag
  Arg [1]    : optional - 0 or 1
  Returntype : boolean
  Exceptions : none
  Caller     : self
  Status     : Medium - Most recovery now dynamic using status table

=cut

sub recovery{
  my $self = shift;
  $self->{'recover'} = shift if(@_);
  return $self->{'recover'};
}

=head2 description
  
  Example    : $imp->description("Human chrX H3 Lys 9 methlyation");
  Description: Getter/Setter for the experiment element
  Arg [1]    : optional - experiment description
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description{
  my $self = shift;

  if (@_) {
	$self->{'description'} = shift;
  }

  return $self->{'description'};
}

=head2 format
  
  Example    : $imp->format("Tiled");
  Description: Getter/Setter for the array format
  Arg [1]    : optional - array format
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub format{
  my ($self) = shift;	
  $self->{'format'} = shift if(@_);
  return $self->{'format'};
}

=head2 experiment
  
  Example    : my $exp = $imp->experiment();
  Description: Getter/Setter for the Experiment element
  Arg [1]    : optional - Bio::EnsEMBL::Funcgen::Experiment
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : throws if arg is not an Experiment
  Caller     : general
  Status     : Stable

=cut

sub experiment{
  my ($self) = shift;	

  if (@_) {
	
    if (! $_[0]->isa('Bio::EnsEMBL::Funcgen::Experiment')) {
      throw("Must pass a Bio::EnsEMBL::Funcgen::Experiment object");
    }

    $self->{'experiment'} = shift;
  }

  return $self->{'experiment'};
}

=head2 db
  
  Example    : $imp->db($funcgen_db);
  Description: Getter/Setter for the db element
  Arg [1]    : optional - Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : throws if arg is not an DBAdaptor
  Caller     : general
  Status     : Stable

=cut

sub db{
  my $self = shift;

  if (defined $_[0] && $_[0]->isa("Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor")) {
    $self->{'db'} = shift;
  } elsif (defined $_[0]) {
    throw("Need to pass a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor");
  }
  
  return $self->{'db'};
}

=head2 pass
  
  Example    : $imp->pass("password");
  Description: Getter/Setter for the db password
  Arg [1]    : optional - db password
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub pass{
  my $self = shift;
  $self->{'pass'} = shift if(@_);
  return $self->{'pass'};
}

=head2 pass
  
  Example    : $imp->host("hoastname");
  Description: Getter/Setter for the db hostname
  Arg [1]    : optional - db hostname
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub host{
  my $self = shift;
  $self->{'host'} = shift if(@_);
  return $self->{'host'};
}

=head2 port
  
  Example    : $imp->port(3306);
  Description: Getter/Setter for the db port number
  Arg [1]    : optional - db port number
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub port{
  my $self = shift;
  $self->{'port'} = shift if(@_);
  return $self->{'port'};
}

=head2 user
  
  Example    : $imp->user("user_name");
  Description: Getter/Setter for the db user name
  Arg [1]    : optional - db user name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub user{
  my $self = shift;
  $self->{'user'} = shift if(@_);
  return $self->{'user'};
}

=head2 dump_fasta
  
  Example    : if($self->dump_fasta()){...do fasta dump...}
  Description: Getter/Setter for the dump_fasta flag
  Arg [1]    : optional - 0 or 1
  Returntype : boolean
  Exceptions : none
  Caller     : self
  Status     : Stable

=cut


sub dump_fasta{
  my $self = shift;
  $self->{'_dump_fasta'} = shift if @_;
  return $self->{'_dump_fasta'};
}



=head2 species
  
  Example    : $imp->species("homo_sapiens");
  Description: Getter/Setter for species
  Arg [1]    : optional - species name(alias?)
  Returntype : string
  Exceptions : none ? throw if no alias found?
  Caller     : general
  Status     : Medium - may move reg alias look up to this method

=cut

sub species{
  my $self = shift;

  #should we do reg alias look up here?

  $self->{'species'} = shift if(@_);
	
  return $self->{'species'};
}

=head2 get_dir
  
  Example    : $imp->get_dir("import");
  Description: Retrieves full path for given directory
  Arg [1]    : mandatory - dir name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : at risk - move to Helper?

=cut

sub get_dir{
  my ($self, $dirname) = @_;
  return $self->get_data("${dirname}_dir");
}

=head2 norm_method
  
  Example    : my $norm_method = $imp->norm_method()
  Description: Getter/Setter for normalisation method
  Arg [1]    : mandatory - method name
  Returntype : string
  Exceptions : none ? throw if no analysis with logic name
  Caller     : general
  Status     : At risk - restrict to logic_name and validate against DB, allow multiple

=cut


sub norm_method{
  my $self = shift;

  if (@_) {
    $self->{'norm_method'} = shift;
  } elsif (! defined  $self->{'norm_method'}) {
    $self->{'norm_method'}= $self->get_config('norm_method');
  }

  return $self->{'norm_method'};
}


=head2 get_config

  Arg [1]    : mandatory - name of the data element to retrieve from the config hash
  Example    : %dye_freqs = %{$imp->get_config('dye_freqs')};
  Description: returns data from the definitions hash
  Returntype : various
  Exceptions : none
  Caller     : Importer
  Status     : at risk - replace with direct calls in the inherited Defs class?

=cut


sub get_config{
  my ($self, $data_name) = @_;
  return $self->get_data('config', $data_name); #will this cause undefs?
}





=head2 register_experiment
  
  Example    : $imp->register_experiment()
  Description: General control method, performs all data import and normalisations
  Arg [1]    : optional - dnadb DBAdaptor
  Returntype : none
  Exceptions : throws if arg is not Bio::EnsEMBL::DBSQL::DBAdaptor
  Caller     : general
  Status     : Medium

=cut


sub register_experiment{
  my ($self) = shift;
  
  #Need to check for dnadb passed with adaptor to contructor
  if (@_) { 

    if ( ! $_[0]->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")) {
      throw("You need to pass a valid dnadb adaptor to register the experiment");
    }
    $self->db->dnadb($_[0]);
  } elsif ( ! $self->db) {
    throw("You need to pass/set a DBAdaptor with a DNADB attached of the relevant data version");
  }
  
  #This could still be the default core db for the current version
  #warn here if not passed DB?
  #These should be vendor independent, only the read methods should need specific order?
  
  $self->init_experiment_import();
  #can we just have init here instead?
    
  if($self->{'write_mage'} || $self->{'no_mage'}){
	$self->read_data("array");

	if(! $self->{'no_mage'}){
	  $self->log("Please check and edit autogenerated tab2mage file:\t".$self->get_config('tab2mage_file'));
	  exit;
	}
  }elsif(! $self->{'no_mage'}){#This should be a no_channel flag, set dependent on import mode(gff_chip, gff_chan)
	#Need to accomodate chip level imports in validate!!
	$self->validate_mage() if (! $self->{'skip_validate'});
  }


  $self->read_data("probe");
  $self->read_data("results");


  

  #warn("we need to access the default method for the vendor here, or override using the gff option");
  #my $tmp_logic_name = ($self->vendor() eq "SANGER") ? "SangerPCR" : "RawValue";
  #$self->import_results("raw", $tmp_logic_name);
  
  #Need to be able to run this separately, so we can normalise previously imported sets with different methods
  #should be able t do this without raw data files e.g. retrieve info from DB

  my $norm_method = $self->norm_method();
  
  if (defined $norm_method) {
	$self->R_norm($norm_method);
    #$self->import_results("norm", $norm_method);
  }
  
  
  return;
}


=head2 validate_mage
  
  Example    : $imp->validate_mage() if(! $imp->{'write_mage'};
  Description: Validates auto-generated and manually edited mage against
               Experiment information, aswell as checking replicate defitions.
               Updates mage_xml table and replicate information accordingly.
               Any other differences are logged or an error is thrown if the 
               difference is deemed critical.
  Returntype : none
  Exceptions : throws if ...?
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : At risk

=cut


#THis is hardcoded for a channel level import at present
#Validation may fail on channels for Sanger?

sub validate_mage(){
  my ($self, $mage_xml, $update) = @_;

  $self->log("Validating mage file:\t".$self->get_config('mage_xml_file'));


  my (%echips, @log);
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name('RawValue');

  #need to change this to default
  my $vsn_anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name('VSN_GLOG');

  my $chan_rset = $self->get_import_ResultSet($anal, 'channel');
  my $rset =  $self->get_import_ResultSet($vsn_anal, 'experimental_chip');


  #doesn't really matter whether we call channel or experimental_chip?
  #yes it does, if we get a channel level set, we're going to resuse the channel cc_ids
  #for the exp_chip result sets, causing the roll back bug

  if(! $rset){
	my @tmp = @{$self->db->get_ResultSetAdaptor->fetch_all_by_name_Analysis($self->name()."_IMPORT", $anal)};

	if(scalar(@tmp) > 1){
	  throw('Found more than one IMPORT ResultSet for '.$self->name().'_IMPORT with analysis '.$anal->logic_name());
	}
	$rset = shift @tmp;
  }
  
  if(! $rset){
	throw('Cannot find ResultSet, are you trying to import a new experiment which already has a tab2mage file present?  Try removing the file, or specifying the -write_mage flag to parse_and_import.pl');
  }

  if(! -l $self->get_dir('output').'/MAGE-ML.dtd'){
	system('ln -s '.$ENV{'EFG_DATA'}.'/MAGE-ML.dtd '.$self->get_dir('output').'/MAGE-ML.dtd') == 0 || 
	  throw('Failed to link MAGE-ML.dtd');
  }
 
  $self->log('VALIDATING MAGE XML');
  my $reader = Bio::MAGE::XML::Reader->new();
  $mage_xml ||= $self->get_config('mage_xml_file');
  $self->{'mage'} = $reader->read($mage_xml);
  
  #this should only ever return 1 for an import
  foreach my $mage_exp(@{$self->{'mage'}->getExperiment_package->getExperiment_list()}){	

	if($mage_exp->getName() ne $self->name()){
	  $self->log('MAGE experiment name ('.$mage_exp->getName().') does not match import name ('.$self->name().')');
	}
	
	#add more experiment level validation here?

	foreach my $assay (@{$mage_exp->getBioAssays()}){

	  if($assay->isa('Bio::MAGE::BioAssay::PhysicalBioAssay')){#channel
		$self->log('Validating PhysicalBioAssay "'.$assay->getName()."'\n");#hyb name(this is the file name for measured assays

		my $bioassc = $assay->getBioAssayCreation();#This is a Hybridisation
		my $array = $bioassc->getArray();#this is an ArrayChip
		my $design_id = $array->getArrayDesign->getIdentifier();
		my $chip_uid = $array->getArrayIdentifier();
	

		foreach my $echip(@{$rset->get_ExperimentalChips()}){
		
		  if($echip->unique_id() eq $chip_uid){
			$self->log("Found ExperimentalChip:\t".$chip_uid);

			if(! exists $echips{$chip_uid}){
			  $echips{$chip_uid} = {(
									 total_biorep            => undef,
									 total_biotechrep        => undef,
									 experimental_biorep     => undef,
									 experimental_biotechrep => undef,
									 total_dye               => undef,
									 experimental_dye        => undef,
									 cell_type               => undef,
									 feature_type            => undef,
									)};
			}

			#Validate ArrayChip
			my ($achip) = @{$self->db->get_ArrayChipAdaptor->fetch_all_by_ExperimentalChips([$echip])};

			if($achip->design_id() ne $design_id){
			  push @log, "ArrayDesign Identifier (${design_id}) does not match ArrayChip design ID (".
				$achip->design_id().")\n\tSkipping channel and replicate validation";
			  #skip the channel/replicate validation here?	  
			} 
			else {			#validate channels and replicate names
					
			  foreach my $src_biomat (@{$bioassc->getSourceBioMaterialMeasurements()}) { #Channel materials(X1)?
				my $biomat = $src_biomat->getBioMaterial();	#LabelledExtract (IP/Control)
				#we could sub this passing $echip and biomat?
				#messy to pass regexs and populate correct echip hash attrs
				#also messy to populate log
				#keeping nested loop also prevents further obfuscation
				#do we need to do all the defined checks, or maybe just the first one?
				#Then we can skip all following warning?

				foreach my $treat (@{$biomat->getTreatments()}) {
				  #As there is effectively one more level of material extraction for the IP channel
				  #this loop will returns materials an iteration out of sync for each channel

				  foreach my $ssrc_biomat (@{$treat->getSourceBioMaterialMeasurements()}) {	#Channel measurement(x1)
					my $sbiomat = $ssrc_biomat->getBioMaterial();
					#This will either be techrep name for control of IP name for experimental channel
					#SOM0035_BR1_TR2 IP  #Immunoprecicpitate
					#SOM0035_BR1_TR2     #Extract
				
					if ($sbiomat->getName() =~ /BR[0-9]+_TR[0-9]+$/) { #Total

					  if (! defined $echips{$chip_uid}{'total_biotechrep'}) {
						$echips{$chip_uid}{'total_biotechrep'} = $sbiomat->getName();
					  }
					  else{
						push @log, "Found two TOTAL Channels on same chip with biotechreps:\t".$sbiomat->getName().
						  " and ".$echips{$chip_uid}{'total_biotechrep'};
					  }
					}else{#Experimental

					  #get feature type from assay
					  my @factor_values = @{$assay->getBioAssayFactorValues()};
					  my ($feature_type);
					  
					  foreach my $fvalue(@factor_values){
						
						if($fvalue->getValue()->getCategory() eq 'Immunoprecipitate'){
						  $feature_type = $fvalue->getName();
						  $feature_type =~ s/anti\s*-\s*//;
						  $feature_type =~ s/\s*antibody\s*//;
						}
					  }
					  $echips{$chip_uid}{'feature_type'} = $feature_type;
					}

					foreach my $ttreat (@{$sbiomat->getTreatments()}) {
									  
					  foreach my $tsrc_biomat (@{$ttreat->getSourceBioMaterialMeasurements()}) {
						my $tbiomat = $tsrc_biomat->getBioMaterial();
						#SOM0035_BR1_TR2     #Extract (exp)
						#SOM0035_BR1              #Sample (total)
									   	
						if ($tbiomat->getName() =~ /BR[0-9]+_TR[0-9]+$/) { #experimental
						  
						  if (! defined $echips{$chip_uid}{'experimental_biotechrep'}) {
							$echips{$chip_uid}{'experimental_biotechrep'} = $tbiomat->getName();
						  }
						  else{
							push @log, "Found two EXPERIMENTAL Channels on same chip with biotechreps:\t".$tbiomat->getName().
							  " and ".$echips{$chip_uid}{'experimental_biotechrep'};
						  }
						
						  my $dye = $biomat->getLabels()->[0]->getName();
							
						  foreach my $chan (@{$echip->get_Channels()}) {
							  
							if ($chan->type() eq 'EXPERIMENTAL') {
								
							  if (uc($dye) ne uc($chan->dye())) {
								push @log, "EXPERIMENTAL channel dye mismatch:\tMAGE = ".uc($dye).' vs DB '.uc($chan->dye);
							  } else {
								$echips{$chip_uid}{'experimental_dye'} = uc($dye);
							  }
							}
						  }
						} 
						else { #control
											  
						  if (! defined $echips{$chip_uid}{'total_biorep'}) {
							$echips{$chip_uid}{'total_biorep'} = $tbiomat->getName();
						  }
						  else{
							push @log, "Found two TOTAL Channels on same chip with biotechreps:\t".$tbiomat->getName().
							  " and ".$echips{$chip_uid}{'total_biorep'};
						  }
						  
						  my $dye = $biomat->getLabels()->[0]->getName();
						
						  foreach my $chan (@{$echip->get_Channels()}) {
						  
							if ($chan->type() eq 'TOTAL') {
							
							  if (uc($dye) ne uc($chan->dye())) {
								push @log, "TOTAL channel dye mismatch:\tMAGE = ".uc($dye).' vs DB '.uc($chan->dye);
							  } 
							  else {
								$echips{$chip_uid}{'total_dye'} = uc($dye);
							  }
							}
						  }
						}
						#could do one more iteration and get Source and FeatureType?
						#we should really extend this, and then update the EC cell_type and feature_types
						#these features might not be biotmats tho...need to check


						foreach my $ftreat (@{$tbiomat->getTreatments()}) {
									  
						  foreach my $fsrc_biomat (@{$ftreat->getSourceBioMaterialMeasurements()}) {
							my $fbiomat = $fsrc_biomat->getBioMaterial();
							#EXPERIMENTAL - biorep 
							#TOTAL        - source/cell type
							my $cell_type;

							if($fbiomat->getName() =~ /BR[0-9]+$/){#EXPERIMETNAL
							
							  if(! defined $echips{$chip_uid}{'experimental_biorep'}){
								$echips{$chip_uid}{'experimental_biorep'} = $fbiomat->getName();
							  }
							  else{
								push @log, "Found two Experimental Channels on same chip with bioreps:\t".$fbiomat->getName().
								  " and ".$echips{$chip_uid}{'experimental_biorep'};
							  }


							  #last treatment/measurement/biomat level should go here
							  #as TOTAL channel does not have another level and will fail
							  foreach my $xtreat (@{$fbiomat->getTreatments()}) {
								
								foreach my $xsrc_biomat (@{$xtreat->getSourceBioMaterialMeasurements()}) {
								  my $xbiomat = $xsrc_biomat->getBioMaterial();
								  
								  foreach my $char(@{$xbiomat->getCharacteristics()}){
									$cell_type = $char->getValue() if($char->getCategory() eq 'CellType');
								  }
								}
							  }

							}else{#this should be BioSource
							  #which should have CellType as characteristic
							  #we could change tab2mage and have this as a factor value, 
							  #but don't want to start messing with "standard" format
						
							  foreach my $char(@{$fbiomat->getCharacteristics()}){
								$cell_type = $char->getValue() if($char->getCategory() eq 'CellType');
							  }
							}
						
							#can have cell_type validation here
							if(! defined $echips{$chip_uid}{'cell_type'}){
							  $echips{$chip_uid}{'cell_type'} = $cell_type;
							}
							elsif( $echips{$chip_uid}{'cell_type'} ne $cell_type){
							  push @log, "Found Channels on same chip (${chip_uid}) with different cell types:\t".
								$cell_type." and ".$echips{$chip_uid}{'cell_type'};
							}
						  }
						}
					  }
					}
				  }
				}
			  }
			}
		  }						#end of echip
		}						#end of foreach echip
	  }							#end of physbioassay	
	}							#end of foreach assay
  }								#end of foreach exp



  #we should fail here with log before we update the result sets

   #we need to build rep names
  #we're currently using sample labels, in the tab2mage file
  #altho' previous sets have been using exp name
  #these have been manually patched afterwards

  #More desirable to have exp name as rset name, but no way of doing BR validation
  #based on sample label, if we don't have it in the tab2mage
  #if we change it in the DB then we need to update the tab2mage

  #no way to do this when generating tab2mage as the user hasn't yet defined the reps
  #we could just make reps based on sample labels
  #then we just assume that alterations made by the user are correct
  #as we can no longer validate using sample labels
  #can still validate using cell/feature type

  #no longer need vendor specific validation as this will be done in tab2mage generation


  #We need to validate reps here
  #the update ec records as appropriate and then create rsets

  my (%bio_reps, %tech_reps);
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
 
  foreach my $echip (@{$rset->get_ExperimentalChips()}) {

	my ($biorep, $biotechrep);

	if (! exists $echips{$echip->unique_id()}) {
	  push @log, "No MAGE entry found for ExperimentalChip:\t".$echip->unique_id();
	} 
	else {

	  foreach my $chan_type('total', 'experimental'){
		
		$biorep = $echips{$echip->unique_id()}{$chan_type.'_biorep'};
		$biotechrep = $echips{$echip->unique_id()}{$chan_type.'_biotechrep'};

		if (! defined $biotechrep) {
		  push @log, 'ExperimentalChip('.$echip->unique_id().') Extract field do not meet naming convention(SAMPLE_BRN_TRN)';
		}							#! defined biorep? will never occur at present
		elsif ($biotechrep !~ /$biorep/) {
		  push @log, "Found Extract(techrep) vs Sample(biorep) naming mismatch\t${biotechrep}\tvs$biorep";
		} 
		
		if ( ! $echips{$echip->unique_id()}{$chan_type.'_dye'}) {
		  push @log, "No ".uc($chan_type)." channel found for ExperimentalChip:\t".$echip->unique_id();
		}

	  }

	  #Is this is really implicit in the test above
	  if($echips{$echip->unique_id()}{'experimental_biorep'} ne $echips{$echip->unique_id()}{'total_biorep'}){
		push @log, "Found biorep mismatch between channels of ExperimentalChip ".$echip->unique_id().":\n".
		  "\tEXPERIMENTAL\t".$echips{$echip->unique_id()}{'experimental_biorep'}."\tTOTAL\t".
			$echips{$echip->unique_id()}{'total_biorep'};
	  }

	  #Is this is really implicit in the test above
	  if($echips{$echip->unique_id()}{'experimental_biotechrep'} ne $echips{$echip->unique_id()}{'total_biotechrep'}){
		push @log, "Found biotechrep mismatch between channels of ExperimentalChip ".$echip->unique_id().":\n".
		  "\tEXPERIMENTAL\t".$echips{$echip->unique_id()}{'experimental_biotechrep'}."\tTOTAL\t".
			$echips{$echip->unique_id()}{'total_biotechrep'};
	  }

		   
	}


	#Now we need to validate ec has same feature/cell type as other ecs in this br
	#this does not handle import sets which ARE allowed to have same name but different types
	
	if(exists $bio_reps{$biorep}){


	  if(! defined $bio_reps{$biorep}{'cell_type'}){
		push @log, "Found undefined CellType for biorep $biorep";
	  }
	  elsif($bio_reps{$biorep}{'cell_type'}->name() ne  $echips{$echip->unique_id()}{'cell_type'}){
		push @log, "Found CellType mismatch between $biorep and ExperimentalChip ".$echip->unique_id();
	  }
	  
	  
	  if(! defined $bio_reps{$biorep}{'feature_type'}){
		push @log, "Found undefined FeatureType for biorep $biorep";
	  }
	  elsif($bio_reps{$biorep}{'feature_type'}->name() ne  $echips{$echip->unique_id()}{'feature_type'}){
		push @log, "Found FeatureType mismatch between $biorep and ExperimentalChip ".$echip->unique_id();
	  } 

	}else{

	  if(defined $echips{$echip->unique_id()}{'cell_type'}){

		my $cell_type = $ct_adaptor->fetch_by_name($echips{$echip->unique_id()}{'cell_type'});

		if(! defined $cell_type){
		  push @log, 'CellType '.$echips{$echip->unique_id()}{'cell_type'}.' does not exist in the database, please use the import_type.pl script';
		}else{
		  $bio_reps{$biorep}{'cell_type'} = $cell_type;
		  $tech_reps{$biotechrep}{'cell_type'} = $cell_type;
		}
	  }else{
		warn "No CellType specified for ExperimentalChip:\t".$echip->unique_id()."\n";
	  }


	  if(defined $echips{$echip->unique_id()}{'feature_type'}){
		my $feature_type = $ft_adaptor->fetch_by_name($echips{$echip->unique_id()}{'feature_type'});

		if(! defined $feature_type){
		  push @log, 'FeatureType '.$echips{$echip->unique_id()}{'feature_type'}.' does not exist in the database, please use the import_type.pl script';
		}
		else{
		  $bio_reps{$biorep}{'feature_type'} = $feature_type;
		  $tech_reps{$biotechrep}{'feature_type'} = $feature_type;
		}
	  }else{
		warn "No FeatureType specified for ExperimentalChip:\t".$echip->unique_id()."\n";
	  }
	}

	push @{$tech_reps{$biotechrep}{'echips'}}, $echip->unique_id();
	push @{$bio_reps{$biorep}{'echips'}}, $echip->unique_id();	
  }
  



  if (@log) {
	$self->log("MAGE VALIATION REPORT\n\t".join("\n::\t", @log));
	throw("MAGE VALIDATION FAILED\nPlease correct tab2mage file and try again:\t".$self->get_config('tab2mage_file'));
  } else {
	$self->log('MAGE VALDIATION SUCCEEDED');
  }


  #we also need to build the tech rep results sets(not displayable)
  #do we need to have result sets for each biorep too?
  #update ExperimentalChip replicate info
  my (%rsets);
  my %types = (
			   feature => {},
			   cell    => {},
			  );


  warn "We need to protect against duplicating these replicate result sets";
  #fetch by exp and rset name?
  #any set with a NULL value can be duplicated


  #This needs to update and split the import/top level sets so they are of same types
  #update ec type here as we have ec context
  #careful not to update multiple times, just once for each ec
 
  my $eca = $self->db->get_ExperimentalChipAdaptor();

  foreach my $echip (@{$rset->get_ExperimentalChips()}) {
	my ($cell_type, $feature_type);

	#Set biorep info and rset
	foreach my $biorep (keys %bio_reps){

	  foreach my $chip_uid(@{$bio_reps{$biorep}{'echips'}}){

		if($chip_uid eq $echip->unique_id()){
		  $echip->biological_replicate($biorep);
		  $cell_type = $bio_reps{$biorep}{'cell_type'};
		  $feature_type = $bio_reps{$biorep}{'feature_type'};

		  if(! defined $rsets{$biorep}){
			
			$rsets{$biorep} = Bio::EnsEMBL::Funcgen::ResultSet->new
			  (
			   -NAME         => $biorep,#this may not be unique, prepend with exp name? Force method to use Experiment_and_name?
			   -ANALYSIS     => $rset->analysis(),
			   -TABLE_NAME   => 'experimental_chip',
			   -FEATURE_TYPE => $feature_type,
			   -CELL_TYPE    => $cell_type,
			  );

			#record cell and feature types
			$types{'feature'}{$feature_type->name()} = $feature_type;
			$types{'cell'}{$cell_type->name()} = $cell_type;
		  }
		  
		  $rsets{$biorep}->add_table_id($echip->dbID(), $rset->get_chip_channel_id($echip->dbID()));
		}
	  }		
	}

	#reset echip types
	$echip->feature_type($feature_type);
	$echip->cell_type($cell_type);

	
	#set tech rep info and rset
	foreach my $techrep(keys %tech_reps){
	  
	  foreach my $chip_uid(@{$tech_reps{$techrep}{'echips'}}){
		
		if($chip_uid eq $echip->unique_id()){
		  $echip->technical_replicate($techrep);

		  if(! defined $rsets{$techrep}){
			$rsets{$techrep} = Bio::EnsEMBL::Funcgen::ResultSet->new
			  (
			   -NAME       => $techrep,#this may not be unique, prepend with exp name? Force method to use Experiment_and_name?
			   -ANALYSIS   => $rset->analysis(),
			   -TABLE_NAME => 'experimental_chip',
			   -FEATURE_TYPE => $tech_reps{$techrep}{'feature_type'},
			   -CELL_TYPE    => $tech_reps{$techrep}{'cell_type'},
			  );
		  }
		  $rsets{$techrep}->add_table_id($echip->dbID(), $rset->get_chip_channel_id($echip->dbID()));
		}
	  }
	}

	$echip->adaptor->update_replicate_types($echip);#store rep info
  }

  $self->log("Created replicate ResultSets:\n::\t\t\t\t\t\t".join("\n::\t\t\t\t\t\t", map $_->name, values %rsets));


  ### Reset/Update/Clean import sets type fields
  my $sql;

  if(scalar keys %{$types{'feature'}} >1){
	$self->log('Resetting IMPORT FeatureType to NULL for multi-FeatureType Experiment');
	$sql = "UPDATE result_set set feature_type_id='NULL' where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';

  }else{
	my ($ftype) = values %{$types{'feature'}};

	if(! defined $rset->feature_type()){
	  $self->log('Updating IMPORT FeatureType to '.$ftype->name());
	  $sql = "UPDATE result_set set feature_type_id=".$ftype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
	elsif($rset->feature_type->dbID ne $ftype->dbID()){
	  warn 'FeatureType mismatch between IMPORT sets('.$rset->feature_type->name().') vs meta sets('.$ftype->name.
	  "\nUpdating to IMPORT to match meta";
	  $self->log('WARNING: FeatureType mismatch.  Updating IMPORT FeatureType('.$rset->feature_type->name().') to match meta('.$ftype->name.')');
	  $sql = "UPDATE result_set set feature_type_id=".$ftype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';

	}
  }

  $self->db->dbc->do($sql) if $sql;

  undef $sql;

  if(scalar keys %{$types{'cell'}} >1){
	$self->log('Resetting IMPORT CellType to NULL for multi-CellType Experiment');
	my $sql = "UPDATE result_set set cell_type_id='NULL' where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
  }else{
	my ($ctype) = values %{$types{'cell'}};

	if(! defined $rset->cell_type()){
	  $self->log('Updating IMPORT CellType to '.$ctype->name());
	  $sql = "UPDATE result_set set cell_type_id=".$ctype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
	elsif($rset->cell_type->dbID ne $ctype->dbID()){
	  warn 'CellType mismatch between IMPORT sets('.$rset->cell_type->name().') vs meta sets('.$ctype->name.
	  "\nUpdating to IMPORT to match meta";
	  $self->log('WARNING: FeatureType mismatch.  Updating IMPORT CellType('.$rset->cell_type->name().') to match meta('.$ctype->name.')');
	  $sql = "UPDATE result_set set cell_type_id=".$ctype->dbID()." where result_set_id in (".$rset->dbID().', '.$chan_rset->dbID().')';
	}
  }

  $self->db->dbc->do($sql) if $sql;

  ### Generate new top level sets here based on br type combos
  #we risk duplicating sets here if import set is set to one cell/featuretype
  #duplicate anyway, as import is really just for easy tracking of all chips during import

  my %toplevel_sets;
  my $toplevel_cnt = 1;
  #could tidy up toplevel_sets implmentation

  foreach my $new_rset(values %rsets){
	
	my $ftype_name = (defined $new_rset->{'feature_type'}) ? $new_rset->{'feature_type'}->name() : undef;
	my $ctype_name = (defined $new_rset->{'cell_type'}) ? $new_rset->{'cell_type'}->name() : undef;

	if(! exists $toplevel_sets{$ftype_name}){
	  $toplevel_sets{$ftype_name} = {};
	  $toplevel_sets{$ftype_name}{'feature_type'} = $new_rset->{'feature_type'};
	}



	if(! exists $toplevel_sets{$ftype_name}{$ctype_name}){
	  $toplevel_sets{$ftype_name}{$ctype_name}{'cell_type'} = $new_rset->{'cell_type'};
	  $toplevel_sets{$ftype_name}{$ctype_name}{'rsets'} = [$new_rset];
	}else{
	  push @{$toplevel_sets{$ftype_name}{$ctype_name}{'rsets'}}, $new_rset;
	}
  }



  #build toplevel sets for each feature/cell type combo using constituent rsets
  foreach my $ftype_name(keys %toplevel_sets){
	
	foreach my $ctype_name(keys %{$toplevel_sets{$ftype_name}}){
	  
	  next if $ctype_name eq 'feature_type';#skip feature type

	  $self->log("Creating toplevel ResultSet for:\t".$self->experiment->name()."\t${ftype_name}\t${ctype_name}");

	  #we need to give these a different key so we're not overwriting in the rset hash
	  $rsets{$self->experiment->name().'_'.$toplevel_cnt} = Bio::EnsEMBL::Funcgen::ResultSet->new
		(
		 -NAME       => $self->experiment->name(),
		 -ANALYSIS   => $rset->analysis(),
		 -TABLE_NAME => 'experimental_chip',
		 -FEATURE_TYPE => $toplevel_sets{$ftype_name}{'feature_type'},
		 -CELL_TYPE    => $toplevel_sets{$ftype_name}{$ctype_name}{'cell_type'},
		);

	  #add consituent table ids
	  foreach my $new_rset(@{$toplevel_sets{$ftype_name}{$ctype_name}{'rsets'}}){
		
		foreach my $ec_id(@{$new_rset->table_ids()}){

		  #Only add it if it has not already been added
		  if(!  $rsets{$self->experiment->name().'_'.$toplevel_cnt}->get_chip_channel_id($ec_id)){
			$rsets{$self->experiment->name().'_'.$toplevel_cnt}->add_table_id($ec_id, $new_rset->get_chip_channel_id($ec_id));
		  }
		}
	  }
	  $toplevel_cnt++;
	}
  }

  $self->log('Storing ResultSets');
  #Store new tech, biol and toplevel type rsets
  foreach my $new_rset(values %rsets){
	$new_rset->add_status('DAS_DISPLAYABLE');
	$rset->adaptor->store($new_rset);
  }

  my $xml_file = open_file($self->get_config('mage_xml_file'));

  #slurp in changing separator to null so we get it all in one string.
  $self->experiment->mage_xml(do{ local ($/); <$xml_file>});
  close($xml_file);

  $self->experiment($self->db->get_ExperimentAdaptor->update_mage_xml_by_Experiment($self->experiment()));

  return;
}


=head2 store_set_probes_features

  Arg [1]    : mandatory - array chip id
  Arg [2]    : optional - Bio::EnsEMBL::Funcgen::ProbeSet
  Arg [3]    : mandatory - hashref of keys probe id, values are 
               hash of probe/features with values 
               Bio::EnsEMBL::Funcgen::Probe/Features for a given 
               probe set if defined.
  Example    : $self->store_set_probes_features($ac->dbID(), $ops, \%pfs);
  Description: Stores probe set, probes and probe features 
  Returntype : none
  Exceptions : none
  Caller     : self
  Status     : Medium

=cut

sub store_set_probes_features{
  my ($self, $ac_id, $pf_hash, $ops) = @_;
  
  ### Deal with ProbeSets
  if ($ops) {
    $ops->size(scalar(keys %$pf_hash));
    ($ops) = $self->db->get_ProbeSetAdaptor->store($ops);
  }

  #If we're going to validate fully, we need to check for probes in this probeset on this array chip
  #Update size if we have any new probes
  #Overkill? Only do on recover? Do not read if array chip is IMPORTED 
  #This does not make any attempt to validate probes/set vs previously stored data
   
  for my $probe_id (keys %$pf_hash) {
    
    #set probeset in probe and store
    #the process corresponding feature
    my $probe = $pf_hash->{$probe_id}->{'probe'};
    $probe->probeset($ops) if $ops;
    ($probe) = @{$self->db->get_ProbeAdaptor->store($probe)};

	#Can't use get_all_Arrays here as we can't guarantee this will only ever be the array we've generated
    #Might dynamically load array if non-present
    #This is allowing multiple dbIDs per probe???  Is this wrong?
    #$self->cache_probe_info($probe->get_probename(), $probe->dbID());###########Remove as we're now importing all then resolving

        
    foreach my $feature (@{$pf_hash->{$probe_id}->{'features'}}) {
      $feature->probe($probe);
      ($feature) = @{$self->db->get_ProbeFeatureAdaptor->store($feature)};
    }
  }

  undef $ops;					#Will this persist in the caller?
  undef %{$pf_hash};

  return;
}


=head2 cache_slice

  Arg [0]    : string - region_name e.g. X
  Arg [1]    : optional - coordinate system name e.g. supercontig, defaults to chromosome
  Example    : my $slice = $self->cache_slice(12);
  Description: Caches or retrieves from cache a given slice
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throws f no region name specified
  Caller     : self
  Status     : At risk

=cut

sub cache_slice{
  my ($self, $region_name, $cs_name) = @_;

  throw("Need to define a region_name to cache a slice from") if ! $region_name;

  $cs_name ||= 'chromosome';
  $self->{'slice_cache'} ||= {};
  $region_name =~ s/chr//;
  $region_name = "MT" if $region_name eq "M";
  
  #can we handle UN/random chromosomes here?
  
  
  if (! exists $self->{'slice_cache'}->{$region_name}) {
    $self->{'slice_cache'}->{$region_name} = $self->slice_adaptor->fetch_by_region($cs_name, $region_name);
    warn("-- Could not generate a slice for ${cs_name}:$region_name\n") if ! defined $self->{'slice_cache'}->{$region_name};
  }
  
  return $self->{'slice_cache'}->{$region_name};
}


=head2 cache_probe_info

  Arg [0]    : mandatory - probe name
  Arg [1]    : mandatory - probe dbID
  Arg [2]    : optioanl int - x coord of probe on array
  Arg [3]    : optional int - y coord of probe on array
  Example    : $self->cache_probe_info("Probe1", $probe->dbID());
               Or for result files which do not have X & Y, we need to cache 
               X & Y from the design files: $self->cache_probe_info('Probe2', $probe->dbID(), $x, $y);
  Description: Setter for probe cache values
  Returntype : none
  Exceptions : throws is cache conflict encountered
  Caller     : self
  Status     : At risk - merge with following?

=cut

sub cache_probe_info{
  my ($self, $pname, $pid, $x, $y) = @_;

  throw('Deprecated, too memory expensive, now resolving DB duplicates and using Tied File cache');
  throw("Must provide a probe name and id") if (! defined $pname || ! defined $pid);


  #do we need to loop through the file here?
  #if we make sure we're testing for a previous dbID before storing probes then we don't need to do this
  #we can catch the error when we get the probe id as we can check for >1 id for the same probe name
  #if (defined $self->{'_probe_cache'}->{$pname} && ($self->{'_probe_cache'}->{$pname}->[0] != $pid)) {
  #  throw("Found two differing dbIDs for $pname, need to sort out redundant probe entries");
  #}
  
  $self->{'_probe_cache'}->{$pname} = (defined $x && defined $y) ? [$pid, $x, $y] : [$pid];
  
  return;
}


=head2 get_probe_id_by_name_Array

  Arg [1]    : mandatory - probe name
  Example    : $pid = $self->get_probe_id_by_name($pname);
  Description: Getter for probe cache values
  Returntype : int
  Exceptions : none
  Caller     : self
  Status     : At risk - merge with previous, move to importer?

=cut

sub get_probe_id_by_name_Array{
  my ($self, $name, $array) = @_;
  
  #this is only ever called for fully imported ArrayChips, as will be deleted if recovering
  $self->resolve_probe_data() if(! exists $self->{'_probe_cache'}{$array->name()});
  
  #we want to cycle through the given cache starting from the last position or 0.
  #we don't want to have to test for the size of the cache each time as this will be quite expensive
  #so we should store sizes separately along with current position


  my ($pid, $line);
 
  #check current line
  if($line = $self->{'_probe_cache'}{$array->name()}{'current_line'}){
	if($line =~ /^${name}\t/){
	  $pid = (split/\t/o, $line)[1];
	}
  }


  if(! $pid){
	while($line = $self->{'_probe_cache'}{$array->name()}{'handle'}->getline()){
	  
	  if($line =~ /^${name}\t/){
		$pid = (split/\t/o, $line)[1];
		$self->{'_probe_cache'}{$array->name()}{'current_line'} = $line;
		last;
	  }
	}
  }

  #do not remove this
  if(! $pid){
	throw("Did not find probe name ($name) in cache, cache may need rebuilding, results may need sorting, or do you have an anomolaous probe?")
  }else{
	chomp $pid;
  }

  return $pid;
}

=head2 get_probe_cache_by_Array

  Arg[1]     : Bio::EnsEMBL::Funcgen::Array
  Arg[2]     : boolean - from db flag, only to be used by Importer->resolve_probe_data !
  Example    : $self->get_probe_cache_by_Array();
  Description: Gets the probe info cache which is an array tied to a file
  Returntype : Boolean - True if cache has been generated and set successfully
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

#from db flag should only be used by importer
#this is because there is no guarantee that it will be resolved unless
#called by resolve_probe_data
#which then renames the file and resets the handle
#can we clean this up and protect/hide this functionality?
#can we check the cache file name in the get methods and throw if it contains unresolved?
#make this private?

sub get_probe_cache_by_Array{
  my ($self, $array, $from_db) = @_;

  my $msg = "Getting probe cache for ".$array->name();
  $msg .= " from DB" if $from_db;
  $self->log($msg);#, 1);

  if(! ($array && $array->isa('Bio::EnsEMBL::Funcgen::Array') && $array->dbID())){
	throw('Must provide a valid stored Bio::EnsEMBL::Funcgen::Array object');
  }

  my $set = 0;
  my $cache_file = $self->get_dir('cache').'/'.$array->name().'.probe_cache';

  ### Generate and resolve fresh cache from DB
  if($from_db){

	$cache_file .= '.unresolved';

	if(exists $self->{'_probe_cache'}{$array->name()}){
	  $self->log('Rebuilding probe_cache from DB for '.$array->name(), 1);


	  #untie @{$self->{'_probe_cache'}{$array->name()}{'entries'}};
	  #close($self->{'_probe_cache'}{$array->name()}{'handle'});#do we need to do this?
	  delete $self->{'_probe_cache'}{$array->name()};#implicitly closes
	  $self->log('Deleted old cache', 1);
	}else{
	  $self->log('Building probe_cache from DB for '.$array->name(), 1);
	}
	
	#Move this to ProbeAdaptor?
	#This is where we'd set the unique key for a vendor and resolves duplicates based on the key
	my $cmd = 'SELECT name, probe_id from probe WHERE array_chip_id IN ('.join(',', @{$array->get_array_chip_ids()}).') ORDER by name, probe_id';
	$cmd = 'mysql '.$self->db->connect_string()." -e \"$cmd\" >".$cache_file;
	run_system_cmd($cmd);
	
  }
 
  ### Set cache
  if(-f $cache_file){ 
	$self->log('MD5 check here?',1);
	$self->{'_probe_cache'}{$array->name()}{'current_line'} = undef;
	$self->{'_probe_cache'}{$array->name()}{'handle'} = open_file($cache_file);

	#can we do a select count instead? and do this instead of the MD5?
	#$cmd = "wc -l $cache_file";
	#my $size = `$cmd`;

	$set = 1;
  }
  else{
	warn 'Failed to get probe cache for array:'.$array->name();
  }
 
  return $set;
}


#should reorganise these emthods to split reading the array data, and the actual data
#currently:
#meta reads array and chip data
#probe reads probe_set, probes, which should definitely be in array, probe_feature? and results
#native data format may not map to these methods directly, so may need to call previous method if required data not defined


=head2 read_data
  
  Example    : $self->read_data("probe")
  Description: Calls each method in data_type array from config hash
  Arg [1]    : mandatory - data type
  Returntype : none
  Exceptions : none
  Caller     : self
  Status     : At risk

=cut

sub read_data{
  my($self, $data_type) = @_;
  map {my $method = "read_${_}_data"; $self->$method()} @{$self->get_config("${data_type}_data")};
  return;
}


=head2 design_type
  
  Example    : $self->design_type("binding_site_identification")
  Description: Getter/Setter for experimental design type
  Arg [1]    : optional - design type
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub design_type{
  my $self = shift;
  return $self->{'design_type'};
}


=head2 get_chr_seq_region_id
  
  Example    : $seq_region_id = $self->get_seq_region_id('X');
  Description: Calls each method in data_type array from config hash
  Arg [1]    : mandatory - chromosome name
  Arg [2]    : optional - start value
  Arg [3]    : optional - end value
  Returntype : int
  Exceptions : none
  Caller     : self
  Status     : At risk

=cut

#convinience wrapper method
#could we use the seq region cache instead?
#this seems like a lot of overhead for getting the id
sub get_chr_seq_region_id{
  my ($self, $chr, $start, $end) = @_;
  #what about strand info?

  #do we need the start and stop?

  #use start and stop to prevent problems with scaffodl assemblies, i.e. >1 seq_region_id
  #my $slice = $self->slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end);
  #we could pass the slice back to the slice adaptor for this, to avoid dbid problems betwen DBs
  
  ###would need to implement other cs's here

  return $self->slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end)->get_seq_region_id();
}

#convinience method

sub slice_adaptor{
  my $self = shift;

  if (! defined $self->{'slice_adaptor'}) {
	$self->{'slice_adaptor'} =  $self->db->get_SliceAdaptor();
  }

  return $self->{'slice_adaptor'};
}


=head2 vsn_norm
  
  Example    : $self->vsn_norm();
  Description: Convinience/Wrapper method for vsn R normalisation
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

#Have Norm class or contain methods in importer?
#Need to have analysis set up script for all standard analyses.

sub vsn_norm{
  my $self = shift;
  return $self->R_norm("VSN_GLOG");
}


=head2 farm
  
  Arg [1]    : Boolean
  Example    : $importer->farm(1);
  Description: Flag to turn farm submission on
  Returntype : Boolean
  Exceptions : Throws is argument not a boolean
  Caller     : general
  Status     : At risk

=cut


sub farm{
  my ($self, $farm) = @_;

  $self->{'farm'} ||= undef;#define farm

  if (defined $farm) {
    throw("Argument to farm must be a boolean 1 or 0")  if(! ($farm == 1 || $farm == 0));
    $self->{'farm'} = $farm;
  }

  return $self->{'farm'};

}


=head2 R_norm
  
  Example    : $self->R_norm(@logic_names);
  Description: Performs R normalisations for given logic names
  Returntype : none
  Exceptions : Throws if R exits with error code or if data not not valid for analysis
  Caller     : general
  Status     : At risk

=cut

sub R_norm{
  my ($self, @logic_names) = @_;
  #This currently normalises a single two colour chip at a time
  #rather than normalising across a set of chip
  #also does in sets of analyses
  #good for keeping data separate, but not efficient in terms of querying
  #convert to use one script which only queries for each echip once, then does each anal


  my $aa = $self->db->get_AnalysisAdaptor();
  my $rset_adaptor = $self->db->get_ResultSetAdaptor();
  my $ra_id = $aa->fetch_by_logic_name("RawValue")->dbID();
  my %r_config = (
				  "VSN_GLOG"      => {( libs         => ['vsn'],
										#norm_method => 'vsn',
									  )},
				  
				  "T.Biweight" => {(
									libs => ['affy'],
									#norm_method => 'tukey.biweight',
								   )},
				 );
  

  foreach my $logic_name (@logic_names) {
    #throw("Not yet implemented TukeyBiweight") if $logic_name eq "Tukey_Biweight";

	#this has already been chcecked and set as the norm_analysis
	#need to resolve this multi method approach
    my $norm_anal = $aa->fetch_by_logic_name($logic_name);


	
	#This only creates an RSet for the IMPORT set
	#So if we re-run with a different analysis
	#tab2mage will have already been validated
	#So RSet generation will be skipped
	#We need to recreate the each non-import RSet for this norm analysis
	
	#This also means the RSets are being created before the data has been imported
	#This avoids having to parse tab2mage each time but means we have an uncertain status of these Rsets

    my $rset = $self->get_import_ResultSet($norm_anal, 'experimental_chip');
	
	my @chips = ();
  
    if (! $rset) {
      $self->log("All ExperimentalChips already have status:\tIMPORTED_${logic_name}");
    } else {					#Got data to normalise and import
      my @dbids;
      my $R_file = $self->get_dir("norm")."/${logic_name}.R";
      my $job_name = $self->experiment->name()."_${logic_name}";
      my $outfile = $self->get_dir("norm")."/result.${logic_name}.txt";
      my $errfile = $self->get_dir("norm")."/${logic_name}.out";

      my $cmdline = "$ENV{'R_PATH'} --no-save < $R_file >$errfile 2>&1";
      my $bsub = "bsub -K -J $job_name ".$ENV{'R_BSUB_OPTIONS'}.
		" -e $errfile $ENV{'R_FARM_PATH'} CMD BATCH $R_file"; #--no-save?

      my $r_cmd = (! $self->farm()) ? $cmdline : $bsub;

      $self->backup_file($outfile);	#Need to do this as we're appending in the loop
  
      #setup qurey
      #warn "Need to add host and port here";
      #Set up DB, defaults and libs for each logic name
      my $query = "options(scipen=20);library(RMySQL);"; #scipen is to prevent probe_ids being converted to exponents
      
      #foreach my $ln(@logic_names){
	
      foreach my $lib (@{$r_config{$logic_name}{'libs'}}) {
		$query .= "library($lib);";
      }
      #}
      
      $query .= "con<-dbConnect(dbDriver(\"MySQL\"), host=\"".$self->host()."\", port=\"".$self->port()."\", dbname=\"".$self->db->dbc->dbname()."\", user=\"".$self->user()."\"";
      $query .= (defined $self->pass()) ? ", pass=\"".$self->pass()."\")\n" : ")\n";
      
      
      #Build queries for each chip
      foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {
    

		#should implement logic name here?
		#can't as we need seperate ResultSet for each

  
		if ($echip->has_status("IMPORTED_".$logic_name)) {
		  $self->log("ExperimentalChip ".$echip->unique_id()." already has status:\tIMPORTED_".$logic_name);
		} else {
	  
		  #warn "Need to roll back here if recovery, as norm import may fail halfway through";
		  
		  push @chips, $echip;
		  my $cc_id = $rset->get_chip_channel_id($echip->dbID());

		  #if ($self->recovery()){
		#	$self->log('Rolling back results for ExperimentalChip('.$echip->dbID().") $logic_name");
		#	$self->db->rollback_results($cc_id) if $self->recovery();							  
		#  }
		  
		  $self->log("Building $logic_name R cmd for ".$echip->unique_id());
		  @dbids = ();
	  
		  foreach my $chan (@{$echip->get_Channels()}) {

			if ($chan->type() eq "EXPERIMENTAL") {
			  push @dbids, $chan->dbID();
			} else {
			  unshift @dbids, $chan->dbID();
			}
		  }
	  
		  throw("vsn does not accomodate more than 2 channels") if (scalar(@dbids > 2) && $logic_name eq "VSN_GLOG");
	  
		  #should do some of this with maps?
		  #HARDCODED metric ID for raw data as one
	
		  #Need to get total and experimental here and set db_id accordingly
		  #can probably do this directly into one df
	  
		  $query .= "c1<-dbGetQuery(con, 'select r.probe_id as PROBE_ID, r.score as CONTROL_score, r.X, r.Y from result r, chip_channel c, result_set rs where c.table_name=\"channel\" and c.table_id=${dbids[0]} and c.result_set_id=rs.result_set_id and rs.analysis_id=${ra_id} and c.chip_channel_id=r.chip_channel_id')\n";

		  $query .= "c2<-dbGetQuery(con, 'select r.probe_id as PROBE_ID, r.score as EXPERIMENTAL_score, r.X, r.Y from result r, chip_channel c, result_set rs where c.table_name=\"channel\" and c.table_id=${dbids[1]} and c.result_set_id=rs.result_set_id and rs.analysis_id=${ra_id} and c.chip_channel_id=r.chip_channel_id')\n";
	  
	  
	
		  #MvA plot here before doing norm


		  if($logic_name eq 'T.Biweight'){

			#log2 ratios
			$query .= 'lr_df<-cbind((log(c2["EXPERIMENTAL_score"], base=exp(2)) - log(c1["CONTROL_score"], base=exp(2))))'."\n";
			
			#Adjust using tukey.biweight weighted average
			#inherits first col name
			$query .= 'norm_df<-(lr_df["EXPERIMENTAL_score"]-tukey.biweight(as.matrix(lr_df)))'."\n";
			$query .= 'formatted_df<-cbind(rep("0", length(c1["PROBE_ID"])), c1["PROBE_ID"], sprintf("%.3f", norm_df[,1]), rep("'.$cc_id.'", length(c1["PROBE_ID"])), c1["X"], c1["Y"])'."\n";
		  
		}
		elsif($logic_name eq 'VSN_GLOG'){
		  #could do this directly
		  $query .= "raw_df<-cbind(c1[\"CONTROL_score\"], c2[\"EXPERIMENTAL_score\"])\n";
		  #variance stabilise
		  $query .= "norm_df<-vsn(raw_df)\n";
    
	  
		  #do some more calcs here and print report?
		  #fold change exponentiate? See VSN docs
		  #should do someplot's of raw and glog and save here?
		  #set log func and params
		  #$query .= "par(mfrow = c(1, 2)); log.na = function(x) log(ifelse(x > 0, x, NA));";
		  #plot
		  #$query .= "plot(exprs(glog_df), main = \"vsn\", pch = \".\");". 
		  #  "plot(log.na(exprs(raw_df)), main = \"raw\", pch = \".\");"; 
		  #FAILS ON RAW PLOT!!
		  #par(mfrow = c(1, 2)) 
		  #> meanSdPlot(nkid, ranks = TRUE) 
		  #> meanSdPlot(nkid, ranks = FALSE) 
		
	
		  #Now create table structure with glog values(diffs)
		  #3 sig dec places on scores(doesn't work?!)
		  $query .= "formatted_df<-cbind(rep(\"0\", length(c1[\"PROBE_ID\"])), c1[\"PROBE_ID\"], sprintf(\"%.3f\", (exprs(norm_df[,2]) - exprs(norm_df[,1]))), rep(\"".$cc_id."\", length(c1[\"PROBE_ID\"])), c1[\"X\"], c1[\"Y\"])\n";
		  
		}
		#load back into DB
		#c3results<-cbind(rep("", length(c3["probe_id"])), c3["probe_id"], c3["c3_score"], rep(1, length(c3["probe_id"])), rep(1, length(c3["probe_id"])))
		#may want to use safe.write here
		#dbWriteTable(con, "result", c3results, append=TRUE)
		#dbWriteTable returns true but does not load any data into table!!!
		
		$query .= "write.table(formatted_df, file=\"${outfile}\", sep=\"\\t\", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)\n";
		
		#tidy up here?? 
		}
      }
  
      $query .= "q();";
  
      open(RFILE, ">$R_file") || throw("Cannot open $R_file for writing");
      print RFILE $query;
      close(RFILE);
 

	  $self->log("Submitting $logic_name job to farm:\t".localtime());
	  run_system_cmd($r_cmd);
	  $self->log("Finished $logic_name job:\t".localtime());

	  #Now load file and update status
	  #Import directly here to avoid having to reparse all results if we crash!!!!
	  $self->log("Importing:\t$outfile");
	  $self->db->load_table_data("result",  $outfile);
	  $self->log("Finishing importing:\t$outfile");
	  

	  foreach my $echip(@chips){
		$echip->adaptor->store_status("IMPORTED_${logic_name}", $echip);
	  }

	  #Recreate all non-import RSets for analysis if not already present
	  #

	  my $rset_a = $self->db->get_ResultSetAdaptor();
	  my %seen_rsets;
	 
	  foreach my $anal_rset(@{$rset_a->fetch_all_by_Experiment($self->experiment)}){
		next if($anal_rset->name =~ /_IMPORT$/);
		next if(exists $seen_rsets{$anal_rset->name});
		next if $anal_rset->analysis->logic_name eq $norm_anal->logic_name;
		$seen_rsets{$rset->name} = 1;
		$anal_rset->analysis($norm_anal);
		$anal_rset->{'dbID'} = undef;
		$anal_rset->{'adaptor'} = undef;
		
		#add the chip_channel_ids from the new anal IMPORT set
		foreach my $table_id(@{$anal_rset->table_ids}){
		  $anal_rset->{'table_id_hash'}{$table_id} = $rset->get_chip_channel_id($table_id);
		}

		$self->log('Adding new ResultSet '.$anal_rset->name.' with analysis '.$norm_anal->logic_name);
		$rset_a->store($anal_rset);
	  }



    }
  }

  return;
}

#can we sub this? args: table_name, logic_name
#also use result_set_name
#would also clean all data for result set if recovery
#return would be result_set

sub get_import_ResultSet{
  my ($self, $anal, $table_name) = @_;

  if (!($anal && $anal->isa("Bio::EnsEMBL::Analysis") && $anal->dbID())) {
    throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
  }

  $self->log("Getting import $table_name ResultSet for analysis:\t".$anal->logic_name());

  my ($rset, @new_chip_channels);
  my $result_adaptor = $self->db->get_ResultSetAdaptor();
  my $logic_name = ($anal->logic_name() eq "RawValue") ? "" : "_".$anal->logic_name();

  if(($anal->logic_name()) eq 'RawValue' && ($table_name eq 'experimental_chip')){
	throw("Cannot have an ExperimentalChip ResultSet with a RawValue analysis, either specify 'channel' or another analysis");
  }

  my $status = "IMPORTED${logic_name}";
  
  #could drop the table name here and use analysis hash?
  #warn("Need to implement chip_sets here");
  
  foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {

    #clean chip import and generate rset
	
	if($table_name eq 'experimental_chip'){
	  
	  if ($echip->has_status($status)) { #this translates to each channel have the IMPORTED_RawValue status
		$self->log("ExperimentalChip(".$echip->unique_id().") already has status:\t".$status);
	  } 
	  else {
		$self->log("Found ExperimentalChip(".$echip->unique_id().") without status $status");

		push @new_chip_channels, $echip;
	  }

	}else{#channel
	  
	  foreach my $chan(@{$echip->get_Channels()}){

		if ($chan->has_status($status)) { #this translates to each channel have the IMPORTED_RawValue status
		  $self->log("Channel(".$echip->unique_id()."_".$self->get_config('dye_freqs')->{$chan->dye()}.") already has status:\t".$status);
		} 
		else {
		  $self->log("Found Channel(".$echip->unique_id()."_".$self->get_config('dye_freqs')->{$chan->dye()}.") without status $status");
		  push @new_chip_channels, $chan;
		}
	  }
	}
  
	if (( ! $rset) && @new_chip_channels) {
	  my(@tmp) = @{$result_adaptor->fetch_all_by_name_Analysis($self->name()."_IMPORT", $anal)};

	  if(scalar(@tmp) > 1){
		throw('Found more than one IMPORT ResultSet for '.$self->name().'_IMPORT with analysis '.$anal->logic_name());
	  }

	  $rset = shift @tmp;


	  #do we need to throw here if not recovery?
	  #what if we want the import result set elsewhere during the first import?
	  
	  #if ($self->recovery()) {
	  
		#fetch by anal and experiment_id
	  #Need to change this to result_set.name!
	  #	warn("add chip set handling here");
	  
		#my @tmp = @{$result_adaptor->fetch_all_by_Experiment_Analysis($self->experiment(), $anal)};
	  #throw("Found more than one ResultSet for Experiment:\t".$self->experiment->name()."\tAnalysis:\t".$anal->logic_name().')' if (scalar(@tmp) >1);
	  #$rset = $tmp[0];
	  
	  #warn "fetching rset with ".$self->name()."_IMPORT ". $anal->logic_name;
	  
	  #$rset = $result_adaptor->fetch_by_name_Analysis($self->name()."_IMPORT", $anal);
	  warn("Warning: Could not find recovery ResultSet for analysis ".$anal->logic_name()) if ! $rset;
	  #}
	  
	  if (! $rset) {
		$self->log("Generating new ResultSet for analysis ".$anal->logic_name());
		#this is throwing if feature_type not defined
		#but we are handling multiple feature types in mage?
		#feature type does not get set even if we only see one feature_type in tab2mage
		#we need to relax this param in Set.pm and implement in FeatureSet only

		$rset = Bio::EnsEMBL::Funcgen::ResultSet->new
		  (
		   -analysis   => $anal,
		   -table_name => $table_name,
		   -name       => $self->name()."_IMPORT",
		   -feature_type => $self->feature_type(),
		   -cell_type    => $self->cell_type(),
		  );
		
		($rset) = @{$result_adaptor->store($rset)};
	  }
	}
  }

  #do we need this here as we're rolling back in the read methods?
  #we only want to roll back those chips/channels which have not been registered
  
  if ($self->recovery()) {

	my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
	
	foreach my $cc(@new_chip_channels){
	  
	  #only roll back if already part of import set
	  #Not previously registered if not 
	  if($rset->contains($cc) && $rset->get_chip_channel_id($cc->dbID())){
		
		if($table_name eq 'channel'){
		  my $chan_name = $ec_adaptor->fetch_by_dbID($cc->experimental_chip_id())->unique_id()."_".
			$self->get_config('dye_freqs')->{$cc->dye()};
		  $self->log("Rolling back results for $table_name:\t".$chan_name);
		  
		}else{
		  $self->log("Rolling back results for $table_name:\t".$cc->unique_id);
		}
		
		$self->db->rollback_results($rset->get_chip_channel_id($cc->dbID()));
	  }
	}
  }

  
  #check whether it is present in the ResultSet and add if not
  if ($rset) {
	#ids will already be present if not rset i.e. already imported

	foreach my $cc(@new_chip_channels){
	  $rset->add_table_id($cc->dbID()) if(! $rset->contains($cc));
	}
  }


  
  if ($rset) {
    $result_adaptor->store_chip_channels($rset);
  } else {
    $self->log("All ExperimentalChips have status:\t$status");
  }
  
  #this only returns a result set if there is new data to import
  return $rset;
}



=head2 resolve_probe_data

  Example    : $self->resolve_probe_data();
  Description: Resolves DB probe duplicates and builds local probe cache
  Returntype : none
  Exceptions : ????
  Caller     : general
  Status     : At risk

=cut

sub resolve_probe_data{
  my $self = shift;

  $self->log("Resolving probe data", 1);

  warn "Probe cache resolution needs to accomodate probesets too!";

  foreach my $array(@{$self->arrays()}){
	my $resolve = 0;
	
	if($self->get_probe_cache_by_Array($array)){#cache already generated

	  #check if we have any new unresolved array chips to add to the cache
	  foreach my $achip(@{$array->get_ArrayChips()}){
		
		if($achip->has_status('RESOLVED')){
		  $self->log("ArrayChip has RESOLVED status:\t".$achip->design_id());#, 1);
		  next;
		}else{
		  $self->log("Found un-RESOLVED ArrayChip:\t".$achip->design_id());
		  $resolve = 1;
		  last;
		}
	  }
	}else{#no cache file
	  $resolve = 1;
	  $self->log('No probe cache found for array '.$array->name());
	}

	if($resolve){
	  $self->log('Resolving array duplicates('.$array->name().') and rebuilding probe cache.', 1);
	  $self->get_probe_cache_by_Array($array, 1);#get from DB

	  #we need ot make sure we mark cache as unresolved, so we don't use it by mistake.
	  


	  my ($line, $name, $pid, @pids);
	  #my $index = 0;
	  my $tmp_name = '';
	  my $tmp_id = '';

	  #miss the header

	  while ($line = $self->{'_probe_cache'}{$array->name}{'handle'}->getline()){
		($name, $pid) = split/\t/o, $line;

		if($name eq $tmp_name){

		  if($pid != $tmp_id){
			push @pids, $pid;
			#should reset to pid here if we have x y data else undef
			#ignore this and force result to have x y
		  }

		  #can't do this naymore unless we figure out how to move the line pointer
		  #would still need to sed the file anyway, better to regen from DB?
		  #undef $self->{'_probe_cache'}{$array->name}{'entries'}->[$i];#delete true or to be resolved duplicate
		}
		elsif($name ne $tmp_name){#new probe
		  $self->tidy_duplicates(\@pids) if(scalar(@pids) > 1);
		  $tmp_name = $name;
		  $tmp_id = $pid;
		  @pids = ($pid);
		  #$index = $i + 1;
		}
	  }

	  $self->tidy_duplicates(\@pids) if(scalar(@pids) > 1);

	  #remove empty lines from cache
	  #only sensible way to do this is to sed the file
	  #greping !/^\n/ doubles the memory
	  #splicing may work but would have to manage indexes inside loop
	  #and will be very inefficient
	  #do sed or recreate from DB?
	  #my $cache_file = $self->get_dir('cache').'/'.$array->name().'.probe_cache';
	  #my $cmd = "sed '/^\$/d' $cache_file > ${cache_file}.tmp";
	  #$self->run_system_cmd($cmd);
	  #$cmd = "mv ${cache_file}.tmp $cache_file";
	  #$self->run_system_cmd($cmd);
	  
	  #$self->log("Deleting gappy cache here", 1);
	  #untie @{$self->{'_probe_cache'}{$array->name()}{'entries'}};
	  #delete $self->{'_probe_cache'}{$array->name()};#needs to be delete to regenerate correctly

	  $self->get_probe_cache_by_Array($array, 1);	#refresh cache from DB
	  #rename resovled cache and reset cache handle
	  my $cmd = 'mv '.$self->get_dir('cache').'/'.$array->name().'.probe_cache.unresolved '.
		$self->get_dir('cache').'/'.$array->name().'.probe_cache';
	  run_system_cmd($cmd);
	  $self->get_probe_cache_by_Array($array);


	  warn "Only generate MD5 here, as this is guranteed to be correct";
	  
	  foreach my $achip(@{$array->get_ArrayChips()}){
		
		if(! $achip->has_status('RESOLVED')){
		  $self->log("Updating ArrayChip to RESOLVED status:\t".$achip->design_id());
		  $achip->adaptor->store_status('RESOLVED', $achip);
		}
	  }
	  
	  $self->log('Finished building probe cache for '.$array->name(), 1);
	}
  }

  $self->log('Finished resolving probe data', 1);

  return;
}


sub tidy_duplicates{
  my ($self, $pids) = @_;

  my $pfa = $self->db->get_ProbeFeatureAdaptor();
  my ($feature, %features);
			
  foreach my $dup_id(@$pids){
			  
	foreach $feature(@{$pfa->fetch_all_by_Probe_id($dup_id)}){
	  #can we safely assume end will be same too?
	  push @{$features{$feature->seq_region_name().':'.$feature->start()}}, $feature;
	}
  }
			
  my (@reassign_ids, @delete_ids);
			
  foreach my $seq_start_key(keys %features){
	my $reassign_features = 1;
	
	foreach $feature(@{$features{$seq_start_key}}){
	  
	  if($feature->probe_id() == $pids->[0]){
		$reassign_features = 0;
	  }else{
		push @delete_ids, $feature->dbID();
	  }
	}
	
	#This assumes that we actually have at least one element to every seq_start_key array
	if($reassign_features){
	  my $new_fid = pop @delete_ids;
	  push @reassign_ids, $new_fid;
	}
  }
			
  #resolve features first so we don't get any orphaned features if we crash.
  $pfa->reassign_features_to_probe(\@reassign_ids, $pids->[0]) if @reassign_ids;
  $pfa->delete_features(\@delete_ids) if @delete_ids;

  return;
}

1;
