
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

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Importer;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_date);
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::Defs::DesignDefs;
use Bio::EnsEMBL::Funcgen::Defs::SangerDefs;
use Bio::EnsEMBL::Funcgen::Defs::NimblegenDefs;
use Bio::EnsEMBL::Funcgen::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use strict;
use vars qw(@ISA);




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
                    -norm_method  Normalisation method (default = vsn_norm, put defaults in Defs?)
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
	  $farm, $ssh, $fasta, $recover, $reg_config, $mage_tab)
	= rearrange(['NAME', 'FORMAT', 'VENDOR', 'GROUP', 'LOCATION', 'CONTACT', 'SPECIES', 
				 'ARRAY_NAME', 'ARRAY_SET', 'ARRAY_FILE', 'DATA_DIR', 'RESULT_FILES',
				 'FEATURE_TYPE_NAME', 'CELL_TYPE_NAME', 'EXPERIMENT_DATE', 'DESCRIPTION',
				 'USER', 'HOST', 'PORT', 'PASS', 'DBNAME', 'DB', 'DATA_VERSION', 'DESIGN_TYPE',
				 'OUTPUT_DIR', 'INPUT_DIR',	#to allow override of defaults
				 'FARM', 'SSH', 'DUMP_FASTA', 'RECOVER', 'REG_CONFIG', 'MAGE_TAB', 'WRITE_MAGE'], @_);
  #add mail flag
  #add user defined norm methods!!!!!!!!!!!!!!!!!!!!!!!!!
  #would have to make sure GroupDefs is inherited first so we can set some mandatory params
  #before checking in ArrayDefs and here

  #define parent defs class based on vendor
  throw("Mandatory argument -vendor not defined") if ! defined $vendor;
  my $defs_type = lc($vendor);
  $defs_type = ucfirst($defs_type)."Defs";
  unshift @ISA, 'Bio::EnsEMBL::Funcgen::Defs::'.$defs_type;

  #would need to set up output dir here to avoid errors on instatiation of Helper and log/debug files
  #Create object from parent class
  my $self = $class->SUPER::new(@_);

  warn "called super new on $class @ISA";

    
  #Set vars and test minimum mandatory params for any import type
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
  $self->{'feature_type_name'} = $ftype_name if $ftype_name;
  $self->{'cell_type_name'} = $ctype_name if $ctype_name;
  $self->experiment_date($exp_date) if $exp_date;
  $self->description($desc) if $desc;
  $self->{'user'} = $user || throw('Mandatory param -user not met');
  $self->{'host'} = $host || 'localhost';
  $self->{'port'} = $port || '3306';
  $self->{'pass'} = $pass || throw('Mandatory param -pass not met');
  $self->dbname($dbname) if $dbname; #overrides autogeneration of dbname
  $self->db($db) if $db;		#predefined db
  $self->{'data_version'} = $data_version || throw('Mandatory param -data_version not met');
  $self->{'design_type'} = $design_type || 'binding_site_identification'; #remove default?
  $self->output_dir($output_dir) if $output_dir; #defs default override
  $self->input_dir($input_dir) if $input_dir; #defs default override
  $self->{'farm'} = $farm || 1;
  $self->{'ssh'} = $ssh || 0;
  $self->{'fasta'} = $fasta || 0;
  $self->{'recover'} = $recover || 0;
  #check for ~/.ensembl_init to mirror general EnsEMBL behaviour
  $self->{'reg_config'} = $reg_config || ((-f "$ENV{'HOME'}/.ensembl_init") ? "$ENV{'HOME'}/.ensembl_init" : undef);
    

  #do this in init_import
  #if we write Exp info in init then a reimport will fail unless recover is set?


  if (defined $mage_tab) {
	$self->{'mage_tab'} = $mage_tab;

	#Experiment section	
	#domain	ebi.ac.uk
	#accession	
	#quality_control	biological_replicate
	#experiment_design_type	binding_site_identification_design
	#name	
	#release_date	2007-03-21
	#submission_date	2006-04-21
	#submitter	Nathan Johnson
	#submitter_email	njohnson@ebi.ac.uk
	#investigator	Nathan Johnson
	#investigator_email	njohnson@ebi.ac.uk
	#organization European Bioinformatics Institute
	#address	European Bioinformatics Institute; Hinxton; Cambridge; CB10 1SD
	#publication_title	
	#authors	
	#journal	
	#volume		
	#issue		
	#pages		
	#year
	#pubmed_id		
  }

  #Set vendor specific attr dependent vars
  $self->set_defs();


  my $host_ip = '127.0.0.1';

  ### LOAD AND RE-CONFIG REGISTRY ###
  if (! defined $self->{'_reg_config'} && ! %Bio::EnsEMBL::Registry::registry_register) {
	
	#current ensembl DBs
	$reg->load_registry_from_db(
								-host => "ensembldb.ensembl.org",
								-user => "anonymous",
								-verbose => $self->verbose(),
							   );
      
      
	#Get standard FGDB
	$self->db($reg->get_DBAdaptor($self->species(), 'funcgen'));
      
	#reset species to standard alias to allow dbname generation
	$self->species($reg->get_alias($self->species()));
      
	#configure dnadb
	#should use meta container here for schem_build/data_version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
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
  }


  $self->debug(2, "Importer class instance created.");
  $self->debug_hash(3, \$self);
    
  return ($self);
}


#Kept separate from new as it is not necessary to have native format raw data
#change name as need dir struc for processing aswell as import, may have imported in a different way
#Need to separate this further as we need still need to set the Experiment object if we're doing a re-normalise/analyse
#Move exeriment/probe/raw result import tests to register experiment?
#Make all other register methods private, so we don't bypass the previously imported exp check


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
  $self->{'input_dir'} ||= $self->get_dir("data")."/input/".$self->vendor()."/".$self->name();
  throw("input_dir is not defined or does not exist (".$self->get_dir("input").")") if(! -d $self->get_dir("input")); #Helper would fail first on log/debug files
  

  #this is now done in control script as the log file is created before the dir is :?
  if (! defined $self->get_dir("output")) {
    $self->{'output_dir'} = $self->get_dir("data")."/output/".$self->vendor()."/".$self->name();
    mkdir $self->get_dir("output") if(! -d $self->get_dir("output"));
    chmod 0744, $self->get_dir("output");
  }
  
  $self->create_output_dirs("raw", "norm");


  if (@{$self->result_files()}) {
    $self->log("Found result files arguments:\n".join("\n", @{$self->result_files()}));
  }


  if ($self->{'feature_type_name'}) {
    my $ftype = $self->db->get_FeatureTypeAdaptor->fetch_by_name($self->{'feature_type_name'});
    
    if (! $ftype) {
      throw("FeatureType '".$self->{'feature_type_name'}."' is not valid or is not present in the DB\n".
			"Please create your FeatureType before importing the experiment");
    }
    
    $self->feature_type($ftype);
    
  }

  if ($self->{'cell_type_name'}) {
    my $ctype = $self->db->get_CellTypeAdaptor->fetch_by_name($self->{'cell_type_name'});
    throw("CellType '".$self->{'cell_type_name'}."' is not valid or is not present in the DB\n".
		  "Please create your CellType before importing the experiment");
    
    $self->cell_type($ctype);
  }


  #check for cell||feature and warn if no met file supplied?


  
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
  
  return;
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
	
  #throw here

  foreach my $name (@dirnames) {
    $self->{"${name}_dir"} = $self->get_dir("output")."/${name}" if(! defined $self->{"${name}_dir"});
    mkdir $self->get_dir($name) if(! -d $self->get_dir($name));
    chmod 0744, $self->get_dir($name);
  }
  
  return;
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
  $self->{'vendor'} = shift if(@_);
  $self->{'vendor'} = uc($self->{'vendor'});
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
  Status     : at risk

=cut


sub output_dir{
  my ($self) = shift;
  $self->{'output_dir'} = shift if(@_);


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
    $self->{'experiment_date'} = &get_date("date", $self->get_def("chip_file")),
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
  
  if (@_) {
    $self->{'dbname'} = shift;
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
  $self->{'fasta'} = shift if(@_);
  return $self->{'fasta'};
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
    $self->{'norm_method'}= $self->get_def('norm_method');
  }

  return $self->{'norm_method'};
}


=head2 get_def

  Arg [1]    : mandatory - name of the data element to retrieve from the defs hash
  Example    : %dye_freqs = %{$imp->get_def('dye_freqs')};
  Description: returns data from the definitions hash
  Returntype : various
  Exceptions : none
  Caller     : Importer
  Status     : at risk - replace with direct calls in the inherited Defs class?

=cut


sub get_def{
  my ($self, $data_name) = @_;
  return $self->get_data('defs', $data_name); #will this cause undefs?
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
  } elsif ( ! $self->db()) {
    throw("You need to pass/set a DBAdaptor with a DNADB attached of the relevant data version");
  }
  
  #This could still be the default core db for the current version
  #warn here if not passed DB?
  #These should be vendor independent, only the read methods should need specific order?
  #Need to totally separate parse/read from import, so we can just do one method if required, i.e. normalise
  #Also need to move id generation to import methods, which would check that a previous import has been done, or check the DB for the relevant info?
  
  $self->init_experiment_import();

  #check here is exp already stored?  Will this work properly?
  #then throw if not recovery
  #else do following, but change to _private type methods
  #as bypassing this register_experiment method and calling directly will cause problems with duplication of data
  #we need to make this more stringent, maybe we can do a caller in each method to make sure it is register experiment calling each method
  
  $self->read_data("array");
  $self->read_data("probe");
  $self->read_data("results");

  warn("we need to access the default method for the vendor here, or override using the gff option");
  my $tmp_logic_name = ($self->vendor() eq "SANGER") ? "SangerPCR" : "RawValue";
  $self->import_results("raw", $tmp_logic_name);
  
  #Need to be able to run this separately, so we can normalise previously imported sets with different methods
  #should be able t do this without raw data files e.g. retrieve info from DB

  my $norm_method = $self->norm_method();
  
  if (defined $norm_method) {
    #$self->$norm_method;
    $self->R_norm($norm_method);
    $self->import_results("norm", $norm_method);
  }
  
  
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
  #my ($self, $ac_id, $ops, $probes, $features) = @_;

  my ($self, $ac_id, $pf_hash, $ops) = @_;
  

  #just call these directly rather than setting each time?
  #my $op_a = $self->db->get_ProbeAdaptor();
  #my $opf_a = $self->db->get_ProbeFeatureAdaptor();


  #if(scalar(@$probes) != scalar(@$features)){
  #  my $t = scalar(@$probes);
  #  my $f = scalar(@$features);
  #  throw("Does not accomodate multiple features($f) per probe $t");
  #}

  #Maybe do some validation here against all probes for probeset and ac_id? 

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
    $self->cache_probe_info($probe->get_probename(), $probe->dbID());

        
    foreach my $feature (@{$pf_hash->{$probe_id}->{'features'}}) {
      $feature->probe($probe);
      ($feature) = @{$self->db->get_ProbeFeatureAdaptor->store($feature)};
    }
  }

  undef $ops;					#Will this persist in the caller?
  undef %{$pf_hash};

  # undef @$probes;
  #  undef @$features,

  return;
}

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

  throw("Must provide a probe name and id") if (! defined $pname || ! defined $pid);


  if (defined $self->{'_probe_cache'}->{$pname} && ($self->{'_probe_cache'}->{$pname}->[0] != $pid)) {
    throw("Found two differing dbIDs for $pname, need to sort out redundant probe entries");
  }
  
  $self->{'_probe_cache'}->{$pname} = (defined $x && defined $y) ? [$pid, $x, $y] : [$pid];
  

 
  return;
}

=head2 get_probe_id_by_name

  Arg [1]    : mandatory - probe name
  Example    : $pid = $self->get_probe_id_by_name($pname);
  Description: Getter for probe cache values
  Returntype : int
  Exceptions : none
  Caller     : self
  Status     : At risk - merge with previous, move to importer?

=cut


sub get_probe_id_by_name{
  my ($self, $name) = @_;
  
  #this is only ever called for fully imported ArrayChips, as will be deleted if recovering
  $self->build_probe_info_cache() if(! defined $self->{'_probe_cache'});  

  return (exists $self->{'_probe_cache'}->{$name}) ? $self->{'_probe_cache'}->{$name}->[0] : undef;
}

=head2 build_probe_info_cache

  Example    : $self->build_probe_info_cache();
  Description: Build the probe info cache from the DB
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : At risk 

=cut


sub build_probe_info_cache{
  my $self = shift;

  $self->log('Building probe_info_cache from DB');

  if (defined $self->{'_probe_cache'}) {
	throw('If you really want to over-write the probe info cache, you must explicitly undef the cache first');
  }

  my %tmp = %{$self->db->get_ProbeAdaptor->fetch_probe_cache_by_Experiment($self->experiment())};

  #map to the cache here 
  map $self->{'_probe_cache'}{$_} = [ $tmp{$_} ] , keys %tmp;

  return;
}


=head2 get_probe_x_y_by_name

  Arg [1]    : mandatory - probe name
  Example    : ($x, $y) = @{$self->get_probe_x_y_by_name($pname)};
  Description: Getter for probe x y cache values
  Returntype : LISTREF
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut


sub get_probe_x_y_by_name{
  my ($self, $name) = @_;
  
  my @xy;

  #This does not test for x & y as we don't want to slow it down for every single probe call
  #slightly odd behaviour have it still build the x/yless cache from the DB even tho it will return undef
  #This is because we don't want it ti be testable for rebuilding the x y cache

  if (! defined $self->{'_probe_cache'}) {
	$self->build_probe_info_cache();
  }

  if (exists  $self->{'_probe_cache'}->{$name} &&
	  scalar(@{$self->{'_probe_cache'}->{$name}}) == 3) {
	@xy = ($self->{'_probe_cache'}->{$name}->[1], $self->{'_probe_cache'}->{$name}->[2]);
  }

  return \@xy;
}


#the generic read_methods should go in here too?
#should reorganise these emthods to split reading the array data, and the actual data
#currently:
#meta reads array and chip data
#probe reads probe_set, probes, which should definitely be in array, probe_feature? and results
#native data format may not map to these methods directly, so may need to call previous method if required data not defined

=head2 import_results
  
  Example    : $self->import_results()
  Description: Imports results into DB from file
  Arg [1]    : mandatory - results dir
  Returntype : none
  Exceptions : throws if R
  Caller     : general
  Status     : Medium

=cut

sub import_results{
  my ($self, $dir, @logic_names) = @_;


  foreach my $logic_name (@logic_names) {
    $self->log("Importing $logic_name data");


    my $loaded = 0;
    my $status = ($logic_name eq "RawValue") ? "IMPORTED" : "IMPORTED_${logic_name}";
  
    
    #should we split this into vendor specific load methods?
    #This loop behaves wierdly due to the naming/array design of each platform
    #Nimblegen loads all experimental chips from an arraychip(can be All_pair if only one Achip)
    #Sanger has individual Echip files
    
    foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {
      
      if ( ! $echip->has_status($status)) {	#must have some results

		if ($dir ne "norm") {
	  
		  if ($self->vendor() eq "NIMBLEGEN") {
	      
			#if(! $loaded){

			#  foreach my $array(@{$self->arrays()}){#load data for all echips
		
			#foreach my $design_id(@{$array->get_design_ids()}){
			#  my $achip = $array->get_ArrayChip_by_design_id($design_id);
			foreach my $chan (@{$echip->get_Channels()}) {
			  my $chan_name = $echip->unique_id()."_".$self->get_def('dye_freqs')->{$chan->dye()};
			  $self->log("Loading $logic_name data for Channel:\t$chan_name");
			  $self->db->load_table_data("result",  $self->get_dir($dir)."/result.${chan_name}.txt");
			  $self->log("Finishing loading $logic_name data for Channel:\t${chan_name}");
			}
			# }
			#}
		  } elsif ($self->vendor() eq "SANGER") { #this is norm data, but we're importing norm, rather than processing.
			warn ("Need to fix the sanger methods so imports as norm, not raw");
			#IMPORTED is not correct status for these Echips, well, should have
			#IMPORTED?  and IMPORTED_NORM_METHOD
	    
	    
			$self->log("Loading $logic_name data for ExperimentalChip:\t".$echip->unique_id());
			$self->db->load_table_data("result",  $self->get_dir('norm')."/result.${logic_name}.".$echip->unique_id().".txt");
			$self->log("Finishing loading $logic_name data for ExperimentalChip:\t".$echip->unique_id());
		  }
		} elsif (! $loaded) {	#norm data, imports all echip/set data as one
		  #should add logic_name here
	  
		  $self->log("Loading $logic_name norm data from ".$self->get_dir($dir)."/result.${logic_name}.txt");
		  $self->db->load_table_data("result",  $self->get_dir($dir)."/result.${logic_name}.txt");
		  $self->log("Finished loading $logic_name norm data");
		}
      
	
		$echip->adaptor->set_status($status, $echip);
		$loaded = 1;	
      }
    }

    $self->log("Finished loading $logic_name $dir results");

  }
  
  
  $self->log("Finished loading results ");
  return;
}

=head2 read_data
  
  Example    : $self->read_data("probe")
  Description: Calls each method in data_type array from defs hash
  Arg [1]    : mandatory - data type
  Returntype : none
  Exceptions : none
  Caller     : self
  Status     : At risk

=cut

sub read_data{
  my($self, $data_type) = @_;
  map {my $method = "read_${_}_data"; $self->$method()} @{$self->get_def("${data_type}_data")};
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
  Description: Calls each method in data_type array from defs hash
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

  $self->{'farm'} ||= undef;

  if ($farm) {

    throw("Argument to farm must be a boolean 1 or 0")  if($farm != 1 || $farm != 0);

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
  my %r_libs = (
				"VSN_GLOG"      => ['vsn'],
				"TukeyBiweight" => ['affy'],
			   );




  foreach my $logic_name (@logic_names) {
    throw("Not yet implemented TukeyBiweight") if $logic_name eq "TukeyBiweight";
    my $norm_anal = $aa->fetch_by_logic_name($logic_name);
    my $rset = $self->get_import_ResultSet('experimental_chip', $norm_anal);
  
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
	
      foreach my $lib (@{$r_libs{$logic_name}}) {
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
	  
		  my $cc_id = $rset->get_chip_channel_id($echip->dbID());
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
	
	  
		  $query .= "c1<-dbGetQuery(con, 'select r.probe_id, r.score as ${dbids[0]}_score from result r, chip_channel c, result_set rs where c.table_name=\"channel\" and c.table_id=${dbids[0]} and c.result_set_id=rs.result_set_id and rs.analysis_id=${ra_id} and c.chip_channel_id=r.chip_channel_id')\n";
		  $query .= "c2<-dbGetQuery(con, 'select r.probe_id, r.score as ${dbids[1]}_score from result r, chip_channel c, result_set rs where c.table_name=\"channel\" and c.table_id=${dbids[1]} and c.result_set_id=rs.result_set_id and rs.analysis_id=${ra_id} and c.chip_channel_id=r.chip_channel_id')\n";
	  
	  
	
		  #should do some sorting here?  Probes are in same order anyway
		  #does this affect how vsn works?  if not then don't bother and just load the correct probe_ids for each set
		  $query .= "raw_df<-cbind(c1[\"${dbids[0]}_score\"], c2[\"${dbids[1]}_score\"])\n";		
		  #variance stabilise
		  $query .= "vsn_df<-vsn(raw_df)\n";
    
	  
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
		  $query .= "glog_df<-cbind(rep(\"\", length(c1[\"probe_id\"])), c1[\"probe_id\"], sprintf(\"%.3f\", (exprs(vsn_df[,2]) - exprs(vsn_df[,1]))), rep(\"".$cc_id."\", length(c1[\"probe_id\"])))\n";
	  
	  
		  #load back into DB
		  #c3results<-cbind(rep("", length(c3["probe_id"])), c3["probe_id"], c3["c3_score"], rep(1, length(c3["probe_id"])), rep(1, length(c3["probe_id"])))
		  #may want to use safe.write here
		  #dbWriteTable(con, "result", c3results, append=TRUE)
		  #dbWriteTable returns true but does not load any data into table!!!
	  
		  $query .= "write.table(glog_df, file=\"${outfile}\", sep=\"\\t\", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)\n";
		
		  #tidy up here?? 
		}
      }
  
      $query .= "q();";
  
      open(RFILE, ">$R_file") || die("Cannot open $R_file for writing");
      print RFILE $query;
      close(RFILE);
 
      system($r_cmd) == 0 or throw("R $logic_name normalisation failed with error code $? ($R_file)");
    }
  }

  return;
}

#can we sub this? args: table_name, logic_name
#also use result_set_name
#would also clean all data for result set if recovery
#return would be result_set

sub get_import_ResultSet{
  my ($self, $table_name, $anal) = @_;
  
  if (!($anal && $anal->isa("Bio::EnsEMBL::Analysis") && $anal->dbID())) {
    throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
  }

  $self->log("Getting import ResultSet for analysis:\t".$anal->logic_name());

  my ($rset);
  my $result_adaptor = $self->db->get_ResultSetAdaptor();
  my $logic_name = ($anal->logic_name() eq "RawValue") ? "" : "_".$anal->logic_name();
  my $status = "IMPORTED${logic_name}";
  
  #could drop the table name here and use analysis hash?
  warn("Need to implement chip_sets here");
  
  foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {

    #clean chip import and generate rset
    if ($echip->has_status($status)) { #this translates to each channel have the IMPORTED_RawValue status
      $self->log("ExperimentalChip(".$echip->unique_id().") already has status:\t".$status);
    } else {

      $self->log("Found ExperimentalChip(".$echip->unique_id().") without status $status");

      if ( ! $rset) {
	
		if ($self->recovery()) {
		  warn ("Should use ResultSet name here, currently retrieving using the analysis and experiment id"); #experiment name?
		  #fetch by anal and experiment_id
		  #Need to change this to result_set.name!
		  warn("add chip set handling here");
	  
		  my @tmp = @{$result_adaptor->fetch_all_by_Experiment_Analysis($self->experiment(), $anal)};
		  throw("Found more than one ResultSet for Experiment:\t".$self->experiment->name()."\tAnalysis:\t$logic_name") if (scalar(@tmp) >1);
		  $rset = $tmp[0];

		  warn("Warning: Could not find recovery ResultSet for analysis ".$anal->logic_name()) if ! $rset;
		}
	
		if (! $rset) {
		  warn "Generating new ResultSet for analysis ".$anal->logic_name();
		  $self->log("Generating new ResultSet for analysis ".$anal->logic_name());
	  
		  $rset = Bio::EnsEMBL::Funcgen::ResultSet->new
			(
			 -analysis   => $anal,
			 -table_name => $table_name,
			);
	  
		  $result_adaptor->store($rset);
		}
      }
      
      if ($self->recovery()) {
	
		my (@cc_ids);
	
		if ($table_name eq 'experimental_chip') {
	  
		  push @cc_ids, $rset->get_chip_channel_id($echip->dbID()) if($rset->contains($echip));
		} else {
		  foreach my $chan (@{$echip->get_Channels()}) {
			push @cc_ids, $rset->get_chip_channel_id($chan->dbID()) if($rset->contains($chan));
		  }
		}

		$self->db->rollback_results(@cc_ids) if (@cc_ids);	
      }
    }
    
    #check whether it is present in the ResultSet and add if not
    if ($rset) {
      #ids will already be present if not rset i.e. already imported
      if ($table_name eq 'channel') {
	
		foreach my $chan (@{$echip->get_Channels()}) {
		  $rset->add_table_id($chan->dbID()) if(! $rset->contains($chan));
		}
      } else {					#experimental_chip
		$rset->add_table_id($echip->dbID()) if(! $rset->contains($echip));
      }
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


sub backup_file{
  my ($self, $file_path) = @_;

  throw("Must define a file path to backup") if(! $file_path);

  if (-f $file_path) {
    $self->log("Backing up:\t$file_path");
    system ("mv ${file_path} ${file_path}.".`date '+%T'`);
  }

  return;

}



1;
