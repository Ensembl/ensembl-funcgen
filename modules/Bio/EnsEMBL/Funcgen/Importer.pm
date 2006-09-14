
=head1 NAME

Bio::EnsEMBL::Funcgen::Importer
  
=head1 SYNOPSIS





=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Importer;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_date);
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Funcgen::Helper Bio::EnsEMBL::Funcgen::ArrayDefs);

use Bio::EnsEMBL::Funcgen::Experiment;#
use Bio::EnsEMBL::Funcgen::ArrayDefs;#rename FormatDefs?
use Bio::EnsEMBL::Funcgen::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;#eventually add this to Registry
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::ConfigRegistry;
my $reg = "Bio::EnsEMBL::Registry";



################################################################################

=head2 new

 Description : 
               

 Arg  [1]    : hash containing optional attributes :-
               
 ReturnType  : Experiment

 Example     : my $Exp = Bio::EnsEMBL::Importer->new(
                                                      
                                                     );

 Exceptions  : 

=cut

################################################################################

sub new{
    my ($caller, %args) = @_;

    my ($self, %attrdata, $attrname, $argname, $db);
    my $reg = "Bio::EnsEMBL::Registry";

    my $class = ref($caller) || $caller;

	#Create object from parent class
	$self = $class->SUPER::new(%args);

    # objects private data and default values
    %attrdata = (
				 #User defined/built 
				 name        => undef,
				 format      => 'Tiled',
				 vendor      => undef,
				 group       => undef,
				 species     => undef,
				 data_version => undef,
				 recover     => 0,
				 location    => undef,
				 contact     => undef,
				 data_dir    => $ENV{'EFG_DATA'},#?
				 dump_fasta  => 0,
				 norm_method => "vsn_norm",
				 description => undef,
				 #DBDefs, have ability to override here, or restrict to DBDefs.pm?
				 pass       => undef,
				 host       => undef,
				 user       => undef,
				 port       => undef,


		 #vars to handle array chip sets
		 #no methods for these as we're replacing with a meta file or something
		 array_set => 0,
		 array_name => undef,
				 

				 #Need to separate pipeline vars/methods from true Experiment methods?
				 #Separate Pipeline object(or just control scripts? Handing step/dir validation?
				 output_dir => undef,

				 #ArrayDefs defined
				 input_dir  => undef,#Can pass this to over-ride ArrayDefs default?
				 array_defs => undef,
				 import_dir => undef,#parsed native data for import
				 norm_dir   => undef,
								 

				 #Data defined
				 #_group_dbid      => undef,
				 #_experiment_id   => undef,
				 echips          => {},
				 arrays          => [],
				 achips          => undef,
				 channels        => {},#?

		 #Other
		 db    => undef,#this should really be an ExperimentAdaptor, but it is the db at the moment?
		 dbname => undef,#to over-ride autogeneration of eFG dbname
				 #check for ~/.ensembl_init to mirror general EnsEMBL behaviour
				 reg_config    => (-f "$ENV{'HOME'}/.ensembl_init") ? "$ENV{'HOME'}/.ensembl_init" : undef,
			
		
				 #HARDCODED
				 #Need to handle a lot more user defined info here which may not be caught by the data files
				 design_type  => "binding_site_identification",#Hard coded MGED type term for now, should have option to enable other array techs?
				);

    # set each class attribute using passed value or default value
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname} && defined $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }


	
	#Can some of these be set in ArrayDefs or "Vendor"Defs?
	#pass?

	foreach my $tmp("name", "vendor", "format", "group", "data_dir", "data_version", "species", "host", "user"){
		$self->throw("Mandatory arg $tmp not been defined") if (! defined $self->{$tmp});
	}

	#Set vendor specific vars/methods
	$self->set_defs();

	
	### LOAD AND RE-CONFIG REGISTRY ###
	if(! defined $self->{'_reg_config'} && ! %Bio::EnsEMBL::Registry::registry_register){
	
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
	
		if(! $self->db() || ($self->data_version() ne $self->db->_get_schema_build())){
		  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
							    -host => 'ensembldb.ensembl.org',
							    -user => 'anonymous',
							    -dbname => $self->species()."_core_".$self->data_version(),
							    -species => $self->species(),
							   );
		}else{
		  $db = $self->db->dnadb();
		}


		$self->{'dbname'} ||= $self->species()."_core_".$self->data_version();

		#generate and register DB with local connection settings
		$db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
								   -user => $self->user(),
								   -host => $self->host(),
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

	}else{#from config
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

sub init_import{
  my ($self) = shift;

  
  #Need to import to egroup here if not present and name, location & contact specified
  $self->validate_group();



	
  #fetch experiment
  #if recovery and ! experiment throw 
  #else if ! experiment new and store
  
  #rename instance to name? Do we need this composite fetch?
  my $exp_adaptor = $self->db->get_ExperimentAdaptor();
  
  #print "XXXXXXXXXX featching by name ".$self->name()."\n";

  my $exp = $exp_adaptor->fetch_by_name($self->name());#, $self->group());
  #should we not just do store here, as this will return the experiment if it has already been stored?

  if ($self->recovery() && (! $exp)){
    warn("No previously stored experiment defined with recovery mode, Importing as normal"); 
  }

  if((! $self->recovery()) && $exp){
    throw("Your experiment name is already registered in the database, please choose a different \"name\", this will require renaming you input directory, or specify -recover if you are working with a failed import. Or specify recovery?");
    #can we skip this and store, and then check in register experiment if it is already stored then throw if not recovery
  }
  else{#niether or both?? or recover and exp
 
    
    $exp = Bio::EnsEMBL::Funcgen::Experiment->new(
						  -GROUP => $self->group(),
						  -NAME  => $self->name(),
						  -DATE  => &get_date("date", $self->get_def("chip_file")),
						  -PRIMARY_DESIGN_TYPE => $self->design_type(),
						  -DESCRIPTION => $self->description(),
						  -ADAPTOR => $self->db->get_ExperimentAdaptor(),
						 );
    
    ($exp) =  @{$exp_adaptor->store($exp)};	#skip this bit?	
  }

  
  $self->experiment($exp);
  
  #Should we separate path on group here too, so we can have a dev/test group?
  
  #Set and validate input dir
  $self->{'input_dir'} = $self->get_def('input_dir') if(! defined $self->get_dir("input"));
  $self->throw("input_dir is not defined or does not exist") if(! -d $self->get_dir("input"));#Helper would fail first on log/debug files
  
  if(! defined $self->get_dir("output")){
    $self->{'output_dir'} = $self->get_dir("data")."/".$self->vendor()."/".$self->name();
    mkdir $self->get_dir("output") if(! -d $self->get_dir("output"));
  }
  
  $self->create_output_dirs("import", "norm");
  
  #remove and add specific report, this is catchig some Root stuff
  $self->log("Initiated efg import with following parameters:\n".Data::Dumper::Dumper(\$self));
  
  return;
}




sub validate_group{
	my ($self) = shift;

	my $group_ref = $self->db->fetch_group_details($self->group());

	if (! $group_ref){
		if($self->location() && $self->contact()){
			$self->db->import_group($self->group(), $self->location, $self->contact());
		}else{
			throw("Group ".$self->group()." does not exist, please specify a location and contact to register the group");
		}
	}

	return;
}

sub create_output_dirs{
	my ($self, @dirnames) = @_;
	
	foreach my $name(@dirnames){
		$self->{"${name}_dir"} = $self->get_dir("output")."/${name}" if(! defined $self->{"${name}_dir"});
		mkdir $self->get_dir($name) if(! -d $self->get_dir($name));
	}

	return;
}

### GENERIC ACCESSOR METHODS ###

sub vendor{
	my ($self) = shift;

	if(@_){
		$self->{'vendor'} = shift;
	}
	
	return $self->{'vendor'};
}

sub location{
	my ($self) = shift;

	if(@_){
		$self->{'location'} = shift;
	}
	
	return $self->{'location'};
}

sub contact{
	my ($self) = shift;

	if(@_){
		$self->{'contact'} = shift;
	}
	
	return $self->{'contact'};
}



sub name{
	my ($self) = shift;	

	if(@_){
		$self->{'name'} = shift;
	}

	return $self->{'name'};
}

sub verbose{
	my ($self) = shift;	

	if(@_){
		$self->{'verbose'} = shift;
	}

	return $self->{'verbose'};
}

sub data_version{
	my ($self) = shift;	

	if(@_){
		$self->{'data_version'} = shift;
		#have reset_dnadb here?
		#Can only do this if we set data_version directly in new
		#rather than calling this method
		#as reset_dnadb assumes db is set
	}

	return $self->{'data_version'};
}




sub group{
	my ($self) = shift;	

	if(@_){
		$self->{'group'} = shift;
	}

	return $self->{'group'};
}

sub dbname{
  my ($self) = shift;	
  
  if(@_){
    $self->{'dbname'} = shift;
  }
  
  return $self->{'dbname'};
}

sub recovery{
	my $self = shift;

	if(@_){
		$self->{'recover'} = shift;
	}

	return $self->{'recover'};
}


sub description{
	my $self = shift;

	if(@_){
		$self->{'description'} = shift;
	}

	return $self->{'description'};
}


sub format{
	my ($self) = shift;	

	if(@_){
		$self->{'format'} = shift;
	}

	return $self->{'format'};
}



sub experiment{
	my ($self) = shift;	

	if(@_){
		$self->{'experiment'} = shift;
	}

	return $self->{'experiment'};
}


sub db{
	my $self = shift;

	my $db_adaptor = "Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor";

	if(defined $_[0] && $_[0]->isa($db_adaptor)){
		$self->{'db'} = shift;
	}elsif(defined $_[0]){
		throw("Need to pass a valid $db_adaptor");
	}

	return $self->{'db'};
}

sub pass{
	my $self = shift;

	$self->{'pass'} = shift if(@_);
	return $self->{'pass'};
}

sub host{
	my $self = shift;

	$self->{'host'} = shift if(@_);
	return $self->{'host'};
}

sub port{
	my $self = shift;

	$self->{'port'} = shift if(@_);
	return $self->{'port'};
}

sub user{
	my $self = shift;

	$self->{'user'} = shift if(@_);
	return $self->{'user'};
}




sub dump_fasta{
	my $self = shift;

	$self->{'dump_fasta'} = shift if(@_);

	return $self->{'dump_fasta'};
}


#convinience wrapper methods, put in helper?
sub get_id{
	my ($self, $id_name) = @_;

	return $self->get_data("${id_name}_id");
}


#used for generating dbadaptors
sub species{
	my $self = shift;

	$self->{'species'} = shift if(@_);
	
	return $self->{'species'};
}

sub get_dir{
	my ($self, $dirname) = @_;

	return $self->get_data("${dirname}_dir");
}

sub norm_method{
	my $self = shift;

	if(@_){
		$self->{'norm_method'} = shift;
	}

	return $self->{'norm_method'};
}


sub register_experiment{
	my ($self) = shift;

	#Need to check for dnadb passed with adaptor to contructor
	if(@_){ 
		if( ! $_[0]->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
			throw("You need to pass a valid dnadb adaptor to register the experiment");
		}
		$self->db->dnadb($_[0]);
	}
	elsif( ! $self->db()){
		throw("You need to pass/set a DBAdaptor with a DNADB attached of the relevant data version");
	}

	#This could still be the default core db for the current version
	#warn here if not passed DB?




	#These should be vendor independent, only the read methods should need specific order?
	#Need to totally separate parse/read from import, so we can just do one method if required, i.e. normalise
	#Also need to move id generation to import methods, which would check that a previous import has been done, or check the DB for the relevant info?

	$self->init_import();

	#check here is exp already stored?  Will this work properly?
	#then throw if not recovery
	#else do following, but change to _private type methods
	#as bypassing this register_experiment method and calling directly will cause problems with duplication of data
	#we need to make this more stringent, maybe we can do a caller in each method to make sure it is register experiment calling each method

	$self->read_data("array");#rename this or next?

	#$self->import_experiment();#imports experiment and chip data
	$self->read_data("probe");
	$self->read_data("results");
	#$self->import_probe_data();
	$self->import_results("import");

	#Need to be able to run this separately, so we can normalise previously imported sets with different methods
	#should be able t do this without raw data files e.g. retrieve info from DB
	my $norm_method = $self->norm_method();
	$self->$norm_method;
	$self->import_results("norm");

	return;
}


#the generic read_methods should go in here too?
#should reorganise these emthods to split reading the array data, and the actual data
#currently:
#meta reads array and chip data
#probe reads probe_set, probes, which should definitely be in array, probe_feature? and results
#native data format may not map to these methods directly, so may need to call previous method if required data not defined

sub import_probe_data{
  my $self = shift;
  
  #This has been genericised for multiple chips, currently only one/data set with Nimblegen
  #would need to write different files for each chip and store the file names
  
  #Need to cycle through array chip getting status of achip e.g. REGISTERED, IMPORTED?
  #This is a array:array_chip i.e. 1 to 1 relationship at the moment
  #Cannot use array->array_chips, as this may bring back more than is valid for this experiment
  
  
  
  #foreach my $achip_id(keys %{$self->get_data('achips')}){
  foreach my $design_id(@{$self->arrays->[0]->get_design_ids()}){
    
    #would need to get chip specific paths here e.g. probe_set.chip_design_id.txt
    #how are we going to validate the probe information?
    #Currently assuming data is correct and complete if we have a pre-registered array_chip
    
    if(! $self->arrays()->[0]->get_achip_status($design_id, "IMPORTED")){
      $self->db->load_table_data("oligo_probe_set",  $self->get_dir("import")."/probe_set.txt");
      $self->db->load_table_data("oligo_probe",  $self->get_dir("import")."/probe.txt");
      $self->db->load_table_data("oligo_feature",  $self->get_dir("import")."/probe_feature.txt");
      
      #error catch here!! DO not want to set imported if we've not imported properly
      $self->db->set_status('array_chip', $self->achip_data($design_id, 'dbID'), "IMPORTED");
    }else{
      $self->log("Probe(set & feature) data for array_chip already IMPORTED");
    }
  }
  
  return;
}

sub import_results{
  my ($self, $source_dir) = @_;
  
  foreach my $array(@{$self->arrays()}){

    foreach my $design_id(@{$array->get_design_ids()}){
      my %ac = %{$array->get_array_chip_by_design_id($design_id)};
      print "Loading reults for ".$ac{'name'}."\n";
      $self->db->load_table_data("result",  $self->get_dir($source_dir)."/result.".$ac{'name'}.".txt");
    }
  }

  return;
}

sub read_data{
	my($self, $data_type) = @_;
	
	map {my $method = "read_${_}_data"; $self->$method()} @{$self->get_def("${data_type}_data")};

	return;
}




sub design_type{
	my $self = shift;
	return $self->{'design_type'};
}





#These are direct references to table names, which should be separated buried in the DBAdaptor, via DBDefs.pm?
#Some of these may also be different for different vendors > ArrayDefs?
#All dbids should be retrieved this way, as they are lazy loaded




#nimblegen specific!
sub get_channel_dbid{
	my ($self, $chan_uid) = @_;

	my ($chip_uid);

	#chan_uid is ${chip_uid}_${dye_freq}
	#This is not stored in the DB, so has to be retrieved 
	#This only works if each channel on a chip uses a different dye

	if( ! $self->channel_data($chan_uid, 'dbID')){
		($chip_uid = $chan_uid) =~ s/_.*//;
		$self->channel_data($chan_uid, 'dbid', $self->db->fetch_channel_dbid_by_echip_dye($self->get_echip($chip_uid)->dbID(),
																								  ${$self->get_channel($chan_uid)}{'dye'}));
	}
	
	return $self->channel_data($chan_uid, 'dbid');
}







#could we use the seq region cache instead?
#this seems like a lot of overhead for getting the id
sub get_chr_seq_region_id{
	my ($self, $chr, $start, $end) = @_;
	#what about strand info?

	#use start and stop to prevent problems with scaffodl assemblies, i.e. >1 seq_region_id
	#my $slice = $self->slice_adaptor->fetch_by_region("chromosome", $chr, $start, $end);
	#we could pass the slice back to the slice adaptor for this, to avoid dbid problems betwen DBs

	return $self->db->get_SliceAdaptor->fetch_by_region("chromosome", $chr, $start, $end)->get_seq_region_id();
}


#We need to enable Importer for previous imports
#currently dies if experiment name already registered
#This will enable re-imports of raw data via different analyses/normalisations
#Have Norm class or contain methods in importer?
#Need to have analysis set up script for all standard analyses.

sub vsn_norm{
  my $self = shift;
  #This currently normalises a single two colour array at a time
	
  my @dbids;
  my $aa = $self->db->get_AnalysisAdaptor();
  my $ra_id = $aa->fetch_by_logic_name("RawValue")->dbID();
  my $va_id = $aa->fetch_by_logic_name("VSN_GLOG")->dbID();
  my $R_file = $self->get_dir("norm")."/norm.R";
  my $outfile = $self->get_dir("norm")."/result.txt";
  my $r_cmd = "R --no-save < $R_file >".$self->get_dir("norm")."/R.out 2>&1";
  
  unlink($outfile);#Need to do this as we're appending in the loop
  
  #setup qurey
  #scipen is to prevent probe_ids being converted to exponents
  my $query = "options(scipen=20);library(vsn);library(RMySQL);".
    "con<-dbConnect(dbDriver(\"MySQL\"), dbname=\"".$self->db->dbc->dbname()."\", user=\"".$self->user()."\"";
  $query .= (defined $self->pass()) ? ", pass=\"".$self->pass()."\")\n" : ")\n";
  
  #currently having separate session fo
  
  
  
  #This is now retrieving a ExperimentalChip obj
  
  foreach my $echip(values %{$self->get_data("echips")}){
    print "Performing VSN for ".$echip->unique_id()."\n";
    @dbids = ();

    foreach my $chan(@{$echip->get_Channels()}){
      
      if($chan->type() eq "EXPERIMENTAL"){
	push @dbids, $chan->dbID();
      }else{
	unshift @dbids, $chan->dbID();
      }
    }
    
    
    throw("vsn_norm does not accomodate more than 2 channels") if scalar(@dbids > 2);
    
    #should do some of this with maps?
    #HARDCODED metric ID for raw data as one
    
    #Need to get total and experimental here and set db_id accordingly
    
    
    $query .= "c1<-dbGetQuery(con, 'select oligo_probe_id, score as ${dbids[0]}_score from result where table_name=\"channel\" and table_id=${dbids[0]} and analysis_id=${ra_id}')\n";
    $query .= "c2<-dbGetQuery(con, 'select oligo_probe_id, score as ${dbids[1]}_score from result where table_name=\"channel\" and table_id=${dbids[1]} and analysis_id=${ra_id}')\n";
    
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
    $query .= "glog_df<-cbind(rep(\"\", length(c1[\"oligo_probe_id\"])), c1[\"oligo_probe_id\"], format(exprs(vsn_df[,2]) - exprs(vsn_df[,1]), nsmall=3), rep(\"${va_id}\", length(c1[\"oligo_probe_id\"])), rep(\"".$echip->dbID()."\", length(c1[\"oligo_probe_id\"])),   rep(\"experimental_chip\", length(c1[\"oligo_probe_id\"])))\n";
    
    
    #load back into DB
    #c3results<-cbind(rep("", length(c3["probe_id"])), c3["probe_id"], c3["c3_score"], rep(1, length(c3["probe_id"])), rep(1, length(c3["probe_id"])))
    #may want to use safe.write here
    #dbWriteTable(con, "result", c3results, append=TRUE)
    #dbWriteTable returns true but does not load any data into table!!!
    
    $query .= "write.table(glog_df, file=\"${outfile}\", sep=\"\\t\", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)\n";
    
    warn("Need to implement R DB import\n");
    
    #tidy up here?? 
  }
  
  
  #or here, specified no save so no data will be dumped
  $query .= "q();";
  
  
  #This is giving duplicates for probe_ids 2 & 3 for metric_id =2 i.e. vsn'd data.
  #duplicates are not present in import file!!!!!!!!!!!!!!!!!!!!
  
  open(RFILE, ">$R_file") || die("Cannot open $R_file for writing");
  print RFILE $query;
  close(RFILE);
  
  system($r_cmd) == 0 or throw("R normalisation failed with error code $? ($R_file)");
  
  return;
}







1;

