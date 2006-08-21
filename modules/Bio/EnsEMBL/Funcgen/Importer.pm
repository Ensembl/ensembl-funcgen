
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

    my ($self, %attrdata, $attrname, $argname);
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
				 recover     => 0,
				 location    => undef,
				 contact     => undef,
				 data_dir    => $ENV{'EFG_DATA'},#?
				 dump_fasta  => 0,
				 norm_method => "vsn_norm",
				 description => undef,
				 #DBDefs, have ability to override here, or restrict to DBDefs.pm?
				 pass       => undef,

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
				 channels        => {},

				 #Other
				 db    => undef,#this should really be an ExperimentAdaptor, ut it is the db at the moment
				 #check for ~/.ensembl_init to mirror general EnsEMBL behaviour
				 reg_config    => (-f "$ENV{'HOME'}/.ensembl_init") ? "$ENV{'HOME'}/.ensembl_init" : undef,
			
		
				 #HARDCODED
				 #Need to handle a lot more user defined info here which may not be caught by the data files
				 design_type  => "binding_site_identification",#Hard coded MGED type term for now, should have option to enable other array techs?

				 #Need to add all aliases in Importer?Just check get_alias and reset to Registry defined standard
				 #This is used to retrieve the efg and dnadb adaptors
				 #Need to validate this somehow?
				 species      => "Human",
				 #description
				 #Associated design types (separate table) timecourse etc?
				 #epi_feature_name (separate table)
				);

    # set each class attribute using passed value or default value
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname} && defined $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }


	#Set vars and defaults
	#$self->



	
	#Not all of these are required if using ArrayDefs?
	#Not required if just fetching the data

	foreach my $tmp("name", "vendor", "format", "group", "data_dir"){
		$self->throw("Mandatory arg $tmp not been defined") if (! defined $self->{$tmp});
	}

	#Set vendor specific vars/methods
	$self->set_defs();


	#load registry
	if(! defined $self->{'_reg_config'} && ! %Bio::EnsEMBL::Registry::registry_register){
		#current ensembl DBs
		$reg->load_registry_from_db(
							   -host => "ensembldb.ensembl.org",
							   -user => "anonymous",
							   #-verbose => 1,
							  );
	}else{#from config
		$reg->load_all($self->{'_reg_config'}, 1);
	}

	#add EFG DB to reg here!! Will be done automatically when incoporated into Regsitry
	#Merely calling the DBAdaptor configures the registry


	$self->debug(2, "Importer class instance created.");
	$self->debug_hash(3, \$self);

    return ($self);
}


#Kept separate from new as it is not necessary to have native format raw data
#change name as need dir struc for processing aswell as import, may have imported in a different way
sub init_import{
	my ($self) = shift;


	#Need to import to egroup here if not present and name, location & contact specified
	$self->validate_group();



	
	#fetch experiment
	#if recovery and ! experiment throw 
	#else if ! experiment new and store

	#rename instance to name? Do we need this composite fetch?
	my $exp_adaptor = $self->db->get_ExperimentAdaptor();
	my $exp = $exp_adaptor->fetch_by_name($self->name());#, $self->group());

	#should we not just do store here, as this will return the experiment if it has already been stored?



	if ($self->recovery() && (! $exp)){
		throw("No previously stored experiment defined with recovery mode, Need to implement recovery mode"); 
	}
	elsif((! $self->recovery()) && $exp){
		throw("Your experiment name is already registered in the database, please choose a different \"name\", this will require renaming you input directory, or specify -recover if you are working with a failed import. Or specify recovery?");	
	}
	else{
		$exp = Bio::EnsEMBL::Funcgen::Experiment->new(
													  -GROUP => $self->group(),
													  -NAME  => $self->name(),
													  -DATE  => &get_date("date", $self->get_def("chip_file")),
													  -PRIMARY_DESIGN_TYPE => $self->design_type(),
													  -DESCRIPTION => $self->description(),
													  -ADAPTOR => $self->db->get_ExperimentAdaptor(),
													 );

		($exp) =  @{$exp_adaptor->store($exp)};		
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
		$self->{'vendor'} = shift;
	}

	return $self->{'name'};
}

sub group{
	my ($self) = shift;	

	if(@_){
		$self->{'group'} = shift;
	}

	return $self->{'group'};
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

	if(@_ && $_[0]->isa($db_adaptor)){
		$self->{'db'} = shift;
	}else{
		throw("Need to pass a valid $db_adaptor") if(! $_[0]->isa($db_adaptor));
	}

	return $self->{'db'};
}

sub pass{
	my $self = shift;

	if(@_){
		$self->{'pass'} = shift;
	}

	return $self->{'pass'};
}

sub dump_fasta{
	my $self = shift;

	if(@_){
		$self->{'dump_fasta'} = shift;
	}

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

	if(@_){
		$self->{'species'} = shift;
	}

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
	if( ! @_ || ! $_[0]->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")){
		throw("You need to pass a valid dnadb adaptor to register the experiment");
	}
	else{
		$self->db->dnadb($_[0]);
	}
		




	#These should be vendor independent, only the read methods should need specific order?
	#Need to totally separate parse/read from import, so we can just do one method if required, i.e. normalise
	#Also need to move id generation to import methods, which would check that a previous import has been done, or check the DB for the relevant info?

	$self->init_import();
	$self->read_data("array");#rename this or next?
	$self->import_experiment();#imports experiment and chip data
	$self->read_data("probe");
	$self->read_data("results");
	$self->import_probe_data();
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
	
	


	foreach my $achip_id(keys %{$self->get_data('achips')}){
		#would need to get chip specific paths here e.g. probe_set.chip_design_id.txt
		#how are we going to validate the probe information?
		#Currently assuming data is correct and complete if we have a pre-registered array_chip

		

		if(! $self->arrays()->[0]->get_achip_status($achip_id, "IMPORTED")){
			$self->db->load_table_data("oligo_probe_set",  $self->get_dir("import")."/probe_set.txt");
			$self->db->load_table_data("oligo_probe",  $self->get_dir("import")."/probe.txt");
			$self->db->load_table_data("oligo_feature",  $self->get_dir("import")."/probe_feature.txt");

			#error catch here!! DO not want to set imported if we've not imported properly
			$self->db->set_status('array_chip', $self->achip_data($achip_id, 'dbID'), "IMPORTED");
			
			
		}else{
			$self->log("Probe(set & feature) data for array_chip already IMPORTED");
		}
	}

	return;
}

sub import_results{
	my ($self, $source_dir) = @_;

	$self->db->load_table_data("result",  $self->get_dir($source_dir)."/result.txt");
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

#imports experiment and chip info
sub import_experiment{
	my ($self) = shift;

	#This is all now done in read methods

	return;

	###CHIPS
	my ($alter, $chip_uid, $chan_uid);
	
	foreach $chip_uid(keys %{$self->get_data('echips')}){
		$alter = 0;
		
		#Get array_chip_id for given design_id
		#print $chipf "\t".$chip_uid."\t".$self->get_experiment_dbid()."\t".$self->get_achip_dbid($self->echip_data($chip_uid, 'design_id'))."\n";
		#Should we check for unique chip_unique_id here? We may have duplications between vendors

		
		foreach $chan_uid(keys %{$self->get_channels($chip_uid)}){
			$alter = 0;
		
			#Validate/Alter/Delete channels and results
			if($self->get_channel_dbid($chan_uid)){
				
				if($self->recovery()){
					#set alter/overwrite flag here for register_entry
					#delete/validate associated results
					$alter = 1;
					
					#should we set a flag to validate last result vs. last result in new results?
					#$self->db_adaptor->delete_results_by_channel($self->get_channel_dbid($chan_uid));
				}else{
					#Do we need this as we should have thrown already in int with the experiment_dbid
					#admin may have delete an experiment and not it's associated chips/results etc
				}
			}

			$self->db->register_entry("channel", $self->get_channel_dbid($chan_uid), $alter, $self->get_echip($chip_uid)->dbid(), 
											  $self->channel_data($chan_uid, 'sample_label'), $self->channel_data($chan_uid, 'dye'), 
											  $self->channel_data($chan_uid, 'type'));

			#Would do other experimental_variants here, inc sample species
			#how are we going to do this from text files?
			#manually define and write script to update chips, this could be predefined and warn if no extra exp_vars defined
			#This would be failry hefty manual work for large arrays :(	
		}
	}
	return;
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

	return $self->db->get_DNADBSliceAdaptor->fetch_by_region("chromosome", $chr, $start, $end)->get_seq_region_id();

}

1;

