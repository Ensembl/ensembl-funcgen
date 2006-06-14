
=head1 NAME

Bio::EnsEMBL::Funcgen::Experiment
  
=head1 SYNOPSIS




=back

=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Experiment;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_date);
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Funcgen::Helper Bio::EnsEMBL::Funcgen::ArrayDefs);

use Bio::EnsEMBL::Funcgen::ArrayDefs;#rename FormatDefs?
use Bio::EnsEMBL::Funcgen::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

#should also hold location of files aswell as methods

################################################################################

=head2 new

 Description : 
               

 Arg  [1]    : hash containing optional attributes :-
               
 ReturnType  : Experiment

 Example     : my $Exp = Bio::EnsEMBL::Experiment->new(
                                                      
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

	#change this to re-arrange?
    # objects private data and default values
    %attrdata = (
				 #User defined/built 
				 _instance    => undef,
				 _format      => undef,
				 _vendor      => undef,
				 _group       => undef,
				 _recover     => 0,
				 _data_dir    => $ENV{'EFG_DATA'},#?
				 _dump_fasta  => 0,
				 _norm_method => "vsn_norm",


				 #DBDefs, have ability to override here, or restrict to DBDefs.pm?
				 _pass       => undef,

				 #Need to separate pipeline vars/methods from true Experiment methods?
				 #Separate Pipeline object(or just control scripts? Handing step/dir validation?
				 _output_dir => undef,

				 #ArrayDefs defined
				 _input_dir  => undef,#Can pass this to over-ride ArrayDefs default?
				 _array_defs => undef,
				 _import_dir => undef,#parsed native data for import
				 _norm_dir   => undef,
								 

				 #Data defined
				 _group_dbid      => undef,
				 _experiment_id   => undef,
				 _echips          => {},
				 _array           => undef,
				 _achips          => undef,
				 _experiment_desc => "",

				 #Other
				 _db_adaptor    => undef,
				 _achip_reg     => {},


				 #HARDCODED
				 #Need to handle a lot more user defined info here which may not be caught by the data files
				 _design_type  => "binding_site_identification",#Hard coded MGED type term for now, should have option to enable other array techs?
				 _species      => "Homo sapiens",
	#description
	#Associated design types (separate table) timecourse etc?
	#epi_feature_name (separate table)
				);

    # set each class attribute using passed value or default value
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname} && defined $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }

	
	#Not all of these are required if using ArrayDefs?
	#Not required if just fetching the data

	foreach my $tmp("instance", "vendor", "format", "group", "data_dir"){
		$self->throw("Mandatory arg $tmp not been defined") if (! defined $self->{"_${tmp}"});
	}
	

	$self->set_defs();

	$self->debug(2,"Experiment class instance created.");
	$self->debug_hash(3, \$self);

    return ($self);
}


#Kept separate from new as it is not necessary to have native format raw data
#change name as need dir struc for processing aswell as import, may have imported in a different way
sub init_import{
	my ($self) = shift;

	#Check group and instance are valid

	if(! $self->get_group_dbid()){
		$self->throw("Your group(".$self->group().
					 ") is not valid, please use a valid group or register a new one in egroup table");
	}


	if( ! $self->recover() && $self->get_experiment_dbid()){
		$self->throw("Your experiment name is already registered in the database, please choose a different \"instance\", this will require renaming you input directory, or specify -recover if you are working with a failed import");
		#Done here to avoid creating incorrect directories
	}elsif($self->recover()){
		$self->throw("Need to implement recovery properly, clean/alter all data associated with this experiment dbid (experiment/experimental_chip/channel/variable/result/design_type/target) and validate array_chip/probe/probe_set/probe_feature to ensure that all were imported correctly");
	}


	#Should we separate path on group here too, so we can have a dev/test group?

	#Set and validate input dir
	$self->{'_input_dir'} = $self->get_def('input_dir') if(! defined $self->get_dir("input"));
	$self->throw("input_dir is not defined or does not exist") if(! -d $self->get_dir("input"));

	if(! defined $self->get_dir("output")){
		$self->{'_output_dir'} = $self->get_dir("data")."/".$self->vendor()."/".$self->instance();
		mkdir $self->get_dir("output") if(! -d $self->get_dir("output"));
	}
  
	$self->create_output_dirs("import", "norm");

	#remove and add specific report, this is catchig some Root stuff
	$self->log("Initiated efg import with following parameters:\n".Data::Dumper::Dumper(\$self));
	

	return;
}


sub create_output_dirs{
	my ($self, @dirnames) = @_;
	
	foreach my $name(@dirnames){
		$self->{"_${name}_dir"} = $self->get_dir("output")."/${name}" if(! defined $self->{"_${name}_dir"});
		mkdir $self->get_dir($name) if(! -d $self->get_dir($name));
	}

	return;
}

### GENERIC ACCESSOR METHODS ###

sub vendor{
	my ($self) = shift;	
	return $self->{'_vendor'};
}

sub instance{
	my ($self) = shift;	
	return $self->{'_instance'};
}

sub group{
	my ($self) = shift;	
	return $self->{'_group'};
}

sub recover{
	my $self = shift;
	return $self->{'_recover'};
}

sub experiment_desc{
	my $self = shift;
	return $self->{'_experiment_desc'};
}


sub format{
	my ($self) = shift;	
	return $self->{'_format'};
}


sub db_adaptor{
	my ($self) = @_;

	
	if(! defined $self->{'_db_adaptor'}){
		$self->{'_db_adaptor'} = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
																	   #have Group/DBDefs.pm override?
																	   pass => $self->pass(),
																	  );
	}

	return $self->{'_db_adaptor'};
}

sub pass{
	my $self = shift;
	return $self->{'_pass'};
}

#convinience wrapper methods, put in helper?
sub get_id{
	my ($self, $id_name) = @_;

	return $self->get_data("${id_name}_id");
}


sub species{
	my $self = shift;
	return $self->{'_species'};
}

sub get_dir{
	my ($self, $dirname) = @_;

	return $self->get_data("${dirname}_dir");
}

sub norm_method{
	my $self = shift;
	return $self->{'_norm_method'};
}


sub register_experiment{
	my ($self) = shift;

	#These should be vendor independent, only the read methods should need specific order?
	#Need to totally separate parse/read from import, so we can just do one method if required, i.e. normalise
	#Also need to move id generation to import methods, which would check that a previous import has been done, or check the DB for the relevant info?

	$self->init_import();
	$self->read_data("array");#rename this or next?
	$self->import_array();#Should only work if not present or if full recovery.  (Shouldn't really want to recover array, so separate into different script which checks all Experiment associated with array before removing from DB)
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

	foreach my $achip_id(keys %{$self->get_data('achips')}){
		#would need to get chip specific paths here e.g. probe_set.chip_design_id.txt

		#how are we going to validate the probe information?
		#Currently assuming data is correct and complete if we have a pre-registered array_chip

		if(! $self->is_registered($achip_id)){
			$self->db_adaptor->load_table_data("probe_set",  $self->get_dir("import")."/probe_set.txt");
			$self->db_adaptor->load_table_data("probe",  $self->get_dir("import")."/probe.txt");
			$self->db_adaptor->load_table_data("probe_feature",  $self->get_dir("import")."/probe_feature.txt");
		}else{
			$self->log("Probe(set & feature) data for array_chip already registered, not validated");
		}
	}

	return;
}

sub import_results{
	my ($self, $source_dir) = @_;

	$self->db_adaptor->load_table_data("result",  $self->get_dir($source_dir)."/result.txt");
	return;
}

sub read_data{
	my($self, $data_type) = @_;
	
	map {my $method = "read_${_}_data"; $self->$method()} @{$self->get_def("${data_type}_data")};

	return;
}



#Import methods go here as they should be vendor independent



#do we need DBDefs too?
#This would contain order of table field names for printing

#Hardcoded for nimblegen at the mo, would need to go through each individual chip if >1 design/DVD
sub import_array{
	my $self = shift;


	#IMPORT ARRAY
	#Setting size to one for nimblegen
	#coudl set alter value, currently 0 to recover?
	$self->db_adaptor->register_entry("array", $self->get_array_dbid(), 0, $self->array_data('design_name'), $self->format(), "1", $self->species(), $self->vendor(), $self->array_data('description'));

	#IMPORT ARRAY CHIPS
	#Currently only one for nimblegen
	foreach my $achip_id(keys %{$self->get_data('achips')}){
		#Need to record which ones are already present for validation and importing extra features
		$self->reg_achip($achip_id) if($self->get_achip_dbid($achip_id));	
		$self->db_adaptor->register_entry("array_chip", $self->get_achip_dbid($achip_id), 0, $achip_id, $self->get_array_dbid(),
										  $self->achip_data($achip_id, 'design_name'), $self->achip_data($achip_id, 'description'));				
	}

	return;
}

sub reg_achip{
	my ($self, $achip_id) = @_;
	$self->{'_achip_reg'}{$achip_id} = 1 if(defined $achip_id);
	return;
}

sub is_registered{
	my ($self, $achip_id) = @_;
	return (exists $self->{'_achip_reg'}{$achip_id}) ? $self->{'_achip_reg'}{$achip_id} : 0;
}

#primary
sub design_type{
	my $self = shift;
	return $self->{'_design_type'};
}

#imports experiment and chip info
sub import_experiment{
	my ($self) = shift;

	###EXPERIMENT
	#Pre-reg check already done in init_import, this accounts for recovery
	#Do we need to catch whether this returns 0 and clean up experiment table
	#add flag to alter table if any inconsistencies found?
	#change 0 to recover?
	$self->db_adaptor->register_entry("experiment", $self->get_experiment_dbid(), 0, $self->instance(), $self->get_group_dbid(),
									  get_date("date", $self->get_def("chip_file")), $self->design_type(), $self->experiment_desc());

	###CHIPS
	my ($alter, $chip_uid, $chan_uid);
	
	foreach $chip_uid(keys %{$self->get_data('echips')}){
		$alter = 0;
		
		#Get array_chip_id for given design_id
		#print $chipf "\t".$chip_uid."\t".$self->get_experiment_dbid()."\t".$self->get_achip_dbid($self->echip_data($chip_uid, 'design_id'))."\n";
		#Should we check for unique chip_unique_id here? We may have duplications between vendors
		#sanity check for now
		$self->throw("experimental_chip.chip_unique_id $chip_uid already exists, already imported or vendor id clash(implement vendor look up) Use recover?") if($self->get_echip_dbid($chip_uid));#&& ! $self->recover()
		
		$self->db_adaptor->register_entry("experimental_chip", $self->get_echip_dbid($chip_uid), $alter, $chip_uid, 
										  $self->get_experiment_dbid(), $self->get_achip_dbid($self->echip_data($chip_uid, 'design_id')));

		
		foreach $chan_uid(keys %{$self->get_channels($chip_uid)}){
			$alter = 0;
		
			#Validate/Alter/Delete channels and results
			if($self->get_channel_dbid($chan_uid)){
				
				if($self->recover()){
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

			$self->db_adaptor->register_entry("channel", $self->get_channel_dbid($chan_uid), $alter, $self->get_echip_dbid($chip_uid), 
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

sub get_group_dbid{
	my ($self) = @_;

	$self->{'_group_dbid'} ||= $self->db_adaptor->fetch_dbid_by_table_field("egroup", $self->group());
	return $self->{'_group_dbid'};
}

sub get_experiment_dbid{
	my ($self) = @_;

	$self->{'_experiment_dbid'} ||= $self->db_adaptor->fetch_dbid_by_table_field("experiment", $self->instance());
	return 	$self->{'_experiment_dbid'};
}

sub get_array_dbid{
	my ($self) = @_;
	
	#should this also query on vendor?
	$self->{'_array'}{'dbid'} ||= $self->db_adaptor->fetch_dbid_by_table_field("array", $self->array_data('design_name'));

	return 	$self->{'_array'}{'dbid'};
}


#We may need to add a vendor clause to these queries as it's not guaranteed that
#ids will be unique between vendors

sub get_achip_dbid{
	my ($self, $achip_id) = @_;

	${$self->get_data('achips', $achip_id)}{'dbid'} ||= $self->db_adaptor->fetch_dbid_by_table_field("array_chip", $achip_id, 'design_id');	
	return $self->achip_data($achip_id, 'dbid');
}

sub get_echip_dbid{
	my ($self, $echip_id) = @_;

	${$self->get_data('echips', $echip_id)}{'dbid'} ||= $self->db_adaptor->fetch_dbid_by_table_field("experimental_chip", $echip_id, 'chip_unique_id');	
	return $self->echip_data($echip_id, 'dbid');
}

#nimblegen specific!
sub get_channel_dbid{
	my ($self, $chan_uid) = @_;

	my ($chip_uid);

	#chan_uid is ${chip_uid}_${dye_freq}
	#This is not stored in the DB, so has to be retrieved 
	#This only works if each channel on a chip uses a different dye

	if( ! $self->channel_data($chan_uid, 'dbid')){
		($chip_uid = $chan_uid) =~ s/_.*//;
		$self->channel_data($chan_uid, 'dbid', $self->db_adaptor->fetch_channel_dbid_by_echip_dye($self->get_echip_dbid($chip_uid),
																								  ${$self->get_channel($chan_uid)}{'dye'}));
	}
	
	return $self->channel_data($chan_uid, 'dbid');
}



sub vsn_norm{
	my $self = shift;
	#This currently normalises a single two colour array at a time
	
	#HACK!
	my $metric_id = 2;#for glog tranformed data

	my @dbids;
	my $R_file = $self->get_dir("norm")."/norm.R";
	my $outfile = $self->get_dir("norm")."/result.txt";
	my $r_cmd = "R --no-save < $R_file >".$self->get_dir("norm")."/R.out 2>&1";

	unlink($outfile);#Need to do this as we're appending in the loop

	#setup qurey
	#scipen is to prevent probe_ids being converted to exponents
	my $query = "options(scipen=20);library(vsn);library(RMySQL);".
		  "con<-dbConnect(dbDriver(\"MySQL\"), dbname=\"efg_test\", user=\"ensadmin\", pass=\"ensembl\")\n";
	#currently having separate session fo

	foreach my $chip_id(keys %{$self->get_data("echips")}){		
		@dbids = ();

		foreach my $chan_id(keys %{$self->echip_data($chip_id, "channels")}){
			push @dbids, $self->get_channel_dbid($chan_id);
		}

		#should do some of this with maps?
		#HARDCODED metric ID for raw data as one
		$query .= "c1<-dbGetQuery(con, \"select probe_id, score as ${dbids[0]}_score from result where channel_id=${dbids[0]} and metric_id=1\")\n";
		$query .= "c2<-dbGetQuery(con, \"select probe_id, score as ${dbids[1]}_score from result where channel_id=${dbids[1]} and metric_id=1\")\n";

		#should do some sorting here?  Probes are in same order anyway
		#does this affect how vsn works?  if not then don't bother and just load the correct probe_ids for each set
		$query .= "raw_df<-cbind(c1[\"${dbids[0]}_score\"], c2[\"${dbids[1]}_score\"])\n";
		
		#glog tranform(vsn)
		$query .= "glog_df<-vsn(raw_df)\n";
		
		#should do someplot's of raw and glog and save here
		
		#set log func and params
		#$query .= "par(mfrow = c(1, 2)); log.na = function(x) log(ifelse(x > 0, x, NA));";
		#plot
		#$query .= "plot(exprs(glog_df), main = \"vsn\", pch = \".\");". 
		#  "plot(log.na(exprs(raw_df)), main = \"raw\", pch = \".\");"; 
		#FAILS ON RAW PLOT!!
		
		
		#par(mfrow = c(1, 2)) 
		#> meanSdPlot(nkid, ranks = TRUE) 
		#> meanSdPlot(nkid, ranks = FALSE) 
		
		#reform dataframe/matrix for each channel
		#3 sig dec places on scores
		$query .= "c1r<-cbind(rep(\"\", length(c1[\"probe_id\"])), c1[\"probe_id\"], format(exprs(glog_df[,1]), nsmall=3), rep(${metric_id}, length(c1[\"probe_id\"])), rep(${dbids[0]}, length(c1[\"probe_id\"])))\n";
		$query .= "c2r<-cbind(rep(\"\", length(c2[\"probe_id\"])), c2[\"probe_id\"], format(exprs(glog_df[,2]), nsmall=3), rep(${metric_id}, length(c2[\"probe_id\"])), rep(${dbids[1]}, length(c2[\"probe_id\"])))\n";
		
		#load back into DB
		#c3results<-cbind(rep("", length(c3["probe_id"])), c3["probe_id"], c3["c3_score"], rep(1, length(c3["probe_id"])), rep(1, length(c3["probe_id"])))
		#may want to use safe.write here
		#dbWriteTable(con, "result", c3results, append=TRUE)
		#dbWriteTable returns true but does not load any data into table!!!

		$query .= "write.table(c1r, file=\"${outfile}\", sep=\"\\t\", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)\n";
		$query .= "write.table(c2r, file=\"${outfile}\", sep=\"\\t\", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)\n";
		#tidy up here?? 
	}


	#or here, specified no save so no data will be dumped
	$query .= "q();";


	#This is giving duplicates for probe_ids 2 & 3 for metric_id =2 i.e. vsn'd data.
	#duplicates are not present in import file!!!!!!!!!!!!!!!!!!!!

	open(RFILE, ">$R_file") || die("Cannot open $R_file for writing");
	print RFILE $query;
	close(RFILE);
	
	system($r_cmd);

	return;
}

1;

