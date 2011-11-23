
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

Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  
=head1 SYNOPSIS

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => "ensembldb.ensembl.org",
   -dbname => "mus_musculus_funcgen_41_36b",
   -species => "Mus_musculus",
   -user => "anonymous",
   -dnadb => $mouse_core_db,
   -port => '3307',
  );

my $experiment_adaptor = $db->get_ExperimentAdaptor();

=back

=head1 DESCRIPTION

This is a wrapper method for Bio::EnsEMBL::DBAdaptor, providing Funcgen
specific methods.

=cut

################################################################################

#To do
#1 Remove need for dnadb to be set, so we can do non feature imports/queries without setting the dnadb?

package Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);# Bio::EnsEMBL::Funcgen::Helper);

use DBI;

use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
my $reg = "Bio::EnsEMBL::Registry";


=head2 new

  Arg [-DNADB]: (optional) Bio::EnsEMBL::DBSQL::DBAdaptor DNADB 
               DNADB will be automatically select using the given parameters,
               the current registry dnadb host or ensembldb.
  Arg [-NO_CACHE]: (optional) int 1
               This option will turn off caching for slice features, so, 
               every time a set of features is retrieved, they will come from
               the database instead of the cache. This option is only recommended
               for advanced users, specially if you need to store and retrieve
               features. It might reduce performance when querying the database if 
               not used properly. If in doubt, do not use it or ask in ensembl-dev               
  Arg [..]   : Other args are passed to superclass
               Bio::EnsEMBL::DBSQL::DBConnection
  Example    : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql' );
  Exmaple2   : $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                    -species => 'Homo_sapiens',
                                                    -group   => 'core'
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql');
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $class, @args ) = @_;

  #Can't do this here yet due to auto-dnadb setting
  #$group ||= 'funcgen';#could pass production here, but would also require species to be passed and is_multi!?

  #Force group to be funcgen as this is the only valid group.
  my $self = $class->SUPER::new(@args, '-group', 'funcgen');
 
  #Currently only uses dnadb params to auto select
  #If attached dnadb or default reg dnadb does not match given assembly version.

  if($self->species eq 'DEFAULT'){  #Auto set species if not set
	
	#Can't do list_value_by_key as this depends on species, so we get a circular reference
	#This has already been set in the registry in SUPER::new above as DEFAULT!
	#So we need to reset this in the registry here?

	$self->{'_species'} = ${$self->get_MetaContainer->list_value_by_key('species.production_name')}[0] || 'DEFAULT';

	if($self->species ne 'DEFAULT'){ #Reset this in the registry to the correct species
	  $self = Bio::EnsEMBL::Utils::ConfigRegistry::gen_load($self);
	  #This causes duplicate software vs DB release warnings
	}
  }
  #else should we redefine the species as the standard alias if it does not match?
  #This would prevent external_db species name testing
  #Maybe the solution is to make external_db multi_species

  my ( $dnadb_host, $dnadb_user, $dnadb_port, $dnadb_pass, $dnadb_assm)
    = rearrange( [ 'DNADB_HOST', 'DNADB_USER',
                   'DNADB_PORT', 'DNADB_PASS',
				   'DNADB_ASSEMBLY',
                 ],
                 @args );

  my $default_dnadb = $self->SUPER::dnadb;
  my ($default_host, $default_port, $default_user, $default_pass, $default_assm, $efg_assm, $dnadb_defined);


  if($default_dnadb->group eq 'core'){
	#This means you have loaded a registry or pass a dnadb to the efg DBAdaptor
	
	$default_host = $default_dnadb->dbc->host;
	$default_port = $default_dnadb->dbc->port;
	$default_user = $default_dnadb->dbc->username;
	$default_pass = $default_dnadb->dbc->password;
	($default_assm = (split/_/, $self->_get_schema_build($default_dnadb))[1]) =~ s/[a-z]//;
	$dnadb_defined = 1;
  }


  #Defaults now set in dnadb as we want to test registry first;
  #These will pick default_dnadb values if params not explicitly set
  $self->{'dnadb_host'} = $dnadb_host || $default_host || 'ensembldb.ensembl.org';

  #Do we need this to be dnadb_ports 3306 5306 for different mysqls?
  #This is not correct for ensembldb, but we over-ride this in set_dnadb_by_assembly_version.
  $self->{'dnadb_port'} = $dnadb_port || $default_port || 3306;
  $self->{'dnadb_user'} = $dnadb_user || $default_user || 'anonymous';
  $self->{'dnadb_pass'} = $dnadb_pass || $default_pass || undef;
  ($efg_assm = (split/_/, $self->_get_schema_build($self))[1]) =~ s/[a-z]//;
  $dnadb_assm ||= $default_assm || $efg_assm;
  $self->{'dnadb_assm'} = $dnadb_assm;


  #Now we only want to reset the dnadb if it does not match the dnadb_assm
  #Use dnadb method here as this will either return a predefined dnadb(attached or reg) or auto select
  #Can we change this so that it only does this when we call dnadb?
  #This resulted in circular reference, so we need to be careful about changing this

  if($self->_get_schema_build($self->dnadb()) !~ /[0-9]+_${dnadb_assm}[a-z]*$/){
	#Do we need to consider reg_config here?
	#We could potentially have two version of the core DB in the config
	#But we would expect the user to handle predefining the dnadb correctly in this case
	warn ':: WARNING: Unable to match assembly version between the dnadb name ('.
	  $self->dnadb->dbc->dbname.') and the specified -dnadb_assm '.$self->dnadb_assembly.
		"\nMaybe you need to rename your DBs according to the Ensembl naming convention e.g. myprefix_homo_sapiens_55_37";

	if($dnadb_defined && $dnadb_host){
	  warn ":: Over-riding pre-defined dnadb host values(reg/-dnadb arg) with dnadb params:\t".
		$self->dnadb_user.'@'.$self->dnadb_host.':'.$self->dnadb_port;
	}
	$self->set_dnadb_by_assembly_version($self->dnadb_assembly);
  }

  return $self;
} 


#Move these stored method to BaseAdaptor?

=head2 is_stored_and_valid

  Arg [1]    : string - class namespace
  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable e.g. ResultSet etc.
  Example    : $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  DESCRIPTION: Validates object class and stored status
  Returntype : none
  Exceptions : Throws if Storable is not valid or stored
  Caller     : general - Adaptors, objects will probably be better off implementing in situ.
               This is to avoid having to test for the adaptor for every object which could slow things down
  Status     : At risk

=cut


#This has to be in the DBAdaptor rather than Storable as we're 
#calling isa on self otherwise which we don't know whether we can

sub is_stored_and_valid{
  my ($self, $class, $obj) = @_;

  if(! (ref($obj) && $obj->isa($class) && $obj->is_stored($self))){
	#is_stored checks adaptor params and dbID, but not whether the adaptor matches the class
	throw('Must provide a valid stored '.$class."\nParameter provided was:\t$obj");
  }

  return;
}

=head2 are_stored_and_valid

  Arg [1]    : string - class namespace
  Arg [2]    : ARRAYREF os Bio::EnsEMBL::Funcgen::Storable objects e.g. ResultSet 
  Arg [3]    : String : return value method name
  Example    : $db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', \@rsets);
  DESCRIPTION: Wrapper for is_stored_and_valid. Will optionally return array of values
               defined by calling method name arg on each object passed 
  Returntype : ARRAYREF - contents defined by optional method name arg
  Exceptions : Throws if object list is not an ARRAY with at least one element
  Caller     : general 
  Status     : At risk

=cut

#Add method params?

sub are_stored_and_valid{
  my ($self, $class, $obj_list, $method_name) = @_;

  my @return_vals;

  if( (ref($obj_list) ne 'ARRAY') ||
	  (scalar(@$obj_list) <=0) ){
	throw('You must provide an ARRAYREF of objects to validate');
  }

  foreach my $obj(@$obj_list){
	$self->is_stored_and_valid($class, $obj);

	if($method_name){
	  #test can method here?
	  push @return_vals, $obj->$method_name;
	}
  }

  return \@return_vals;
}



#Move these to Helper.pm! Check method dependencies first!

=head2 load_table_data

  Arg [1]    : string - table name
  Arg [1]    : string - file path for file to load
  Example    : $db->load_table_data("result",  $self->get_dir($results_dir)."/result.txt");
  DESCRIPTION: Generic method to load a file into a specified table
  Returntype : none
  Exceptions : Throws if argument not supplied
  Caller     : general
  Status     : At risk - only used by for results at present, to be removed

=cut

sub load_table_data{
  my ($self, $table, $file, $ssh) = @_;

  chmod 0755, $file;

  #  warn("Importing $table data from $file");
  #if this gives an Errcode: 2, then your mysql instance cannot see the file.
  #This could be due to a soft link on a visible directory to an unmounted filesystem
  #change this to use the mysqlimport?



  #This is failing as ssh is not set up to login silently without password prompt
  #Need to defined ssh keys?

  #(my $tmp_file = $file) =~ s/.*\///;
  #$tmp_file = '/tmp/'.$tmp_file;

  #my $scp = 'scp $(hostname):'.$file." ".$self->dbc->host().":${tmp_file}";
  #my $sql = "load data infile '$tmp_file' into table $table";
  #$self->dbc->do($sql);
  #remove tmp file via ssh if load successful

  my $cmd = 'mysqlimport -L '.$self->connect_string().' '.$file;
  system($cmd) == 0 || throw("Failed to load data from $file\nExit code:\t".($?>>8)."\n$!");
  
  return;
}



=head2 get_available_adaptors

  Example    : my %pairs = %{$dba->get_available_adaptors()};
  Description: gets a hash of the available adaptors
  ReturnType : reference to a hash
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::ConfigRegistry
  Status     : Stable

=cut


#will adding SliceAdaptor here use the dna DB? i.e. the core DB rather than the efg DB?

sub get_available_adaptors{
  my ($self) = shift;
  
  my %pairs = (
			   'Channel'            => 'Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor',
			   'ExperimentalChip'   => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor',
			   'ArrayChip'          => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor',
			   'Array'              => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor',
			   'ProbeSet'           => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor',
			   'Probe'              => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor',
			   'ProbeFeature'       => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor',
			   'AnnotatedFeature'   => 'Bio::EnsEMBL::Funcgen::DBSQL::AnnotatedFeatureAdaptor',
			   'RegulatoryFeature'  => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor',
			   'Experiment'         => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
			   'ExperimentalGroup'  => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor',
			   'DataSet'            => 'Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor',
			   'FeatureType'        => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor',
			   'FGCoordSystem'      => 'Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor',#prepended FG to override core adaptor?
			   'MetaCoordContainer' => 'Bio::EnsEMBL::Funcgen::DBSQL::MetaCoordContainer',
			   'FeatureSet'         => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor',
			   'ResultSet'          => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor',
			   'DataSet'            => 'Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor',
			   'InputSet'           => 'Bio::EnsEMBL::Funcgen::DBSQL::InputSetAdaptor',
			   'ExternalFeature'    => 'Bio::EnsEMBL::Funcgen::DBSQL::ExternalFeatureAdaptor',
			   'CellType'           => 'Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor',
			   'DBEntry'            => 'Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor',
			   'Slice'              => 'Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor',
			   'ResultFeature'      => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor',
			   'MotifFeature'       => 'Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor',
			   'BindingMatrix'      => 'Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor',
			   'SegmentationFeature'=> 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationFeatureAdaptor',
			   
			 
			   #add required EnsEMBL(core) adaptors here
			   #Should write/retrieve from efg not dna db
			   'UnmappedObject'     => 'Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor',
			   'Analysis'           => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
			   "MetaContainer"      => 'Bio::EnsEMBL::DBSQL::MetaContainer',
	      );
  
  return (\%pairs);
}

=head2 _get_schema_build

  Arg [1]    : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor or Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $shema_build = $db->_get_schema_build($slice->adaptor->db());
  DESCRIPTION: 
  Returntype : string
  Exceptions : Throws if argument not supplied
  Caller     : general
  Status     : At risk - replace with MetaContainer method

=cut


#Slightly hacky convinience method to get the data/schema.version/build from a feature slice

sub _get_schema_build{
  my ($self, $db) = @_;

  #Have to explicitly pass self->db to this method if required, this highlights which db is being tested 
  throw("Need to define a DBAdaptor to retrieve the schema_build from") if (! $db);
  #avoided using dnadb by default to avoid obfuscation of behaviour
  
  my @dbname = split/_/, $db->dbc->dbname();



  my $schema_build = pop @dbname;
  $schema_build = pop(@dbname).'_'.$schema_build;
  return $schema_build;
}

#=head2 get_SliceAdaptor
#
#  Arg [1]    : (optional) int - coord_system_id
#  Example    : my $slice_adaptor = $db->get_SliceAdaptor($cs->dbID());
#  DESCRIPTION: Retrieves a slice adaptor from the dnadb corresponding 
#               to the coord_system_id, or retrieves from the default dnadb
#  Returntype : Bio::EnsEMBL::DBSQL::SLiceAdaptor
#  Exceptions : Throws if arguments not supplied
#  Caller     : general
#  Status     : At risk - remove and add this to BaseFeatureAdaptor->fetch_all_by_Slice_constraint
#
#=cut

#Funcgen specific, get's Adaptor from dnadb, or validates/autogenerates from coord_system_id
#Only imlpmented in _obj_from_sth, rely on feature_slice elsewhere

#Not in registry as get_adaptor will not take $cs_id arg

#Move all this dnadb specif stuff to dnadb, to ensure all dnadb derived object are from correct DB
#All dnadb centric methods should then either use the default or pass a new coordsysid to redefine the dnadb
#Should we make this mandatory to ensure dnadb is redefined, this would avoid getting data from wrong db, but maybe a pain in the butt
#also, changing dnadb would work, which isn't pretty

#Are all dnadb(feature) data retrievals mediated by a Slice?
#ProbeFeatureADaptor has by probe/probeset queries which would retrieve for all DBs/coord systems,
#any further dnadb derived methods on the objects would have to resolve coord system issue and use correct dnadb
#or should we only retrieve for current dnadb?

#rename this DNADB|FGSliceAdaptor?
#as this works differently to normal method
#the problem arises when we get features from the DB by none Slice methods, these may not refer to the current dnadb
#so we have to implement checks in non slice based feature calls to make sure we nest the correct dnadb adaptor

#sub get_SliceAdaptor{
#  my ($self, $cs_id) = @_;

  #$cs_id is only used in ProbeFeatureAdaptor, no longer used
  #but is this correct?


  #Need to add check if current cs_id refers to current dnadb
  
  #extract this to a "validate_dnadb" method
  #This will be called for each noon Slice based fetch method for each feature returned
  #or should we group the fetch statements by coord system id and try and do it more efficiently
  
  #is this "validate_coordsystem"?
  
  #Can we cache the DNA DBAdaptors against the FG csis rather than doing this everytime?
  #will this be too much memory overhead? Registry is already a cache, can we just reference the registry?
 


 
#  if($cs_id){
#    my $csa = $self->get_FGCoordSystemAdaptor();
#    my $fg_cs = $csa->fetch_by_dbID($cs_id);
    #my $schema_build = $fg_cs->schema_build();
    #Get species here too
    
#	if(! $fg_cs->contains_schema_build($self->_get_schema_build($self->dnadb()))){
    #if($schema_build ne $self->_get_schema_build($self->dnadb())){
#	  my $lspecies = $reg->get_alias($self->species());
      #warn "Generating dnadb schema_build is $schema_build and dnadb is ".$self->_get_schema_build($self->dnadb())."\n";

      #get from cs_id
      #can we return direct from registry for older versions?
      #best to generate directl as we may have only loaded the current DBs
      #set dnadb here and return after block

	  #should we really permanently set this here
	  #what is we were on ens-livemirror?
	  #we would then lose that association
	  #should we change dnadb to be totally dynamic anyway and only set it for the current default?

#	  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
#		(						
#		 -host => "ensembldb.ensembl.org",
#		 -user => "anonymous",
#		 -species => $lspecies,
#		 -dbname => $lspecies.'_core_'.$fg_cs->get_latest_schema_build(),
#		 -group => 'core',
#		 #-port  => 5306,
#		);


	  #This new port only has from 48 onwards!!!
  
      
#      $self->dnadb($dnadb);
      
#    }
#  }
  
#  return $self->dnadb->get_SliceAdaptor();#this causes circular reference if dnadb not set i.e if this is generated from scratch without a dnadb rather than from the reg?????
#}


sub dnadb_host{
  my ($self, $host) = @_;
  $self->{'dnadb_host'} = $host if $host;
  return $self->{'dnadb_host'};
}

sub dnadb_port{
  my ($self, $port) = @_;
  $self->{'dnadb_port'} = $port if $port;
  return $self->{'dnadb_port'};
}
 
sub dnadb_pass{
 my ($self, $pass) = @_;
 $self->{'dnadb_pass'} = $pass if $pass;
 return $self->{'dnadb_pass'};
}

sub dnadb_user{
  my ($self, $user) = @_;
  $self->{'dnadb_user'} = $user if $user;
  return $self->{'dnadb_user'};
}

sub dnadb_assembly{
  my ($self, $assm) = @_;
  $self->{'dnadb_assm'} = $assm if $assm;
  return $self->{'dnadb_assm'};
}

#Redefine dbadb here to add coordsystem

=head2 dnadb

  Arg [1]:     Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]:     string - coord_system name e.g. chromosome
  Usage :      my $dnadb = $db->dnadb(); 
  Description: returns the database adaptor where the dna lives i.e. the core db for a given species
               There are at least 2 cases where you need to set this explicitly
               1.  If you want to retrieve features on an assembly which is not the default in 
               the correspeonding core DB with matching schema_build
               2.  If the corresponding core DB is not available on the default ensembl DB 
               server(ensembldb/ens-livemirror) i.e. before a new release.
  Status :     At risk. - Might remove validation of CS

=cut

#This is not taking account of the registry which may have already been loaded
#So we may be setting the dnadb correctly here
#But it won't be the default core db in the registry, it will be cached as species1 or something?

sub dnadb { 
  my ($self, $dnadb, $cs_name) = @_; 

  #super dnadb automatically sets the current DBAdaptor as the dnadb
  #this is the only way of checking whether it has been defined properly.
 

  #This needs to use the same host as the registry, which may differ from the efg host!
  #What if we pass a reg config file? with multiple hosts?
  #We probably always want to use the host of the current dnadb?

  #Automatically selects the current core DB from the reg and may be wrong!!
  #How do we want to deal with this?
  
  #So we want a method in DBAdaptor to reset efg DB and dna DB in registry?
  #Or do we just let this happen

  if($dnadb || $self->SUPER::dnadb->group() ne 'core'){

	if(! $dnadb){#Guess dnadb
	  #We are only guessing the dnadb if we have not already defined one
	  #Either by attaching to the efg db
	  #Or by loading the reg
	  #The Importer will automatically check if the dnadb matches and set_by_assembly_version
	  #Do we want this functionality in the efg db->new?
	  #The problem is redefining the db in the registry will not take the extra params required?
	  #Maybe the method will take this. Then we can implement this check and set_dnadb_by_assembly_version
	  #in new.
	  
	  

	  return $self->set_dnadb_by_assembly_version($self->dnadb_assembly);
	}
	
	
	$self->SUPER::dnadb($dnadb); 

	#set default coordsystem here
	#there might not be a chromosome level if we just have a scaffold assembly
	#supercontig is already loaded as we use toplevel?
	#This should really get all the default CSs from the core CSAdaptor?
	#We should also do this during the update?
	#This will also enable people to query using clone/contig/supercontig level slices
	#How will this work? Where will the mapping between CSs be done?
	
	my @cs_names;
	@cs_names = ($cs_name) if $cs_name;
	#$cs_name ||= 'chromosome';

	if(! $cs_name){
  
	  foreach my $cs(@{$dnadb->get_CoordSystemAdaptor->fetch_all_by_attrib('default_version')}){
		push @cs_names, $cs->name;
	  }
	}


	foreach my $cs_name(@cs_names){

	  my $cs;
	  eval { $cs = $dnadb->get_CoordSystemAdaptor->fetch_by_name($cs_name)};
	  my $error = $@;
	  
	  if($error){
		my ($schema, $build) = split/_/, $self->_get_schema_build($dnadb);
		$build =~ s/[a-z]//;
		throw("It appears that the schema of ".$dnadb->dbc->dbname.
			  ' is incompatible with your current core API version('.$reg->software_version.
			  ").  You could try using the $schema version of the core API, or alternatively try specifying ".
			  "different -dnadb/registry_host parameters to point to a make recent version containing build $build\n");
		
	  }
	  
	  
	  #this will only add the default assembly for this DB, if we're generating on another we need to add it separately.
	  #or shall we fetch/add all by name?
	  
	  #This is a non-obious store behaviour!!!!!!!!!!!!!!!!!
	  #This can result in coord_system entries being written
	  #unknowingly if you are using the efg DB with a write user/pass
	  $self->get_FGCoordSystemAdaptor->validate_and_store_coord_system($cs);
	}
  }

  return $self->SUPER::dnadb();#never pass @_ here!
} 


=head2 set_dnadb_by_assembly_version

  Arg [1]:     string - Assembly version e.g. for homo_sapiens_core_49_36k it would be 36
  Usage :      $efgdb->set_dnadb_by_assembly_version('36'); 
  Description: Sets the dnadb to the latest version given the assembly version
  Exceptions:  Throws if no assembly version provided or cannot for appropriate dnadb on ensembldb
  Status :     At risk

=cut



sub set_dnadb_by_assembly_version{
  my ($self, $assm_ver) = @_;

  throw('Must provide and assembly version to set the dnadb') if ! defined $assm_ver;

  my $reg_lspecies = $reg->get_alias($self->species());
   #The registry has incremented the species as we have recreated the efg DB
  #possibly using a different schema_build
  #This set true lspecies to allow dnadb detection
  #in multi DB environments e.g. DAS server
  my $lspecies = $reg_lspecies;
  $lspecies =~ s/[0-9]+$// if($lspecies =~ /[0-9]$/);


  throw('Either specify a species parameter or set species.production_name in the meta table to set dnadb automatically, alternatively pass a dnadb parameter') if $lspecies eq 'default';
	
  #So we use params first
  #else registry params
  #else ensembldb
 
  my @ports = ($self->dnadb_port);

  #Start with lastest MySQL instances
  #We are over-riding specified port here, only for known hosts
  #we should really account for this and make it nr
  if($self->dnadb_host eq 'ensdb-archive'){#
	@ports = (5304, 3304);
  } 
  elsif($self->dnadb_host eq 'ensembldb.ensembl.org'){
	@ports = (5306, 3306);
  }

  
  #We should probably allow for non-ensembldb core DBs here too
  #Do we need to account for other ports, staging etc for release?
  #These should run fine on 3306.
  #my $current_port = $self->dnadb->dbc->port;
  #This is assuming dnadb port will only ever be one or the other
  #This assumption is restricted to ensembldb in the port loop
  #my $tmp_port = ($current_port == 3306) ? 5306 : 3306;
  
  #Do we need to get species from reg here?
  #And test the species has been set?

  my $sql = 'show databases like "'.$lspecies.'_core_%_'.$assm_ver.'%"';
  my ($dbh, @dbnames, $port, $host_port);

  foreach $port(@ports){	
	#This is probably duplicating connections and over-riding any connection 
	#pooling going on in the base DBConnection if we are using the same host port
	#as the registry connection

	$dbh = DBI->connect('DBI:mysql:host='.$self->dnadb_host.";port=${port}",
						$self->dnadb_user, 
						$self->dnadb_pass, 
						{'RaiseError' => 1});

	eval { @dbnames = map {$_ = "@$_"} @{$dbh->selectall_arrayref($sql)}; };

	if($@){
	  throw('Failed to fetch dna DB names from '.$self->dnadb_host.":${port}"."\n$@");
	}



	#sort and filter out non-core DBs
	#This will always take the latest release, not the latest genebuild version
	#Which is probably what we want anyway
	@dbnames = grep(/core_[0-9]/, sort @dbnames);
  
	if(scalar(@dbnames)==0){
	  warn(':: Failed to find '.$self->species.' core DB for assembly version '.$assm_ver.' using '
		   .$self->dnadb_user.'@'.$self->dnadb_host.':'.$port);
	}
	else{
	  $host_port = $port;
	  last;
	}
  }

  throw("Failed to find dnadb with assembly version $assm_ver. Maybe you want to set -dnadb_host?") if(scalar(@dbnames)==0);


  warn ":: Auto-selecting build $assm_ver core DB as:\t".
	$self->dnadb_user.'@'.$dbnames[$#dbnames].':'.$self->dnadb_host.':'.$host_port."\n";


  my $db = $reg->reset_DBAdaptor($reg_lspecies, 'core', $dbnames[$#dbnames], $self->dnadb_host, $host_port, $self->dnadb_user, $self->dnadb_pass);
  
  $self->dnadb($db);
  return $db;
}





#General Status methods
#will Move to Bio::EnsEMBL::Funcgen::DBSQL::Status

=head2 fetch_all_states

  Arg [1]    : string - table name
  Arg [2]    : int - table id
  Example    : my @states = @{$db->fetch_all_states('channel', 1)};
  Description: Retrieves all states associated with the given table record
  Returntype : Listref
  Exceptions : Throws if arguments not supplied
  Caller     : general
  Status     : At risk - Move to Status

=cut

sub fetch_all_states{
	my ($self, $table, $id) = @_;


	throw("DBAdaptor::fetch_all_states is deprecated");


	throw("Need to specifiy a table and an id to retrieve status") if (! $table || ! $id);


	my $sql = "SELECT state FROM status WHERE table_name=\"$table\" AND table_id=\"$id\"";

	my @states = map{ $_ = "@$_"} @{$self->dbc->db_handle->selectall_arrayref($sql)};

	return \@states;
}


=head2 fetch_status_by_name

  Arg [1]    : string - table name
  Arg [2]    : int - table id
  Arg [3]    : string - status
  Example    : if($db->fetch_status_by_name('channel', 1, 'IMPORTED'){ ... };
  Description: Retrieves given state associated with the table record
  Returntype : ARRAYREF
  Exceptions : Throws if arguments not supplied
  Caller     : general
  Status     : At risk - Move to Stasus

=cut



sub fetch_status_by_name{
	my ($self, $table, $id, $state) = @_;

	throw("DBAdaptor::fetch_status_by_name is deprecated");

	throw("Need to specify a table and an id to retrieve status") if (! $table || ! $id || ! $state);

	#should we enum the state?


	my $sql = "SELECT state FROM status WHERE table_name=\"$table\" AND table_id=\"$id\" AND state=\"$state\"";
	return $self->dbc->db_handle->selectrow_array($sql);
}


=head2 set_status

  Arg [1]    : string - table name
  Arg [2]    : int - table id
  Arg [3]    : string - status
  Example    : $db->set_status('channel', 1, 'IMPORTED');
  DESCRIPTION: RETRIEVES GIVEN STATE ASSOCIATED WITH THE table record
  Returntype : ARRAYREF
  Exceptions : Throws if arguments not supplied
  Caller     : general
  Status     : At risk - Move to Status

=cut


sub set_status{
	my ($self, $table, $id, $state) = @_;

	throw("DBAdaptor::set_status is deprecated");

	throw("Need to supply a table, dbid and a valid status") if (!($table && $id && $state));

	my $sql = "INSERT INTO status(table_id, table_name, state) VALUES(\"$id\", \"$table\", \"$state\")";
	$self->dbc->do($sql);

	return;
}


sub stable_id_prefix{
  my $self = shift;

  if(! defined $self->{'stable_id_prefix'}){
	($self->{'stable_id_prefix'}) = @{$self->dnadb->get_MetaContainer->list_value_by_key('species.stable_id_prefix')};

	#Only add R if it is defined
	$self->{'stable_id_prefix'} .= 'R' if $self->{'stable_id_prefix'};
  }

  return $self->{'stable_id_prefix'}
}


=head2 connect_string

  Example    : my $import_cmd = 'mysqlimport '.$db->connect_string()." $table_file";
  Description: Retrieves the mysql cmdline connection string
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub connect_string{
  my $self = shift;

  return '-h'.$self->dbc->host().' -u'.$self->dbc->username().' -p'.$self->dbc->password()
	.' -P'.$self->dbc->port().' '.$self->dbc->dbname();
}



# DEPRECATED METHODS #

sub fetch_group_details{
  my ($self, $gname) = @_;

  deprecate("Please use ExperimentalGroupAdaptor");

  throw("Need to specify a group name") if ! $gname;
  my $sql = "SELECT * from experimental_group where name=\"$gname\"";
  return $self->dbc->db_handle->selectrow_array($sql);
}

sub import_group{
  my ($self, $gname, $loc, $contact) = @_;

  deprecate('Please use ExperimentalGroup/Adaptor to import experimental_group info');
  throw('import_group no longer supported');
}


1;

