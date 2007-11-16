=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::Helper
  
=head1 SYNOPSIS


 e.g. 


 my $object = Bio::EnsEMBL::Object->new
 (
     logging     => 1,
     log_file    => "/tmp/Misc.log",
     debug_level => 2,
     debug_file  => "/tmp/Misc.dbg",
 );

 $object->log("This is a log message.");
 $object->debug(1,"This is a debug message.");
 $object->system("rmdir /tmp/test");


 ----------------------------------------------------------------------------


=head1 OPTIONS

=over 8


=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> performs several health check and update methods prior to release.

=cut

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::HealthChecker;

use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Analysis;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

################################################################################

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = {};

  bless $self, $class;

  #validate and set type, analysis and feature_set here
  my ($db, $builds) = rearrange(['DB', 'BUILDS'], @_);
  
  throw('You must define a type of external_feature to import') if(! defined $type);

  if (! ($db && ref($db) &&
		 $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  #test connection
  $db->dbc->db_handle;
  
  $self->{'db'} = $db;
  $self->{'builds'} = (@$builds) ? @$builds : ('DEFAULT');
  
  return $self;
}

sub db{
  my ($self) = @_;

  return $self->{'db'};
}


#wrapper method

=hea2

  Arg[0]     : boolean - 
  Arg[1]     : boolean -
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : StableArg[0] - boolean force

=cut

sub update_db_for_release{
  my ($self, @args) = @_;

  if(@args)

  #do seq_region_update to validate dnadb first
  #hence avoiding redoing longer methods
  $self->validate_new_seq_regions($force_srs);
  $self->
  

  $self->log(' :: Finished updating '.$self->db->dbhost.':'.$self->db->dbc->dbname.' for release');
}

sub validate_new_seq_regions{
  my ($self, $force) = @_;

  

  #do we need to add the none default levels here?
  #or are we only bothered about those which constitute the toplevel?

  #To make sure we have all the correct levels in eFG we need to get all the names.
  #then get all by name from the core db and set them as the dnadb.
  # we also need to get all the toplevel seq_regions and store them in the seq_region table
  #use BaseFeatureAdaptor::_pre_store with and array of pseudo feature on each top level slice

  #Validate the efgdb and dnadb schema version are the same first
  
  if(! $force){
	my $efgdb_sm = join('_', @{$self->get_schema_and_build($self->db->dbc->dbname)});
	my $dnadb_sm = join('_', @{$self->get_schema_and_build($self->db->dnadb->dbc->dbname)});
						
	if($efgdb_sm ne $dnadb_sm){
	  $self->log("WARNING Skipped validate_new_seq_regions as schema_versions are mismatched:\t".
				 "efgdb $efgdb_sm\tdnadb $dnadb_sm");
	  return 0;
	}
  }
  
  my $pf_adaptor = $efg_db->get_ProbeFeatureAdaptor();
  my $slice_adaptor = $efg_db->get_SliceAdaptor();
  
  $self->log(" :: Validating new coord_systems/seq_regions :: ::");

  foreach my $build(@builds){
	
	$self->log("Importing seq_region/coord_system info for build:\t".$build);
	
	foreach my $slice(@{$slice_adaptor->fetch_all('toplevel', $build)}){
	  
	  if($slice->start() != 1){
		#we must have some sort of PAR linked region i.e. Y
		$slice = $slice_adaptor->fetch_by_region($slice->coord_system_name(), $slice->seq_region_name());
	  }
	  
	
	  #we need test if it needs doing first?
	  #we would need to test for the coord_systems outside of this loop
	  #and then for each seq_region inside the loop if the coord_system is present

	  $self->log("_pre_storing seq_region info for slice:\t".$slice->name());
	  
	  my $pseudo_feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
		(
		 -slice => $slice,
		 -start => $slice->start(),
		 -end   => $slice->end(),
		 -strand => 0,
		);
	  
	  $pf_adaptor->_pre_store($pseudo_feature);
	  
	}
  }

  $self->log("Finished validating seq_regions\n");

  return;
}



sub update_meta_schema_version{
  my ($self) = @_;

  my $schema_version = $self->get_schema_and_build($self->db->dbname)->[0];

  my $sql = "UPDATE meta set meta_value='.$schema_version.' where meta_key='schema_version'";
  $efg_db->dbc->db_handle->do($sql);

  $self->log(" :: Updated meta.schema_version to $schema_version :: ::\n");

}


sub update_meta_coord{
  my ($self, @table_names) = @_;

  my $sql = 'UPDATE meta set meta_value=48 where meta_key="schema_version"';
  $efg_db->dbc->db_handle->do($sql);

  #set default table_name
  if(! @table_names || scalar(@table_names) == 0){

	@table_names = qw(
					  regulatory_feature
					  probe_feature
					  external_feature
					  annotated_feature
					 );
  }

  #backup meta coord
  if(system("mysql -h$host -P$port -u$user -p$pass -N "
			. "-e 'SELECT * FROM meta_coord' ${species}_funcgen_${schema_build}"
			. "> ${species}_funcgen_${schema_build}.meta_coord.backup"
		   ) != 0 ){
	throw("Can't dump the original meta_coord for back up");#will this get copied to log?
  } 
  else {
	$self->log(" :: Updating meta_coord table. Original meta_coord table backed up in "
			   . "${species}_funcgen_${schema_build}.meta_coord.backup :: ::\n");
  }


  #Update each max_length for table_name and coord_system
  foreach my $table_name(@table_names){
	my $sql1 = "select distinct(cs.name), mc.coord_system_id, cs.version, mc.max_length from coord_system cs, meta_coord mc where mc.table_name='$table_name' and mc.coord_system_id=cs.coord_system_id";
	
	$self->log("Updating meta_coord max_length for $table_name:\nname\tcoord_system_id\tversion\tmax_length");
	
	#can we test for emtpy array here? Then skip delete.
	
	my @info = @{$efg_db->dbc->db_handle->selectall_arrayref($sql1)};
	
	map {print join("\t", @{$_})."\n"} @info;
	
	# Clean old entries
	$self->log("Deleting old meta_coord entries");
	my $sql = "DELETE FROM meta_coord WHERE table_name ='$table_name'";
	$efg_db->dbc->db_handle->do($sql);
	
	# Generate new max_lengths
	$self->log("Generating new max_lengths");
	
	#Is this query running for each redundant cs_id?
	#would it be more efficient to retrieve the NR cs_ids first and loop the query for each cs_id?

	$sql =
	  "INSERT INTO meta_coord "
		. "SELECT '$table_name', s.coord_system_id, "
		  . "MAX( t.seq_region_end - t.seq_region_start + 1 ) "
			. "FROM $table_name t, seq_region s "
			  . "WHERE t.seq_region_id = s.seq_region_id "
				. "GROUP BY s.coord_system_id";
  
	$efg_db->dbc->db_handle->do($sql);

	$self->log("New max_lengths for $table_name are:");
  
	@info = @{$efg_db->dbc->db_handle->selectall_arrayref($sql1)};

	map {$self->log(join("\t", @{$_})."\n")} @info;
  }  

  $self->log("Finished updating meta_coord max_lengths\n");

  return;
}

#Check displayable data_sets

sub check_displayable_sets{
  my $self = shift;
  
  $self->log(" :: Checking DataSets :: ::\n");

  my @dsets = @{$efg_db->get_DataSetAdaptor->fetch_all()};
  
  
  foreach my $dset(@dsets){
	print ":: Found DataSet ".$dset->display_name.'('.$dset->name.")\n";
	
	if($dset->is_displayable){
	  
	  #my $text = "
	  
	}
	
  }
  return;
}




### Check for regulatory meta entries for all regulatory feature_sets


1;
