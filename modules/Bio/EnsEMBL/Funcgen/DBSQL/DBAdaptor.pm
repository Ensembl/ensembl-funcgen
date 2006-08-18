
=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  
=head1 SYNOPSIS




=back

=head1 DESCRIPTION

B<This program> 

=cut

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);# Bio::EnsEMBL::Funcgen::Helper);

use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
#use Bio::EnsEMBL::Funcgen::Helper;
my $reg = "Bio::EnsEMBL::Registry";

################################################################################

=head2 new

 Description : 
               

 Arg  [1]    : hash containing optional attributes :-
               
 ReturnType  : Experiment

 Example     : my $Exp = Bio::EnsEMBL::Funcgen::DBAdaptor->new(
                                                      
                                                     );

 Exceptions  : 

=cut

################################################################################



#these should be removed if overlap with core DBAdaptor or moved to ImportAdaptor?
# or use store methods on each object

sub fetch_dbid_by_table_field{
	my ($self, $table, $name, $field) = @_;

	$field ||= "name";

	my $sql = "select ${table}_id from $table where $field =\"$name\"";

	#print "sql is $sql\n";

	return $self->dbc->db_handle->selectrow_array($sql);
}	

sub fetch_channel_dbid_by_echip_dye{
	my ($self, $chip_dbid, $dye) = @_;

	my $sql = "select channel_id from channel where experimental_chip_id =\"$chip_dbid\" and dye = \"$dye\"";

	return $self->dbc->db_handle->selectrow_array($sql);
}	


sub insert_table_row{
	my ($self, $table, @values) = @_;

	#This assumes first field is an auto increment dbid field
	my $sql = "insert into $table values (\"\", \"". join("\", \"", @values)."\")";

	#print "sql is $sql\n";
	#we need to catch this and stack trace?
	return $self->dbc->do($sql);#do we really want to return? Or capture fail here and warn
}


sub get_last_table_id{
	my ($self, $table) = @_;
	
	#can't use multiline query :(
	#my $sql = "insert into $table (${table}_id) values (\"\"); select LAST_INSERT_ID()";
	
	#my $sql = "insert into $table (${table}_id) values (\"\")";
	#$self->dbc->do($sql);
	#Don't need to check if if any higher as insert using autoincrement will select highest unused id
	#return $self->dbc->db_handle->last_insert_id(undef, undef, undef, undef);	

	my $sql = "select ${table}_id from $table order by ${table}_id desc limit 1";
	return $self->dbc->db_handle->selectrow_array($sql);
}

sub load_table_data{
	my ($self, $table, $file) = @_;

	#$self->log("Loading $table data from $file");
	my $sql = "load data infile \"$file\" into table $table";
	$self->dbc->do($sql);
	#$self->log("Finished loading $table data");
	return;
}

#Only validates if already present
#add flag to alter table if any inconsistencies found?
#This could be heavily utilised in the recovery method to avoid having to delete entries for incomplete imports
sub register_entry{
	my ($self, $table, $dbid, $alter, @data) = @_;

	my ($sql, @db_entry);
	my $valid = 1;#just throw here instead? or warn as may be able
	$self->throw("Alter entry function not yet implemented") if ($alter);

	if(! defined $dbid){
		$self->insert_table_row($table, @data);
	}else{
		#$self->log("$table $dbid already registered, validating");
		$sql = "select * from $table where ${table}_id = \"$dbid\"";

		@db_entry = $self->dbc->db_handle->selectrow_array($sql);
		shift @db_entry;

		for my $i(0..$#db_entry){

			if($db_entry[$i] ne $data[$i]){
				$valid = 0;
				#should really get field names here too, and store if doing repeat validations to reduce no. of queries
				warn("Validating $table entry $dbid found the following mismatch:\nData\t$data[$i]\nDB\t$db_entry[$i]");
			}
		}
	}
	
	return $valid;#return value to allow caller to throw
}


#Used by ConfigRegistry to make adaptors available
#will adding SliceAdaptor here use the dna DB? i.e. the core DB rather than the efg DB?

sub get_available_adaptors{
  my ($self) = shift;
  
  my %pairs = (
	       'Channel'            => 'Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor',
	       'Experimental'       => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
	       'ExperimentalChip'   => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor',
	       'OligoArray'         => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor',
	       'OligoProbeSet'      => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoProbeSetAdaptor',
	       'OligoProbe'         => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoProbeAdaptor',
	       'OligoFeature'       => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoFeatureAdaptor',
	       'PredictedFeature'   => 'Bio::EnsEMBL::Funcgen::DBSQL::PredictedFeatureAdaptor',
	       'Experiment'         => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
	       'ResultSet'          => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor',
	       'FGCoordSystem'      => 'Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor',#prepended FG o override core  adaptor
	       'MetaCoordContainer' => 'Bio::EnsEMBL::Funcgen::DBSQL::MetaCoordContainer',
	       
	       
	       #add required EnsEMBL(core) adaptors here
	       #Should write/retrieve from efg not dna db
	       'Analysis'           => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
	       "MetaContainer"      => "Bio::EnsEMBL::DBSQL::MetaContainer",
	      );
  
  return (\%pairs);
}


#Hacky convinience method to get the data/schema.version/build from a feature slice

sub _get_schema_build{
	my ($self, $adaptor) = @_;


	throw("Need to define a DBAdaptor to retrieve the schema_build from") if (! $adaptor);
	
	my $schema_build = $adaptor->dbc->dbname();
	$schema_build =~ s/[a-zA-Z_]*//;

	return $schema_build;
}

#Funcgen specific, get's Adaptor from dnadb, or validates/autogenerates from coord_system_id
#Only imlpmented in _obj_from_sth, rely on feature_slice elsewhere

#Not in registry as get_adaptor will not take $cs_id arg

sub get_SliceAdaptor{
	my ($self, $cs_id) = @_;

	if(! $cs_id && ! $self->dnadb()){
		throw("Need to set a valid dnadb or pass a Funcgen coord_system_id");
	}

	if($cs_id){
		my $csa = $self->get_FGCoordSystemAdaptor();
		my $fg_cs = $csa->fetch_by_dbID($cs_id);
		my $schema_build = $fg_cs->schema_build();
		
		if($self->dnadb()){#validity check here, but let next block return
			if($schema_build ne $self->_get_schema_build($self->dnadb())){
				throw("Supplied coord_system_id does not map to schema_build of dnadb");
			}
		}
		else{#get from cs_id
			#can we return direct from registry for older versions?
			#best to generate directl as we may have only loaded the current DBs
			#set dnadb here and return after block

			#Need to get latin name for DB species from somewhere, meta?
			my $species = "homo_sapiens";


			warn("This needs to handle querying local or ensembl DB");

			my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
			  (						
			   -host => "ensembldb.ensembl.org",
			   -user => "anonymous",
			   -dbname => "${species}_core_${schema_build}",
			  );

			$self->dnadb($dnadb);

		}
	}

	return $self->dnadb->get_SliceAdaptor();
}


#Redefine dbadb here to add coordsystem

=head2 dnadb

  Title :   dnadb 
  Usage :   my $dnadb = $db->dnadb(); 
  Function: returns the database adaptor where the dna lives Useful if you only want 
            to keep one copy of the dna on disk but have other databases with genes and 
            features in Returns : dna database adaptor 
  Args :    Bio::EnsEMBL::DBSQL::BaseAdaptor
  Status : Medium Risk. 
         : Use the Registry method add_DNAAdaptor/get_DNAAdaptor instead??  Not for eFG

=cut


sub dnadb { 
	my $self = shift; 

	if(@_) { 

		$reg->add_DNAAdaptor($self->species(),$self->group(),$_[0]->species(),$_[0]->group()); 

		#set default coordsystem here
		my $cs = $_[0]->get_CoordSystemAdaptor->fetch_by_name('chromosome');
		#this will only add the default assembly for this DB, if we're generating on another we need to add it separately.
		#or shall we fetch/add all by name?

		$self->get_FGCoordSystemAdaptor->validate_coord_system($cs);
	}

	return $self->SUPER::dnadb(@_);
} 


#Group methods, as not adaptor/class for Group(Incorporated in Experiment)

sub fetch_group_details{
	my ($self, $gname) = @_;

	throw("Need to specify a group name") if ! $gname;
	my $sql = "SELECT * from egroup where name=\"$gname\"";
	return $self->dbc->db_handle->selectrow_array($sql);
}



sub import_group{
	my ($self, $gname, $loc, $contact) = @_;

	throw("Need to supply a group name, location and contact") if (!($gname && $loc && $contact));


	my $sql = "INSERT INTO egroup(name, location, contact) VALUES(\"$gname\", \"$loc\", \"$contact\")";
	$self->dbc->do($sql);


	#$self->dbc->db_handle->last_insert_id(undef, undef, undef, undef);	
	return;#return last insert id here?
}


sub fetch_all_states{
	my ($self, $table, $id) = @_;

	throw("Need to specifiy a table and an id to retrieve status") if (! $table || ! $id);


	my $sql = "SELECT state FROM status WHERE table=\"$table\" AND table_id=\"$id\"";
	return $self->dbc->db_handle->selectall_arrayref($sql);
}


sub fetch_status_by_name{
	my ($self, $table, $id, $state) = @_;

	throw("Need to specifiy a table and an id to retrieve status") if (! $table || ! $id || ! $state);


	my $sql = "SELECT state FROM status WHERE table_name=\"$table\" AND table_id=\"$id\" AND state=\"$state\"";
	return $self->dbc->db_handle->selectrow_array($sql);
}


sub set_status{
	my ($self, $table, $id, $state) = @_;

	throw("Need to supply a table, dbid and a valid status") if (!($table && $id && $state));

	my $sql = "INSERT INTO status(table_id, table_name, state) VALUES(\"$id\", \"$table\", \"$state\")";
	$self->dbc->do($sql);

	return;
}

	       
1;

