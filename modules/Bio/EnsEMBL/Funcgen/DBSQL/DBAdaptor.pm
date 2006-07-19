
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
#use Bio::EnsEMBL::Funcgen::Helper;


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

#need to set group here somewhere
#only required when added to Registry/ConfigRegistry?


#sub new{
#    my ($caller, %args) = @_;

#    my ($self, %attrdata, $attrname, $argname);
#    my $class = ref($caller) || $caller;


    #define defaults here before call to DBAdaptor
    #change this to re-arrange?
    # objects private data and default values
    #none of these are set here, all passed to DBAdaptor
#    %attrdata = (
#		 _-dbname => "efg_test",
#		 #_-host   => "localhost",
#		 #_-port   => undef,#default handled in dbc
#		 _-user   => "ensadmin",
#		 _-pass   => undef,
#		 #_-dbc    => undef,#Pass existing connection to avoid multiple connections?
#		 _-group  => 'funcgen',#?
#		 _-species => 'Multi',#??????????????????
#		 #_-dnadb =>#?????????		 
#		);


#    # set each class attribute using passed value or default value
#    foreach $attrname (keys %attrdata){
#        ($argname = $attrname) =~ s/^_//; # remove leading underscore
#	$args{$argname} = (exists $args{$argname}  && defined $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
#    }


 #   #Create object from parent class
 #   $self = $class->SUPER::new(%args);


	
    #Do we need these? Yes!, but make all above undef(port?) and use GroupDefs to define
    
    #foreach my $tmp("dbname", "user", "host"){
    #	$self->throw("Mandatory arg $tmp not been defined") if (! defined $self->{"_${tmp}"});
    #}
    
    #if(! $self->dbc()){
    #	$self->dbc(new Bio::EnsEMBL::DBSQL::DBConnection(
    #-user   => $self->user(),
    #													 -host   => $self->host(),
    #													 -dbname => $self->dbname(),
    #													 -host   => $self->host(),
    #													 -pass   => $self->pass(),
    #													));
    #}
    
    
    
    
#    $self->debug(2,"DBAdaptor class instance created.");
#    $self->debug_hash(3, \$self);
    
#    return ($self);
#  }


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
				 'OligoArray'    => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor',
				 'OligoProbeSet' => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoProbeSetAdaptor',
				 'OligoProbe'    => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoProbeAdaptor',
				 'OligoFeature'  => 'Bio::EnsEMBL::Funcgen::DBSQL::OligoFeatureAdaptor',
				 'Experiment'    => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
				 'ResultSet'     => 'Bio::EnsEMBL::Funcgen::DBSQL::ResultSetAdaptor',
				 'FGCoordSystem'   => 'Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor',#prepended FG o override core  adaptor
				 #why doesn't this work as above?
				 'MetaCoordContainer' => 'Bio::EnsEMBL::Funcgen::DBSQL::MetaCoordContainer',


				 #add required EnsEMBL(core) adaptors here
				 #Should write/retrieve from efg not dna db
				 'Analysis'           => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
				 "MetaContainer"      => "Bio::EnsEMBL::DBSQL::MetaContainer",
				);
	
	return (\%pairs);
}	       
1;

