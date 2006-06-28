
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

@ISA = qw(Bio::EnsEMBL::Funcgen::Helper);# Bio::EnsEMBL::GroupDefs);

use Bio::EnsEMBL::DBSQL::DBConnection;
#use Bio::EnsEMBL::GroupDefs;#Here on pass from Experiment?
use Bio::EnsEMBL::Funcgen::Helper;


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

sub new{
    my ($caller, %args) = @_;

    my ($self, %attrdata, $attrname, $argname);
    my $class = ref($caller) || $caller;

	#Create object from parent class
	$self = $class->SUPER::new(%args);

	#change this to re-arrange?
    # objects private data and default values
    %attrdata = (
				 _dbname => "efg_test",
				 _host   => "localhost",
				 _port   => undef,#default handled in dbc
				 _user   => "ensadmin",
				 _pass   => undef,
				 _dbc    => undef,#Pass existing connection to avoid multiple connections?
				);

    # set each class attribute using passed value or default value
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }

	
	#Do we need these? Yes!, but make all above undef(port?) and use GroupDefs to define

	foreach my $tmp("dbname", "user", "host"){
		$self->throw("Mandatory arg $tmp not been defined") if (! defined $self->{"_${tmp}"});
	}
	
	if(! $self->dbc()){
		$self->dbc(new Bio::EnsEMBL::DBSQL::DBConnection(
														 -user   => $self->user(),
														 -host   => $self->host(),
														 -dbname => $self->dbname(),
														 -host   => $self->host(),
														 -pass   => $self->pass(),
														));
	}



	$self->debug(2,"DBAdaptor class instance created.");
	$self->debug_hash(3, \$self);

    return ($self);
}


sub dbc{
	my ($self, $dbc) = @_;

	if(defined $dbc){

		if(! $dbc->isa('Bio::EnsEMBL::DBSQL::DBConnection')){
			$self->throw("$dbc is no a DBConnection\n");
		}
	
		$self->{'_dbc'} = $dbc;
	}

	return 	$self->{'_dbc'};
}


#genericise?
sub user{
	my ($self, $user) = @_;
	$self->{'_user'} = $user if($user);
	return 	$self->{'_user'};
}

sub pass{
	my ($self, $pass) = @_;
	$self->{'_pass'} = $pass if($pass);
	return 	$self->{'_pass'};
}


sub dbname{
	my ($self, $dbname) = @_;
	$self->{'_dbname'} = $dbname if($dbname);
	return 	$self->{'_dbname'};
}

sub port{
	my ($self, $port) = @_;
	$self->{'_port'} = $port if($port);
	return 	$self->{'_port'};
}

sub host{
	my ($self, $host) = @_;
	$self->{'_host'} = $host if($host);
	return 	$self->{'_host'};
}



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

	$self->log("Loading $table data from $file");
	my $sql = "load data infile \"$file\" into table $table";
	$self->dbc->do($sql);
	$self->log("Finished loading $table data");
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
		$self->log("$table $dbid already registered, validating");
		$sql = "select * from $table where ${table}_id = \"$dbid\"";
		@db_entry = $self->dbc->db_handle->selectrow_array($sql);
		shift @db_entry;

		for my $i(0..$#db_entry){

			if($db_entry[$i] ne $data[$i]){
				$valid = 0;
				#should really get field names here too, and store if doing repeat validations to reduce no. of queries
				$self->warn("Validating $table entry $dbid found the following mismatch:\nData\t$data[$i]\nDB\t$db_entry[$i]");
			}
		}
	}
	
	return $valid;#return value to allow caller to throw
}

1;

