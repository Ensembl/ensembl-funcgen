#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor - A database adaptor for fetching and
storing Funcgen OligoArray objects.

=head1 SYNOPSIS

my $oaa = $db->get_OligoArrayAdaptor();

my $array = $oaa->fetch_by_name('Array-1');
my @arrays = @{$oaa->fetch_all()};

=head1 DESCRIPTION

The OligoArrayAdaptor is a database adaptor for storing and retrieving
Funcgen OligoArray objects.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
OligoArrayAdaptor modules written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::OligoArrayAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::OligoArray;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_array_chip_dbID

  Arg [1]    : int - dbID of array_chip
  Example    : my $array = $oaa->fetch_by_array_chip_dbID($ac_dbid);
  Description: Retrieves a named OligoArray object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::OligoArray
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_array_chip_dbID {
    my $self = shift;
    my $ac_dbid = shift;
	my $sth = $self->prepare("
		SELECT a.array_id
		FROM array a, array_chip ac
		WHERE a.array_id = ac.array_id
		AND ac.array_chip_id = $ac_dbid
	");


	$sth->execute();
	my ($array_id) = $sth->fetchrow();

	return $self->fetch_by_dbID($array_id);
}



=head2 fetch_by_name

  Arg [1]    : string - name of an array
  Example    : my $array = $oaa->fetch_by_name('Array-1');
  Description: Retrieves a named OligoArray object from the database.
  Returntype : Bio::EnsEMBL::Funcgen::OligoArray
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;
    
    my $result = $self->generic_fetch("a.name = '$name'");
	
    if (scalar @$result > 1) {
		warning("Array $name is not unique in the database, but only one result has been returned");
    } 


    #should have fetch by name vendor, to provide uniqueness?
    #should check for this on import!
    return $result->[0];
}

#=head2 fetch_all_by_type
#
#  Arg [1]    : List of strings - type(s) (e.g. AFFY or OLIGO)
#  Example    : my @arrays = @{$oaa->fetch_all_by_type('OLIGO')};
#  Description: Fetch all arrays of a particular type.
#  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoArray objects
#  Exceptions : Throws if no type is provided
#  Caller     : General
#  Status     : Medium Risk

#=cut

#sub fetch_all_by_type {
#	my ($self, @types) = @_;
	
#	throw('Need type as parameter') if !@types;
	
#	my $constraint;
#	if (scalar @types == 1) {
#		$constraint = qq( oa.type = '$types[0]' );
#	} else {
#		$constraint = join q(','), @types;
#		$constraint = qq( oa.type IN ('$constraint') );
#	}

#	return $self->generic_fetch($constraint);
#}

=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoArray - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty OligoArray objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::OligoArray getters
  Status     : Medium Risk

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _tables {
	my $self = shift;
	
	return ['array', 'a'];
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw( a.array_id a.name a.format a.vendor a.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoArray objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $array_id, $name, $format, $vendor, $description);
	
	$sth->bind_columns(\$array_id, \$name, \$format, \$vendor, \$description);
	
	while ( $sth->fetch() ) {
		my $array = Bio::EnsEMBL::Funcgen::OligoArray->new(
														   -dbID        => $array_id,
														   -adaptor     => $self,
														   -name        => $name,
														   -format      => $format,
														   -vendor      => $vendor,
														   -description => $description,
														  );

		push @result, $array;

		#if ($parent_id) {
		#	my $parent_array = Bio::EnsEMBL::Funcgen::OligoArray->new(
		#		-dbID    => $parent_id,
		#		-adaptor => $self,
		#	);
		#	$array->superset($parent_array);
		#}
	}
	return \@result;
}


=head2 _fetch_array_chips_by_array_dbID

  Arg [1]    : int - dbID or array
  Example    : None
  Description: Retrieves array_chip hashes guven an array id
  Returntype : Listref of array_chip hashes
  Exceptions : Throws if no array dbID specified
  Caller     : Internal
  Status     : At risk - rplace with ArrayChipAdaptor

=cut

sub _fetch_array_chips_by_array_dbID {
	my ($self, $array_dbID) = @_;

	throw("Must specifiy an array_id to retrieve array_chips") if (! $array_dbID);
	my ($array_chip_id, $design_id, $name, %ac_tmp);

	my $sth = $self->prepare("select array_chip_id, design_id, name from array_chip where array_id = $array_dbID");	
	$sth->execute();

	$sth->bind_columns(\$array_chip_id, \$design_id, \$name);
	
	while ( $sth->fetch() ) {
		$ac_tmp{$design_id} = {(
					array_id => $array_dbID,
					dbID => $array_chip_id,
					name => $name,
				       )};
	}

	return \%ac_tmp;
}





=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::OligoArray objects
  Example    : $oaa->store($array1, $array2, $array3);
  Description: Stores given OligoArray objects in the database. This
               method checks for arrays previously stored and updates 
               and new array_chips accordingly.
  Returntype : Listref of stored OligoArray objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


#This works slightly differently as arary_chip are not stored as objects, 
#yet we need to retrieve a dbID for the array before we know about all the array_chips

sub store {
  my $self = shift;
  my @args = @_;
  
  my ($sarray);
  
  my $sth = $self->prepare("
			INSERT INTO array
			(name, format, vendor, description)
			VALUES (?, ?, ?, ?)");
  
    
  foreach my $array (@args) {
    if ( !$array->isa('Bio::EnsEMBL::Funcgen::OligoArray') ) {
      warning('Can only store OligoArray objects, skipping $array');
      next;
    }
    
    # Has array already been stored?
    # What about array_chips? Still need to check they've been imported
    if (!( $array->dbID() && $array->adaptor() == $self )){
      #try and fetch array here and set to array if valid
      $sarray = $self->fetch_by_name($array->name());#this should be name_vendor?
      
      if( ! $sarray){
	warn("Storing array\n");
	
	$sth->bind_param(1, $array->name(),         SQL_VARCHAR);
	$sth->bind_param(2, $array->format(),       SQL_VARCHAR);
	$sth->bind_param(3, $array->vendor(),       SQL_VARCHAR);
	$sth->bind_param(4, $array->description(),  SQL_VARCHAR);
	
	$sth->execute();
	my $dbID = $sth->{'mysql_insertid'};
	$array->dbID($dbID);
	$array->adaptor($self);
      }
      else{
	warn("Array already stored, using previously stored array");# validating array_chips\n");
	$array = $sarray;
      }
    }
    
    #$array = $self->store_array_chips($array, $sarray);
    
    #need to make sure we have the full object here, so query again or can we use array or sarray?
    #there may be problems with repopulating array with full achip complement if experiment only uses some of the achips?
    #but we're reging the achips in the importer, so this shouldn't matter?
    
  }
  
  
  #will this be the original or the stored objects?
  
  return \@args;
  #return $array;
}



#This method takes an array, and optinally a stored sarray else it retrieves one
#Missing array_chips in sarray are stored and set in array
#then the array size is updated if the array_chips hash is larger the the size attribute
#do we need this size attribute?  Can we not remove it and do a dynamic keys on the array_chips hash?


#sub store_array_chip{
#  my ($self, $array) = @_; #, $sarray) = @_;
#  
#  
#  #need to check for array here?
#  
#  my $ac_sth = $self->prepare("INSERT INTO array_chip
#                                 (design_id, array_id, name)
#                                 VALUES (?, ?, ?)");
#  
#  #$sarray  = $self->fetch_by_name($array->name()) if(! $sarray);
#  
#  foreach my $design_id(keys %{$array->array_chips()}){
#    my $sdesign_id = 0;
#    
#    if($sarray && $sarray->array_chips()){
#      ($sdesign_id) = grep(/$design_id/, keys %{$sarray->array_chips});
#    }
#    
#    
#    #do we need to set some flag here and in _obj_from_sth
#    #registered flag, true for retrieved, false for just stored achips, 
#    #this will be used by importer to run probe import etc
#    #Use status?
#    
#    #Change these warns to logs?
#
#    if(! $sdesign_id){
#      warn("Storing array chip:\t".$array->array_chips->{$design_id}{'name'}."\n");
#      
#      #if this is throwing because of no dbID, this is because you've deleted the array_chips and not the array
#
#      $ac_sth->bind_param(1, $design_id,                                SQL_VARCHAR);
#      $ac_sth->bind_param(2, $array->dbID(),                            SQL_INTEGER);
#      $ac_sth->bind_param(3, $array->array_chips->{$design_id}{'name'}, SQL_VARCHAR);
#      $ac_sth->execute();
#   
#      $array->array_chips->{$design_id}{'dbID'} = $ac_sth->{'mysql_insertid'};
#      #$self->db->set_status('array_chip', $array->array_chips->{$design_id}{'dbID'}, 'STORED');
#      $sarray->add_array_chip($design_id, $array->array_chips->{$design_id}) if ($sarray);
#    }else{
#      warn("Array chip already stored:\t".$array->array_chips->{$design_id}{'name'}."\n");
#
#    }
#    
#  }
#  
#  #Set new achips
#  #if($sarray){
#  #  $array->array_chips($sarray->array_chips()) if ($sarray);#Can we not just set array to sarray?
#  #}
#  
#  return $sarray || $array;
#}


=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('array');
}


#New Funcgen methods
#fetch_all_by_group?
#fetch_by_channel_id?
#fetch_by_chip_id?
#fetch_by_probe_id
#fetch_by_probe_set_id


1;

