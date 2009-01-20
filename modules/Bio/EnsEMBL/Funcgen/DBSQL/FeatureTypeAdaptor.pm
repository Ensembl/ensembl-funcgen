#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor - A database adaptor for fetching and
storing Funcgen FeatureType objects.

=head1 SYNOPSIS

my $ft_adaptor = $db->get_FeatureTypeAdaptor();

my $feature_type = $ft_adaptor->fetch_by_name("H3K4me3");


=head1 DESCRIPTION

The FeatureTypeAdaptor is a database adaptor for storing and retrieving
Funcgen FeatureType objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string - name of FeatureType
  Arg [1]    : optional string - class of FeatureType
  Example    : my $ft = $ft_adaptor->fetch_by_name('H3K4me2');
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType object
  Exceptions : Throws if more than one FeatureType for a given name found
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name{
  my ($self, $name, $class) = @_;

  throw("Must specify a FeatureType name") if(! $name);


  my $constraint = " name = ?";

  $constraint .= " AND class = ?" if $class;


  $self->bind_param_generic_fetch($name,   SQL_VARCHAR);
  $self->bind_param_generic_fetch($class,  SQL_VARCHAR) if $class;

  my @fts = @{$self->generic_fetch($constraint)};
  

  #This should never happen?
  if(scalar @fts >1){
    $class ||= "";
    throw("Found more than one FeatureType:$class $name");
  }


  return $fts[0];
}

=head2 fetch_all_by_class

  Arg [1]    : string - class of FeatureType
  Example    : my $ft = $ft_adaptor->fetch_all_by_class('Histone');
  Description: Fetches all FeatureTypes of a given class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_class{
  my ($self, $class) = @_;

  throw("Must specify a FeatureType class") if(! defined $class);


  my $constraint = " class = ? ";


  #Use bind param method to avoid injection
  $self->bind_param_generic_fetch($class, SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_associated_SetFeature

  Arg [1]    : Bio::EnsEMBL:SetFeature
  Example    : my $assoc_ftypes = $ft_adaptor->fetch_all_by_associated_SetFeature($ext_feature);
  Description: Fetches all associated FeatureTypes for a given SetFeature. Note: The main FeatureType for
               a SetFeature is accessible via the feature_type method.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if arg is not valid or stored
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_associated_SetFeature{
  my ($self, $sfeat) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::SetFeature', $sfeat);
  my $constraint;

  #Need to protect against SQL injection here?
  #Not so vital as this never takes user supplied data

  my $sql = 'SELECT feature_type_id from associated_feature_type where feature_table="'.$sfeat->feature_set->type.'" and feature_id='.$sfeat->dbID;



  my @ft_ids = map {$_ = "@$_"} @{$self->dbc->db_handle->selectall_arrayref($sql)};
  $constraint = ' feature_type_id in ('.join(',',@ft_ids).') ' if @ft_ids;

  return ($constraint) ? $self->generic_fetch($constraint) : [];
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
	
  return (
	  ['feature_type', 'ft'],
	 );
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
	
  return qw( ft.feature_type_id ft.name ft.class ft.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $ft_id, $name, $class, $desc);
	
	$sth->bind_columns(\$ft_id, \$name, \$class, \$desc);
	
	while ( $sth->fetch() ) {
	  my $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
							     -dbID        => $ft_id,
							     -NAME        => $name,
							     -CLASS       => $class,
							     -DESCRIPTION => $desc,
							     -ADAPTOR     => $self,
							    );
	  
	  push @result, $ftype;
	  
	}
	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::FeatureType objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores given Channel objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;
  
  
  my $sth = $self->prepare("
			INSERT INTO feature_type
			(name, class, description)
			VALUES (?, ?, ?)");
    
  
  
  foreach my $ft (@args) {
    if ( ! $ft->isa('Bio::EnsEMBL::Funcgen::FeatureType') ) {
      warning('Can only store FeatureType objects, skipping $ft');
      next;
    }
    
    if (!( $ft->dbID() && $ft->adaptor() == $self )){
      
      #Check for previously stored FeatureType
      my $s_ft = $self->fetch_by_name($ft->name(), $ft->class());
	
      if(! $s_ft){
	$sth->bind_param(1, $ft->name(),           SQL_VARCHAR);
	$sth->bind_param(2, $ft->class(),          SQL_VARCHAR);
	$sth->bind_param(3, $ft->description(),    SQL_VARCHAR);
	
	$sth->execute();
	my $dbID = $sth->{'mysql_insertid'};
	$ft->dbID($dbID);
	$ft->adaptor($self);
      }
      else{
		$ft = $s_ft;
		warn("Using previously stored FeatureType:\t".$ft->name()."\n"); 
      }
    }
  }

  return \@args;
}


=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$ec_a->list_dbIDs()};
  Description: Gets an array of internal IDs for all ExperimentalChip objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : At risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('feature_type');
}



1;

