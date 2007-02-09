#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor - A database adaptor for fetching and
storing Funcgen ExperimentalChip objects.

=head1 SYNOPSIS

my $ec_a = $db->get_ExperimentalChipAdaptor();

my @ecs = @{$ec_a->fetch_all_by_Experiment($exp)};


=head1 DESCRIPTION

The ExperimentalChipAdaptor is a database adaptor for storing and retrieving
Funcgen ExperimentalChip objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalChipAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_all_by_experiment_dbID

  Arg [1]    : int - dbID of Experiment
  Example    : my @ecs = @{$ec_a->fetch_all_by_experiment_dbID($ac_dbid);
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_experiment_dbID {
    my $self = shift;
    my $e_dbid = shift;

	my ($ec_id, @results);

	throw("Must specify an experiemntal dbID") if(! $e_dbid);


	my $sth = $self->prepare("
		SELECT ec.experimental_chip_id
		FROM experimental_chip ec, experiment e
		WHERE ec.experiment_id = e.experiment_id
        AND e.experiment_id = $e_dbid
	");



	#can we do a generic fetch here?


	$sth->execute();


	while ($ec_id = $sth->fetchrow()){
	  #warn("got ec id $ec_id\n");
	  push @results, $self->fetch_by_dbID($ec_id);
	}

	return \@results;
}

=head2 fetch_contigsets_by_experiment_dbID

  Arg [1]    : int - dbID of Experiment
  Example    : my @track_sets = @{$ec_a->fetch_contigsets_experiment_dbID($ac_dbid);
  Description: returns a list of track sets which are each an arrayref containing 
               the track name as the first element, with the rest of the array being
	       the contiguous chip making up the trackset i.e. chips to be displayed on
	       the same track
  Returntype : Listref of mixed types
  Exceptions : None
  Caller     : General
  Status     : At Risk - hardcoded for v41 release

=cut

#contiguous/subsets/tracksets?  i.e. want to display on same track
sub fetch_contigsets_by_experiment_dbID {
    my $self = shift;
    my $e_dbid = shift;
    
    throw("deprecated, use ResultSet");

    my (@tracksets, @hack1);
    #my @hack1 = ("H3K9ac - Human Bone Osteosarcoma Epithelial Cells (U2OS)");
    #my @hack2 = ("H3kgac-2");
    #46092 + 46078; 46082 + 46075
    my $HeLa = "Human Epithelial Carcinoma Cells (HeLa)";
    my $GM06990 = "Human B-Lymphocyte Cells (GM06990)";

    #what are we going to return? arrayref to list of arrays of echips?
    #where do we get set name from?
    #first element should be set name
    #hashref to key = set name values = array of echips

    #differentiating purely on chip uid at present
    
    #This is currently a hack!!
    #Need ti implement contig_set_id in experimental_chip

    
    #HACKY HACKY!! NEED TO IMPLEMENT CONTIG SET IN DB AND API

    foreach my $echip (@{$self->fetch_all_by_experiment_dbID($e_dbid)}){
         
      if($self->db->species() =~ /homo/i){

	
	#Hacky set control needs handling in EC adaptor using the status tables
	#Also need to build/return name(Tissue/feature type name/desc) & FeatureType class > how to render
	

	if($echip->unique_id() eq "46092" || $echip->unique_id() eq "46078"){
	  push @hack1, "H3K9ac - Human Bone Osteosarcoma Epithelial Cells (U2OS)" if (scalar(@hack1) == 0);
	  push @hack1, $echip;
	}
	elsif($e_dbid == 7){
	  push @hack1, "H3K4me1 - $GM06990" if (scalar(@hack1) == 0);
	  push @hack1, $echip;
	}
	elsif($e_dbid == 6){
	  push @hack1, "H3K4me1 - $HeLa" if (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 8){
	  push @hack1, "H3K4me2 - $GM06990" if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 5){
	   push @hack1, "H3K4me2 - $HeLa" if  (scalar(@hack1) == 0); 
	   push @hack1, $echip;
	}
	elsif($e_dbid == 4){
	  push @hack1, "H3K4me3 - $HeLa" if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 9){
	  push @hack1, "H3ac - $GM06990" if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 3){
	  push @hack1, "H3ac - $HeLa"  if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
       	}
	elsif($e_dbid == 10){
	  push @hack1, "H4ac - $GM06990"  if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 2){
	  push @hack1, "H4ac - $HeLa" if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}
	elsif($e_dbid == 11){
	  push @hack1, "H3K4me3 - $GM06990" if  (scalar(@hack1) == 0); 
	  push @hack1, $echip;
	}

      }elsif($self->db->species() =~ /mus/i){
	
	next if ($echip->unique_id() != "48317");
	#&& $echip->unique_id() != "48316" &&
	#  $echip->unique_id() != "48320" && $echip->unique_id() != "65797");				    
	my @tmp = ("H3K4me3 - Mouse embyronic fibroblast (MEFf)", $echip);
	push @tracksets, \@tmp;
      }
      else{
	warn "No ExperimentalChip set hacks for species other than human or mouse";
      }
    }

    if($self->db->species() =~ /homo/i){
      @tracksets = (\@hack1);#, \@hack2);
    }


    return \@tracksets;
}

=head2 fetch_by_unique_and_experiment_id

  Arg [2]    : int - unique_id
  Arg [1]    : int - dbID of Experiment
  Example    : my $ec = ec_a->fetch_by_unique_and_experiment_id($c_uid, $exp_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_unique_and_experiment_id {
  my ($self, $c_uid, $e_dbid) = @_;
    
  throw("Must provide and unique_id and and experiment_id") if(! $c_uid || ! $e_dbid);

  my $sth = $self->prepare("
		SELECT ec.experimental_chip_id
		FROM experimental_chip ec
		WHERE ec.unique_id ='$c_uid'
        AND ec.experiment_id = $e_dbid
	");
  

  $sth->execute();
  my ($ec_id) = $sth->fetchrow();
  
	
  return $self->fetch_by_dbID($ec_id) if $ec_id;
}



=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return ['experimental_chip', 'ec'];
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw( ec.experimental_chip_id  ec.unique_id 
		   ec.experiment_id         ec.array_chip_id 
		   ec.feature_type_id       ec.cell_type_id
		   ec.replicate );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
	
  my (@result, $ec_id, $c_uid, $exp_id, $ac_id, $ftype_id, $ctype_id, $rep, $ftype, $ctype);

  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();
  
  
  $sth->bind_columns(\$ec_id, \$c_uid, \$exp_id, \$ac_id, \$ftype_id, \$ctype_id, \$rep);
  
  while ( $sth->fetch() ) {

    $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
    $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;
    
    my $array = Bio::EnsEMBL::Funcgen::ExperimentalChip->new(
							     -dbID           => $ec_id,
							     -unique_id      => $c_uid,
							     -experiment_id  => $exp_id,
							     -array_chip_id  => $ac_id,
							     -feature_type   => $ftype,
							     -cell_type      => $ctype,
							     -replicate      => $rep,
							     -adaptor        => $self,
							    );
	  
    push @result, $array;
    
  }
  return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Example    : $oaa->store($ec1, $ec2, $ec3);
  Description: Stores given ExperimentalChip objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : ARRAYREF
  Exceptions : Throws if passed non-ExperimentalChip arg or if ExperimentalChip already stored but arg has no dbID
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;
  
  my ($sarray);
  
  my $sth = $self->prepare("
			INSERT INTO experimental_chip
			(unique_id, experiment_id, array_chip_id, feature_type_id, cell_type_id, replicate)
			VALUES (?, ?, ?, ?, ?, ?)");
  
    
  
  foreach my $ec (@args) {
    throw('Can only store ExperimentalChip objects') if ( ! $ec->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip') );
    
    if (!( $ec->dbID() && $ec->adaptor() == $self )){
      
      my $s_ec = $self->fetch_by_unique_and_experiment_id($ec->unique_id(), $ec->experiment_id());
      throw("ExperimentalChip already exists in the database with dbID:".$s_ec->dbID().
	    "\nTo reuse/update this ExperimentalChip you must retrieve it using the ExperimentalChipAdaptor".
	    "\nMaybe you want to use the -recover option?") if $s_ec;
      
      my $ftype_id = (defined $ec->feature_type()) ? $ec->feature_type->dbID() : undef;
      my $ctype_id = (defined $ec->cell_type()) ? $ec->cell_type->dbID() : undef;
      
      $sth->bind_param(1, $ec->unique_id(),      SQL_VARCHAR);
      $sth->bind_param(2, $ec->experiment_id(),  SQL_VARCHAR);
      $sth->bind_param(3, $ec->array_chip_id(),  SQL_VARCHAR);
      $sth->bind_param(4, $ftype_id,             SQL_INTEGER);
      $sth->bind_param(4, $ctype_id,             SQL_INTEGER);
      $sth->bind_param(4, $ec->replicate(),      SQL_VARCHAR);
      
      $sth->execute();
      my $dbID = $sth->{'mysql_insertid'};
      $ec->dbID($dbID);
      $ec->adaptor($self);
      
      #}
      #else{
      #	  $ec = $s_ec;
      
      #my @states = @{$self->db->fetch_all_states('experimental_chip', $ec->dbID())};
      #	  my @states = @{$self->db->get_StatusAdaptor->fetch_all_states($ec)};
      #	  warn("Using previously stored ExperimentalChip (".$ec->unique_id().") with states\t@states\n");
      #  }
    }else{
      #assume we want to update the states
      warn('You may want to use $exp_chip->adaptor->store_states($exp_chip)');
      $self->store_states($ec);
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
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;	
    return $self->_list_dbIDs('experimental_chip');
}



1;

