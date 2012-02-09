#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor
#

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

Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor - A database adaptor for fetching and
storing Funcgen FeatureType objects.

=head1 SYNOPSIS

my $ft_adaptor = $db->get_FeatureTypeAdaptor();

my $feature_type = $ft_adaptor->fetch_by_name("H3K4me3");


=head1 DESCRIPTION

The FeatureTypeAdaptor is a database adaptor for storing and retrieving
Funcgen FeatureType objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::FeatureType

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

#Exported from BaseAdaptor
$true_tables{feature_type} = [['feature_type', 'ft']];
#deref so we don't modify the true tables
@{$tables{feature_type}}   = @{$true_tables{feature_type}};



#Regulatory evidence feature type information

my $core_label = 'DNase1 & TFBS';
my $attr_label = 'Hists & Pols';

my %regulatory_evidence_labels = 
  (
   $core_label => {
				   name      => 'Open chromatin & TFBS',
				   long_name => 'Open chromatin & Transcription factor binding sites',
				   label     => $core_label,
				  },
   
   
   $attr_label  => {
					name      => 'Histones & polymerases',
					long_name => 'Histone modifications & RNA polymerases',
					label     => $attr_label,
				   }
  );


my %regulatory_evidence_info = 
  (
   'Transcription Factor'         => $regulatory_evidence_labels{$core_label},
   'Transcription Factor Complex' => $regulatory_evidence_labels{$core_label},
   'Open Chromatin'               => $regulatory_evidence_labels{$core_label},
   'Polymerase'                   => $regulatory_evidence_labels{$attr_label},
   'Histone'                      => $regulatory_evidence_labels{$attr_label},
  );

#


=head2 fetch_by_name

  Arg [1]    : string - name of FeatureType
  Arg [2]    : optional string - class of FeatureType
  Arg [3]    : optional Bio::EnsEMBL::Analysis - Analysis used to generate FeatureType
  Example    : my $ft = $ft_adaptor->fetch_by_name('H3K4me2');
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType object (or ARRAY if called in ARRAY context)
  Exceptions : Throws if more than one FeatureType for a given name found.
               Throws if Analysis is defined but not valid.
  Caller     : General
  Status     : At risk

=cut

#Should really change this to fetch_all_by_name
#incorporating that functionality here will conditioanlly change the return type!!


sub fetch_by_name{
  my ($self, $name, $class, $analysis) = @_;

  throw("Must specify a FeatureType name") if(! $name);

  my $constraint = ' name = ? ';

  $constraint .= ' AND class = ? ' if $class;

  if($analysis){
	
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
	$constraint .= ' AND analysis_id = ? ';
  }



  $self->bind_param_generic_fetch($name,           SQL_VARCHAR);
  $self->bind_param_generic_fetch($class,          SQL_VARCHAR) if $class;
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER) if $analysis;
  

  my @fts = @{$self->generic_fetch($constraint)};
  

  #This can happen if using a redundant name between classes or analyses
  #remove?
  if( wantarray && (scalar @fts >1) ){
    $class ||= "";
    throw("Found more than one FeatureType:$class $name");
  }

  return (wantarray) ? @fts : $fts[0];
}


=head2 fetch_all_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Analysis
  Example    : my @fts = @{$ft_adaptor->fetch_all_by_Analysis($analysis);
  Description: Fetches all FeatureTypes for a given Analysis.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if Analysis not valid
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_Analysis{
  my ($self, $analysis) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
  my $constraint = ' analysis_id = ? ';

  #Use bind param method to avoid injection
  $self->bind_param_generic_fetch($analysis->dbID, SQL_INTEGER);

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_class

  Arg [1]    : string - class of FeatureType
  Example    : my @fts = @{$ft_adaptor->fetch_all_by_class('Histone')};
  Description: Fetches all FeatureTypes of a given class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if class not defined
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





sub fetch_all_by_associated_SetFeature{
  my ($self, $sfeat) = @_;
  deprecate('Please use the more generic fetch_all_by_association method');
  return $self->fetch_all_by_association($sfeat);
}

=head2 fetch_all_by_association

  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable
  Example    : my $assoc_ftypes = $ft_adaptor->fetch_all_by_association($ext_feature);
  Description: Fetches all associated FeatureTypes for a given Storable. 
               Note: Where appropriate, the main FeatureType for a Storable is 
               accessible via the feature_type method.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::FeatureType objects
  Exceptions : Throws if arg is not valid or stored
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_association{
  my ($self, $storable) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Storable', $storable);
  
  push @{$tables{feature_type}}, ['associated_feature_type', 'aft'];

  my $table_name = $storable->adaptor->_main_table->[0];

  my $constraint = 'aft.feature_type_id=ft.feature_type_id AND aft.table_name="'.$table_name.
	'" AND aft.table_id='.$storable->dbID;

  my $feature_types =  $self->generic_fetch($constraint);
  #Reset tables
  @{$tables{feature_type}} = @{$true_tables{feature_type}}; 

  return $feature_types;
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
	
  return @{$tables{feature_type}};
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
	
  return qw( ft.feature_type_id ft.name ft.class ft.analysis_id ft.description ft.so_accession ft.so_name);
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
	
	my (@result, $ft_id, $name, $class, $anal_id, $desc, $so_acc, $so_name, %analysis_hash);
	my $anal_a = $self->db->get_AnalysisAdaptor;
	
	$sth->bind_columns(\$ft_id, \$name, \$class, \$anal_id, \$desc, \$so_acc, \$so_name);
	
	$analysis_hash{0} = undef;

	  while ( $sth->fetch() ) {

	  $anal_id ||= 0;
	  
	  if (! exists $analysis_hash{$anal_id}){
		$analysis_hash{$anal_id} = $anal_a->fetch_by_dbID($anal_id);
	  }
	  
	  my $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
		(
		 -dbID         => $ft_id,
		 -NAME         => $name,
		 -CLASS        => $class,
		 -ANALYSIS     => $analysis_hash{$anal_id},
		 -DESCRIPTION  => $desc,
		 -SO_ACCESSION => $so_acc,
		 -SO_NAME      => $so_name,
		 -ADAPTOR      => $self
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
  Exceptions : Throws if FeatureType not valid
               Throws if Analysis defined but not valid
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;
  
  
  my $sth = $self->prepare("
			INSERT INTO feature_type
			(name, class, analysis_id, description, so_accession, so_name)
			VALUES (?, ?, ?, ?, ?, ?)");
    
  
  
  foreach my $ft (@args) {

    if ( ! (ref($ft) && $ft->isa('Bio::EnsEMBL::Funcgen::FeatureType') )) {
      throw('Can only store FeatureType objects, skipping $ft');
    }

    my $anal_id;
	my $analysis = $ft->analysis;
	
	if($analysis){	
	  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis);
	  $anal_id = $analysis->dbID;
	}



    if (! ( $ft->dbID && 
			( $ft->adaptor == $self) )
	   ){
      
      #Check for previously stored FeatureType
      my $s_ft = $self->fetch_by_name($ft->name, $ft->class, $ft->analysis);
	
      if(! $s_ft){
		$sth->bind_param(1, $ft->name,          SQL_VARCHAR);
		$sth->bind_param(2, $ft->class,         SQL_VARCHAR);
		$sth->bind_param(3, $anal_id,           SQL_INTEGER);
		$sth->bind_param(4, $ft->description,   SQL_VARCHAR);
		$sth->bind_param(5, $ft->so_accession,  SQL_VARCHAR);
		$sth->bind_param(6, $ft->so_name,       SQL_VARCHAR);
		
		
		$sth->execute();
		my $dbID = $sth->{'mysql_insertid'};
		$ft->dbID($dbID);
		$ft->adaptor($self);
      }
      else{
		$ft = $s_ft;
		warn("Using previously stored FeatureType:\t".$ft->name."\n"); 
      }
    }
  }

  return \@args;
}

=head2 list_regulatory_evidence_classes

  Args       : None
  Example    : 
  Description: 
  Returntype : Array of Strings
  Exceptions : None
  Caller     : web code
  Status     : At risk - remove in favour of get_regulatory_evidence_info

=cut

sub list_regulatory_evidence_classes {
    my ($self) = @_;
	
	#change this to return the whole hash

    return keys(%regulatory_evidence_info);
}

=head2 get_regulatory_evidence_labels

  Args       : None
  Example    : 
  Description: 
  Returntype : HASREF
  Exceptions : None
  Caller     : web code
  Status     : At risk

=cut

sub get_regulatory_evidence_labels{
  return \%regulatory_evidence_labels;
}

=head2 get_regulatory_evidence_info

  Args       : None
  Example    : 
  Description: Returns all regulatory evidence info keyed on feature type class
  Returntype : HASHREF
  Exceptions : None
  Caller     : web code
  Status     : At risk

=cut

sub get_regulatory_evidence_info{
  return \%regulatory_evidence_info;
}


#Deprecated

#list_dbIDs now uses inherited BaseAdaptor::list_dbIDs


1;

