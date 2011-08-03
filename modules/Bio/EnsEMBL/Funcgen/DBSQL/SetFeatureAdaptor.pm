#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor
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

Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor - A base database adaptor for SetFeature adaptors.
storing SetFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_SetFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The SetFeatureAdaptor is a base adaptor for all SetFeature adaptors.
e.g. SetFeature, AnnotatedFeature etc.  It provides common methods
across all feature types.



=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::SetFeature

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::SetFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA @EXPORT);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}}, '%tables', '%true_tables');
#required for child adaptor's store and _obj_from_sth methods

=head2 fetch_all_by_FeatureType_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) hashref - params hash
                                      associated => 1, #Also return feature which have the associated FeatureType
                                      logic_name => 'analysis.logic_name'
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice with main or associated FeatureType. 
               This is mainly used by external FeatureSets which can sometimes have more 
               than one associated FeatureType. NOTE: This is intended to work for FeatureTypes at the 
               feature level, not the more generic FeatureSet level FeatureTypes.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType_FeatureSets {
  my ($self, $ftype, $fsets, $params) = @_;

  if ($self->_feature_class eq 'annotated'){
	throw("This method is not appropriate for FeatureSets with feature_class 'annotated', please use standard fetch_all_by_FeatureSets or get_all_Features on an individual FeatureSet");


	#Why is this not appropriate? #As there is not feature_type_id in annotated_feature?
	

  }

  #How do we validate the parameters?
  #We can't catch typos as we don't know what might need passing on to other methods
  #Set all these so we can use them without exists
  $params->{'logic_name'} ||= undef;
  $params->{'associated'} ||= undef;
  
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

 
  #Need to restrict to the current default cs's
  #This requires separate query or query extension
  #Need to make @tables and @true_tables exported?
  #Or just attrs?
  #Set table to true_tables in new?
  my ($table_name, $table_syn) = @{$self->_main_table};

  my $constraint = $table_syn.'.feature_set_id = fs.feature_set_id AND '.
	$table_syn.'.feature_type_id='.$ftype->dbID.' AND '.$table_syn.".feature_set_id ".$self->_generate_feature_set_id_clause($fsets);

  #We should really pass the params hash on here
  $constraint = $self->_logic_name_to_constraint($constraint, $params->{logic_name});


  my @features = @{$self->generic_fetch($constraint)};


  #This is an interim solution and really needs changing!
  #Can we genericise this as a lazy loader function?
  #Isn't there already something liek this in core?
  if ($params->{'associated'}){

	for my $fset(@{$fsets}){
	  #We want to bring back features from the same set
	  #We are not bringing back different class of feature i.e. AnnotatedFeatures and ExternalFeatures
	  #Let's not catch this, but let it silently fail so peaople.

	  next if $table_name ne lc($fset->feature_class).'_feature';

	  my $sql = 'SELECT table_id from associated_feature_type where table_name="'.$table_name.'" and feature_type_id='.$ftype->dbID;

	  my @dbids = map $_ = "@$_", @{$self->dbc->db_handle->selectall_arrayref($sql)};

	  if(@dbids){
		$constraint = " $table_syn.${table_name}_id in (".join(',',@dbids).') ' if @dbids;
		push @features, @{$self->generic_fetch($constraint, $params->{logic_name})};
	  }
	}
  }

  return \@features;
  #return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_Feature_associated_feature_types

  Arg [1]    : Bio::EnsEMBL::SetFeature
  Arg [2]    : (optional) hashref - params hash, all entries optional
                                      -logic_name => 'analysis.logic_name'
                                      -include_direct_links => 1, #Also return feature which are linked by Feature->feature_type
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $sfa->fetch_all_by_Feature_associated_feature_types($feature);
  Description: Retrieves a list of all features linked via the associated FeatureTypes of 
               the given Feature in the same FeatureSet.  This is mainly used by external
               FeatureSets which can sometimes have more than one associated FeatureType.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::SetFeature objects
  Exceptions : Throws if SetFeature not stored and valid
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Feature_associated_feature_types {
  my ($self, $feature, $params) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::SetFeature', $feature);

  #We always want associated!
  #Otherwise it would be just a normal fetch_all_by_FeatureType_FeatureSets query 
  #Which is more efficient.


  #How do we validate the parameters?
  #We can't catch typos as we don't know what might need passing on to other methods
  #Set all these so we can use them without exists
  $params->{'include_direct_links'} ||= undef;
  $params->{'logic_name'}           ||= undef;

  my %dbIDs;#Use a hash to reduce dbID redundancy
  my ($table_name, $table_syn) = @{$self->_main_table};
  my $fset = $feature->feature_set;
  my ($sql, $constraint, @features);

  #Direct FeatureType
  if($params->{'include_direct_links'}){
	#Just grab dbIDs here rather than use generic fetch
	$sql = "SELECT ${table_name}_id from $table_name where feature_type_id=".$feature->feature_type->dbID.' and feature_set_id='.$fset->dbID;

	#This just sets each value to a key with an undef value
	map {$dbIDs{"@$_"} = undef} @{$self->dbc->db_handle->selectall_arrayref($sql)};
  }

 
    
  my @assoc_ftypes = @{$feature->associated_feature_types};

  if(@assoc_ftypes){

	my $ftype_ids = join(', ', (map $_->dbID, @assoc_ftypes));
  
	#Now we need to restrict the Features based on the FeatureSet of the original Feature, we could make this optional.
	$sql = "SELECT table_id from associated_feature_type aft, $table_name $table_syn where aft.table_name='".$fset->feature_class."_feature' and aft.feature_type_id in ($ftype_ids) and aft.table_id=${table_syn}.${table_name}_id and ${table_syn}.feature_set_id=".$fset->dbID;
	#This just sets each value to a key with an undef value
	
	map {$dbIDs{"@$_"} = undef} @{$self->dbc->db_handle->selectall_arrayref($sql)};
  }



  if(keys %dbIDs){
	$constraint = " $table_syn.${table_name}_id in (".join(',', keys %dbIDs).')';
	push @features, @{$self->generic_fetch($constraint, $params->{logic_name})};
  }


  return \@features;
}





=head2 fetch_all_by_Slice_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureType {
  my ($self, $slice, $type, $logic_name) = @_;
	
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $type);
  
  my $ft_id = $type->dbID();
  
  my $constraint = $self->_main_table->[1].".feature_set_id = fs.feature_set_id AND ".
	"fs.feature_type_id = '$ft_id'";

  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_by_FeatureSets

  Arg [1]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [2]    : optional - analysis.logic_name
  Example    : my $features = $set_feature_adaptor->fetch_all_by_FeatureSets(@fsets);
  Description: Retrieves a list of features specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk - Maybe not sensible to fetch all data in one call.

=cut

sub fetch_all_by_FeatureSets {
  my ($self, $fsets, $logic_name) = @_;
	
  my $constraint = $self->_main_table->[1].'.feature_set_id '.$self->_generate_feature_set_id_clause($fsets); 
  #could have individual logic_names for each annotated feature here?
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->generic_fetch($constraint);
}


=head2 _generate_feature_set_id_clause

  Arg [1]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Example    : my $fset_d_clause = $self->_generate_feature_set_id_clause($fsets);
  Description: Build feature_set id clause for FeatureSet methods
  Returntype : string
  Exceptions : Throws if FeatureSets are passed
               Throws if FeatureSet feature_class does not match adaptor feature_class
               Throws if FeatureSet is not valid 
  Caller     : Bio::EnsEMBL::DBSQL::SetFeatureAdaptor
  Status     : At Risk

=cut

sub _generate_feature_set_id_clause{
  my ($self, $fsets) = @_;

  my @fs_ids;

  if(! ( (ref($fsets) eq 'ARRAY') &&
		 scalar(@{$fsets}) > 0) ){
	throw('Must provide a list of Bio::EnsEMBL::FeatureSet objects');
  }

  my $fclass = $self->_feature_class;
  
  foreach my $fset (@{$fsets}) {
	
	if (! (ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"))){
	  throw('Not a FeatureSet object');
	}

	if($fset->feature_class ne $fclass){
	  throw('FeatureSet feature_class \''.$fclass.'\' does not match adaptor feature_class \''.$fset->feature_class.'\'');
	}

	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);

	push (@fs_ids, $fset->dbID());
  }

  return scalar(@fs_ids >1) ? 'IN ('.join(',', @fs_ids).')' : '= '.$fs_ids[0];
}



=head2 fetch_all_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [3]    : optional - analysis.logic_name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_FeatureSets($slice, \@fsets);
  Description: Retrieves a list of features on a given slice, specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureSets {
  my ($self, $slice, $fsets, $logic_name) = @_;
	
  my $constraint = $self->_main_table->[1].'.feature_set_id '.
	$self->_generate_feature_set_id_clause($fsets); 

  #could have individual logic_names for each annotated feature here?
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_Iterator_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Arrayref of Bio::EnsEMBL::FeatureSet objects
  Arg [3]    : optional - analysis.logic_name
  Arg [4]    : optional - iterator chunk length. Default is 1MB
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $iter = $ofa->fetch_Iterator_by_Slice_FeatureSets($slice, \@fsets);

           
  Description: Simple Iterator wrapper method for fetch_all_by_Slice_FeatureSets
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut


sub fetch_Iterator_by_Slice_FeatureSets{
  my ($self, $slice, $fsets, $logic_name, $chunk_length) = @_;

  return $self->fetch_Iterator_by_Slice_method
	($self->can('fetch_all_by_Slice_FeatureSets'),
	 [$slice, $fsets, $logic_name],
	 0,#Slice idx
	 $chunk_length #default is 1000000
	);
}




# Redefine BaseFeatureAdaptor method as analysis now abstracted to feature_set
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
# DO WE HAVE TO CALL THIS EXPLICITLY in this adaptor, or will BaseAdaptor use redefined method?


sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if (!$logic_name);

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if(! $an) {
	  #warn or throw?
	  warn("No analysis associated with logic_name $logic_name");
	  return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " fs.analysis_id = $an_id";
  return $constraint;
}



=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _default_where_clause {
  my $self = shift;
	
  return $self->_main_table->[1].'.feature_set_id = fs.feature_set_id';
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an GROUP/ORDER BY clause. 
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut


sub _final_clause {
  my $self = shift;
  #return '';
  return ' ORDER BY '.$self->_main_table->[1].'.seq_region_id, '.$self->_main_table->[1].'.seq_region_start';
}


#re-written for non-standard feature table

=head2 fetch_all_by_logic_name 

  Arg [1]     : string $logic_name the logic name of the analysis of features to obtain 
  Example     : $fs = $a->fetch_all_by_logic_name('foobar'); 
  Description : Returns an arrayref of features created from the database. only 
                features with an analysis of type $logic_name will be returned. 
  Returntype  : arrayref of Bio::EnsEMBL::SetFeatures
  Exceptions  : thrown if $logic_name not defined
  Caller      : General 
  Status      : At risk

=cut 

sub fetch_all_by_logic_name { 
  my $self = shift; 
  my $logic_name = shift || throw( "Need a logic_name" ); 

  my $constraint;
  my $an = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name); 
  $constraint = ' '.$self->_main_table->[1].'.feature_set_id=fs.feature_set_id and fs.analysis_id = '.$an->dbID() if($an);

  return (defined $constraint) ? $self->generic_fetch($constraint) : undef;
} 


=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all SetFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs($self->_main_table->[0]);
}


=head2 _feature_class

  Example    : if($self->feature_class ne $fset->feature_class){ throw('some error'); }
  Description: Convenience method to retrieve the feature class for this adaptor
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub _feature_class{
  my $self = shift;

  #use the first word of the table name as the class
  my $fclass;
  ($fclass = $self->_main_table->[0]) =~ s/_.*//;#use the first word of the table name as the class

  return $fclass;
}

1;
