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

=cut

package Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use vars qw(@ISA);



@ISA = ('Bio::EnsEMBL::Funcgen::Utils::Helper');

# Base functionality for external_feature parsers

#Make this inherit from Helper?
#Then change all the prints to logs

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  #validate and set type, analysis and feature_set here
  my ($type, $db, $clobber, $archive, $import_fsets) = rearrange(['TYPE', 'DB', 'CLOBBER', 'ARCHIVE', 'IMPORT_SETS'], @_);
  
  throw('You must define a type of external_feature to import') if(! defined $type);

  if (! ($db && ref($db) &&
		 $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  throw('You can only specify either -clobber or -archive, but not both') if($clobber && $archive);

  $self->{'display_name_cache'} = {};
  $self->{'db'} = $db;
  $self->{type} = $type;
  $self->{'clobber'} = $clobber if defined $clobber;
  $self->{'archive'} = $archive if defined $archive;


  #This is not fully implemented yet and need to be validated against the config feature_set
  #pass something like set1,set2 and split and validate each.
  #Or do this in the calling script?

  throw('-import_sets not fully implemented yet') if defined $import_fsets;
  $self->{'import_sets'} = (defined $import_fsets) ? @{$import_fsets} : undef;
  
  $self->log("Parsing and loading $type ExternalFeatures");

  return $self;

}


=head2 db

  Args       : None
  Example    : my $feature_set_adaptor = $seld->db->get_FeatureSetAdaptor
  Description: Getter for the DBAdaptor.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub db{
  my $self = shift;

  return $self->{'db'};
}

=head2 import_sets

  Args       : None
  Example    : foreach my $import_set_name(@{$self->import_sets}){ ... do the import ... }
  Description: Getter for the list of import feature set names, defaults to all in parser config.
  Returntype : Arrayref of import feature_set names
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub import_sets{
  my $self = shift;

  return $self->{'import_sets'} || [keys %{$self->{'feature_sets'}}];
}


=head2 set_feature_sets

  Args       : None
  Example    : $self->set_feature_sets;
  Description: Imports feature sets defined by import_sets.
  Returntype : None
  Exceptions : Throws if feature set already present and clobber or archive not set
  Caller     : General
  Status     : Medium Risk

=cut

#This is done after validate and store feature_types

sub set_feature_sets{
  my $self = shift;

  throw('Must provide a set feature_set config hash') if ! defined $self->{'feature_sets'};


  my $fset_adaptor = $self->db->get_FeatureSetAdaptor;
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
  
  foreach my $fset_name(@{$self->import_sets}){

	$self->log("Defining FeatureSet:\t$fset_name");  
	my $fset = $fset_adaptor->fetch_by_name($fset_name);

	#we don't need data sets for external_feature sets!

	if(defined $fset){
	  $self->log("Found previous FeatureSet $fset_name");

	  if($self->{'clobber'}){

		$self->rollback_FeatureSet($fset);#Need to pass \@slices here?
	  }
	  elsif($self->{'archive'}){
		my $archive_fset =  $fset_adaptor->fetch_by_name($fset_name."_v".$self->{'archive'});

		if(defined $archive_fset){
		  throw("You are trying to create an archive external feature_set which already exists:\t${fset_name}_v".$self->{archive});
		}

		my $sql = "UPDATE feature_set set name='$fset_name}_v".$self->{archive}."' where name='$fset_name'";
		$self->db->dbc->do($sql);
		undef $fset;
	  }else{
		throw("You are trying to create an external feature_set which already exists:\t$fset_name\nMaybe to want to clobber or archive?");
	  }
	}

	if(! defined $fset){
	  #don't need to use RNAFeatureType here as this is the setwide generic feature_type
	  #or do we have separate tables for external_feature and external_rna_feature?

	  #validate analysis first
	  my $analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'});

	  if(! defined $analysis){
		
		$self->log('Analysis '.$self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'}.
		  " not found, storing from config hash");		
		$analysis_adaptor->store(Bio::EnsEMBL::Analysis->new(%{$self->{'feature_sets'}{$fset_name}{analysis}}));
		$analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{analysis}{-logic_name});
		warn "fetched sotre analysis $analysis";
	  }

	  #replace hash config with object
	  $self->{'feature_sets'}{$fset_name}{'analysis'} = $analysis;

	  warn "analysis is $analysis ".$analysis->dbID;


	  my $display_name = (exists $self->{'feature_sets'}{$fset_name}{'display_label'}) ? $self->{'feature_sets'}{$fset_name}{'display_label'} : $fset_name;

	  

	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													 -name         => $fset_name,
													 -feature_class=> 'external',
													 -analysis     => $self->{'feature_sets'}{$fset_name}{'analysis'},
													 -feature_type => ${$self->{'feature_sets'}{$fset_name}{'feature_type'}},
													 -display_label => $display_name,
													);

	  ($fset) = @{$self->db->get_FeatureSetAdaptor->store($fset)};
	}

	#Now replace config hash with object
	$self->{feature_sets}{$fset_name} = $fset;
  }

  return;
}


=head2 validate_and_store_feature_types

  Args       : None
  Example    : $self->validate_and_store_feature_types;
  Description: Imports feature types defined by import_sets.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Medium Risk - Move this to (Base)Importer.pm

=cut

sub validate_and_store_feature_types{
  my $self = shift;

  # LBL enhancers have both positive and negative varieties

  #feature type class/naming is logical but not intuitive
  #+-----------------+-------------------------+----------+----------------------------------------------------+
  #| feature_type_id | name                    | class    | description                                        |
  #+-----------------+-------------------------+----------+----------------------------------------------------+
  #|          398680 | VISTA Enhancer          | Enhancer | Enhancer identified by positive VISTA assay        |
  #|          398681 | VISTA Target - Negative | Region   | Enhancer negative region identified by VISTA assay |
  #+-----------------+-------------------------+----------+----------------------------------------------------+


  #if (lc($type) eq 'VISTA') {
  #  return (validate_type($db_adaptor, 'VISTA Enhancer') && validate_type($db_adaptor, 'VISTA Target - Negative'));
  #}

  #my $sth = $self->db->dbc->prepare("SELECT analysis_id FROM analysis WHERE logic_name=?");

  

  #remove lc as mysql doesn't care about case
  #$sth->execute($type);
  #if ($sth->fetchrow_array()) {
  #  print "Type $type is valid\n";
  #} else {
  #  print "Type $type is not valid - is there an entry for $type in the analysis table?\n";
  #  return 0;
  #}

  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;

  foreach my $import_set(@{$self->import_sets}){

	my $ftype_config = ${$self->{'feature_sets'}{$import_set}{'feature_type'}};
	my $ftype = $ftype_adaptor->fetch_by_name($ftype_config->{'name'});

	$self->log("Validating $import_set FeatureType:\t".$ftype_config->{'name'});

	if(! defined $ftype){
	  $self->log("FeatureType '".$ftype_config->{'name'}."' for external feature_set ".$self->{'type'}." not present");
	  $self->log("Storing using type hash definitions");
	
	  $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
													   -name => $ftype_config->{'name'},
													   -class => $ftype_config->{'class'},
													   -description => $ftype_config->{'description'},
													  );
	  ($ftype) = @{$ftype_adaptor->store($ftype)};
	}

	#Replace hash config with object
	$self->{'feature_types'}{$ftype_config->{'name'}} = $ftype;
  }

  return;
}



=head2 project_feature

  Args [0]   : Bio::EnsEMBL::Feature
  Args [1]   : string - Assembly e.g. NCBI37
  Example    : $self->project($feature, $new_assembly);
  Description: Projects a feature to a new assembly via the AssemblyMapper
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk - is this in core API? Move to Utils::Helper?

=cut



# --------------------------------------------------------------------------------
# Project a feature from one slice to another
sub project_feature {
  my ($self, $feat, $new_assembly) = @_;

  # project feature to new assembly
  my $feat_slice = $feat->feature_Slice;


  if(! $feat_slice){
	throw('Cannot get Feature Slice for '.$feat->start.':'.$feat->end.':'.$feat->strand.' on seq_region '.$feat->slice->name);
  }

  my @segments = @{ $feat_slice->project('chromosome', $new_assembly) };

  if(! @segments){
	$self->log("Failed to project feature:\t".$feat->display_label);
	return;
  }
  elsif(scalar(@segments) >1){
	$self->log("Failed to project feature to distinct location:\t".$feat->display_label);
	return;
  }

  my $proj_slice = $segments[0]->to_Slice;
  
  if($feat_slice->length != $proj_slice->length){
	$self->log("Failed to project feature to comparable length region:\t".$feat->display_label);
	return;
  }


  # everything looks fine, so adjust the coords of the feature
  $feat->start($proj_slice->start);
  $feat->end($proj_slice->end);
  $feat->strand($proj_slice->strand);
  my $slice_new_asm = $self->slice_adaptor->fetch_by_region('chromosome', $proj_slice->seq_region_name, undef, undef, undef, $new_assembly);
  $feat->slice($slice_new_asm);

  return $feat;

}

sub slice_adaptor{
  my $self = shift;

  if(! defined $self->{'slice_adaptor'}){
	$self->{'slice_adaptor'} = $self->db->get_SliceAdaptor;
  }
  
  return $self->{'slice_adaptor'};
}


1;
