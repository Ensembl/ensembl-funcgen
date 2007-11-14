package Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Analysis;

# Base functionality for external_feature parsers

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = {};

  #my $self = {(
#			  feature_sets => {(
#								vista   => {('VISTA enhancer set' => undef)},
#								cisred  => {(
#											 'cisRED search regions' => undef,
#											 'cisRED group motifs'   => undef,
#											)},
#								miranda => {('miRanda miRNA' => undef)},
#							   )},
#			  )};
  bless $self, $class;

  #validate and set type, analysis and feature_set here
  my ($type, $db, $clobber, $archive) = rearrange(['TYPE', 'DB', 'CLOBBER', 'ARCHIVE'], @_);
  
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
  
  print ":: Parsing and loading $type ExternalFeatures\n";

  return $self;

}

sub db{
  my ($self, $db) = @_;

  if($db){

	if(! $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor')){
	  throw('You must prove a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
	}
	
	$self->{'db'} = $db;
  }

  return $self->{'db'};
}


sub set_feature_sets{
  my $self = shift;

  throw('Must provide a set feature_set config hash') if ! defined $self->{'feature_sets'};


  my $fset_adaptor = $self->db->get_FeatureSetAdaptor;
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
  
  foreach my $fset_name(keys %{$self->{feature_sets}}){
	
	#check feature_type is validated?

	my $fset = $fset_adaptor->fetch_by_name($fset_name);

	#what about data sets?
	#we don't need data sets for external_sets!

	if(defined $fset){

	  if($self->{'clobber'}){
		#Need to clobber any DBEntries first!!!

		if($self->{'type'} eq 'cisred'){
		  my @ext_feat_ids = 	map @{$_}, @{$self->db->dbc->db_handle->selectall_arrayref('select external_feature_id from external_feature where feature_set_id='.$fset->dbID)};
		  
		  if(@ext_feat_ids){

			my ($core_ext_dbid) = $self->db->dbc->db_handle->selectrow_array('select external_db_id from external_db where db_name="core"');

			if($core_ext_dbid){
			  #double table delete?
			  my $sql = "delete x, ox from object_xref ox, xref x where ox.ensembl_object_type='ExternalFeature' and x.external_db_id=$core_ext_dbid and ox.xref_id=x.xref_id and ensembl_id in(".join(', ', @ext_feat_ids).')';
			  print ":: Clobbering xrefs for $fset_name\n";

			  #warn "sql is $sql";

			  $self->db->dbc->do($sql);
			}
		  }
	
		  print ":: Clobbering old features for external feature_set:\t$fset_name\n";
		  my $sql = 'delete from external_feature where feature_set_id='.$fset->dbID;
		  $self->db->dbc->do($sql);
		}

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

	if(!defined $fset){
	  #don't need to use RNAFeatureType here as this is the setwide generic feature_type
	  #or do we have separate tables for external_feature and external_rna_feature?

	  #validate analysis first
	  my $analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'});
	  if(! defined $analysis){
		
		print ':: Analysis '.$self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'}.
		  " not found, storing from config hash\n";
		

		$analysis_adaptor->store(Bio::EnsEMBL::Analysis->new
								 (
								  %{$self->{'feature_sets'}{$fset_name}{'analysis'}}
								  #-logic_name    => $self->{'feature_sets'}{$fset_name}{'analysis'}{'-logic_name'},
								  #-description   => $self->{'feature_sets'}{$fset_name}{'analysis'}{'-description'},
								  #-display_label => $self->{'feature_sets'}{$fset_name}{'analysis'}{'-display_label'},
								  #-diplayable    => $self->{'feature_sets'}{$fset_name}{'analysis'}{'-displayable'},
								 )
								);
		
		$analysis = $analysis_adaptor->fetch_by_logic_name($self->{'feature_sets'}{$fset_name}{'analysis'});
	  }

	  #replace hash config with object
	  $self->{'feature_sets'}{$fset_name}{'analysis'} = $analysis;

	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													 -name         => $fset_name,
													 -type         => 'external',
													 -analysis     => $self->{'feature_sets'}{$fset_name}{'analysis'},
													 -feature_type => ${$self->{'feature_sets'}{$fset_name}{'feature_type'}},
													);

	  ($fset) = @{$self->db->get_FeatureSetAdaptor->store($fset)};
	}

	#Now replace config hash with object
	$self->{feature_sets}{$fset_name} = $fset;
  }

  return;
}




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

  foreach my $ftype_name(keys %{$self->{'feature_types'}}){

	my $ftype = $ftype_adaptor->fetch_by_name($ftype_name);


	if(! defined $ftype){
	  print ":: FeatureType '".$ftype_name."' for external feature_set ".$self->{'type'}." not present\n".
		":: Storing using type hash definitions\n";
	
	  $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
													   -name => $ftype_name,
													   -class => $self->{'feature_types'}{$ftype_name}{'class'},
													   -description => $self->{'feature_types'}{$ftype_name}{'description'},
													  );
	  ($ftype) = @{$ftype_adaptor->store($ftype)};
	}

	#Replace hash config with object
	$self->{'feature_types'}{$ftype_name} = $ftype;
  }

  return;
}



# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
# Type is always all lower case


sub get_display_name_by_stable_id{
  my ($self, $stable_id, $type) = @_;



  if($type !~ /(gene|transcript|translation)/){
	warn "Cannot get display_name for stable_id $stable_id with type $type\n";
	return;
  }
  
  if(! exists $self->{'display_name_cache'}->{$stable_id}){
	($self->{'display_name_cache'}->{$stable_id}) = $self->db->dnadb->dbc->db_handle->selectrow_array("SELECT x.display_label FROM ${type}_stable_id s, $type t, xref x where t.display_xref_id=x.xref_id and s.${type}_id=t.gene_id and s.stable_id='${stable_id}'");
  }

  return $self->{'display_name_cache'}->{$stable_id};
}

# --------------------------------------------------------------------------------
# Project a feature from one slice to another
sub project_feature {
  my ($self, $feat, $new_assembly) = @_;

  # just use a SimpleFeature for convenience
  #my $feat = Bio::EnsEMBL::SimpleFeature->new
  #  (-start    => $start,
  #   -end      => $end,
  #   -strand   => $strand,
  #   -slice    => $slice,
  #   -analysis => $analysis,
  #   -display_label => $label,
  #   -score   => 0);

  # project feature to new assembly
  my $feat_slice = $feat->feature_Slice;
  my @segments = @{ $feat_slice->project('chromosome', $new_assembly) };


  #next what?
  #next unless (@segments);
  #next if (scalar(@segments) > 1);

  if(! @segments){
	print "Failed to project feature:\t".$feat->display_label."\n";
  }
  elsif(scalar(@segments) >1){
	print "Failed to project feature to distinct location:\t".$feat->display_label."\n";
	return;
  }

  my $proj_slice = $segments[0]->to_Slice;

  if($feat_slice->length != $proj_slice->length){
	print "Failed to project feature to comparable length region:\t".$feat->display_label."\n";
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
