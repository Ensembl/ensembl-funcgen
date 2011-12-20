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

Bio::EnsEMBL::Funcgen::Parsers::BaseImporter

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the start of merging the Importer with the  InputSet & BaseExternal Parsers.
This class holds all generic methods used for importing data into the funcgen schema.
Move all generic methods from Importer to here, and move format specific methods to new parsers.
Then remove Importer completely.

=head1 SEE ALSO


=cut

package Bio::EnsEMBL::Funcgen::Parsers::BaseImporter;

use strict;


use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::FeatureType;
#use Bio::EnsEMBL::Funcgen::Utils::Helper;
#use vars qw(@ISA);
#@ISA = ('Bio::EnsEMBL::Funcgen::Utils::Helper');

use base qw(Bio::EnsEMBL::Funcgen::Utils::Helper); #@ISA change to parent with perl 5.10


#new
#edb_release over-ride, to enable loading of old data.



=head2 validate_and_store_config

  Args       : None
  Example    : $self->validate_and_store_config;
  Description: Imports feature types defined by import_sets.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk 

=cut

#Taken from InputSet
#Needs to support 'static' config from external parsers

#Need to get the config has format working with this and set_feature_sets?
#Or just migrate define_and_validate_sets here now?

#Now takes arrayref, need to change in other callers

#All analyses and ftypes now stored here(inc fset only defined), so don't need to do this in the define sets method

#define/set sets method should write to user config first, call validate_store, before defining sets?

#Should we allow empty hashes to fetch from DB?

sub validate_and_store_config{
  my ($self, $fset_names) = @_;

  if( (ref($fset_names) ne 'ARRAY') ||
	  (scalar(@$fset_names) == 0) ){
	throw('Must pass FeatureSet names to validate_and_store_config');
  }

  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;


  #Do we want to enable no config? Just don't call this if we have no config!

  my ($static_config, $user_config, $fset_config);

  #Set here to avoid auto-vivying in tests below
  $user_config   = $self->{user_config} if exists ${$self}{user_config};
  $static_config = $self->{static_config} if exists ${$self}{static_config};
  my $config = $user_config || $static_config;

  if(! $config){
	throw('No user or static config found');
  }
  elsif($user_config && $static_config){
	throw('BaseImporter does not yet support overriding static config with user config');
	#Would need to over-ride on a key by key basis, account for extra config in either static or user config?
  }
  
  #Store config for each feature set
  #inc associated feature_types and analyses
  #add cell_types in here?

  foreach my $import_set(@$fset_names){
	#warn "validating config for $import_set";

	
	if(exists ${$config}{feature_sets}{$import_set}){
	  $fset_config = $config->{feature_sets}{$import_set};
	}
	
	if(! $fset_config){
	  throw("Could not find config for:\t$import_set");
	}

	$self->log("Validating and storing config for:\t$import_set");

	  
	#If we grab the ftype and analysis from the feature set first
	#Then we can remove the reundancy in the config
	#would need to handle case i.e. -ANALYSIS -analysis
	#use rearrange!

	#Set in analyses and feature_types config, should use same ref if already exist
	#else we have a duplicated entry using the same name, which may have different attrs!
	
	#This is going to be a problem as we are assigning a new value to the key?
	#Do we need to use references in feature_sets?

	if(exists ${$fset_config}{feature_set}){
	  #my ($fset_analysis, $fset_ftype) = rearrange(['ANALYSIS', 'FEATURE_TYPE'], %{$fset_config->{feature_set}});

	  #Can't use rearrange for key we are setting as we need to now the case
	  #This is grabbing top level config keys, not hashes
	  #config names refer to top level keys and may not match logic_name or ftype name
	  my $fset_params          = $fset_config->{feature_set};
	  my $fset_analysis_key    = (exists ${$fset_params}{-analysis})            ? '-analysis'                        : '-ANALYSIS';
	  my $analysis_config_name = (exists  ${$fset_params}{$fset_analysis_key} ) ? $fset_params->{$fset_analysis_key} : undef;
	  my $fset_ftype_key       = (exists ${$fset_params}{-feature_type})        ? '-feature_type'                    : '-FEATURE_TYPE';
	  my $ftype_config_name    = (exists ${$fset_params}{$fset_ftype_key})      ? $fset_params->{$fset_ftype_key}    : undef;
	  	  

	  #Check we have these in feature_sets analyses/feature_types already?
	  #Should that be a hash with empty {} values
	  #Then we can use the top level analyses hash if required
	  #Setting feature_set analysis/feature_type only if not present i.e. we don't over-write the top level
	  #defining hash
	  #Postponing checking here will lose context if we need to throw
	  #i.e. we won't know where the error has come from if the name is not present in the top level hash
	  #Can we call direct from here instead?

  

	  #Leave to store to catch these mandatory attrs

	  if($analysis_config_name){
		
		if(ref(\$analysis_config_name) ne 'SCALAR'){
		  throw("You must set feature_set($import_set) -analysis config as a string value i.e. a key referencing the top level analyses config");
		}
		
		if(! exists ${$fset_config}{analyses}{$analysis_config_name}){ #Set in analyses to validate/store below
		  $fset_config->{analyses}{$analysis_config_name}     = {}; #Don't have to set a ref to top level here, as these will all get set to the same obj ref in validate/store
		}

		#Only use refs in the feature_set
		#top level and feature_set analyses get set to same obj at same time		  
		$fset_config->{feature_set}{$fset_analysis_key} = \$fset_config->{analyses}{$analysis_config_name};
		  
	  }

	  
	  if($ftype_config_name){		
	
		if(ref(\$ftype_config_name) ne 'SCALAR'){
		  throw("You must set feature_set($import_set) -feature_type config as a string value i.e. a key referencing the top level feature_types config");
		}

		if(! exists ${$fset_config}{feature_types}{$ftype_config_name}){ #Set in analyses to validate/store below
		  $fset_config->{feature_types}{$ftype_config_name} = {};
		}
		
		$fset_config->{feature_set}{$fset_ftype_key}  = \$fset_config->{feature_types}{$ftype_config_name};

	  }
	}

	
	#Can self ref user config if 'do' will work with %config, specified as the last line
	#Merge these two loops?
	#Need to account for additional config keys in user/static config

	
	if(exists ${$fset_config}{'analyses'}){
	  
	  foreach my $logic_name(keys %{$fset_config->{'analyses'}}){
		$fset_config->{'analyses'}{$logic_name} =
		  $self->validate_and_store_analysis($logic_name, $fset_config->{'analyses'}{$logic_name});
	  }
	}

	if(exists ${$fset_config}{'feature_types'}){
	  
	  foreach my $ftype_key(keys %{$fset_config->{'feature_types'}}){
		$fset_config->{'feature_types'}{$ftype_key} = 
		  $self->validate_and_store_feature_type($ftype_key, $fset_config->{'feature_types'}{$ftype_key});
	  }
	}
	
	#if(exists ${$fset_config}{'feature_set'}){
	#  #Just ignore this for now, as set_feature_set and define_and_validate_sets deal with this repectively
	#  #Solution is to extend define_and_validate_set to support external feature static config
	#  #Then remove set_feature_sets
	# Also need to consider InputSet::define_sets
	#}
  }

  return;
}


#Need some method to get and add analyses and ftypes to reduce hardcoding of key strings
#Can we add Class->new to config to provide some of this?
#Would proliferate constructor calls into config
#May also prevent some logic imposed by validate/set methods i.e. defaults?
#Has to be one or other as we can't assume we can use methods if we are dealing with a hash

#Should we allow empty hashes to default to DB entry?
#Move HASH test into these methods

#Need to use rearrange in these methods for case safety
#rearrange is quite slow, fine for import

#Currently can't use empty hash to default to DB, as this gives us no key in the config
#Can we set this to the string required for the fetch method in the referenced hash(analyses/feature_types)?
# eq test would work, but would have to skip HASH test

#Maybe we just set the analysis in the feature_set by the logic_name string, rather than a ref
#Then we test for presence in the relevant hash
#Would have to set obj and ref when validating/storing
#Would also have to do this in the feature_set feature_types/analyses i.e. we would be able to ref the whole hash
#would have to be an arrayrefs

#Do we want to support no entries in top level definition if we can find it in the DB?


#Could merge these two using a config hash to define adaptor method, params etc?
#Would be easy to add more types to store e.g. cell_type

sub validate_and_store_analysis{
  my ($self, $logic_name, $analysis_params) = @_; 
  my $analysis;
  eval {$self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $analysis_params)};

  if(! $@){ #Can assume we have already set this valid Analysis
	$analysis = $analysis_params
  }
  else{
	my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
	
	### Validate config entries
	#Catches undef or empty {} at feature_sets and top levels, defaults to DB
	my $invalid_entry = 1;
	my $obj_logic_name;
	my $got_config = 0;


  if( (! defined $analysis_params ) ||
	  (ref($analysis_params) eq 'HASH') ){ #We have a valid entry
	
	if(  (ref($analysis_params) eq 'HASH') && 
		 (%{$analysis_params}) ){ #number of keys
	  $got_config = 1;
	  $invalid_entry = 0;
	}
	else{ #empty feature_set analyses hash -> check top level analyses first
	  #$invalid_entry = 1;

	  if(exists ${$self->{static_config}{analyses}}{$logic_name}){
		$analysis_params = $self->{static_config}{analyses}{$logic_name};
		
		if( (! defined $analysis_params ) ||
			(ref($analysis_params) eq 'HASH') ){
		  $invalid_entry = 0;
		  
		  if( (ref($analysis_params) eq 'HASH') &&
			  (%{$analysis_params}) ){ #number of keys
			$got_config = 1;
		  }#else is empty top level {}

		}#else is invalid
	  }
	  else{ #No top level config, assume we want to use the DB
		$invalid_entry = 0;
	  }
	}
  }
  
  
  if($invalid_entry){
	throw("You have defined a none HASH value in your config for analysis:\t$logic_name\n".
		  "Please define config as HASH, or use empty HASH or undef to use existing default config or DB");
  }
	
  if($got_config){
	($obj_logic_name) = rearrange(['LOGIC_NAME'], %{$analysis_params});
  }else{
	$obj_logic_name   = $logic_name;
  }
  	
  
  if($logic_name ne $obj_logic_name){ #Not a show stopper as this is just the config key
	warn "Found analysis key name - logic_name mismatch in config:\t$logic_name vs $obj_logic_name\n";
  }
 
  $analysis         = $analysis_adaptor->fetch_by_logic_name($obj_logic_name);
  
  

  if($got_config){

	my $config_anal      = Bio::EnsEMBL::Analysis->new(%{$analysis_params});
	
	if(! defined $analysis){	
	  $self->log('Analysis '.$obj_logic_name." not found in DB, storing from config");		
	  $analysis_adaptor->store($config_anal);
	  $analysis = $analysis_adaptor->fetch_by_logic_name($obj_logic_name);	
	}
	else{
	  
	  my $not_same = $analysis->compare($config_anal);
	  #Analysis::compare returns the opposite of what you expect!
	  
	  if($not_same){
		throw('There is a param mismatch between the '.$obj_logic_name.
			  ' Analysis in the DB and config. Please rectify or define a new logic_name');
	  }
	}
  }
  elsif(! defined $analysis){
	throw("Cannot fetch $obj_logic_name analysis from DB, please check your config key or define new top level analyses config");
  }
}

  return $analysis;
}

sub validate_and_store_feature_type{
  my ($self, $ftype_name, $ftype_params) = @_;
  my $ftype;
  eval {$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype_params)};
  
  if(! $@){ #Can assume we have already set this valid Analysis
	$ftype = $ftype_params
  }
  else{
	my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
	
	#Validate config entries
	#Catches undef or empty {} at feature_sets and top levels, defaults to DB
	my $invalid_entry = 1;
	my ($name, $class);
  
  my $got_config = 0;

  if( (! defined $ftype_params ) ||
	  (ref($ftype_params) eq 'HASH') ){ #We have a valid entry
	
	if(  (ref($ftype_params) eq 'HASH') && 
		 (%{$ftype_params}) ){ #number of keys
	  $got_config = 1;
	  $invalid_entry = 0;
	}
	else{ #empty feature_set ftype hash -> check top level ftype first
	  #$invalid_entry = 1;
	  
	  if(exists ${$self->{static_config}{feature_types}}{$ftype_name}){
		$ftype_params = $self->{static_config}{feature_types}{$ftype_name};
		
		if( (! defined $ftype_params ) ||
			(ref($ftype_params) eq 'HASH') ){
		  $invalid_entry = 0;
		  
		  if( (ref($ftype_params) eq 'HASH') &&
			  (%{$ftype_params}) ){ #number of keys
			$got_config = 1;
		  }#else is empty top level {}

		}#else is invalid
	  }
	  else{ #No top level config, assume we want to use the DB
		$invalid_entry = 0;
	  }
	}
  }
  
 
  if($invalid_entry){
	throw("You have defined a none HASH value in your config for feature_type:\t$ftype_name\n".
		  "Please define config as HASH, or use empty HASH or undef to use existing default config or DB");
  }
	
  if($got_config){
	($name, $class) = rearrange(['NAME', 'CLASS'], %{$ftype_params});
  }else{
	$name   = $ftype_name;
  }
  	

   
  #Can't use rearrange for key we are setting as we need to now the case
  my $analysis_key = (exists ${$ftype_params}{-analysis}) ? '-analysis' : '-ANALYSIS';
  my $analysis;

  if(exists ${$ftype_params}{$analysis_key}){
	#This is slightly redundant as we may have already validated this analysis
	my ($lname) = rearrange(['LOGIC_NAME'], %{$ftype_params->{$analysis_key}});
	$ftype_params->{$analysis_key} = $self->validate_and_store_analysis($lname, $ftype_params->{$analysis_key});
	$analysis = $ftype_params->{$analysis_key};
  }
  
  my @ftypes = $ftype_adaptor->fetch_by_name($name, $class, $analysis);
   
  if(scalar(@ftypes) > 1){
	throw("Unable to fetch unique feature_type $name. Please specify top level config to define class (and analysis)");
  }
  else{
	$ftype = $ftypes[0];#can be undef
  }
  
  
  if($got_config){
	my $config_ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(%{$ftype_params});

	if($ftype){
	  
	  if(! $ftype->compare($config_ftype)){
		my $label = $name."($class";
		$label .= (defined $analysis) ? ' '.$analysis->logic_name.')' : ')';
		
		throw('There is a param mismatch between the '.$name.
			  ' FeatureType in the DB and config. Please rectify in the config.');
	  }
	}
	else{
	  $self->log('FeatureType '.$name." not found in DB, storing from config");		
	  ($ftype) = @{$ftype_adaptor->store($config_ftype)};
	}
  }
  elsif(! defined $ftype){
	throw("Cannot fetch $name feature_type from DB, please check your config key or define new top level feature_types config");
  }
}

  return $ftype;
}


1;
