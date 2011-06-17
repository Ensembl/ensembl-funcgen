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

Bio::EnsEMBL::Funcgen::Utils::HealthChecker
  
=head1 SYNOPSIS

=head1 DESCRIPTION

B<This program> performs several health check and update methods prior to release.

=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::HealthChecker;

use strict;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Analysis;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);



#TO DO
# 1 DONE Print all fails and warnings in summary at end of script.
# 2 validate_RegulatoryFeature_Sets
# 3 Some of these can be migrated or mirrored in java HCs for safety


################################################################################

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
    
  #validate and set type, analysis and feature_set here
  my ($db, $builds, $skip_mc, $check_displayable, $skip_analyse, $meta_coord_tables, $skip_xrefs) = 
	rearrange(['DB', 'BUILDS', 'SKIP_META_COORD', 'CHECK_DISPLAYABLE', 'SKIP_ANALYSE', 'META_COORD_TABLES', 'SKIP_XREF_CLEANUP'], @_);
  
  
  if (! ($db && ref($db) &&
		 $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('You must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }
  
  #test connection
  $db->dbc->db_handle;
  
  $self->{'db'} = $db;
  $self->{'mysql_connect_string'} = 'mysql -h'.$db->dbc->host.' -u'.$db->dbc->username.' -p'
	.$db->dbc->password.' '.$db->dbc->dbname.' -P'.$db->dbc->port;
  $self->{'dbname'} = $db->dbc->dbname;
  $self->{'builds'} = (scalar(@$builds)>0) ? $builds : [];
  $self->{'skip_meta_coord'} = $skip_mc;
  $self->{'skip_xrefs'} = $skip_xrefs;
  $self->{'skip_analyse'} = $skip_analyse;
  $self->{'check_displayable'} = $check_displayable;
  
  if(defined $meta_coord_tables){

	throw('-skip_meta_coord is set, Cannot build meta_coord entries for tables '.join(', ', @$meta_coord_tables));

	if(! ref($meta_coord_tables) eq 'ARRAY'){
	  throw('-meta_coord_tables parameter must be an array ref');
	}

	@{$self->{'meta_coord_tables'}} = @$meta_coord_tables;
  }

  return $self;
}

sub db{
  my ($self) = @_;
  
  return $self->{'db'};
}



=head2 update_db_for_release

  Arg[0]     : 
  Example    : 
  Description: Wrapper method to perform all common update functions
  Returntype : 
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub update_db_for_release{
  my ($self, @args) = @_;
  
  if(@args){
  }
	
  #do seq_region_update to validate dnadb first
  #hence avoiding redoing longer methods
  $self->validate_new_seq_regions;#($force_srs);
  $self->update_meta_schema_version;
  $self->check_meta_strings;
  $self->check_meta_species;
  $self->set_current_coord_system;
  $self->update_meta_coord;
  $self->clean_xrefs;
  $self->validate_DataSets;
  $self->check_stable_ids;
    $self->log_data_sets();
  $self->analyse_and_optimise_tables;#ALWAYS LAST!!

  $self->log_header('??? Have you dumped/copied GFF dumps ???');

  #Log footer? Pass optional counts hash?
  $self->log('Finished updating '.$self->{'dbname'}." for release\n\n");
}

sub validate_new_seq_regions{
  my ($self, $force) = @_;
  

  #We need to add some functionality to handle non-standard schema_build progression here
  #This should be used before any data is loaded
  #It should also warn if there any duplicates if it is not run before coord_systems are duplicated
  #Should this just be handled in the BaseFeatureAdaptor/CoordSystemAdaptor?
  

  
  
  #do we need to add the none default levels here?
  #or are we only bothered about those which constitute the toplevel?

  #To make sure we have all the correct levels in eFG we need to get all the names.
  #then get all by name from the core db and set them as the dnadb.
  # we also need to get all the toplevel seq_regions and store them in the seq_region table
  #use BaseFeatureAdaptor::_pre_store with and array of pseudo feature on each top level slice

  #Validate the efgdb and dnadb schema version are the same first
  #This is because we expect the schem_build to be the same for a release
  #However, this will autoset the dnadb if no defined, so will always match!

  if(! $force){
	my $efgdb_sm = join('_', @{$self->get_schema_and_build($self->{'dbname'})});
	my $dnadb_sm = join('_', @{$self->get_schema_and_build($self->db->dnadb->dbc->dbname)});
	
	if($efgdb_sm ne $dnadb_sm){
	  $self->report("WARNING Skipped validate_new_seq_regions as schema_versions are mismatched:\t".
				 "efgdb $efgdb_sm\tdnadb $dnadb_sm");
	  return 0;
	}
  }
  
  my $pf_adaptor    = $self->db->get_ProbeFeatureAdaptor();
  my $slice_adaptor = $self->db->dnadb->get_SliceAdaptor();
  my $dnadb_csa     = $self->db->dnadb->get_CoordSystemAdaptor;
  
  $self->log_header('Validating new coord_systems/seq_regions');
  
  my @slices;
  my %versioned_levels;
  my $default_version;
  
  #Grab unversioned top level slices and versioned levels
  #May miss some old versioned level if the new assembly no longer has them
  foreach my $slice(@{$slice_adaptor->fetch_all('toplevel', undef, 1)}){

	if (! $slice->coord_system->version){
	  push @slices, $slice;
	}
	else{

	  if($default_version &&
		  ($default_version ne $slice->coord_system->version)){
		throw("Found more than one default CoordSystem version:\t${default_version}\t".$slice->coord_system->version);
	  }
	  else{
		$default_version = $slice->coord_system->version;
	  }
	}
  }
  
 
  #Get all versioned levels for all builds
  foreach my $cs(@{$dnadb_csa->fetch_all}){
	
	if($cs->version){
	  $versioned_levels{$cs->version} ||= [];	
	  push @{$versioned_levels{$cs->version}}, $cs->name;
	}
  }
  
  push @{$self->{'builds'}}, $default_version if scalar(@{$self->{'builds'}}) == 0;
  


  #Grab slices for each versioned level
  foreach my $build(@{$self->{'builds'}}){
	
	if(! exists $versioned_levels{$build}){
	  throw("CoordSystem version $build does not exist in the dnadb ".$self->db->dnadb->dbc->dbname);
	}
	
	foreach my $level(@{$versioned_levels{$build}}){
	  $self->log("Getting slices for $level $build");
	  push @slices, @{$slice_adaptor->fetch_all($level, $build)};
	}
  }
  
  $self->log("Importing seq_region/coord_system info for builds:\t".join(',', @{$self->{'builds'}}));
  
  foreach my $slice(@slices){

	if($slice->start() != 1){
	  $self->log("Reslicing slice:\t".$slice->name());
	  #we must have some sort of PAR linked region i.e. Y
	  $slice = $slice_adaptor->fetch_by_region($slice->coord_system_name(), $slice->seq_region_name());
	}
	
	
	#we need test if it needs doing first?
	#we would need to test for the coord_systems outside of this loop
	#and then for each seq_region inside the loop if the coord_system is present
	
	$self->log("_pre_storing seq_region info for slice:\t".$slice->name());
	
	my $pseudo_feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
	  (
	   -slice => $slice,
	   -start => 0,
	   -end   => 0,
	   -strand => 0,
		);
	
	$pf_adaptor->_pre_store($pseudo_feature);
	#This will create a meta_coord entry of max_length 1 for features which have an absent meta_coord entry
	
  }

  
  $self->log("Finished validating seq_regions\n");
  
  return;
}



sub update_meta_schema_version{
  my ($self) = @_;
  
  my $schema_version = $self->get_schema_and_build($self->{'dbname'})->[0];
  
  

  my $sql = 'DELETE from meta where meta_key="schema_version"';
  $self->db->dbc->db_handle->do($sql);
  $sql = "INSERT into meta(meta_key, species_id, meta_value) values ('schema_version', NULL, '$schema_version')";
 
  $self->db->dbc->db_handle->do($sql);

  $self->log_header("Updated meta.schema_version to $schema_version");

}




sub update_meta_coord{
  my ($self, @table_names) = @_;
  
  my $species_id = $self->db()->species_id();
  
  if($self->{'skip_meta_coord'}){
	$self->log("Skipping meta_coord update\n");
	return;
  }
  

  $self->log_header('Updating meta_coord table');


  #set default table_name
  if(! @table_names || scalar(@table_names) == 0){

	#Can we do this via DBAdaptor and get all available adaptors which are BaseFeatureAdaptors then grab the first table name

	if(defined $self->{'meta_coord_tables'}){
	  @table_names = @{$self->{'meta_coord_tables'}};
	}
	else{#default

	  @table_names = qw(
						regulatory_feature
						probe_feature
						external_feature
						annotated_feature
						result_feature
					 );
	}
  }

  #backup meta coord
  #if(system("mysql -h$host -P$port -u$user -p$pass -N "
#			. "-e 'SELECT * FROM meta_coord' ${species}_funcgen_${schema_build}"
#			. "> ${species}_funcgen_${schema_build}.meta_coord.backup"
#		   ) != 0 ){

  if(system($self->{'mysql_connect_string'}." -e 'SELECT * FROM meta_coord'"
			. '> '.$self->{'dbname'}.'meta_coord.backup'
		   ) != 0 ){

	throw("Can't dump the original meta_coord for back up");#will this get copied to log?
  } 
  else {
	$self->log('Original meta_coord table backed up in '. $self->{'dbname'}.'.meta_coord.backup');
  }
  

  #Update each max_length for table_name and coord_system

  my $sql;

  foreach my $table_name(@table_names){
  
	$sql = "select distinct(cs.name), mc.coord_system_id, cs.version, mc.max_length from coord_system cs, meta_coord mc where mc.table_name='$table_name' and mc.coord_system_id=cs.coord_system_id and cs.species_id = $species_id";
	
	$self->log('');
	$self->log("Updating meta_coord max_length for $table_name:");
	$self->log("name\tcoord_system_id\tversion\tmax_length");
	
	#can we test for emtpy array here? Then skip delete.
	
	my @info = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
	
	#log this
	map {$self->log(join("\t", @{$_}))} @info;
	
	# Clean old entries
	$self->log("Deleting old meta_coord entries");
	$sql = "DELETE mc FROM meta_coord mc, coord_system cs WHERE mc.table_name ='$table_name' and mc.coord_system_id = cs.coord_system_id and cs.species_id = $species_id";
	$self->db->dbc->db_handle->do($sql);
	
	# Generate new max_lengths
	$self->log("Generating new max_lengths");
	
	#Is this query running for each redundant cs_id?
	#would it be more efficient to retrieve the NR cs_ids first and loop the query for each cs_id?
	
	#Can we get the dbID of the largest feature for ease of checking?
	#This won't work as we're grouping by coord_system
	#would need to select distinct coord_system_id for table first
	#This may well slow down quite a bit doing it this way
	
	$sql = "select distinct s.coord_system_id from coord_system cs, seq_region s, $table_name t WHERE t.seq_region_id = s.seq_region_id and s.coord_system_id = cs.coord_system_id and cs.species_id = $species_id";
	my @cs_ids = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
	#Convert single element arrayrefs to scalars
	map $_ = ${$_}[0], @cs_ids;

	$self->log("New max_lengths for $table_name are:");
	
	$self->log(join("\t", ('coord_system_id', 'max_length', 'longest record dbID')));

	foreach my $cs_id(@cs_ids){
	  #This will always give a length of 1 even if there are no features present

	  my @cs_lengths;

	  #The probe_feature table is now too big to do this in one go
	  #We need to break this down into sr_ids

	  $sql = "SELECT distinct t.seq_region_id from $table_name t, seq_region sr where t.seq_region_id=sr.seq_region_id and sr.coord_system_id=$cs_id";
	  
	  my @sr_ids = @{$self->db->dbc->db_handle->selectcol_arrayref($sql)};
	  

	  #Get longest feature for all seq_regions
	  foreach my $sr_id(@sr_ids){
		$sql = "SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.${table_name}_id "
		  . "FROM $table_name t "
			. "WHERE t.seq_region_id = $sr_id ";
		$sql .= ' and t.window_size=0' if $table_name eq 'result_feature';
		$sql .= " order by max desc limit 1";
		

		#Problem here is that DBs without 0 wsize result_feture entries will not get a meta_coord entry
		#We need to implement this in the _pre_store method too?


		my ($cs_length, $table_id);
		($cs_length, $table_id) = $self->db->dbc->db_handle->selectrow_array($sql);
		push @cs_lengths, [$cs_length, $table_id] if $cs_length;
	  }


	  if(@cs_lengths){
		#This will now contain a list of arrays refs contain the max length and feature id for 
		#each seq_region in this coord_system
		#Now sort to get the longest
		#Can't sort on 2 day array in the normal way
		#One list list lists, the comparatee is no longer a list but a reference
		@cs_lengths = sort { $b->[0] <=> $a->[0] } @cs_lengths;
		$self->log(join("\t\t", ($cs_id, @{$cs_lengths[0]})));
		$sql = "INSERT INTO meta_coord values('${table_name}', $cs_id, ".$cs_lengths[0][0].')';
		$self->db->dbc->db_handle->do($sql);
	  }
	}
  }  
  
  $self->log("Finished updating meta_coord max_lengths\n");
  
  return;
}




sub check_meta_species{
  my ($self) = @_;

  $self->log_header('Checking meta species.production_name');

  my $mc = $self->db->get_MetaContainer;
  my @latin_names = @{$mc->list_value_by_key('species.production_name')};
  
  my $dbname_species = $self->db->dbc->dbname;
  $dbname_species =~ s/_funcgen_.*//;
  
  if(scalar(@latin_names) > 1){
	$self->report("FAIL:\tFound more than one species.production_name in meta:\t".join(", ", @latin_names));
  }
  elsif(scalar(@latin_names) == 1 && ($latin_names[0] ne $dbname_species)){
	$self->report("FAIL:\tFound mismatch between meta species.production_name and dbname:\t".$latin_names[0]." vs $dbname_species");
  }
  elsif(scalar(@latin_names) == 0){
	$self->report("WARNING:\tFound no meta species.production_name setting as:\t$dbname_species");
	$self->db->dbc->db_handle->do("INSERT into meta(species_id, meta_key, meta_value) values(1, 'species.production_name', '$dbname_species')");
  }

  #else is okay

  return;
}



#Move to Java HC? Or update if update flag specified
#Using same code used by build_reg_feats!

sub check_meta_strings{
  my ($self, $update) = @_;
  

  $self->log_header('Checking meta strings');

  warn "Need to check/update regbuild.version and regbuild.initial_release_date regbuild.last_annotation_update";

  #update flag?

  my @regf_fsets;
  my $passed = 1;
  my $fset_a = $self->db->get_FeatureSetAdaptor;
  my $mc = $self->db->get_MetaContainer;
  my $regf_a = $self->db->get_RegulatoryFeatureAdaptor;
  #We now want to chek all build
  @regf_fsets = @{$fset_a->fetch_all_by_type('regulatory')};
 
  
  if(scalar(@regf_fsets) == 0){
	$self->report("WARNING: Found no regulatory FeatureSets for check_meta_strings");
  }
  else{
	
	warn "Need to check/update regbuild.version and regbuild.initial_release_date regbuild.last_annotation_update";
  

	#How do we validate this?
	#Check all feature_sets exist
	#Pull back some features from a test slice and check the number of bits match.
	#Check the feature_type string exists and matches else create.
	
	
	foreach my $fset(@regf_fsets){
	  $self->log_header("Validating meta entries for FeatureSets:\t".$fset->name);
		
	  #Fail for old versions as we want to remove these
	  if( $fset->name =~ /_v[0-9]+$/){
		$self->report("FAIL:\t".$fset->name." is an old RegulatoryFeature set, please remove!");
	  }

	  my $cell_type = (defined $fset->cell_type) ? $fset->cell_type->name : 'core';
	  my $fset_string_key  = "regbuild.${cell_type}.feature_set_ids";
	  my $ftype_string_key = "regbuild.${cell_type}.feature_type_ids";
	  my $fset_string = $mc->list_value_by_key($fset_string_key)->[0];
	  my $ftype_string = $mc->list_value_by_key($ftype_string_key)->[0];

	  $self->log('Validating '.$fset->name.":\n\t$fset_string_key($fset_string) vs $ftype_string_key($ftype_string)");

	  #Test fset vs ftype string
	  if(! defined $fset_string && ! defined $ftype_string){
		$self->report("FAIL:\tNo $fset_string_key or $ftype_string_key found in meta table");
	  }
	  elsif(! defined $fset_string){
		$self->report("FAIL:\tNo $fset_string_key found in meta table");
	  }
	  else{
		my @fset_ids = split/,/, $fset_string;
		my @ftype_ids;
		my @new_ftype_ids;
		my $ftype_fail = 0;
		
		if(defined $ftype_string){
		  @ftype_ids = split/,/, $ftype_string;
		}
		else{
		  $self->report("WARNING:\tNo $ftype_string_key found in meta table, will update using $fset_string_key");
		}
	  
		
		#Now need to work backwards through ftypes to remove pseudo ftypes before validating
		#New string should be A,A,A;S,S,S,S,S,S;P,P,P
		#Where A is and Anchor/Seed set
		#S is a supporting set
		#P is a pseudo feature type e.g. TSS proximal
		

		if(scalar(@fset_ids) != scalar(@ftype_ids)){
		  $self->report("FAIL:\tLength mismatch between $fset_string_key and $ftype_string_key");
		}
		
		foreach my $i(0..$#fset_ids){
		  my $supporting_set_id = $fset_ids[$i];
		  my $sset = $fset_a->fetch_by_dbID($supporting_set_id);
		  
		  if(! defined $sset){
			$self->report("FAIL:\t$fset_string_key $supporting_set_id does not exist in the DB");
		  }
		  else{
			#test/build ftype string
		  
			if(defined $ftype_string){
			  
			  if($sset->feature_type->dbID != $ftype_ids[$i]){
				$ftype_fail = 1;
				$self->report("FAIL:\t$fset_string_key $supporting_set_id(".$sset->name.") FeatureType(".$sset->feature_type->name.") does not match $ftype_string_key $ftype_ids[$i]");
			  }
			}
	
			push @new_ftype_ids, $sset->feature_type->dbID;
	
		  }
		}


		#Set ftype_string
		#This will not account for pseudo ftypes?  Remove!!!?
		my $new_ftype_string = join(',', @new_ftype_ids);

		if(! defined $ftype_string){
		  $self->log("Updating $ftype_string_key to:\t$new_ftype_string");
		  $self->db->dbc->db_handle->do("INSERT into meta(species_id, meta_key, meta_value) values(1, '$ftype_string_key', '$new_ftype_string')");
		}
		elsif($ftype_fail){
		  $self->report("FAIL:\t$ftype_string_key($ftype_string) does not match $fset_string_key types($new_ftype_string)");
		}


		#Finally validate versus a reg feat
		#Need to change this to ftype string rather than fset string?
		my $id_row_ref = $self->db->dbc->db_handle->selectrow_arrayref('select regulatory_feature_id from regulatory_feature where feature_set_id='.$fset->dbID.' limit 1');

		if(! defined $id_row_ref){
		  $self->report("FAIL:\tNo RegulatoryFeatures found for FeatureSet ".$fset->name);
		}
		else{
		  my ($regf_dbID) = @$id_row_ref;
		  my $rf_string = $regf_a->fetch_by_dbID($regf_dbID)->binary_string;
		  
		  if(length($rf_string) != scalar(@fset_ids)){
			$self->report("FAIL:\tRegulatory string length mismatch between RegulatoryFeature($regf_dbID) and $fset_string_key:\n$rf_string(".length($rf_string).")\n$fset_string(".scalar(@fset_ids).")");
		  }
		}
	  }
	}
  }

  return;
}


#Change this to log sets and incorporate RegFeat FeatureSet as standard
#Grab all reg fsets
#grab all displayable data sets which aren't reg sets?

sub log_data_sets{
  my $self = shift;
  
  my ($status);
  my $txt = 'Checking ';

  $status = 'DISPLAYABLE' if($self->{'check_displayable'});

	
  my @dsets;
  my $dsets = $self->db->get_DataSetAdaptor->fetch_all($status);
  @dsets = @$dsets if defined $dsets;

  $txt .= scalar(@dsets).' ';
  $txt.= $status.' ' if($self->{'check_displayable'});

  $txt .= 'DataSets';

  $self->log_header($txt);
  
  
  foreach my $dset(@dsets){
	$self->log_set("DataSet:\t\t", $dset) ;

	my $fset = $dset->product_FeatureSet;
	$self->log_set("Product FeatureSet:\t", $fset) if $fset;
	
	my @supporting_sets = @{$dset->get_supporting_sets};
	
	$self->log('Found '.scalar(@supporting_sets).' supporting sets:');

	if(my @supporting_sets = @{$dset->get_supporting_sets}){
	  #type here could be result, experimental or feature
	  #and feature could be annotated or experimental
	  #Move this to log set?

	  map { my $type = ref($_);
			$type =~ s/.*://;			
			$type .= '('.$_->feature_class.')' if($type eq 'FeatureSet');
			#Need to sprintf $type here to fixed width
			$self->log_set($type.":\t", $_)} @supporting_sets;
	}
	$self->log();
  }

  return;
}

sub log_set{
  my ($self, $text, $set) = @_;
  
  $text .= $set->display_label.'('.$set->name.')';
  $text .= "\tDISPLAYABLE" if($set->is_displayable);
  $self->log($text);

  return;
}


sub check_stable_ids{
  my ($self, @slices) = @_;
  
  my $species_id = $self->db()->species_id();

  $self->log_header('Checking stable IDs');

  my $fset_a = $self->db->get_FeatureSetAdaptor;

  my @regf_fsets = @{$fset_a->fetch_all_by_type('regulatory')};

  if(!@regf_fsets){
	$self->report('WARNING: No regulatory FeatureSets found');
  }
  else{

	foreach my $fset(@regf_fsets){

	  if($fset->name =~ /_v[0-9]$/){
		$self->log("Skipping stable_id test on archived set:\t".$fset->name);
		next;
	  }
	  
	  #Can't count NULL field, so have to count regulatory_feature_id!!!

	  #getting SR product here!!
	  my $sql = "select count(rf.regulatory_feature_id) from regulatory_feature rf, seq_region sr, coord_system cs where rf.stable_id is NULL and rf.seq_region_id = sr.seq_region_id and sr.coord_system_id = cs.coord_system_id and cs.species_id = $species_id and rf.feature_set_id=".$fset->dbID;
	  
	  
	  
	  my ($null_sids) = @{$self->db->dbc->db_handle->selectrow_arrayref($sql)};
	  
	  if($null_sids){
		$self->report("FAIL: Found a total of $null_sids NULL stable IDs for ".$fset->name);
		
		my $slice_a = $self->db->get_SliceAdaptor;
		
		if(! @slices){
		  @slices = @{$slice_a->fetch_all('toplevel', 1)};
		}
		
		foreach my $slice(@slices){
		  my $sr_name=$slice->seq_region_name;
		  $sql = 'select count(rf.stable_id) from regulatory_feature rf, seq_region sr, coord_system cs where rf.seq_region_id=sr.seq_region_id and sr.name="'.$sr_name.'" and sr.coord_system_id = cs.coord_system_id and cs.species_id = '.$species_id.' and rf.stable_id is NULL and rf.feature_set_id='.$fset->dbID;
		  ($null_sids) = @{$self->db->dbc->db_handle->selectrow_arrayref($sql)};
		  
		  #This is not reporting properly.

		  $self->log($fset->name.":\t$null_sids NULL stable IDs on ".$slice->name) if $null_sids;
		}
	  }
	  else{
		$self->log($fset->name.":\tNo NULL stable IDs found");
	  }
	}
  }

  return;

}


#This is for mart to enable them to join to the seq_region table without
#creating a product from all the reundant seq_region entries for each schema_build

sub set_current_coord_system{
  my ($self) = @_;



  my $schema_build = $self->db->_get_schema_build($self->db->dnadb);
  $self->log_header("Setting current coord_system on $schema_build");

  my $sql = "update coord_system set is_current=False where schema_build !='$schema_build'";
  $self->db->dbc->do($sql);
  $sql = 'update coord_system set is_current=True where schema_build ="'.$schema_build.'" and attrib like "%default_version%"';
  $self->db->dbc->do($sql);

  return;
}

sub validate_DataSets{
  my $self = shift;
  
  $self->log_header('Validating DataSets');


  #checks regualtory feature data and supporting sets
  #links between DataSet and FeatureSet, i.e. correct naming, not linking to old set
  #Displayable, DAS_DISPLAYABLE, IMPORTED_ASSM, MART_DISPLAYABLE
  #Naming of result_set should match data/feature_set
  #warn non-attr feature/result sets which are DISPLAYABL
  #warn about DISPLAYABLE sets which do not have displayable set in analysis_description.web_data


  my $fset_a = $self->db->get_FeatureSetAdaptor;
  my $dset_a = $self->db->get_DataSetAdaptor;
  my ($dset_states, $rset_states, $fset_states) = $self->get_regbuild_set_states($self->db);

  my %rf_fsets;
  my %set_states;
  my $sql;

 RF_FSET: foreach my $rf_fset(@{$fset_a->fetch_all_by_type('regulatory')}){
	my $rf_fset_name = $rf_fset->name;
	

	$self->log("Validating $rf_fset_name");


	$rf_fsets{$rf_fset_name} = $rf_fset;#Do we only need the name for checking the dsets independantly?

	if($rf_fset_name =~ /_v[0-9]+$/){
	  $self->report("FAIL:\tFound archived regulatory FeatureSet:\t$rf_fset_name");
	  next RF_FSET;
	}
	
	foreach my $state(@$fset_states){
	  
	  if(! $rf_fset->has_status($state)){
		$self->report("WARNING:\tUpdating FeatureSet $rf_fset_name with status $state");
		
		$sql = 'INSERT into status select '. $rf_fset->dbID.
		  ", 'feature_set', status_name_id from status_name where name='$state'";

		$self->db->dbc->db_handle->do($sql);
	  }
	}

	#Do we need to warn about other states?
		


	my $rf_dset = $dset_a->fetch_by_product_FeatureSet($rf_fset);

	if(! $rf_dset){
	    $self->report("FAIL:\tNo DataSet for FeatureSet:\t$rf_fset_name");
	  next RF_FSET;
	}

	

	if($rf_fset_name ne $rf_dset->name){
	  $self->report("FAIL:\tFound Feature/DataSet name mismatch:\t$rf_fset_name vs ".$rf_dset->name);
	  next RF_FSET;
	}

	
	foreach my $state(@$dset_states){
	  
	  if(! $rf_dset->has_status($state)){
		#Can we just update and warn here?
		#Or do this separately in case we want some control over this?
		$self->report("WARNING:\tUpdating DataSet $rf_fset_name with status $state");

		$sql = 'INSERT into status select '.$rf_dset->dbID.
		  ", 'data_set', status_name_id from status_name where name='$state'";
		$self->db->dbc->db_handle->do($sql);
	  }
	}



	#We should check fset ctype matches all attr_set ctypes?
	#May have problems if we want to merge two lines into the same ctype build



	foreach my $ra_fset(@{$rf_dset->get_supporting_sets}){
	
	  foreach my $state(@$fset_states){
	  
		if(! $ra_fset->has_status($state)){
		  #Can we just update and warn here?
		  #Or do this separately in case we want some control over this?
		  $self->report("WARNING:\tUpdating FeatureSet ".$ra_fset->name." with status $state");

		  $sql = 'INSERT into status select '.$ra_fset->dbID.
			", 'feature_set', status_name_id from status_name where name='$state'";
		  $self->db->dbc->db_handle->do($sql);

		}
	  }



	  my $ra_dset = $dset_a->fetch_by_product_FeatureSet($ra_fset);


	  my @sset = @{$ra_dset->get_displayable_supporting_sets('result')};
	  
	  if(scalar(@sset) > 1){#There should only be one
		$self->report("FAIL:\tThere should only be one DISPLAYABLE supporting ResultSet for DataSet:\t".$ra_dset->name);
	  }
	  elsif(scalar(@sset) == 0){
		$self->report("FAIL:\tThere are no DISPLAYABLE supporting ResultSet for DataSet:\t".$ra_dset->name);
		next; #$ra_fset
	  }

	  my $ra_rset = $sset[0];
	
	  foreach my $state(@$rset_states){
	  
		if(! $ra_rset->has_status($state)){
		  #Can we just update and warn here?
		  #Or do this separately in case we want some control over this?
		  $self->report("WARNING:\tUpdating ResultSet ".$ra_rset->name." with status $state");

		  $sql = 'INSERT into status select '.$ra_rset->dbID.
			", 'result_set', status_name_id from status_name where name='$state'";
		  $self->db->dbc->db_handle->do($sql);
		}
	  }
	}
  }

  return;
} # End of validate_DataSets






sub analyse_and_optimise_tables{
  my $self = shift;

  #myisamchk --analyze. or analyze statement

  if($self->{'skip_analyse'}){
	$self->log_header('Skipping analyse/optimise tables');
	return;
  }

  $self->log_header("Analysing and optimising tables");
  

  my $sql = 'show tables;';
  my @tables = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};	  
  map $_ = "@{$_}", @tables;
  my  $analyse_sql = 'analyze table ';
  my $optimise_sql = 'optimize table ';   



  foreach my $table(@tables){
	$self->log("Analysing and optimising $table");



	#Remove analyse as optimise does everything this does
	my @anal_info = @{$self->db->dbc->db_handle->selectall_arrayref($analyse_sql.$table)};

	foreach my $line_ref(@anal_info){
	  my $status = $line_ref->[3];
	  $self->report("FAIL: analyse $table status $status") if (!($status eq 'OK' || $status eq 'Table is already up to date'));
	}

	my @opt_info = @{$self->db->dbc->db_handle->selectall_arrayref($optimise_sql.$table)};

	foreach my $line_ref(@opt_info){
	  
	  my $status = $line_ref->[3];
	  $self->report("FAIL: optimise $table status $status") if (!( $status eq 'OK' || $status eq 'Table is already up to date'));
	}
	
  }

  return;
}# end of analyse_and_optimise_tables


sub clean_xrefs{
  my ($self) = @_;

  if($self->{'skip_xrefs'}){
	$self->log_header('Skipping clean_xrefs');
	return;
  }

  $self->log_header("Cleaning unlinked xref records");
 
  my $sql = 'DELETE x FROM xref x LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id WHERE ox.xref_id IS NULL';
  #Should this also take accoumt of unmapped_objects?
  #No, as unmapped_object doesn't use xref, but probably should

   my $row_cnt = $self->db->dbc->do($sql);

  $self->reset_table_autoinc('xref', 'xref_id', $self->db);
  $row_cnt = 0 if $row_cnt eq '0E0';
  $self->log("Deleted $row_cnt unlinked xref records");


  #Now remove old edbs
  $self->log_header("Cleaning unlinked external_db records");

  #Need to account for xref and unmapped_object here
  $sql = 'DELETE edb FROM external_db edb '.
	'LEFT JOIN xref x ON x.external_db_id = edb.external_db_id '.
	  'LEFT JOIN  unmapped_object uo ON uo.external_db_id=edb.external_db_id '.
		'WHERE x.external_db_id IS NULL and uo.external_db_id is NULL';
  $row_cnt = $self->db->dbc->do($sql);

  $self->reset_table_autoinc('external_db', 'external_db_id', $self->db);
  $row_cnt = 0 if $row_cnt eq '0E0';
  $self->log("Deleted $row_cnt unlinked external_db records");


  #Shouldn't clean orphaned oxs here as this means a rollback been done underneath the ox data
  #or we have xref_id=0!
  #Leave this to HC?

  



  return;
}


1;
