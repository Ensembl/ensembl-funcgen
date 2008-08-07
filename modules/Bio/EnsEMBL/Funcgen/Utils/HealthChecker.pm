=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::Helper
  
=head1 SYNOPSIS


 e.g. 


 my $object = Bio::EnsEMBL::Object->new
 (
     logging     => 1,
     log_file    => "/tmp/Misc.log",
     debug_level => 2,
     debug_file  => "/tmp/Misc.dbg",
 );

 $object->log("This is a log message.");
 $object->debug(1,"This is a debug message.");
 $object->system("rmdir /tmp/test");


 ----------------------------------------------------------------------------


=head1 OPTIONS

=over 8


=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> performs several health check and update methods prior to release.

=cut

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


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

################################################################################

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
    
  #validate and set type, analysis and feature_set here
  my ($db, $builds, $skip_mc, $check_displayable, $skip_analyse, $meta_coord_tables) = 
	rearrange(['DB', 'BUILDS', 'SKIP_META_COORD', 'CHECK_DISPLAYABLE', 'SKIP_ANALYSE', 'META_COORD_TABLES'], @_);
  
  
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
  $self->{'builds'} = (scalar(@$builds)>0) ? $builds : ['DEFAULT'];
  $self->{'skip_meta_coord'} = $skip_mc;

  if(defined $meta_coord_tables){

	throw('-skip_meta_coord is set, Cannot build meta_coord entries for tables '.join(', ', @$meta_coord_tables));

	if(! ref($meta_coord_tables) eq 'ARRAY'){
	  throw('-meta_coord_tables parameter must be an array ref');
	}

	@{$self->{'meta_coord_tables'}} = @$meta_coord_tables;
  }

  $self->{'skip_analyse'} = $skip_analyse;
  $self->{'check_displayable'} = $check_displayable;
  
  return $self;
}

sub db{
  my ($self) = @_;
  
  return $self->{'db'};
}


#wrapper method

=head2

  Arg[0]     : boolean - 
  Arg[1]     : boolean -
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : StableArg[0] - boolean force

=cut

sub update_db_for_release{
  my ($self, @args) = @_;
  
  if(@args){
  }
	
  #do seq_region_update to validate dnadb first
  #hence avoiding redoing longer methods
  $self->validate_new_seq_regions;#($force_srs);
  $self->update_meta_schema_version;
  $self->update_meta_coord;
  $self->check_meta_strings;
  $self->analyse_tables;
  

  $self->log('??? Have you dumped/copied GFF dumps ???');
  $self->log("??? Have you diff'd the sql for each species vs. a fresh schema ???");
  $self->log('Need to implement check meta string check');

  #Log footer? Pass optional counts hash?
  $self->log('Finished updating '.$self->{'dbname'}." for release\n\n");
}

sub validate_new_seq_regions{
  my ($self, $force) = @_;
  
  
  
  #do we need to add the none default levels here?
  #or are we only bothered about those which constitute the toplevel?

  #To make sure we have all the correct levels in eFG we need to get all the names.
  #then get all by name from the core db and set them as the dnadb.
  # we also need to get all the toplevel seq_regions and store them in the seq_region table
  #use BaseFeatureAdaptor::_pre_store with and array of pseudo feature on each top level slice

  #Validate the efgdb and dnadb schema version are the same first
  
  if(! $force){
	my $efgdb_sm = join('_', @{$self->get_schema_and_build($self->{'dbname'})});
	my $dnadb_sm = join('_', @{$self->get_schema_and_build($self->{'dbname'})});
	
	if($efgdb_sm ne $dnadb_sm){
	  $self->log("WARNING Skipped validate_new_seq_regions as schema_versions are mismatched:\t".
				 "efgdb $efgdb_sm\tdnadb $dnadb_sm");
	  return 0;
	}
  }
  
  my $pf_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $slice_adaptor = $self->db->get_SliceAdaptor();
  
  $self->log_header('Validating new coord_systems/seq_regions');
  
  foreach my $build(@{$self->{'builds'}}){
	
	$self->log("Importing seq_region/coord_system info for build:\t".$build);
	
	foreach my $slice(@{$slice_adaptor->fetch_all('toplevel', $build, 1)}){
	  #1 is non-reference flag, essential for haplotype regions

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
  }

  $self->log("Finished validating seq_regions\n");
  
  return;
}



sub update_meta_schema_version{
  my ($self) = @_;
  
  my $schema_version = $self->get_schema_and_build($self->{'dbname'})->[0];
  
  

  my $sql = 'DELETE from meta where meta_key="schema_version"';
  $self->db->dbc->db_handle->do($sql);
  $sql = "INSERT into meta(meta_key, meta_value) values ('schema_version', '$schema_version')";
 
  $self->db->dbc->db_handle->do($sql);

  $self->log_header("Updated meta.schema_version to $schema_version");

}




sub update_meta_coord{
  my ($self, @table_names) = @_;
  
  if($self->{'skip_meta_coord'}){
	$self->log("Skipping meta_coord update\n");
	return;
  }
  
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
	$self->log_header('Updating meta_coord table. Original meta_coord table backed up in '
			   . $self->{'dbname'}.'.meta_coord.backup');
  }
  

  #Update each max_length for table_name and coord_system
  foreach my $table_name(@table_names){
	my $sql1 = "select distinct(cs.name), mc.coord_system_id, cs.version, mc.max_length from coord_system cs, meta_coord mc where mc.table_name='$table_name' and mc.coord_system_id=cs.coord_system_id";
	
	$self->log('');
	$self->log("Updating meta_coord max_length for $table_name:\n\tname\tcoord_system_id\tversion\tmax_length");
	
	#can we test for emtpy array here? Then skip delete.
	
	my @info = @{$self->db->dbc->db_handle->selectall_arrayref($sql1)};
	
	#log this
	map {$self->log("\t".join("\t", @{$_})."\n")} @info;
	
	# Clean old entries
	$self->log("Deleting old meta_coord entries");
	my $sql = "DELETE FROM meta_coord WHERE table_name ='$table_name'";
	$self->db->dbc->db_handle->do($sql);
	
	# Generate new max_lengths
	$self->log("Generating new max_lengths");
	
	#Is this query running for each redundant cs_id?
	#would it be more efficient to retrieve the NR cs_ids first and loop the query for each cs_id?
	
	#Can we get the dbID of the largest feature for ease of checking?
	#This won't work as we're grouping by coord_system
	#would need to select distinct coord_system_id for table first
	#This may well slow down quite a bit doing it this way
	
	$sql = "select distinct s.coord_system_id from seq_region s, $table_name t WHERE t.seq_region_id = s.seq_region_id";
	my @cs_ids = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
	#Convert single element arrayrefs to scalars
	map $_ = ${$_}[0], @cs_ids;

	$self->log("New max_lengths for $table_name are:");
	
	#wtf?
	#map {$self->log(join("\t", @{$_}))} ['coord_system_id', 'max_length', 'longest record dbID'];
	$self->log(join("\t", ('coord_system_id', 'max_length', 'longest record dbID')));


	foreach my $cs_id(@cs_ids){
	  #This will always give a length of 1 even if there are no features present

	  $sql = "SELECT s.coord_system_id, (t.seq_region_end - t.seq_region_start + 1 ) as max, t.${table_name}_id "
			. "FROM $table_name t, seq_region s "
			  . "WHERE t.seq_region_id = s.seq_region_id "
				. "and s.coord_system_id=${cs_id} order by max desc limit 1";


	  @info = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
	  #Convert one multi element array_ref into array
	  @info = @{$info[0]};
	  $self->log(join("\t\t", @info));

	  $sql = "INSERT INTO meta_coord values(\"${table_name}\", \"${cs_id}\", \"$info[1]\")";


	  #$sql = "INSERT INTO meta_coord "
	  #. "SELECT '$table_name', s.coord_system_id, "
	  #	. "MAX( t.seq_region_end - t.seq_region_start + 1 ) "
	  #	  . "FROM $table_name t, seq_region s "
	  #		. "WHERE t.seq_region_id = s.seq_region_id "
	  #		  . "GROUP BY s.coord_system_id";
	  
		$self->db->dbc->db_handle->do($sql);
	}

	#$self->log("New max_lengths for $table_name are:");
  
	#@info = @{$self->db->dbc->db_handle->selectall_arrayref($sql1)};
	
	#map {$self->log(join("\t", @{$_})."\n")} @info;
  }  
  
  $self->log("Finished updating meta_coord max_lengths\n");
  
  return;
}


sub check_meta_strings{
  my ($self, $all_builds) = @_;
  
  #update flag?

  my @regf_fsets;
  my $passed = 1;
  my $fset_a = $self->db->get_featureSetAdaptor;
  my $mc = $self->db->get_MetaContainer;
  my $regf_a = $self->db->get_RegulatoryFeatureAdaptor;

  if($all_builds){
	@regf_fsets = @{$fset_a->fetch_all_by_type('regulatory')};
  }else{
	my $fset = $fset_a->fetch_by_name('RegulatoryFeatures');
	push @regf_fsets, $fset if defined $fset;
  }
  
  my @meta_keys = ('regulatory_string_feature_set_id', 'regulatory_string_feature_type_id');

  #What about anchor/seed sets?

  if(scalar(@regf_fsets) == 0){
	$self->log_header("Found no RegulatoryFeature sets for meta_string check");
  }
  else{

	$self->log_header("Validating meta entries for FeatureSets:\t".join("\t", (map $_->name, @regf_fsets)));
	
	
	#How do we validate this?
	#Check all feature_sets exist
	#Pull back some features from a test slice and check the number of bits match.
  #Check the feature_type string exists and matches else create.
	
	
	foreach my $fset(@regf_fsets){
	#get version number of build
	  my (undef, $build_version) = split/v/, $fset->name;
	$build_version = (defined $build_version) ? '_v'.$build_version : '';
	my $fset_string_key = 'regulatory_string_feature_set_id'.$build_version;
	my $ftype_string_key = 'regulatory_string_feature_type_id'.$build_version;
	my $fset_string = $mc->list_value_by_key($fset_string_key)->[0];
	my $ftype_string = $mc->list_value_by_key($ftype_string_key)->[0];

	$self->log('Validating '.$fset->name.":\n\t$fset_string_key($fset_string) vs $ftype_string_key($ftype_string)");

	#Test fset vs ftype string
	if(! defined $fset_string && ! defined $ftype_string){
	  $self->log("FAIL:\tNo $fset_string_key or $ftype_string_key found in meta table");
	}
	elsif(! defined $fset_string){
	  $self->log("FAIL:\tNo $fset_string_key found in meta table");
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
		$self->log("WARNING:\tNo $ftype_string_key found in meta table, will update using $fset_string_key");
	  }
	  

	  #Now need to work backwards through ftypes to remove pseudo ftypes before validating
	  #New string should be A,A,A;S,S,S,S,S,S;P,P,P
	  #Where A is and Anchor/Seed set
	  #S is a supporting set
	  #P is a pseudo feature type e.g. TSS proximal
	  

	  if(scalar(@fset_ids) != scalar(@ftype_ids)){
		$self->log("FAIL:\tLength mismatch between $fset_string_key and $ftype_string_key");
	  }

	  foreach my $i(0..$#fset_ids){
		my $supporting_set_id = $fset_ids[$i];
		my $sset = $fset_a->fetch_by_dbID($supporting_set_id);

		if(! defined $sset){
		  $self->log("FAIL:\t$fset_string_key $supporting_set_id does not exist in the DB");
		}
		else{
		  #test/build ftype string
		  
		  if(defined $ftype_string){
			
			if($sset->feature_type->dbID != $ftype_ids[$i]){
			  $ftype_fail = 1;
			  $self->log("FAIL:\t$fset_string_key $supporting_set_id(".$sset->name.") FeatureType(".$sset->feature_type->name.") does not match $ftype_string_key $ftype_ids[$i]");
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
		$self->db->dbc->db_handle->do("INSERT into meta values(NULL, '$ftype_string_key', '$new_ftype_string')");
	  }
	  elsif($ftype_fail){
		$self->log("FAIL:\t$ftype_string_key($ftype_string) does not match $fset_string_key types($new_ftype_string)");
	  }


	  #Finally validate versus a reg feat
	  #Need to change this to ftype string rather than fset string?

	  my ($regf_dbID) = @{$self->db->dbc->db_handle->selectrow_arrayref('select regulatory_feature_id from regulatory_feature where feature_set_id='.$fset->dbID.' limit 1')};
	  
	  if(! defined $regf_dbID){
		$self->log("FAIL:\tNo RegulatoryFeatures found for FeatureSet ".$fset->name);
	  }
	  else{
		my $rf_string = $regf_a->fetch_by_dbID($regf_dbID)->{'display_label'};#Direct access to avoid feature type
			
		if(length($rf_string) != scalar(@fset_ids)){
		  $self->log("FAIL:\tRegulatory string length mismatch between RegulatoryFeature($regf_dbID) and $fset_string_key:\n$rf_string(".length($rf_string).")\n$fset_string(".scalar(@fset_ids).")");
		}
	  }
	}
  }
  }
  return;
}


#Change this to log sets and incorporate RegFeat FeatureSet as standard


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
	$self->log_set("Found DataSet:\t\t", $dset) ;

	my $fset = $dset->product_FeatureSet;
	$self->log_set("Product FeatureSet:\t", $fset) if $fset;
	
	#my @supporting_sets = @{$dset->get_supporting_sets};

	if(my @supporting_sets = @{$dset->get_supporting_sets}){
	  map $self->log_set("SupportingSet:\t\t", $_), @supporting_sets;
	}
  }

  $self->log_header('Checking Regulatory FeatureSets');

  
  
  

  return;
}

sub log_set{
  my ($self, $text, $set) = @_;
  
  #if(! $set->isa('Bio::EnsEMBL::Funcgen::DataSet')){
#	$text .= $set->set_type.":\t";
#  }
 
  $text .= $set->display_label.'('.$set->name.')';
  $text .= "\tDISPLAYABLE" if($set->is_displayable);
  $self->log($text);

  return;
}

sub analyse_tables{
  my $self = shift;

  #myisamchk --analyze. or analyze statement

  if($self->{'skip_analyse'}){
	$self->log('Skipping analyse tables');
	return;
  }


  

  my $sql = 'show tables;';
  my @tables = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};	  
  map $_ = "@{$_}", @tables;
  $sql = 'analyze table ';


  foreach my $table(@tables){
	$self->log("Analysing table $table:");

	my @anal_info = @{$self->db->dbc->db_handle->selectall_arrayref($sql.$table)};

	foreach my $line_ref(@anal_info){
	  $self->log(join("\t", @$line_ref));
	}
  }

  return;
}


### Check for regulatory meta entries for all regulatory feature_sets


1;
