=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::HealthChecker

=head1 SYNOPSIS

=head1 DESCRIPTION

B<This program> provides several methods to health check and update tables prior to
release. Using the updte_DB_for_release method runs the following:

  validate_new_seq_regions    - _pre_stores seq_region & coord_system info from new core DB
  check_meta_species_version  - Validates meta species and version wrt dbname
  set_current_coord_system    - Updates coord_system.is_current to 1 for current schema_build (required for mart)
  update_meta_coord           - Regenerates meta_coord.max_length values (required for Slice range queries)
  clean_xrefs                 - Removes old unused xref and external_db records
  log_data_sets               - Logs all DataSets
  analyse_and_optimise_tables - Does what is says


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::HealthChecker;

use strict;
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::ProbeFeature;


use base qw( Bio::EnsEMBL::Funcgen::Utils::Helper );


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
  my ($db, $builds, $skip_mc, $check_displayable, $skip_analyse, $meta_coord_tables, $skip_xrefs, $fix) =
	rearrange(['DB', 'BUILDS', 'SKIP_META_COORD', 'CHECK_DISPLAYABLE', 'SKIP_ANALYSE', 'META_COORD_TABLES', 'SKIP_XREF_CLEANUP', 'FIX'], @_);


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
  $self->{fix} = $fix;

  if(defined $meta_coord_tables){

	throw('-skip_meta_coord is set, Cannot build meta_coord entries for tables '.join(', ', @$meta_coord_tables));

	if(! ref($meta_coord_tables) eq 'ARRAY'){
	  throw('-meta_coord_tables parameter must be an array ref');
	}

	@{$self->{'meta_coord_tables'}} = @$meta_coord_tables;
  }

  return $self;
}

sub db{  return $_[0]->{db};  }

sub fix{ return $_[0]->{fix}; }

=head2 update_db_for_release

  Arg[0]     :
  Example    :
  Description: Wrapper method to perform all common update functions
  Returntype :
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub update_db_for_release {
  my ($self, @args) = @_;

  if(@args){
  }

  #do seq_region_update to validate dnadb first
  #hence avoiding redoing longer methods

  $self->validate_new_seq_regions;
  $self->set_current_coord_system;
  $self->update_meta_coord;
#  $self->clean_xrefs;
#  $self->log_data_sets();
  $self->analyse_and_optimise_tables;#ALWAYS LAST!!

  $self->log_header('??? Have you dumped/copied GFF dumps ???');

  #Log footer? Pass optional counts hash?
  $self->log('Finished updating '.$self->{'dbname'}." for release\n\n");
}

sub validate_new_seq_regions {
  my $self = shift;

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

    if($slice->start != 1){
      $self->log("Reslicing slice:\t".$slice->name);
	    #we must have some sort of PAR linked region i.e. Y
      $slice = $slice_adaptor->fetch_by_region($slice->coord_system_name(), $slice->seq_region_name());
    }


	  #we need test if it needs doing first?
	  #we would need to test for the coord_systems outside of this loop
	  #and then for each seq_region inside the loop if the coord_system is present

    $self->log("_pre_storing seq_region info for slice:\t".$slice->name());

    $pf_adaptor->_pre_store(
      Bio::EnsEMBL::Funcgen::ProbeFeature->new(-slice => $slice,
                                               -start => 0,
                                               -end   => 0,
                                               -strand => 0     ) );
	  #This will create a meta_coord entry of max_length 1 for features which have an absent meta_coord entry
  }

  $self->log("Finished validating seq_regions\n");
  return;
}


sub update_meta_coord{
  my ($self, @table_names) = @_;

  if($self->{'skip_meta_coord'}){
    $self->log("Skipping meta_coord update\n");
    return;
  }

  $self->log_header('Updating meta_coord table');

  my $species_id           = $self->db->species_id;

  #set default table_name
  if(! @table_names || scalar(@table_names) == 0){
	  #Can we do this via DBAdaptor and get all available adaptors
    #which are BaseFeatureAdaptors then grab the first table name

    if(defined $self->{'meta_coord_tables'}){
      @table_names = @{$self->{'meta_coord_tables'}};
    }
    else{#default
      @table_names = qw(
	regulatory_feature
	probe_feature
	external_feature
	annotated_feature
	mirna_target_feature
      );
    }
  }

  #backup meta coord
  my $cmd = $self->{'mysql_connect_string'}.' -e "SELECT * FROM meta_coord" > '.$self->{'dbname'}.'meta_coord.backup';
  run_system_cmd($cmd);
	$self->log('Original meta_coord table backed up in '. $self->{'dbname'}.'.meta_coord.backup');

  #Update each max_length for table_name and coord_system
  my $sql;

  foreach my $table_name(@table_names){
    $sql = 'select distinct(cs.name), mc.coord_system_id, cs.version, mc.max_length from '.
    "coord_system cs, meta_coord mc where mc.table_name='$table_name' and ".
    "mc.coord_system_id=cs.coord_system_id and cs.species_id = $species_id";

    $self->log('');
    $self->log("Updating meta_coord max_length for $table_name:");
    $self->log("name\tcoord_system_id\tversion\tmax_length");

  	#can we test for emtpy array here? Then skip delete.
    my @info = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
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

#    $sql = 'select distinct s.coord_system_id from coord_system cs, seq_region s, '.
#     "$table_name t WHERE t.seq_region_id = s.seq_region_id and s.coord_system_id = ".
#     "cs.coord_system_id and cs.species_id = $species_id";
    $sql = "
      SELECT DISTINCT
        coord_system_id
      FROM
        seq_region
      WHERE
        seq_region_id
      IN (
          SELECT DISTINCT
            seq_region_id
          FROM
            $table_name
          )";

    my @cs_ids = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
	  #Convert single element arrayrefs to scalars
    map $_ = ${$_}[0], @cs_ids;

    #This should probably fail if features are present but don't match any seq_region_ids in the seq_region_table
    #i.e. they have been incorrectly imported directly into mysql
    $self->log("New max_lengths for $table_name are:");
    $self->log(join("\t", ('coord_system_id', 'max_length', 'longest record dbID')));

    foreach my $cs_id(@cs_ids){
      #This will always give a length of 1 even if there are no features present
      my @cs_lengths;

      #The probe_feature table is now too big to do this in one go
      #We need to break this down into sr_ids
      $sql = "SELECT distinct t.seq_region_id from $table_name t, seq_region sr ".
       "where t.seq_region_id=sr.seq_region_id and sr.coord_system_id=$cs_id";
      my @sr_ids = @{$self->db->dbc->db_handle->selectcol_arrayref($sql)};

      #Get longest feature for all seq_regions
      foreach my $sr_id(@sr_ids){
        $sql = "SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.${table_name}_id ".
         "FROM $table_name t WHERE t.seq_region_id = $sr_id ";
        $sql .= ' and t.window_size=0' if $table_name eq 'result_feature';
        $sql .= ' order by max desc limit 1';

        #Problem here is that DBs without 0 wsize result_feture entries will not get a meta_coord entry
        #We need to implement this in the _pre_store method too?
        my ($cs_length, $table_id);
        ($cs_length, $table_id) = $self->db->dbc->db_handle->selectrow_array($sql);
        push @cs_lengths, [$cs_length, $table_id] if $cs_length;
      }

      if(@cs_lengths){
        #This will now contain a list of arrays refs contain the max length and feature id for
		    #each seq_region in this coord_system. Now sort to get the longest
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


sub analyse_and_optimise_tables {
  my $self = shift;

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

  foreach my $table(@tables) {
    $self->log("Analysing and optimising $table");

    #Remove analyse as optimise does everything this does
    my @anal_info = @{$self->db->dbc->db_handle->selectall_arrayref($analyse_sql.$table)};

    foreach my $line_ref(@anal_info){
    my $status = $line_ref->[3];
    $self->report_fail("FAIL: analyse $table status $status") if (!($status eq 'OK' || $status eq 'Table is already up to date'));
    }

    my @opt_info = @{$self->db->dbc->db_handle->selectall_arrayref($optimise_sql.$table)};

    foreach my $line_ref(@opt_info){

    my $status = $line_ref->[3];
      $self->report_fail("FAIL: optimise $table status $status") if (!( $status eq 'OK' || $status eq 'Table is already up to date'));
    }
  }

  return;
}# end of analyse_and_optimise_tables


sub clean_xrefs {
  my ($self) = @_;

  if($self->{'skip_xrefs'}){
	$self->log_header('Skipping clean_xrefs');
	return;
  }

  $self->log_header("Cleaning unlinked xref records");

  # Deletes all entries from the xref table that don't link to an object xref.
  #
  my $sql = 'DELETE x FROM xref x LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id WHERE ox.xref_id IS NULL';

   my $row_cnt = $self->db->dbc->do($sql);

  $self->db->reset_table_autoinc('xref', 'xref_id');
  $row_cnt = 0 if $row_cnt eq '0E0';
  $self->log("Deleted $row_cnt unlinked xref records");

  #Now remove old edbs
  $self->log_header("Cleaning unlinked external_db records");

  # Seems to be delete all entries from the external_db table that neither
  # link to an xref nor to an unmapped object.
  #
  $sql = 'DELETE edb FROM external_db edb '.
	'LEFT JOIN xref x ON x.external_db_id = edb.external_db_id '.
	  'LEFT JOIN  unmapped_object uo ON uo.external_db_id=edb.external_db_id '.
		'WHERE x.external_db_id IS NULL and uo.external_db_id is NULL';
  $row_cnt = $self->db->dbc->do($sql);

  $self->db->reset_table_autoinc('external_db', 'external_db_id');
  $row_cnt = 0 if $row_cnt eq '0E0';
  $self->log("Deleted $row_cnt unlinked external_db records");

  return;
}


1;
