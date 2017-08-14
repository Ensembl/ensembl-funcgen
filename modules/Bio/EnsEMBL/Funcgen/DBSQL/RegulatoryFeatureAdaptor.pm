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

  Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryFeatureAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'RegulatoryFeature');
  my $regulatory_feature = $regulatory_feature_adaptor->fetch_by_stable_id('ENSR00000000011');

  print 'Stable id:        ' . $regulatory_feature->stable_id                              . "\n";
  print 'Analysis:         ' . $regulatory_feature->analysis->logic_name                   . "\n";
  print 'Feature type:     ' . $regulatory_feature->feature_type->name                     . "\n";
  print 'Epigenome count:  ' . $regulatory_feature->epigenome_count                        . "\n";
  print 'Slice name:       ' . $regulatory_feature->slice->name                            . "\n";
  print 'Coordinates:      ' . $regulatory_feature->start .' - '. $regulatory_feature->end . "\n";
  print 'Regulatory build: ' . $regulatory_feature->get_regulatory_build->name             . "\n";

=head1 DESCRIPTION

  The RegulatoryFeatureAdaptor is a database adaptor for storing and retrieving
  RegulatoryFeature objects. The FeatureSet class provides convenient wrapper
  methods to the Slice functionality within this adaptor.

=cut
package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Data::Dumper;
use DBI qw(:sql_types);

# One day:
# use base 'Bio::EnsEMBL::Feature';
use Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;
use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor);


=head2 fetch_by_stable_id

  Arg [1]    : String $stable_id - The stable id of the regulatory feature to retrieve
  Arg [2]    : optional - Bio::EnsEMBL::FeatureSet
  Example    : my $rf = $rf_adaptor->fetch_by_stable_id('ENSR00000309301');
  Description: Retrieves a regulatory feature via its stable id.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
    my $self      = shift;
    my $stable_id = shift;
    
    my $constraint = "rf.stable_id = ? and rb.is_current = true";
    $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
    
    my ($regulatory_feature) = @{$self->generic_fetch($constraint)};
    return $regulatory_feature;
}

sub fetch_by_stable_id_RegulatoryBuild {
    my $self             = shift;
    my $stable_id        = shift;
    my $regulatory_build = shift;

    return $self->_fetch_by_stable_id_regulatory_build_id($stable_id, $regulatory_build->dbID);
}

sub _fetch_by_stable_id_regulatory_build_id {
    my $self      = shift;
    my $stable_id = shift;
    my $regulatory_build_id = shift;

    my $constraint = "rf.stable_id = ? and rf.regulatory_build_id = ?";
    $self->bind_param_generic_fetch($stable_id,           SQL_VARCHAR);
    $self->bind_param_generic_fetch($regulatory_build_id, SQL_VARCHAR);
    
    my ($regulatory_feature) = @{$self->generic_fetch($constraint)};
    return $regulatory_feature;
}

sub fetch_Iterator {
    my $self      = shift;
    my $constraint = "regulatory_build.is_current = true";
   return $self->_generic_fetch_Iterator(10, $constraint);
}

sub fetch_Iterator_by_RegulatoryBuild {
    my $self             = shift;
    my $regulatory_build = shift;
    my $constraint = "regulatory_feature.regulatory_build_id = " . $regulatory_build->dbID;
    return $self->_generic_fetch_Iterator(10, $constraint);
}

# Nicked from Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor and customised.
#
sub _generic_fetch_Iterator {

    my ($self, $cache_size, $full_constraint) = @_;

    # prepare and execute a query to fetch all dbIDs
    my $sth = $self->prepare(qq{
        SELECT      regulatory_feature_id
        FROM        regulatory_feature JOIN regulatory_build using (regulatory_build_id)
        WHERE       $full_constraint
    });
    $sth->execute;

    my $regulatory_feature_id;
    $sth->bind_columns(\$regulatory_feature_id);

    $cache_size ||= 1000;
    
    my @cache;

    my $items_to_fetch = 1;

    return Bio::EnsEMBL::Utils::Iterator->new(sub{

        if (@cache == 0 && $items_to_fetch) {
            
            # our cache is empty, and there are still items to fetch, so
            # fetch the next chunk of dbIDs and create objects from them
	    #
            my @dbIDs;

            my $item_count = 0;

            while( $sth->fetch ) {

                push @dbIDs, $regulatory_feature_id;
                if (++$item_count == $cache_size) {
                    # we have fetched a cache's worth of dbIDs, so flag that
                    # there are still items to fetch and last out of the loop
                    $items_to_fetch = 1;
                    last;
                }
                # if this is the last row, this flag will be 0 outside the loop
                $items_to_fetch = 0;
            }
            $sth->finish unless $items_to_fetch;
            @cache = @{ $self->fetch_all_by_dbID_list(\@dbIDs) } if @dbIDs;
        }
        return shift @cache;
    });
}

sub _fake_multicell_activity {

  my $self = shift;
  my $actual_regulatory_activity = shift;
  
  my $multicell_regulatory_activity = Bio::EnsEMBL::Funcgen::RegulatoryActivity->new;
  $multicell_regulatory_activity->activity('ACTIVE');
  $multicell_regulatory_activity->_epigenome_id(undef);
  $multicell_regulatory_activity->_is_multicell(1);

  my $multicell_regulatory_evidence = Bio::EnsEMBL::Funcgen::RegulatoryEvidence->new;
  $multicell_regulatory_evidence->db($self->db);

  foreach my $current_regulatory_activity (@$actual_regulatory_activity) {
  
    my $regulatory_evidence = $current_regulatory_activity->regulatory_evidence;
  
    $multicell_regulatory_evidence->add_supporting_annotated_feature_id(
      $regulatory_evidence->supporting_annotated_feature_ids
    );
    $multicell_regulatory_evidence->add_supporting_motif_feature_id(
      $regulatory_evidence->supporting_motif_feature_ids
    );
  }
  $multicell_regulatory_activity->regulatory_evidence($multicell_regulatory_evidence);
  return $multicell_regulatory_activity;
}

sub _true_tables {
  return (
    [ 'regulatory_feature',  'rf' ],
    [ 'regulatory_build',    'rb' ],
  );
}

sub _columns {
  my $self = shift;

  return qw(
    rf.regulatory_feature_id
    rf.seq_region_id
    rf.seq_region_start
    rf.seq_region_end
    rf.seq_region_strand
    rf.bound_start_length
    rf.bound_end_length
    rf.feature_type_id
    rf.stable_id
    rf.epigenome_count
    rb.regulatory_build_id
    rb.is_current
    rb.analysis_id
  );
}

# Prevent superclass from adding a clause that won't work.
sub _default_where_clause {
  return 'rf.regulatory_build_id = rb.regulatory_build_id';
}

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;
  
  my $sa = ($dest_slice) ? $dest_slice->adaptor->db->get_SliceAdaptor() : $self->db->dnadb->get_SliceAdaptor;

  my $feature_type_adaptor     = $self->db->get_FeatureTypeAdaptor();
  my $analysis_adaptor         = $self->db->get_AnalysisAdaptor();  
  my $regulatory_build_adaptor = $self->db->get_RegulatoryBuildAdaptor();
  
  my @feature_from_sth;
  my (%fset_hash, %slice_hash, %sr_name_hash, %sr_cs_hash, %ftype_hash);

  my (
    $sth_fetched_dbID,
    $seq_region_id,
    $sth_fetched_seq_region_start,
    $sth_fetched_seq_region_end,
    $sth_fetched_seq_region_strand,
    $sth_fetched_bound_start_length,
    $sth_fetched_bound_end_length,
    $sth_fetched_feature_type_id,
    $sth_fetched_stable_id,
    $sth_fetched_epigenome_count,
    $sth_fetched_rb_dbid,
    $sth_fetched_rb_is_current,
    $sth_fetched_analysis_id,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$seq_region_id,
    \$sth_fetched_seq_region_start,
    \$sth_fetched_seq_region_end,
    \$sth_fetched_seq_region_strand,
    \$sth_fetched_bound_start_length,
    \$sth_fetched_bound_end_length,
    \$sth_fetched_feature_type_id,
    \$sth_fetched_stable_id,
    \$sth_fetched_epigenome_count,
    \$sth_fetched_rb_dbid,
    \$sth_fetched_rb_is_current,
    \$sth_fetched_analysis_id,
  );

  my ($dest_slice_start, $dest_slice_end);
  my ($dest_slice_strand, $dest_slice_length, $dest_slice_sr_name);

  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
    $dest_slice_length  = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
  }
  
  my $project_slice_coordinates_to_destination_slice = sub {
  
    # If the destination slice starts at 1 and is forward strand, nothing needs doing
    unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {

      if ($dest_slice_strand == 1) {
        $sth_fetched_seq_region_start    = $sth_fetched_seq_region_start - $dest_slice_start + 1;
        $sth_fetched_seq_region_end      = $sth_fetched_seq_region_end   - $dest_slice_start + 1;
      } else {
        my $tmp_seq_region_start = $sth_fetched_seq_region_start;
        $sth_fetched_seq_region_start    = $dest_slice_end - $sth_fetched_seq_region_end       + 1;
        $sth_fetched_seq_region_end      = $dest_slice_end - $tmp_seq_region_start + 1;
        $sth_fetched_seq_region_strand   *= -1;
      }
    }
  };

  # The current regulatory feature that is being constructed
  my $regulatory_feature_under_construction;
  
  my $fetch_slice_with_cache = sub {
    my $seq_region_id = shift;
    
    my $slice = $slice_hash{'ID:'.$seq_region_id};

    if (!$slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{'ID:'.$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }
    my $seq_region_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};
    return ($slice, $seq_region_name, $sr_cs);
  };

  # Flag to indicate that this feature should be skipped. This can happen 
  # when a feature is not on the destination slice.
  #
  my $current_feature_not_on_destination_slice = undef;
  my %regulatory_build_cache;
  
  ROW: while ( $sth->fetch() ) {
  
    $ftype_hash{$sth_fetched_feature_type_id} = $feature_type_adaptor->fetch_by_dbID($sth_fetched_feature_type_id) 
      if ! exists $ftype_hash{$sth_fetched_feature_type_id};

    # Get the slice object
    my ($slice, $seq_region_name) = $fetch_slice_with_cache->($seq_region_id);
    
    if ($mapper) {
    
      # If we are here, that means that there is a feature on a seq region 
      # that is not toplevel.
      #
      # This is a data issue and should never happen.
      #s
      throw("There are features in the database that haven't been mapped to toplevel!");
    }
    
    # If a destination slice was provided convert the coords
    if ($dest_slice) {
    
      $project_slice_coordinates_to_destination_slice->();

      my $current_feature_not_on_destination_slice = 
        $sth_fetched_seq_region_end < 1 
        || $sth_fetched_seq_region_start > $dest_slice_length
        || ( $dest_slice_sr_name ne $seq_region_name );

      next ROW
      if ($current_feature_not_on_destination_slice);

      $slice = $dest_slice;
    }

    my $regulatory_feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast({
	'start'             => $sth_fetched_seq_region_start,
	'end'               => $sth_fetched_seq_region_end,
	'_bound_lengths'    => [$sth_fetched_bound_start_length, $sth_fetched_bound_end_length],
	'strand'            => $sth_fetched_seq_region_strand,
	'slice'             => $slice,
	'_analysis_id'      => $sth_fetched_analysis_id,
	'adaptor'           => $self,
	'dbID'              => $sth_fetched_dbID,
	'feature_type'      => $ftype_hash{$sth_fetched_feature_type_id},
	'stable_id'         => $sth_fetched_stable_id,
	'epigenome_count'   => $sth_fetched_epigenome_count,
	'regulatory_build_id' => $sth_fetched_rb_dbid,
	
	});
    push @feature_from_sth, $regulatory_feature;
  }
  return \@feature_from_sth;
}

=head2 store

  Args       : Array of Bio::EnsEMBL::Funcgen::RegulatoryFeature objects
  Example    : $regulatory_feature_adaptor->store(@regulatory_features);
  Description: Stores given RegulatoryFeature objects in the database. Sets 
		dbID and adaptor on the objects that it stores.
  Returntype : Listref of stored RegulatoryFeatures
  Exceptions : Throws, if a list of RegulatoryFeature objects is not provided or if
               the Analysis, Epigenome and FeatureType objects are not attached or stored.
               Throws, if analysis of set and feature do not match
               Warns if RegulatoryFeature already stored in DB and skips store.
  Caller     : Regulatory Build
  Status     : Stable

=cut

sub store {
  my ($self, @regulatory_feature) = @_;

  if (scalar(@regulatory_feature) == 0) {
    throw('Must call store with a list of RegulatoryFeature objects');
  }
  foreach my $rf (@regulatory_feature) {
    if( ! ref $rf || ! $rf->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature') ) {
      throw('Feature must be an RegulatoryFeature object');
    }
  }
  
  my $sth_store_regulatory_feature = $self->prepare("
    INSERT INTO regulatory_feature (
      seq_region_id,
      seq_region_start,
      seq_region_end,
      bound_start_length,
      bound_end_length,
      seq_region_strand,
      feature_type_id,
      stable_id,
      epigenome_count,
      regulatory_build_id
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
  );
  # When loading the regulatory build, this would lead to many error messages 
  # about duplicate entries being printed to screen. Errors are either handled
  # or rethrown.
  #
  $sth_store_regulatory_feature->{PrintError} = 0;

  my $sth_store_regulatory_evidence = $self->prepare("
    INSERT INTO regulatory_evidence (
      regulatory_feature_id, 
      attribute_feature_id, 
      attribute_feature_table
    ) VALUES (?, ?, ?)"
  );
  
  my $sth_regulatory_activity = $self->prepare("
    INSERT INTO regulatory_activity (
      regulatory_feature_id,
      epigenome_id,
      activity
    ) VALUES (?,?,?);
  ");

  my $db = $self->db();

  foreach my $current_regulatory_feature (@regulatory_feature) {

    my $seq_region_id;
    ($current_regulatory_feature, $seq_region_id) = $self->_pre_store($current_regulatory_feature);
    $current_regulatory_feature->adaptor($self);

    $sth_store_regulatory_feature->bind_param( 1, $seq_region_id,                                   SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 2, $current_regulatory_feature->start,               SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 3, $current_regulatory_feature->end,                 SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 4, $current_regulatory_feature->bound_start_length,  SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 5, $current_regulatory_feature->bound_end_length,    SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 6, $current_regulatory_feature->strand,              SQL_TINYINT);
    $sth_store_regulatory_feature->bind_param( 7, $current_regulatory_feature->feature_type->dbID,  SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param( 8, $current_regulatory_feature->stable_id,           SQL_VARCHAR);
    $sth_store_regulatory_feature->bind_param( 9, $current_regulatory_feature->epigenome_count,     SQL_INTEGER);
    $sth_store_regulatory_feature->bind_param(10, $current_regulatory_feature->regulatory_build_id, SQL_INTEGER);
    
    eval {
      # Store and set dbID
      $sth_store_regulatory_feature->execute;
      $current_regulatory_feature->dbID( $self->last_insert_id );
    };
    if ($@) {
      my $error_message = $@;
      
      # If the regulatory feature already exists, the error message will look 
      # like this:
      #
      # DBD::mysql::st execute failed: Duplicate entry 
      # '179363-1579-0-2-200-ENSR00000000001-0-0' for key 
      # 'uniqueness_constraint_idx' at [..]/ensembl-funcgen/modules/Bio/EnsEMBL/Funcgen/DBSQL/RegulatoryFeatureAdaptor.pm line 539, <$fh> line 10001.
      #
      # It would be possible to check, whether the regulatory feature already 
      # exists in the database, but doing so would make storing slower. (I guess)
      #
      my $regulatory_feature_already_exists = $error_message =~ /uniqueness_constraint_idx/;
      
      # The uniqueness constraint is in place to avoid duplicate entries in 
      # the regulatory feature table. 
      #
      # The regulatory build script however creates a new regulatory feature
      # for every possible activity of a regulatory feature and stores that as
      # a new feature.
      #
      # The insertion of a duplicate regulatory feature is caught here. Then 
      # the existing regulatory feature is retrieved and the new activity is 
      # added to it.
      #
      if (! $regulatory_feature_already_exists) {
      
	# If the error message is about something else, then rethrow.
	#
	throw($error_message);
      }
      my $existing_regulatory_feature = $self->_fetch_by_stable_id_regulatory_build_id( 
	$current_regulatory_feature->stable_id, 
	$current_regulatory_feature->regulatory_build_id 
      );
      
      # This can happen during the regulatory build, when there are features,
      # but not stable ids yet. And the script is being rerun.
      #
      if (! defined $existing_regulatory_feature) {
        throw($error_message);
      }
      
      # Set the database id so the attributes and activities can be linked to this.
      $current_regulatory_feature->dbID( $existing_regulatory_feature->dbID );
    }
    
    if (! defined $current_regulatory_feature->regulatory_activity) {
      throw('Feature has no regulatory activity.');
    }
    if (ref $current_regulatory_feature->regulatory_activity ne 'ARRAY') {
      throw('Regulatory activity must be an array.');
    }
    
    use Data::Dumper;
    if (! defined $current_regulatory_feature->dbID) {
      throw(
	"Error storing the regulatory feature: "
	. Dumper($current_regulatory_feature)
      );
    }
    
    # Store the activities of the current regulatory feature in the various feature sets.
    #
    REGULATORY_ACTIVITY:
    foreach my $current_regulatory_activity (@{$current_regulatory_feature->regulatory_activity}) {

      next REGULATORY_ACTIVITY if ($current_regulatory_activity->_is_multicell);

      $sth_regulatory_activity->bind_param(1,  $current_regulatory_feature->dbID,           SQL_INTEGER);
      $sth_regulatory_activity->bind_param(2,  $current_regulatory_activity->_epigenome_id, SQL_INTEGER);
      $sth_regulatory_activity->bind_param(3,  $current_regulatory_activity->activity);

      eval {
	$sth_regulatory_activity->execute();
      };
      if ($@) {
	use Carp;
	$Data::Dumper::Maxdepth = 3;
	confess(
	  Dumper({
	    error => $@,
	    regulatory_activity => $current_regulatory_activity,
	    regulatory_feature => $current_regulatory_feature,
	  })
	);
      }
      
#       # Store the regulatory_evidence
#       #
#       # Note that the regulatory build script bypasses the api for loading 
#       # regulatory attributes, so this probably never gets called.
#       #
#       # That is a good thing, because this code links the attributes to the 
#       # regulatory features. In the new schema (v85 and above) regulatory
#       # attributes are linked to regulatory_feature_feature_sets.
#       #
#       my $regulatory_evidence = $current_regulatory_activity->get_RegulatoryEvidence;
# 
#       foreach my $id (@{$regulatory_evidence->supporting_motif_feature_ids}) {
# 
#         $sth_store_regulatory_evidence->bind_param(1, $current_regulatory_feature->dbID, SQL_INTEGER);
#         $sth_store_regulatory_evidence->bind_param(2, $id,  SQL_INTEGER);
#         $sth_store_regulatory_evidence->bind_param(3, 'motif', SQL_VARCHAR);
#         
#         $sth_store_regulatory_evidence->execute();
#         
#       }
#       foreach my $id (@{$regulatory_evidence->supporting_annotated_feature_ids}) {
# 
#         $sth_store_regulatory_evidence->bind_param(1, $current_regulatory_feature->dbID, SQL_INTEGER);
#         $sth_store_regulatory_evidence->bind_param(2, $id,  SQL_INTEGER);
#         $sth_store_regulatory_evidence->bind_param(3, 'annotated', SQL_VARCHAR);
#         
#         $sth_store_regulatory_evidence->execute();
#         
#       }
    }
  }
  return @regulatory_feature;
}

sub valid_activities {
  return ('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');
}

sub valid_activities_as_string {
  return join ', ', valid_activities;
}

sub is_valid_activity {
  my $self = shift;
  my $activity_to_test = shift;
  
  my @valid_activity = valid_activities;  
  foreach my $current_valid_activity (@valid_activity) {
    return 1 if ($activity_to_test eq $current_valid_activity);
  }
  return;
}

sub _make_arrayref_if_not_arrayref {
  my $obj = shift;
  my $obj_as_arrayref;
  
  if (ref $obj eq 'ARRAY') {
    $obj_as_arrayref = $obj;
  } else {
    $obj_as_arrayref = [ $obj ];
  }
  return $obj_as_arrayref;
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $regf_adaptor->fetch_all_by_Slice($slice);
  Description: Retrieves a list of features on a given slice, specific for the current
               default RegulatoryFeature set.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice {
  my ($self, $slice) = @_;
  
  return $self->_fetch_all_by_Slice_Epigenomes_Activity_RegulatoryBuild(
    $slice, undef, undef, undef
  );
}
sub fetch_all_by_Slice_RegulatoryBuild {
  my ($self, $slice, $regulatory_build) = @_;
  
  return $self->_fetch_all_by_Slice_Epigenomes_Activity_RegulatoryBuild(
    $slice, undef, undef, $regulatory_build
  );
}

sub fetch_all_by_Slice_FeatureSets {
  die("Regulatory features are no longer linked ot feature sets. Use fetch_all_by_Slice_Epigenomes instead.");
}

sub _fetch_all_by_Slice_Epigenomes_Activity_RegulatoryBuild {
  my ($self, $slice, $epigenomes, $activity, $selected_regulatory_build) = @_;

  if (defined $activity) {
    if (! $self->is_valid_activity($activity)) {
      die(
	qq(\"$activity\"is not a valid activity. Valid activities are: ) . valid_activities_as_string
      );
    }
  }
  
  #explicit super call, just in case we ever re-implement in here
#   my $all_regulatory_features = $self->SUPER::fetch_all_by_Slice($slice);
  my $all_regulatory_features = $self->SUPER::fetch_all_by_Slice_constraint($slice);
  
  
  if (defined $selected_regulatory_build) {
    #
    # Discard regulatory features that are not part of the selected regulatory build.
    #
    my $filtered_regulatory_features;
    REGULATORY_FEATURE: foreach my $current_regulatory_feature (@$all_regulatory_features) {
      if ($current_regulatory_feature->get_regulatory_build->dbID == $selected_regulatory_build->dbID) {
	push @$filtered_regulatory_features, $current_regulatory_feature;
	next REGULATORY_FEATURE;
      }
    }
    $all_regulatory_features = $filtered_regulatory_features;
  } else {
    #
    # Discard regulatory features that are not part of the current regulatory build.
    #
    my $filtered_regulatory_features;
    REGULATORY_FEATURE: foreach my $current_regulatory_feature (@$all_regulatory_features) {
      if ($current_regulatory_feature->get_regulatory_build->is_current) {
	push @$filtered_regulatory_features, $current_regulatory_feature;
	next REGULATORY_FEATURE;
      }
    }
    $all_regulatory_features = $filtered_regulatory_features;
  }

  if (defined $epigenomes) {
  
    $epigenomes = _make_arrayref_if_not_arrayref($epigenomes);

    my $filtered_regulatory_features;

    REGULATORY_FEATURE: foreach my $current_regulatory_feature (@$all_regulatory_features) {
      foreach my $current_epigenome (@$epigenomes) {
	if ($current_regulatory_feature->has_activity_in($current_epigenome)) {
	  push @$filtered_regulatory_features, $current_regulatory_feature;
	  next REGULATORY_FEATURE;
	}
      }
    }
    $all_regulatory_features = $filtered_regulatory_features;
  }
  
  if (defined $activity) {
  
    my $filtered_regulatory_features;
    
    REGULATORY_FEATURE: foreach my $current_regulatory_feature (@$all_regulatory_features) {
      if ($current_regulatory_feature->has_epigenomes_with_activity($activity)) {
	push @$filtered_regulatory_features, $current_regulatory_feature;
	next REGULATORY_FEATURE;
      }
    }
    $all_regulatory_features = $filtered_regulatory_features;
  }
  return $all_regulatory_features;
}

# sub _default_where_clause {
#   return 'ra.regulatory_feature_id = rf.regulatory_feature_id';
# }


=head2 fetch_all_by_attribute_feature

  Arg [1]    : Bio::Ensembl::Funcgen::AnnotatedFeature or MotifFeature
  Example    : my @regfs = @{$regf_adaptor->fetch_all_by_attribute_feature($motif_feature)};
  Description: Retrieves a list of RegulatoryFeatures which contain the given attribute feature.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature objects
  Exceptions : Throws is argument not valid
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_attribute_feature {
  my ($self, $attr_feat) = @_;

  my $attr_class = ref($attr_feat);
  
#   use Carp;
#   confess("This is never used.");
# 
  my %valid_attribute_features = (
    'Bio::EnsEMBL::Funcgen::MotifFeature'     => 'motif',
    'Bio::EnsEMBL::Funcgen::AnnotatedFeature' => 'annotated',
  );
  
  if(! exists $valid_attribute_features{$attr_class}) {
    throw("Attribute feature must be one of:\n\t".join("\n\t", keys(%valid_attribute_features)));
  }

  $self->db->is_stored_and_valid($attr_class, $attr_feat);
  my $attr_feat_table = $valid_attribute_features{$attr_class};

  my $rf_ids = $self->db->dbc->db_handle->selectall_arrayref(
    "select regulatory_feature_id "
    . "from regulatory_evidence join regulatory_activity using (regulatory_activity_id) join regulatory_feature using (regulatory_feature_id) join regulatory_build using (regulatory_build_id) "
    . "where attribute_feature_table='${attr_feat_table}' and attribute_feature_id=".$attr_feat->dbID. " and regulatory_build.is_current=1"
  );
  
  my @rf_ids_flattened = map { @$_ } @$rf_ids;
  my @results = $self->_fetch_by_dbID_list(@rf_ids_flattened);
  return \@results;
}

sub _fetch_by_dbID_list {
  my $self = shift;
  my @db_id = @_;
  
  my @fetched_objects;  
  foreach my $current_db_id (@db_id) {  
    my $object = $self->fetch_by_dbID($current_db_id);
    push @fetched_objects, $object;
  }
  return @fetched_objects;
}

1;
