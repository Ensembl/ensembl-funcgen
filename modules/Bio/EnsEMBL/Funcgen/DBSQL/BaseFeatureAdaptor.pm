#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
#
# Copyright (c) 2006 Ensembl
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor - An Base class for all
Funcgen FeatureAdaptors, redefines some methods to use the Funcgen DB

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for Funcgen feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods
common to all Funcgen feature adaptors.

=head1 CONTACT

Contact Ensembl development list for info: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA @EXPORT);
use strict;


use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
#use Bio::EnsEMBL::DBSQL::BaseAdaptor;#?
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

our $SLICE_FEATURE_CACHE_SIZE = 4;
our $MAX_SPLIT_QUERY_SEQ_REGIONS = 3;





#Added schema_build to feature cache key
#Do not remove

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
  my($self, $slice, $constraint, $logic_name) = @_;

  my @result;

  if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
    throw("Bio::EnsEMBL::Slice argument expected.");
  }

  $constraint ||= '';
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  #if the logic name was invalid, undef was returned
  return [] if(!defined($constraint));

  #check the cache and return if we have already done this query

  #Added schema_build to feature cache for EFG
  my $key = uc(join(':', $slice->name, $constraint, $self->db->_get_schema_build($slice->adaptor())));


  if(exists($self->{'_slice_feature_cache'}->{$key})) {
    return $self->{'_slice_feature_cache'}->{$key};
  }

  my $sa = $slice->adaptor();

  # Hap/PAR support: retrieve normalized 'non-symlinked' slices
  my @proj = @{$sa->fetch_normalized_slice_projection($slice)};

  if(@proj == 0) {
    throw('Could not retrieve normalized Slices. Database contains ' .
          'incorrect assembly_exception information.');
  }

  # Want to get features on the FULL original slice
  # as well as any symlinked slices

  # Filter out partial slices from projection that are on
  # same seq_region as original slice

  my $sr_id = $slice->get_seq_region_id();

  @proj = grep { $_->to_Slice->get_seq_region_id() != $sr_id } @proj;

  my $segment = bless([1,$slice->length(),$slice ],
                      'Bio::EnsEMBL::ProjectionSegment');
  push( @proj, $segment );


  # construct list of Hap/PAR boundaries for entire seq region
  my @bounds;
  my $ent_slice = $sa->fetch_by_seq_region_id($sr_id);
  $ent_slice    = $ent_slice->invert() if($slice->strand == -1);
  my @ent_proj  = @{$sa->fetch_normalized_slice_projection($ent_slice)};

  shift @ent_proj; # skip first
  @bounds = map {$_->from_start - $slice->start() + 1} @ent_proj;


  # fetch features for the primary slice AND all symlinked slices
  foreach my $seg (@proj) {
    my $offset = $seg->from_start();
    my $seg_slice  = $seg->to_Slice();

    my $features = $self->_slice_fetch($seg_slice, $constraint); ## NO RESULTS

    # if this was a symlinked slice offset the feature coordinates as needed
    if($seg_slice->name() ne $slice->name()) {

    FEATURE:
      foreach my $f (@$features) {
        if($offset != 1) {

          $f->{'start'} += $offset-1;
          $f->{'end'}   += $offset-1;
        }

        # discard boundary crossing features from symlinked regions
        foreach my $bound (@bounds) {
          if($f->{'start'} < $bound && $f->{'end'} >= $bound) {
            next FEATURE;
          }
        }

        $f->{'slice'} = $slice;
        push @result, $f;
      }
    }
    else {
      push @result, @$features;
    }
  }

  $self->{'_slice_feature_cache'}->{$key} = \@result;


  return \@result;
}




####Redefined for Funcgen to trust the slice the feature is generated on
####and store on a non-core/dna i.e. slice enabled DB
####Also handles multiple coord systems, generating new coord_system_ids as appropriate
####Need to think about retrieval, i.e. mapping coord_system_id back to the core/dnadb
####version may not mean same seq_region_ids?  So need to use name somehow, name may not be unique!!
####Using the schema version may guarantee stability of seq_region_ids

#
# Helper function containing some common feature storing functionality
#
# Given a Feature this will return a copy (or the same feature if no changes 
# to the feature are needed) of the feature which is relative to the start
# of the seq_region it is on. The seq_region_id of the seq_region it is on
# is also returned.
#
# This method will also ensure that the database knows which coordinate
# systems that this feature is stored in.
#

sub _pre_store {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand());


  my $db = $self->db();
  my $slice = $feature->slice();

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region
  if($slice->start != 1 || $slice->strand != 1) {
	  throw("You must generate your feature on a slice starting at 1 with strand 1");
	  #move feature onto a slice of the entire seq_region
	  #$slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
	  #                                         $slice->seq_region_name(),
	  #                                         undef, #start
	  #                                         undef, #end
	  #                                         undef, #strand
	  #                                         $slice->coord_system->version());
	  #$feature = $feature->transfer($slice);
	  #if(!$feature) {
	  #  throw('Could not transfer Feature to slice of ' .
	  #        'entire seq_region prior to storing');
	  #}
  }
  
  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;#from core/dnadb
 
  #retrieve corresponding Funcgen coord_system and set id in feature
  my $csa = $self->db->get_FGCoordSystemAdaptor();#had to call it FG as we were getting the core adaptor
  my $fg_cs = $csa->validate_coord_system($cs);
  $fg_cs = $csa->fetch_by_name_schema_build_version($cs->name(), $self->db->_get_schema_build($slice->adaptor()), $cs->version());
  $feature->coord_system_id($fg_cs->dbID());

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];


  #Need to do this for Funcgen DB
  my $mcc = $db->get_MetaCoordContainer();
  $mcc->add_feature_type($fg_cs, $tabname, $feature->length);

 # my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
  my $seq_region_id = $slice->get_seq_region_id();

  #would never get called as we're not validating against a different core DB
  #if(!$seq_region_id) {
  #  throw('Feature is associated with seq_region which is not in this dnadb.');
  #}

  return ($feature, $seq_region_id);
}





1;


