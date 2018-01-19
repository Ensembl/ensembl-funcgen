#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor - Funcgen Feature Adaptor base class

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for Funcgen feature adaptors. This base class is simply a way
to redefine some methods to use with the Funcgen DB.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

use vars qw(@EXPORT); #require Exporter is done in BaseAdaptor
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

### To do ###
#Add Translation method for fetch by FeatureSets methods?
#Test and implement feature mapping between coord systems
#Correct/Document methods!!!
#Implement externaldb db_name version registry

my %warnings;

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Wrapper method for core BaseAdaptor, build seq_region cache for features
  Returntype : ARRAYREF of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : FeatureAdaptor classes
  Status     : at risk

=cut

sub generic_fetch {
  my $self = shift;

  #need to wrap _generic_fetch to always generate the
  #seq_region_cache other wise non-slice based fetch methods fail

  #build seq_region cache here once for entire query
  #Using default schema_build here
  #So would need to set dnadb appropriately
  #This is cleaning tmp cache values such that
  #nested feature fetches cause failure e.g. when regfeats retrieve their reg_attrs
  #We need a way of always generating the tmp cache, or having it persist?
  #This is because we haven't built the tmp cache in non_slice based methods i.e. we haven't run get_seq_region_id_by_Slice
  #This is the same for all non-Slice based methods
  #And the is no way around it as we are not providing that info about the new DB by passing a slice!
  #The only way to get around this is to make the tmp_cache persistant

#   $self->build_seq_region_cache();

  return $self->SUPER::generic_fetch(@_);
}


# TODO compare this to core methods and remove new_assembly support?
# 

sub _pre_store {
  my ($self, $feature, $new_assembly) = @_;

  if (!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand());

  my $db = $self->db();
  my $slice = $feature->slice();

  if (!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region
  if ($slice->start != 1 || $slice->strand != 1) {

    #move feature onto a slice of the entire seq_region
    $slice = $slice->adaptor->fetch_by_region($slice->coord_system->name(),
                                              $slice->seq_region_name(),
                                              undef, #start
                                              undef, #end
                                              undef, #strand
                                              $slice->coord_system->version());
    $feature = $feature->transfer($slice);

    if (!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }
  
  my $seq_region_id = $slice->adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}

# May need to clone _get_by_Slice in here
# if coord system mapping is required




sub fetch_all_by_Gene_FeatureSets{
  my ($self, $gene, $fsets, $dblinks) = @_;

  if (! ( ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))) {
    throw("You must pass a valid Bio::EnsEMBL::Gene object");
  }

  my @features = @{$self->fetch_all_by_stable_Storable_FeatureSets($gene, $fsets)};

  if ($dblinks) {

    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      push @features, @{$self->fetch_all_by_Transcript_FeatureSets($transcript, $fsets, $dblinks)};
    }
  }

  return \@features;
}

sub fetch_all_by_Transcript_FeatureSets{
  my ($self, $transc, $fsets, $dblinks) = @_;

  if (! ( ref($transc) && $transc->isa('Bio::EnsEMBL::Transcript'))) {
    throw("You must pass a valid Bio::EnsEMBL::Transcript object");
  }


  my @features = @{$self->fetch_all_by_stable_Storable_FeatureSets($transc, $fsets)};

  if ($dblinks) {
    my $translation = $transc->translation;
    push @features, @{$self->fetch_all_by_stable_Storable_FeatureSets($translation, $fsets)} if $translation;
  }

  return \@features;
}



=head2 fetch_all_by_display_label

  Arg [1]    : String $label - display label of feature to fetch
  Example    : my $feat = $adaptor->fetch_all_by_display_label('abd-A_dpp:REDFLY:TF000092');
  Description: Returns the features with the given display label.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Feature objects
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

#There is no result_feature.display_label attribute.
#Move this to individual feature adaptors to avoid over-riding?


sub fetch_all_by_display_label {
  my ($self, $label) = @_;

  throw('You must defined a display_label argument') if ! defined $label;

  my $table_syn = $self->_main_table->[1];
  my $constraint = "${table_syn}.display_label = ?";
  $self->bind_param_generic_fetch($label, SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}

=head2 _count_features_by_field_id

  Arg [1]    : string     - table field to count
  Arg [2]    : string/int - id to count
  Example    : my $probe_feature_count = $pfa->count_feature_by_field_id('probe_id', $probe_id);
  Description: Returns a count of the features with the acciated field id
  Returntype : string/int - count of features
  Exceptions : Throws args are not set
  Caller     : FeatureAdaptors
  Status     : At risk

=cut

#This does not assume one record per feature
#But does assume primary key if ${table_name}_id
#Can't move to core due to cs issues, but could mirror implementation.

sub _count_features_by_field_id {
  my ($self, $count_field, $count_id) = @_;
  #Any other params here?

  if (! ($count_field && $count_id)) {
    throw('You must provide a count name and a count id to count by');
  }

  my ($table_name, $table_syn) = @{$self->_main_table};
  my $table_id   = "${table_name}_id";
  my @cs_ids     = @{$self->_get_coord_system_ids};
  my $sql = "SELECT count(distinct($table_id)) from $table_name $table_syn, seq_region sr where ${table_syn}.${count_field}=? and ${table_syn}.seq_region_id in(select distinct(seq_region_id) from seq_region where coord_system_id in(".join(',', @cs_ids).'))';
  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $count_id, SQL_INTEGER);
  $sth->execute;

  return $sth->fetchrow_array;
}



=head2 force_reslice

  Arg [1]    : Optional - Boolean
  Example    : if($self->force_reslice){
                    # Reslice features past ends of destination Slice
               }
  Description: Sets/Returns force_reslice boolean
  Returntype : Boolean
  Exceptions : None
  Caller     : FeatureAdaptors::_objs_from_sth
  Status     : At risk

=cut

sub force_reslice{
  my $self  = shift;
  my $force = shift;

  if (defined $force) {
    $self->{force_reslice} = $force;
  }

  return $self->{force_reslice};
}


# Currently over-riding this method in SetFeatureAdaptor
# But might be nice to remove that and genericise this by getting 
# the analysis_id table via a method which can be over-ridden e.g. 
# SetFeatureAdaptor::_analysis_id_table_syn which would return the relevant set table name and syn
# BaseFeatureAdaptor::_analysis_id_table_syn would simply return the _main_table syn
# 

#Should this use the cache via AnalysisAdpator::fetch_by_logic_name?
#As _logic_name_to_constriant does?
#_logic_name_to_constraint returns undef if logic_name is not known
# allowing caller to return [] early.
# How would this work in compose_constraint?
# simply delete this first if exists
# How does this work wrt redundancy between constraints and optional_constraints

# For generic methods, there is a slight redundancy in calling
# the _main_table method
# Can we call this in compose_constraint and pass to contrain method

sub _constrain_logic_names {
  my $self        = shift;
  my $logic_names = shift;
  assert_ref($logic_names, 'ARRAY');
  
  if(! @$logic_names){
    throw('Must pass an Arrayref of logic_name strings to contrain by');  
  }

  for my $lname(@$logic_names){
    if(! defined $lname){
      throw('Found undefined logic_name value');  
    }  
  }

  my $syn = $self->_main_table->[1];
  $self->_tables([['analysis', 'a']]);

  my $constraint = $syn.'.analysis_id = a.analysis_id and a.logic_name in ("'.
    join('", "', @$logic_names).'")';
    
  return ($constraint, {});   #{} = not further constraint conf
}

1;


