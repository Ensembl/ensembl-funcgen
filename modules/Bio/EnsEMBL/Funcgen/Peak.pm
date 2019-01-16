=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Peak - A module to represent an enriched feature mapping i.e. a peak call.

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Peak;

use strict;
use warnings;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID              => 'dbID',
    db                => 'db',
    peak_calling_id   => 'peak_calling_id',
    summit            => 'summit',
    score             => 'score',
    slice             => 'slice',
    seq_region_id     => 'seq_region_id',
    seq_region_start  => 'seq_region_start',
    seq_region_end    => 'seq_region_end',
    seq_region_strand => 'seq_region_strand',
    peak_calling      => 'set_PeakCalling',
  };
}

use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

sub dbID              { return shift->_generic_get_or_set('dbID',              @_);}
sub db                { return shift->_generic_get_or_set('db',                @_);}
sub adaptor           { return shift->_generic_get_or_set('db',                @_);}
sub peak_calling_id   { return shift->_generic_get_or_set('peak_calling_id',   @_);}

=head2 summit

  Example    : 
  Description: Accessor for the summit attribute. This is the base at which 
               the peak had the highest level of enrichment.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub summit            { return shift->_generic_get_or_set('summit',            @_);}

=head2 score

  Example    : 
  Description: Accessor for the score. This is the score assigned by the peak 
               caller.
  Returntype : Float
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub score             { return shift->_generic_get_or_set('score',             @_);}

=head2 start

  Example    : 
  Description: Accessor for the start. This is the start of the peak on the 
               slice.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub start             { return shift->_generic_get_or_set('start',             @_);}

=head2 end

  Example    : 
  Description: Accessor for the end. This is the end of the peak on the slice.
  Returntype : Int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub end               { return shift->_generic_get_or_set('end',               @_);}
sub seq_region_id     { return shift->_generic_get_or_set('seq_region_id',     @_);}
sub seq_region_start  { return shift->_generic_get_or_set('seq_region_start',  @_);}
sub seq_region_end    { return shift->_generic_get_or_set('seq_region_end',    @_);}
sub seq_region_strand { return shift->_generic_get_or_set('seq_region_strand', @_);}
sub strand            { return shift->_generic_get_or_set('strand',            @_);}

=head2 slice

  Example    : 
  Description: Accessor for the slice attribute.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub slice { 
  my $self  = shift;
  my $slice = shift;
  
  if ($slice) {
    $self->seq_region_id($slice->get_seq_region_id);
  }
  return $self->_generic_get_or_set('slice', $slice);
}

=head2 fetch_PeakCalling

  Example    : 
  Description: Fetches the peak calling object representing the peak 
               calling that generated this peak.
  Returntype : Bio::EnsEMBL::Funcgen::PeakCalling
  Exceptions : None
  Caller     : general
  Status     : Deprecated

=cut

sub fetch_PeakCalling {
  my $deprecation_message = 'It will be removed in release 100. Please use '
      . 'Bio::EnsEMBL::Funcgen::Peak::get_PeakCalling instead.';
  deprecate($deprecation_message);
  return shift->_generic_fetch('peak_calling', 'get_PeakCallingAdaptor', 'peak_calling_id');
}

=head2 get_PeakCalling

  Example    :
  Description: Gets the peak calling object representing the peak
               calling that generated this peak.
  Returntype : Bio::EnsEMBL::Funcgen::PeakCalling
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub get_PeakCalling {
    return shift->_generic_fetch('peak_calling', 'get_PeakCallingAdaptor', 'peak_calling_id');
}

=head2 set_PeakCalling

  Args       : Object of type Bio::EnsEMBL::Funcgen::PeakCalling
  Example    : 
  Description: Setter for the peak calling object representing the peak 
               calling that generated this peak.
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub set_PeakCalling {
  my $self = shift;
  my $peak_calling = shift;
  
  $self->_generic_set('peak_calling', 'Bio::EnsEMBL::Funcgen::PeakCalling', $peak_calling);
  $self->peak_calling_id($peak_calling->dbID);
  return;
}

=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label {
    my $self = shift;

    #auto generate here if not set in table
    #need to go with one or other, or can we have both, split into diplay_name and display_label?
    
    if(! $self->{'display_label'}  && $self->adaptor){
      $self->{'display_label'} = $self->fetch_FeatureType->name()." -";
      $self->{'display_label'} .= " ".$self->fetch_Epigenome->display_label();
      $self->{'display_label'} .= " Enriched Site";
    }
	
    return $self->{'display_label'};
}

=head2 display_id

  Example    : my $label = $feature->display_id;
  Description: Getter for the display_id of this feature. This was created 
               for generating the display id used in big bed files. Converting
               from bed to bigbed causes problems, if 
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_id {
    my $self = shift;
    my $peak_calling = $self->fetch_PeakCalling;

    if(! $self->{'display_id'}  && $self->adaptor){
      $self->{'display_id'} = join '_', 
        $peak_calling->fetch_FeatureType->name(),
        $peak_calling->fetch_Epigenome->production_name(),
        "_Enriched_Site";
    }
    return $self->{'display_id'};
}

=head2 get_underlying_structure

  Example    : my @loci = @{ $af->get_underlying_structure() };
  Description: Returns and array of loci consisting of:
                  (start, (motif_feature_start, motif_feature_end)*, end)
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : Deprecated

=cut

sub get_underlying_structure{
  my $deprecation_message = 'It will be removed in release 100';
  deprecate($deprecation_message);
  return [];
}

=head2 fetch_all_MotifFeatures

  Example    : my $overlapping_motif_features = $peak->fetch_all_MotifFeatures
  Description: Returns all MotifFeatures that overlap with a Peak
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::MotifFeature objects
  Exceptions : none
  Caller     : General
  Status     : Deprecated

=cut

sub fetch_all_MotifFeatures {
    my $self = shift;
    my $deprecation_message = 'It will be removed in release 100. Please use '
        . 'Bio::EnsEMBL::Funcgen::Peak::get_all_MotifFeatures instead.';
    deprecate($deprecation_message);
    return $self->get_all_MotifFeatures();
}

=head2 get_all_MotifFeatures

  Example    : my $overlapping_motif_features = $peak->get_all_MotifFeatures
  Description: Returns all MotifFeatures that overlap with a Peak
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::MotifFeature objects
  Exceptions : none
  Caller     : General
  Status     : Stable

=cut

sub get_all_MotifFeatures {
    my $self = shift;

    my $motif_features
        = $self->adaptor()->_fetch_overlapping_MotifFeatures( $self );

    return $motif_features;
}

sub seq_region_name {
  my $self = shift;
  my $slice = $self->slice;
  return $slice->seq_region_name;
}

sub feature_so_acc {
  my $self = shift;
  return $self->fetch_PeakCalling->fetch_FeatureType->so_accession;
}

=head2 summary_as_hash

  Example       : $segf_summary = $annotf->summary_as_hash;
  Description   : Retrieves a textual summary of this Peak.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;
  my $peak_calling = $self->fetch_PeakCalling;
  my $slice = $self->slice;
  
  return
    {
      id               => $self->dbID,
      feature_type     => $peak_calling->fetch_FeatureType->name,
      epigenome        => $peak_calling->fetch_Epigenome->name,
      source           => $peak_calling->fetch_Analysis->logic_name,
      seq_region_name  => $slice->seq_region_name,
      start            => $slice->seq_region_start,
      end              => $slice->seq_region_end,
      description      => $peak_calling->display_label,
      strand           => $self->strand,
      summit           => $self->summit,
      score            => $self->score,
    };
}
1;

