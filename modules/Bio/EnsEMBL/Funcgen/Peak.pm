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

use Bio::EnsEMBL::Utils::Exception qw( throw );

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

=head2 score

  Example    : my $score = $feature->score;
  Description: Getter for the score attribute for this feature. 
  Returntype : String (float)
  Exceptions : None
  Caller     : General
  Status     : Stable

=head2 summit

  Arg [1]    : (optional) int - summit postition
  Example    : my $peak_summit = $feature->summit;
  Description: Getter for the summit attribute for this feature. 
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub dbID              { return shift->_generic_get_or_set('dbID',              @_);}
sub db                { return shift->_generic_get_or_set('db',                @_);}
sub adaptor           { return shift->_generic_get_or_set('db',                @_);}
sub peak_calling_id   { return shift->_generic_get_or_set('peak_calling_id',   @_);}
sub summit            { return shift->_generic_get_or_set('summit',            @_);}
sub score             { return shift->_generic_get_or_set('score',             @_);}
sub start             { return shift->_generic_get_or_set('start',             @_);}
sub end               { return shift->_generic_get_or_set('end',               @_);}
sub seq_region_id     { return shift->_generic_get_or_set('seq_region_id',     @_);}
sub seq_region_start  { return shift->_generic_get_or_set('seq_region_start',  @_);}
sub seq_region_end    { return shift->_generic_get_or_set('seq_region_end',    @_);}
sub seq_region_strand { return shift->_generic_get_or_set('seq_region_strand', @_);}
sub strand            { return shift->_generic_get_or_set('strand',            @_);}
sub slice             { return shift->_generic_get_or_set('slice',             @_);}

sub fetch_PeakCalling {
  return shift->_generic_fetch('peak_calling', 'get_PeakCallingAdaptor', 'peak_calling_id');
}

sub set_PeakCalling {
  my $self = shift;
  my $obj  = shift;
  return shift->_generic_set('peak_calling', 'Bio::EnsEMBL::Funcgen::PeakCalling', $obj);
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
      $self->{'display_label'} = $self->feature_type->name()." -";
      $self->{'display_label'} .= " ".$self->epigenome->display_label();
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

    if(! $self->{'display_id'}  && $self->adaptor){
      $self->{'display_id'} = join '_', 
        $self->feature_type->name(),
        $self->epigenome->production_name(),
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
  Status     : At Risk - This is TFBS specific and could move to TranscriptionFactorFeature

=cut

sub get_underlying_structure{
  return [];
}

=head2 get_associated_MotifFeatures

  Example    : my @assoc_mfs = @{ $af->get_associated_MotifFeatures };
  Description: There are none in the database, so this always returns undef.
  Returntype : undef
  Exceptions : None
  Caller     : General
  Status     : stable

=cut

sub get_associated_MotifFeatures{
  return [];
}

sub SO_term {
  my $self = shift;
  return $self->feature_type->so_accession;
}

=head2 summary_as_hash

  Example       : $segf_summary = $annotf->summary_as_hash;
  Description   : Retrieves a textual summary of this Peak.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;
  my $feature_set = $self->feature_set;

  return
    {
      feature_type     => $self->feature_type->name,
      epigenome        => $self->epigenome->name,
      source           => $feature_set->analysis->logic_name,
      seq_region_name  => $self->seq_region_name,
      start            => $self->seq_region_start,
      end              => $self->seq_region_end,
      description      => $feature_set->display_label,
      strand           => $self->strand,
      summit           => $self->summit,
      score            => $self->score,
    };
}
1;

