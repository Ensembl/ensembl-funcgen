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

Bio::EnsEMBL::Peak - An object to represent a peak.

=head1 SYNOPSIS

  use strict;
  use warnings;
  use Bio::EnsEMBL::Registry;
  use List::Util qw( min );

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org',
      -user => 'anonymous'
  );

  my $peak_adaptor  = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'Peak');
  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'core',    'Slice');

  # Fetch a slice
  my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '17', 63_992_802, 64_038_237);

  # Fetch all peaks on the slice
  my $peaks_on_slice = $peak_adaptor->fetch_all_by_Slice($slice);

  my $number_of_peaks_on_slice = @$peaks_on_slice;
  print "There are $number_of_peaks_on_slice peaks on the slice.\n";

  # Print the first ten
  my $max_features_to_print = 10;

  # This prints:
  # 
  # There are 697 peaks on the slice.
  # H3K36me3 - Fetal Stomach Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K27me3 - EM CD8+ ab T cell (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K9me3 - CD14+CD16- monocyte (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K36me3 - H1-trophoblast Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K36me3 - naive B cell (To) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K36me3 - CD14+CD16- monocyte (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K36me3 - naive B cell (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K27me3 - neutrophil (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K27me3 - M2 macrophage (VB) Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  # H3K36me3 - Fetal Stomach Enriched Site (chromosome:GRCh37:17:63992802:64038237:1)
  #
  for my $i ( 1.. min($max_features_to_print, $number_of_peaks_on_slice) ) {
    print_feature($peaks_on_slice->[$i]);
  }

  sub print_feature {
    my $feature = shift;
    print 
      $feature->display_label
      . " (" 
      . $feature->slice->name 
      . ")\n";
  }


=head1 DESCRIPTION

Peaks are the enriched regions that were found during a ChIP-seq or a similar
high throughput experiment.

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
sub slice             { return shift->_generic_get_or_set('slice',             @_);}

=head2 fetch_PeakCalling

  Example    : 
  Description: Fetches the peak calling object representing the peak 
               calling that generated this peak.
  Returntype : Bio::EnsEMBL::Funcgen::PeakCalling
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub fetch_PeakCalling {
  return shift->_generic_fetch('peak_calling', 'get_PeakCallingAdaptor', 'peak_calling_id');
}

=head2 fetch_FeatureType

  Example    : my $feature_type = $peak->fetch_FeatureType;
  Description: Fetches the type of feature that this peak represent on the 
               sequence. E.g.: 
                  - Cfos,
                  - Cjun,
                  - DNase1,
                  - Gabp,
                  - H3K23me2 or
                  - Rad21.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub fetch_FeatureType {

  my $self = shift;
  my $peak_calling = $self->fetch_PeakCalling;
  my $feature_type = $peak_calling->fetch_FeatureType;
  return $feature_type;
}

=head2 fetch_Epigenome

  Example    : my $epigenome = $peak->fetch_Epigenome;
  Description: Fetches the epigenome on which the high throughput experiment 
               was performed that yielded this peak.
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub fetch_Epigenome {

  my $self = shift;
  my $peak_calling = $self->fetch_PeakCalling;
  my $epigenome = $peak_calling->fetch_Epigenome;
  return $epigenome;
}

=head2 set_PeakCalling

  Args       : Object of type Bio::EnsEMBL::Funcgen::PeakCalling
  Example    : $peak->set_PeakCalling($peak_calling);
  Description: Setter for the peak calling object representing the peak 
               calling that generated this peak.
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

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

    if(! $self->{'display_id'}  && $self->adaptor){
      $self->{'display_id'} = join '_', 
        $self->fetch_FeatureType->name(),
        $self->fetch_Epigenome->production_name(),
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
  return $self->fetch_FeatureType->so_accession;
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
      feature_type     => $self->fetch_FeatureType->name,
      epigenome        => $self->fetch_Epigenome->name,
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

