#
# Ensembl module for Bio::EnsEMBL::Funcgen::MotifFeature
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

Bio::EnsEMBL::MotifFeature - A module to represent a feature mapping as based
on a binding matrix e.g position weight matrix

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MotifFeature;

my $feature = Bio::EnsEMBL::Funcgen::MotifFeature->new
 (
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
    -DISPLAY_LABEL => $text,
    -SCORE         => $score,
    -FEATURE_TYPE  => $ftype,
    -STABLE_ID     => 1,
 );

=head1 DESCRIPTION

A MotifFeature object represents the genomic placement of a sequence motif.
For example a transcription factor binding site motif associated with a
position weight matrix. These are generally associated with AnnotatedFeatures
of the corresponding FeatureType.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Funcgen::MotifFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

use base qw(Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable);

use constant SEQUENCE_ONTOLOGY => {
  acc  => 'SO:0000235',
  term => 'TF_binding_site',
};


=head2 new

  Arg [-SCORE]          : (optional) int - Score given by the motif mapper.
  Arg [-SLICE]          : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-BINDING_MATRIX] : Bio::EnsEMBL::Funcgen::BindingMatrix - Binding Matrix associated to this feature.
  Arg [-START]          : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]            : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]  : string - Display label for this feature
  Arg [-STRAND]         : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]           : (optional) int - Internal database ID.
  Arg [-ADAPTOR]        : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::MotifFeature->new(
                                									  -SLICE          => $chr_1_slice,
								                                	  -START          => 1_000_000,
                                									  -END            => 1_000_024,
								                                	  -STRAND         => -1,
                                									  -BINDING_MATRIX => $bm,
                                  									-SCORE          => $score,
                                                    -STABLE_ID     => 1 );

  Description: Constructor for MotifFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Exceptions : Throws if BindingMatrix not valid
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self   = $class->SUPER::new(@_);

    my ($score, $binding_matrix, $stable_id) =
        rearrange([ 'SCORE', 'BINDING_MATRIX', 'STABLE_ID' ], @_);

    throw('Must supply a -score parameter') if !defined $score;
    throw('Must supply a -binding_matrix parameter') if !defined $binding_matrix;

    assert_ref($binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix');

    $self->{score}                         = $score;
    $self->{binding_matrix}                = $binding_matrix;
    $self->{stable_id}                     = $stable_id if $stable_id;
    $self->{overlapping_Peaks}             = undef;
    $self->{overlapping_RegulatoryFeature} = undef;

    return $self;
}


=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the calling code
               is very disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast { return bless ($_[1], $_[0]); }

=head2 get_BindingMatrix

  Example    : my $binding_matrix = $bmf->get_BindingMatrix();
  Description: Getter for the BindingMatrix object
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_BindingMatrix { return shift->{binding_matrix}; }

=head2 score

  Example    : my $score = $feature->score();
  Description: Getter for the score attribute for this feature. This is the
               score given by the motif mapping software.
  Returntype : Scalar (double)
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub score { return shift->{score}; }

=head2 get_all_overlapping_Peaks

  Example    : my $peaks = $motif_feature->get_all_overlapping_Peaks;
  Description: Gets all Peaks that overlap with this motif feature
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Peak objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub get_all_overlapping_Peaks {
    my $self = shift;

    if ( !$self->{overlapping_Peaks} ) {
        $self->{overlapping_Peaks} =
            $self->adaptor()->_fetch_all_overlapping_Peaks($self);
    }

    return $self->{overlapping_Peaks};
}

=head2 get_all_overlapping_Peaks_by_Epigenome

  Arg [1]    : Bio::EnsEMBL::Funcgen::Epigenome object
  Example    : my $peak =
             :    $motif_feature->fetch_overlapping_Peak_by_Epigenome($epigenome);
  Description: Fetches all overlapping Peaks for a particular Epigenome
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Peak objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub get_all_overlapping_Peaks_by_Epigenome {
    my ($self, $epigenome) = @_;

    my $peaks =
      $self->adaptor->_fetch_all_overlapping_Peaks_by_Epigenome($self,
                                                                $epigenome);

    return $peaks;
}

=head2 is_experimentally_verified_in_Epigenome

  Arg [1]    : Bio::EnsEMBL::Funcgen::Epigenome object
  Example    : my $is_verified =
             :    $motif_feature->is_experimentally_verified_in_Epigenome($epigenome);
  Description: Returns true if the motif feature is experimentally verified
             : (i.e. has an overlapping Peak) in the given epigenome.
  Returntype : Boolean
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub is_experimentally_verified_in_Epigenome {
    my ($self, $epigenome) = @_;

    my $is_experimentally_verified = 0;

    my $peaks = $self->get_all_overlapping_Peaks_by_Epigenome($epigenome);

    if (scalar @{$peaks} > 0){
        $is_experimentally_verified = 1;
    }

    return $is_experimentally_verified;
}

=head2 get_all_Epigenomes_with_experimental_evidence

  Example    : my $epigenomes =
             :    $motif_feature->get_all_Epigenomes_with_experimental_evidence;
  Description: Returns a list of Epigenomes where the motif feature has been
               experimentally verified
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::Epigenome objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub get_all_Epigenomes_with_experimental_evidence {
  my $self = shift;
  my $peaks = $self->get_all_overlapping_Peaks;
  my %epigenome_dbIDs;
  for my $peak (@{$peaks}){
    $epigenome_dbIDs{$peak->get_PeakCalling->epigenome_id} = 1;
  }

  my $epigenome_adaptor = $self->adaptor->db->get_adaptor('Epigenome');
  my @dbID_list = keys %epigenome_dbIDs;
  my $epigenomes = $epigenome_adaptor->fetch_all_by_dbID_list(\@dbID_list);
  return $epigenomes;
}

=head2 is_position_informative

  Arg [1]    : Scalar - 1-based integer position within the motif
  Example    : $mf->is_position_informative($pos);
  Description: Indicates if a given position within the motif is highly informative
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : At High risk

=cut

sub is_position_informative {
    my $self     = shift;
    my $position = shift;

    return $self->get_BindingMatrix->is_position_informative($position);
}


=head2 infer_variation_consequence

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature
  Arg [2]    : Boolean - returns result in linear scale (default is log scale)
  Example    : my $vfs = $vf_adaptor->fetch_all_by_Slice($slice_adaptor->fetch_by_region('toplevel',$mf->seq_region_name,$mf->start,$mf->end,$mf->strand));
               foreach my $vf (@{$vfs}){
                   print $mf->infer_variation_consequence($vf)."\n";
               }

  Description: Calculates the potential influence of a given variation in a motif feature.
               Returns a value between -100% (lost) and +100% (gain) indicating the difference
               in strength between the motif in the reference and after the variation.

  Returntype : Scalar (numeric) or undef
  Exceptions : Throws if argument is not a Bio::EnsEMBL::Variation::VariationFeature
               Warns if the VariationFeature is not contained within the MotifFeature
  Caller     : General
  Status     : At High risk

=cut

sub infer_variation_consequence{
  my $self   = shift;
  my $vf     = shift;
  my $linear = shift;
  assert_ref($vf, 'Bio::EnsEMBL::Variation::VariationFeature');
  my $vf_sr_start = $vf->seq_region_start;
  my $sr_start    = $self->seq_region_start;
  my $allele      = $vf->allele_string(undef, $self->seq_region_strand);

  if($allele !~ /^[ACTG]\/[ACTG]$/){
    throw("Unsupported variation allele:\t".$allele."\nCurrently only SNPs supported");
  }

  # From now on, assumes variation is a SNP
  if( ! (($self->seq_region_name eq $vf->seq_region_name) &&
         ($sr_start <= $vf_sr_start) &&
         ($self->seq_region_end >= $vf_sr_start))){

    warn('VariationFeature('.$vf->variation_name.
      ") is not contained within MotifFeature:\t".$self->slice->name."\n");
    return 0; # 0 as we need a consequence delta
  }

  my $ref_seq = $self->seq; # Get the stranded seq
  $allele =~ s/^.*\/\s*//;  # Get the strand specific non-ref allele

  my $vf_idx = $vf_sr_start - $sr_start; # 0 based to avoid unecessary -1 in substr

  my $var_seq = substr($ref_seq, 0, $vf_idx).$allele.
    substr($ref_seq, $vf_idx + 1);  # + length($variant));

  # relative affinity only works with strand matched seq
  # in 5'->3' orientation. We already have the strand seq
  # so just need to reverse if -1

  if($self->seq_region_strand == -1){
    $var_seq = reverse($var_seq);#tr/ACGT/TGCA/;
    $ref_seq = reverse($ref_seq);
  }

  my $bm     = $self->get_BindingMatrix;
  my $var_ra = $bm->relative_sequence_similarity_score($var_seq, $linear);
  my $ref_ra = $bm->relative_sequence_similarity_score($ref_seq, $linear);

  return (defined $var_ra && defined $ref_ra ) ? (100 * ($var_ra - $ref_ra)) : undef;
}


=head2 stable_id

  Example    : my $stable_id = $feature->stable_id();
  Description: Getter for the stable_id attribute for this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub stable_id { return shift->{stable_id}; }

=head2 summary_as_hash
  Example       : $motif_feature_summary = $motif_feature->summary_as_hash;
  Description   : Retrieves a textual summary of this MotifFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)
=cut

sub summary_as_hash {
    my $self = shift;

    my $summary = {
        binding_matrix_stable_id     => $self->get_BindingMatrix->stable_id,
        start                        => $self->seq_region_start,
        end                          => $self->seq_region_end,
        strand                       => $self->strand,
        seq_region_name              => $self->seq_region_name,
        stable_id                    => $self->stable_id,
        score                        => $self->score,
        transcription_factor_complex => join( ',',
            @{ $self->get_BindingMatrix->get_TranscriptionFactorComplex_names } ),
    };

    my $epigenomes = $self->get_all_Epigenomes_with_experimental_evidence;

    if (scalar @{$epigenomes} > 0) {
        my @epigenome_names_list;
        for my $epigenome ( @{$epigenomes} ) {
            push @epigenome_names_list, $epigenome->short_name();
        }
        $summary->{epigenomes_with_experimental_evidence} = join ',',
          @epigenome_names_list;
    }

    return $summary;
}

1;
