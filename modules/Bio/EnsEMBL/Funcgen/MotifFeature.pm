#
# Ensembl module for Bio::EnsEMBL::Funcgen::MotifFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::MotifFeature - A module to represent a feature mapping as based
on a binding matrix e.g position weight matrix

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MotifFeature;

my $feature = Bio::EnsEMBL::Funcgen::MotifFeature->new(
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
    -DISPLAY_LABEL => $text,
    -SCORE         => $score,
    -FEATURE_TYPE  => $ftype,
); 



=head1 DESCRIPTION

A MotifFeature object represents the genomic placement of a sequence motif.
For example a transcription factor binding site motif associated with a 
position weight matrix. These are generally associated with AnnotatedFeatures
of the corresponding FeatureType.


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor


=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::MotifFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Funcgen::Storable;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable);


=head2 new

 
  Arg [-SCORE]          : (optional) int - Score assigned by analysis pipeline
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
									  -DISPLAY_LABEL  => $text,
									  -SCORE          => $score,
                                                                         );


  Description: Constructor for MotifFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Exceptions : Throws if BindingMatrix not valid
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  my ($score, $bmatrix, $dlabel) = rearrange(['SCORE', 'BINDING_MATRIX', 'DISPLAY_LABEL'], @_);
    

  if(! (ref($bmatrix) && $bmatrix->isa('Bio::EnsEMBL::Funcgen::BindingMatrix'))){
	throw('You must pass be a valid Bio::EnsEMBL::Funcgen::BindingMatrix');
  }

  $self->{'binding_matrix'} = $bmatrix;
  $self->{'score'}          = $score  if $score;
  $self->display_label($dlabel) if $dlabel;

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

sub new_fast {
  return bless ($_[1], $_[0]);
}



=head2 binding_matrix

  Example    : my $bmatrix_name = $mfeat->binding_matrix->name;
  Description: Getter for the BindingMatrix attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub binding_matrix{
  my $self = shift;

  return $self->{'binding_matrix'};
}

=head2 score

  Arg [1]    : (optional) int - score
  Example    : my $score = $feature->score();
  Description: Getter/Setter for the score attribute for this feature. 
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub score {
    my $self = shift;
	
    $self->{'score'} = shift if @_;
		
    return $self->{'score'};
}


=head2 display_label

  Arg [1]    : string - display label
  Example    : my $label = $feature->display_label();
  Description: Getter/Setter for the display label of this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label {
  my ($self, $dlabel) = @_;
	
  #auto generate here if not set in table
  if ($dlabel){
	$self->{'display_label'} = $dlabel;
  }
  elsif(! defined $self->{'display_label'}){
	$self->{'display_label'} = $self->binding_matrix->feature_type->name.':'
	  .$self->binding_matrix->name();
  }
  
  return $self->{'display_label'};
}



=head2 associated_annotated_features

  Example    : my @associated_afs = @{$feature->associated_annotated_features()};
  Description: Getter/Setter for associated AnntoatedFeatures.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen:AnnotatedFeature objects
  Exceptions : None
  Caller     : General
  Status     : At risk - may change to associated_transcript_factor_features

=cut

sub associated_annotated_features{
  my ($self, $afs) = @_;
  
  #Lazy load as we don't want to have to do a join on all features when most will not have any

 
  if(defined $afs){

	if(ref($afs) eq 'ARRAY'){

	  foreach my $af(@$afs){
	
		if( ! $af->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature') ){
		  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::AnnotatedFeature objects');
		}
		#test is stored in adaptor
	  }

	  if(defined $self->{'associated_annotated_features'}){
		warn('You are overwriting associated_annotated_features for the MotifFeature');
		#we could simply add the new ones and make them NR.
	  }

	  $self->{'associated_annotated_features'} = $afs;
	}
	else{
	  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::AnnotatedFeature objects');
	}
  }


  if(! defined $self->{'associated_annotated_features'}){

	if(defined $self->adaptor){
	  $self->{'associated_annotated_features'} = 
		$self->adaptor->db->get_AnnotatedFeatureAdaptor->fetch_all_by_associated_MotifFeature($self);
	}
  }
  
  #This has the potential to return undef, or an arrayref which may be empty.
  return $self->{'associated_annotated_features'};
}


=head2 infer_variation_consequence

  Arg [1]    : Bio::EnsEMBL::Variation::VariationFeature
  Arg [2]    : boolean - 1 if results in linear scale (default is log scale)
  Example    : my $vfs = $vf_adaptor->fetch_all_by_Slice($slice_adaptor->fetch_by_region('toplevel',$mf->seq_region_name,$mf->start,$mf->end,$mf->strand));
               foreach my $vf (@{$vfs}){
                   print $mf->infer_variation_consequence($vf)."\n";
               }

  Description: Calculates the potential influence of a given variation in a motif feature.
               Returns a value between -100% (lost) and +100% (gain) indicating the difference 
               in strength between the motif in the reference and after the variation.

               The variation feature slice needs to be the motif feature, including the strand
  Returntype : float
  Exceptions : throws if argument is not a  Bio::EnsEMBL::Variation::VariationFeature
               throws if the variation feature is not contained in the motif feature
  Caller     : General
  Status     : At risk

=cut

sub infer_variation_consequence{
  my ($self, $vf, $linear) = (shift, shift, shift);

  if(! $vf->isa('Bio::EnsEMBL::Variation::VariationFeature')){
    throw "We expect a Bio::EnsEMBL::Variation::VariationFeature object, not a ".$vf->class;
  }

  #See if these checks are required or if there are more efficient ways to do the checks...
  #if(($self->slice->seq_region_name ne $vf->slice->seq_region_name) ||
  #   ($self->slice->start != $vf->slice->start) || 
  #   ($self->slice->end != $vf->slice->end) ){
  #  throw "Variation and Motif are on distinct slices";
  #}
  #if(!(($vf->start >= $self->start) && ($vf->end <= $self->end ))){
  #  throw "Variation should be entirely contained in the Motif";
  #}

  if( ($vf->start < 1) || ($vf->end > $self->binding_matrix->length)){ throw "Variation not entirely contained in the motif feature"; }

  if(!($vf->allele_string =~ /^[ACTG]\/[ACTG]$/)){ throw "Currently only SNPs are supported"; }

  my $ref_seq = $self->seq;

  my $variant = $vf->allele_string;
  $variant =~ s/^.*\///;
  $variant =~ s/\s*$//;

  my ($vf_start,$vf_end) = ($vf->start, $vf->end);
  if($vf->strand == -1){
    #Needed for insertions
    $variant = reverse($variant);
    $variant =~ tr/ACGT/TGCA/;
  }
  my $var_seq = substr($ref_seq,0, $vf_start - 1).$variant.substr($ref_seq, $vf_start+length($variant)-1);

  my $bm = $self->{'binding_matrix'};
  return 100 * ($bm->relative_affinity($var_seq,$linear) - $bm->relative_affinity($ref_seq,$linear));
  
}

1;

