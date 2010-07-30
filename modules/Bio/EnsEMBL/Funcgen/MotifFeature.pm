#
# Ensembl module for Bio::EnsEMBL::Funcgen::MotifFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::MotifFeature - A module to represent a feature mapping as based
on a binding matrix e.g position weight matrix

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MotifFeature;

my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
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

 
  Arg [-SCORE]        : (optional) int - Score assigned by analysis pipeline
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-START]        : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]: string - Display label for this feature
  Arg [-STRAND]       : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-FEATURE_TYPE] : Bio::EnsEMBL::Funcgen::FeatureType - The feature type corresponding to the
                        gene associated with thie given motif.
  Arg [-dbID]         : (optional) int - Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.

  Example    : my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
										                                  -SLICE         => $chr_1_slice,
									                                      -START         => 1_000_000,
									                                      -END           => 1_000_024,
									                                      -STRAND        => -1,
									                                      -DISPLAY_LABEL => $text,
									                                      -SCORE         => $score,
                                                                          -FEATURE_TYPE  => $ftype,
                                                                         );


  Description: Constructor for AnnotatedFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::AnnotatedFeature
  Exceptions : Throws if FeatureType not valid
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
	$self->{'display_label'} = $self->binding_matrix->name()." motif";
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



1;

