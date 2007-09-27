#
# Ensembl module for Bio::EnsEMBL::Funcgen::ExternalFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::ExternalFeature - A module to represent an externally curated feature 
mapping from an external_db.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ExternalFeature;

my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
    -DISPLAY_LABEL => $text,
    -FEATURE_SET   => $fset,
    -FEATURE_TYPE  => $ftype,
);



=head1 DESCRIPTION

An ExternalFeature object represents the genomic placement of an externally curated
feature from and DB external to Ensembl.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ExternalFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::SetFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::SetFeature);


=head2 new

 
  Arg [-FEATURE_SET]  : Bio::EnsEMBL::Funcgen::FeatureSet
  Arg [-FEATURE_TYPE] : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-START]        : int - The start coordinate of this feature relative to the start of the slice
		                it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int -The end coordinate of this feature relative to the start of the slice
	                    it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-DISPLAY_LABEL]: string - Display label for this feature
  Arg [-STRAND]       : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int - Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.
  Example             : my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
										                                  -SLICE         => $chr_1_slice,
									                                      -START         => 1_000_000,
                      				                                      -END           => 1_000_024,
									                                      -STRAND        => -1,
									                                      -DISPLAY_LABEL => $text,
									                                      -FEATURE_SET   => $fset,
                                                                          -FEATURE_TYPE  => $ftpe,
                                                                         );


  Description: Constructor for ExternalFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::ExternalFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  my ($ftype) = rearrange(['FEATURE_TYPE'], @_);
  
  $self->feature_type($ftype) if $ftype;
		
  return $self;
}


=head2 display_label

  Arg [1]    : string - display label
  Example    : my $label = $feature->display_label();
  Description: Getter and setter for the display label of this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label {
  my $self = shift;
	
  $self->{'display_label'} = shift if @_;
    
  if(! $self->{'display_label'}  && $self->adaptor()){
	
	$self->{'display_label'} = $self->feature_set->feature_type->name().' - ';
	$self->{'display_label'} .= $self->cell_type->name() if $self->cell_type();
	$self->{'display_label'} .= $self->feature_type->name() if(defined $self->{'feature_type'});
  }
	
  return $self->{'display_label'};
}


=head2 feature_type

  Example    : my $ft_name = $efeature->feature_type()->name();
  Description: Getter/Setter for the feature_type attribute for this feature.
               Defaults to FeatureSet FeatureType if absent?
  Returntype : Bio::EnsEMBL::Funcgen:FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub feature_type{
  my ($self, $ftype) = @_;
  
  if(defined $ftype && ! (ref($ftype) && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType"))){
	throw("Must pass valid Bio::EnsEMBL::Funcgen::FeatureType object");
  }else{
	$self->{'feature_type'} = $ftype;
  }

  return (defined $self->{'feature_type'}) ? $self->{'feature_type'} : $self->feature_set->feature_type();
}

#=head2 external_db_name

#  Example    : my $edb_name = $efeature->external_db_name();
#  Description: Getter for the external_db_name  attribute for this feature.
#  Returntype : string - external_db_name
#  Exceptions : None
#  Caller     : General
#  Status     : At risk

#=cut

#sub external_db_name{
#  my $self = shift;
  
#  return $self->feature_set->external_db_name();
#}

1;

