#
# Ensembl module for Bio::EnsEMBL::Funcgen::OligoFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::OligoFeature - A module to represent an oligonucleotide probe's
genomic mapping.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::OligoFeature;

my $feature = Bio::EnsEMBL::Funcgen::OligoFeature->new(
	-PROBE         => $probe,
	-MISMATCHCOUNT => 0,
	-SLICE         => $chr_1_slice,
	-START         => 1_000_000,
	-END           => 1_000_024,
	-STRAND        => -1,
); 

#build_id/version, seq_region_id, seq_region_strand -  from slice?
#cigar_line?


=head1 DESCRIPTION

An OligoFeature object represents the genomic placement of an OligoProbe
object. The data are stored in the oligo_feature table.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
AffyFeature module written by Ian Sealy and Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::OligoFeature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Feature;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [-PROBE]        : Bio::EnsEMBL::Funcgen::OligoProbe - probe
        An OligoFeature must have a probe. This probe must already be stored if
		you plan to store the feature.
  Arg [-MISMATCHCOUNT]: int
        Number of mismatches over the length of the probe. 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice
        The slice on which this feature is.
  Arg [-START]        : int
        The start coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int
        The end coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-STRAND]       : int
        The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int
        Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
        Database adaptor.
  Example    : my $feature = Bio::EnsEMBL::Funcgen::OligoFeature->new(
				   -PROBE         => $probe,
				   -MISMATCHCOUNT => 0,
				   -SLICE         => $chr_1_slice,
				   -START         => 1_000_000,
				   -END           => 1_000_024,
				   -STRAND        => -1,
#-ANALYSIS?
			   ); 
  Description: Constructor for OligoFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::OligoFeature
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;
	
	my $class = ref($caller) || $caller;
	
	my $self = $class->SUPER::new(@_);
	
	my ($probe, $mismatchcount, $coord_sys_id )
		= rearrange(['PROBE', 'MISMATCHCOUNT', 'COORD_SYSTEM_ID'], @_);
	
	$self->probe($probe);
	$self->mismatchcount($mismatchcount);

	#do we need to validate this against the db?  Grab from slice and create new if not present?  Will this be from the dnadb?
	
	$self->coord_system_id($coord_sys_id);
	
	return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::OligoFeature
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new_fast {
   my ($class, $hashref)  = @_;


   return bless ($hashref, $class);
}

=head2 probeset

  Arg [1]    : (optional) string - probeset
  Example    : my $probeset = $feature->probeset();
  Description: Getter and setter for the probeset for this feature. Shortcut
               for $feature->probe->probeset(), which should be used instead.
			   Probeset is not persisted if set with this method.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->probeset() because this may be removed

=cut

sub probeset {
    my $self = shift;
	
    $self->{'probeset'} = shift if @_;
	
    if ($self->{'probe'}) {
		$self->{'probeset'} = $self->probe()->probeset();
    }
	
    return $self->{'probeset'};
}

=head2 mismatchcount

  Arg [1]    : int - number of mismatches
  Example    : my $mismatches = $feature->mismatchcount();
  Description: Getter and setter for number of mismatches for this feature.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub mismatchcount {
    my $self = shift;
	
    $self->{'mismatchcount'} = shift if @_;
	
    return $self->{'mismatchcount'};
}


=head2 coord_system_id

  Arg [1]    : int - dbID of corresponding coord_system for DB of origin
  Example    : $feature->coord_system_id($cs_id);
  Description: Getter and setter for the coord system id for this feature.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub coord_system_id {
    my $self = shift;
	
    $self->{'coord_system_id'} = shift if @_;
	
    return $self->{'coord_system_id'};
}



=head2 probelength

  Args       : None 
  Example    : my $probelength = $feature->probelength();
  Description: Getter for the length of the probe. Shortcut for
               $feature->probe->length(), which should be used instead.
			   Originally, this method returned the length of the feature,
			   which was often, but not always, the same as the length of the
			   probe.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->length() because this may be removed

=cut

sub probelength {
    my $self = shift;
	
    return $self->probe->length();
}

=head2 probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoProbe - probe
  Example    : my $probe = $feature->probe();
  Description: Getter, setter and lazy loader of probe attribute for
               OligoFeature objects. Features are retrieved from the database
			   without attached probes, so retrieving probe information for a
			   feature will involve another query.
  Returntype : Bio::EnsEMBL::Funcgen::OligoProbe
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probe {
    my $self = shift;
	my $probe = shift;
    if ($probe) {
		if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::Funcgen::OligoProbe') ) {
			throw('Probe must be a Bio::EnsEMBL::Funcgen::OligoProbe object');
		}
		$self->{'probe'} = $probe;
    }
	if ( !defined $self->{'probe'} && $self->dbID() && $self->adaptor() ) {
	    $self->{'probe'} = $self->adaptor()->db()->get_OligoProbeAdaptor()->fetch_by_OligoFeature($self);
	}
    return $self->{'probe'};
}

=head2 get_results_by_channel

  Arg [1]    : int - channel_id (mandatory)
  Arg [2]    : string - Analysis name e.g. RawValue, VSN (optional)
  Example    : my @results = $feature->results();
  Description: Getter, setter and lazy loader of results attribute for
               OligoFeature objects.
  Returntype : List ref to arrays containing ('score', 'Analysis logic_name');
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_results_by_channel {
    my $self = shift;
	my $channel_id = shift;
	my $anal_name = shift;

	#$self->{'results'} ||= {};
	$self->{'results_complete'} ||= 0;
	
	if(! $self->{'results'} || ($anal_name && ! exists $self->{'results'}{$anal_name})){
		#fetch all, set complete set flag
		$self->{'results_complete'} ||= 1 	if(! $anal_name);

		foreach my $results_ref(@{$self->adaptor->fetch_results_by_channel_analysis($self->probe->dbID(), 
																				  $channel_id, $anal_name)}){

			$self->{'results'}{$$results_ref[1]} = $$results_ref[0];
		}
	}


	
    return $self->{'results'}
}


#Will this be too slow, can we not do one query across all tables


1;

