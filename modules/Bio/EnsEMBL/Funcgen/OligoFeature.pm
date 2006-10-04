#
# Ensembl module for Bio::EnsEMBL::Funcgen::OligoFeature
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::OligoFeature - A module to represent an oligonucleotide probe
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
core OligoFeature module written by Ian Sealy.

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
				   -MISMATCHCOUNT => 0,#remove?
				   -SLICE         => $chr_1_slice,
				   -START         => 1_000_000,
				   -END           => 1_000_024,
				   -STRAND        => -1,
								      #-ANALYSIS?CIGARLINE!!!!!!!!!!
			   ); 
  Description: Constructor for OligoFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::OligoFeature
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($probe, $mismatchcount, $coord_sys_id )
    = rearrange(['PROBE', 'MISMATCHCOUNT', 'COORD_SYSTEM_ID'], @_);

  #Need to add analysis/cigar_line(remove mismatch?)

	
  $self->probe($probe);
  $self->mismatchcount($mismatchcount);
  
  #do we need to validate this against the db?  Grab from slice and create new if not present?  Will this be from the dnadb?
  
  #do we need this coordsys id if we're passing a slice?  We should have the method but not in here?

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
  Status     : High Risk

=cut

sub mismatchcount {
    my $self = shift;
	

    #replace with dynamic check of cigarline?

    $self->{'mismatchcount'} = shift if @_;
	
    return $self->{'mismatchcount'};
}


=head2 cigar_line

  Arg [1]    : str - Cigar line alignment annotation
  Example    : my $cg = $feature->cigar_line();
  Description: Getter and setter for number of the cigar line attribute for this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub cigar_line {
  my $self = shift;
  
  $self->{'cigar_line'} = shift if @_;
	
  return $self->{'cigar_line'};
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

=head2 get_results_by_channel_id

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

sub get_results_by_channel_id {
    my $self = shift;
    my $channel_id = shift;
    my $anal_name = shift;

    warn "This method not fully implemented, remove/deprecate?";

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


#The experiment/al chip specificity has already been done by the ofa->fetch_all_by_Slice_Experiment
#This may be called with no preceding Experiment specificity
#this would return results for all experiments
#do we need to set a default Experiment?


#THis should return both Chip and Channel based results
#just Chip for now
#maybe retrieve and hash all if not Analysis object passed?  Then return what?  


=head2 get_result_by_Analysis_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : listref - Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my $result = $feature->get_result_by_Analysis_ExperimentalChips($anal, \@echips);
  Description: Getter of results attribute for a given Analysis and set of ExperimentalChips
  Returntype : float
  Exceptions : Throws is no Analysis or ExperimentalChips are not passed?
  Caller     : General
  Status     : High Risk

=cut


#make ExperimentalChips optional?

#or have ResultSetAdaptor?  Do we need a ResultSet?
#may not have ExperimentalChip, so would need to return ec dbID aswell


######This will break/return anomalous if
#ECs are passed from different experiments
#ECs are passed from different Arrays


sub get_result_by_Analysis_ExperimentalChips{
    my ($self, $anal, $exp_chips) = @_;

    throw("Need to pass listref of ExperiemntalChips") if(scalar(@$exp_chips) == 0);
    throw("Need to pass a valid Bio::EnsEMBL::Analysis") if ! $anal->isa("Bio::EnsEMBL::Analysis");

    my (%query_ids, %all_ids);
    my $anal_name = $anal->logic_name();
    
    foreach my $ec(@$exp_chips){
				
      throw("Need to pass a listref of Bio::EnsEMBL::Funcgen::ExperimenalChip objects") 
	if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");

		#my $tmp_id = $self->adaptor->db->get_OligoArrayAdaptor->fetch_by_array_chip_dbID($ec->array_chip_id())->dbID();
      
		#$array_id ||= $tmp_id;
      
      #throw("You have passed ExperimentalChips from different if($array_id != $tmp_id)
      
      if(exists  $all_ids{$ec->array_chip_id()}){
	throw("Multiple chip query only works with contiguous chips within an array, rather than duplicates");
      }
      

      $all_ids{$ec->dbID()} = 1;
      $query_ids{$ec->dbID()} = 1 if(! exists $self->{'results'}{$anal_name}{$ec->dbID()});
      
    }
    
    
    my @ec_ids = keys %query_ids;
    my @all_ids = keys %all_ids;
    
    
    #warn "ec ids @ec_ids\n";
    #warn "all ids @all_ids\n";
    
    #$self->{'results'} ||= {};
    #$self->{'results_complete'} ||= 0;#do we need this now?
    
    if((scalar(@all_ids) - scalar(@ec_ids))> 1){
		throw("DATA ERROR - There is more than one result stored for the following ExperimentalChip ids: @all_ids");
	}		
	elsif(! $self->{'results'} || (($anal_name && scalar(@ec_ids) > 0) && scalar(@all_ids) == scalar(@ec_ids))){
		#fetch all, set complete set flag
		#$self->{'results_complete'} ||= 1 	if(! $anal_name);

		#would need to look up chip and channel analyses here and call relevant fetch
		#or pass the chip and then build the query as = or IN dependent on context of logic name

		#if there are multiple results, last one will overwrite others
		my @result_refs = @{$self->adaptor->fetch_results_by_probe_experimental_chips_analysis($self->probe->dbID(), \@ec_ids, $anal_name)};
	
		#could do foreach here to deal with retrieving all i.e. no logic name
		throw("Fetched more than one result for this OligoFeature, Analysis and ExperimentalChips") if (scalar(@result_refs) >1);
		#Can supply mutliple chips, but probe ids "should" be unique(in the DB at least) amongst contiguous array_chips
		#build the cache based on logic name and table_id
		#cahce key??  should we cat the ec_ids together?
		$self->{'results'}{$anal_name}{":".join(":", @ec_ids).":"} = $result_refs[0]->[0];
	}

	#do we return the ec ids here to, or do we trust that the user will know to only pass contiguous rather than duplicate chips

	#how are we going to retrieve the result for one of many possible ec id keys?
	#options, cat ec dbids as key, and grep them to find full key, then return result
	#this may hide the duplicate chip problem
	#If a query has already been made and cached,another query with one differing ID(duplicate result) may never be queried as we already have a cahced result
	#We shoulld pick up duplicates before this happens
	#If we try and mix ExperimentalChips from different experiments, then this would also cause multiple results, and hence hide some data
	
	my @keys;
	foreach my $id(@all_ids){
	  my @tmp = grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});
	  #Hacky needs sorting, quick fix for release!!

	  if(@tmp){
	    push @keys, grep(/:${id}:/, keys %{$self->{'results'}{$anal_name}});

	    last;
	  }

	}

	throw("Got more than one key for the results cache") if scalar(@keys) > 1;

    return $self->{'results'}{$anal_name}{$keys[0]};
}


#Will this be too slow, can we not do one query across all tables


1;

