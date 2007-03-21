#
# Ensembl module for Bio::EnsEMBL::Funcgen::Probe
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::Probe - A module to represent a nucleotide probe.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Probe;

#

my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
        -PROBE_SET     => $probe_set,
	    -NAME          => 'Probe-1',
        -ARRAY         => $array,
        -ARRAY_CHIP_ID => $ac_dbid,
	    -CLASS         => "EXPERIMENTAL",
);

=head1 DESCRIPTION

An Probe object represents an probe on a microarray. The data (currently the 
name, probe_set_id, length, pair_index and class) are stored
in the oligo_probe table. 

For Affy arrays, a probe can be part of more than one array, but only part of
one probeset. On each Affy array the probe has a slightly different name. For
example, two different complete names for the same probe might be
DrosGenome1:AFFX-LysX-5_at:535:35; and Drosophila_2:AFFX-LysX-5_at:460:51;. In
the database, these two probes will have the same oligo_probe_id. Thus the same
Affy probe can have a number of different names and complete names depending on
which array it is on.

=head1 AUTHOR

This module was created by Nathan Johnson, but is almost entirely based on the
Affy/OligoProbe module written by Arne Stabenau and Ian Sealy.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Probe;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-NAME]          : string - probe name
        Used when the probe is on one array.
  Arg [-NAMES]         : Listref of strings - probe names
        Used when the probe is on multiple arrays.
  Arg [-ARRAY]          : Bio::EnsEMBL::Funcgen::Array
        Used when the probe is on one array.
  Arg [-ARRAYS]         : Listref of Bio::EnsEMBL::Funcgen::Array
        Used when the probe is on multiple arrays.
  Arg [-ARRAY_CHIP_ID]  : int - array_chip db ID
        Used when the probe is on one array.
  Arg [-NAMES]          : Listref of ints - arary_chip db IDs
        Used when the probe is on multiple arrays.
  Arg [-PROBE_SET]      : Bio::EnsEMBL::ProbeSet
        Each probe is part of one(and only one) probeset, if not probe set
        then probeset = probe i.e. probe_set size = 1
  Arg [-LENGTH]         : int - probe length
        Will obviously be the same for all probes if same probe
		is on multiple arrays.
  Arg [-CLASS]          : string - probe class e.g. CONTROL, EXPERIMENTAL
        Will be the same for all probes if same probe is on
		multiple arrays.


  Example    : my $probe = Bio::EnsEMBL::Probe->new(
                   -NAME          => 'Probe-1',
				   -PROBE_SET     => $probe_set,
                   -ARRAY         => $array,
                   -ARRAY_CHIP_ID => $array_chip_id,
				   -LENGTH        => 25,
                   -CLASS         => 'EXPERIMENTAL',
      
               );
  Description: Creates a new Bio::EnsEMBL::Probe object.
  Returntype : Bio::EnsEMBL::Probe
  Exceptions : Throws if not supplied with probe name(s) and array(s)
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;
  
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
  
  my (
      $names,          $name,
      $array_chip_ids, $array_chip_id,
      $arrays,         $array,
      $probeset,       $aclass,
      $length
     ) = rearrange([
		    'NAMES',          'NAME',
		    'ARRAY_CHIP_IDS', 'ARRAY_CHIP_ID',
		    'ARRAYS',         'ARRAY',
		    'PROBE_SET',      'CLASS',
		    'LENGTH'
		   ], @_);
  
	
  @$names = ($name) if(ref($names) ne "ARRAY");
  @$array_chip_ids = ($array_chip_id) if (ref($array_chip_ids) ne "ARRAY");
  @$arrays = ($array) if (ref($arrays) ne "ARRAY");
  
  #We need to record duplicates for each probe_set i.e. each array.
  #the relationship is really array_chip to name, as everything else stays the same
  #can't have same probe_set_id as this wouldn't maintain relationship
  #need unique ps id's or array_chip_id in probe table?
  #Then we can miss probeset id's out totally if required
  #or should we just duplicate everything with unique db IDs
  
  
  if (defined $$names[0]) {
    
    if(scalar(@$names) != scalar(@$array_chip_ids)){
      throw("You have not specified valid name:array_chip_id pairs\nYou need a probe name for each Array");
    }
    
    if(defined $$arrays[0]){ 
      if(scalar(@$names) != scalar(@$arrays)){
	throw("You have not specified valid name:Array pairs\nYou need a probe name for each Array\n");
      }
    }
    else{
      warn("You have not specified and Array objects, this will result in multiple/redundant queries based on the array_chip_id\nYou should pass Array objects to speed up this process");
    }
    
    # Probe(s) have been specified
    # Different names reflect different array
    
    for my $i(0..$#{$names}){
      $self->add_array_chip_probename($$array_chip_ids[$i], $$names[$i], $$arrays[$i]);
    }
  } else {
    throw('You need to provide a probe name (or names) to create an Probe');
  }
  
  $self->probeset($probeset) if defined $probeset;
  $self->class($aclass)      if defined $aclass;
  $self->length($length)     if defined $length;
  return $self;
}

=head2 add_array_chip_probename

  Arg [1]    : int - db ID of array_chip
  Arg [2]    : string - probe name
  Arg [3]    : Bio::EnsEMBL::Funcgen::Array - optional, will get obj from DBAdaptor if not passed
  Example    : $probe->add_array_chip_probename($ac_dbid, $probename);
  Description: Adds a probe name / array pair to a probe, allowing incremental
               generation of a probe.
  Returntype : None
  Exceptions : None
  Caller     : General,
               Probe->new(),
               ProbeAdaptor->_obj_from_sth(),
			   AffyProbeAdaptor->_obj_from_sth()
  Status     : Medium Risk

=cut

sub add_array_chip_probename {
    my $self = shift;
    my ($ac_dbid, $probename, $array) = @_;
    $self->{ 'arrays'     } ||= {};
    $self->{ 'probenames' } ||= {};

    #mass redundancy here, possibility of fetching same array over and over!!!!!!!!!!!!!!
    if(! defined $array){
      $array = $self->adaptor()->db()->get_ArrayAdaptor()->fetch_by_array_chip_id($ac_dbid);
    }
    
    #mapping between probename and ac_dbid is conserved through array name between hashes
    #only easily linked from arrays to probenames,as would have to do foreach on array name
    
    $self->{ 'arrays'     }->{$ac_dbid} = $array;
    $self->{ 'probenames' }->{$array->name()} = $probename;

    return;
}


=head2 get_all_ProbeFeatures

  Args       : None
  Example    : my $features = $probe->get_all_ProbeFeatures();
  Description: Get all features produced by this probe. The probe needs to be
               database persistent.
  Returntype : Listref of Bio::EnsEMBL:Funcgen::ProbeFeature objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_ProbeFeatures {
	my $self = shift;
	if ( $self->adaptor() && $self->dbID() ) {
		return $self->adaptor()->db()->get_ProbeFeatureAdaptor()->fetch_all_by_Probe($self);
	} else {
		warning('Need database connection to retrieve Features');
		return [];
	}    
}

=head2 get_all_Arrays

  Args       : None
  Example    : my $arrays = $probe->get_all_Arrays();
  Description: Returns all arrays that this probe is part of. Only works if the
               probe was retrieved from the database or created using
			   add_Array_probename (rather than add_arrayname_probename).
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_Arrays {
    my $self = shift;
	# Do we have Array objects for this probe?
    if (defined $self->{'arrays'}) {
		return [ values %{$self->{'arrays'}} ];
    } elsif ( $self->adaptor() && $self->dbID() ) { 
		# Only have names for arrays, so need to retrieve arrays from database
		warning('Not yet implemented');
		return [];
    } else {
		warning('Need database connection to get Arrays by name');
		return [];
    }
}

=head2 get_all_probenames

  Args       : None
  Example    : my @probenames = @{$probe->get_all_probenames()};
  Description: Retrieves all names for this probe. Only makes sense for probes
               that are part of a probeset (i.e. Affy probes), in which case
			   get_all_complete_names() would be more appropriate.
  Returntype : Listref of strings
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_probenames {
    my $self = shift;
    return [ values %{$self->{'probenames'}} ];
}

=head2 get_probename

  Arg [1]    : string - array name
  Example    : my $probename = $probe->get_probename('Array-1');
  Description: For a given array, retrieve the name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is not known for this probe
  Caller     : General
  Status     : Medium Risk

=cut


#we can have dulplicate probes on same array for Nimblegen
#what defines and unique probe?
#If we have a duplicate on the same array or even on the same array_chip, then we can still return the same name
#Needs more work

sub get_probename {
    my ($self, $arrayname) = @_;
    
    if (! $arrayname){
      
      #Sanity check that there is only one non-AFFY array
      my @ac_ids = keys %{$self->{'arrays'}};

      if((scalar @ac_ids == 1) && ($self->get_all_Arrays()->[0]->vendor() ne "AFFY")){
	$arrayname = $self->get_all_Arrays()->[0]->name();
      }
      else{
	throw("Cannot retrieve probename without arrayname if more than 1 array chip(@ac_ids) and not NIMBELGEN(".$self->get_all_Arrays()->[0]->vendor().")\n");
      }
    }

	
    my $probename = $self->{'probenames'}->{$arrayname};

    if (!defined $probename) {
		throw('Unknown array name');
    }
	
    return $probename;
}

=head2 get_all_complete_names

  Args       : None
  Example    : my @compnames = @{$probe->get_all_complete_names()};
  Description: Retrieves all complete names for this probe. The complete name
               is a concatenation of the array name, the probeset name and the
			   probe name.
  Returntype : Listref of strings
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_complete_names {
    my $self = shift;
	
    my @result = ();
	
    my $probeset = $self->probeset()->name();
    $probeset .= ':' if $probeset;
	

    warn "For Nimblegen this need to be Container:Seqid::probeid?";

    #while ( my ($arrayname, $probename) = each %{$self->{'probenames'}} ) {
    while ( my (undef, $array) = each %{$self->{'arrays'}} ) {
      #would have to put test in here for $self->arrays()->vendor()
      #if($array->vendor() eq "AFFY"){
      
      push @result, $array->name().":$probeset".$self->{'probenames'}{$array->name()};
      #}
      #else{
      #	push @result, $self->{'probenames'}{$array->name()};
      #}
    }
    
    return \@result;
}



#For affy this matters as name will be different, but not for Nimblegen
#Need to consolidate this
#have get name method which throws if there is more than one array
#detects array vendor and does appropriate method

=head2 get_complete_name

  Arg [1]    : string - array name
  Example    : my $compname = $probe->get_complete_name('Array-1');
  Description: For a given array, retrieve the complete name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is not known for this probe
  Caller     : General
  Status     : Medium Risk

=cut

sub get_complete_name {
    my $self = shift;
    my $arrayname = shift;

    my $probename = $self->{'probenames'}->{$arrayname};
    if (!defined $probename) {
		throw('Unknown array name');
    }
	
	my $probeset = $self->probeset()->name();
	$probeset .= ':' if $probeset;

	

	
    return "$arrayname:$probeset$probename";
}

=head2 probeset

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::ProbeSet
  Example    : my $probe_set = $probe->probeset();
  Description: Getter and setter of probe_set attribute for Probe objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probeset {
    my $self = shift;

	#This is nesting a parent object in a child object....circular reference?
	#Need to make sure we code around this when generating ProbeSets

    $self->{'probe_set'} = shift if @_;
    return $self->{'probe_set'};
}

=head2 class

  Arg [1]    : (optional) string - class
  Example    : my $class = $probe->class();
  Description: Getter and setter of class attribute for Probe
               objects e.g. CONTROL, EXPERIMENTAL
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub class {
    my $self = shift;
    $self->{'class'} = shift if @_;
    return $self->{'class'};
}

=head2 length

  Arg [1]    : (optional) int - probe length
  Example    : my $probelength = $probe->length();
  Description: Getter and setter of length attribute for Probe
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub length {
    my $self = shift;
    $self->{'length'} = shift if @_;
    return $self->{'length'};
}


### ARRAY DESIGN SPECIFIC METHODS
### Sould really split this into BaseProbeDesign and eFGProbeDesign wrapper class?
### Quick implementation for now
### These analyses are just the core Analysis classes representing
### the output from a RunnableDB ensembl-analysis.
### Note scores added to previously stored probes may never make it into the DB?
### May also get duplicate record for previously stored probes?
### This will only happen if we're Collapsing or we fully implemnt the nr probe id to seq theory

=head2 add_Analysis_score

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_analysis_score($analysis, $score);
  Description: Setter for probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub add_Analysis_score{
    my ($self, $anal, $score) = @_;

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    throw("Must provide a score to add to the probe") if ! defined $score;

    $self->{'analysis'}{$anal->dbID()} = $score;

    return;
}

=head2 add_Analysis_CoordSystem_score

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : Bio::EnsEMBL::CoordSystem
  Arg [3]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_analysis_coord_system_score($analysis, $coord_sys, $score);
  Description: Setter for coord system dependant probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub add_Analysis_CoordSystem_score{
    my ($self, $anal, $cs, $score) = @_;

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    if(! ($cs && $cs->dbID() && $cs->isa("Bio::EnsEMBL::Funcgen::CoordSystem"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::CoordSystem");
    }

    throw("Must provide a score to add to the probe") if ! defined $score;

    $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()} = $score;

    return;
}

=head2 get_score_by_Analysis

  Arg [1]    : Bio::EnsEMBL::Analysis
  Example    : my $anal_score = $probe->get_analysis_score($analysis);
  Description: Setter for probe analysis attributes from an array design
  Returntype : string
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub get_score_by_Analysis{
  my ($self, $anal) = @_;
  
  $self->get_all_design_scores() if ! defined $self->{'analysis'};

  if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
    throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
  }


  return (exists $self->{'analysis'}{$anal->dbID()}) ? $self->{'analysis'}{$anal->dbID()} : undef;
}

=head2 get_score_by_Analysis_CoordSystem

  Arg [1]    : Bio::EnsEMBL::Analysis
  Arg [2]    : Bio::EnsEMBL::CoordSystem
  Arg [3]    : string - analysis score (as string a precision may differ between analyses)??
  Example    : $probe->add_analysis($analysis, $coord_sys, $score);
  Description: Setter for coord system dependant probe analysis attributes from an array design
  Returntype : None
  Exceptions : throws if args are not met or valid
  Caller     : General
  Status     : at risk

=cut

sub get_score_by_Analysis_CoordSystem{
    my ($self, $anal, $cs) = @_;

    $self->get_all_design_scores() if ! defined $self->{'analysis_coord_system'};

    if(! ($anal && $anal->dbID() && $anal->isa("Bio::EnsEMBL::Analysis"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Analysis");
    }

    if(! ($cs && $cs->dbID() && $cs->isa("Bio::EnsEMBL::Funcgen::CoordSystem"))){
      throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::CoordSystem");
    }

    my $score = undef;

    if(exists $self->{'analysis_coord_system'}{$anal->dbID()} &&
       exists $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()}){
      $score = $self->{'analysis_coord_system'}{$anal->dbID()}{$cs->dbID()};
    }

    return $score;
}


=head2 get_all_design_scores

  Arg [1]    : Boolean - No fetch flag, to fetch design scores from DB, used in adaptor
               To avoid testing DB for each probe when no design scores have been added.
  Example    : my @anal_score_coordsets = @{$probe->get_all_design_scores()};
  Description: Gets all design scores as analysis_id, score and optionally coord_system_id
  Returntype : ARRAYREF
  Exceptions : throws if no fetch flag is not defined and adaptor or probe is not defined and or stored.
  Caller     : General
  Status     : at risk

=cut

#not named get_all_Analysis_scores as this would imply only non-cs dependent scores
#hence named after table, as this returns simple table records

sub get_all_design_scores{
  my ($self, $no_fetch) = @_;

  my ($analysis_id, $cs_id, $score, @design_scores);

  if(! $no_fetch){#can assume we have none stored already due to implementation of add methods
    
    throw("Probe must have and adaptor to fetch design scores from the DB") if(! $self->adaptor());
    
    foreach my $probe_analysis(@{$self->adaptor->fetch_all_design_records($self)}){
      #we can't use the add methods here as would be cyclical
      #nor do we need extra validation
      
      ($analysis_id, $cs_id, $score) = @$probe_analysis;
      
      if($cs_id){
	$self->{'analysis_coord_system'}{$analysis_id}{$cs_id} = $score;
      }else{
	$self->{'analysis'}{$analysis_id} = $score;
      }
    }
  }
  
  #populate array from attrs
  if(exists $self->{'analysis_coord_system'}){
    
    foreach $analysis_id(keys %{$self->{'analysis_coord_system'}}){
      
      foreach $cs_id(keys %{$self->{'analysis_coord_system'}{$analysis_id}}){
	push @design_scores, [$analysis_id, $self->{'analysis_coord_system'}{$analysis_id}{$cs_id}, $cs_id];
      }
    }
  }
  
  if(exists $self->{'analysis'}){
    
    foreach $analysis_id(keys %{$self->{'analysis'}}){
      
      push @design_scores, [$analysis_id, $self->{'analysis'}{$analysis_id}];
    }
  }
  

  return \@design_scores;

}


#do we need get_all methods for Analysis and Analysis_CoordSystem?
#maybe if we split into another Class and Adaptor

1;

