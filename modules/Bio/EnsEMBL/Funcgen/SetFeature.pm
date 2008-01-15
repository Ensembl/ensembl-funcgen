# Ensembl module for Bio::EnsEMBL::Funcgen::SetFeature
#
# Copyright (c) 2007 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Funcgen::SetFeature - Ensembl specific set feature.

=head1 SYNOPSIS

    my $feat = new Bio::EnsEMBL::Feature(-start         => 100,
                                         -end           => 220,
                                         -strand        => -1,
                                         -slice         => $slice,
                                         -feature_set   => $fset,
                                         -display_label => $label,
                                      );

    my $start  = $feat->start;
    my $end    = $feat->end;
    my $strand = $feat->strand;

    #move the feature to the chromosomal coordinate system
    $feature = $feature->transform('chromosome');

    #move the feature to a different slice (possibly on another coord system)
    $feature = $feature->transfer($new_slice);

    #project the feature onto another coordinate system possibly across
    #boundaries:
    @projection = @{$feature->project('contig')};

    #change the start, end, and strand of the feature in place
    $feature->move($new_start, $new_end, $new_strand);

=head1 DESCRIPTION

This is a simple wrapper method for the core Feature class to contain common generic
Funcgen SetFeature methods.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::SetFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [-FEATURE_SET]  : Bio::EnsEMBL::Funcgen::FeatureSet
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
  Example             : my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
										                                  -SLICE         => $chr_1_slice,
									                                      -START         => 1_000_000,
                      				                                      -END           => 1_000_024,
									                                      -STRAND        => -1,
									                                      -DISPLAY_LABEL => $text,
									                                      -FEATURE_SET   => $fset,
                                                                                  );


  Description: Constructor for SetFeature objects. Should never be called directly, only by children.
  Returntype : Bio::EnsEMBL::Funcgen::SetFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);
  
  my ($display_label, $fset, $ftype)
    = rearrange(['DISPLAY_LABEL', 'FEATURE_SET', 'FEATURE_TYPE'], @_);


  #should have feature_type & cell_type here too?

  
  $self->display_label($display_label) if $display_label;

  if($ftype){

	if(! (ref($ftype) && $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType'))){
	  throw('feature_type param must be a valid Bio::EnsEMBL::Funcgen::FeatureType');
	}

	$self->{'feature_type'} = $ftype;
  }

  if(! (ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"))){
	throw("Must pass valid Bio::EnsEMBL::Funcgen::FeatureSet object");
  }
  $self->{'feature_set'}= $fset;

		
  return $self;
}


=head2 feature_set

  Arg [1]    : (optional) Bio::EnsEMBL::FeatureSet 
  Example    : $efeature->feature_set($fset);
  Description: Getter for the FeatureSet attribute for this feature. 
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub feature_set {
  my $self = shift;

  return $self->{'feature_set'};
}




=head2 get_all_DBEntries

  Example    : my @dbentries = @{ $gene->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this transcript.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranscriptAdaptor::store
  Status     : Stable

=cut

sub get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  #Need to add tests for valid objects for xrefs

  # if not cached, retrieve all of the xrefs for this gene

  my @tables = $self->adaptor->_tables;

  if(!defined $self->{$cache_name} && $self->adaptor()) {

	my @tables = $self->adaptor->_tables;
	@tables = split/_/, $tables[0]->[0];
	my $object_type = join('', (map ucfirst($_), @tables));
	
    $self->{$cache_name} = 
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntery->new(...);
               $transcript->add_DBEntry($dbe);
  Description: Associates a DBEntry with this transcript. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               storable (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the calling code 
               is very disciplined.
  Returntype : Bio::EnsEMBL::Funcgen::SetFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  return bless ($_[1], $_[0]);
}


=head2 cell_type

  Example    : my $cell_name = $efeature->cell_type()->name();
  Description: Getter for the cell_type attribute for this feature.
               May not always be set for ExternalFeatures.
  Returntype : Bio::EnsEMBL::Funcgen:CellType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub cell_type{
	my $self = shift;

	return $self->feature_set->cell_type();
}

=head2 feature_type

  Example    : my $ft_name = $efeature->feature_type()->name();
  Description: Getter for the feature_type attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen:FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub feature_type{
  my $self = shift;
  
  #why is this not a setter?
  #this should only be set in new

  return (defined $self->{'feature_type'}) ?  $self->{'feature_type'} : $self->feature_set->feature_type();
}


#redefined from Feature to use FeatureSet, unless specifically set

=head2 analysis

  Example    : my $analysis = $efeature->feature_type()->name();
  Description: Getter for the type attribute for this feature.
  Returntype : Bio:EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if analysis passed is not a valid Bio::EnsEMBL::Analysis
  Caller     : General
  Status     : At risk

=cut

sub analysis{
  my $self = shift;


  #this is to allow multi analysis sets, but the adaptor currently  throws if they are not the same on store
  if(@_){

    if($_[0]->isa("Bio::EnsEMBL::Analysis")){
      $self->{'analysis'} = $_[0];
    }else{
      throw("Must pass a valid Bio::EnsEMBL::Analysis");
    }
    
  }

  return (defined $self->{'analysis'}) ? $self->{'analysis'} : $self->feature_set->analysis();
}


1;
