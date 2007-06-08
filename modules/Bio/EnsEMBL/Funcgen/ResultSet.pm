#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultSet
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::ResultSet - A module to represent ResultSet.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultSet;

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new(

); 



=head1 DESCRIPTION

A ResultSet object provides access to a set raw results from an Experiment. A set will be one or more 
contiguous chips to be treated as one set, with the same analysis. Duplicate sets will form a separate
result set, as will the same raw data analysed or normalised in a different manner.


=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ResultSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-ANALYSIS] :



  Example    : my $feature = Bio::EnsEMBL::Funcgen::ResultSet->new(
                                                                   -dbid        => $dbid,
                                                                   -analysis    => $analysis,
                                                                   -table_name  => 'experimental_chip',
                                                                   -table_id    => $ec_id,
			                                          ); 
  Description: Constructor for ResultSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if no experiment_id defined
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($analysis, $table_name, $table_id, $ftype, $ctype, $name)
    = rearrange(['ANALYSIS', 'TABLE_NAME', 'TABLE_ID', 'FEATURE_TYPE', 'CELL_TYPE', 'NAME'], @_);


  $self->{'table_id_hash'} = {};

  #maybe don't need tha analysis args as mandatory as we're testing in the adaptor store method
  if (! $table_name){
    throw("Need to pass the following arg:\t-table_name");
  }

 
  
  #do we need some control of creating new objects with dbID and adding result_groups/feature_sets and them storing/updating them
  #potential for someone to create one from new using a duplicate dbID and then linking incorrect data to a pre-existing ResultGroup
  #we need to verify that each table_name/id in the set is from the same experiment


  $self->analysis($analysis) if $analysis;
  $self->table_name($table_name);
  $self->add_table_id($table_id) if $table_id;
  $self->feature_type($ftype) if $ftype;
  $self->cell_type($ctype) if $ctype;
  $self->name($name) if $name;

  return $self;
}



#methods
#set wide display label(predicted_feature) + more wordy label for wiggle tracks?
#defined by experiment type i.e. time course would require timepoint in display label
#deal with this dynamically or have display_label in table
#Need call on type, or fetch all would

#_get_ec_ids or contigsets?
#this should now be an intrinsic part of this class/adaptor

#cell line
#feature_type
#displayable...should have one for the whole set and one for each raw and predicted?

#have analysis as arg? Or do we get all analysis sets?
#we need to be able to set analyses for ResultSets dynamically from DB
#pick up all ResultSets 
#displayable field in ResultSets also?

#If we have mixed types in the same experiment then we could get promoter features and histone wiggle tracks displayed togeter
#Not v.good for display purposes?  We may want to separate the promoter and histone tracks, or we may want ll the experiment data together but of mixed types.
#We need to be able to pull back the experiment type for each set, therefore this needs setting on an ec level, not an experiment level.
#This is also v.reliant on putting contig set info in place, otherwise we may get mixed chip types in same set.

#get_raw_analysis_name
#get_predicted_feature_analysis_name
#set ResultFeatures and PredictedFeatures in hash keyed by analysis_name?



=head2 table_name

  Arg [1]    : (optional) string - table_name (experimental_chip or channel)
  Example    : $result_set->experiment_id($exp_id);
  Description: Getter and setter for the table_name for this ResultSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub table_name{
    my $self = shift;

    if (@_){
      
      if($self->{'table_name'} && ($self->{'table_name'} ne $_[0])){
	throw("Cannot mix  table name/types of a ResultSet");
      }
	
      $self->{'table_name'} = $_[0];
    }

    return $self->{'table_name'};
}


=head2 name

  Arg [1]    : (optional) string - name of the result set
  Example    : $result_set->name($name);
  Description: Getter and setter for the name attribute for this ResultSet.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub name{
    my $self = shift;

    if (@_){
      $self->{'name'} = shift;
	}

    return $self->{'name'};
}



=head2 analysis

  Arg [1]    : (optional) - Bio::EnsEMBL::Analysis
  Example    : $anal_id = $rset->analysis->dbID();
  Description: Getter and setter for the analysis attribute for this ResultSet.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub analysis {
  my $self = shift;
	
  if(@_){
	  throw("Must pass a valid Bio::EnsEMBL::Analysis object") if (! $_[0]->isa("Bio::EnsEMBL::Analysis"));
	  $self->{'analysis'} = shift;
  }
		
  return $self->{'analysis'};
}


=head2 feature_type

  Arg [1]    : (optional) - Bio::EnsEMBL::Funcgen::FeatureType
  Example    : $fname = $rset->feature_type->name();
  Description: Getter and setter for the feature_type attribute for this ResultSet.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : Throws if arg is not a valid Bio::EnsEMBL::Funcgen::FeatureType
  Caller     : General
  Status     : At Risk

=cut


sub feature_type {
  my $self = shift;
	
  if(@_){
    throw("Must pass a valid Bio::EnsEMBL::Funcgen::FeatureType object") if (! $_[0]->isa("Bio::EnsEMBL::Funcgen::FeatureType"));
    $self->{'feature_type'} = shift;
  }
		
  return $self->{'feature_type'};
}



=head2 cell_type

  Arg [1]    : (optional) - Bio::EnsEMBL::Funcgen::CellType
  Example    : $cname = $rset->cell_type->name();
  Description: Getter and setter for the cell_type attribute for this ResultSet.
  Returntype : Bio::EnsEMBL::Funcgen::CellType
  Exceptions : Throws if arg is not a valid Bio::EnsEMBL::Funcgen::CellType
  Caller     : General
  Status     : At Risk

=cut


sub cell_type {
  my $self = shift;
	
  if(@_){
    throw("Must pass a valid Bio::EnsEMBL::Funcgen::CellType object") if (! $_[0]->isa("Bio::EnsEMBL::Funcgen::CellType"));
    $self->{'cell_type'} = shift;
  }
		
  return $self->{'cell_type'};
}




=head2 add_table_id

  Example    : $result_set->add_table_id($ec_id, $cc_id);
  Description: Caches table_id chip_channel_id to the ResultSet.
               The unique chip_channel_id is used to key into the result table,
               it also reduces redundancy and enable mapping of results to chips
               rather than just the ResultSet.  This enables result retrieval
               based on chips in the same set which  have a differing status.
  Returntype : None
  Exceptions : Throws if no table_id defined
  Caller     : General
  Status     : At Risk

=cut

sub add_table_id {
  my ($self, $table_id, $cc_id) = @_;
  
  if (! defined $table_id){	
    throw("Need to pass a table_id");
  }else{
    
    if((exists $self->{'table_id_hash'}->{$table_id}) && (defined $self->{'table_id_hash'}->{$table_id})){
      throw("You are attempting to redefine a chip_channel_id which is already defined");
    }
    
    $self->{'table_id_hash'}->{$table_id} = $cc_id;    
    
  }

  return;
}


=head2 table_ids

  Example    : $result_set->feature_group_id($fg_id);
  Description: Getter and setter for the feature_group_id for this ResultSet.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub table_ids {
  my $self = shift;
  
  return [ keys %{$self->{'table_id_hash'}} ];
}

=head2 chip_channel_ids

  Example    : my @rset_cc_ids = @{$result_set->chip_channel_ids()};
  Description: Getter for the chip channel ids for this ResultSet.
  Returntype : listref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub chip_channel_ids {
  my $self = shift;
  
  return [ values %{$self->{'table_id_hash'}} ];
}

=head2 contains

  Example    : if($result_set->contains($chip_or_channel)){...do some chip or channel erpartions here...};
  Description: Returns true if the given Channel or ExperimentalChip is part of this ResultSet
  Returntype : boolean
  Exceptions : warns if ResultSet table name is not of argument type
  Caller     : General
  Status     : At Risk

=cut


sub contains{
  my ($self, $chip_channel) = @_;

  my $contains = 0;
  my @tables = $chip_channel->adaptor->_tables();
  my ($table_name, undef) = @{$tables[0]};

  if($table_name ne $self->table_name()){
    warn("ResultSet(".$self->table_name().") cannot contain ${table_name}s");
  }else{
    $contains = 1 if (exists $self->{'table_id_hash'}->{$chip_channel->dbID()});
  }
  
  return $contains;
}

=head2 get_chip_channel_id

  Arg [1]    : int - ExperimentalChip dbID
  Example    : $result_set->get_chip_channel_id($ec_id);
  Description: Retrieves a chip_channel_id from the cahce given an ExperimentalChip dbID
  Returntype : int
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub get_chip_channel_id{
  my ($self, $table_id) = @_;
  
  return (exists $self->{'table_id_hash'}->{$table_id}) ?  $self->{'table_id_hash'}->{$table_id} : undef;
}

=head2 get_ExperimentalChips

  Example    : my @ecs = @{$result_set->get_ExperimentalChips()};
  Description: Retrieves a chip_channel_id from the cahce given an ExperimentalChip dbID
  Returntype : Listref of ExperimentalChip object
  Exceptions : warns is not an experimental_chip ResultSet
  Caller     : General
  Status     : At Risk

=cut

sub get_ExperimentalChips{
  my $self = shift;
  
  if(! defined $self->{'experimental_chips'}){
    my $ec_adaptor = $self->adaptor->db->get_ExperimentalChipAdaptor();
    
	if($self->table_name() eq "experimental_chip"){

	  foreach my $ec_id(@{$self->table_ids()}){
#		  warn "Getting ec with id $ec_id";
		push @{$self->{'experimental_chips'}}, $ec_adaptor->fetch_by_dbID($ec_id);
		#should this be hashed on chip_channel_id?
	  }
	}else{
	  #warn("Retrieving ExperimentalChips for a Channel ResultSet");
	  
	  my %echips;
	  my $chan_adaptor = $self->adaptor->db->get_ChannelAdaptor(); 
		  
	  foreach my $chan_id(@{$self->table_ids()}){
		my $chan = $chan_adaptor->fetch_by_dbID($chan_id);
		$echips{$chan->experimental_chip_id} ||= $ec_adaptor->fetch_by_dbID($chan->experimental_chip_id);
	  }
	  
	  @{$self->{'experimental_chips'}} = values %echips;
	}
  }

  return $self->{'experimental_chips'};
}



=head2 get_replicate_set_by_chip_channel_id

  Arg[0]     : int - chip_channel_id
  Example    : my $rep_set_name = $result_set->get_replicate_set_by_chip_channel_id($cc_id);
  Description: Retrieves the replicate set name defined by the corresponding ExperimentalChip
  Returntype : String - replicate set name
  Exceptions : 
  Caller     : General
  Status     : At Risk - implement for Channels?

=cut


sub get_replicate_set_by_chip_channel_id{
  my ($self, $cc_id) = @_;

  if( ! defined $self->{'_replicate_cache'}){

	warn "Generating replicate cache!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";


	foreach my $ec (@{$self->get_ExperimentalChips()}){
	  
	  $self->{'_replicate_cache'}{$self->get_chip_channel_id($ec->dbID())} = $ec->replicate();
	  

	}
  }


  #warn here of absent replicate info?

  return (exists $self->{'_replicate_cache'}{$cc_id}) ?  $self->{'_replicate_cache'}{$cc_id} : undef;

}


=head2 display_label

  Example    : print $rset->display_label();
  Description: Getter for the display_label attribute for this ResultSet.
               This is more appropriate for teh predicted_features of the set.
               Use the individual display_labels for each raw result set.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub display_label {
  my $self = shift;
  
  if(! $self->{'display_label'}){
    
    #This should display some info about the chip set/duplicte set if there is more than one set of data for a feature_set!!!!!!!!!!!!!!!
    
    #Some tomfoolery here to accomdate sets which we do not know the feature or cell type for.
    #should we make cell_type and feature_type mandatory?

    if(defined $self->feature_type()){
      $self->{'display_label'} = $self->feature_type->name()." - ";
    }else{
       $self->{'display_label'} = "FEATURE TYPE NOT KNOWN - ";
    }

    if(defined $self->cell_type()){
      $self->{'display_label'} .= ($self->cell_type()->display_label()) ? $self->cell_type->display_label() : $self->cell_type->name();
    }else{
      $self->{'display_label'} .= "CELL TYPE NOT KNOWN";
    }

    $self->{'display_label'} .= " Enriched Sites";
  }
  
  return $self->{'display_label'};
}


=head2 get_displayable_ResultFeatures_by_Slice

  Example    : my @results = @{$ResultSet->get_all_displayable_results()};
  Description: wrapper to get_all_results with displayable flag passed
  Returntype : List ref to array refs containing ($display_label, @result_features);
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_displayable_ResultFeatures_by_Slice{
  my ($self, $slice) = @_;
  return $self->get_ResultFeatures_by_Slice($slice, 'DISPLAYABLE');
}


#Is it possible to to have one chip displayable and another not within the same set.
#one may fail validation, but would just omit from result_set?
#Or maybe we only want certain chips displayed, locational context is not split evenly over chips so not that useful :?
#If we handle this here then all we have to do is filter table_ids first before passing them to the fetch.

sub get_ResultFeatures_by_Slice{
  my ($self, $slice, $status) = @_;


  #this does not store the ResultFeatures anywhere so we need to be mindful that calling this will always result in a DB query
  #Can we generate a slice hash for this?
  #This would eat more memory but may be faster for web display if we can support some sort of dynamic querying with respect to expanded/shifted slice and cached ResultFeatures.
  #would have to cater for 49 bp overhang which may result in duplicate ResultFeatures which overlap ends of adjacent/overlapping slices
  #49bp ? should be altered to cope with overlap properly.


  #do we need to set the replicate_cache here first, or should this be done in new?

  return $self->adaptor->fetch_ResultFeatures_by_Slice_ResultSet($slice, $self, $status);
}



1;
