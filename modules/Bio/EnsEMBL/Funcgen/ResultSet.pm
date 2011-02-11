#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultSet
#

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::ResultSet - A module to represent ResultSet.
 

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultSet;

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new(
                                                       -dbid        => $dbid,
                                                       -analysis    => $analysis,
                                                       -table_name  => 'experimental_chip',
                                                       -table_id    => $ec_id,
); 



=head1 DESCRIPTION

A ResultSet object provides access to a set raw results from an Experiment. A set will be one or more 
contiguous chips to be treated as one set, with the same analysis. Duplicate sets will form a separate
result set, as will the same raw data analysed or normalised in a different manner.

=cut

#To do
#Change add_table_id to add_ExperimentalChip_Channel?


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::ResultSet;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate);
use Bio::EnsEMBL::Funcgen::Set;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);


=head2 new

  Arg [-ANALYSIS] :



  Example    : my $feature = Bio::EnsEMBL::Funcgen::ResultSet->new(
                                                                   -dbid        => $dbid,
                                                                   -analysis    => $analysis,
                                                                   -table_name  => 'experimental_chip',
                                                                   -table_id       => $ec_id,
                                                                   -result_feature_set => 1,
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
	
  my $self = $class->SUPER::new(@_, ('-feature_class' => 'result'));
	
  my ($table_name, $table_id, $rf_set, $dbfile_data_dir)
    = rearrange(['TABLE_NAME', 'TABLE_ID', 'RESULT_FEATURE_SET', 'DBFILE_DATA_DIR'], @_);

  $self->{'table_id_hash'} = {};

  #maybe don't need tha analysis args as mandatory as we're testing in the adaptor store method
  if (! $table_name){
    throw("Need to pass the following arg:\t-table_name");
  }
  
  #do we need some control of creating new objects with dbID and adding result_groups/feature_sets and them storing/updating them
  #potential for someone to create one from new using a duplicate dbID and then linking incorrect data to a pre-existing ResultGroup
  #we need to verify that each table_name/id in the set is from the same experiment
   $self->table_name($table_name);
  $self->add_table_id($table_id) if $table_id;
  $self->result_feature_set($rf_set) if $rf_set;
  $self->dbfile_data_dir($dbfile_data_dir) if $dbfile_data_dir;

  return $self;
}


#These are CollectionContainer? methods
#For a core track the would probably be in the Analysis
#All other collection methods are in ResultFeatureAdaptor(and parents)

=head2 get_dbfile_path_by_window_size

  Arg[1]     : int - window size
  Arg[2]     : OPTIONAL Bio::EnsEMBL::Slice Used when generating individual seq_region Collections
  Example    : my $filepath = $self->get_dbfile_path_by_ResultSet_window_size($rset, $wsize);
  Description: Generates the default dbfile path for a given ResultSet and window_size
  Returntype : string
  Exceptions : Throws if Slice is not valid
  Caller     : general
  Status     : At risk

=cut

sub get_dbfile_path_by_window_size{
  my ($self, $window_size, $slice) = @_;

  if($slice){
   
	if(! (ref($slice) && $slice->isa("Bio::EnsEMBL::Slice"))){
	  throw('You must provide a valid Bio::EnsEMBL::Slice');
	}

	$window_size .= '.'.$slice->seq_region_name;
  }

  return $self->dbfile_data_dir.'/result_features.'.$self->name.'.'.$window_size.'.col';
}


=head2 dbfile_data_dir

  Arg[1]     : OPTIONAL string - data directory for this ResultSet
  Example    : my $dbfile_data_dir = $self->dbfile_data_dir;
  Description: Getter/Setter for the root dbfile data directory for this ResultSet
  Returntype : string
  Exceptions : None
  Caller     : self
  Status     : at risk

=cut



sub dbfile_data_dir{
  my ($self, $data_dir) = @_;

  $self->{'dbfile_data_dir'} = $data_dir if defined $data_dir;

  return $self->{'dbfile_data_dir'};
}



=head2 result_feature_set

  Arg [1]    : optional - boolean 0 or 1.
  Example    : if($rset->result_feature_set){ ...use result_feature table ...};
  Description: Getter and setter for the result_feature_set attribute.
  Returntype : boolean
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub result_feature_set{
  my $self = shift;

  $self->{'result_feature_set'} = shift if @_;;
  return $self->{'result_feature_set'};
}


=head2 table_name

  Arg [1]    : (optional) string - table_name (experimental_chip, channel or input_set)
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



=head2 add_table_id

  Example    : $result_set->add_table_id($ec_id, $cc_id);
  Description: Caches table_id result_set_input_id to the ResultSet. In the case of an 
               array ResultSet, the unique result_set_input_id is used to key into the 
               result table, it also reduces redundancy and enable mapping of results to chips
               rather than just the ResultSet.  This enables result retrieval
               based on chips in the same set which  have a differing status.
               In the case of a sequencing ResultSet, this simply refers to the InputSet ids.
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
      throw("You are attempting to redefine a result_set_input_id which is already defined");
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


sub chip_channel_ids {
  my $self = shift;

  deprecate('ResultSet::chip_channel_ids is deprecated, please use result_set_input_ids');
  
  return $self->result_set_input_ids;
}

=head2 result_set_input_ids

  Example    : my @rset_rsi_ids = @{$result_set->result_set_input_ids()};
  Description: Getter for the input ids for this ResultSet.
  Returntype : arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub result_set_input_ids {
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

=head2 get_result_set_input_id

  Arg [1]    : int - dbID (experimental_chip, channel or input_set)
  Example    : $result_set->get_result_set_input_id($ec_id);
  Description: Retrieves a result_set_input_id from the cache given a dbID
  Returntype : int
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub get_result_set_input_id{
  my ($self, $table_id) = @_;
  
  return (exists $self->{'table_id_hash'}->{$table_id}) ?  $self->{'table_id_hash'}->{$table_id} : undef;
}


sub get_chip_channel_id{
  my ($self, $table_id) = @_;
  
  deprecate('ResultSet::get_chip_channel_ids is dperecated, please us get_result_set_input_id');
  return $self->get_result_set_input_id($table_id);
}



=head2 get_InputSets

  Example    : my @ecs = @{$result_set->get_ExperimentalChips()};
  Description: Retrieves a chip_channel_id from the cahce given an ExperimentalChip dbID
  Returntype : Listref of ExperimentalChip object
  Exceptions : warns is not an experimental_chip ResultSet
  Caller     : General
  Status     : At Risk

=cut

sub get_InputSets{
  my $self = shift;
  
  if($self->table_name ne 'input_set'){
	warn 'Cannot get_InputSets for an array based ResultSet';
	return;
  }

  

  if(! defined $self->{'input_sets'}){
    my $is_adaptor = $self->adaptor->db->get_InputSetAdaptor();
    
	foreach my $is_id(@{$self->table_ids()}){	
	  push @{$self->{'input_sets'}}, $is_adaptor->fetch_by_dbID($is_id);
	}
  }

  return $self->{'input_sets'};
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
  
  if($self->table_name eq 'input_set'){
	warn 'Cannot get_ExperimentalChips for an InputSet ResultSet';
	return;
  }

  if(! defined $self->{'experimental_chips'}){
    my $ec_adaptor = $self->adaptor->db->get_ExperimentalChipAdaptor();
    
	if($self->table_name() eq "experimental_chip"){

	  foreach my $ec_id(@{$self->table_ids()}){
		#warn "Getting ec with id $ec_id";
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



=head2 get_replicate_set_by_result_set_input_id

  Arg[0]     : int - chip_channel_id
  Example    : my $rep_set_name = $result_set->get_replicate_set_by_result_set_input_id($cc_id);
  Description: Retrieves the replicate set name defined by the corresponding ExperimentalChip
  Returntype : String - replicate set name
  Exceptions : 
  Caller     : General
  Status     : At Risk - implement for Channels?

=cut

#Where is this used?

sub get_replicate_set_by_result_set_input_id{
  my ($self, $cc_id) = @_;

  if( ! defined $self->{'_replicate_cache'}){

	warn "Generating replicate cache!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";


	foreach my $ec (@{$self->get_ExperimentalChips()}){
	  
	  $self->{'_replicate_cache'}{$self->get_result_set_input_id($ec->dbID())} = $ec->replicate();
	  

	}
  }


  #warn here of absent replicate info?

  return (exists $self->{'_replicate_cache'}{$cc_id}) ?  $self->{'_replicate_cache'}{$cc_id} : undef;

}

sub get_replicate_set_by_chip_channel_id{
  my ($self, $cc_id) = @_;

  deprecate('Please use get_replicate_set_by_result_set_input_id instead');
  return $self->get_replicate_set_by_result_set_input_id($cc_id);
}


=head2 get_displayable_ResultFeatures_by_Slice

  Arg[1]     : Bio::EnsEMBL::Slice
  Arg[2]     : Boolean - with probe flag, will nest Probe object in ResultFeature 
  Example    : my @results = @{$ResultSet->get_all_displayable_ResultFeatures_by_Slice($slice)};
  Description: Simple wrapper method for ResultFeatureAdaptor::fetch_all_by_Slice_ResultSet
  Returntype : Arrayref of ResultFeatures
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut


sub get_displayable_ResultFeatures_by_Slice{
  my ($self, $slice, $with_probe, $max_bins, $window_size, $constraint) = @_;
  return $self->adaptor->fetch_ResultFeatures_by_Slice_ResultSet($slice, $self, 'DISPLAYABLE', $with_probe, $max_bins, $window_size, $constraint);
}




=head2 get_ResultFeatures_by_Slice

  Arg[1]     : Bio::EnsEMBL::Slice
  Arg[2]     : string - Status name e.g. 'DISPLAYABLE'
  Arg[3]     : Boolean - with probe flag, will nest Probe object in ResultFeature 
  Arg[4]     : int - Max bins i.e. pixel width of display
  Arg[5]     : int - window_size
  Arg[6]     : string - constraint
  Example    : my @rfs_with_rpobe = @{$ResultSet->get_all_ResultFeatures_by_Slice($slice, undef, 1)};
  Description: Simple wrapper method for ResultFeatureAdaptor::fetch_all_by_Slice_ResultSet
  Returntype : Arrayref of ResultFeatures
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_ResultFeatures_by_Slice{
  my ($self, $slice, $status, $with_probe, $max_bins, $window_size, $constraint) = @_;
  return $self->adaptor->db->get_ResultFeatureAdaptor->fetch_all_by_Slice_ResultSet($slice, $self, $status, $with_probe, $max_bins, $window_size, $constraint);
}



#Floats unpack inaccurately so need 3 sigfiging
#This should match the format in which they are originally stored
#This is dependant on ResultSet type i.e. reads or intensity?
#No format for reads!
#Should this be set in the ResultSet instead?
#It may be more efficient for the caller to test for format first rather than blindly printf'ing
#even if there is no format?
#This needs setting in new, so we don't have to eval for every score.

sub score_format{
  return '%.3f';
}




=head2 log_label

  Example    : print $rset->log_label();
  Description: Get a string of the unique key fields for logging purposes
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub log_label {
  my $self = shift;

  my $label;
  
  if(defined $self->feature_type()){
	$label = $self->feature_type->name.":";
  }else{
	$label = "Unknown FeatureType:";
  }

  if(defined $self->cell_type()){
	$label .= $self->cell_type->name;
  }else{
	$label .= "Uknown CellType";
  }

  return $self->name.":".$self->analysis->logic_name.":".$label;
}



1;
