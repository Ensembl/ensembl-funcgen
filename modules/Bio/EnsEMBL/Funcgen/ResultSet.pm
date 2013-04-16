#
# Ensembl module for Bio::EnsEMBL::Funcgen::ResultSet
#

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::ResultSet - Represents a set of experimental results (signal).


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ResultSet;

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new
  (
   -dbid        => $dbid,
   -analysis    => $analysis,
   -table_name  => 'input_set',
   -table_id    => $input_set_id,
   -feature_class => 'dna_methylation',
  );



=head1 DESCRIPTION

A ResultSet object provides access to a set raw or processed signal values from an Experiment,
and can consist of any number of inputs which can either be InputSets or ExperimentalChips.
A ResultSet can represent a single or multiple merged replicates of raw or processed(normalised)
signal data.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Set

=cut

package Bio::EnsEMBL::Funcgen::ResultSet;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::Set;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Set);

#Valid enum field hashes to prevent loading NULLs

#result_set.feature_class
my %valid_classes =
  (
   result => undef,
   dna_methylation => undef,
  );

#result_set_input.table_name
my %valid_table_names =
  (
   experimental_chip => undef,
   input_set         => undef,
   channel           => undef,
  );


=head2 new

  Arg [-ANALYSIS] :



  Example    : my $feature = Bio::EnsEMBL::Funcgen::ResultSet->new
                 (
                  -dbid               => $dbid,
                  -analysis           => $analysis,
                  -support            => \@support_objects,
                  -result_feature_set => 1,
                  -feature_class      => 'result',
                  -feature_type       => $ftype,
	             );
  Description: Constructor for ResultSet objects.
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : Throws if feature class is not valid.
               Throws if table name not defined.
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  my ($table_name, $table_id, $rf_set, $dbfile_data_dir, $support)
    = rearrange(['TABLE_NAME', 'TABLE_ID', 'RESULT_FEATURE_SET', 'DBFILE_DATA_DIR', 'SUPPORT'], @_);
  # TEST MANDATORY PARAMS

  #explicit type check here to avoid invalid types being imported as NULL
  #and subsequently throwing errors on retrieval
  my $type = $self->feature_class;
  $self->{table_id_hash}      = {};

  if ( !( $type && exists $valid_classes{$type} ) ) {
    throw( 'You must define a valid ResultSet type e.g. ' .
           join( ', ', keys %valid_classes ) );
  }

  #todo fully deprecate add_table_id and -table_name and -table_id
  #There is some redundancy here between add_table_id and -support methods
  #and the add_table_id method currently allows empty ResultSets
  #this may have been to support an import where we don't have access to the support objects? 
  #check before removal
        
  if(defined $support){ 
    $self->_set_support($support, $table_name);     
  }
  else{
    deprecate('The -table_name, -table_id params and the add_table_id methods are now deprecated. Please use -support');
    
    if (! (defined $table_name &&
           exists $valid_table_names{$table_name}) ){
       throw('You need to pass a valid -table_name e.g. '.
            join(', ', keys %valid_table_names));
    }
  
    #This will currently allow a ResultSet with no support/table_ids   
    $self->_add_table_id($table_id)           if defined $table_id;
  }
  
  $self->{table_name}         = $table_name;
  $self->{result_feature_set} = (defined $rf_set) ? 1 : 0;
  $self->dbfile_data_dir($dbfile_data_dir) if defined $dbfile_data_dir;

  return $self;
}

=head2 reset_relational_attributes

  Arg[1] : Hashref containing the following mandatory parameters:
            -analysis     => Bio::EnsEMBL::Analysis,
            -feature_type => Bio::EnsEMBL::Funcgen::FeatureType,
            -cell_type    => Bio::EnsEMBL::Funcgen::CellType,
            -support      => Arrayref of valid support objectse.g. InputSet

  Description: Resets all the relational attributes of a given ResultSet. 
               Useful when creating a cloned object for migration beween DBs 
  Returntype : None
  Exceptions : Throws if any of the parameters are not defined or invalid.
  Caller     : Migration code
  Status     : At risk

=cut

sub reset_relational_attributes{
  my ($self, $params_hash, $no_db_reset) = shift;
  
  #This uses new support param, rather that add_table_id/table_name/table_id
   
  my ($support, $analysis, $feature_type, $cell_type)
    = rearrange(['SUPPORT', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE'], @$params_hash);
  
  #flush table ID cache and add support
  $self->{table_id_hash} = undef;
  $self->_add_support($support);
  
  #is_stored (in corresponding db) checks will be done in store method
 
  if(! (defined $analysis &&
        ref($analysis) eq 'Bio::EnsEMBL::Analysis') ){
    throw('You must pass a valid Bio::EnsEMBL::Analysis');
  }
  
  if(! (defined $feature_type &&
        ref($feature_type) eq 'Bio::EnsEMBL::Funcgen::FeatureType') ){
    throw('You must pass a valid Bio::EnsEMBL::FeatureType');
  }
  
  
  if(! (defined $cell_type &&
        ref($cell_type) eq 'Bio::EnsEMBL::CellType') ){
    throw('You must pass a valid Bio::EnsEMBL::CellType');
  }
  
  $self->{cell_type}    = $cell_type;
  $self->{feature_type} = $feature_type;
  $self->{analysis}     = $analysis;
  
  #Finally undef the dbID and adaptor by default
  if(! $no_db_reset){
    $self->{adaptor} = undef;
    $self->{dbID}    = undef;
  }
  
  return;
}

sub add_support{
  my($self, $support, $table_name) = @_;

  if(! (defined $support &&
        (ref($support) eq 'ARRAY')) ){
    throw('You must pass an Arrayref of support objects');
  }
  
  #Check passed table name matches objects
  #remove this once we remove add_table_ids
  #although we still need to check all support object are the same
  my $tmp_table_name = $table_name; 
  my $sdba = $support->[0]->adaptor;  
  (my $class_name = ref($support->[0])) =~ s/.*://g; 
  
  if(! $sdba){
    throw("No $class_name adaptor found. All ResultSet -support must be stored");
  }

  $table_name = $sdba->_main_table->[0];    

  if(defined $tmp_table_name &&
    ( $tmp_table_name ne $table_name) ){
     throw('Specified -table_name does not match -support table name. '.
          'Please omit -table_name and ensure all -support objects are of the same class'); 
  }
  
  #Validate table name
  if(! exists $valid_table_names{$table_name}){
    throw("The -table_name or -support objects do not refer to a valid".
          " support table:\t$table_name\nValid table names:\t".
          join(', ', keys %valid_table_names));
  }
  
  #Add the table_ids and set the support 
  foreach my $sset(@$support){
    
    if (! ( defined $sset && 
            ref($sset) && 
            $sset->isa('Bio::EnsEMBL::Funcgen::'.$class_name)) ){
         throw(ref($sset)." is not a valid support type for this $class_name ResultSet\n");
    }
    
    #we don't have access to is_stored_and_valid here
    #so this may throw an error        
    $self->_add_table_id($sset->dbID);   
  }
  
  $self->{support} = $support;  
  
  return;
}


=head2 experimental_group

  Example    : my $rset_exp_group = $rset->experimental_group;
  Description: Convenience method to get the unique experimental group name
               for this ResutlSet. Only works for ResultSets supported by InputSets.
  Returntype : String or undef (if no unique name found).
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

#Need a better way for handling source refs/track info for tracks which can share the same analysis
#Is this the analysis description right place to be putting source refs?
#Yes if we are only ever going to have one source/signal track
#ChIP-Seq ResultFeature tracks have merged sources and only have source refs corresponding peak zmenus
#This will also appear in the track info and config

sub experimental_group{
  my $self = shift;

  if(! exists $self->{experimental_group}){
    #Not undef check as undef is a valid value
    #for mixed project ResultSets

    if($self->table_name ne 'input_set'){
      throw('Cannot currently get ExperimentalGroup for a ResultSet with non-InputSet support'); 
    }

    my $exp_group;
    my @isets = @{$self->get_InputSets};

    if(@isets){
      $exp_group = $isets[0]->get_Experiment->get_ExperimentalGroup->name;

      foreach my $iset(@isets){

        if($exp_group ne $iset->get_Experiment->get_ExperimentalGroup->name){
          #Mixed experimental_group ResultSet
          $exp_group = undef;
          last;
        }
      }
    }

    $self->{experimental_group} = $exp_group;
  }

  return $self->{experimental_group};
}


=head2 display_label

  Example    : print $set->display_label;
  Description: Getter for the display_label attribute for this Set.
  Returntype : String
  Exceptions : Warns if feature_class not explicitly handled
  Caller     : General
  Status     : At Risk

=cut

#Add project name in here
#Where is this used in the webcode? Track names?
#Not for ChIP-Seq DataSets as these use the merged class descriptions
#e.g. Open Chromatin & TFBS
#Just for dna_methylation, as these are always show on one track.

sub display_label {
  my $self = shift;

  if(! defined $self->{display_label}){

    if($self->feature_class eq 'result'){ # We have a ChIP/DNase1/FAIRE signal set
      $self->{display_label} = $self->feature_type->name.' '.$self->cell_type->name.' signal';
    }
    elsif($self->feature_class eq 'dna_methylation'){

      my $project = $self->project || '';

      if($project){
        $project = ' '.$project;
      }

      $self->{display_label} = $self->analysis->display_label.' '.$self->cell_type->name.$project;
    }
    else{
      warn 'ResultSet feature class('.$self->feature_class.') not explicitly handled in display_label';
      $self->{display_label} = $self->feature_type->name.' '.$self->cell_type->name;
    }
  }

  return $self->{display_label};
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

  Arg[1]     : String (optional) - data directory for this ResultSet
  Example    : my $dbfile_data_dir = $self->dbfile_data_dir;
  Description: Getter/Setter for the root dbfile data directory for this ResultSet
  Returntype : String
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer and new
  Status     : at risk

=cut

# currently using setter functionality in importer

sub dbfile_data_dir{
  my ($self, $data_dir) = @_;
  $self->{dbfile_data_dir} = $data_dir if defined $data_dir;
  return $self->{dbfile_data_dir};
}


=head2 result_feature_set

  Example    : if($rset->result_feature_set){ ...use result_feature table ...};
  Description: Getter for the result_feature_set attribute.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub result_feature_set{ return $_[0]->{result_feature_set}; }


=head2 table_name

  Example    : if($result_set->table_name eq 'input_set'){
                 #Do something here
               }
  Description: Getter for the table_name of a ResultSet.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub table_name{ return $_[0]->{table_name}; }


=head2 _add_table_id

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

sub _add_table_id {
  my ($self, $table_id, $cc_id) = @_;

  if (! defined $table_id){
    throw('Need to pass a table_id');
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

sub table_ids { return [ keys %{$_[0]->{'table_id_hash'}} ]; }


=head2 result_set_input_ids

  Example    : my @rset_rsi_ids = @{$result_set->result_set_input_ids()};
  Description: Getter for the input ids for this ResultSet.
  Returntype : arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub result_set_input_ids { return [ values %{$_[0]->{'table_id_hash'}} ]; }


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
  my $table_name = $chip_channel->adaptor->_main_table->[0];

  if($table_name ne $self->table_name){
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


=head2 get_support

  Example    : my @support = @{$result_set->get_support};
  Description: Retrieves a list of objects supporting the ResultSet
  Returntype : Arrayref of objects e.g. InputSet, ExperimentalChip or Channel
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub get_support {
  my $self = $_[0];
  
  if(! defined $self->{support}){
    my $adaptor_method = 'get_'.
                         join('',  (map ucfirst($_), split(/_/, $self->table_name))). #export from EFGUtils?
                         'Adaptor'; 
                          
    #Adaptor may be absent if we have an unstored object and not used -support in new
    #This will disappear once -table_id and -table_name are removed from new               
    my $supp_adaptor = $self->adaptor->db->$adaptor_method;  
  
    foreach my $id(@{$self->table_ids()}){
      push @{$self->{support}}, $supp_adaptor->fetch_by_dbID($id);
    }   
  }
    
  return $self->{support};
}




=head2 get_ExperimentalChips

  Example    : my @ecs = @{$result_set->get_ExperimentalChips()};
  Description: Retrieves a list of ExperimentalChip objects
  Returntype : Listref of ExperimentalChip object
  Exceptions : Warns if ResultSet is an InputSet
  Caller     : General
  Status     : At Risk

=cut

sub get_ExperimentalChips{
  my $self = shift;

  my $echips;

  if( $self->table_name ne 'experimental_chip' ||
      $self->table_name ne 'channel' ){
	warn 'Cannot get_ExperimentalChips for an ResultSet with table_name '.$self->table_name;
  }
  elsif($self->table_name eq 'experimental_chip'){
    $echips = $self->get_support;
  }
  else{  #table_name is channel 
    
    if(! defined $self->{experimental_chips}){
      my %echips;
      my $chan_adaptor = $self->adaptor->db->get_ChannelAdaptor;
      my $ec_adaptor   = $self->adaptor->db->get_ExperimentalChipAdaptor;
  
      foreach my $chan_id( @{$self->table_ids} ){
        my $chan = $chan_adaptor->fetch_by_dbID($chan_id);
        $echips{$chan->experimental_chip_id} ||= $ec_adaptor->fetch_by_dbID($chan->experimental_chip_id);
      }
      
      @{$self->{'experimental_chips'}} = values %echips;
    }
    
    $echips = $self->{'experimental_chips'};
  }

  return $echips;
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

  if(defined $self->feature_type){
    $label = $self->feature_type->name.":";
  }else{
    $label = "Unknown FeatureType:";
  }

  if(defined $self->cell_type){
    $label .= $self->cell_type->name;
  }else{
    $label .= "Unknown CellType";
  }

  return $self->name.":".$self->analysis->logic_name.":".$label;
}


=head2 compare_to

  Args[1]    : Bio::EnsEMBL::Funcgen::ResultSet (mandatory)
  Args[2]    : Boolean - Shallow flag, omits nested object comparisons which require
                dbID and is_stored checks.
               i.e. assumes all primary key components are the same.
  Example    : my %diffs = %{$rset->compare_to($other_rset, 1)};
  Description: Compare this ResultSet to another. Does not compare support
  Returntype : Hashref of key attribute name keys and value which differ.
  Exceptions : Throws if arg is not a valid ResultSet
  Caller     : General
  Status     : At Risk

=cut

#%diffs
#keys define the attribute/method/test
#If the key is a string it is a simple warning
#If it is an array ref, it shows the differences between the attributes tested
#if it a hash ref, it is a nest %diffs has from a nested object
#identity of this ResultSet handled in caller, not in diffs hash.

#self (or rset) would have to be stored to work with DBAdaptor check in shallow_compare_Storables


#slightly odd as setting shallow omits calling the shallow_compare_Storables method

sub compare_to {
  my ($self, $rset, $diffs, $shallow) = @_;
    
  if(! (defined $rset &&
        ref($rset) &&
        $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet')) ){
      throw('You must pass a valid Bio::EnsEMBL::Funcgen::ResultSet to compare');
  }
  
  if(defined $diffs &&
    (ref($diffs) ne 'HASH') ){
    throw('Diffs hash mush be passed as Hashref');  
  }
  else{  
    $diffs = {};
  }
  
  $self->compare_string_methods($rset, [ qw(name table_name feature_class get_all_states) ], $diffs);
  
  #We know table_ids are the same, but are they from the same db?
  #Test InputSets from one with DBAdaptor::is_stored from the other
  #InputSets have to be stored otherwise we would have thrown in new or add_table_id
  #Is this the same for other inputs/support i.e. ExperimentalChips

  #Now deal with support
  #this could be subbed out to compare_arrays 
  #args would be string/numerical and sort flag
  my @support       = sort {$a->dbID <=> $b->dbID} @{$self->get_support};
  my @other_support = sort {$a->dbID <=> $b->dbID} @{$rset->get_support};
        
  if(scalar(@support) != scalar(@other_support)){
    $diffs->{'ResultSet::get_support - size'} = [join(',', map ($_->dbID, @support)), 
                                                 join(',', map ($_->dbID, @other_support))];
  }
   
  if(! $shallow){
  
    if(scalar(@support) == scalar(@other_support)){
     
      for my $i(0..$#support){
        my %support_diffs = %{$self->shallow_stored_Storables($support[$i], $other_support[$i])};
   
        if(keys %support_diffs){
          $diffs->{'ResultSet - shallow_compare support'} = \%support_diffs;
          last; #Bailing out here may miss further mismatches
        }
      } 
    }
    
    foreach my $obj_method( qw(feature_type cell_type analysis) ){
      $diffs->{'ResultSet::'.$obj_method} = $self->compare_stored_Storables($self->obj_method, $rset->$obj_method);
    }   
  }
    
  return $diffs;
}

  
### DEPRECATED ###

sub get_InputSets{ #DEPRECATED IN v72
  deprecate('get_InputSets is now deprecated, please use get_support.');
  return $_[0]->get_support;
}


sub add_table_id { #DEPRECATED IN v72
   deprecate('The add_table_id method is no deprecated, please use the -support param in new');
   return $_[0]->_add_table_id($_[1]);
}

1;
