=head1 Bio::EnsEMBL::Funcgen::Hive::DefineResultSets
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineResultSets;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects validate_checksum);
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( merge_bams );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  $self->init_branching_by_analysis;
  
  my $input_subset_ids = $self->get_param_method('input_subset_ids',  'required');
  assert_ref($input_subset_ids, 'HASH', 'InputSubset dbIDs');
  
  $self->set_param_method('controls', delete $input_subset_ids->{controls});
  
  if(scalar(keys %$input_subset_ids) == 0){
    throw('No (signal) input_subset_ids defined. Must pass input_subset_ids hash of (optional) \'controls\''.
      ' and ResultSet name keys and input_subset_id arrayref values');
  }

  $self->get_param_method('alignment_analysis', 'required');

  # Undef in analysis DefineResultSets
  # 1 in DefineMergedReplicateResultSet
  #
  # This is how this module knows where it is in the ersa pipeline and 
  # changes its behaviour accordingly.
  #
  my $merge_idr_replicates = $self->get_param_method('merge_idr_replicates', 'silent');

  return;
}

sub _are_controls {
  my $control_input_subsets        = shift; 
  my $all_controls = 1;
  
  foreach my $ctrl(@$control_input_subsets){
    
    if(! $ctrl->is_control){
      $all_controls = 0;
      last;
    }
  }
  
  return $all_controls;
}


sub _are_signals {
  my $signal_input_subsets     = shift; 
  my $all_signal_input_subsets = 1;
  
  foreach my $sig(@$signal_input_subsets) {
    
    if($sig->is_control){
      $all_signal_input_subsets = 0;
      last;
    }
  }
  
  return $all_signal_input_subsets;
} 


sub run {

  my $self                  = shift;
  my $helper                = $self->helper;
  my $result_set_adaptor    = $self->out_db->get_ResultSetAdaptor;
  my $input_subset_ids      = $self->input_subset_ids;
  my $control_input_subsets = [];
  
  my $alignment_analysis = $self->alignment_analysis;
  my $alignment_analysis_object = &scalars_to_objects(
    $self->out_db, 'Analysis', 'fetch_by_logic_name', [$alignment_analysis]
  )->[0];

  my $controls = $self->controls;

  # This is run in the DefineResultSets analysis
  #
  if ( $controls && assert_ref($controls, 'ARRAY', 'controls') && @$controls) {
  
    $control_input_subsets = &scalars_to_objects($self->out_db, 'InputSubset',
                                                'fetch_by_dbID',
                                                $controls);
    
    if(! &_are_controls($control_input_subsets)){
      throw("Found unexpected non-control InputSubsets specified as controls\n\t".
        join("\n\t", map($_->name, @$control_input_subsets)));
    }
    
    my %exps;
    
    foreach my $ctrl(@$control_input_subsets){
      $exps{$ctrl->experiment->name} = $ctrl->experiment;
    }
    
    if( scalar(keys(%exps)) != 1 ){
      throw("Failed to identify a unique control Experiment for :\n".
        join(' ', keys(%exps)));  
    }
    my ($exp) = values(%exps);#We only have one
  }
 
  my (%result_sets, %replicate_bam_files);

  # merge_idr_replicates is undef in DefineResultSets analysis
  #
  my $merge_idr_replicates = $self->merge_idr_replicates;
  my $foo_bar = 'thisisdeprecatedandshouldntappearanywhereifitdoespleaseletmeknow';
  
  # Looks like this:
  #
  # input_subset_ids => {'K562:hist:BR1_H3K27me3_3526' => [3219,3245,3429],'controls' => [3458]}
  #
  foreach my $set_name (keys %$input_subset_ids) {
  
    # $set_name is 'K562:hist:BR1_H3K27me3_3526'
    #
    # $parent_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
    #
    my $parent_set_name = $set_name.'_'.$alignment_analysis_object->logic_name;
    
    # $signal_input_subsets is set to the input sub set objects of
    #
    # [3219,3245,3429]
    #
    my $signal_input_subsets = &scalars_to_objects(
      $self->out_db, 'InputSubset', 'fetch_by_dbID', $input_subset_ids->{$set_name}
    );
    
    if (! &_are_signals($signal_input_subsets)) {
      throw("Found unexpected controls in signal InputSubsets\n\t".
      join("\n\t", map($_->name, @$signal_input_subsets)));
    }

    # The object representing H3K27me3
    #
    my $feature_type        = $signal_input_subsets->[0]->feature_type;
    
    # False, because H3K27me3 is a broad feature type.
    #
    my $is_idr_feature_type = $self->is_idr_FeatureType($feature_type);
    
    # Cell type object for 'K562:hist:BR1'
    #
    my $cell_type           = $signal_input_subsets->[0]->cell_type;
    
    #Define a single rep set with all of the InputSubsets 
    #i.e. non-IDR merged or post-IDR merged     
    
    # This is an array with one element. The one element is an array 
    # reference to $signal_input_subsets.
    #
    my @rep_sets = ($signal_input_subsets);

    # Evaluates to 1 for our example dataset.
    #
    my $has_signal_replicates = (scalar(@$signal_input_subsets) > 1) ? 1 : 0;
    
    my $current_set_to_be_processed_with_idr = $is_idr_feature_type && $has_signal_replicates;

    if($current_set_to_be_processed_with_idr) {
        # Split single rep sets
        @rep_sets = map {[$_]} @$signal_input_subsets;  
    }

    foreach my $rep_set (@rep_sets) {
    
      my $result_set_name;
      if($current_set_to_be_processed_with_idr) {
      
        # There will be only 1 in the $rep_set
        $result_set_name = $parent_set_name . '_TR'.$rep_set->[0]->replicate;

      } else {
	# Reminder: $parent_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
	#
	$result_set_name = $parent_set_name;
      }

      my $result_set_constructor_parameters = {
	-RESULT_SET_NAME     => $result_set_name,
	-SUPPORTING_SETS     => [@$rep_set, @$control_input_subsets],
	-DBADAPTOR           => $self->out_db,
	-RESULT_SET_ANALYSIS => $alignment_analysis_object,
	-ROLLBACK            => $self->param_silent('rollback'),
	-RECOVER             => $self->param_silent('recover'),
	-FULL_DELETE         => $self->param_silent('full_delete'),
	-CELL_TYPE           => $cell_type,
	-FEATURE_TYPE        => $feature_type
      };

      if ($current_set_to_be_processed_with_idr) {

	# This goes to the branch with replicate signals
      
	$result_sets{$foo_bar}->{$parent_set_name} ||= [];  
	push @{$result_sets{$foo_bar}->{$parent_set_name}}, $result_set_constructor_parameters
	
      } else {
	
	# This goes to the branch for merged signals
	
	$result_sets{$foo_bar}{merged} ||= [];
	push @{$result_sets{$foo_bar}{merged}}, $result_set_constructor_parameters
      }

    }
  }
  
  my %batch_params = %{$self->batch_params};
  my %branch_sets;
  my $tracking_adaptor = $self->tracking_adaptor;

  foreach my $result_set_group_name (keys %{$result_sets{$foo_bar}}) {
  
    my $result_set_group = $result_sets{$foo_bar}->{$result_set_group_name};

    foreach my $result_set (@$result_set_group) {
    
      $result_set = $helper->define_ResultSet(%{$result_set});

      $self->helper->debug(1, "Caching $foo_bar branch jobs for ".$result_set->name);

      $branch_sets{$foo_bar}{$result_set_group_name}{set_names} ||= [];
      $branch_sets{$foo_bar}{$result_set_group_name}{dbIDs}     ||= [];
      push @{$branch_sets{$foo_bar}{$result_set_group_name}{set_names}}, $result_set->name;
      push @{$branch_sets{$foo_bar}{$result_set_group_name}{dbIDs}},     $result_set->dbID;
    }
  }
  
  $self->helper->debug(1, "Processing cached branch $foo_bar");

  # Pick an arbitrary set for access to the controls 
  my ($any_group) = keys(%{$branch_sets{$foo_bar}});

  # result_set_groups has all the information necessary for 
  # JobFactorySignalProcessing to create jobs for running the bwa 
  # analysis triplet on the signals.
  #
  $self->branch_job_group(
    2,
    [
      {
	%batch_params,
	result_set_groups => $branch_sets{$foo_bar},
	set_type  => 'ResultSet',
	dbID      => $branch_sets{$foo_bar}{$any_group}{dbIDs}->[0],
	set_name  => $branch_sets{$foo_bar}{$any_group}{set_names}->[0]
      }
    ]
  );
  
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;