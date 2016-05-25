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
use Data::Dumper;

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  #$self->init_branching_by_analysis;
  
  my $input_subset_ids = $self->get_param_method('input_subset_ids',  'required');
  assert_ref($input_subset_ids, 'HASH', 'InputSubset dbIDs');
  
  $self->set_param_method('controls', delete $input_subset_ids->{controls});
  
  if(scalar(keys %$input_subset_ids) == 0){
    throw('No (signal) input_subset_ids defined. Must pass input_subset_ids hash of (optional) \'controls\''.
      ' and ResultSet name keys and input_subset_id arrayref values');
  }

  $self->get_param_method('alignment_analysis', 'required');

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

sub _throw_if_control_parameters_are_invalid {
  my $self = shift;
  my $controls = shift;
  
  my $control_parameters_are_valid = $controls && assert_ref($controls, 'ARRAY', 'controls') && @$controls;
  
  if ( !$control_parameters_are_valid ) {
    throw("Invalid controls!");
  }
}

sub _throw_if_control_input_subsets_are_invalid {
  my $self = shift;
  my $control_input_subsets = shift;
  if(! &_are_controls($control_input_subsets)){
    throw("Found unexpected non-control InputSubsets specified as controls\n\t".
      join("\n\t", map($_->name, @$control_input_subsets)));
  }
  
  my $control_experiments = $self->_fetch_control_experiments($control_input_subsets);
  
  if( scalar(keys(%$control_experiments)) != 1 ) {
    throw("Failed to identify a unique control Experiment for :\n".
      join(' ', keys(%$control_experiments)));
  }
}

sub _fetch_control_experiment {
  my $self = shift;
  my $control_input_subsets = shift;
  
  my $control_experiments = $self->_fetch_control_experiments($control_input_subsets);
  
  # There should only ever be one.
  #
  my ($exp) = values(%$control_experiments);
  return $exp;
}

sub _fetch_control_experiments {
  my $self = shift;
  my $control_input_subsets = shift;
  
  my %exps;
  foreach my $ctrl(@$control_input_subsets) {
    $exps{$ctrl->experiment->name} = $ctrl->experiment;
  }
  return \%exps;
}

sub run {

  my $self                  = shift;
  my $helper                = $self->helper;
  #my $result_set_adaptor    = $self->out_db->get_ResultSetAdaptor;
  my $input_subset_ids      = $self->input_subset_ids;
  my $control_input_subsets = [];
  
  my $alignment_analysis = $self->alignment_analysis;
#   my $alignment_analysis = 'alignment_analysis';
  
  my $alignment_analysis_object = &scalars_to_objects(
    $self->out_db, 'Analysis', 'fetch_by_logic_name', [$alignment_analysis]
  )->[0];

  my $controls = $self->controls;
  
  $self->_throw_if_control_parameters_are_invalid($controls);
  
  $control_input_subsets = &scalars_to_objects(
    $self->out_db, 
    'InputSubset',
    'fetch_by_dbID',
    $controls
  );
  
  $self->_throw_if_control_input_subsets_are_invalid($control_input_subsets);
  
# ------------------------------------------------------------------------------

  my $control_experiment = $self->_fetch_control_experiment($control_input_subsets);
  
  my $result_set_name = $control_experiment->name . '_control_alignment';

  my $control_result_set_constructor_parameters = {
    -RESULT_SET_NAME      => $result_set_name,
#     -RESULT_SET_REPLICATE => undef,
    -SUPPORTING_SETS      => [@$control_input_subsets],
    -DBADAPTOR            => $self->out_db,
    -RESULT_SET_ANALYSIS  => $alignment_analysis_object,
    -ROLLBACK             => $self->param_silent('rollback'),
    -RECOVER              => $self->param_silent('recover'),
    -FULL_DELETE          => $self->param_silent('full_delete'),
    -EPIGENOME            => $control_experiment->epigenome,
    -FEATURE_TYPE         => $control_experiment->feature_type
  };
  
  my $control_result_set = $helper->define_ResultSet(%$control_result_set_constructor_parameters);

# ------------------------------------------------------------------------------
  
  my (%result_sets, %replicate_bam_files);

  # Looks like this:
  #
  # input_subset_ids => {'K562:hist:BR1_H3K27me3_3526' => [3219,3245,3429],'controls' => [3458]}
  #
  foreach my $experiment_name (keys %$input_subset_ids) {
  
    # $experiment_name is 'K562:hist:BR1_H3K27me3_3526'
    #
    # $parent_result_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
    #
    my $parent_result_set_name = $experiment_name . '_' . $alignment_analysis_object->logic_name;
    
    # $signal_input_subsets is set to the input sub set objects of
    #
    # [3219,3245,3429]
    #
    my $signal_input_subsets = &scalars_to_objects(
      $self->out_db, 'InputSubset', 'fetch_by_dbID', $input_subset_ids->{$experiment_name}
    );
    
    if (! &_are_signals($signal_input_subsets)) {
      throw("Found unexpected controls in signal InputSubsets\n\t".
      join("\n\t", map($_->name, @$signal_input_subsets)));
    }

    (
      my $key, 
      my $constructor_parameters
    ) = $self->create_result_set_constructor_parameters(
      {
	signal_input_subsets      => $signal_input_subsets,
	control_input_subsets     => $control_input_subsets,
	alignment_analysis_object => $alignment_analysis_object,
	parent_result_set_name    => $parent_result_set_name,
      }
    );
    $result_sets{$key} = $constructor_parameters;
  }
  my %batch_params = %{$self->batch_params};
  my %branch_sets;

  $Data::Dumper::Maxdepth = 3;
  print Dumper(\%result_sets);
  die();
  
  my @hive_jobs_fix_experiment_id;
  foreach my $result_set_group_name (keys %result_sets) {
  
    my $result_set_group = $result_sets{$result_set_group_name};

    foreach my $current_result_set_constructor_parameter (@$result_set_group) {
    
      my $result_set = $helper->define_ResultSet(%{$current_result_set_constructor_parameter});
      
      push @hive_jobs_fix_experiment_id, { 
	dbID      => $result_set->dbID, 
	#technical_replicate => $current_result_set_constructor_parameter->{-RESULT_SET_REPLICATE},
      };

      $self->helper->debug(1, "Caching branch jobs for ".$result_set->name);

      $branch_sets{$result_set_group_name}{set_names} ||= [];
      $branch_sets{$result_set_group_name}{dbIDs}     ||= [];
      push @{$branch_sets{$result_set_group_name}{set_names}}, $result_set->name;
      push @{$branch_sets{$result_set_group_name}{dbIDs}},     $result_set->dbID;
    }
  }
  
  push @hive_jobs_fix_experiment_id, { 
    dbID      => $control_result_set->dbID,
  };

  $self->helper->debug(1, "Processing cached branch");

  # result_set_groups has all the information necessary for 
  # JobFactorySignalProcessing to create jobs for running the bwa 
  # analysis triplet on the signals.
  #
  $self->branch_job_group(
    2,
    \@hive_jobs_fix_experiment_id,
    3,
    [
      {
	%batch_params,
	result_set_groups => \%branch_sets,
	set_type  => 'ResultSet',
	dbID      => $control_result_set->dbID,
	set_name  => $control_result_set->name,
      }
    ]
  );
  
  return;
}

sub create_result_set_constructor_parameters_with_idr {
  my $self  = shift;
  my $param = shift;
  
  my $signal_input_subsets      = $param->{signal_input_subsets};
  my $control_input_subsets     = $param->{control_input_subsets};
  my $alignment_analysis_object = $param->{alignment_analysis_object};
  my $parent_result_set_name    = $param->{parent_result_set_name};

  my @constructor_parameters;
  foreach my $replicate_input_subset (@$signal_input_subsets) {

      my $technical_replicate_number = $replicate_input_subset->technical_replicate;
      my $result_set_name = $parent_result_set_name . '_TR' . $technical_replicate_number;
  
      my $result_set_constructor_parameters = {
	-RESULT_SET_NAME      => $result_set_name,
	-SUPPORTING_SETS      => [$replicate_input_subset, @$control_input_subsets],
	-DBADAPTOR            => $self->out_db,
	-RESULT_SET_ANALYSIS  => $alignment_analysis_object,
	-ROLLBACK             => $self->param_silent('rollback'),
	-RECOVER              => $self->param_silent('recover'),
	-FULL_DELETE          => $self->param_silent('full_delete'),
	-EPIGENOME            => $replicate_input_subset->epigenome,
	-FEATURE_TYPE         => $replicate_input_subset->feature_type
      };
      push @constructor_parameters, $result_set_constructor_parameters;
    }
    return ($parent_result_set_name, \@constructor_parameters);
}

sub create_result_set_constructor_parameters_merged {

  my $self  = shift;
  my $param = shift;
  
  my $signal_input_subsets      = $param->{signal_input_subsets};
  my $control_input_subsets     = $param->{control_input_subsets};
  my $alignment_analysis_object = $param->{alignment_analysis_object};
  my $parent_result_set_name    = $param->{parent_result_set_name};

  my $result_set_constructor_parameters = {
    -RESULT_SET_NAME      => $parent_result_set_name,
    -SUPPORTING_SETS      => [@$signal_input_subsets, @$control_input_subsets],
    -DBADAPTOR            => $self->out_db,
    -RESULT_SET_ANALYSIS  => $alignment_analysis_object,
    -ROLLBACK             => $self->param_silent('rollback'),
    -RECOVER              => $self->param_silent('recover'),
    -FULL_DELETE          => $self->param_silent('full_delete'),
    -EPIGENOME            => $signal_input_subsets->[0]->epigenome,
    -FEATURE_TYPE         => $signal_input_subsets->[0]->feature_type
  };
  return ('merged', [ $result_set_constructor_parameters ]);
}

sub create_result_set_constructor_parameters {

  my $self  = shift;
  my $param = shift;
  
  my $signal_input_subsets = $param->{signal_input_subsets};
  my $has_signal_replicates = (scalar(@$signal_input_subsets) > 1) ? 1 : 0;
  my $is_idr_feature_type = $self->is_idr_FeatureType(
    $signal_input_subsets->[0]->feature_type
  );

  my $current_set_to_be_processed_with_idr = $is_idr_feature_type && $has_signal_replicates;

  if($current_set_to_be_processed_with_idr) {
    return $self->create_result_set_constructor_parameters_with_idr($param);
  }
  return $self->create_result_set_constructor_parameters_merged($param);
}


sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;