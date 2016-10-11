
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::MergeIDRReplicateAlignments

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeIDRReplicateAlignments;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects validate_checksum);
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( merge_bams );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );


sub fetch_input {
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

  # Undef in analysis MergeIDRReplicateAlignments
  # 1 in MergeIDRReplicateAlignments
  #
  # This is how this module knows where it is in the ersa pipeline and 
  # changes its behaviour accordingly.
  #
  my $merge_idr_replicates = $self->get_param_method('merge_idr_replicates', 'silent');

  if($merge_idr_replicates) {
    
    if( scalar(keys %$input_subset_ids) != 1 ) {
      $self->throw_no_retry('Cannot currently specify > 1 input_subsets_ids group in merge_idr_replicate');
    }

    # Dataflowed from PostprocessIDR
    $self->get_param_method('max_peaks', 'required'); 

    # Batch flown
    my $permissive_peak_logic_name = $self->param_required('permissive_peaks');

    $self->set_param_method(
      'idr_peak_analysis_id',
      scalars_to_objects($self->out_db, 'Analysis', 'fetch_by_logic_name', $permissive_peak_logic_name)->[0]->dbID
    );
  }
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
  
  # HACK
  my $db = $self->param_required('out_db');
  my $result_set_adaptor = $db->get_ResultSetAdaptor;
  $result_set_adaptor->{file_type} = 'BAM';
  
#   my $result_set_adaptor    = $self->out_db->get_ResultSetAdaptor;
  my $input_subset_ids      = $self->input_subset_ids;
  my $branch;
  
  my $alignment_analysis = $self->alignment_analysis;
  my $alignment_analysis_object = &scalars_to_objects(
    $self->out_db, 'Analysis', 'fetch_by_logic_name', [$alignment_analysis]
  )->[0];

  my $control_branch = 'Preprocess_'.$alignment_analysis.'_control';
  my $controls = $self->controls;
  my $control_input_subsets = [];
  
  if(
    $controls 
    && assert_ref($controls, 'ARRAY', 'controls') 
    && @$controls) {
    $control_input_subsets = &scalars_to_objects(
      $self->out_db, 'InputSubset', 'fetch_by_dbID', $controls
    );
    if(! &_are_controls($control_input_subsets)) {
      throw("Found unexpected non-control InputSubsets specified as controls\n\t"
      . join("\n\t", map($_->name, @$control_input_subsets)));
    }
  }

#   use Data::Dumper;
#   print Dumper($control_input_subsets);
#   die();
  
  my (%result_sets, %replicate_bam_files);

  # merge_idr_replicates is undef in MergeIDRReplicateAlignments analysis
  #
  my $merge_idr_replicates = $self->merge_idr_replicates;
  
  # Looks like this:
  #
  # input_subset_ids => {'K562:hist:BR1_H3K27me3_3526' => [3219,3245,3429],'controls' => [3458]}
  #
  foreach my $set_name (keys %$input_subset_ids) {
  
#     die;
  
    # $set_name is 'K562:hist:BR1_H3K27me3_3526'
    #
    # $parent_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
    #
#     my $parent_set_name = $set_name.'_'.$alignment_analysis_object->logic_name;
    
    
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
    my $parent_set_name = $signal_input_subsets->[0]->experiment->name . '_' . $alignment_analysis_object->logic_name;

    # The object representing H3K27me3
    #
    my $feature_type        = $signal_input_subsets->[0]->feature_type;
    
    # False, because H3K27me3 is a broad feature type.
    #
    my $is_idr_feature_type = $self->is_idr_FeatureType($feature_type);
    
    # Epigenome object for 'K562:hist:BR1'
    #
    my $epigenome           = $signal_input_subsets->[0]->epigenome;
    
    #Define a single rep set with all of the InputSubsets 
    #i.e. non-IDR merged or post-IDR merged     
    
    # This is an array with one element. The one element is an array 
    # reference to $signal_input_subsets.
    #
    my @rep_sets = ($signal_input_subsets);

    # Evaluates to 1 for our example dataset.
    #
    my $has_signal_replicates = (scalar(@$signal_input_subsets) > 1) ? 1 : 0;

    $branch = 'DefineMergedDataSet';

    $replicate_bam_files{$parent_set_name}{rep_bams} = [];
    
    my $biological_replicates_present = $self->signal_input_subsets_have_biological_replicates($signal_input_subsets);
    
    foreach my $current_replicate (@$signal_input_subsets) {

      my $rep_result_set_name;
      
      if ($biological_replicates_present) {
	$rep_result_set_name = $self->create_result_set_name_for_biological_replicate({
	  parent_result_set_name => $parent_set_name,
	  replicate_number       => $current_replicate->biological_replicate
	});
      } else {
	$rep_result_set_name = $self->create_result_set_name_for_technical_replicate({
	  parent_result_set_name => $parent_set_name,
	  replicate_number       => $current_replicate->technical_replicate
	});
      }
      
      my $result_set = $result_set_adaptor->fetch_by_name($rep_result_set_name);
	
      if(! defined $result_set) {
	$self->throw_no_retry(
	  "Could not find ResultSet for post-IDR merged ResultSet: "
	 . $rep_result_set_name
	);
      }
      
      my $current_bam_file = $self->db_output_dir
      . '/' . $self->get_alignment_files_by_ResultSet_formats($result_set);
      
      push @{$replicate_bam_files{$parent_set_name}{rep_bams}},	$current_bam_file;
    }

#     use Data::Dumper;
#     die(Dumper(\%replicate_bam_files));
    
    foreach my $rep_set (@rep_sets) {

      # Reminder: $parent_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
      #
      my $result_set_name = $parent_set_name;

      my $result_set_constructor_parameters = {
	-RESULT_SET_NAME     => $result_set_name,
	-SUPPORTING_SETS     => [@$rep_set, @$control_input_subsets],
	-DBADAPTOR           => $self->out_db,
	-RESULT_SET_ANALYSIS => $alignment_analysis_object,
	-ROLLBACK            => $self->param_silent('rollback'),
	-RECOVER             => $self->param_silent('recover'),
	-FULL_DELETE         => $self->param_silent('full_delete'),
	-EPIGENOME           => $epigenome,
	-FEATURE_TYPE        => $feature_type
      };

      $result_sets{$branch}{merged} ||= [];
      push @{$result_sets{$branch}{merged}}, $result_set_constructor_parameters

    }
  }
  
  #Now do the actual ResultSet generation and cache the output_ids on the correct branch
  my %batch_params = %{$self->batch_params};
  my %branch_sets;
  my $tracking_adaptor = $self->tracking_adaptor;

  foreach my $branch(keys %result_sets) {
  
    foreach my $result_set_group_name (keys %{$result_sets{$branch}}){
      my $result_set_group = $result_sets{$branch}->{$result_set_group_name};

      foreach my $result_set (@$result_set_group) {
        $result_set = $helper->define_ResultSet(%{$result_set});

	$tracking_adaptor->store_tracking_info($result_set, {
	    allow_update => 1,
	    info         => {
	      idr_max_peaks        => $self->max_peaks,
	      idr_peak_analysis_id => $self->idr_peak_analysis_id,
	    }
	  }
	);

	my $merged_file = $self->get_alignment_path_prefix_by_ResultSet($result_set).'.bam'; 

	merge_bams({
	  input_bams => $replicate_bam_files{$result_set->name}{rep_bams}, 
	  output_bam => $merged_file,
	  debug      => $self->debug,
	});
	
	$result_set->adaptor->dbfile_data_root($self->db_output_dir);
	$result_set->dbfile_path($merged_file);
	$result_set->adaptor->store_dbfile_path($result_set, 'BAM');

	$self->branch_job_group(
	  2,
	  [
	    {
	      %batch_params,
	      dbID       => $result_set->dbID, 
	      set_name   => $result_set->name,
	      set_type   => 'ResultSet',
	      bam_file   => $merged_file,
	    }
	  ]
	);

      }
    }
  }
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;