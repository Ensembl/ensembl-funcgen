
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineResultSets

=head1 DESCRIPTION

The modules creates and stores ResultSet objects given a list of input_subset_ids. Dynamic data flow is performed 
via the 'branching_by_analysis' mechanism which allows branch names to be used to data flow. This allows dataflow 
to be defined by the inputs of an analysis wrt the pre-defined analysis config.

It functions in two modes:

  1 DefineResultSets analysis
  This runs prior to any alignment analyses. For experiments which are destined for IDR analysis, 
  it creates individual replicate ResultSets, else it creates a merged replicate ResultSet. Branching by
  analysis data flows to:
    Preprocess_${alignment_analysis}_(control|merged|replicate) & PreprocessIDR (dependant on 'replciate')
    
  2 DefineMergedReplicateResultSet analysis
  This runs post IDR analysis and also merges the replicate alignments. This always dataflows to DefineMergedDataSet

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
  #
  my $merge_idr_replicates = $self->get_param_method('merge_idr_replicates', 'silent');
  
  print "\n\n--------> $merge_idr_replicates \n\n";

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
  my $result_set_adaptor    = $self->out_db->get_ResultSetAdaptor;
  my $input_subset_ids      = $self->input_subset_ids;
  my $control_input_subsets = [];
  my $branch;
  
  my $alignment_analysis = $self->alignment_analysis;
  my $alignment_analysis_object = &scalars_to_objects(
    $self->out_db, 'Analysis', 'fetch_by_logic_name', [$alignment_analysis]
  )->[0];

  my $control_branch = 'Preprocess_'.$alignment_analysis.'_control';
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
  
  print "\n----- merge_idr_replicates ----------------------------------------\n";
  print Dumper($merge_idr_replicates);
  print "\n-------------------------------------------------------------\n";
   
  foreach my $set_name(keys %$input_subset_ids) {
    my $parent_set_name = $set_name.'_'.$alignment_analysis_object->logic_name;
    
    my $signal_input_subsets = &scalars_to_objects($self->out_db, 'InputSubset',
                                                  'fetch_by_dbID',
                                                  $input_subset_ids->{$set_name});
    if(! &_are_signals($signal_input_subsets)){
      throw("Found unexpected controls in signal InputSubsets\n\t".
      join("\n\t", map($_->name, @$signal_input_subsets)));
    }

    my $feature_type        = $signal_input_subsets->[0]->feature_type;
    my $is_idr_feature_type = $self->is_idr_FeatureType($feature_type);
    my $cell_type           = $signal_input_subsets->[0]->cell_type;
    
    #Define a single rep set with all of the InputSubsets 
    #i.e. non-IDR merged or post-IDR merged     
    my @rep_sets = ($signal_input_subsets);
    my $has_signal_replicates = (scalar(@$signal_input_subsets) >1) ? 1 : 0;    

    # Only run in DefineResultSets analysis
    #
    if(!$is_idr_feature_type || !$has_signal_replicates) {
      # single rep ID or multi-rep non-IDR
      $branch = 'Preprocess_'.$alignment_analysis.'_merged';
    }
    
    if($merge_idr_replicates && $is_idr_feature_type && $has_signal_replicates) {
    
        $branch = 'DefineMergedDataSet';

        $replicate_bam_files{$parent_set_name}{rep_bams} = [];

        foreach my $rep (@$signal_input_subsets) {
          #Pull back the rep Rset to validate and get the alignment file for merging
          my $rep_result_set_name = $parent_set_name.'_TR'.$rep->replicate;
          my $result_set = $result_set_adaptor->fetch_by_name($rep_result_set_name);            
          #Could also fetch them with $result_set_a->fetch_all_by_supporting_Sets($rep).
           
          if(! defined $result_set){
            $self->throw_no_retry("Could not find ResultSet for post-IDR merged ResultSet:\t".
              $rep_result_set_name); 
          }
      
          #todo validate controls are the same
          #This should already have been done in PreprocessIDR
          #but probably a good idea to do here too
          #As we may get here by means other than PreprocessIDR?
                    
          push @{$replicate_bam_files{$parent_set_name}{rep_bams}}, 
            $self->get_alignment_files_by_ResultSet_formats($result_set, ['bam'])->{bam};
        }
    }
    
    if(! $merge_idr_replicates && $is_idr_feature_type && $has_signal_replicates) {
        #RunIDR semaphore handled implicitly later 
        $branch = 'Preprocess_'.$alignment_analysis.'_replicate';
        
        # Split single rep sets
        @rep_sets = map {[$_]} @$signal_input_subsets;  
    }

    # Don't use control branch, if merge_idr_replicates is set. In this case this 
    # module is being used in the DefineMergedReplicateResultSet analysis. 
    # This module also backs the DefineResultSets analysis. In the future, 
    # this module should be split into two modules for each of the two 
    # analyses. Until then, setting merge_idr_replicates serves to tell the module, 
    # which analysis it is currently running.
    # 
    if (!$merge_idr_replicates) {
      $branch = $control_branch;
    }
    
    
    foreach my $rep_set (@rep_sets) {
      my $result_set_name = $parent_set_name;#.'_'.$alignment_analysis_object->logic_name;
    
      if($is_idr_feature_type && $has_signal_replicates && ! $merge_idr_replicates){
        #there will be only 1 in the $rep_set 
        $result_set_name .= '_TR'.$rep_set->[0]->replicate;
      }

      my $cache_ref;

      if(! $merge_idr_replicates) {
	# branch can be replicate(no control) or control(with reps)  
        $result_sets{$branch}->{$parent_set_name} ||= [];  
        $cache_ref                            = $result_sets{$branch}->{$parent_set_name};
      }

      if($merge_idr_replicates) {

        #branches can be
        # merged (no controls)
        # control (single rep IDR set or merged) 
        # or DefineMergedDataSet i.e. this is the IDR analysis is DefineMergedReplicateResultSet
        #is used of merged key here correct for DefineMergedDataSet?

        $result_sets{$branch}{merged} ||= [];
        $cache_ref                = $result_sets{$branch}{merged};
      }

      push @$cache_ref, {
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
    }
  }
  
  use Data::Dumper;
  $Data::Dumper::Maxdepth = 3;
  print Dumper(\%result_sets);
  
  
  #Now do the actual ResultSet generation and cache the output_ids on the correct branch
  my %batch_params = %{$self->batch_params};
  my %branch_sets;
  my $tracking_adaptor = $self->tracking_adaptor;

  #We need to test for CONTROL_ALIGNED here
  #but this is currently only set on the ResultSets
  #and these maybe entirely new result sets
  #so this has to be set on the experiment too!
  #we can test that above and set the branch accordingly

  foreach my $branch(keys %result_sets){
  
    foreach my $result_set_group_name (keys %{$result_sets{$branch}}){
      my $result_set_group = $result_sets{$branch}->{$result_set_group_name};

      foreach my $result_set (@$result_set_group) {
        $result_set = $helper->define_ResultSet(%{$result_set});

        if($merge_idr_replicates) {

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
        }

        if($branch =~ /(_merged$|^DefineMergedDataSet$)/){# (no control) job will only ever have 1 result_set
          $self->branch_job_group($branch, [{%batch_params,
                                             dbID       => $result_set->dbID, 
                                             set_name   => $result_set->name,
                                             set_type    => 'ResultSet',
                                             }]);
        }
        elsif($branch =~ /(_control$|_replicate$)/){
          
          $self->helper->debug(1, "Cacheing $branch branch jobs for ".$result_set->name);
          
          $branch_sets{$branch}{$result_set_group_name}{set_names} ||= [];
          $branch_sets{$branch}{$result_set_group_name}{dbIDs}     ||= [];
          push @{$branch_sets{$branch}{$result_set_group_name}{set_names}}, $result_set->name;
          push @{$branch_sets{$branch}{$result_set_group_name}{dbIDs}},     $result_set->dbID;
                    
          if($branch =~ /_replicate$/o){
            #Just create the individual fan output_ids
            $branch_sets{$branch}{$result_set_group_name}{output_ids} ||= [];
            push @{$branch_sets{$branch}{$result_set_group_name}{output_ids}}, {%batch_params,
                                                                          dbID     => $result_set->dbID,
                                                                          set_name => $result_set->name,
                                                                          set_type => 'ResultSet'};    
          }
        }
        else{ #Sanity check
          #branch is always defined within this module
          $self->throw_no_retry("$branch is not supported by DefineResultSets");  
        }
      }
    }
  }
  
  
  #Now flow job_groups of result_sets to control job/replicate & IDR 
  foreach my $branch(keys %branch_sets){            
    $self->helper->debug(1, "Processing cached branch $branch");
    #what about other branches in here
  
    if($branch =~ /_replicate$/){
      
      foreach my $rep_set(keys %{$branch_sets{$branch}}){ 
        #Add semaphore RunIDR job
        $self->branch_job_group($branch, $branch_sets{$branch}{$rep_set}{output_ids}, 'PreprocessIDR', 
                                [{%batch_params,
                                 dbIDs     => $branch_sets{$branch}{$rep_set}{dbIDs},
                                 set_names => $branch_sets{$branch}{$rep_set}{set_names},
                                 set_type  => 'ResultSet'}]);
      }           
    }       
    elsif($branch =~ /_control$/){    
      #Pick an arbitrary set for access to the controls 
      my ($any_group) = keys(%{$branch_sets{$branch}});
   
      #result_set_groups here is used by MergeControlAlignments_and_QC to flow correctly
      $self->branch_job_group($branch, 
                              [{%batch_params,
                               result_set_groups => $branch_sets{$branch},
                               set_type  => 'ResultSet',
                               dbID      => $branch_sets{$branch}{$any_group}{dbIDs}->[0],
                               set_name  => $branch_sets{$branch}{$any_group}{set_names}->[0]}]);   
    }
  }
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;