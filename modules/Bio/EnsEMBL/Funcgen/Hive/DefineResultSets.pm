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
  #my $result_set_adaptor    = $self->out_db->get_ResultSetAdaptor;
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
    
      my $replicate = $rep_set->[0]->replicate;
    
      my $result_set_name;
      if($current_set_to_be_processed_with_idr) {
      
        # There will be only 1 in the $rep_set
        $result_set_name = $parent_set_name . '_TR'.$replicate;

      } else {
	# Reminder: $parent_set_name is 'K562:hist:BR1_H3K27me3_3526_bwa_samse'
	#
	$result_set_name = $parent_set_name;
      }

      my $result_set_constructor_parameters = {
	-RESULT_SET_NAME      => $result_set_name,
	-RESULT_SET_REPLICATE => $replicate,
	-SUPPORTING_SETS      => [@$rep_set, @$control_input_subsets],
	-DBADAPTOR            => $self->out_db,
	-RESULT_SET_ANALYSIS  => $alignment_analysis_object,
	-ROLLBACK             => $self->param_silent('rollback'),
	-RECOVER              => $self->param_silent('recover'),
	-FULL_DELETE          => $self->param_silent('full_delete'),
	-CELL_TYPE            => $cell_type,
	-FEATURE_TYPE         => $feature_type
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
  
  my @hive_jobs_fix_experiment_id;
  foreach my $result_set_group_name (keys %{$result_sets{$foo_bar}}) {
  
    my $result_set_group = $result_sets{$foo_bar}->{$result_set_group_name};

    foreach my $current_result_set_constructor_parameter (@$result_set_group) {
    
=head1 A note on storing of result sets and their supporting sets

How supporting sets for result sets are stored and get retrieved:

They get stored when the result set is created like this:

```
 $result_set = $helper->define_ResultSet(%{$result_set});
```

This is the next command below.

At the end of the define_ResultSet method there is a call to 
_validate_rollback_Set like this:

```
 return $self->_validate_rollback_Set($stored_rset, $rset, 'result_set', $rollback_level,
$rset_adaptor, $slices, $recover, $full_delete,
$rset_mode);
```

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/Utils/Helper.pm#L893

This is the method in which the result set object is stored in the database. 
The result set is stored towards the end of _validate_rollback_Set, if the 
full_delete flag has been set:

```
 # REDEFINE STORED SET AFTER FULL DELETE
if($full_delete){
($stored_set) = @{ $adaptor->store($new_set) };
} 
```

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/Utils/Helper.pm#L1122

When looking at the store method used above, it might not be immediately 
obvious where the store method in the ResultSetAdaptor stores the supporting 
sets:

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/DBSQL/ResultSetAdaptor.pm#L461

The line that triggers storing of the supporting sets is this one:

```
$self->store_chip_channels($rset);
```

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/DBSQL/ResultSetAdaptor.pm#L499

Note that in the store_chip_channels method the supporting data is not stored in the 
supporting_set table as the name might suggest, but in result_set_input:


```
 my $sth = $self->prepare('INSERT INTO result_set_input (result_set_id, table_id, table_name)'.
' VALUES (?, ?, ?)');
my $sth1 = $self->prepare('INSERT INTO result_set_input '.
'(result_set_input_id, result_set_id, table_id, table_name) VALUES (?, ?, ?, ?)');
```

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/DBSQL/ResultSetAdaptor.pm#L661

The supporting set table is used for objects that are supported by result 
sets, not for storing support for a result set.

For manual inspection the supporting sets for a result set can be selected 
from the database like this:

```
select * from result_set join result_set_input using (result_set_id) join input_subset on (table_id=input_subset_id) where result_set.name = "K562:hist:BR2_H3K27ac_3526_bwa_samse"
```

Result sets missing controls can be found like this:

```
select result_set.name, sum(input_subset.is_control=1) as c from result_set 
join result_set_input using (result_set_id) 
join input_subset on (table_id=input_subset_id) 
join experiment on (result_set.experiment_id = experiment.experiment_id) 
join experimental_group using (experimental_group_id) 
where experimental_group.name = "CTTV020"
group by result_set.name
having c=0

```

Also note that when ResultSet objects are created, they don't have the 
supporting sets added to them when they are constructed:

https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/DBSQL/ResultSetAdaptor.pm#L427

They are lazy loaded when the support is requested using the table_ids.

=cut
      my $result_set = $helper->define_ResultSet(%{$current_result_set_constructor_parameter});
      
      push @hive_jobs_fix_experiment_id, { 
	dbID      => $result_set->dbID, 
	replicate => $current_result_set_constructor_parameter->{-RESULT_SET_REPLICATE},
      };

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
    \@hive_jobs_fix_experiment_id,
    3,
    [
      {
	%batch_params,
	result_set_groups => $branch_sets{$foo_bar},
	set_type  => 'ResultSet',
	dbID      => $branch_sets{$foo_bar}{$any_group}{dbIDs}->[0],
	set_name  => $branch_sets{$foo_bar}{$any_group}{set_names}->[0],
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