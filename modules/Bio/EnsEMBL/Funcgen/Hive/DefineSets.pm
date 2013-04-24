=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::DefineSets

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineSets;

use base ('Bio::EnsEMBL::Funcgen::Hive::Base');#?



#Bio::EnsEMBL::Funcgen::Hive(::Config)
#We don't need to discriminate between Runnables and RunnableDBs anymore
#Just name the modules accordingly!


use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names 
                                               strip_param_flags run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Importer;
#use Data::Dumper;

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;


#We're assuming here that the Experiment and InputSet/SubSets have all been previously registered, 
#and now we want simply to define/fetch the data set, feature and result set based on these data.

#Do we want to add support here to add a collection result set to an existing data set which has no result_set
#and vice versa wrt peak set?

#Naming issues here as peaks set will always have analysis in name.
#Do we actually need to create a result_set only data set?

#todo -slice_import_status?

sub fetch_input {   # fetch parameters...
  my $self = shift @_;

  
  #why is this needed? Can we set this in the config as the src_root?  
  if(!defined($ENV{EFG_SRC})){ throw "NEED to define EFG_SRC"; }

  
  $self->param('name',$self->param('input_set'));
  $self->param('output_dir',$self->param('output_dir')."/".$self->param('input_set'));

  $ENV{EFG_DATA} = $self->param('output_dir');
  #Either test $ENV{EFG_DATA} or use a changed version of Importer.pm 
  #Are these ENV vars required by collection code?

  $self->set_Importer_from_params;

  return;
}

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;
  
  my $Imp = $self->param('importer');
  
  
  #This needs to call a generic method which will use
  #optioanl command line filters (ftype.ctype, project, name etc) and poll the input_set_tracking table
  #based on the expectations of this analysis and the status set for a given input_set
  #Need to handle some form of rollback here, such that we can fetch input_ids and rollback an IMPORTED analysis
  #rollback will fail if there are depedancies which already exist (IMPORTED or not)
  
  #This method can be re-used by other confs to define their inputs, if we don't dataflow between conf linking analyses
  
  #Then for each input_set_id
  #we fetch the input_set and define_DataSet based on the 


  foreach my $input_set_id(@input_set_ids){
    
    #This does not allow for us to create just the ResultSet
    #we should handle that here
    
    #do all these come from the importer?
    #can we start to separate these out a little bit
   
   
       
    
    $Imp->define_DataSet
      (
        -NAME => $set_name,
        'DESCRIPTION', 'FEATURE_CLASS', 'DISPLAY_LABEL', 'FEATURE_SET_ANALYSIS'
    'RESULT_SET_ANALYSIS', 'RESULT_SET_MODE'
    'DBADAPTOR', 'ROLLBACK', 'SLICES', 'SUPPORTING_SETS', 
      'FEATURE_TYPE', 'CELL_TYPE', 'FEATURE_CLASS'
    #dataflow here as we now use data_set_ids rather than input_set_ids 
  }
  


  return 1;
}


sub write_output {  # Create the relevant jobs
  my $self = shift @_;
  
  #Don't need explicit data flow here as it will be done via branch 1.
   
  return;
}



#Here we need to handle pre-existing data and rollback!
#Or just do this separately for now?
#Would need to use compare methods here?
#


#We will have pre-defined InputSets here
#What will input parameters actually be? built from dir names or manually specified
#or picked up from status in tracking table? 
#We could then specify filters in conjunction with tracking states to
#run subsets
#This would allow use to pre-config the tracking DB.
#This would require the status of the input_set to be updated
#else we would risk re-running a job that has already be run.
#Is this duplication of states between hive and tracking sensible?

#This should really run as one job, rather than a job for each set?
#This should take filter arguments to fetch the required input_set based on the status
#or a list of input_set names
#This in turn will create the reset of the ids to flow onto the other analyses
#what about manually creation of these input_ids, if we want to rerun just one analysis(peaks/signal)
#simply re_run this analysis, in recovery mode, with the appropriate config, such that it dataflows
#correctly
#need to explicitly flow the IDs as we don't have the correct ids at this point.
#
#This is not possible for all the other pipelines, as they won't share this analysis
#but is this the utility of the tracking status
#this will allow use to re-generate the input_ids, regardless of whether we are 
#re-running/overwriting

#Do we want to rollback separately, or allow the pipeline to do this?
#probably allow the pipeline to do it based on the rollback_level!
#rollback level can be set by an analysis, as this will be defined by the context of the analysis

#If we don't update the status of the input_set, then e can use the name to identify which sets have been IMPORTED
#is this easy?
#will this give us tracking of all aspects of the pipeline? MotifFeature association on FeatureSets?

#rollback
#Should we allow the pipeline to perform rollback of entire IMPORTED data set?
#Or just single analysis which we are rerunning i.e. features for a peak analysis
#what about the case where the collections have fallen over?
#Currently we can't rollback a result_set which has a data set.
#if this where we want to recover? i.e. allow re-use of an existing result_set(with an associated data_set)
#how does the rollback for this currently work?
#
#We need a -recover option, which will simply overwrite exisint feature data given that the set definition is unchanged
#This needs unset IMPORTED status


sub define_Sets {
  my ($self, $analysis, $input_subset, $fset_name) = @_;
  #Todo make it more generic and accept multiple input_subsets
  #Also maybe pass parameters as hash list...
  
  #aren't vars above set in params and available via *private* methods?
  #$self->_analysis(), $self->_input_file(), $self->_feature_set_name()
  
  #Global parameters set in Funcgen->fetch_input 
  my $efgdba = $self->_efgdba();
  my $set_name  =  $self->_set_name();
  my $group = $self->_group();
  my $cell_type =  $self->_cell_type();
  my $feature_type =  $self->_feature_type();
  
  
  my $iset_name = $set_name;
  my $dset_name = $fset_name;

  # set experiment: Reuse if already exists? (This comes from result sets)
  my $ea = $efgdba->get_ExperimentAdaptor;
  my $exp = $ea->fetch_by_name($set_name);

  my @date = (localtime)[5,4,3];
  $date[0] += 1900; $date[1]++;
  
  if (! defined $exp) {
    
    #Group needs to be set manually, like Cell_Type and Feature_Type
    #Do not create Group on the fly here, as it will cause concurrency issues...    
    $exp = Bio::EnsEMBL::Funcgen::Experiment->new
      (
       -NAME => $set_name,
       -EXPERIMENTAL_GROUP => $group,
       -DATE => join('-', @date),
       -PRIMARY_DESIGN_TYPE => 'binding_site_identification',
       -ADAPTOR => $ea,
      );

    ($exp) =  @{$ea->store($exp)};

  }
  throw("Can't create experiment $set_name ") unless $exp;

  my $isa = $efgdba->get_InputSetAdaptor();
  my $iset = $isa->fetch_by_name($iset_name);

  if (! defined $iset){
    
    $iset = Bio::EnsEMBL::Funcgen::InputSet->new
      (
       -name         => $iset_name,
       -experiment   => $exp,
       -feature_type => $feature_type,
       -cell_type    => $cell_type,
       -vendor       => 'SOLEXA',
       -format       => 'SEQUENCING',
       -feature_class => 'result'
       # Analysis is not being used??
       #-analysis     => $self->feature_analysis,
      );
    warn "Storing new InputSet:\t$iset_name\n";
    ($iset)  = @{$isa->store($iset)};
    
    $iset->add_new_subset($input_subset);
    $iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
  } else {

    #We only expect one subset here (? why??)... 
    #shouldn't we be adding the control file also when used?? But this is SWEmbl-specific...
    #And it should be the same file name...
    #Maybe do some file checking  here???
    warn "InputSet already exists:\t$iset_name\n";
    my @issets = @{$iset->get_InputSubsets};
    
    #if(scalar(@issets) > 1){
    #  throw("InputSet $iset_name has more than one InputSubset:\t".join("\t", (map $_->name, @issets)));
    #} elsif((scalar(@issets) == 1) && ($issets[0]->name ne $self->param('input_file'))){
    #  throw("InputSet $iset_name already has an InputSubset(".$issets[0]->name.") which does not match ".$self->param('input_file'));
    #} elsif(scalar(@issets) == 0){ #we can just add this InputSubset
    #  $iset->add_new_subset($self->input_id);
    #  $iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
    #}

    if(scalar(@issets)==0){
      #we can just add this InputSubset. Add an extra 'input:' as prefix?
      $iset->add_new_subset($input_subset);
      $iset->adaptor->store_InputSubsets($iset->get_InputSubsets);
    } else {
      #warn("Need to uncomment this section!! - it was commented just for testing purposes!!");
      #we just need to check if our file(s) is(are) already here...
      if(!$iset->get_subset_by_name($input_subset)){ 
      	#throw("InputSet $iset_name has InputSubsets(".join("\t", (map $_->name, @issets)).") which do not match ".$input_subset);
	#warn("InputSet $iset_name has InputSubsets(".join("\t", (map $_->name, @issets)).") which do not match ".$input_subset);
      }
    } 
  }

  my $fsa = $efgdba->get_FeatureSetAdaptor();
  my $fset = $fsa->fetch_by_name($fset_name);
  
  if ( ! defined $fset ) {
    
    $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
      (
       -analysis      => $analysis,
       -feature_type  => $feature_type,
       -cell_type     => $cell_type,
       -name          => $fset_name,
       -feature_class => 'annotated',
       -experiment_id => $exp->dbID,
       #The adaptor is needed to store!
       -adaptor       => $fsa

      );
    
    warn "Storing new FeatureSet:\t$fset_name\n";
    ($fset) = @{$fsa->store($fset)};
		
  } 
  else {
    warn "FeatureSet already exists:\t$fset_name\n";

    if(@{$efgdba->get_AnnotatedFeatureAdaptor->fetch_all_by_FeatureSets([$fset])}){
      throw "Feature Set $set_name already contains data. Please rollback before rerunning";
    }

  }

  my $dsa = $efgdba->get_DataSetAdaptor;
  my $dset = $dsa->fetch_by_name($dset_name);
  
  
    if ( ! defined $dset ) {
      
      $dset = Bio::EnsEMBL::Funcgen::DataSet->new
        (
         -SUPPORTING_SETS     => [$iset],
         -FEATURE_SET         => $fset,
         -DISPLAYABLE         => 1,
         -NAME                => $dset_name,
         -SUPPORTING_SET_TYPE => 'input',
        );
      
      warn "Storing new DataSet:\t$dset_name\n";
      ($dset) = @{$dsa->store($dset)}
    } 
  else {
    
    warn "DataSet already exists:\t$dset_name\n";
    
    # need to check whether InputSets and supporting_sets are the same and 
    # possibly add InputSet to supporting_sets
    
    my $ssets = $dset->get_supporting_sets();
    
    my %ssets_dbIDs = ();
    map { $ssets_dbIDs{$_->dbID}='' } (@{$ssets});
    $dset->add_supporting_sets([ $iset ]) if (! exists $ssets_dbIDs{$iset->dbID}); 
	  
  }
  
}








1;
