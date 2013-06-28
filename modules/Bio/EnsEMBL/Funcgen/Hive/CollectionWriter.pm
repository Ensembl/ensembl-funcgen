=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::CollectionWriter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::CollectionWriter;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseImporter');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names
                                               get_feature_file ); 
                                               #strip_param_args strip_param_flags run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw(throw);

#use Data::Dumper;

#global values for the Helper... maybe pass as parameters...
$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;

sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  $self->SUPER::fetch_input;  
  $self->helper->debug(1, "CollectionWriter::fetch_input after SUPER::fetch_input");

  my $set_type       = $self->param_required('set_type');
  my $adaptor_method = 'get_'.$set_type.'Adaptor'; 
  my $dbid     = $self->param_required('dbID');
  my $db       = $self->param_required('out_db');
  my $set_name = $self->param_required('set_name');
  
  my $set     = $self->out_db->$adaptor_method->fetch_by_dbID($dbid);
  
  if(! defined $set){
    throw("Could not fetch $set_type with dbID $dbid ($set_name)"); 
  }
  
  #$self->param('data_set', $dset);
  my $rset;
  
  if($set_type eq 'DataSet'){
    my @rsets = @{$set->get_supporting_sets('result')};
                    
    if(scalar @rsets != 1){
      throw("Expected 1 ResultSet, found:\n".$self->helper->dump(\@rsets));
    }
    
    $rset = $rsets[0];
  }
  elsif($set_type eq 'ResultSet'){
    $rset = $set;
  }
  else{
   throw($set_type.' set_type not supported. Must be DataSet or ResultSet'); 
  }
  
  
  $self->param('result_set', $rset);


  #Need to validate slices here, mandatory for fanned jobs and wrap up


  
  $self->helper->debug(1, "CollectionWriter::fetch_input got ResultSet");
  

  #Define current default Importer params for this analysis
     
   
  
  #where will the bed file be generated?
  #file_format should dictate output and input file types
  #i.e. lowest common denominator until we implement a bam parser
  #for peak import and collection generation
  my $input_sets = $rset->get_support;  
  my $align_dir  = $self->alignment_dir($input_sets);
  #This will have validated that input_set are all part of the same experiment
  #(and they have the same alignment logic_name, when this is implemented)
  
  my $suffix = lc($self->param_required('feature_file_format'));
    
  
  #Is this true, doesn't bam have samse too?
  if ($suffix ne 'bam'){
    $suffix = 'samse.'.$suffix;
  }
    
  #Currently don't support > 1 result_file, but that is caught in the Importer somewhere
  #Move this logic to the Importer and make input_files/data_dir optional
  my @result_files;
  $self->param('cell_type', $rset->cell_type); 
  my $bam_params = {header => $self->sam_header}; #Just in case we need to convert 

  for my $iset(@$input_sets){
    push @result_files, get_feature_file($align_dir.'/'.$iset->name.".${suffix}", $bam_params);    
  }
  
  $self->helper->debug(1, "CollectionWriter::fetch_input setting new_importer_params");
  
  
    #Arguably, some of the should be set in default_importer_params
    #but this does the same jobs, and prevents having to flow two separate hashes
  
  
  
  $self->param('output_dir', 
               $self->param('output_dir').'/result_feature/'.$rset->name);
  
  
  
  $self->new_Importer_params( 
   {
    -input_feature_class => 'result',
    -format              => 'SEQUENCING',#This needs changing to different types of seq
    -recover             => 1, 
    -force               => 1, #for store_window_bins_by_Slice_Parser
    #-feature_analysis  => #Don't need this now as we already have a reseult set
    #and is not mandatory
       
    -output_set  => $rset,
    -input_files => \@result_files,
	 
    #Can't we set this based on the output_set and feature_class?
    -output_dir  => $self->param('output_dir'),
    #we really want to send the intermediate files to work_dir
    #but we don't have support for this in the Importer
    #-name        => $input_sets->[0]->name,
   });
    
  
  #todo why do we have to set EFG_DATA?

  $ENV{EFG_DATA} = $self->param('output_dir');
    
  
  my $imp = $self->get_Importer;
  #This uses the params from the input_id dataflowed from the prepare step
 
 
  if((! $self->param_silent('wrap_up')) &&
    ($imp->prepared)){
    $self->param_required('slices');     
  }
 

  return 1;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;
  my $Imp = $self->param('importer');
  
  my @slices = @{&generate_slices_from_names($Imp->slice_adaptor, 
                                             $self->param_silent('slices'), 
                                             $self->param_silent('skip_slices'), 
                                             'toplevel', 
                                             0, 1)};#nonref, incdups
  $Imp->slices(\@slices);
  
  
  
  if(! $self->param_silent('wrap_up')){ #Prepare or write slice col
  
    
  
    if( ! $Imp->prepared ){
      #$Imp->init_experiment_import; #Define and store sets once here before setting off parallel farm jobs.
      #We already have the experiment and inputset/subsets defined
      #a lot of other stuff in here was array specific (dirs, noram analysis etc)
      #we do however need to makesure result_files are defined
      


      #Preparing data...
      $Imp->read_and_import_data('prepare');

      #This should be stricted to those slices which are seen to avoid errors
      #when trying to open col files!

      my %seen_slices = %{$Imp->slice_cache};

      #Dataflow only jobs for slices that the imported has seen
      foreach my $slice (values %seen_slices){
        #my $new_input_id = eval($self->input_id); 
        #wtf? This would return nothing?! is this missing some params from the original input_id?
        #and would likely be resetting val in the original input_id
        #These were cell_type, feature_type, input_set and result_file
        #and are not in pipeline_wide_parameters    
       #$new_input_id->{slice} = $slice->seq_region_name;
       #$new_input_id->{result_file} = $Imp->output_file;
       #$new_input_id->{total_features} = $Imp->counts('total_features');
       #push @rep_out_ids, $new_input_id;
    
    
       #Dataflow explicitly here for each each slice
    
       $self->dataflow_output_id(
        {dbID           => $self->param('dbID'),
         set_name       => $self->param('set_name'),#mainly for readability
         set_type       => $self->param('set_type'),
         slices         => [$slice->seq_region_name],
         new_importer_params => 
          {
           -input_files    => [$Imp->output_file],
           -total_features => $Imp->counts('total_features'),
           -batch_job      => 1, 
	       -prepared       => 1,  
	        #These can be replaced with rset dbID
            #cell_type      => $self->param('cell_type'),
            #feature_type   => 
            #input_set      =>
          }
        },  
        2); #branch 2 for fan of slice jobs
      }
    
      
      #We should have thrown already if the preparation has failed
      #i.e. we don't have the correct rollback defined
      
      #Do we want jobs to fail if the don't have the right states, or silently overlook them
      #e.g. we have a bunch of inputs defined, but only some need their collections generating
      #We don't really want failed jobs blocking the pipeline
      #we just want to run the jobs it needs to and ignore the rest
      #unless there is a real failure, where we would throw in the Importer
      #so we want to 'recover' i.e. simply use the available ResultSets if they have the appropriate 
      #states.
      
      #So let's say we have some IMPORTED sets, and some unimported sets, how do we differentiate
      #with the rollback mechanism?
      #Current rollback modes are full(delete everything) or recover(allow th
      #
      
            
      #Can only dataflow here if we have successfully managed to prepare
      
      $self->dataflow_output_id({dbID     => $self->param('dbID'),
                                 set_name => $self->param('set_name'),#mainly for readability
                                 set_type => $self->param('set_type'),
                                 slices   => [keys %seen_slices]}, 1);
    } 
    else { #These are the fanned slice jobs  
      #Test for new_import_params here?
      #Or might these already be specified?
    
       
      #The params from the input_id are picked up in set_Importer_from_params called in fetch_input
      #This generates a slice .col file
      #$Imp->register_experiment;
      $Imp->read_and_import_data(undef);
      
      #Need to pass slices names to wrap up step from here or from
      #the prepare step?
      
      
    } 
  }
  else{ #Wrap up
    #set name and dbID will always be result_set here


    # run the merge script here...
    
    #todo put this bulk of this in a package such that we can call directly
    #and remove cmdline call from here
    #can optionally maintain script
    my $rset = $self->param_required('result_set');
    
    #don't have acces to these vars anymore as they have been overwritten by the 
    #out_db DBAdaptor
    #move to module? Collector
    #How is this accessed via the importer?
    #update the script to use that module
    
    
    warn "HARDCODING force_overwrite for merge_and_index_slice_collections!";
    $Imp->feature_adaptor->merge_and_index_slice_collections($Imp->output_set, 
                                                             $Imp->slices, 
                                                             $self->param('output_dir'),
                                                             1);
    
    
    #my $cmd = $ENV{EFG_SRC}.'/scripts/import/merge_and_index_collections.pl '.
    # ' -dbhost '.$self->param('host').
	# ' -dbport '.$self->param('port').
	# ' -dbuser '.$self->param('user').
	# ' -dbpass '. $self->param('pass').
	# ' -dbname '. $self->param('dbname').
	# ' -data_dir '. $self->param('output_dir').
	# ' -result_set_name '.$rset->name;
	#todo, change this to dbID
		  
    #run_system_cmd($cmd);
 
    #set appropriate states...
    $rset->adaptor->set_imported_states_by_Set($rset);
  }

  return;
}

# No more dataflow required here
sub write_output { return; }

1;
