
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::RunPeaks

=head1 DESCRIPTION

This class runs a PeakCaller analysis and loads the results into the DB

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunPeaks;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception              qw( throw );
use Bio::EnsEMBL::Utils::Scalar                 qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils      qw( scalars_to_objects );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( _init_peak_caller 
                                                    _run_peak_caller );
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift;
  
  # This is used in slice_objects in BaseDB, who knows what it does
  #
  $self->param('include_slice_duplicates',1);
  
  $self->SUPER::fetch_input;

  my $set_type = $self->param_required('set_type');  
  my $pftypes  = $self->get_param_method('process_file_types', 'silent', [undef]);
  assert_ref($pftypes, 'ARRAY', 'process_file_types');
  
  if(scalar(@$pftypes) != 1){
    #Single undef value, to enable single iteration of loops
    #below where we want default value of $peak_caller->out_file($file_type);
    #See write_output for limitations of current support
    $self->throw_no_retry("RunPeaks currently only supports a single process_file_type:\n\t".
      'process_file_type => ["'.join('", "', @$pftypes).'"]');  
  }
  
  #filter_max_peaks ensures validation of max_peaks variables for known IDR analyses
  #max_peaks also needs grabbing from result_set_tracking and dataflow from IdentifySetInputs
  
  my $max_peaks = $self->get_param_method('max_peaks', 'silent');
  
  if($self->param_silent('filter_max_peaks') &&
     (! defined $max_peaks)){
    $self->throw_no_retry('The filter_max_peaks param has been set, but no max_peaks param has been set');
  }
  
  # HACK
  my $db = $self->param_required('out_db');
  my $result_set_adaptor = $db->get_ResultSetAdaptor;
  $result_set_adaptor->{file_type} = 'BAM';

  my ($fset, $result_set, $analysis);

  if($set_type eq 'ResultSet'){
    
    # This is run in run_SWEmbl_R0005_replicate 
    #die('This is never run!');
  
    $result_set = $self->fetch_Set_input('ResultSet'); 
    
    #This is likely permissive peaks for pre_IDR rep 
    
    
    #Can we not auto detect based on run_idr and is_idr_feature_type
    #No this is breaking the link between setting the analysis and the branching
    #which is based on the FeatureSet logic_name
    #This branch config is not loaded in IdentifyReplicateResultSet
    #How are we goign to do this, such that we don't risk 
    #passing the wrong analysis! i.e. It has to match the hive analysis name
    #
    
    my $peak_analysis = $self->param_required('peak_analysis');
    $analysis = scalars_to_objects($self->out_db, 'Analysis',
                                                  'fetch_by_logic_name',
                                                  [$peak_analysis])->[0];
    if(! defined $analysis){
      $self->throw_no_retry("Could not find peak_analysis in DB:\t".$peak_analysis);  
    }                            
  }
  else{ 
    $fset     = $self->fetch_Set_input('FeatureSet');
    $analysis = $fset->analysis;
    $result_set     = $self->ResultSet; 
  }
  
#   my $signal_alignment = $result_set->dbfile_path;
  my $signal_alignment  = $self->db_output_dir . '/' . $self->get_alignment_files_by_ResultSet_formats($result_set,  undef);
  my $control_alignment = $self->db_output_dir . '/' . $self->get_alignment_files_by_ResultSet_formats($result_set, 1);
#   die($control_alignment);

  # The code for building this path is duplicated in PreProcessIDR. This 
  # should be configured somewhere centrally.
  #
  $self->get_output_work_dir_methods(
    $self->peaks_output_dir 
    . '/' . $result_set->experiment->name
    . '/' . $analysis->logic_name 
  );

#   my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($result_set, undef, 1);#validate aligned flag 
#   my $control_prefix = $self->get_alignment_path_prefix_by_ResultSet($result_set, 1, 1);#and control flag 
  
  my $sam_ref_fai = $self->sam_ref_fai($result_set->epigenome->gender);  #Just in case we need to convert

  #These maybe things like extra input/reference files
  #where we don't want to store the filepath in the DB.
#   my $sensitive_caller_params = $self->param_silent($analysis->program.'_parameters') || {};
  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $species          = $self->param('species');
  my $assembly         = $self->param('assembly');
  my $epigenome_gender = $result_set->epigenome->gender;

  my $chromosome_lengths_relative = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $epigenome_gender,
    assembly         => $assembly,
    file_type        => 'chromosome_lengths_by_species_assembly',
  });
  my $reference_data_root_dir = $self->param('reference_data_root_dir');
  
  my $chromosome_length_file = $reference_data_root_dir . '/' . $chromosome_lengths_relative;
  
#   die($chromosome_length_file);

  my $sensitive_caller_params = {
    # This is used by CCAT
    -chr_file => $chromosome_length_file
  };
  
  my $pfile_path = ( defined $self->bin_dir ) ?
    $self->bin_dir.'/'.$analysis->program_file : $analysis->program_file;
    
  my $dirname = join '/', (
    $self->peaks_output_dir,
    $result_set->experiment->name,
    $analysis->logic_name
  );
  my $file_prefix    = $result_set->name.'.'.$analysis->logic_name;
  
  my $file_extension = '.txt';
  
  my $peaks_file = join '/', (
    $dirname,
    $file_prefix . '.' . $file_extension
  );
  
  # HACK CCAT can't use bam files, but the previous analysis created appropriate bed files.
  if ($analysis->logic_name eq 'ccat_histone') {
    $signal_alignment  =~ s/.bam$/.bed/;
    $control_alignment =~ s/.bam$/.bed/;
  }
  
  my @init_peak_caller_args = (
    -analysis       => $analysis,
#     -align_prefix   => $align_prefix,
#     -control_prefix => $control_prefix,
    
    -signal_alignment  => $signal_alignment,
    -control_alignment => $control_alignment,
    
    -sam_ref_fai    => $sam_ref_fai,
    -debug          => $self->debug,
    -peak_module_params =>
     {%$sensitive_caller_params,
      -program_file      => $pfile_path,
      -out_file_prefix   => $result_set->name.'.'.$analysis->logic_name,
      -out_dir           => $self->output_dir,
      -output_file       => $peaks_file,
      -convert_half_open => 1, # Before loading into DB
      #-is_half_open      => $self->param_silent('is_half_open') || 0}, # Now defined by PeakCaller or subclass defined on out/input formats   
    });
  
  my $peak_runnable = _init_peak_caller(@init_peak_caller_args);
   
  $self->set_param_method( 'peak_runnable', $peak_runnable );

  return;
}  # end sub fetch_input



#out_file may need a file_type spec here
#There is only partial support here for multiple output
#process_file_types array is specified in the analysis config
#all would be post processed uniformly, although this is only specified 
#for CCAT at present which has no $max_peaks specified
#write_output also lacks the ability to specify the extra FeatureSets required
#for the outher output files. Hence we only expect 1 process_file_type there

#mv max peak filtering to write_output?
#This seems more appropriate now, given that we need to do file/header manipulation
#for filtering. And sort if format specific
#Do we want to be able filter in run, without 'writing' output

sub run {
  my $self      = shift;
  my $max_peaks = $self->max_peaks;
   
  if( ! $self->param_silent('reload')){
    my $pcaller = $self->peak_runnable;

    if(! eval { _run_peak_caller(-peak_caller => $pcaller, 
                                 -max_peaks   => $self->max_peaks, 
                                 -file_types  => $self->process_file_types,
                                 -debug       => $self->debug); 1 }) {
      $self->throw_no_retry('Failed to call run on '.ref($pcaller)."\n$@");
    }
  }
  
  return;
}

sub write_output {
  my $self = shift;
  
  # When loading results disconnecting when inactive can lead to very poor 
  # performance. Hundreds of thousands or even millions of rows may have 
  # to be inserted into the annotated_feature table. The default setting of 
  # disconnecting after every statement would mean that a new connection has 
  # to be established for every insert statement and it would be disconnected
  # thereafter. Therefore this behaviour is deactivated here. 
  # 
  my $db = $self->get_param_method('out_db', 'required');
  $db->dbc->disconnect_when_inactive(0);
  
  if (!$self->can('FeatureSet')) {
    warn "Skipping load features as no FeatureSet is defined";
    return;
  }
  
  my $fset = $self->FeatureSet;
  
  if (!$fset) {
    die "The FeatureSet has been set, but is undefined!";
    return;
  }
  
#   $self->helper->rollback_FeatureSet($fset);

  my $af_adaptor = $fset->adaptor->db->get_AnnotatedFeatureAdaptor;    
  my $params     = {
    -file_type      => $self->process_file_types->[0],
    -processor_ref  => $self->can('store_AnnotatedFeature'),
    -processor_args => [
	$self,#Need to pass self as this calls a code ref
	$fset,
	$af_adaptor
      ]
  };

  #Need list context otherwise would get last value of list returned
  #i.e. $retvals, which we don't need here
  my ($feature_cnt) = $self->peak_runnable->process_features($params);
  
  #As CCAT is very oddly calling duplicate regions we need to skip over these in 
  #store_AnnotatedFeature
  #Hence we need to validate the feature_cnt returned matches a the DB count 
  #else throw here, rather than in store_AnnotatedFeature
  #so we can finish the load.
  #Luckily there is now flow from the CCAT analysis, so we can review these and accept them
  #throw before or after imported states?
		      
  $self->helper->debug(1, "Processed $feature_cnt features for ".$fset->name.'('.$fset->dbID.')');
  my $stored_features = $af_adaptor->generic_count('af.feature_set_id='.$fset->dbID);
  
  #Arguably this should be after the follwing test
#   $fset->adaptor->set_imported_states_by_Set($fset);
  
  if($feature_cnt != $stored_features) {
    # Could change to a normal throw if we move above the satus setting
    $self->throw_no_retry('Processed feature count does not match stored '.
      "feature count:\t$feature_cnt vs $stored_features");
  }
  return;
} ## end sub write_output

#Should we attempt to skip features which overhang end of slices?
#This may be a 'feature' of the peak caller?

# Does not handle 1/2 open here

sub store_AnnotatedFeature {
  my ( $self, $fset, $af_adaptor, $fhash ) = @_;
  my $err;
  
  my $slice;
  
  eval {  
    $slice = $self->get_Slice( $fhash->{-seq_region} );
  };

  if ( $slice ) {
    delete ${$fhash}{-seq_region};
    
    if(! eval {$af_adaptor->store( Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
                                    (%$fhash, -slice => $slice, -feature_set => $fset)); 1}){
      $err = 'Could not create and store ' .
        $fset->name . " AnnotatedFeature with attributes:\n\t" .
        join( "\n\t", ( map { "$_ => " . $fhash->{$_} } keys %$fhash ) ).
        "\n\t-slice => $slice\n\t-feature_set => $fset\n$@";
    }
  }
  else{
    $err = "Failed to get slice ".$fhash->{-seq_region};
  }

  warn $err if $err;
  return $err;
}

1;


