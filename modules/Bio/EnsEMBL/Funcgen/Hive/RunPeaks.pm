
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::RunPeaks

=head1 DESCRIPTION

This class runs a PeakCaller analysis and loads the results into the DB

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunPeaks;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;


use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

#todo Merge in IDR QC here??? this can be done on the bed files as part of run?
#actually replicates will need to run separately, so IDR QC will have to be
#in a separate analysis, which means we will likely have to load the features before
#we run IDR on the replicate bed files.

#todo rename logic_names to be generic e.g. SWEmbl_tight
#then we can have different parameter sets between species

#todo Delegate to the peak analysis runnable, such that it is entirely independant of the hive
#this will require a standard interface
#and moving a lot of the fetch_input stuff in here
#rename this RunPeaks

#input_dir and work_dir can be separate
#such that we can't point to remote fastqs
#input_dir should never be altered
#where as work_dir can be written to and cleaned up afterwards
#is this fully supported?

#params
#-output dir should never be dataflowed! ANd is really only safe when running one analysis
#as it will send all output there
#-reload

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;

  my $fset     = $self->fetch_Set_input('feature');
  my $analysis = $fset->analysis;
  my $dset     = $self->data_set; 
  my $input_sets = $dset->get_supporting_sets('input');

  if ( scalar(@$input_sets) != 1 ) {
    throw( 'Expected 1 InputSet, found ' .
           scalar(@$input_sets) . ' for DataSet ' . $dset->name );
  }

#This is required for getting files, move this to get_alignment_file_by_InputSets?
  $self->set_param_method( 'cell_type', $fset->cell_type, 'required' );

  #do we need both experiment name and logic name in here?

  $self->get_output_work_dir_methods( $self->db_output_dir . '/peaks/' .
      $fset->get_InputSet->get_Experiment->name . '/' . $analysis->logic_name );

  my $peak_module = $self->validate_package_from_path($analysis->module);
  my $formats = $peak_module->input_formats;
  my $filter_format = $self->param_silent('bam_filtered') ? undef : 'bam';    
  my $control_file;
  my $align_file = $self->get_alignment_file_by_InputSet_formats($input_sets->[0], 
                                                                 $formats,
                                                                 undef,  # control flag
                                                                 undef,  # all_formats flag
                                                                 $filter_format
                                                                );

  if ( ! $self->get_param_method( 'skip_control', 'silent' ) ) {
    #This throws if not found
    $control_file = $self->get_alignment_file_by_InputSet_format($input_sets->[0], 
                                                                 $formats, 
                                                                 1,     # control flag
                                                                 undef, # all_formats flag
                                                                 $filter_format); 
  }

  #align and control file could potentially be different formats here
  #shall we let the peak caller handle that or test here?


  #work dir is now based on the peaks and we get the input alignment file
  #above, so we don't need this dir

#my $work_dir = $self->workdir.'/'.join('/', ($self->workdir,
#                                             'alignments',
#                                             $self->species,
#                                             $self->assembly)
#                                             $experiment_name);
#workdir is currently not used!!!!
#This should be used for the samtools sort/merge stuff
#so needs defining in BaseDB?
#or just Base.pm
#Should we mirror all the subdirs from the workdir root?
#Then we can set workdir dynamically by just subing the output_dir with the work_dir

  #output_dir method should also set work_dir
  #input_dir should nevr be written to, unless it is also the output_dir

  #my $input_dir = $self->param_silent('input_dir') || $work_dir;
  #$self->set_dir_param_method('input_dir', $input_dir);

  #todo create SWEmbl runnable which does not inherit from hive
  #pass through self for access to debug/log etc?

#We must flow separately for each analysis! from the CollctionWriter preprocess job

  #This is hardcoding for packages in Bio::EnsEMBL::Funcgen::Hive
  #use full package name in analysis!

 

  #validate program_file isn't already path?

  my $pfile_path =
    ( defined $self->bin_dir ) ?
    $self->bin_dir.'/'.$analysis->program_file :
    $analysis->program_file;

  my $peak_runnable = $peak_module->new(
    -program_file      => $pfile_path,
    -parameters        => $analysis->parameters,
    -align_file        => $align_file,
    -control_file      => $control_file,
    -out_file_prefix   => $fset->name . '.' . $analysis->logic_name,
    -out_dir           => $self->output_dir,
    -convert_half_open => 1,

    #todo change this to separate flags, reload will take priority over retry?
    -is_half_open => $self->param_silent('is_half_open') ||
      0,    #default to closed coords
  );

  #How are we going to support filtered and unfiltered feature/data_set here?
  #should we keep the peak and the alignment in the same DBagnostic
  #directory?
  #No we should at least have then in analysis_logic_name subdirs
  #and these are DB specific.

  $self->set_param_method( 'peak_runnable', $peak_runnable );

  return 1;
} ## end sub fetch_input

sub run {
  my $self = shift;

  if( ! ( $self->param_silent('reload') && 
          -e $self->peak_runnable->out_file ) ){
    $self->peak_runnable->run;
  }

  return;
}

sub write_output {
  my $self = shift;
  my $fset = $self->feature_set;

  if ( $fset->has_status('IMPORTED') ) {
    throw( "Cannot imported feature into a \'IMPORTED\' FeatureSet:\t" .
           $fset->name );

    #rollback should be outside of this module!
  }

  $self->peak_runnable->store_features( $self->can('store_AnnotatedFeature'),
                              $fset,
                              $fset->adaptor->db->get_AnnotatedFeatureAdaptor );

  $fset->adaptor->set_imported_states_by_Set($fset);

  my $batch_params = $self->get_batch_params;

  #Log counts here?

#No data flow to PeaksQC here as this is done via semaphore from PreprocessAlignment

  #todo update tracking states?

  return;
} ## end sub write_output


sub store_AnnotatedFeature {
  my ( $self, $fset, $af_adaptor, $fhash ) = @_;

  if ( my $slice = $self->get_Slice( $fhash->{-seq_region} ) ) {
    delete ${$fhash}{-seq_region};

    eval {
      $af_adaptor->store( Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
                                 %$fhash, -slice => $slice -feature_set => $fset
                          ) );
    };

    if ($@) {
      throw( 'Could not create and store ' .
             $fset->name . " AnnotatedFeature with attributes:\n\t" .
             join( "\n\t", ( map { "$_ => " . $fhash->{$_} } keys %$fhash ) ) );
    }
  }

  # else this can only happen with -slices/skip_slices

  return;
}

1;

# Private function only to be called by subclasses of this class
# gets the number of reads in a sam or bed file
#sub _get_number_of_reads {
#   my ($self, $file, $file_type) = (shift, shift, shift);
#   if(($file_type ne "bed") && ($file_type ne "sam")){ throw "Only bed and sam file types supported"; }
#   my $nbr_reads = 0;
#   #If needed, add an option to check if is zipped or not...
#   my $open_cmd = "gzip -dc $file |";
#   open(FILE,$open_cmd);
#   while(<FILE>){
#     if($file_type eq "sam"){
#       next if /^\@SQ/;
#     }else {
#       next if /track name=/o;
#     }
#     $nbr_reads++;
#   }
#   close FILE;
#   return $nbr_reads;
#}

# Private function only to be called by subclasses of this class
# gets the number of reads in a sam or bed file
#sub _get_slices {
#  #NOT DONE!!
#   my ($self, $file, $file_type) = (shift, shift, shift);
#   if(($file_type ne "bed") && ($file_type ne "sam")){ throw "Only bed and sam file types supported"; }
#   my $nbr_reads = 0;
#   #If needed, add an option to check if is zipped or not...
#   my $open_cmd = "gzip -dc $file |";
#   open(FILE,$open_cmd);
#   while(<FILE>){
#     if($file_type eq "sam"){
#       next if /^@SQ/;
#     }else {
#       next if /track name=/o;
#     }
#     $nbr_reads++;
#   }
#   close FILE;
#   return $nbr_reads;
#}

=item input_dir 

=item skip_control

=item align_file

=item control_file

#is now injected as 'align_file'
#sub _input_file {
#  return $_[0]->_getter_setter('input_file',$_[1]);
#}

=cut

