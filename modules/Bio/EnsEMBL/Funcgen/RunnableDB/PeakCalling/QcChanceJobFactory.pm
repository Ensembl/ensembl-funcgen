package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcChanceJobFactory;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
    my $self = shift;
    
    my $species         = $self->param_required('species');
    #my $experiment_name = $self->param_required('experiment');
    my $alignment_name  = $self->param_required('alignment');
    my $data_root_dir   = $self->param_required('data_root_dir');
    my $tempdir         = $self->param_required('tempdir_peak_calling');
    
    my $alignment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Alignment');

    if ($alignment_name eq NO_CONTROL_FLAG) {
      $self->say_with_header("Skipping this alignment, because it is a hack to indicate a fake experiment. (Should be done properly one day.)", 1);
      return;
    }
    
    my $signal_alignment  = $alignment_adaptor->fetch_by_name($alignment_name);
    
   if (! defined $signal_alignment) {
      $self->throw("Can't find signal alignment with name: $alignment_name");
   }
    
    my $signal_experiment = $signal_alignment->fetch_Experiment;
    
    if ($signal_experiment->is_control) {
        $self->warning("Chance is not run on controls, no jobs will be generated.");
        return;
    }
    
    my $control_experiment = $signal_experiment->get_control;
    
    if (! defined $control_experiment) {
      $self->say_with_header("This experiment has no control, so chance will be skipped.", 1);
      return;
    }
    
    my $control_alignment = $alignment_adaptor->fetch_complete_deduplicated_by_Experiment($control_experiment);
    
    if (! defined $control_alignment) {
        die("Can't find control alignment for " . $signal_alignment->name);
    }

    my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');

    my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
    my $default_assembly = $default_chromosome_coordsystem->version;
    my $epigenome = $signal_experiment->epigenome;

    my $epigenome_gender          = $epigenome->gender;
    my $epigenome_production_name = $epigenome->production_name;

    # epigenome_gender is used as a substitution parameter in another ehive 
    # parameter. If it is undefined, it will lead to an error.
    #
    if (! defined $epigenome_gender) {
        $epigenome_gender = '';
    }

    my $chance_bin_file = 'chance_bin_file_for_'.$epigenome_production_name.'_'.$epigenome_gender.'.bed';

    my $signal_bam_data_file = $signal_alignment->fetch_bam_DataFile;
    my $signal_bam_file_name = $signal_bam_data_file->path;

    my $control_bam_data_file = $control_alignment->fetch_bam_DataFile;
    my $control_bam_file_name = $control_bam_data_file->path;

    my $signal_bam_file  = $data_root_dir . '/' . $species . '/' . $default_assembly . '/' . $signal_bam_file_name;
    my $control_bam_file = $data_root_dir . '/' . $species . '/' . $default_assembly . '/' . $control_bam_file_name;

    if (! -e $signal_bam_file) {
        die("$signal_bam_file doesn't exist!");
    }
    if (! -e $control_bam_file) {
        die("$control_bam_file doesn't exist!");
    }

    use File::Basename;
    (my $signal_bam_file_base_name,  my $signal_bam_directory)  = fileparse($signal_bam_file);
    (my $control_bam_file_base_name, my $control_bam_directory) = fileparse($control_bam_file);
    
    my $experiment_name = $signal_experiment->name;
    my $experiment_id   = $signal_experiment->dbID;
    
    my $overridden_tempdir = "$tempdir/$species/qc_chance/$experiment_id";

    my $input_id = {

        # This is for JobFactoryArgenrich to generate the jobs for indexing the
        # pairs of signal and control bam files.
        #
        # signal_alignment_id is not used, it is included, so for jobs with shared controls ehive doesn't skip them.
        # The output is expected in the accumulator
        #
        column_names => [ 'kind', 'file', 'sourcedir', 'chance_tempdir', 'signal_alignment_id' ],
        inputlist    => [
            [ 'signal',  $signal_bam_file_base_name,   $signal_bam_directory,  $overridden_tempdir, $signal_alignment->dbID ],
            [ 'control', $control_bam_file_base_name,  $control_bam_directory, $overridden_tempdir, $signal_alignment->dbID ],
        ],
        alignment_has_duplicates  => $signal_alignment->has_duplicates,
        experiment_name           => $experiment_name,
        epigenome_gender          => $epigenome_gender,
        epigenome_production_name => $epigenome_production_name,
        chance_bin_file           => $chance_bin_file,
        assembly                  => $default_assembly,
        species                   => $species,
        argenrich_outfile         => 'argenrich_outfile.txt',
        
        signal_alignment          => $alignment_name,
        control_alignment         => $control_alignment->name,
        
        chance_tempdir => $overridden_tempdir,

        # Connection details for the db to which the results will be written
        tracking_db_user   => $alignment_adaptor->dbc->user,
        tracking_db_pass   => $alignment_adaptor->dbc->password,
        tracking_db_host   => $alignment_adaptor->dbc->host,
        tracking_db_name   => $alignment_adaptor->dbc->dbname,
        tracking_db_port   => $alignment_adaptor->dbc->port,
    };

  $self->dataflow_output_id($input_id, 2);
  return;
}

1;
