package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CallPeaks;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use JSON::XS qw(encode_json);

use constant {
  BRANCH_STORE_PEAKS => 2,
};

sub run {

  my $self = shift;
  
  my $plan          = $self->param_required('execution_plan');
  my $species       = $self->param_required('species');
  my $data_root_dir = $self->param_required('data_root_dir');
  my $tempdir       = $self->param_required('tempdir');

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
        summarise
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  
  print summarise($plan_expanded);
  
  my $call_peaks_plan = $plan_expanded->{call_peaks};
  
  my $signal_bam_file 
    = $call_peaks_plan
      ->{input}
      ->{signal}
      ->{output}
      ->{real}
  ;
  
  my $control_bam_file;
  
  if ($plan_expanded->{meta_data}->{experiment_has_control} eq FALSE) {
    $control_bam_file = undef;
  } else {
    $control_bam_file 
      = $call_peaks_plan
        ->{input}
        ->{control}
        ->{output}
        ->{real}
    ;
  }
  
  my $chromosome_lengths_by_species_assembly
    = $call_peaks_plan
      ->{chromosome_lengths_by_species_assembly};

  my $name                  = $call_peaks_plan->{name};
  my $output_dir            = $call_peaks_plan->{output_dir};
  my $samtools_fasta_index  = $call_peaks_plan->{samtools_fasta_index};
  my $peak_calling_strategy = $call_peaks_plan->{peak_calling_strategy};
  my $logic_name            = $call_peaks_plan->{analysis};
  
  if (! defined $logic_name) {
    die("No analysis logic name was defined for peak calling!");
  }
  
  my $analysis_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Analysis'
    );
  my $peak_analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);
  
  if (! defined $peak_analysis) {
    die("Can't fetch peak analysis with logic name $logic_name!");
  }
  
  my $reference_data_root_dir = $self->param('reference_data_root_dir');

  use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( 
    _init_peak_caller 
    _run_peak_caller
  );
  
  my $out_dir            = $tempdir . '/' . $output_dir;
  
  use File::Path qw( make_path );
  make_path( $out_dir );

  my $peaks_to_load_file = $out_dir . "/${name}.peaks_to_load.json";
  my $peaks_file         = $out_dir . "/${name}.output.txt";
  my $out_file_prefix    = $name;

  my $peak_caller_binary = $peak_analysis->program;
  
  open my $peak_out, '>' . $peaks_to_load_file || die("Couldn't open file $peaks_to_load_file for writing!");
  
  my $sensitive_caller_params = {
    -chr_file => $reference_data_root_dir . '/' . $chromosome_lengths_by_species_assembly,
    -FDR_THRESHOLD => 1,
  };

  my $control_alignment;
  
  if ($control_bam_file) {
    $control_alignment = $data_root_dir . '/' . $control_bam_file;
  } else {
    $control_alignment = undef;
  }

  # Set to one for ccat
  my $is_half_open = 0;
  
  my @init_peak_caller_args = (
    -analysis           => $peak_analysis,
    -signal_alignment   => $data_root_dir . '/' . $signal_bam_file,
    -control_alignment  => $control_alignment,
    -sam_ref_fai        => $reference_data_root_dir . '/' . $samtools_fasta_index,
    -debug              => $self->debug,
    -is_half_open       => $is_half_open,
    -peak_module_params => {
        %$sensitive_caller_params,
        -program_file      => $peak_caller_binary,
        -out_file_prefix   => $out_file_prefix,
        -out_dir           => $out_dir,
        -output_file       => $peaks_file,
        # Before loading into DB
        -convert_half_open => 1,
    });

  my $peak_runnable = _init_peak_caller(@init_peak_caller_args);
  
  my $store_Peak = sub {
    my ( $self, undef, undef, $feature_hash ) = @_;
    $peak_out->print(encode_json($feature_hash)."\n");
    return;
  };

  my $file_type;
  my @max_peaks = ( 'one_item' );
  
  if ($peak_calling_strategy eq CALL_BROAD_PEAKS) {
    $file_type = 'significant_region';
  }
  if ($peak_calling_strategy eq CALL_NARROW_PEAKS) {
    $file_type = 'txt';
  }
  if ($peak_calling_strategy eq CALL_TIGHT_PEAKS) {
    $file_type = 'txt';
  }
  
  my $idr_strategy = $plan->{idr}->{strategy};
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw( SKIP_IDR );
  
  if ($idr_strategy eq SKIP_IDR) {
    @max_peaks = ();
  }
  
  if ($idr_strategy ne SKIP_IDR) {
  
    my $experiment_name = $plan->{idr}->{name};
    my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'experiment');
    my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
    if (! defined $experiment) {
      confess("Can't fetch experiment with name $experiment_name!");
    }
    my $idr_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'idr');
    my $idr = $idr_adaptor->_fetch_by_experiment_id($experiment->dbID);
    
    if (! defined $idr) {
      confess("Couldn't find idr for experiment $experiment_name!");
    }
    @max_peaks = ( -max_peaks => $idr->max_peaks );
  }

  if (! defined $file_type) {
    $self->throw(
      "file_type wasn't set for peak_calling_strategy=$peak_calling_strategy"
    );
  }
  if (@max_peaks == 1) {
    $self->throw(
      "max_peaks wasn't set for peak_calling_strategy=$peak_calling_strategy"
    );
  }

  my $peak_calling_succeeded = undef;
  my $error_message;

  my $TIMEOUT_IN_SECONDS = 10 * 3600;
  
  my $childPid;
  
  eval {
      local $SIG{ALRM} = sub { die "alarm" };
      
      local $SIG{INT} = sub { 
        $self->warning("Got INT signal!");
        sleep(30);
      };
      
      alarm($TIMEOUT_IN_SECONDS);
      
      if ($childPid = fork()) {
        wait();
      } else {
        _run_peak_caller(
            -peak_caller => $peak_runnable, 
            -debug       => $self->debug,
            @max_peaks
        );
      }

      $peak_calling_succeeded = 1;
      alarm(0);
      
      my $params = {
        -file_type      => $file_type,
        -processor_ref  => $store_Peak,
        -processor_args => [
          undef,
          undef,
          undef
        ]
      };
      $peak_runnable->process_features($params);
      $peak_out->close;
      
      if (! -e $peaks_to_load_file) {
        $self->warning("The file $peaks_to_load_file was not created!");
        sleep(20);
        $self->throw("Can't continue without the peaks.");
      }
  };
  if ($@) {
      $peak_calling_succeeded = 0;
      
      my $killed_for_timeout = $@ =~ /alarm/;
      
      if ($killed_for_timeout) {
        kill 9, $childPid;
        wait;
        $error_message = "Peak calling job " . $self->input_job->dbID . " timed out after $TIMEOUT_IN_SECONDS seconds.";
      }
      
      if (! $killed_for_timeout) {
        #$self->throw($@);
        $error_message = "Peak calling job " . $self->input_job->dbID . " failed with error message: $@";
      }
      print "\n\n------------------------------------------->$error_message";
      $self->warning($error_message);
  }

  $self->dataflow_output_id(
    {
      'species'                => $species,
      'peaks_to_load_file'     => $peaks_to_load_file,
      'execution_plan'         => $plan,
      'peak_calling_succeeded' => $peak_calling_succeeded,
      'error_message'          => $error_message
    }, 
    BRANCH_STORE_PEAKS
  );
}

1;
