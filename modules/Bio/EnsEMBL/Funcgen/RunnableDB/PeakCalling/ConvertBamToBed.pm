package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ConvertBamToBed;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Data::Dumper;
use Carp;

sub run {

  my $self = shift;
  
  my $species             = $self->param_required('species');
  my $data_root_dir       = $self->param_required('data_root_dir');
  my $convert_controls    = $self->param_required('convert_controls');

  my $execution_plan_list;
  
  if ($self->param_is_defined('execution_plan_list')) {
    $execution_plan_list = $self->param_required('execution_plan_list');
  } else {
    $execution_plan_list = [ $self->param_required('execution_plan') ];
  }

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
        summarise
  );
  
  my %seen_alignment;
  my @alignments_to_convert;
  foreach my $execution_plan (@$execution_plan_list) {
  
    my $execution_plan_expanded = resolve_nonterminal_symbols($execution_plan);
    lock_execution_plan($execution_plan_expanded);
    
    print summarise($execution_plan_expanded);
    
    my $alignment_plans = $execution_plan_expanded->{bam_to_bed};
    
    foreach my $alignment_plan (@$alignment_plans) {
    
      my $want_to_convert_this_alignment
        = (($alignment_plan->{is_control} eq TRUE ) &&  $convert_controls)
      ||  (($alignment_plan->{is_control} eq FALSE) && !$convert_controls);
    
      if ($want_to_convert_this_alignment) {
        
        my $alignment_name = $alignment_plan->{name};
        my $already_scheduled = exists $seen_alignment{$alignment_name};
        
        if (! $already_scheduled) {
          push @alignments_to_convert, $alignment_plan;
          $seen_alignment{$alignment_name} = 1;
        }
      }
    }
  }

  $self->say_with_header("Found " . @alignments_to_convert . " alignments to convert to bed format.", 1);
#   print("Found " . @alignments_to_convert . " alignments to convert to bed format.");
#   die;
  my @bam_2_bed_conversion;
  
  foreach my $alignment_to_convert (@alignments_to_convert) {
  
    my $bam_file = $data_root_dir . '/' . $alignment_to_convert->{input}->{output}->{real};
    my $bed_file = $data_root_dir . '/' . $alignment_to_convert->{output}->{real};
    
    print Dumper($alignment_to_convert);
    
    if (! -e $bam_file) {
        confess("Expected bam file $bam_file, but it doesn't exist!");
    }
    
    push @bam_2_bed_conversion, {
        from => $bam_file,
        to   => $bed_file,
    }
  }
  
  print Dumper(\@bam_2_bed_conversion);
  
  my $run_options = {
    use_bash_pipefail => 1
  };

  foreach my $bam_to_bed_conversion (@bam_2_bed_conversion) {
  
    my $bam_file = $bam_to_bed_conversion->{from};
    my $bed_file = $bam_to_bed_conversion->{to};
    
    my $cmd = "bamToBed -i $bam_file > $bed_file";
    
    use File::Basename qw( dirname );
    my $dirname = dirname($bed_file);

    use File::Path qw( make_path );
    make_path( $dirname );
    
    $self->say_with_header("Running $cmd", 1);
    
    my $has_failed = $self->run_system_command($cmd, $run_options);
    if ($has_failed) {
      $self->throw("The following command failed:\n" . $cmd)
    }
  }
  
}

1;
