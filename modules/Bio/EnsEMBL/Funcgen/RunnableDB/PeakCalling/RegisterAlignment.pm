package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterAlignment;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

sub run {

  my $self = shift;
  
  my $species = $self->param_required('species');
  my $plan    = $self->param_required('execution_plan');
  
  lock_execution_plan($plan);
  
  print Dumper($plan);
  
  my $align_plan = $plan
    ->{input}
  ;
  
  my $alignment_name   = $align_plan->{name};
  my $read_files       = $align_plan->{input}->{read_files};
  my $to_gender        = $align_plan->{to_gender};
  my $bam_file         = $align_plan->{output}->{stored};
  my $experiment_name  = $align_plan->{from_experiment};
  my $ensembl_analysis = $align_plan->{ensembl_analysis};
  
  my $remove_duplicates_ensembl_analysis = $plan->{analysis};
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
  my $is_complete     = $align_plan->{is_complete} eq TRUE;
  my $is_control      = $align_plan->{is_control}  eq TRUE;
  
  my $experiment_adaptor = Bio::EnsEMBL::Registry
  ->get_adaptor(
      $species, 
      'funcgen', 
      'Experiment'
  );
  
  # For single end reads, the name is in the "name", for paired end it is in
  # "1" and "2":
  #
  my $read_names = [ 
    map { 
      $_->{type} eq SINGLE_END 
      ? 
        (
          $_->{name} 
        )
      : 
        ( 
          $_->{1}, 
          $_->{2}
        )
    } @$read_files 
  ];
#   
#   #my $read_names = 
#   die(Dumper($read_names));
  
  my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
  
  my $alignment_with_duplicates_id
    = $self->register_alignment(
    {
      bam_file_path       => $bam_file,
      read_names          => $read_names,
      species             => $species,
      alignment_name      => $alignment_name,
      experiment          => $experiment,
      has_duplicates      => 1,
      is_control          => $is_control,
      logic_name          => $ensembl_analysis,
      to_gender           => $to_gender,
      is_complete         => $is_complete,
    }
  );
  my $remove_duplicates_plan = $plan;
  
  my $alignment_name = $remove_duplicates_plan->{name};
  my $read_names     = $read_names;
  my $bam_file       = $remove_duplicates_plan->{output}->{stored};
  
  my $alignment_no_duplicates_id
    = $self->register_alignment(
    {
      bam_file_path       => $bam_file,
      read_names          => $read_names,
      species             => $species,
      alignment_name      => $alignment_name,
      experiment          => $experiment,
      has_duplicates      => 0,
      is_control          => $is_control,
      logic_name          => $remove_duplicates_ensembl_analysis,
      to_gender           => $to_gender,
      is_complete         => $is_complete,
    }
  );
  
  my $alignment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Alignment'
    );

  my $alignment_with_duplicates = $alignment_adaptor->fetch_by_dbID($alignment_with_duplicates_id);
  my $alignment_no_duplicates   = $alignment_adaptor->fetch_by_dbID($alignment_no_duplicates_id);
  
  $alignment_with_duplicates -> deduplicated_alignment_id($alignment_no_duplicates_id);
  $alignment_no_duplicates   -> source_alignment_id($alignment_with_duplicates_id);
  
  $alignment_adaptor->update($alignment_with_duplicates);
  $alignment_adaptor->update($alignment_no_duplicates);

  return;
}

sub register_alignment {

  my $self  = shift;
  my $param = shift;
  
  my $bam_file_path  = $param->{bam_file_path};
  my $read_names     = $param->{read_names};
  my $species        = $param->{species};
  my $alignment_name = $param->{alignment_name};
  my $experiment     = $param->{experiment};
  my $has_duplicates = $param->{has_duplicates};
  my $is_control     = $param->{is_control};
  my $logic_name     = $param->{logic_name};
  my $to_gender      = $param->{to_gender};
  my $is_complete    = $param->{is_complete};

  my $read_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'ReadFile'
      );
  
  my @read_file_ids;
  foreach my $read_name (@$read_names) {
    my $read_file = $read_file_adaptor->fetch_by_name($read_name);
    
    if (! defined $read_file) {
      $self->throw("Can't fetch read file with name " . Dumper($read_name));
    }
    
    push @read_file_ids, $read_file->dbID;
  }
  
  my $data_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'DataFile'
    );
  use Bio::EnsEMBL::Funcgen::DataFile;
  my $data_file = $data_file_adaptor->fetch_by_path($bam_file_path);
  
  if (! defined $data_file) {
    $data_file = Bio::EnsEMBL::Funcgen::DataFile
      ->new(
        -table_id     => 0,
        -table_name   => 'alignment',
        -path         => $bam_file_path,
        -file_type    => 'BAM',
        -md5sum       => undef,
      );
    $data_file_adaptor->store($data_file);
  }

  my $analysis_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Analysis'
    );
  my $alignment_analysis = $analysis_adaptor->fetch_by_logic_name($logic_name);

  my $alignment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Alignment'
    );

  if (! defined $alignment_analysis) {
    die("No analysis $logic_name in the database!");
  }

  my $alignment = $alignment_adaptor->fetch_by_name($alignment_name);
  
  if (! defined $alignment) {
    $alignment = Bio::EnsEMBL::Funcgen::Alignment->new(
        -name           => $alignment_name,
        -analysis_id    => $alignment_analysis->dbID,
        -read_file_ids  => \@read_file_ids,
        -bam_file_id    => $data_file->dbID,
        -experiment_id  => $experiment->dbID,
        -has_duplicates => $has_duplicates,
        -is_control     => $is_control,
        -to_gender      => $to_gender,
        -is_complete    => $is_complete,
    );
    $alignment_adaptor->store($alignment);
  }
  return $alignment->dbID
}

1;
