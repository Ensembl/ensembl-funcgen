package Bio::EnsEMBL::Funcgen::Report::FastQC;

use strict;
use base 'Bio::EnsEMBL::Funcgen::Report::Generator';

sub template {
  my $self = shift;
  my $template = $self->template_dir . '/quality_checks/fastqc.html';
  return $template;
}

sub _static_content {
  return {
    title    => 'FastQC',
  };
}

sub _dynamic_content {

  my $self = shift;
  
  my $species = $self->species;
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $read_file_adaptor = $funcgen_adaptor->get_ReadFileAdaptor;
  my $all_read_files = $read_file_adaptor->fetch_all;

# my $save_file = 'compute_fastqc_statistics';

# if (1) {

  my $datasets = $self->compute_datasets($all_read_files);
  
#   open my $out, '>', $save_file;
#   $out->print(Dumper($datasets));
#   $out->close;
# 
# } else {
# 
#   local $/;
#   open my $in, $save_file;
#   
#   my $dumped = <$in>;
#   
#   no strict;
#   $datasets = eval $dumped;
#   use strict;
#   
#   $in->close;
# }

  return {
    datasets => $datasets,
    dbc      => $funcgen_adaptor->dbc,
    species  => $species,

  };
}

sub create_read_file_dataset {

  my $self = shift;
  
  my $read_files                = shift;
  my $read_file_filter_callback = shift;

  my $filtered_peak_callings = [ grep { $read_file_filter_callback->($_) } @$read_files ];
  
  my $final_dataset = {
      fastqc_outcome => $self->fetch_fastqc_counts($filtered_peak_callings),
  };
  
  return $final_dataset;
}

sub fetch_fastqc_counts {

  my $self = shift;

  my $read_files = shift;
  
  my $pass_fail_warn_hash_factory = sub {
    return {
      'PASS' => 0,
      'FAIL' => 0,
      'WARN' => 0,
    };
  };
  
  my %fastqc_outcome_count = (

    basic_statistics             => $pass_fail_warn_hash_factory->(),
    per_base_sequence_quality    => $pass_fail_warn_hash_factory->(),
    per_tile_sequence_quality    => $pass_fail_warn_hash_factory->(),
    per_sequence_quality_scores  => $pass_fail_warn_hash_factory->(),
    per_base_sequence_content    => $pass_fail_warn_hash_factory->(),
    per_sequence_gc_content      => $pass_fail_warn_hash_factory->(),
    per_base_n_content           => $pass_fail_warn_hash_factory->(),
    sequence_length_distribution => $pass_fail_warn_hash_factory->(),
    sequence_duplication_levels  => $pass_fail_warn_hash_factory->(),
    overrepresented_sequences    => $pass_fail_warn_hash_factory->(),
    adapter_content              => $pass_fail_warn_hash_factory->(),
    kmer_content                 => $pass_fail_warn_hash_factory->(),

    run_failed                   => 0,

  );

  $self->process_read_file_qc_values_from_peak_callings(
  
    $read_files,
    
    sub {
    
      my $fastqc = shift;
      
      my $basic_statistics              = $fastqc->basic_statistics;
      my $per_base_sequence_quality     = $fastqc->per_base_sequence_quality;
      my $per_tile_sequence_quality     = $fastqc->per_tile_sequence_quality;
      my $per_sequence_quality_scores   = $fastqc->per_sequence_quality_scores;
      my $per_base_sequence_content     = $fastqc->per_base_sequence_content;
      my $per_sequence_gc_content       = $fastqc->per_sequence_gc_content;
      my $per_base_n_content            = $fastqc->per_base_n_content;
      my $sequence_length_distribution  = $fastqc->sequence_length_distribution;
      my $sequence_duplication_levels   = $fastqc->sequence_duplication_levels;
      my $overrepresented_sequences     = $fastqc->overrepresented_sequences;
      my $adapter_content               = $fastqc->adapter_content;
      my $kmer_content                  = $fastqc->kmer_content;
      my $run_failed                    = $fastqc->run_failed;

      $fastqc_outcome_count{basic_statistics}              {$basic_statistics}++;
      $fastqc_outcome_count{per_base_sequence_quality}     {$per_base_sequence_quality}++;
      $fastqc_outcome_count{per_tile_sequence_quality}     {$per_tile_sequence_quality}++;
      $fastqc_outcome_count{per_sequence_quality_scores}   {$per_sequence_quality_scores}++;
      $fastqc_outcome_count{per_base_sequence_content}     {$per_base_sequence_content}++;
      $fastqc_outcome_count{per_sequence_gc_content}       {$per_sequence_gc_content}++;
      $fastqc_outcome_count{per_base_n_content}            {$per_base_n_content}++;  
      $fastqc_outcome_count{sequence_length_distribution}  {$sequence_length_distribution}++;
      $fastqc_outcome_count{sequence_duplication_levels}   {$sequence_duplication_levels}++;
      $fastqc_outcome_count{overrepresented_sequences}     {$overrepresented_sequences}++;
      $fastqc_outcome_count{adapter_content}               {$adapter_content}++;
      $fastqc_outcome_count{kmer_content}                  {$kmer_content}++;
      
      if ($run_failed) {
        $fastqc_outcome_count{run_failed}++;
      }
    },
    
    sub {
      $fastqc_outcome_count{run_failed}++;
    }
  );
  return \%fastqc_outcome_count
}

sub process_read_file_qc_values_from_peak_callings {

  my $self = shift;

  my $read_files                    = shift;
  my $read_file_callback            = shift;
  my $read_file_run_failed_callback = shift;

  PEAK_CALLING:
  foreach my $read_file (@$read_files) {

    my $fastqc = $read_file->fetch_FastQC;
    
    if ($fastqc) {
      $read_file_callback->($fastqc);
      next PEAK_CALLING;
    }
    $read_file_run_failed_callback->($fastqc);
  }
  return;
}

sub compute_datasets {

  my $self = shift;
  
  my $all_read_files = shift;
  my @datasets;

  push 
    @datasets, 
    {
      title => 'All consortia combined',
      all => $self->create_read_file_dataset(
        $all_read_files,
        sub { return 1; }
      ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Blueprint',
      all => $self->create_read_file_dataset(
        $all_read_files,
        sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'BLUEPRINT' }
      ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'ENCODE',
      all => $self->create_read_file_dataset(
        $all_read_files,
        sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'ENCODE' }
      ),
    }
  ;

  push 
    @datasets, 
    {
      title => 'Roadmap Epigenomics',
      all => $self->create_read_file_dataset(
        $all_read_files,
        sub { shift->fetch_ReadFileExperimentalConfiguration->get_Experiment->get_ExperimentalGroup->name eq 'Roadmap Epigenomics' }
      ),
    }
  ;
  return \@datasets;
}

1;
