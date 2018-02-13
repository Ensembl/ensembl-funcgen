package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SplitFastq;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::Utils::Fastq::Processor;
use Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser;

use constant {
  BRANCH_ALIGN => 2,
  BRANCH_MERGE => 3,
};

sub run {

  my $self = shift;
  my $species      = $self->param_required('species');
  my $plan         = $self->param_required('execution_plan');
  my $tempdir      = $self->param_required('tempdir');
  
  my $num_records_per_split_fastq;
  if ($self->param_is_defined('num_records_per_split_fastq')) {
    $num_records_per_split_fastq = $self->param('num_records_per_split_fastq');
  }
  
  if (! defined $num_records_per_split_fastq) {
    $num_records_per_split_fastq = 1_000_000;
  }
  
  print Dumper($plan);
  
  my $align_plan = $plan
    ->{input}
  ;
  
  my $alignment_name   = $align_plan->{name};
  my $read_files       = $align_plan->{input}->{read_files};
  my $to_gender        = $align_plan->{to_gender};
  my $to_assembly      = $align_plan->{to_assembly};
  my $bam_file         = $align_plan->{output}->{real};
  my $ensembl_analysis = $align_plan->{ensembl_analysis};

  $self->say_with_header("Creating alignment $alignment_name");
  
  my $read_file_adaptor
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFile'
  );
  
  my @chunks_to_be_merged;
  
  foreach my $read_file (@$read_files) {

    my $read_file_names;

    if ($read_file->{type} eq 'paired_end') {
      $read_file_names = [
        $read_file->{1},
        $read_file->{2},
      ];

      my $read_file_objects = [ map { $read_file_adaptor->fetch_by_name($_) } @$read_file_names ];
      
      if ($read_file_objects->[0]->file_size != $read_file_objects->[1]->file_size) {
        die("The read files have different lengths!");
      }
    }
    
    if ($read_file->{type} eq 'single_end') {
      $read_file_names = [
        $read_file->{name},
      ];
    }
    
    my $directory_name = join '_', @$read_file_names;
    
    my $alignment_read_file_name = $alignment_name . '/' . $directory_name;
    
    my $dataflow_fastq_chunk = sub {

        my $param = shift;
        
        my $absolute_chunk_file_name = $param->{absolute_chunk_file_name};
        my $current_file_number      = $param->{current_file_number};
        my $current_temp_dir         = $param->{current_temp_dir};
        
        $self->say_with_header("Created fastq chunk $absolute_chunk_file_name");
        
        my $chunk_bam_file = $current_temp_dir . '/' . $alignment_name . '_' . $directory_name . '_' . $current_file_number . '.bam';
        
        push @chunks_to_be_merged, $chunk_bam_file;
        
        my $dataflow_output_id
          = {
            'fastq_file'       => $absolute_chunk_file_name,
            'tempdir'          => $current_temp_dir,
            'bam_file'         => $chunk_bam_file,
            'species'          => $species,
            'to_gender'        => $to_gender,
            'to_assembly'      => $to_assembly,
            'ensembl_analysis' => $ensembl_analysis,
          };
          $self->dataflow_output_id(
            $dataflow_output_id, 
            BRANCH_ALIGN
          );
#         print Dumper(
#           [
#             $dataflow_output_id,
#             $chunk_bam_file
#           ]
#         );
    };

    my $fastq_record_processor = Bio::EnsEMBL::Funcgen::Utils::Fastq::Processor->new(
      -tempdir                     => $tempdir,
      -species                     => $species,
      -alignment_name              => $alignment_read_file_name,
      -num_records_per_split_fastq => $num_records_per_split_fastq,
      -chunk_created_callback      => $dataflow_fastq_chunk,
    );

    split_read_file({
      read_file_names        => $read_file_names,
      dataflow_fastq_chunk   => $dataflow_fastq_chunk,
      fastq_record_processor => $fastq_record_processor,
      read_file_adaptor      => $read_file_adaptor,
    });
  }
  
  $self->dataflow_output_id(
    {
      'species' => $species,
      'chunks'  => \@chunks_to_be_merged,
      'plan'    => $plan,
    }, 
    BRANCH_MERGE
  );
  return;
}

sub split_read_file {

  my $param = shift;
  
  my $read_file_names        = $param->{read_file_names};
  my $dataflow_fastq_chunk   = $param->{dataflow_fastq_chunk};
  my $fastq_record_processor = $param->{fastq_record_processor};
  my $read_file_adaptor = $param->{read_file_adaptor};
  
  my $parser = Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser->new;
  
  my @read_file_objects = map { $read_file_adaptor->fetch_by_name($_) } @$read_file_names;
  my @fastq_files       = map { $_->file } @read_file_objects;
  
  use List::Util qw( uniq );
  my $same_number_of_reads_in_every_file 
    = 1 == uniq map { $_->file_size } @read_file_objects;
  
  if (! $same_number_of_reads_in_every_file) {
    confess(
      "These fastq files are paired, but don't have the same number of reads!\n"
      . Dumper(\@fastq_files)
      . "\n"
      . ( join ' vs ', map { $_->file_size } @read_file_objects )
    );
  }
  
#   my $max   = 500;
#   my $count =   0;

  my @cmds          = map { "zcat $_ |"                                     } @fastq_files;
  my @file_handles  = map { open my $fh, $_ or die("Can't execute $_"); $fh } @cmds;

  use List::Util qw( none any );

  my $all_records_read = undef;

  $fastq_record_processor->tags([ 1 .. @fastq_files ]);
  $fastq_record_processor->init;

  while (
      (! $all_records_read)
#       && 
#       ($count<$max)
    ) {
    
    my @fastq_records = map  { $parser->parse_next_record($_) } @file_handles;
    $all_records_read = none { defined $_ } @fastq_records;
    #$all_records_read = any { ! defined $_ } @fastq_records;
    
    $fastq_record_processor->process(@fastq_records);
#     $count++;
  }
  map { $_->close || die "Can't close file!" } @file_handles;
  #map { $_->close } @file_handles;
  $fastq_record_processor->flush;
  return;
}


1;
