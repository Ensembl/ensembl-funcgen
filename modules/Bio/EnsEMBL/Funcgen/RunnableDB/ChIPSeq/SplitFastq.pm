package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SplitFastq;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_ALIGN => 2,
  BRANCH_MERGE => 3,
};

sub run {

  my $self = shift;
  my $species      = $self->param_required('species');
  my $plan         = $self->param_required('plan');
  my $tempdir      = $self->param_required('tempdir');
  my $in_test_mode = $self->param('in_test_mode');
  
  my $num_records_per_split_fastq = $self->param('num_records_per_split_fastq');
  
  my $max_records = undef;
  
  if (! defined $num_records_per_split_fastq) {
    $num_records_per_split_fastq = 5_000_000;
  }
  
  if ($in_test_mode) {
    $max_records = 2_000;
    $num_records_per_split_fastq = 500;
  }
  
  my $align_plan = $plan
    ->{remove_duplicates}
    ->{align}
  ;
  
  my $alignment_name = $align_plan->{name};
  my $read_files     = $align_plan->{read_files};
  my $to_gender      = $align_plan->{to_gender};
  my $to_assembly    = $align_plan->{to_assembly};
  my $bam_file       = $align_plan->{bam_file};
  
  $self->say_with_header("Creating alignment $alignment_name");
  
  my $read_file_adaptor
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFile'
  );
  
  my @read_file_names;
  foreach my $read_file_name (@$read_files) {
    my $read_file = $read_file_adaptor->fetch_by_name($read_file_name);
    push @read_file_names, $read_file->file;
  }
  
  my $file_list = join ' ', @read_file_names;
  my $cmd = "zcat $file_list |";
  
  $self->say_with_header("Reading fastq files by running:");
  $self->say_with_header($cmd);
  
  my @chunks;

  my $dataflow_fastq_chunk = sub {

    my $param = shift;
    
    my $absolute_chunk_file_name = $param->{absolute_chunk_file_name};
    my $current_file_number      = $param->{current_file_number};
    my $current_temp_dir         = $param->{current_temp_dir};
    
    $self->say_with_header("Created fastq chunk $absolute_chunk_file_name");
    
    my $chunk_bam_file = $current_temp_dir . '/' . $alignment_name . '_' . $current_file_number . '.bam';
    
    push @chunks, $chunk_bam_file;
    
    $self->dataflow_output_id( 
      {
        'fastq_file'  => $absolute_chunk_file_name,
        'species'     => $species,
        'tempdir'     => $current_temp_dir,
        'to_gender'   => $to_gender,
        'to_assembly' => $to_assembly,
        'bam_file'    => $chunk_bam_file,
      }, 
      BRANCH_ALIGN
    );
  };

  my $fastq_record_processor = FastqRecordProcessor->new(
    -tempdir                     => $tempdir,
    -species                     => $species,
    -alignment_name              => $alignment_name,
    -num_records_per_split_fastq => $num_records_per_split_fastq,
    -chunk_created_callback      => $dataflow_fastq_chunk 
  );
  
  my $process_fastq_record = sub {
    my $record = shift;
    $fastq_record_processor->process($record);
  };

  eval {
    $self->parse_fastq_from_cmd({
      max_records    => $max_records,
      cmd            => $cmd,
      process_record => $process_fastq_record,
    });
  };
  if ($@) {
    die(
      "Error processing\n" 
      . Dumper(\@read_file_names) . "\n"
      . "Got: $@"
    );
  }
  $self->dataflow_output_id( 
    {
      'species' => $species,
      'chunks'  => \@chunks,
      'plan'    => $plan,
    }, 
    BRANCH_MERGE
  );
  return;
}

sub parse_fastq_from_cmd {

  my $self  = shift;
  my $param = shift;

  my $cmd = $param->{cmd};
  
  open my $fh, $cmd;
  $param->{fh} = $fh;
  
  $self->parse_fastq($param);
  
  $fh->close;
  return;
}

sub parse_fastq {

  my $self  = shift;
  my $param = shift;

  my $fh             = $param->{fh};
  my $process_record = $param->{process_record};
  my $max_records    = $param->{max_records};

  my $num_lines = 0;
  
  FASTQ_RECORD:
  while (defined (my $current_line = <$fh>)) {
  
    if ($current_line !~ /^\@/) {
    
      # Some of Encode's fastq files in .tgz format have a bit of garbage at 
      # the start of the file.
      #
      warn "Invalid fastq record: $current_line";
      next FASTQ_RECORD;
    }
    
    my @current_fastq_record;
    chomp $current_line;
    push @current_fastq_record, $current_line;
    
    foreach my $line_of_record (2, 3, 4) {
      my $line = <$fh>;
      chomp $line;
      push @current_fastq_record, $line;
    }
    
    $process_record->(\@current_fastq_record);
    
    $num_lines++;
    if (
      defined $max_records 
      && $num_lines == $max_records
    ) {
      last FASTQ_RECORD;
    }
  }
  return;
}

1;

package FastqRecordProcessor;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    tempdir                     => 'tempdir',
    species                     => 'species',
    alignment_name              => 'alignment_name',
    num_records_per_split_fastq => 'num_records_per_split_fastq',
    chunk_created_callback      => 'chunk_created_callback',
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub tempdir                     { return shift->_generic_get_or_set('tempdir',                     @_); }
sub species                     { return shift->_generic_get_or_set('species',                     @_); }
sub alignment_name              { return shift->_generic_get_or_set('alignment_name',              @_); }
sub num_records_in_current_file { return shift->_generic_get_or_set('num_records_in_current_file', @_); }
sub num_records_per_split_fastq { return shift->_generic_get_or_set('num_records_per_split_fastq', @_); }
sub current_file_number         { return shift->_generic_get_or_set('current_file_number',         @_); }
sub current_temp_dir            { return shift->_generic_get_or_set('current_temp_dir',            @_); }
sub current_file_name           { return shift->_generic_get_or_set('current_file_name',           @_); }
sub current_file_handle         { return shift->_generic_get_or_set('current_file_handle',         @_); }
sub chunk_created_callback      { return shift->_generic_get_or_set('chunk_created_callback',      @_); }

sub create_record_temp_dir {
  my $self = shift;
  
  my $record_tempdir = join '/',
    $self->tempdir,
    $self->species,
    $self->alignment_name,
    'partial_alignments',
    $self->current_file_number,
  ;

  return $record_tempdir;
}

sub create_chunk_file_basename {
  my $self = shift;
  my $file_basename = 'chunk_' . $self->current_file_number . '.fastq';
  return $file_basename;
}

sub init {
  my $self = shift;
  $self->num_records_in_current_file(0);
  $self->current_file_number(0);
  $self->update_chunk_dependent_variables;
  return;
}

sub update_chunk_dependent_variables {
  my $self = shift;
  
  $self->current_temp_dir(
    $self->create_record_temp_dir
  );
  $self->current_file_name(
    $self->create_chunk_file_basename
  );
  return;
}

sub absolute_chunk_file_name {
  my $self = shift;
  return $self->current_temp_dir . '/' . $self->current_file_name
}

sub create_new_file_handle {
  my $self = shift;
  
  my $record_tempdir = $self->current_temp_dir;

  use File::Path qw(make_path);
  make_path($record_tempdir);
  
  open my $fh, '>', $self->absolute_chunk_file_name;
  return $fh
}

sub process {
  my $self   = shift;
  my $record = shift;
  
  my $fh = $self->current_file_handle;
  
  use Scalar::Util qw ( openhandle );
  
  if (! defined $fh || ! openhandle($fh)) {
    $fh = $self->create_new_file_handle;
    $self->current_file_handle($fh);
  }
  
  foreach my $line (@$record) {
    $fh->print($line);
    $fh->print("\n");
  }
  
  $self->num_records_in_current_file(
    1 + $self->num_records_in_current_file
  );
  
  my $is_time_for_next_chunk 
    = $self->num_records_in_current_file 
      >= $self->num_records_per_split_fastq;
  
  if ($is_time_for_next_chunk) {
  
    $self->current_file_handle->close;
  
    $self->chunk_created_callback->(
      {
        absolute_chunk_file_name => $self->absolute_chunk_file_name,
        current_file_number      => $self->current_file_number,
        current_temp_dir         => $self->current_temp_dir,
      }
    );

    $self->current_file_number(1 + $self->current_file_number);
    $self->num_records_in_current_file(0);
    
    $self->update_chunk_dependent_variables;
  }
  return;
}

1;
