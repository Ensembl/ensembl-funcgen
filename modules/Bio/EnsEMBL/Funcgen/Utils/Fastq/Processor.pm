package Bio::EnsEMBL::Funcgen::Utils::Fastq::Processor;

use strict;
use Data::Dumper;

=head1

  A fastq record are four consecutive lines from a fastq file. 
  
  The FastqRecordProcessor expects a stream of arrays of four lines to the
  process method.
  
  These will be written into chunks. 
  
  When it has written $num_records_per_split_fastq into one file,
  it will:
    - Call $dataflow_fastq_chunk and
    - start a new chunk file.
  
  For every new chunk file it will create a new directory in $tempdir. It will use
  $alignment_name in the temporary directory it creates so that different 
  alignments don't overwrite each other.
  
  Once it
  
  It will create directories and chunks in these directories.

=head2 Example
  
  use Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser;
  use Bio::EnsEMBL::Funcgen::Utils::Fastq::Processor;

  my $tempdir  = '/hps/nobackup/production/ensembl/mnuhn/peak_calling_paired_end/temp_dir/test.deleteme';
  my $species  = 'mus_musculus';
  my $alignment_name = 'test_alignment';
  my $num_records_per_split_fastq = 100;
  my $dataflow_fastq_chunk = sub { 
    my $param = shift;
    print Dumper($param);
  };

  my $parser = Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser->new;

  #   /nfs/production/panda/ensembl/funcgen/source_files/fastq/mus_musculus/encode/v92/ENCFF664IWE.fastq.gz
  #   /nfs/production/panda/ensembl/funcgen/source_files/fastq/mus_musculus/encode/v92/ENCFF228QFU.fastq.gz

  my $fastq_1_file = '/nfs/production/panda/ensembl/funcgen/source_files/fastq/mus_musculus/encode/v92/ENCFF664IWE.fastq.gz';
  my $fastq_2_file = '/nfs/production/panda/ensembl/funcgen/source_files/fastq/mus_musculus/encode/v92/ENCFF228QFU.fastq.gz';

  my $cmd_1 = "zcat $fastq_1_file |";
  my $cmd_2 = "zcat $fastq_2_file |";

  my $fastq_record_processor = Bio::EnsEMBL::Funcgen::Utils::Fastq::Processor->new(
    -tempdir                     => $tempdir,
    -species                     => $species,
    -alignment_name              => $alignment_name,
    -num_records_per_split_fastq => $num_records_per_split_fastq,
    -chunk_created_callback      => $dataflow_fastq_chunk,
    -tags                        => [ 1, 2 ]
  );

  open my $fh_1, $cmd_1 or die("Can't execute $cmd_1");
  open my $fh_2, $cmd_2 or die("Can't execute $cmd_2");

  my $max   = 500;
  my $count =   0;

  while (
      (my $fastq_1_record = $parser->parse_next_record($fh_1)) 
      && 
      (my $fastq_2_record = $parser->parse_next_record($fh_2)) 
      && 
      ($count<$max)
    ) {

    $fastq_record_processor->process($fastq_1_record, $fastq_2_record);
    $count++;
  }

=cut

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
    tags                        => 'tags',
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
sub tags                        { return shift->_generic_get_or_set('tags',                        @_); }

sub create_record_temp_dir {
  my $self = shift;
  
  my $record_tempdir = join '/',
    $self->tempdir,
    $self->alignment_name,
    'partial_alignments',
    $self->current_file_number,
  ;
  return $record_tempdir;
}

sub create_chunk_file_basename {
  my $self = shift;
  my $tags = $self->tags;
  return [ map { 'chunk_number_' . $self->current_file_number . '_tag_' . $_ . '.fastq' } @$tags ];
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
  return [ map { $self->current_temp_dir . '/' . $_ } @{$self->current_file_name} ];
}

sub create_new_file_handle {
  my $self = shift;
  
  my $record_tempdir = $self->current_temp_dir;

  use File::Path qw(make_path);
  make_path($record_tempdir);
  
  my $absolute_chunk_file_names = $self->absolute_chunk_file_name;
  
  my @fhs;
  for my $absolute_chunk_file_name (@$absolute_chunk_file_names) {
    open my $fh, '>', $absolute_chunk_file_name || die("Can't open file $absolute_chunk_file_name !");
    push @fhs, $fh;
  }
  
  return \@fhs;
}

sub flush {
  my $self   = shift;

  my $fhs = $self->current_file_handle;
  
  use Scalar::Util qw ( openhandle );
  
  if (! (defined $fhs->[0] && openhandle($fhs->[0]))) {
    return;
  }
  
  $self->chunk_created_callback->(
    {
      absolute_chunk_file_name => $self->absolute_chunk_file_name,
      current_file_number      => $self->current_file_number,
      current_temp_dir         => $self->current_temp_dir,
    }
  );
  foreach my $fh (@$fhs) {
    $fh->close;
  }
  return;
}

sub process {
  my $self   = shift;
  my @record = @_;
  
  my $fhs = $self->current_file_handle;
  
  use Scalar::Util qw ( openhandle );
  
  if (! defined $fhs->[0] || ! openhandle($fhs->[0])) {
    $fhs = $self->create_new_file_handle;
    $self->current_file_handle($fhs);
  }
  
  my $tags = $self->tags;
  
  for (my $tag_index = 0; $tag_index < scalar @$tags; $tag_index++) {
  
    my $current_record = $record[$tag_index];
    my $current_fh     = $fhs->[$tag_index];
    
    foreach my $line (@$current_record) {
      $current_fh->print($line);
      $current_fh->print("\n");
    }
  }
  
  $self->num_records_in_current_file(
    1 + $self->num_records_in_current_file
  );
  
  my $is_time_for_next_chunk 
    = $self->num_records_in_current_file 
      >= $self->num_records_per_split_fastq;
  
  if ($is_time_for_next_chunk) {
    
    foreach my $fh (@$fhs) {
      $fh->close;
    }
  
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
