package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::MergeBamFiles;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $chunks        = $self->param_required('chunks');
  my $data_root_dir = $self->param_required('data_root_dir');
  my $species       = $self->param_required('species');
  my $plan          = $self->param_required('plan');
  
  my $align_plan = $plan
    ->{remove_duplicates}
    ->{align}
  ;
  my $assembly       = $align_plan->{to_assembly};
  my $bam_file_real  = $align_plan->{bam_file}->{real};

  use File::Basename qw( dirname basename );
  my $dirname = dirname($bam_file_real);
  my $full_path = $data_root_dir . '/' . $dirname;

  use File::Path qw(make_path remove_tree);
  make_path($full_path);
  
  use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# merge_bams 
  merge_bams({
    input_bams => $chunks,
    output_bam => $data_root_dir . '/' . $bam_file_real,
    debug => $self->debug,
  });
  return;
}

1;
