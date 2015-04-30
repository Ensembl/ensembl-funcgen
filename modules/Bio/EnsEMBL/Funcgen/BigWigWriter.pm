=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::CollectionWriter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::BigWigWriter;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseImporter');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names );
use Bio::EnsEMBL::Utils::Exception qw(throw);
use File::Temp;

$main::_debug_level = 0;
$main::_tee = 0;
$main::_no_log = 1;

####################################################
## Preparations
####################################################
sub fetch_input {  
  my $self = shift;
  my $rset = $self->fetch_Set_input('ResultSet');
  #Set some module defaults
  $self->param('disconnect_if_idle', 1);

  $self->SUPER::fetch_input;    
  $self->helper->debug(1, "CollectionWriter::fetch_input after SUPER::fetch_input");  
  $self->helper->debug(1, "CollectionWriter::fetch_input got ResultSet:\t".$rset);
  
  $self->init_branching_by_analysis;  #Set up the branch config
  
  $self->get_output_work_dir_methods($self->db_output_dir.'/result_feature/');
  $self->set_param_method('bam_files', $self->get_alignment_files_by_ResultSet_formats($rset, ['bam']), 'required'); 
  $self->set_param_method('rset_name', $rset->name, 'required'); 
  $self->set_param_method('chrom_lengths', $self->create_chrom_lengths, 'required'); 
  
  #todo why do we have to set EFG_DATA?
  $ENV{EFG_DATA} = $self->output_dir;
}

####################################################
## Core function 
####################################################
sub run {
  my $self = shift;
  my @bams =  @{$self->bam_files};
  my $chrom_lengths = $self->get_chrom_lengths;
  my $output = $self->output_dir.$self->get_rset_name.'.bw';
  # The fixedSummaries within the wigToBigWig code presets the zoom levels within the BigWig file to
  # fit the Ensembl default levels

  # NOTE: We are creating a single wig file and converting it to BigWig
  # This may cost a lot of memory at the wigToBigWig stage.
  # If the problem arises, we can create one bigWig per slice then merge them with bigWigCat
  my $cmd = "wiggletools write - sum " . join(" ", @bams) . " | wigToBigWig -fixedSummaries stdin $chrom_lengths $output";
  system($cmd) || die("Failed when running command:\n$cmd\n");
  unlink $chrom_lengths;
}

####################################################
## Cleaning up
####################################################
sub write_output { 
  my $self = shift;    
  $self->dataflow_job_groups;
  unlink $self->get_chrom_lengths;
  return; 
}

####################################################
## Dumping chromosome lengths into text file 
####################################################
sub create_chrom_lengths {
  my ($self) = @_;
  my ($fh, $name) = tempfile(DIR => $self->work_dir, UNLINK => 0);
  $self->fetch_chrom_lengths($fh);
  close $fh;
  return $name;
}

sub fetch_chrom_lengths {
  my ($self, $fh) = @_;
  print_log("Fetching chromosome lengths from core DB\n");
  my $slice_adaptor = $self->out_db->get_SliceAdaptor();
  my @slices = @{ $slice_adaptor->fetch_all('toplevel', undef, undef, 0) };

  foreach my $slice (@slices) {
    print $fh join("\t", ($slice->seq_region_name(), $slice->end() - $slice->start() + 1)) . "\n";
  }
}

1;
