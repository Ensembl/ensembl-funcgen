
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::RunWiggleTools

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd );
# use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( write_chr_length_file);
# use Bio::EnsEMBL::Utils::Exception         qw( throw );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

sub ResultSet {
    my $self = shift;
    $self->{'_ResultSet'} = shift if @_;
    return $self->{'_ResultSet'};
}


sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;

  if($self->get_param_method('set_type', 'silent')) {
  
    my $db            = $self->param_required('out_db');
    my $type          = $self->param('type');
    my $result_set_id = $self->param('result_set_id');
    
    my $result_set_adaptor = $db->get_ResultSetAdaptor;
    # HACK
    $result_set_adaptor->{file_type} = 'BAM';
# Oh the tediousness of the parameters!
# 
# We need a result set.
# 
# Ideally we have a result set id:
    
    if (defined $result_set_id) {
      my $result_set = $result_set_adaptor->fetch_by_dbID($result_set_id);
      $self->ResultSet($result_set);
    }

# No? Then do we have at least a dbID in this weird parameter 
# passing mess?

    if (! defined $result_set_id && defined $self->param('dbID')) {

# If is is a data set dbID, then we can fetch the result set id from that.

      if ($self->param('set_type') eq 'DataSet') {
      
        my $dbID = $self->param('dbID');
        my $dataset_adaptor = $db->get_DataSetAdaptor;
        my $dataset = $dataset_adaptor->fetch_by_dbID($dbID);
        
        my $rsets = $dataset->get_supporting_sets('result');
        if(scalar @$rsets != 1) {
            throw("Expected 1 ResultSet, found:\n".$self->helper->dump(\@$rsets));
        }
        $self->ResultSet($rsets->[0]);

#         use Data::Dumper;
#         print Dumper($rsets->[0]);
#         die;
      
      }
    
    }
  my $result_set = $self->ResultSet;
#     my $result_set = $self->fetch_Set_input('ResultSet');
    $self->helper->debug(1, "RunWiggleTools::fetch_input got ResultSet:\t".$result_set->name);

    if($self->get_param_method('input_files', 'silent')) {
      $self->throw_no_retry("It is unsafe to specify the input_files are set parameters.\n".
      'Please remove input_files from input_id');
    }
    my $is_control;
    if ($type eq 'control') {
      $is_control = 1;
    }
    $self->input_files([
      $self->db_output_dir
      . '/' . $self->get_alignment_files_by_ResultSet_formats($result_set, $is_control)
    ]);

    my $output_prefix = $self->bigwig_output_dir . '/' . $result_set->name;
    
    warn "output_prefix = $output_prefix";
    
    $self->get_param_method(
      'output_prefix',
      'silent',
      $output_prefix
    );
    $self->hive_run_system_cmd("mkdir -p " . $self->bigwig_output_dir);
  }
  return;
}


####################################################
## Core function
####################################################

sub run {
  my $self       = shift;
  my $result_set = $self->ResultSet;
  
  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $species          = $self->param('species');
  my $assembly         = $self->param('assembly');
  my $epigenome_gender = $result_set->epigenome->gender;

  my $chromosome_lengths_relative = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $epigenome_gender,
    assembly         => $assembly,
    file_type        => 'chromosome_lengths_by_species_assembly',
  });
  my $reference_data_root_dir = $self->param('reference_data_root_dir');
  
  my $chromosome_length_file = $reference_data_root_dir . '/' . $chromosome_lengths_relative;
  
  my $output = $self->output_prefix.'.bw';

  $self->throw_no_retry("bigWigToWig can't cope with colons, these should not be in the experiment names! ($output)")
    if ($output =~ /:/);

  my $cmd_wiggletools = "wiggletools write - ". $self->_build_rpkm_cmd;
  my $cmd_wigToBigWig = 'wigToBigWig -fixedSummaries stdin ' . $chromosome_length_file . ' ' . $output;
  
  my $cmd = $cmd_wiggletools . ' | ' . $cmd_wigToBigWig;
  
  my $cmd = "wiggletools write - ". $self->_build_rpkm_cmd . ' | '.'wigToBigWig -fixedSummaries stdin ' . $chromosome_length_file . ' ' . $output;

  $self->helper->debug(1, "Running:\n\t".$cmd);
  eval {
    $self->hive_run_system_cmd($cmd);
  };
  if ($@) {
    sleep(20);
    # If job was killed for memlimit, allow for the worker to be killed as well.
    die($@);
  }
#   my $temporary_directory_root = $self->temporary_directory_root;
#   
#   my $hostname       = `hostname`;
#   chomp($hostname);
#   my $pid            = $$;
#   my $temporary_file = $hostname . '.' . $pid . '.wig';
#   
#   my $temporary_directory = File::Spec->catfile($temporary_directory_root, 'wiggletools');
#   
#   use File::Path qw(make_path remove_tree);
#   make_path($temporary_directory);
#   
#   my $temporary_file_full_path = "$temporary_directory/$temporary_file";
#   
#   my $cmd_wiggletools_redirected = $cmd_wiggletools . " > " . $temporary_file_full_path;
#   
#   $self->hive_run_system_cmd("rm -f $temporary_file_full_path", undef, 1);
#   $self->hive_run_system_cmd($cmd_wiggletools_redirected, undef, 1);
#   $self->hive_run_system_cmd(
#     "wigToBigWig -fixedSummaries $temporary_file_full_path " . $chromosome_length_file . ' ' . $output, 
#     undef, 1
#   );  
#   $self->hive_run_system_cmd("rm -f $temporary_file_full_path", undef, 1);
#   
#   # If job was killed for memlimit, allow for the worker to be killed as well.
#   sleep 20;

  $result_set->adaptor->dbfile_data_root($self->db_output_dir);
  $result_set->dbfile_path($output);
  $result_set->adaptor->store_dbfile_path($result_set, 'BIGWIG');

  return;
}

sub _build_rpkm_cmd {
  my $self = shift;

  # This assumes bam input with index
  my @one_bam = @{$self->input_files};
  die("There can only ever be one bam file associated with a result set!") unless(@one_bam == 1);
  my $bam_file = $one_bam[0];

  my $cmd = qq(samtools index $bam_file);
  $self->hive_run_system_cmd($cmd, undef, 1);

  my $count_cmd = 'samtools idxstats '.$bam_file.' | awk \'{total = total + $2} END{print total}\'';
  $self->helper->debug(2, "Running:\n\t".$count_cmd);
  my $total_mapped = run_backtick_cmd($count_cmd);

  if(! $total_mapped) {
    $self->throw_no_retry("Failed to get number of mapped reads from index of:\n\t".$bam_file);
  }
  my $wiggle_tools_cmd = 'mean scale '.(10**9 / $total_mapped).' '.$bam_file;

  return $wiggle_tools_cmd;
}


1;
