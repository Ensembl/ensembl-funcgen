
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::RunWiggleTools

=head1 DESCRIPTION

A simple module to convert the bam files associated with a ResultSet into a BigWig

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools;

use warnings;
use strict;
use File::Temp                             qw( tempfile );
use Fatal                                  qw( close );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names run_system_cmd );
use Bio::EnsEMBL::Utils::Exception         qw( throw );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

# Todo
# 1 Move this to SeqTools so it can be re-used in other contexts?
#   Limited utility as it is already a very thin wrapper.
# 2 Add file version support? 

####################################################
## Preparations
####################################################


sub fetch_input {
  my $self = shift;
  $self->param('disconnect_if_idle', 1);  # Set before DB connection via SUPER::fetch_input
  $self->SUPER::fetch_input;
  $self->get_param_method('operator', 'silent', 'write');

  if($self->get_param_method('set_type', 'silent')){  # Set input params
    my $rset = $self->fetch_Set_input('ResultSet');
    $self->helper->debug(1, "RunWiggleTools::fetch_input got ResultSet:\t".$rset->name);

    # check we haven't sepcified input_files
    if($self->get_param_method('input_files', 'silent')){
      $self->throw_no_retry("It is unsafe to specify the input_files are set parameters.\n".
        'Please remove input_files from input_id');
    }

    if($rset->has_status('IMPORTED')){
      # This is not handling assembly specific IMPROTED states!
      # But we have moved to separate DBs for idfferent assemblies
      # So we need to drop that?


      if($self->param_silent('force')){ #force/recover
        $rset->adaptor->revoke_status('IMPORTED', $rset);
      }
      else{
        $self->throw_no_retry("Cannot write bigWig as result set is already marked as IMPORTED\n".
          'To over-ride, run reseed_jobs.pl -append \'{force>1}\'');
        # DebugJob -f will not do this, as this is a beekeeper option
      }
    }

    $self->input_files([$self->get_alignment_files_by_ResultSet_formats($rset, ['bam'])->{bam}]);
    $self->set_param_method('output_format', 'BigWig');
    $self->get_param_method('output_prefix',
                            'silent',
                            $self->db_output_dir.'/result_feature/'.$rset->name);
  }
  else{  # File input params
    my $inp_files = $self->get_param_method('input_files', 'required');  # Let's be forgiving

    if(! ref($inp_files)){
      $self->input_files([$inp_files]);
    }
    else{
      assert_ref($inp_files, 'ARRAY', 'input_files');
    }

    $self->get_param_method('output_format', 'silent', 'wig');  # default to wig
    $self->get_param_method('output_prefix', 'required');

    $self->helper->debug(1, "RunWiggleTools::fetch_input got input_files:\n\t".
      join("\n\t", @{$self->input_files}));
  }


  # TODO
  # Change this directory naming? There is not need for result_feature reference now?
  # Remember we can't change the structure on nfs live easily!
  # Is there any chance that get_alignment_files_by_ResultSet_formats
  # could trigger a filter/sort?
  # This file may already be in use by another analysis i.e. the peak caller

  # This may not actually be used if output_file is specified
  $self->get_output_work_dir_methods($self->db_output_dir.'/result_feature/');
  return;
}


####################################################
## Core function
####################################################

sub run {
  my $self        = shift;
  my @tmpfiles;
  my $chrlen_file = $self->write_chrlen_file;

  # DEFINE THE OUTPUT
  my $write_out = '-';
  my $reformat_cmd = '';
  my $output;

  if(lc($self->output_format) eq 'wig'){
    $write_out = $self->output_prefix.'.wig';
    $output    = $write_out;
  }
  elsif(lc($self->output_format) eq 'bigwig'){
    push @tmpfiles, $self->write_chrlen_file;
    $output       = $self->output_prefix.'.bw';
    $reformat_cmd = ' | wigToBigWig -fixedSummaries stdin '.$tmpfiles[0].' '.$output;
  }
  else{
    $self->throw_no_retry("Output format not supported:\t".$self->output_format);
    # Should really be done in fetch_input, but easier here.
  }

  my $operator = ($self->operator eq 'write') ? '' : $self->operator;

  # BUILD/RUN THE COMMAND
  # Presence of index files is input format specific, so not explicitly tested/generated here
  # Always use write to avoid redirects
  my $cmd = "wiggletools write $write_out $operator ".
   join(' ',@{$self->input_files}).$reformat_cmd;
  run_system_cmd($cmd);
  # Took about 20 mins with ~ 8GB mem for single bam

  # PUT OUT THE TRASH
  if(@tmpfiles){
    unlink @tmpfiles or
     warn "Could not unlinke tmp file:\t$!\n";
   }

  # UPDATE THE DB
  if($self->set_type){  # Update dbfile_registry.path and set ResultSet states
    my $rset = $self->ResultSet;
    $rset->adaptor->dbfile_data_root($self->db_output_dir);
    $rset->dbfile_path($output);
    $rset->adaptor->store_dbfile_path($rset);
    $rset->adaptor->set_imported_states_by_Set($rset);
  }

  return;
}


sub write_output {  return; }  # Nothing to do/flow here?


sub write_chrlen_file {
  my $self = shift;
  my ($fh, $name) = tempfile(DIR => $self->work_dir); #, UNLINK => 0); # is default with function interface

  my $slices = generate_slices_from_names($self->out_db->dnadb->get_SliceAdaptor,
                                          $self->slices,
                                          $self->skip_slices,
                                          'toplevel', 0, 1);  # nonref, incdups
  $self->debug(2, 'Writing lengths for '.scalar(@{$slices}).' toplevel (inc_dups) slices');

  foreach my $slice (@{$slices}) {
    print $fh join("\t", ($slice->seq_region_name, $slice->end - $slice->start + 1)) . "\n";
  }

  # Ideally need to validate this has been written correctly

  close $fh;
  return $name;
}


1;
