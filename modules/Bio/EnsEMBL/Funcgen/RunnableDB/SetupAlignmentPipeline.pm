=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::SetupAlignmentPipeline

=head1 DESCRIPTION

'SetupAlignmentPipeline' Does all the setup before the Alignment is run
Checks for existence of input files, etc...
This Runnable CAN be run multiple times in parallell!

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::SetupAlignmentPipeline;

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(is_gzipped);
#use Data::Dumper;

use base ('Bio::EnsEMBL::Funcgen::RunnableDB::Alignment');

#TODO... Maybe use and update the tracking database...
sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  
  $self->SUPER::fetch_input();

  #Magic default number...
  my $fastq_chunk_size = 8000000;
  if($self->param("fastq_chunk_size")){ $fastq_chunk_size = $self->param("fastq_chunk_size")};
  $self->_fastq_chunk_size($fastq_chunk_size);

  #Sets up the output dir for this experiment_name
  my $output_dir = $self->_output_dir();
  if(! -d $output_dir){ 
    system("mkdir -p $output_dir") && throw("Couldn't create output directory $output_dir");
  }

  my $input_dir = $self->_input_dir();
  if(! -d $input_dir ){ throw " Couldn't find input directory $input_dir"; }

  opendir(DIR,$input_dir);
  my @dirs = grep(/^\d/,readdir(DIR));
  closedir(DIR);

  if(scalar(@dirs)==0){ throw "No replicates found in $input_dir"; } 

  my @replicates;
  foreach my $dir (@dirs){
    #TODO: maybe use some other code for replicates? (e.g. Rep\d )
    if($dir =~ /^(\d)$/){
      my $replicate = $1;

      opendir(DIR,$input_dir."/".$replicate);
      my @files = grep(/.fastq/,readdir(DIR));
      closedir(DIR);

      my @replicate_files;
      foreach my $file (@files){
	#TODO: Check if there are "true" sequence files inside, not just by filename
	if($file =~ /.fastq/){
	  #TODO: Maybe use e.g. Bio::EnsEMBL::Funcgen::Utils::EFGUtils::is_gzipped
	  if($file =~ /^(.*.fastq).gz$/){
	    my $cmd = "gunzip ".$input_dir."/".$replicate."/".$file;
	    $file = $1;    
	    if(system($cmd) != 0){ warn "Problem occurred unzipping $file"; next; }
	  }
	  if($file =~ /^(.*.fastq).bz2$/){
	    my $cmd = "bunzip2 ".$input_dir."/".$replicate."/".$file;
	    $file = $1;    
	    if(system($cmd) != 0){ warn "Problem occurred unzipping $file"; next; }
	  }

	  push(@replicate_files,$file);
	} else { warn "File $file does not seem fastq"; }
      }

      if(scalar(@replicate_files)==0){ throw "No files for replicate $replicate"; } 
      
      push(@replicates, $replicate);

    } else { warn "Invalid replicate $dir ignored";   }
  }

  $self->_replicates(\@replicates);

  return 1;
}

sub run {   
  my $self = shift @_;

  my $fastq_chunk_size = $self->_fastq_chunk_size();

  my %output_ids;
  my $set_name = $self->_set_name();
  foreach my $replicate (@{$self->_replicates()}){
    my $rep_dir = $self->_input_dir()."/".$replicate;
    #TODO: maybe pass the number of lines as parameter...
    
    my $split_cmd="cat ".$rep_dir."/*.fastq | split -d -a 4 -l $fastq_chunk_size - ".
      $self->_output_dir()."/".$set_name.".".$replicate.".";
    if(system($split_cmd) != 0){ throw "Problems running $split_cmd";  }
    
    #Need to remove tmp files from previous runs
    #system('rm -f '.$self->_output_dir()."/${set_name}\.${replicate}\.".'*');

    opendir(DIR,$self->_output_dir());
    my @files = grep(/${set_name}\.${replicate}\./,readdir(DIR));
    closedir(DIR);
    
    my @replicate_input_ids;
    foreach my $file (@files){
      #Need to add the specific file to the input_id...
      my $new_input_id = eval($self->input_id);
      $new_input_id->{"input_file"} = $file;
      push(@replicate_input_ids, $new_input_id);
    }
    $output_ids{$replicate} = \@replicate_input_ids;

    #Rezip the original files
    my $zip_cmd = "gzip ".$rep_dir."/*.fastq";
    if(system($zip_cmd) != 0){ warn "There was a problem zipping back the source files";  }

  }

  $self->_output_ids(\%output_ids);

  return 1;
}


sub write_output {  # Create the relevant job
  my $self = shift @_;

  my %output_ids = %{$self->_output_ids()};
  #TODO add the number of replicates as input so it will behave in a proper fashion!!
  my $new_input_id = eval($self->input_id);
  $new_input_id->{"nbr_replicates"} = scalar(keys %output_ids);
  my ($converge_job_id) = @{ $self->dataflow_output_id($new_input_id, 3, { -semaphore_count => scalar(keys %output_ids) })  };

  foreach my $replicate (keys %output_ids){
    my $new_input_id = eval($self->input_id);
    $new_input_id->{"replicate"} = $replicate;
    my $rep_out_ids = $output_ids{$replicate};
    $new_input_id->{"nbr_subfiles"} = scalar(@$rep_out_ids);
    #Carefull with the id: 1-run_alignmens; 2-wrap_up; 3-converge
    #Add the number of files in each so it will behave appropriately...
    my ($funnel_job_id) = @{ $self->dataflow_output_id($new_input_id, 2, { -semaphore_count => scalar(@$rep_out_ids),  -semaphored_job_id => $converge_job_id })  };
    my $fan_job_ids = $self->dataflow_output_id($rep_out_ids, 1, { -semaphored_job_id => $funnel_job_id } );
  }
  
  return 1;

}

#Private getter / setter to the fastq chunk size
sub _fastq_chunk_size {
  return $_[0]->_getter_setter('fastq_chunk_size',$_[1]);
}

#Private getter / setter to the output_ids list
sub _output_ids {
  return $_[0]->_getter_setter('output_ids',$_[1]);
}

#Private getter / setter to the output_ids list
sub _replicates {
  return $_[0]->_getter_setter('replicates',$_[1]);
}

1;
