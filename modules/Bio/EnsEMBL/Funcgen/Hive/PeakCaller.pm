=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::PeakCaller

=head1 DESCRIPTION

Base class for all peak caller analyses. Provides generic constructor method, 
along with file validation, accessor methods and a simple feature storing.

=cut

package Bio::EnsEMBL::Funcgen::Hive::PeakCaller;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped         gunzip_file 
                                               file_suffix_parse  open_file ); #convert_strand_from_bed


#Move this to utils

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  my ($prog_file, $prog_params, $align_file, $control_file, 
    $out_dir, $out_file_prefix, $reload, $rerun, $no_load, 
    $is_half_open, $convert_half_open, $output_format) =
    rearrange(['PROGRAM_FILE', 'PARAMETERS', 'ALIGN_FILE', 'CONTROL_FILE',
      'OUT_DIR', 'OUT_FILE_PREFIX', 'RELOAD', 'RERUN', 'NO_LOAD', 
      'IS_HALF_OPEN', 'CONVERT_HALF_OPEN', 'OUTPUT_FORMAT'], @_);
      
 
  #may want to have an empty string for the params
 
  if(! ($prog_file && (defined $prog_params) && $align_file &&
        $out_dir && $output_format)){
    #this will give an undef warning
    throw("Some mandatory parameters are not met:\n\t".
      "-PROGRAM_FILE => $prog_file, -ALIGN_FILE => $align_file, ".
      "-PARAMETERS => $prog_params, -OUT_DIR => $out_dir, ".
      "-OUTPUT_FORMAT => $output_format");
  }

  if(! -f $prog_file){
    throw("Program file does not exist or is not a file:\n$prog_file");  
  }


  if(! defined $control_file){
    warn "Running $prog_file without and control file\n"; #omit this warn? 
  }
 

  #pass helper here for debug/log?
  
  
  $self->{control_file}      = $control_file; 
  $self->{prog_file}         = $prog_file;
  $self->{align_file}        = $align_file;
  $self->{parameters}        = defined $prog_params  ? $prog_params : ''; #To avoid warnings
  $self->{out_dir}           = $out_dir;
  $self->{out_file_prefix}   = $out_file_prefix;
  $self->{output_format}     = $output_format;
  $self->{is_half_open}      = $is_half_open;
  $self->{convert_half_open} = $convert_half_open;
  
  #Init files here, so the output file is still available even if we don't run
  #i.e. for reloading
  
  
  #This is an issue as this forces the requirement of the feature files here
  #before we can ask this moduel what format it wants
  
  #change get_alignment_file method to take a list of formats in preference order?
  #This will allow the caller to pick the prefered format, but not fail if it is not available
  
  $self->get_file_info;
   
  return $self;
}

#subtlety here to the gzip status
#SWEmbl will handle both gzipd or both unzipd
#but not mixed gzip/unzipd

#todo handle $unzip, or pass back gzipped status
#what about no unzip defined unzipped files, do we then gzip?

#do we still want a reload mode, where we don't care about the input
#files, so long as we have a valid out_file_prefix

sub init_files {
  my ($self, $unzip) = @_;
  
  if(! defined $self->{file_info}){
  
    my ($align_file, $gzip_align) = gunzip_file($self->align_file);
    my $gzip_control = 0;
    my $control_file = $self->control_file;
  
    if($control_file){
      ($control_file, $gzip_control) = gunzip_file($self->align_file);
    }
  
    my ($file_prefix, $suffix) = file_suffix_parse($align_file);   
   
    if(! defined $self->out_file_prefix){
      $self->{out_file_prefix} = $file_prefix;  
    }
    
    $self->{out_file} = $self->out_dir.'/'.$self->out_file_prefix.'.'.$self->output_format;
  
    $self->{file_info} = [$align_file, $self->{out_file}, $suffix,
                          $control_file, $gzip_align, $gzip_control];
  }
    
  return $self->{file_info};
}

sub file_info         { return $_[0]->{file_info};         }
sub control_file      { return $_[0]->{control_file};      }
sub align_file        { return $_[0]->{align_file};        }
sub program_file      { return $_[0]->{program_file};      }
sub parameters        { return $_[0]->{parameters};        }
sub out_dir           { return $_[0]->{out_dir};           }
sub out_file_prefix   { return $_[0]->{out_file_prefix};   }
sub output_format     { return $_[0]->{output_format};     }
sub out_file          { return $_[0]->{out_file};          }
sub out_file_handle   { return $_[0]->{out_file_handle};   }
sub out_file_header   { return $_[0]->{out_file_header};   } 
sub is_half_open      { return $_[0]->{is_half_open};      }
sub convert_half_open { return $_[0]->{convert_half_open}; }   

  
#File parsing and DB loading may be polluting this PeakCaller code
#but is not mandatory functionality, so we will allow this as is quite useful
#This is completely agnostic toward the type of code ref p
 
sub process_features {
  my ($self, $process_cref, @process_cref_args) = @_;  
  
  assert_ref($process_cref, 'CODE', 'Process feature method');
  
  #parse/init methods can either be implemented directly in subclass as parse_feature
  #without setting output_format, or by setting output_format and relying on generic 
  #output format methods defined in this class, or by re-implementing the later in a subclass 
  my $init_method = 'init_output_file';
  my $parse_method = 'parse_feature';

  if(defined $self->output_format){
    $parse_method = 'parse_'.$self->output_format.'_feature';
    $init_method  = 'init_'.$self->output_format.'_feature';
  } 
  
  if(! $self->can($parse_method)){
    throw(ref($self)."::${parse_method} is not a valid method.".
    ' This class is missing a method or does not have a valid output format');      
  }
   
  if(! $self->can($init_method)){
    throw(ref($self)."::${init_method} is not a valid method.".
    ' This class is missing a method or does not have a valid output format');      
  }
  
  $self->$init_method;
  my $feature_hash;
 
  while( ($feature_hash = $self->parse_feature) &&
         defined $feature_hash ){          
    $process_cref->(@process_cref_args, $feature_hash);       
  }   
}  



sub init_bed_file {
  my $self = $_[0];
    
  $self->{out_file_handle} = open_file($self->out_file);
  #This will validate exists
  
  #check we have a header
  my $prevpos = $self->out_file_handle->getpos;
  my $line;
  
  while(($line = $self->out_file_handle->getline) &&
         defined $line){
    
    #Should we handle # and Region here too?
    
    if($line =~ /^browser/){
      $prevpos = $self->out_file_handle->getpos;
    }
    elsif($line =~ /^track/){
      $prevpos = $self->out_file_handle->getpos;
      $self->{out_file_header} = $line;
    }
    else{
      last;  
    }
  }
  
  $self->out_file_handle->setpos($prevpos);
  
  return $self->out_file_handle;
}


sub parse_bed_feature {
  my $self = $_[0];
  my $line;
  
  eval { $line = $self->out_file_handle->getline; };
  
  if($@){
    throw("Failed to getline from out_file_handle, maybe you need to init_bed_file first.\n$@");
  }
  
  if(defined $line){
    my ($seqid, $start, $end, undef, undef, undef, $score, undef, undef, $summit) = split(/\s+/, $line);
       
    #$strand = convert_strand_from_bed($strand);
    $summit = int($summit + 0.5);#Rounds to nearest integer as int rounds down     

    if($self->is_half_open && $self->convert_half_open){
      $start += 1;
    }
    
    #Handle half_open here, but slice in Caller which knows about the DB

    $line = 
      {-start      => $start,
       -end        => $end,
       -score      => $score,
       -summit     => $summit,
       -seq_region => $seqid   
      };
  }

  return $line;
}


1;

