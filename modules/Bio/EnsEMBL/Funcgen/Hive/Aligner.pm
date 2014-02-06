=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::Aligner

=head1 DESCRIPTION

Base class for all sequence alignement analyses. Provides generic constructor method, 
along with file validation and some accessor methods.

=cut

package Bio::EnsEMBL::Funcgen::Hive::Aligner;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  my ($prog_file, $prog_params, $query_file, $ref_file, $out_dir, $format) =
    rearrange(['PROGRAM_FILE', 'PARAMETERS', 'QUERY_FILE', 
               'REFERENCE_FILE', 'OUTPUT_DIR', 'OUTPUT_FORMAT'], @_);
      
  if(! ($prog_file && $query_file && $ref_file)){
    throw("Some mandatory parameters are not met:\n\t".
      join("\n\t", ("-PROGRAM_FILE   => $prog_file,", 
                    "-QUERY_FILE     => $query_file,",
                    "-REFERENCE_FILE => $ref_file")));
  }

  if(! -f $prog_file){
    throw("Program file does not exist or is not a file:\n\t$prog_file");  
  }

  if(! -f $query_file){
    throw("Query file does not exist or is not a file:\n\t$query_file");  
  }

  if(! -f $ref_file){
    throw("Reference file does not exist or is not a file:\n\t$ref_file");  
  }

  (my $input_dir = $query_file) =~ s/(.*\/)[^\/].*/$1/go;
  $query_file =~ s/$input_dir//;

  if(defined $out_dir){ 
    if(! -d $out_dir){
      throw("Output directory does not exist:\n\t".$out_dir);  
    }
    
    $self->{output_dir}      = $out_dir;
  }
  else{
    $self->{output_dir}      = $input_dir;
  }
  
  #if($format){
  #  
  #  if(! grep(/^${format}$/, @{$self->valid_output_formats}) ){
  #    throw("Output format($format) is not valid. Please specify one of:\t".
  #      join("\t", @{$self->valid_output_formats}));  
  #  }
  #  
  #  #This assumes valid_output_formats are in order of general preference
  #  #i.e. most optimsed/performant/binary first
  #  $self->{output_format} = $self->valid_output_formats->[0]; 
  #  #}
  #else{
  #  $self->{output_format} = $self->valid_input_formats
  #}
  
  $self->{input_dir}    = $input_dir;
  $self->{program_file} = $prog_file;
  $self->{query_file}   = $query_file;
  $self->{parameters}   = defined $prog_params  ? $prog_params : ''; #To avoid warnings
  $self->{ref_file}     = $ref_file;
  
  return $self;
}

sub query_file     { return shift->{query_file};    }
sub reference_file { return shift->{ref_file};      }
sub program_file   { return shift->{program_file};  }
sub parameters     { return shift->{parameters};    }
sub output_dir     { return shift->{output_dir};    }
sub input_dir      { return shift->{input_dir};     }
#sub output_format  { return shift->{output_format}; }

#abstract methods which muct be defined in the the specific aligner sub class  
#sub run
#sub valid_output_formats

1;

