=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::PeakCaller

=head1 DESCRIPTION

Base class for all peak caller analyses. Provides generic constructor method, 
along with file validation, accessor methods and a simple feature storing.

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

  my ($prog_file, $prog_params, $query_file, $ref_file, $out_dir) =
    rearrange(['PROGRAM_FILE', 'PARAMETERS', 'QUERY_FILE', 
               'REFERENCE_FILE', 'OUT_DIR'], @_);
      
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

  $self->{prog_file}  = $prog_file;
  $self->{query_file} = $query_file;
  $self->{parameters} = defined $prog_params  ? $prog_params : ''; #To avoid warnings
  $self->{out_dir}    = $out_dir;
  $self->{ref_file}   = $ref_file;
  
  return $self;
}

sub query_file     { return shift->{query_file};   }
sub reference_file { return shift->{ref_file};     }
sub program_file   { return shift->{program_file}; }
sub parameters     { return shift->{parameters};   }
sub out_dir        { return shift->{out_dir};      }

#abstract class which muct be defined in the the specific aligner sub class  
#sub run{}

1;

