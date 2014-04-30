=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::SWEmbl

=head1 DESCRIPTION

Runs SWEmbl peak caller and optionally parses and processes output.
This is in the Hive namespace, but is not dependnant on any Hive modules and can 
be run as a stand alone job outside of the hive infrastructure.

=cut

package Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::SWEmbl;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use base qw( Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller );#Does not import

my %fswitches = 
 (
  bed   => '-B',
  sam   => '-S',
  maq   => '-M',
  eland => '-E',
  bam   => '-F',
 );

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  #specify -output_format default first, so it can be over-ridden?
  my $self   = $class->SUPER::new(-output_format => 'bed', -is_half_open => 1, @_);
  
  return $self;
}

#todo SWEmbl handles gzipped file, so just use is_gzipped and pass -z
#

#Can call this before new!
#So we don't have to mess about with passing params to run and load methods

#could return keys from fswitches?
#Preference should be done by the caller, not here, so order shouldn't matter
sub input_formats    { return ['bam', 'sam', 'bed']; }

sub requires_control { return 0; }

sub run {   
  my $self = shift;
  my ($align_file, $out_file, $suffix, $control_file, 
    $gzipped_align, $gzipped_control) = @{$self->file_info};
    
  my $format_switch = $fswitches{lc($suffix)};#auto-vivifies
  throw("$suffix is not recognised as a valid SWEmbl input format") if(! $format_switch);
        
        
  #Set -z 
  my $compressed = $gzipped_align ? ' -z ' : '';
        
  if($control_file &&
     ($gzipped_control != $gzipped_align)){
    throw("-i(input_file) and -r(eference/input file cannot have mixed compression states:\n\t".
      $align_file."\n\t".$control_file);  
  }
        
        
  my $cmd = $self->program_file." $format_switch -i $compressed ".
    $self->align_file.' '. $self->parameters.' -o '.$out_file;
  $cmd .= " -r ".$self->control_file if $self->control_file;
  
  warn "Running:\t$cmd\n" if $self->debug;
  
  
  #This did no cause failure when failed
  run_system_cmd($cmd);
     
  return;
}


1;

