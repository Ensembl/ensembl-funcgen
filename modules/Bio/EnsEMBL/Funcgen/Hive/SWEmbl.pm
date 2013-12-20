=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::SWEmbl

=head1 DESCRIPTION

Runs SWEmbl peak caller and optionally parses and processes output.
This is in the Hive namespace, but is not dependnant on any Hive modules and can 
be run as a stand alone job outside of the hive infrastructure.

=cut

package Bio::EnsEMBL::Funcgen::Hive::SWEmbl;

use warnings;
use strict;

#qw the methods even if they are EXPORTED, so we know where they come from
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use base qw( Bio::EnsEMBL::Funcgen::Hive::PeakCaller );#Does not import

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
  my $self   = $class->SUPER::new(-output_format => 'bed', @_);
  #specify -output_format default first, so it can be over-ridden?
  
  #bed is half_open by default
  if(! defined $self->is_half_open){
    $self->{is_half_open} = 1;
  }
 
  return;
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
  
  my ($align_file, $suffix, $out_file, $control_file, 
    $gzip_align, $gzip_control) = @{$self->file_info};

    
  my $format_switch = $fswitches{lc($suffix)};#auto-vivifies
  throw("$suffix is not recognised as a valid SWEmbl input format") if(! $format_switch);
        
  my $command = $self->program_file." $format_switch -z -i ".
    $self->align_file.' '. $self->parameters.' -o '.$out_file;
  
  if($self->control_file){  
    $command .= " -r ".$self->control_file; 
  }
  
  #warn "Running analysis:\t$command";
  if(system($command)) { throw("FAILED to run $command"); }
     
  return 1;
}


1;

