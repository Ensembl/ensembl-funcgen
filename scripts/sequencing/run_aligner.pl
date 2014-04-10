#!/usr/bin/env perl


#This is set up such that it can be required by a module and run as a library 
#without having to fork a system cmd
#However, it is not used in the hive pipeline, as it is easier just to 
#create and run the runnable directly
#This is more useful if you want to run an alignment where you don't have access
#to a DB and the relevant Analysis object.

#ALIGNER: ['PROGRAM_FILE', 'PARAMETERS', 'QUERY_FILE', 
#               'TARGET_FILE', 'OUTPUT_DIR', 'OUTPUT_FORMAT', 'DEBUG']
#BWA: 'SAM_REF_FAI'


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
Getopt::Long::Configure("pass_through"); #Allows unknown options to pass on to @ARGV

&run_aligner();


#renamed run_aligner for clarity when required into a module

sub run_aligner{
  #Assigning to @ARGV allows this script to be required and run from a module 
  #without having to fork via and external system call
  @ARGV = @_ if (! @ARGV) && @_;
  my $aln_pkg;

  GetOptions 
   (#Mandatory
    'aligner=s' => \$aln_pkg,  
    #Optional  
    'help'             => sub { pod2usage(-exitval => 0); }, 
    'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
   );

  my @aligner_params = @ARGV;

#can't do this as we expect the unexpected! <Cue silouette of crazy naked dancing lady> <- 10 points for that reference
# or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
             
 
  if(! defined $aln_pkg){
    pod2usage(-exitval => 1, -message => '-aligner parameter must be defined');
  }
  
  if($aln_pkg !~ /::/){
    warn "Full $aln_pkg namespace was not specified, defaulting to:\tBio::EnsEMBL::Funcgen::Hive::$aln_pkg\n";
    $aln_pkg = "Bio::EnsEMBL::Funcgen::Hive::$aln_pkg";
  }
  
  #Quote so eval treats $aln_pkg as BAREWORD and converts :: to /
  if( ! eval "{require $aln_pkg; 1}"){
    die("Failed to require $aln_pkg\n$@");
  }

  if(! @aligner_params){
    pod2usage(-exitval => 1, 
              -message => "There are no aligner parameters specified\n".
              "For more information please run:\tperldoc $aln_pkg");
  }
 
  my $aligner = $aln_pkg->new(@aligner_params);

  if(! defined $aligner){
    die("Failed to create $aln_pkg from params:\n@aligner_params");
  }

  $aligner->run;

  return;
}

1;


__END__

