#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 DESCRIPTION

Little utility for filtering PWM mappings using the output of pwm_filter_mappings.pl

Input comprises a thresholds file as produced by ensembl-funcgen/scripts/pwm_filter_mappings.pl and a number of original mappings files.

It expects to find the original mappings in the current working directory in files with names like jaspar.version.mappings eg MA0012.1.mappings. Typically these will have been created by using ensembl-funcgen/scripts/pwm_genome_map.pl followed by ensembl-funcgen/scripts/pwm_filter_mappings.pl.

Output files have names like MA0012.1.pwm_map and are written to the directory specified on the command line with -o




=head1 USAGE

The output directory must already exist, this script will not create it.

=head1 EXAMPLES

pwm_thresh_filter.pl -i thresholds  -o filtered/

=head1 SEE ALSO

ensembl-funcgen/scripts/pwm_filter_mappings.pl
ensembl-funcgen/scripts/pwm_genome_map.pl	   
ensembl-funcgen/scripts/pwm_mapping_notes

=cut

#TODO This is much quicker in awk! Currently done in parallelised in a bash script


use strict;
use DBI;
use Env;

use Getopt::Std;
use IO::Handle;
use IO::File;


#$| = 1; #no output buffer

my $out_dir='';
my $infile='';
my $verbose = 0;

my %opt;

if ($ARGV[0]){
&Getopt::Std::getopts('h:o:i:', \%opt) || die ;
}else{
&help_text; 
}

&process_arguments;

unless( -e $out_dir ){die "output directory $out_dir does not exist. Please create it"}


my %facs;
open(IN,$infile) or die "failed to open $infile";
while( my $line = <IN>){
    my @field = split("\t",$line);
    $facs{$field[0]}->{'factor'} = $field[1];
    $facs{$field[0]}->{'thresh'} = $field[2];
    $facs{$field[0]}->{'score_thresh'} = $field[3];
}
close(IN);


foreach my $mat (keys %facs){

  my $source_file = $mat.".mappings";
  &commentary("filtering $source_file\n");

  my $write_file = $out_dir.$mat.".pwm_map";
#    $verbose = 2;

    #TODO This is already calculated in pwm_filter_map.pl
    #This should be written to file/tracking/hive, so it doesn't have to be recalculated
    #This should go in a binding_matrix_analysis table
    #which will contain the max score, and the perc threshold used, and maybe some other stats?
    #in fact we already have the score_thresh, why is this being recalculated?
    
  my $max = &backtick("cut -f 5 $source_file | sort -g -u | tail -1");
  chop $max;
  my $thresh = $facs{$mat}->{'thresh'} * $max/100;
  print "$mat $max $thresh\n" if $verbose;

  my $ifh;
  open($ifh,$source_file) or die "failed to open $source_file";
  my $ofh;
  open($ofh,"> $write_file") or die "failed to open $write_file";
  
  while(my $line = <$ifh>){
    my @field = split("\t",$line);
    
    if($field[4] >= $thresh){
      print $ofh $line or die "failed to print to $write_file";
    }
  }
	    
  close($ifh);
  close($ofh);
}

exit;


 
###################################################################

sub backtick{
    my $command = shift;

    warn "executing $command \n" if ($verbose ==2);

    my $res = `$command`;
    if($?){
        warn "failed to execute $command \n";
        warn "output :-\n $res \n";
	die "exit code $?";
    }

    return $res;
}


# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print  "$_[0]";
}



   
sub err{
    print STDERR "$_[0]\n";
}
  


sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }


    if (exists $opt{o}){
        $out_dir = $opt{o}.'/';
    } else{
	&help_text(" Please supply the name of a directory for the output");
    } 


    if (exists $opt{i}){
        $infile = $opt{i};
    }else{
	&help_text(" Please supply the name of a pwm_filter_mappings.pl output file ");
    }


} 


sub help_text{
    my $msg=shift;
    
    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    pwm_thresh_filter.pl -i <file>

                  [-h] for help

                  [-o] <output directory> - name of a directory for output
                   -i  <input file> - name of a  pwm_filter_mappings.pl output file

END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
