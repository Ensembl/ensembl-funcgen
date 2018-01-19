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

=head1 NAME



=head1 SYNOPSIS

  

=head1 PARAMETERS

  Mandatory:
    -sam_header      <String>
      
  Optional:
    -output_dir <String>  Output directory path
    -help                 Prints the POD SYNOPSIS 
    -man                  Prints this full man page
 
=head1 DESCRIPTION


=cut


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

main() unless(caller);
#This avoids executing run_aligner when including as a libary via require or do
#for general usage or when unit testing


sub main{
  my ($sam_header, $out_file);
  my @tmp_args = @ARGV;

  GetOptions 
   ('sam_header=s' => \$sam_header, 
    'out_file=s'   => \$out_file,
     
    #Optional
    'help'      => sub { pod2usage(-exitval => 0); }, 
    'man|m'     => sub { pod2usage(-exitval => 0, -verbose => 2); },
   ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
  
  if(! (defined $sam_header && -e $sam_header)){
    pod2usage(-exitval => 1, 
              -message => "Must specify a sam_header file parameter");  
  }
  
  if(! defined $out_file){
    ($out_file = $sam_header) =~ s/\.sam//;
    $out_file  =~ s/\.header//;
    $out_file .= '.CCAT_chr_lengths.txt';
    warn "No -out_file paramter specified\n";
  }

  print "Writing output to:\n\t$out_file\n";
  
  

  #Get the size file... similar to the sam header... 
  open(my $sh_fh, '<',   $sam_header);
  open(my $out_fh , '>', $out_file);

  while(<$sh_fh>){
    chomp;
    if(/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/){
      my $slice = $1;
      my $slice_size = $2;
    
      #Mouse Hack!!
      #next if(($self->_species eq 'mus_musculus') && !($slice =~ /chromosome/));
      print $out_fh $slice."\t".$slice_size."\n";
    }else{
      die("Could not process sam header line $_"); 
    }
  }
  
  close($sh_fh);
  close($out_fh);
  print "DONE\n";
}




1;


__END__

