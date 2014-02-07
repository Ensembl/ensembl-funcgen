#!/usr/bin/env/perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;
use Getopt::Long;
##########################INPUT#################
#$pad_start=90745845;
#$pad_end=91644698;
#$chrY="/lustre/scratch109/ensembl/funcgen/ia3/GRCm38/fasta/mus_musculus/GRCm38/fasta/chrY.fa";

my $pad_start='';
my $pad_end='';
my $chr='';
my $help='';
my $assembly='';
usage() if ( @ARGV < 1 or
          ! GetOptions('help|?' => \$help, 'pad_start=i' => \$pad_start, 'pad_end=i' => \$pad_end, 'chromosome' => \$chr, 'assembly' => \$assembly)
          or defined $help );
 
sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: chrY_pad_par.pl [--pad_start start] [--pad_end end] [--chromosome fastafile] [--assembly assembly_name] [--help|-?]\n";
  exit;
}



#////////////////////////////////////////////////

#PAR is inclusive both start and end are included
my $pad_length=($pad_end - $pad_start +1);

#$pad = 'N' x 898854;
my $pad = 'N' x $pad_length;

$pad=[split(//,$pad)];

open (FH, $chr);                                                                                                                                                                       
my @arr = <FH>; 
close FH;
shift @arr;                                                                                                                                                                                 
chomp(@arr);                                                                                                                                                                                
my $line=join("",@arr);
@arr=split(//,$line);
#store chromosome length
$line=length $line;
splice(@arr,($pad_start-1),$pad_length,@{$pad});

open (FHW, ">" . $chr . "_PAR_PADDED.fasta");
print FHW ">chromosome:${assembly}:Y:1:$line:1 chromosome Y\n";
for(my $i=0;$i <= scalar @arr; $i+=70)
{
print FHW join ("",@arr[$i..($i+69)]),"\n";
}
close FHW;
