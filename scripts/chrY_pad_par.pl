#!/usr/bin/env/perl
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
