#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

 bisulphite_qc.pl
  
=head1 SYNOPSIS

 bisulphite_qc.pl /path/to/bed/dnamethylationfeature_bed_file.bed
 
 Quality control check for DNAMethylationFeature bed files 

=head1 OPTIONS

perl bisulphite_qc.pl /path/to/dnamethylationfeature_bed_files/GM12878/GM12878_5MC_ENCODE_HAIB.bed

The only required argument is a DNAMethylationFeature bed file in Ensembl specific DNAMethylation bed format.
The output will be generated in the same directory in which the programme is run.

=head1 DESCRIPTION

This program computes fraction of methy-cytosines at each individual read depth between 1 and 100 and clubs everything
beyond 100 into one group.
The results are plotted into separate line charts for each cytosine context in the bed file using R. 


=head1 output

Text file for each of the cytosine contexts
R code file to generate the plots
Pdf of plots

=cut

use strict;
use warnings;

#INPUT BED FILE
my $file = $ARGV[0];
my ( %hash, %defined_contexts );
open( FH, $file ) or die "cant open $!";
my $line;
while ( chomp( $line = <FH> ) ) {
    my ( @temp, $depth, $context );
    @temp = split( /\s+/, $line );
    ( $context, $depth ) = split "/", $temp[3];
    $defined_contexts{$context} = undef;
    $depth = "gt_100" if $depth > 100;
    
# Here I'm checking methylation level for the cytosine and if it is greater than 0, the cytosine is considered to be methylated.
# However instead of using a 0 value, a better estimate for considering unmethylation would be very useful. This could be probably 
# done through computing a non-conversion/error rate estimation by aligning raw reads against mitochondrial genome which is not
# supposed to contain any methylation and the fraction of mitochondrial cytosines found methylated can be used here for unmethylation.  

    $temp[4] > 0
      ? $hash{ $context . "_" . $depth }{meth}++
      : $hash{ $context . "_" . $depth }{unmeth}++;

}

close(FH);

$file = ( split "/", $file )[-1];
$file =~ s/(.*).bed$/$1/;

open( FO, ">${file}.R" );
print FO "pdf(file=\"${file}.pdf\")\n";

for my $context ( sort keys %defined_contexts )

{
    my $outfile = "${file}_$context";
    open( FHW, ">" . $outfile );

    for my $depth ( 1 .. 100 )

    {
        if (    ( defined $hash{ $context . "_" . $depth }{meth} )
            and ( defined $hash{ $context . "_" . $depth }{unmeth} ) )
        {
            print FHW "$depth\t",
              $hash{ $context . "_" . $depth }{meth} /
              ( $hash{ $context . "_" . $depth }{meth} +
                  $hash{ $context . "_" . $depth }{unmeth} ), "\n";
        }

    }
    print FHW "101\t",
      $hash{ $context . "_gt_100" }{meth} /
      ( $hash{ $context . "_gt_100" }{meth} +
          $hash{ $context . "_gt_100" }{unmeth} );

    close(FHW);

    print FO "data <- read.table(\"$outfile\",sep=\"\t\",header=F)\n";
    print FO
"plot(data,type=\"l\",col=\"blue\",axes=F,xlab=\"Read depth\", ylab=\"methyl-Cs/Covered-Cs\",main=\"$context\")\n";
    print FO "axis(2)\n";
    print FO
"axis(1,las=2,cex.axis=0.8,at=c(1,seq(5,100,by=5)),lab=c(1,seq(5,100,by=5)))\n";
    print FO "box()\n";

}

print FO "dev.off()\n";

close(FO);
system "R CMD BATCH --slave ${file}.R";
