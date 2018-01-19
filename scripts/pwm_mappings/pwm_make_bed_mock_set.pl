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

given a file of genomic regions and a file of fasta header lines from which the length of the chromosomes can be extracted, generates a set of randomly placed regions on the same chromosomes as the originals.

Currently the new may overlap the old and one another.

Currently accepts only pseudo-bed type files with 1 based coords.

=head1 USAGE

ls -1 K562* | xargs -I {} echo make_bed_mock_set.pl~ -g~ /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.id_lines~ -i '~' {} '~' '>' {} _mock | tr -d ' ' | tr '~' ' '

ls -1 ES* | xargs -I {} echo make_bed_mock_set.pl~ -g~ /data/blastdb/Ensembl/funcgen/mus_musculus_male_NCBIM37_unmasked.id_lines~ -i '~' {} '~' '>' {} _mock | tr -d ' ' | tr '~' ' '


=cut



#TO DO
# 1 Strip all of this out and use SeqUtils::randomise_bed_file

use strict;
use DBI;
use Env;

use Getopt::Std;
use IO::Handle;
use IO::File;

$| = 1; #no output buffer

my($user, $password, $driver, $host, $port);
my $outfile='';
my $infile='';
my $id_list;
my $sp;
my $verbose = 0;
my $no_overlap = 0;
my $genome_file;

my %opt;

if ($ARGV[0]){
  &Getopt::Std::getopts('g:h:o:i:n', \%opt) || die ;
}else{
  &help_text; 
}

&process_arguments;

my $ofh;
if($outfile){
     open($ofh,">$outfile") or die "couldn't open file $outfile\n";
}else{
        #or write to STDOUT
     $ofh = new IO::File;
     $ofh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}

open(IN,$genome_file) or die "failed to open $genome_file";
my %chrom;
while(my $line = <IN>){

  #warn $line;
  chomp $line;
  
  #This is new fasta header format as of release 76
  #This needs moving to SeqTools or similar, so we always use the same code for fasta header handling!
  
  my($sr_name, undef, $slice_name) = split(/\s+/, $line);
  $sr_name =~ s/^>//;
  my(undef, undef, undef, $min, $max) = split(/:/,$slice_name);
  $chrom{$sr_name}->{max} = $max;
  $chrom{$sr_name}->{min} = $min;

  if($min != 1){
    die "$sr_name has high start in genome file:\t$genome_file\nScript needs updating to deal with this ";
  }
}
close(IN);


#print $chrom{11}->{max}."\n";

my @orig;
open(IN,$infile) or die "failed to open $infile";
while(my $line = <IN>){
    chop $line;
    #my @field = split("\t",$line);
   # my @feat = @field[0..2];
   # push @orig, \@feat;
  push @orig, [split("\t",$line)];
}
close(IN);

my $wrote = 0;

foreach my $aref (@orig){
  #print join("\t",@$aref)."\n";
  my $len = $aref->[2] - $aref->[1] +1;

  ##### HACK HACK DS Test
  #If the length of feature is greater than genome region where it is we can ignore this one... ?
  # This region is most likely atrefactual and should be eliminated...
  if($len >= $chrom{$aref->[0]}->{max}){
    print STDERR "Length of feature greater than lenght of genome region!\n";
    next;
  }
  ### END OF HACK

  my $new_end = 0;
  my $tries = 0;

  while($new_end < $len){ # should also consider chrom start coord here
    $new_end = int(rand($chrom{$aref->[0]}->{max}));
    $tries ++;
    
    if($tries > 1000){
      die "Data problem :\n".
        "max for ".$aref->[0]." = ".$chrom{$aref->[0]}->{max}."\n".
        join("\t",@$aref)." len = $len\n";
    }
  }

  my $new_start = $new_end - $len +1;
  #my $new_len =  $new_end -  $new_start +1;

  $wrote++;
  print $ofh $aref->[0]."\t$new_start\t$new_end\n";#.join("\t", $aref->[3..$#{$aref}])."\n";
}

close($ofh);

if(! $wrote){
  die("Failed to write any mock peaks to:\t$outfile");
}
else{
  warn "Wrote $wrote mock peaks to:\t$outfile\n";  
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



sub get_names_from_file{
    my $file=shift;

    open(IN, "< $file") or die "couldn't open list file $file";

    my @ret;
    while( <IN> ){
        chop;
	push  @ret, $_ ;
    }

    my $text_list = join(",",@ret);
    #print $text_list;
    return $text_list;

}


sub config{
 
($user =     $ENV{'ENSMARTUSER'}) or return(0); # ecs1dadmin
($password = $ENV{'ENSMARTPWD'}) or return(0); #
($host   =   $ENV{'ENSMARTHOST'}) or return(0); #localhost
($port =     $ENV{'ENSMARTPORT'}) or return(0); #3360
($driver  =  $ENV{'ENSMARTDRIVER'}) or return(0); #mysql
 
}
   
sub err{
    print STDERR "$_[0]\n";
}
  



sub print_mart_results{
    my $res_aref = shift;
    my $ma; #dummy to avoid undef err

    if( $ma->result_type() eq 'sequence'){
	foreach my $aref (@$res_aref){
	    print $aref->[0]."\n";
	    print $aref->[1]."\n";
	}
    }else{
	foreach my $aref (@$res_aref){
	    foreach my $field (@$aref){
	           ($field)? print $field." ":print "unavailable ";
	    }
	    print "\n";
	}
    }
}


sub make_contact{
    my $db = shift;

    unless($driver && $db && $host && $port && $user){
	&err("DB connection parameters not set");
        exit(1);
    }

    # hook up with the server
    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port;mysql_local_infile=1";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 0});
    if ($dbh){
        print STDERR ("connected to $db\n");
    }else{
        &err("failed to connect to database $db");
        exit(1);
 
    }   

            
    return $dbh;
}



sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }


    if (exists $opt{o}){
        $outfile = $opt{o};
    }  


    if (exists $opt{i}){
        $infile = $opt{i};
    }else{
	&help_text(" Please supply the name of a bed file");
    }

    if (exists $opt{n}){
        $no_overlap = 1;
    }



    if  (exists $opt{g}){
        $genome_file = $opt{g};
    }else{
	&help_text("You must provide a file containing fasta ID lines with long sequence identifiers");
    }

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    .pl [-h] for help

                  [-o] <output file> - name of a bed file for output
                   -i  <input file> - name of a bed file 
                  [-n] flag no overlap between old and new set
                   -g  <input file> - file of genome fasta id lines
                  [-] 
                  [-] <> 
                  [-] <> 


END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
