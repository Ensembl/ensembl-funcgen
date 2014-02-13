#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. You may obtain a copy of the License at

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

PWM to genome mapper for routine use by Funcgen.

Uses (a slightly modified) find_pssm_dna C++ program from the MOODS suite as the search engine.

Input comprises a fasta file of DNA sequences and a list of jaspar format pfm files and their associated matrix_list.txt file. The input file is split using the fastaexplode program from the Exonerate suite (http://www.ebi.ac.uk/~guy/exonerate/).

Output is a bed style, tab separated file but with 1-based rather than 0 based coordinates.

=head1 USAGE

/ensembl-funcgen/scripts/pwm_genome_map.pl -g seqs.fasta -o mappings.out 
 JASPAR_CORE_2008/*.pfm

Options 

 -g filename. The file contains 1 or more DNA sequences in fasta format. Typically these would constitute an entire genome, but could also be shorter eg promoter sequences. If the sequences are softmasked the mapping software will convert lower case to upper case internally before doing the mapping. Hence mappings will be provided for softmasked regions.

 -o output_file contains mappings for all PWMs to all sequences. Columns are tab separated. Coordinates are 1 based.

<sequence_name> <start> <end> <motif_name> <score> <strand>



 -a <assembly_version> This option acts as both a flag to abbreviate the chromosome and supercontig names and lets the script know what assembly name occurs in the full length chromosome names.
 
eg if the full chromosome names look like 

>supercontig::HSCHRUN_RANDOM_CTG33:1:41933:1 
>chromosome:GRCh37:Y:1:59373566:1 

then -a GRCh37 will cause them to be shortened to 

HSCHRUN_RANDOM_CTG33
Y

in the output.



The amount of memory required by find_pssm_dna depends on the number of mappings which will be generated. Initially this is unpredictable, so it is wise to allow plenty of memory. Trial and error suggests that 15G is sufficient for the human genome and several hundred PWMs at a probability threshold of 0.001.


Extra blank lines or trailing tabs in pfm files cause the error :-

"Matrix XXX.pfm has wrong alphabet size. Omitting"

=head1 EXAMPLES

/usr/bin/nice -n19 /nfs/users/nfs_d/dkeefe/src/head/ensembl-funcgen/scripts/pwm_genome_map.pl  -g /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.fa examples/data/matrix/JASPAR_CORE_2008/*.pfm

rm -f bsub* ; bsub -q long -o bsub_out -e bsub_err -R 'select[mem>15000] rusage[mem=15000]' -M 15000000 /nfs/users/nfs_d/dkeefe/src/head/ensembl-funcgen/scripts/pwm_genome_map.pl -g /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.fa examples/data/matrix/JASPAR_CORE_2008/*.pfm

pwm_genome_map.pl -g short_seqs.fa -a GRCh37 -o /lustre/scratch103/ensembl/dkeefe/jaspar_pwm_genome_map.out -t transfac examples/data/matrix/transfac32/matrix.dat 

bsub -q normal -o bsub_out_jasp -e bsub_err_jasp -R 'select[mem>15000] rusage[mem=15000]' -M 15000000 pwm_genome_map.pl -g /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.fasta -a GRCh37 -o JASPAR_CORE_v_GRCh37_all.tab -p 0.001 /data/blastdb/Ensembl/funcgen/pwm/JASPAR_CORE_Oct2009/vertebrates/*.pfm

bsub -q normal -o bsub_out_tfac -e bsub_err_tfac -R 'select[mem>15000] rusage[mem=15000]' -M 15000000 pwm_genome_map.pl -g /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.fa -a GRCh37 -o transfac32_v_GRCh37_all.tab -p 0.001  -t transfac /data/blastdb/Ensembl/funcgen/pwm/transfac32/matrix.dat


    bsub -q normal -o bsub_out_tfac -e bsub_err_tfac -R 'select[mem>15000] rusage[mem=15000]' -M 15000000 pwm_genome_map.pl -g /lustre/scratch103/ensembl/funcgen/bwa_indexes/MUS_MUSCULUS/mus_musculus_male_NCBIM37_unmasked.fasta -a NCBIM37 -o transfac32_v_NCBIM37_all.tab  -p 0.001  -t transfac /data/blastdb/Ensembl/funcgen/pwm/transfac32/matrix.dat




=head1 SEE ALSO
    
mysql -u ensro -hens-genomics2 -P3306 -BN -e"select name from feature_type where class in( 'Transcription Factor','Insulator')" dev_homo_sapiens_funcgen_58_37c

mysql -u ensro -hens-genomics1 -P3306 -BN -e"select name from feature_type where class in( 'Transcription Factor','Insulator') " dev_mus_musculus_funcgen_59_37l



# to get peak coords

mysql -u ensro -hens-genomics2 -P3306 -BN -e"select sr.name,af.seq_region_start,af.seq_region_end,af.score,ft.name,ct.name from annotated_feature af,feature_set fs,feature_type ft,seq_region sr,cell_type ct where ft.class in('Transcription Factor','Insulator') and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id and sr.seq_region_id = af.seq_region_id and sr.schema_build = '58_37c' and fs.cell_type_id = ct.cell_type_id and ft.name = 'Srf'" dev_homo_sapiens_funcgen_58_37c 


mysql -u ensro -hens-genomics1 -P3306 -BN -e"select sr.name,af.seq_region_start,af.seq_region_end,af.score,ft.name,ct.name from annotated_feature af,feature_set fs,feature_type ft,seq_region sr,cell_type ct where ft.class in('Transcription Factor','Insulator') and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id and sr.seq_region_id = af.seq_region_id and sr.schema_build = '58_37k' and fs.cell_type_id = ct.cell_type_id and ft.name = 'Cmyc'" dev_mus_musculus_funcgen_59_37l

=head1 TO DO

add perl implementation of fastaexplode

remove unused code

maybe collate output into one file per PWM

=cut

# To do - Some o fthese apply to all the pwm scripts
# 1 Getopt::Long, remove process_arguments
# 2 remove backticks in favour of system. Use EFGUtils::run_system_cmd
# 3 Remove commentary and err subs and flip use STDERR(warn) and STDOUT(print). 
#   We should have nothing on STDERR unless there is a warning/error.
# 4 DONE fastaclean path? Now append path belowb
# 5 Fix backtick which tests do not work
# 6 remove config sub

#Requirements
#Do this in a pre-exec before we process anything?
#test dirs too?
#fastaclean
#$moodsmapper i.e. find_pssm_dna


use strict;
use Env;
use Getopt::Std;
use File::Basename;
use IO::Handle;
use IO::File;
use DBI qw( :sql_types );
#use lib '/nfs/users/nfs_d/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
$ENV{PATH} .= ':/usr/local/ensembl/bin/'; #fastaclean home, probably needs moving to /software/ensembl/bin

my($user, $password, $driver, $host, $port);
my $outfile='';
my $infile='';
my $sp;
my $verbose = 0;
my $pwm_type = 'jaspar';
my $work_dir = "/lustre/scratch103/ensembl/dkeefe/pwm_genome_map_$$/";
#my $work_dir = "/lustre/scratch101/ensembl/ds19/tmp/";

my $genome_file;
#Changed DS
my $moods_mapper = '/software/ensembl/funcgen/find_pssm_dna';
#my $moods_mapper = '~dkeefe/bin/find_pssm_dna';
my $thresh = 0.001;
my $assembly;

my %opt;

if(! @ARGV){
  &help_text;
}

print "$0 @ARGV\n";
&Getopt::Std::getopts('a:s:t:g:h:v:o:i:p:w:', \%opt) || die ;


# get configuration from environment variables
#&config; # this may fail but config can be on command line


&process_arguments;
my @pwm_files = @ARGV;

`rm -rf $work_dir`;
#$verbose = 2;
&backtick("mkdir -p $work_dir"."/genome");


# irrespective of pwm_type we need to create the rev-comp matrices and put all the matrices in a working directory - along with a composite matrix_list.txt file
my %file_max_score;
my $matrix_file;
if($pwm_type eq 'jaspar'){

    # we assume the matrix_list.txt file is in same dir as PWMs
    #$matrix_file = &find_matrix_txt(\@pwm_files) or die
    #    "Unable to find a matrix_list.txt file for the jaspar matrices \n".
    #    join("\n",@pwm_files[0..5])." etc.\n";
    #&backtick("cp $matrix_file $work_dir");

    # we move the PWM files to the working dir and create rev-comp files
    # there too, making additions to the matrix_list in working dir as we go
    foreach my $file (@pwm_files){
        &backtick("cp $file $work_dir"); 
        &rev_comp_matrix($file,$work_dir);
	my($min,$max)= &matrix_min_max($file,$work_dir);
	print basename($file)."\t$min\t$max\n" if $verbose > 1;
	$file_max_score{ basename($file) } = $max;
    }

}
elsif($pwm_type eq 'transfac'){
    # if the pwm_type isn't Jaspar we need to convert them into Jaspar pfm 
    # format and create a matrix_list.txt file. Put new PWM files etc. in 
    # working directory
    # Transfac matrixes all come in a single file 
    if(@ARGV > 1){
	die "Transfac matrices should all be in a single file, not ".
             join("\n",@ARGV)."\n";
    }

    my @pwm_files = &parse_matrixdat_2_pfm($ARGV[0],$work_dir);

    foreach my $file (@pwm_files){
	&rev_comp_matrix($file,$work_dir);
	my($min,$max)= &matrix_min_max($file,$work_dir);
	print basename($file)."\t$min\t$max\n";
	$file_max_score{ basename($file) } = $max;
    }
}
else{
    die "Can only handle Jaspar and Transfac formats at present\n";
}



# we need to explode the genome.fasta file into individual sequences.
# and if abbreviated chr names have been requested change the file names
# as these are used in the output file
# also if the sequences contain characters other than [acgtnACGTN] they need
# to be converted to N otherwise find_pssm_dna outputs the wrong coords.
#$verbose = 2;
my @chr_files = &explode_genome_fasta($genome_file,$work_dir.'genome',$assembly);

my $jdbh;
my $dsn = sprintf( "DBI:%s:%s:host=%s;port=%s",
                     'mysql', 'JASPAR_v5_0', 'ens-genomics1',
                     3306 );

eval {
  $jdbh = DBI->connect( $dsn, 'ensro', undef,
                          { 'RaiseError' => 1 } );
  };

my $error = $@;

  if ( !$jdbh || $error || !$jdbh->ping ) {
    die( "Could not connect to database using [$dsn] as a locator:\n" . $error );
  }
  

# now we map all the PWMs to each genomic sequence

#NJ parallelise this on the farm!!!
#This is currently running to triple the memory footprint of previous runs!
#what is being cached in this loop, and why is there so much more of it for a modest increase 
#in the size of Jaspar
#%file_max_score might be the cuplrit
#Nope this is only accessed in parse_out_2_tab, not populated

foreach my $chr_file (@chr_files){

    my $tab=my $out=$chr_file;
    $tab =~ s/fa/tab/;
    $out =~ s/fa/out/;
    my $command = "$moods_mapper -f  $thresh $chr_file $work_dir"."*.pfm > $out";
    warn $command."\n";
    &backtick("$command");
    &parse_out_2_tab($out,$tab,$jdbh,\%file_max_score);
    #&backtick("rm -f $out");
    #&backtick("rm -f $chr_file");
}


# collate individual .tab files into the specified output file
my $command = "cat $work_dir/genome/*.tab > $outfile ";
&backtick($command);

# count how many mappings we got for each matrix
#$command = "cut -f4 $outfile | sort| uniq -c > $work_dir"."pwm_mapping_counts_".$pwm_type;
#&backtick($command);

# tidy up
`rm -rf $work_dir`;
exit;


 
###################################################################
sub matrix_min_max{
    my($file,$work_dir)=@_;

    open(IN,$file) or die "failed to open $file";

    my %swap =( 0 => 3,
                1 => 2,
                2 => 1,
                3 => 0
	       );

    my @mat;
    my $rows = 0;
    my $cols;
    my $sum; # we assume each col has the same sum, so just calc first
    while(my $line = <IN>){
	chop $line;
	unless($line =~ /[0-9]/){next}
        $line =~ s/^\s+(.*)/$1/; # remove leading whitespace
        $line =~ s/[\[\]]//g; #remove brackets if any
	my @field = split(/\s+/,$line); # split on white space
	#print join("\t",@field)."\n";
	#print join("~",@field)."\n";
	#my @rev= reverse(@field); # to get rev comp
	#$mat[$swap{$rows}] = \@rev;# to get rev comp
	#print join("\t",reverse(@field))."\n";
	$mat[$rows] = \@field;
        $rows++;
	$sum += $field[0];
    }
    if($rows > 4){ die "too many rows in file $file" }
    close(IN); 

    $cols = scalar(@{$mat[0]});
    # convert to probabilities and transpose to facilitate get_max_add()
    my @new_mat;
    for(my $i = 0;$i<$rows;$i++){
	for(my $j=0;$j<$cols;$j++){
	    $new_mat[$j][$i] = ($mat[$i]->[$j])/$sum;
	}
    }


    # jaspar IDs can be MA for core, CN for CNE and PB and PF for PHYLOFACTS 
    # our own PWM IDs are FG
    my($matrix_id) = $file =~ /.*([MPCF][AFNGBHL][0-9]+).[1-9].pfm/;
    print $file."\n".$matrix_id."\n";
    unless($matrix_id){die "problem parsing matrix name $file"}

    my($min_pos,@sums) = &get_max_add(@new_mat);
    return($min_pos,$sums[0]);

}



# this doesn't do exactly the same as find_pssm_dna cos the latter does some
# approximating to allow the use of integer arithmetic.
sub get_max_add{
    my @mat = @_;

    my $len = scalar(@mat);
    #print $len."\n";

    my @maxes;
    my $min_sum;
    for(my $i=0;$i <= $#mat;$i++){
	my @list = sort {$b <=> $a} @{$mat[$i]}; # max is element 0
	#$maxes[$i] = ((($list[0])+0.01)/1.04)/0.25;
	#$maxes[$i] = &log2($maxes[$i]);
	#$min_sum += &log2( ((($list[3])+0.01)/1.04)/0.25 );

	#$maxes[$i] = $list[0]/0.25;
	#$maxes[$i] = &log2($maxes[$i]);
	#$min_sum += &log2( $list[3]/0.25 );

	$maxes[$i] = ((($list[0])+0.002)/1.008)/0.25;
	$maxes[$i] = &log2($maxes[$i]);
	$min_sum += &log2( ((($list[3])+0.002)/1.008)/0.25 );

    }
 
    # last element contains bit score for its own most frequent letter
    # second last contains bit score for its own most frequent letter +last 
    # etc.   
    my @sums;
    $sums[$len] = 0;
    for(my $i=$#mat;$i >= 0;$i--){
        $sums[$i] = $maxes[$i]+$sums[$i+1];
    }

    #print join(" ",@maxes)."\n";
    #print join(" ",@sums)."\n";

    return($min_sum,@sums);
}

sub log2{
    my $n = shift;

    return log($n);#/log(2.0);

}


# converts a single transfac matrix.dat file to multiple jaspar pfm files
sub parse_matrixdat_2_pfm{
    my($tf_file,$work_dir)=@_;

    my $ifh;
    open($ifh,"$tf_file") or die "failed to open file $tf_file";

    my $list_file = $work_dir.'matrix_list.txt';
    my $ofh;
    open($ofh,"> $list_file") or die "failed to open file $list_file";

    my $ac;
    my $name;
    my @files;
    while(my $line = <$ifh>){
	chop $line;

	if($line =~ /^AC/){($ac) = $line =~ /AC  (M[0-9]+)/; }
        if($line =~ /^NA/){($name) = $line =~ /NA  (.+)/; }
	if($line =~ /^P0/){
	my $pfm = &process_matrix($ifh,$ofh,$ac,$name,$work_dir);
	push @files, $pfm;
	}

    }

	close($ifh);
	close($ofh);
    return @files;
}

# reads in a transfac matrix, transposes it and outputs it into a jaspar pfm
# file 
sub process_matrix{
    my($ifh,$ofh,$ac,$name,$work_dir)=@_;

    $ac =~ s/M0/MA/;
    $ac .= '.1'; # all files are give version 1
    # print to list_file
    print $ofh $ac."\t0.0\t ".$name."\n";

    open(OUT,">".$work_dir.$ac.'.pfm') 
        or die "failed to open ".$work_dir.$ac.'.pfm';

    my @mat;
    while(my $line = <$ifh>){
	chop $line;
	if($line =~ /^XX/){
	    my $cols = scalar(@mat);
	    for(my $i=0;$i<4;$i++){
		for(my $j=0;$j<$cols;$j++){
		    print OUT $mat[$j][$i];
		    unless($j == $cols-1){print OUT ' '}
		}
		print OUT "\n";
	    }

	    close(OUT);
	    return($work_dir.$ac.'.pfm');
	}else{
	    my @field = split(/\s+/,$line);
	    my @arr = @field[1..4];
	    push @mat, \@arr;
	}

    }

}




sub parse_out_2_tab{
   #my($out,$tab,$matrix_file,$max_scores_ref)=@_;

  my($out,$tab,$dbh,$max_scores_ref)=@_;
  
  open(IN,$out) or die "couldn't open file $out ";
  open(OUT,">$tab") or die "couldn't open file $tab ";
  
  my $chr = basename $out;
  $chr =~ s/\.out//;
  
  my $strand;
  my $id;
  my $addn;
  my $tf_name;
  my $score_thresh;
  while(my $line = <IN>){
    chop $line;
    if($line eq ''){
	    next;
    }
    elsif($line =~ /^>/){
	    $line =~ s/>//;
	    $id = basename $line;
	    $id =~ s/\.pfm//;
      
      # get TF name from matrix_list file
	    #$tf_name = &backtick("grep '^$id	' $matrix_file |cut -f3");
	    
	    ($tf_name) = $dbh->selectrow_array("SELECT NAME from MATRIX where BASE_ID='$id'");
	    
	    
	    
	    chop($tf_name); #wtf? taking off the version but not the dot?
	    unless($tf_name){
	      warn "Failed to get TF name from the DB:\t$id\n";
             #   warn "Failed to get TF name for $id from $matrix_file";
	    }
      $strand = '1';
	    if($id =~ 'rc'){
        $strand = '-1';
        $id =~ s/rc//;
	    }
	    $score_thresh = $max_scores_ref->{$id.'.pfm'} * 0.7; # 70% of max
	    $score_thresh = 0;
	    &commentary( $max_scores_ref->{$id.'.pfm'} ." = max score for $id $tf_name\n") if $verbose;
      
	    $line = <IN>;
	    ($addn) = $line =~ /.*Length:([0-9]+)/;
    }
    else{
	    my($start,$score)= split("\t",$line);
	    if($score > $score_thresh){
        print OUT join("\t",$chr,($start+1),($start+$addn),$tf_name,
                       $score,$strand,$id)."\n" or 
                         die "failed to print to $tab" ;
	    }
    }
  }
  close(IN);
  close(OUT);
}


# get the individual sequences from the genome file and put them in files
# which have the fasta id as their name and an extension of .fa
# optionally reduce chromosome name to chr_name as in ensembl databases
# returns the list of files produced or dies
sub explode_genome_fasta{
  my($genome_file,$target_dir,$assembly) = @_;
  &commentary("exploding $genome_file to $target_dir") if $verbose;
  
  if (! $genome_file){
    die('Genome file not defined');
  }
  

  warn "GENOME file is $genome_file";
  
	#These tests do not work!
  
  my $res = &backtick("which fastaclean");
  if($res =~ 'not found'){
    die "Need to implement fastaclean functionality here";
  }
  
  
  $res = &backtick("which fastaexplode");
    if($res =~ 'not found'){
      die "Need to implement fastaexplode functionality here";
    }else{
      &backtick("fastaexplode -f $genome_file -d $target_dir");
    }
  
  
  $res = &backtick("ls -1 $target_dir"."/*.fa");
  my @lines = split("\n",$res);
  
  unless(@lines > 0){
    die "ERROR: No genome files produced by splitting".$genome_file;
  }
  
  # now we alter the file names if short chromosome names have been requested
  if($assembly){

    my @short_names;
    foreach my $file_path (@lines){
	    my $name = basename $file_path;
      # sed 's/chromosome:$assembly:/chr/'|sed 's/supercontig:$assembly://' |sed 's/supercontig:://' | sed 's/:[0-9:]*//'
      
	    $name =~ s/chromosome:$assembly://;
	    $name =~ s/supercontig:$assembly://;
	    $name =~ s/supercontig:://;
	    $name =~ s/:[0-9:]*//;
	    my $path = dirname $file_path;
	    my $new = $path.'/'.$name;

      #This is not working as all but the Y chr are already short names
      #temp hack to get it running
      my $command;

      if($file_path =~ /chromosome:GRCh37:Y/){

      

        $command = "fastaclean -a -f $file_path > $new ; rm -f $file_path";
      }
      else{
        my $path = dirname $file_path;
        my $new = $path.'/tmp';
        $command = "fastaclean -a -f $file_path > $new ;".
        " rm -f $file_path ;mv $new $file_path";
      }


        warn "executing: $command";



	    &backtick($command);
	    push @short_names, $new;
    }
    return @short_names;
  }else{
    foreach my $file_path (@lines){
      # just do fastaclean
	    my $path = dirname $file_path;
	    my $new = $path.'/tmp';
	    my $command = "fastaclean -a -f $file_path > $new ;".
        " rm -f $file_path ;mv $new $file_path";

      warn "executing: $command";
	    &backtick($command);
    }
    
  }
  
  
  return @lines;
}



# assumes a jaspar matrix pfm file with 4 rows A C G T
sub rev_comp_matrix{
    my($file,$work_dir)=@_;

    open(IN,$file) or die "failed to open $file";

    my %swap =( 0 => 3,
                1 => 2,
                2 => 1,
                3 => 0
	       );

    my @mat;
    my $rows = 0;
    my $cols;
    while(my $line = <IN>){
	chop $line;
	unless($line =~ /[0-9]/){next}
        $line =~ s/^\s+(.*)/$1/; # remove leading whitespace
        $line =~ s/[\[\]]//g; #remove brackets if any
	my @field = split(/\s+/,$line); # split on white space
	#print join("\t",@field)."\n";
	#print join("~",@field)."\n";
	my @rev= reverse(@field);
	$mat[$swap{$rows}] = \@rev;
	#print join("\t",reverse(@field))."\n";
        $rows++;

    }
    if($rows > 4){ die "too many rows in file $file" }
    close(IN); 

    # jaspar IDs can be MA for core, CN for CNE, PB, PH and PF for PHYLOFACTS
    # our own PWM IDs are FG
    #my($matrix_id) = $file =~ /.*([MPCF][AFNG][0-9]+).pfm/;
    my($matrix_id,$version) = $file =~ /.*([MPCF][AFNGBHL][0-9]+).([0-9]*).pfm/;
    print $file."\n".$matrix_id." $version\n";
    unless($matrix_id){die "problem parsing matrix name $file"}
 

    my $rc_file = $work_dir.$matrix_id.'rc.'.$version.'.pfm';
    print $rc_file."\n";
    open(OUT,"> $rc_file") or die "failed to open file $rc_file";

    for(my $i=0;$i<$rows;$i++){
        print OUT join(" ",@{$mat[$i]})."\n";
    }

    close(OUT);

    #NJ don't need this now, we just need to add support fo DB access where ever this is used
    #was this ever being used again

    # get the relevant line from the matrix_list.txt file
    # by grepping for the id at the start of the line
    # add rc to the ID and append the line to matrix_list.txt
    #my $res = &backtick("grep '^$matrix_id.$version' $work_dir".
 #                       "matrix_list.txt");
 #   print $res;
    #chop($res);
    #my @field = split("\t",$res);
    ##$field[0] .= 'rc';
    #$field[0] = $matrix_id.'rc.'.$version;
    #$res = join("\t",@field);
    ##print $res."\n";
    #open(OUT, ">> $work_dir"."matrix_list.txt") or 
    #    die "failed to open $work_dir"."matrix_list.txt for appending";
    #print OUT $res."\n" or die "failed to write to  $work_dir".
    #                           "matrix_list.txt";
    #close(OUT);


}



# There should be a matrix_list.txt file in the same directory as the pfm files
sub find_matrix_txt{
    my($pwm_files_aref) = @_;

    foreach my $file_path (@$pwm_files_aref){
	my @field = split('/',$file_path);

	if(@field ==1){ 
            # we're working in the matrix directory
            if( -e "./matrix_list.txt"){
		return("./matrix_list.txt");
	    }
	}else{
	    pop @field;
	    my $check = join('/',@field).'/matrix_list.txt';
	    #print "checking $check\n";
	    if( -e $check){return($check)}
	}

    }

    return '';

}


# makes shell execute command and checks for errors
# returns the output of the command unprocessed
sub backtick{
    my $command = shift;

    warn "executing $command \n" if ($verbose ==2);

    my $res = `$command`;
    if($?){
        warn "failed to execute $command\n";
        warn "output:\t$res\n";
	die "exit code $?";
    }

    return $res;
}




# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print STDERR "$_[0]";
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
  




sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }




    if (exists $opt{w}){
        $work_dir = $opt{w}; 
    }

    if (exists $opt{p}){
        $thresh = $opt{p}; 
    }

    if (exists $opt{P}){
        $port = $opt{P}; 
    }

    if (exists $opt{u}){
        $user = $opt{u}; 
    }

    $sp = 'homo_sapiens';
    if (exists $opt{s}){
        $sp = lc($opt{s});
    }


    if (exists $opt{o}){
        $outfile = $opt{o};
    }else{
	&help_text("Please give a name for the output file");
    }


    if (exists $opt{g}){
        $genome_file = $opt{g};
    } 
 
    if (exists $opt{i}){
        $infile = $opt{i};
    }

    if (exists $opt{t}){
        $pwm_type = lc( $opt{t} );
    }

    if (exists $opt{v}){
	$verbose = $opt{v};
    }

    if (exists $opt{a}){
	$assembly = $opt{a};
    }


} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

  pwm_genome_map.pl -g fasta_file -o output_file [options] pwm_file1 [pwm_file2 pwm_file3 ...]

                  [-h] for help
                  [-a] <assembly_version> eg GCRh37 as used in full chr. name
                  [-t] <pwm_type> 'jaspar'=default, 'transfac' 
                  [-p] <probability> threshold for moods mapper default=0.001
                   -g  <file_name> genome fasta file
                   -o  <output file> - name of a file for output
                  [-v] <integer> verbosity level
                  [-s] <species> eg -smus_musculus, default = homo_sapiens  
                  [-w] <dir_name> path of a (spacious) working directory  
                  [-] 
                  [-] <> 
                  [-] <> 


END_OF_TEXT
#                  [-i] <input file> - list of TFs with thresholds

    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}

1;
