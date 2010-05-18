#!/usr/local/ensembl/bin/perl

=head1 DESCRIPTION

PWM to genome mapper for routine use by Funcgen.

Uses (a slightly modified) find_pssm_dna C++ program from the MOODS suite as the search engine.

Input comprises a fasta file of DNA sequences and a list of jaspar format pfm files and their associated matrix_list.txt file.

Output is a bed style, tab separated file but with 1-based rather than 0 based coordinates.


Just under 1h for chr1 with thresh 0.005
561,015,360 mappings in total

5340884 mappings for MA0003
 442082 using my own method

1561940 with thresh 0.001

=head1 AUTHOR(S)

dkeefe@ebi.ac.uk

=head1 USAGE

/ensembl-functgenomics/scripts/pwm_genome_map.pl -g seqs.fasta -o mappings.out 
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



The amount of memory required by find_pssm_dna depends on the number of mappings which will be generated. Initially this is unpredictable, so it is wise to allow plenty of memory. Trial and error suggests that 15G is sufficient for chr1 and several hundred PWMs at a probability threshold of 0.001.


Extra blank lines in pfm files cause the error :-

"Matrix XXX.pfm has wrong alphabet size. Omitting"

=head1 EXAMPLES

/usr/bin/nice -n19 /nfs/users/nfs_d/dkeefe/src/head/ensembl-functgenomics/scripts/pwm_genome_map.pl  -g /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.fa examples/data/matrix/JASPAR_CORE_2008/*.pfm

rm -f bsub* ; bsub -q long -o bsub_out -e bsub_err -R 'select[mem>15000] rusage[mem=15000]' -M 15000000 /nfs/users/nfs_d/dkeefe/src/head/ensembl-functgenomics/scripts/pwm_genome_map.pl -g /data/blastdb/Ensembl/funcgen/human_male_GRCh37_unmasked.fa examples/data/matrix/JASPAR_CORE_2008/*.pfm

=head1 SEE ALSO


=head1 CVS

 $Log: not supported by cvs2svn $

=head1 TO DO

deal with TRANSFAC data

add perl implementation of fastaexplode

add name from matrix_list.txt to bed file

remove unused code

=cut


use strict;
#use DBI;
use Env;
use Getopt::Std;
use File::Basename;
use IO::Handle;
use IO::File;
#use lib '/nfs/users/nfs_d/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';

my($user, $password, $driver, $host, $port);
my $outfile='';
my $infile='';
my $sp;
my $verbose = 0;
my $pwm_type = 'jaspar';
my $work_dir = "/lustre/scratch103/ensembl/dkeefe/pwm_genome_map_$$/";
my $genome_file;
my $moods_mapper = '~dkeefe/bin/find_pssm_dna';
my $thresh = 0.001;
my $assembly;

my %opt;

if ($ARGV[0]){
    &Getopt::Std::getopts('a:s:t:g:h:v:o:i:', \%opt) || die ;
}else{
    &help_text; 
}


# get configuration from environment variables
#&config; # this may fail but config can be on command line


&process_arguments;
my @pwm_files = @ARGV;

`rm -rf $work_dir`;
#$verbose = 2;
&backtick("mkdir -p $work_dir"."/genome");


# irrespective of pwm_type we need to create the rev-comp matrices and put all the matrices in a working directory - along with a composite matrix_list.txt file
my $matrix_file;
if($pwm_type eq 'jaspar'){

    # we assume the matrix_list.txt file is in same dir as PWMs
    $matrix_file = &find_matrix_txt(\@pwm_files) or die
        "Unable to find a matrix_list.txt file for the jaspar matrices \n".
        join("\n",@pwm_files[0..5])." etc.\n";
    &backtick("cp $matrix_file $work_dir");

    # we move the PWM files to the working dir and create rev-comp files
    # there too, making additions to the matrix_file in working dir as we go
    foreach my $file (@pwm_files){
        &backtick("cp $file $work_dir"); 
	&rev_comp_matrix($file,$work_dir,$matrix_file);
    }

}
elsif($pwm_type eq 'transfac'){
    # if the pwm_type isn't Jaspar we need to convert them into Jaspar pfm 
    # format and create a matrix_list.txt file. Put new PWM files etc. in 
    # working directory
    # Transfac matrixes all come in a single file 




    die "Can only handle Jaspar format at present\n";

}
else{
    die "Can only handle Jaspar format at present\n";
}


# we need to explode the genome.fasta file into individual sequences.
# and if abbreviated chr names have been requested change the file names
# as these are used in the output file
$verbose = 2;
my @chr_files = &explode_genome_fasta($genome_file,$work_dir.'genome',$assembly);



# now we map all the PWMs to each genomic sequence
foreach my $chr_file (@chr_files){

    my $tab=my $out=$chr_file;
    $tab =~ s/fa/tab/;
    $out =~ s/fa/out/;
    my $command = "$moods_mapper $thresh $chr_file $work_dir"."*.pfm > $out";
    #print $command."\n";
    &backtick("$command");
    &parse_out_2_tab($out,$tab,$work_dir."matrix_list.txt");
    &backtick("rm -f $out");
    &backtick("rm -f $chr_file");
}


# collate individual .tab files into the specified output file
# optionally reduce chromosome name to chr_name as in ensemble db

my $command = "cat $work_dir/genome/*.tab > $outfile ";
&backtick($command);



#`rm -rf $work_dir`;
exit;


 
###################################################################
sub parse_out_2_tab{
    my($out,$tab,$matrix_file)=@_;

    open(IN,$out) or die "couldn't open file $out ";
    open(OUT,">$tab") or die "couldn't open file $tab ";

    my $chr = basename $out;
    $chr =~ s/\.out//;

    my $strand;
    my $id;
    my $addn;
    my $tf_name;
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
	    $tf_name = &backtick("grep '^$id	' $matrix_file |cut -f3");
	    chop($tf_name);
	    unless($tf_name){
                warn "Failed to get TF name for $id from $matrix_file";
	    }
            $strand = '1';
	    if($id =~ 'rc'){
		$strand = '-1';
		$id =~ s/rc//;
	    }
	    $line = <IN>;
	    ($addn) = $line =~ /.*Length:([0-9]+)/;
	}
#        elsif($line =~ /^H/){
#            # find_pssm_dna outputs 0 based coord
#	    ($addn) = $line =~ /.*Length:([0-9]+)/;
#	    
#	}
        else{
	    my($start,$score)= split("\t",$line);
	    print OUT join("\t",$chr,($start+1),($start+$addn),$tf_name,
                                $score,$strand)."\n";
	}
    }
    close(IN);
    close(OUT);
}


# get the individual sequences from the genome file and put them in files
# which have the fasta id as their name and an extension of .fa
# returns the list of files produced or dies
sub explode_genome_fasta{
    my($genome_file,$target_dir,$assembly) = @_;
    &commentary("exploding $genome_file to $target_dir") if $verbose;

    my $res = &backtick("which fastaexplode");
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
	    my $command = "mv $file_path $new";
	    &backtick($command);
	    push @short_names, $new;
	}
	return @short_names;
    }


    return @lines;
}



# assumes a jaspar matrix pfm file with 4 rows A C G T
sub rev_comp_matrix{
    my($file,$work_dir,$matrix_file)=@_;

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

    # jaspar IDs can be MA for core, CN for CNE and PF for PHYLOFACTS
    # our own PWM IDs are FG
    my($matrix_id) = $file =~ /.*([MPCF][AFNG][0-9]+).pfm/;
    #print $file."\n".$matrix_id."\n";

 

    my $rc_file = $work_dir.$matrix_id.'rc.pfm';
#    print $rc_file."\n";
    open(OUT,"> $rc_file") or die "failed to open file $rc_file";

    for(my $i=0;$i<$rows;$i++){
        print OUT join(" ",@{$mat[$i]})."\n";
    }

    close(OUT);

    # get the relevant line from the matrix_list.txt file
    # by grepping for the id followed by a tab
    # add rc to the ID and append the line to matrix_list.txt
    my $res = &backtick("grep '^$matrix_id	' $work_dir".
                        "matrix_list.txt");
 #   print $res;
    chop($res);
    my @field = split("\t",$res);
    $field[0] .= 'rc';
    $res = join("\t",@field);
#    print $res."\n";
    &backtick("echo '$res' >> $work_dir"."matrix_list.txt");


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
        $password = $opt{p}; 
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
                   -g  <file_name> genome fasta file
                   -o  <output file> - name of a file for output
                  [-i] <input file> - name of a file for input
                  [-n] <gene_name> ensembl stable gene id eg ENSG00000079974
                  [-v] <integer> verbosity level
                  [-s] <species> eg -smus_musculus, default = homo_sapiens  
                  [-w] <dir_name> path of a (spacious) working directory  
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
