=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 DESCRIPTION

Determines the log odds score threshold for filtering PWM mappings to be incorporated into the functional genomics databases.

=head1 USAGE

see

ensembl-functgenomics/scripts/pwm_mapping_notes

grep '>'  /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.fasta > /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines

=head1 EXAMPLES

pwm_filter_mappings.pl -i list -e dev_homo_sapiens_funcgen_58_37c -H ens-genomics1 -g /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines

pwm_filter_mappings.pl -i list -e dev_homo_sapiens_funcgen_60_37e -H ens-genomics1 -P 3306 -u ensadmin -p big_secret -o thresholds 

=head1 SEE ALSO

ensembl-functgenomics/scripts/pwm_filter_mappings.pl
ensembl-functgenomics/scripts/pwm_genome_map.pl
ensembl-functgenomics/scripts/pwm_make_bed_mock_set.pl

=head1 TO DO

add the functionality from pwm_filter_mappings.pl to the end of this script

=head1 CVS

 $Log: not supported by cvs2svn $
 Revision 1.3  2010-11-26 15:56:07  dkeefe
 added a single linkage clustering of the peaks step

 Revision 1.2  2010-11-24 10:41:13  dkeefe
 Generally tidied up, added more POD and converted hard codded
 variables to command line options.

 Revision 1.1  2010-11-23 14:16:43  dkeefe
 Determines the log odds score threshold for filtering PWM mappings to be incorporated into the functional genomics databases.



=cut


use strict;
use DBI;
use Env;
use Getopt::Std;
use IO::Handle;
use IO::File;

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port);
my $outfile='';
my $infile='';
my $verbose = 2;
my $perc_thresh = 5; # FDR
my $genome_descriptor_file;# = "/data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines";
my %opt;

$| = 1; #no output buffer

if ($ARGV[0]){
&Getopt::Std::getopts('u:p:e:H:h:o:i:P:g:', \%opt) || &help_text("Invalid Argument") ;
}else{
&help_text; 
}


# get configuration from environment variables
&config; # this may fail but config can be on command line

my $enc_db = 'dev_homo_sapiens_funcgen_60_37e';# default, can be overridden by args
&process_arguments;

my $ofh;
if($outfile){
     open($ofh,">$outfile") or die "couldn't open file $outfile\n";
}else{
        #or write to STDOUT
     $ofh = new IO::File;
     $ofh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}



# hook up with the server

my $dbh = &make_contact($enc_db);

my %matrix_tf_pairs= &get_list_from_file($infile);

#%matrix_tf_pairs=(
# 'MA0139.1' => 'CTCF',
# 'MA0099.2' => 'Cjun',
# 'MA0099.2' => 'Cfos'

	#	   );

while(my($mat,$tf)=each(%matrix_tf_pairs)){

    &commentary( "$mat\t$tf\n");

    # get just the mappings for this matrix from the all_mappings file
    my $mappings_file = $mat.".mappings";
    unless( -e $mappings_file ){
	my $command = "grep $mat all_mappings.tab > $mappings_file";
	&backtick($command);
    }

    my ($schema_build) = $enc_db =~ /.*funcgen_(.*)/;
    my $q = "select sr.name,af.seq_region_start,af.seq_region_end,af.score,ft.name,ct.name from annotated_feature af,feature_set fs,feature_type ft,seq_region sr,cell_type ct where ft.class in('Transcription Factor','Insulator') and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id and sr.seq_region_id = af.seq_region_id and sr.schema_build = '$schema_build' and fs.cell_type_id = ct.cell_type_id and ft.name = '$tf'";
    &commentary( "executing :-\n$q\n");
    my $aaref = $dbh->selectall_arrayref($q);
    unless(defined $aaref && @$aaref > 0){die "no data returned by query :\n$q"}

    my $peaks_file = $tf.".real_peaks";
    open(OUT,"> $peaks_file") or die "failed to open $peaks_file";
    foreach my $aref (@$aaref){
	if(@$aref->[1] == 0){
	    warn "start coord = 0 \n ".join("\t",@$aref)."\nFiltered out!\n";
	}else{
 	    print OUT join("\t",@$aref)."\n" or die "failed to write to file";
        }
       
    }
    close(OUT);
    # We single linkage cluster the real peaks
    # before we create the mock peaks    
    &backtick("slink.pl -o temp_$$ -U $peaks_file ");
    &backtick("rm -f $peaks_file ");
    &backtick("mv temp_$$ $peaks_file ");

    my $mock_file = $tf.".mock_peaks";
    my $command = "pwm_make_bed_mock_set.pl -g $genome_descriptor_file  -i $peaks_file > $mock_file";
    &backtick($command);

   
    my $real_olaps_file = $tf."_real_peaks_$mat".".olaps";
    $command = "cooccur.pl -o ".$real_olaps_file." $peaks_file $mappings_file " ;

    &backtick($command);

    my $mock_olaps_file = $tf."_mock_peaks_$mat".".olaps";
    $command = "cooccur.pl -o ".$mock_olaps_file."  $mock_file $mappings_file ";
    &backtick($command);

    my $max = &backtick("cut -f 8 $real_olaps_file | sort -g -u | tail -1");
    chop $max;
    my $n_peaks =  &backtick("cat $peaks_file |wc -l");
    chop $n_peaks;


FIND_THRESH:

    my $ifh;
    open($ifh,$mock_olaps_file) or die "failed to open $mock_olaps_file";

    $verbose = 3;	    
       
    my %peaks;
    
    my @mappings; # array of arrays of chr,start,end,log odds score

    while(my $line = <$ifh>){
	chop $line;
	my @field = split("\t",$line);
	my @mapping;
	push @mapping,@field[0..2];
	push @mapping,@field[7];
	
	push @mappings, \@mapping;
    }
    close($ifh);

    &commentary( "$tf\t"); 
    my ($perc_of_max_thresh, $score_thresh) =
               &mock_perc_score_thresh(\@mappings,$n_peaks,$perc_thresh,$max);

    print $ofh "$mat\t$tf\t$perc_of_max_thresh\t $score_thresh\n";

}





$dbh->disconnect;
exit;


 
###################################################################
# $max is the max log odds score for the factor's PWM
# $perc_thresh is the maximum percentage of mock peaks we want to contain a peak# $npeaks is the total number of peaks for the factor
# $aaref contains the mappings to the mock peaks as an array of arrays of
# chr,start,end,log_odds_score
sub mock_perc_score_thresh{
    my($aaref,$n_peaks,$perc_thresh,$max)=@_;

    my $max_allowed = $n_peaks * $perc_thresh /100;
    my $decr = 0.1;
    for(my $perc=100;$perc>0;$perc -= $decr){
        my $thresh = $max * $perc/100;
	my %peaks;
	foreach my $aref (@$aaref){

            if($aref->[3] >= $thresh){
	        my $peak = join('~',$aref->[0],$aref->[1],$aref->[2]);
	#	    print $peak." ".$aref->[3]." $thresh\n";
	  	$peaks{$peak}++;

	    }
	}

        my $peak_count = scalar(keys(%peaks));
	print STDERR "$perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose >2;
        if($peak_count >= $max_allowed){
	    print STDERR "$perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose >1 ;
	    if($peak_count *100/$n_peaks > $perc_thresh+1){
		$perc += $decr;
		$perc = ($perc > 100)? 100:$perc;
		$thresh = $max * $perc/100;
	    }

	    print STDERR "$perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose > 1;


	    return($perc,$thresh);


	}
    }

if($verbose == 99){ 
    die "Either generate more mappings for this factor or remove it from this run";
}
print "PROBLEM: This factor does not reach threshold \n\n";
$verbose = 99; 
# run the procedure with high verbose level to generate diagnostic output
&mock_perc_score_thresh($aaref,$n_peaks,$perc_thresh,$max)

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



# execute
#
#  Arg [1]   : scalar database handle from DBI module
#  Arg [2]   : array of scalars containing lines of text in the format of SQL statements
#  Function  : submits each element of Arg2 to the $dbh->prepare and 
#              $dbh->execute methods and checks for successful execution. 
#              Returns 0 on failure, 1 on success. Emits error messages 
#              via &err.
#  Returntype: int
#  Exceptions: none
#  Example   : execute($dbh, @array);

sub execute{
    my $dbh = shift;
    my (@array)=@_;

    my $sth;

    &commentary("processing SQL.");

    foreach my $query(@array){
	
    	&commentary(".");

	    unless($sth=$dbh->prepare($query)){
               &err("Error: preparation of statement failed on: $query\n");
               &err("database error message: ".$dbh->errstr."\n");
               return(0);
	    }

        unless($sth->execute){ # returns true on success
            &err("Error: statement execution failed on: $query\n");
            &err("statement handle error message:".$sth->errstr."\n");
            return(0);
	    }
    }

    &commentary("\n");

    return(1);
}


# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print  "$_[0]";
}



sub get_list_from_file{
    my $file=shift;

    open(IN, "< $file") or die "couldn't open list file $file";

    my %ret;
    while( <IN> ){
        chop;
	if($_ =~ '#'){next}
	my @field = split("\t", $_) ;
	my @matrices = split(';',$field[1]);
	foreach my $mat (@matrices){
	    if(exists $ret{$mat}){ 
                warn "WARNING same matrix has multiple TFs\n".
                    $mat."\t".$field[0]."\t".
		    $mat."\t".$ret{$mat}." only one TF's peaks will be used\n";
	    }

	    $ret{$mat} = $field[0];
	}
    }
    close(IN);
    return %ret;
   
}


sub config{
 
($user =     $ENV{'ENSMARTUSER'}) or return(0); 
($host   =   $ENV{'ENSMARTHOST'}) or return(0); 
($port =     $ENV{'ENSMARTPORT'}) or return(0); #usually 3360
($driver  =  $ENV{'ENSMARTDRIVER'}) or return(0); # mysql
($password = $ENV{'ENSMARTPWD'}) or return(0); #
 
}
   
sub err{
    print STDERR "$_[0]\n";
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



    if (exists $opt{H}){
        $host = $opt{H}; 
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

    if (exists $opt{o}){
        $outfile = $opt{o};
    }  


    if (exists $opt{i}){
        $infile = $opt{i};
    }  


    if (exists $opt{g}){
	$genome_descriptor_file = $opt{g};
    }

    if  (exists $opt{e}){
        $enc_db = $opt{e};
    } 

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    pwm_filter_mappings.pl [-h] for help
                  [-e] <db_name> ensembl funcgen database name
                  [-H] <host machine> eg ens-genomics1
                  [-u] <database user> 
                  [-o] <output file> - name of a file for output of thresholds
                  [-i] <input file> - list of TF_name=>pfm_file pairings
                  [-p] <mysql password> 
                  [-P] <mysql port> 
                  [-g] <input file> - containing fasta ID lines for the genome 



END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
