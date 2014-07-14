#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Determines the log odds score threshold for filtering PWM mappings to be incorporated into the functional genomics databases.

=head1 USAGE

see

ensembl-funcgen/scripts/pwm_mapping_notes

grep '>'  /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.fasta > /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines

=head1 EXAMPLES

pwm_filter_mappings.pl -i list -e dev_homo_sapiens_funcgen_58_37c -H ens-genomics1 -g /data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines

pwm_filter_mappings.pl -i list -e dev_homo_sapiens_funcgen_60_37e -H ens-genomics1 -P 3306 -u ensadmin -p big_secret -o thresholds 

=head1 SEE ALSO

ensembl-funcgen/scripts/pwm_filter_mappings.pl
ensembl-funcgen/scripts/pwm_genome_map.pl
ensembl-funcgen/scripts/pwm_make_bed_mock_set.pl

=head1 TO DO

add the functionality from pwm_filter_mappings.pl to the end of this script

=cut


use strict;
use DBI;
use Getopt::Std;
use IO::Handle;
use IO::File;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd
                                               run_backtick_cmd );

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port, $mapping_tabfile);
my $outfile='';
my $infile='';
my $verbose = 2;
my $perc_thresh = 5; # 5% background
my $genome_descriptor_file;# = "/data/blastdb/Ensembl/funcgen/homo_sapiens_male_GRCh37_58_37c_unmasked.id_lines";
my $schema_build;
my %opt;

$| = 1; #no output buffer

if ($ARGV[0]){
&Getopt::Std::getopts('t:u:p:e:H:h:o:i:P:g:s:', \%opt) || &help_text("Invalid Argument") ;
}else{
&help_text; 
}



my $enc_db;
&process_arguments;

my $ofh;
if($outfile){
     #open($ofh,">$outfile") or die "couldn't open file $outfile\n";
     
     open($ofh,">>$outfile") or die "couldn't open file $outfile\n";
     
}else{
        #or write to STDOUT
     $ofh = new IO::File;
     $ofh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}



# hook up with the server
#Hardcode this for the moment...
$driver="mysql";
my $dbh = &make_contact($enc_db);

my %matrix_tf_pairs= &get_list_from_file($infile);

#This currently write and read mappings files from cwd, so has to be run from the dir!!??


#We could make this recoverable by skipping if we already have the threshold available in the file


#Need to tidy up peak files!!
#do we even need to sort here? Output from Moods should already be sorted.
#although we need the peaks and the MF sorted in the same order.
#How is the sort handling primes?
#mappings file do not appear to have +ve stranf then -ve strand mapping
#cooccur handles the sorting

MATRIX: while(my($mat,$tf)=each(%matrix_tf_pairs)){

  &commentary("$mat\t$tf\n");
  
  
  if(-e $outfile){
    my $grep_cmd    = "grep -E '^$mat\[\[:space:\]\]' $outfile";
    my $thresh_line = run_backtick_cmd($grep_cmd, 1);#grep returns exit code 1 if not present
  
    if($thresh_line){
      warn "Skipping $mat as is already present in output:\n$thresh_line";
      next MATRIX;  
    }
  }

  # get just the mappings for this matrix from the all_mappings file
  my $mappings_file = $mat.".mappings";
    
  #This grep is VERY SLOOOOOOOOOOWWWWWWWW
  #the pwm_genome_map.pl script could probably write to separate files
  #by bsubbing individual pwm alignment jobs
  #Where are these removed? Are we duplicating footprint again?
  #This wrote an empty file for the first matrix. likely issues with runtime?
  
  if($mat eq 'MA0065.1'){
    warn "Skipping MA0065.1 for now as is rerunning grep";
    next MATRIX;
  }
  elsif(! -e $mappings_file ){
    
    warn "Skipping $mat $tf due to absent mappings file:\$mappings_file";
    next MATRIX;
    
    
    # -F = fixed string
    #LC_ALL changes locale to C, which forces ASCII character set (128 characters) rather than UTF-8's > 110,000 
    my $command = "LC_ALL=C grep -F $mat $mapping_tabfile > $mappings_file";
    run_system_cmd($command);

    #Hack to deal with a special case that messes things up!
    if($mat eq ' PB0182.1'){
      run_system_cmd("awk '\$5 > 9' < PB0182.1.mappings > PB0182.1.mappings.filtered"); #Prime less than 9? wtf, was this supposed to filter score?
      run_system_cmd("mv PB0182.1.mappings.filtered PB0182.1.mappings");
    }
  }
 
  #elsif($mat eq 'MA0341.1'){
  #   warn "Skipping MA0341.1 as does not meet threshold";
  #  next MATRIX;
  #  

  #my ($schema_build) = $enc_db =~ /.*funcgen_(.*)/;
  if(!defined($schema_build)){ die "Need a schema build!"; }
  #schema build better passed as input as db may not hve been "transformed"
    
  #This should use the API!
  #restriction to schema build is killing this query as the DB was copied before the update script was run
  #API would automatically select the current coord_system
  
  #New query
  my $q = 'select sr.name, af.seq_region_start, af.seq_region_end, fs.name, af.score, af.seq_region_strand '.
    'from annotated_feature af, feature_set fs, seq_region sr, feature_type ft '.
    'where ft.class in(\'Transcription Factor\', \'Insulator\', \'Transcription Factor Complex\') '.
    'and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id '.
    'and sr.seq_region_id = af.seq_region_id and sr.schema_build = '.
    "'$schema_build' and ft.name = '$tf' and af.seq_region_start>0";
  #Add sort in here, but needs to be matched to unix sort of MF mappings
  
#select sr.name, af.seq_region_start, af.seq_region_end, fs.name, af.score, af.strand from annotated_feature af, feature_set fs, seq_region sr, feature_type ft where ft.class in("Transcription Factor", "Insulator", "Transcription Factor Complex") and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id and sr.seq_region_id = af.seq_region_id and sr.schema_build ='76_38' and ft.name = 'Egr1'
  
  my $peaks_file = $tf.".real_peaks";
  #temp hack to get this thing shifting
  my $cmd = "mysql --skip-column-names -hens-genomics2 -uensro -e\"$q\" nj1_tracking_homo_sapiens_funcgen_76_38 | sort -k1,1 -k2,2g -k3,3g > $peaks_file";

  #old query  
  #Half of this is being ignored at present as peaks file only has sr, start and end?
  #unecessary joing to celltype?
  #$q = 'select sr.name,af.seq_region_start,af.seq_region_end,af.score,ft.name,ct.name '.
  #  'from annotated_feature af,feature_set fs,feature_type ft,seq_region sr,cell_type ct'.
  #  'where ft.class in("Transcription Factor", "Insulator", "Transcription Factor Complex") '.
  #  'and fs.feature_type_id = ft.feature_type_id and af.feature_set_id = fs.feature_set_id '.
  #  'and sr.seq_region_id = af.seq_region_id and sr.schema_build = '.
  #  "'$schema_build' and fs.cell_type_id = ct.cell_type_id and ft.name = '$tf'";
    
  #Change this to bed format so we can use bedtools!
  #Also need to change mappings file to bed format, in parse_out_to_tab/bed  
    
  #&commentary( "executing :-\n$q\n");
  #warn $cmd;
  run_system_cmd($cmd);

  
=pod

  my $aaref = $dbh->selectall_arrayref($q);
  unless(defined $aaref && @$aaref > 0){die "no data returned by query :\n$q"}


  open(OUT,"> $peaks_file") or die "failed to open $peaks_file";
  
  #put this filter in the query!!! so we can dump straight to a sorted file
  #instead of bloating the memory and iterating
  
  foreach my $aref (@$aaref){
      
	  if($aref->[1] == 0){
	    warn "start coord = 0 \n ".join("\t",@$aref)."\nFiltered out!\n";
    }
    else{
 	    print OUT join("\t",@$aref)."\n" or die "failed to write to file";
    }
  }

  close(OUT);
  
=cut
  
  # We single linkage cluster the real peaks
  # before we create the mock peaks    
  #&backtick("slink.pl -o temp_$$ -U $peaks_file ");
  run_system_cmd("~dz1/.local/bin/bedtools merge -i $peaks_file > temp_$$");
 # &backtick("rm -f $peaks_file ");
  run_system_cmd("mv temp_$$ $peaks_file");
  
  my $mock_file = $tf.".mock_peaks";
  my $command = "pwm_make_bed_mock_set.pl -g $genome_descriptor_file  -i $peaks_file | sort  -k1,1 -k2,2g -k3,3g > $mock_file";
  warn $command;
  run_system_cmd($command);

  #todo pre-sort mappings files, so we can specifiy -sorted in bedtools commands
  

  #These are dying due to using 37GB! due to unsorted input?

  my $real_olaps_file = $tf."_real_peaks_$mat".".olaps";
  #$command = $ENV{SRC}."/ensembl-funcgen/scripts/miscellaneous/cooccur.pl -o ".$real_olaps_file." $peaks_file $mappings_file " ;
  $command = "~dz1/.local/bin/bedtools intersect -sorted -wa -wb -a $peaks_file -b $mappings_file > $real_olaps_file";
  warn $command;
  run_system_cmd($command);

  my $mock_olaps_file = $tf."_mock_peaks_$mat".".olaps";
  #$command = $ENV{SRC}."/ensembl-funcgen/scripts/miscellaneous/cooccur.pl -o ".$mock_olaps_file."  $mock_file $mappings_file ";
  $command = "~dz1/.local/bin/bedtools intersect -sorted -wa -wb -a $mock_file -b $mappings_file > $mock_olaps_file";
  warn $command;
  run_system_cmd($command);


  #Here, need to handle format difference
  #Cut the score, this appear to be cutting the strand, but there was an empty field
  #-u was is redundant here?
  my $max = run_backtick_cmd("cut -f 8 $real_olaps_file | sort -g | tail -1");
  
  #New format
  #my $max = run_backtick_cmd("cut -f 11 $real_olaps_file | sort -g -u | tail -1");
  chop $max;
  my $n_peaks = run_backtick_cmd("wc -l $peaks_file");
  chop $n_peaks;


  #FIND_THRESH:

  my $ifh;
  open($ifh, $mock_olaps_file) or die "failed to open $mock_olaps_file";
  $verbose = 3;	    
  my @mappings; # array of arrays of chr,start,end,log odds score

  while(my $line = <$ifh>){
    chop $line;
    my @field = split("\t",$line);
    push @mappings, [@field[0..2], $field[7]];
    
    #new format
    #push @mappings, [@field[0..2], $field[10]];
  }

  close($ifh);

  &commentary( "$tf\t"); 
  my ($perc_of_max_thresh, $score_thresh) =
   identify_background_threshold(\@mappings,$n_peaks,$perc_thresh,$max);

  print $ofh "$mat\t$tf\t$perc_of_max_thresh\t$score_thresh\n";
}


$dbh->disconnect;


exit;


 
###################################################################
# $max is the max log odds score for the factor's PWM
# $perc_thresh is the maximum percentage of mock peaks we want to contain a peak
# $npeaks is the total number of peaks for the factor
# $aaref contains the mappings to the mock peaks as an array of arrays of chr,start,end,log_odds_score
sub mock_perc_score_thresh{
  my ($aaref, $n_peaks, $perc_thresh, $max) = @_;
  my $max_allowed = $n_peaks * $perc_thresh /100;
  my $decr        = 0.1;
  
  
  
  #Holy chao! This is a 1000 iterations!
  #let's be a bit more clever about this, and either do it binary chop style
  #or at least do initial larger decrements?
  #current percentages are showing anywhere between 33 and 100%
  #are the 100% pfms really meeting the 5% background?
  
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
    
    
    warn "$peak_count >= $max_allowed";
    
    if($peak_count >= $max_allowed){
      print STDERR "CALCULATING FINAL TRESH $perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose >1 ;
	    
      if($peak_count *100/$n_peaks > $perc_thresh+1){
        $perc  += $decr;
        $perc   = ($perc > 100)? 100:$perc;
        $thresh = $max * $perc/100;
	    }

      print STDERR "FINAL $perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose > 1;
      return ($perc,$thresh);
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


sub identify_background_threshold{
  my ($aaref, $n_peaks, $perc_thresh, $max) = @_;
  my $max_allowed = $n_peaks * $perc_thresh /100;
  #my $decr        = 0.1;
  my $decr = 5;
 
  #Holy chao! This was a 1000 iterations!
  #let's be a bit more clever about this, and either do it binary chop style
  #or at least do initial larger decrements?
  #current percentages are showing anywhere between 33 and 100%
  #are the 100% pfms really meeting the 5% background?

  THRESH: for(my $perc=100; $perc>0; $perc -= $decr){
    warn "prec is now $perc";
    my $thresh = $max * $perc/100;
    #my %peaks;
    my $prev_peak = '';
    my $peak_count = 0;

    foreach my $aref (@$aaref){

      if($aref->[3] >= $thresh){
        my $peak = join('~', ($aref->[0],$aref->[1],$aref->[2]) );
        #     print $peak." ".$aref->[3]." $thresh\n";
        #$peaks{$peak}++;
        #These are not NR peaks, as we get get a peak union
        #for each MF overlap
        #we could easily change this to a count by comparing to the previous key
        
        
        if($peak ne $prev_peak){
          $peak_count++;
          $prev_peak = $peak;  
        }
      }
    }

    #my $peak_count = scalar(keys(%peaks));
    #print STDERR "$perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose >2;
    
    
    #warn "$peak_count >= $max_allowed";
    
    if($peak_count >= $max_allowed){
      
      if(($decr != 0.1) && $perc != 100){
        #reset loop from previous iteration with finer grain decrements  
        $perc += $decr;
        $decr = 0.1;
        warn "Resetting: \$perc to $perc and \$decr to $decr"; 
        next THRESH;
      }
         
      print STDERR "CALCULATING FINAL TRESH $perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose >1 ;

      if($peak_count *100/$n_peaks > $perc_thresh+1){
        $perc  += $decr;
        $perc   = ($perc > 100) ? 100 : $perc;
        $thresh = $max * $perc/100;
      }

      print STDERR "FINAL $perc $max_allowed $peak_count ".$peak_count *100/$n_peaks." $thresh\n" if $verbose > 1;
      return ($perc,$thresh);
    }
  }

  if($verbose == 99){ 
    die "Either generate more mappings for this factor or remove it from this run";
  }
  
  print "PROBLEM: This factor does not reach threshold \n\n";
  $verbose = 99; 
  # run the procedure with high verbose level to generate diagnostic output
  identify_background_threshold($aaref,$n_peaks,$perc_thresh,$max)
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

    if (exists $opt{t}){
        $mapping_tabfile = $opt{t}; 
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

    if  (exists $opt{s}){
        $schema_build = $opt{s};
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
