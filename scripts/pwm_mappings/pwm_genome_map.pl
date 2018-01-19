#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

# To do - Some of these apply to all the pwm scripts
# 1 Getopt::Long, remove process_arguments
# 2 remove backticks in favour of system. Use EFGUtils::run_system_cmd
# 3 Remove commentary and err subs and flip use STDERR(warn) and STDOUT(print). 
#   We should have nothing on STDERR unless there is a warning/error.
# 4 DONE fastaclean path? Now append path belowb
# 5 Fix backtick which tests do not work
# 6 remove config sub
# 7 Reduce memory footprint! Currently using > 30GB with Jaspar 5
# 8 Parallelise & make recoverable
# 9 Restrict to motifs we are interested
# 10 which_path for run_moods.pl is checking /usr/local/ensembl/bin which is not apparently in the $PATH?
# 11 Move all batch management/merge to Hive and simplify this script to iterate over the target files
#    or fail. Move as much as poss to MotifTools

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

use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools   qw( explode_fasta_file );
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( read_matrix_file get_revcomp_file_path );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( run_backtick_cmd run_system_cmd which_path );

$ENV{PATH} .= ':/usr/local/ensembl/bin/'; #fastaclean home, probably needs moving to /software/ensembl/bin

my ($user, $password, $driver, $host, $port, $out_dir);
my $outfile='';
my $infile='';
my $sp;
my $verbose = 0;
my $score_thresh_perc = 0.7; # 70% of max
my $pwm_type = 'jaspar';

my $genome_file;
my $thresh = 0.001;
my $assembly;

my %opt;

if(! @ARGV){
  help_text();
}

print "$0 @ARGV\n";
Getopt::Std::getopts('a:s:t:g:h:v:o:i:p:', \%opt) || die ;


# get configuration from environment variables
#&config; # this may fail but config can be on command line


process_arguments();
 
my @pwm_files    = @ARGV;
my $work_dir     = $out_dir.'/tmp_results';
my $lsf_out_dir  = $work_dir.'/lsf_out';
#which_path as run_moods.pl is likely in a subdir of the $PATH
#Do it here, as it's a bit slow to do it in a loop
print "Finding run_moods.pl path...\n";
my $moods_path   = which_path('run_moods.pl');
#Now allowing recovery
#To force, simply delete the relevant output manually?
#if($force){
#  run_system_cmd("rm -rf $work_dir" if -e $work_dir;
#}

run_system_cmd("mkdir -p $lsf_out_dir") if ! -e $lsf_out_dir;

#This will be validated before we run moods, 
#but not before we preprocess the matrices and submit the jobs
#So leave here for now, until we add a DBI or DBConnection helper method to DBAdaptorHelper?
#based on params hash. 


#my $jdbh;
#my $dsn = sprintf( "DBI:%s:%s:host=%s;port=%s",
#                     'mysql', $jdb_name, $jdb_host,
#                     3306 );

#eval {
#  $jdbh = DBI->connect( $dsn, $jdb_user, undef,
#                          { 'RaiseError' => 1 } );
#  };

#my $error = $@;

#if ( ! $jdbh || $error || ! $jdbh->ping ) {
#  die( "Could not connect to database using [$dsn] as a locator:\n" . $error );
#}
  

# we need to explode the genome.fasta file into individual sequences.
# and if abbreviated chr names have been requested change the file names
# as these are used in the output file
# also if the sequences contain characters other than [acgtnACGTN] they need
# to be converted to N otherwise find_pssm_dna outputs the wrong coords.
#$verbose = 2;

#Move genome files outside of workdir
#As these are reference files which maybe re-used.
my @chr_files = @{explode_fasta_file($genome_file, $out_dir.'/genome', $assembly, 1)};


# irrespective of pwm_type we need to create the rev-comp matrices and put all the matrices in a working directory - along with a composite matrix_list.txt file
my %file_max_score;
my $pfm_job_info = {};

if($pwm_type eq 'jaspar'){

  foreach my $file (@pwm_files){
    run_matrix_alignments($file, $pfm_job_info);
  }
}
#elsif($pwm_type eq 'transfac'){
#  # if the pwm_type isn't Jaspar we need to convert them into Jaspar pfm 
#  # format and create a matrix_list.txt file. Put new PWM files etc. in 
#  # working directory
#  # Transfac matrixes all come in a single file 
#  if(@ARGV > 1){
#    die "Transfac matrices should all be in a single file, not ".join("\n",@ARGV)."\n";
#  }
#
#  foreach my $file (parse_matrixdat_2_pfm($ARGV[0],$work_dir)){
#    run_matrix_alignments($file, $pfm_job_info);
#  }
#
#  #This should have really been done upstream of here!
#  #as the jaspar pmf files are not generated
#}
else{
  die "Can only handle Jaspar and Transfac formats at present\n";
}


### Wait for jobs arrays and merge files ###

my $sleep_interval = 10; #seconds
#Start low as most will have been completed by the time all the 
#jobs have been submitted, and the loop itself will take some more time
my $submitted = grep { $_->{job_array_id} } values(%$pfm_job_info);
my $num_still_running = grep { ! $_->{processed} } values(%$pfm_job_info);
my $merged    = grep { $_->{processed} }    values(%$pfm_job_info);



#Also need to count those run but not merged yet
my $seen_failures = 0;


#warn "num still running is $num_still_running";


#This will loop indefinitely if we encounter a non-trnasient job state which is
#not DONE or EXIT e.g. UNKWN

#TODO parallelise these merges by fanning these jobs in the hive

while($num_still_running){

  foreach my $pfm_file(keys %$pfm_job_info){
    my $p_info = $pfm_job_info->{$pfm_file};

    if(! $p_info->{processed}){

      if(defined $p_info->{job_array_id}){
        print "Checking batch jobs for $pfm_file (Job ID = ".$p_info->{job_array_id}.")\n";
        my ($finished, $success, $j_info) = bjobs_poll($p_info->{job_array_id});
        print "Finished $finished/".$p_info->{array_size}.". Success for $success\n";

        if($finished == $p_info->{array_size}){
          

          if($success != $finished){ #We have some EXITs
            #Remove the complete jobs
            delete $j_info->{DONE};
            warn 'Found '.($finished - $success)." failed jobs:\t$pfm_file\nSkipping merge...\n";
            $p_info->{failed_jobs} = $j_info;
            $seen_failures++;
          }
          else{ #Submit merge job
            print "Job array completed successfully:\t$pfm_file\nProceeding with merge...\n";
            merge_files($p_info);
            $merged++;
          }

         #We are starting to see arrays complete
          #So merge time is probably long enough to wait before next iteration
          $sleep_interval = 0;
          $submitted--;
          $p_info->{processed} = 1;
          $num_still_running--;
        }
        else{
          $sleep_interval = 60;#assume we have to wait for something now 
          #but this could be an uncaught state.
        }
      }
      else{ #Must have all files, but not merged yet
        $sleep_interval = 0;
        print "Running merge for previously run $pfm_file batch jobs\n";
        merge_files($p_info);
        $p_info->{processed} = 1;
        $merged++;
        $num_still_running--;
      }
    }
    else{ #Is already processed
      next;
    }

    print "$merged merged. $submitted submitted. $num_still_running in total to process.".
      " Sleeping for $sleep_interval seconds...\n";

    #Having sleep here instead of at start means we may poll bjobs before
    #the batch daemon has had chance to process the bsub, and may return nothing
    #This will result in sleep being set to 60 from beginning
    sleep $sleep_interval;   
  }
}

if($seen_failures){ # Summarise failures and die
  warn "FAILURE SUMMARY\n";

  foreach my $pfm_file(keys %$pfm_job_info){
    my $p_info = $pfm_job_info->{$pfm_file};

    if($p_info->{failed_jobs}){
      warn $pfm_file."\n";

      map { warn $_.":\t".join(' ', @{$p_info->{failed_jobs}->{$_}})."\n" } 
       keys %{$p_info->{failed_jobs}};
    }
  }
  die("$seen_failures/".scalar(%$pfm_job_info).' PFMs had failures');
}
else{
  #Tidyup only if we have seen no failures to allow recovery
  #unlink(@chr_files);#We always run fastaexplode, so may aswell
  run_system_cmd("rm -rf $work_dir");
}

sub merge_files{
  my $p_info = shift;

  #cating inline shouldn't take too much time.
  #and data sets will complete at different times
  #so little time lost here
  
  #Sort here as the rev comp mapping will simply be listed after all the normal mappings

  my $cat_cmd = 'cat '.join(' ', @{$p_info->{output_files}}).' | sort -k1,1 -k2,2n -k3,3n > '.$p_info->{merged_file}.'.tmp';
  print "Merging files to:\t".$p_info->{merged_file}."\n";

  #Wrap this in eval too
  run_system_cmd('rm -f '.$p_info->{merged_file}.'.tmp') if -e $p_info->{merged_file}.'.tmp';

  if(my $exit_code = run_system_cmd($cat_cmd, 1)){
    #use no exit flag, so we can process the rest before failing
    my $err = "Failed to run merge command:\n$cat_cmd\nExit code:\t$exit_code";
    $seen_failures++;
    $p_info->{failed_jobs} = {'MERGE_FAIL' => [ $err ]};    
    warn $err."\n";
  }
  else{
    run_system_cmd('mv '.$p_info->{merged_file}.'.tmp '.$p_info->{merged_file});
    run_system_cmd("rm -f ".join(' ', @{$p_info->{output_files}}));
  }

  return;
}


sub bjobs_poll{
  my $job_id       = shift || die("Must provide an LSF job ID to poll");
  my @bjobs_output = run_backtick_cmd('bjobs '.$job_id);
  my %job_info     = (DONE => [], EXIT => []);

  if($bjobs_output[0] !~ /^JOBID/o){
    #Throw here?
    warn "Unexpected output from bjobs:\n".join("\n", @bjobs_output);
  }
  else{
    #CHECK for DONE/EXIT jobs, else assume we are in some
    #sort of RUN or suspended status
    #what about UNKWN?
    shift @bjobs_output; #Remove the header line
    my ($state, $job_name);

    foreach my $job_line(@bjobs_output){
      (undef, undef, $state, undef, undef, undef, $job_name) = split(/\s+/, $job_line);
      #Get and append array index from job name
      my $batch_job_id = $job_id;

      $batch_job_id .= $1 if $job_name =~ /(\[[0-9]+\])$/o;
      $job_info{$state}||= [];
      push @{$job_info{$state}}, $batch_job_id;
    }
  }
  #We need to return total number 'complete'
  #what about unknown jobs? And other states?
  #Also need to return those which have status EXIT or weren't DONE
  return (scalar(@{$job_info{EXIT}}) + scalar(@{$job_info{DONE}}), 
          (scalar(@{$job_info{DONE}})), 
          \%job_info);
}


#is the rc file being passed as an argument to find_pssm_dna

sub run_matrix_alignments{
  my $pfm_file_path = shift;
  my $pfm_job_info  = shift;
  #my $file_max_score = shift;
  #my $work_dir = shift; #global

  my $pfm_file_name = basename($pfm_file_path);  

  #warn "HARDCODED PB0111.1.pfm only for debug";
  #return if $pfm_file_name ne 'PB0111.1.pfm';

  print "Preparing to submit alignment jobs for:\t$pfm_file_name\n";

  #This has already been created in run_pwm_mapping_pipeline.pl
  my $rc_file       = get_revcomp_file_path($pfm_file_path);

  #This is actually no longer used here?
  #my ($min, $max)                 = matrix_min_max($pfm_file_name, $work_dir);
  #$file_max_score{$pfm_file_name} = $max;
  #print $pfm_file_name."\t$min\t$max\n" if $verbose > 1;
  my $format = 'bed';
  (my $all_mappings_file = $pfm_file_name) =~ s/\.pfm$/.unfiltered.${format}/o;
  $all_mappings_file = $out_dir.'/'.$all_mappings_file;

  $pfm_job_info->{$pfm_file_name} = 
   {job_array_id => undef,
    array_size   => undef,
    output_files => [],
    merged_file  => $all_mappings_file}; 
  #This block currently auto-recovers complete files
  #To force a re-run requires deleting the output
  my $seen = 0;

  # TODO Remove job array submission for slices
  # Biggest job only seems to take ~90secs?

  if(! -e $all_mappings_file){
    my @target_files;
    (my $out_prefix = $work_dir.'/'.$pfm_file_name) =~ s/\.pfm$/.unfiltered/o;

    foreach my $chr_file (@chr_files){
      #explode_genome_fasta makes file named after the seq_region_name e.g 5.fa
      (my $chr = basename($chr_file)) =~ s/\.fa(sta)*$//o;
 
      #This has to match what MotifTools::run_moods creates
      my $chr_out_file = $out_prefix.'.'.$chr.'.'.$format;
      push @{$pfm_job_info->{$pfm_file_name}{output_files}}, $chr_out_file;
      #This does not yet support recovery from out file
      #This needs to be implemented in run_moods.pl!

      if(! -e $chr_out_file){     
        push @target_files, $chr_file;
      }
      #else{
      #  #warn "Already seen $chr_out_file\n";
      #  $seen++;
      #}
    }

    my $asize = scalar(@target_files);
    $pfm_job_info->{$pfm_file_name}{array_size}   = $asize;

    if($asize){

      if($asize != scalar(@{$pfm_job_info->{$pfm_file_name}{output_files}})){
        warn 'Already seen '.(scalar(@{$pfm_job_info->{$pfm_file_name}{output_files}}) - $asize).
          '/'.scalar(@{$pfm_job_info->{$pfm_file_name}{output_files}})." expected output files\n";
      }
    
      my $command = $moods_path." --format $format --job_array ".
       "--out_prefix $out_prefix --queries $pfm_file_path $rc_file --targets ".
       join(' ', @target_files);
      #--host $jdb_host --dbname $jdb_name --user $jdb_user ".
      #default --p_thresh is 0.001

      print "Submitting $asize jobs for:\t$pfm_file_name\n";
      $pfm_job_info->{$pfm_file_name}{job_array_id} = bsub_job_array($pfm_file_name, $asize, "$command");
    }
    else{
      warn 'Skipping run_moods.pl submission. Already seen all expected '.
       "intermediate output for ${pfm_file_name}\n";#.join(' ', @{$pfm_job_info->{$pfm_file_name}{output_files}});
    }
  }
  else{
    $pfm_job_info->{$pfm_file_name}{processed} = 1;
    print "Skipping job submission. Already seen final merged output file:\t$all_mappings_file\n";
  }
  
  return $pfm_job_info;
}


sub bsub_job_array{
  my $job_name    = shift;
  my $array_size  = shift;
  my $job_cmd     = shift;
 
  my $bsub_params = '-M 6000 -R"select[mem>6000] rusage[mem=6000]" '.
   "-o $lsf_out_dir/$job_name.\%J.\%I.out -e $lsf_out_dir/$job_name.\%J.\%I.err";
  my $cmd         = "bsub -J \"${job_name}[1-${array_size}]\" $bsub_params '$job_cmd'";
  my $bsub_output = run_backtick_cmd($cmd);
  #print $bsub_output."\n";

  if($bsub_output =~ /Job \<([0-9]+)\> is submitted to /){
    $bsub_output = $1;
  }
  else{  
    #This is slightly brittle, we could do some retry here?
    throw("Failed to submit job:\n$cmd\n$bsub_output");
  }

  return $bsub_output;
}
 
###################################################################

#Move some of these to MotifTools.pm

#Make this use read_matrix_file, or integrate it into read_matrix_file

sub matrix_min_max{
  my $file = shift;

  open(IN,$file) or die "failed to open $file";


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

    my($min_pos,@sums) = get_max_add(@new_mat);
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
	$maxes[$i] = log2($maxes[$i]);
	$min_sum  += log2( ((($list[3])+0.002)/1.008)/0.25 );

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


#What? This is logn not log2!
sub log2{
    my $n = shift;

    return log($n);#/log(2.0);

}


#TODO 
#Move these transfac methods to MotifTools and integrate into read_matrix_file

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
	my $pfm = process_matrix($ifh,$ofh,$ac,$name,$work_dir);
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




sub process_arguments{
    if ( exists $opt{'h'} ){ 
        help_text();
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
        $out_dir = $opt{o};
    }else{
	   help_text("Please give a path for the output directory");
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
      warn "\n".$msg."\n";
    }

    warn <<"END_OF_TEXT";

  pwm_genome_map.pl -g fasta_file -o output_dir [options] pwm_file1 [pwm_file2 pwm_file3 ...]

                  [-h] for help
                  [-a] <assembly_version> eg GCRh37 as used in full chr. name
                  [-t] <pwm_type> 'jaspar'=default, 'transfac' 
                  [-p] <probability> threshold for moods mapper default=0.001
                   -g  <file_name> genome fasta file
                   -o  <output dir> path of a (spacious) output directory 
                  [-v] <integer> verbosity level
                  [-s] <species> eg -smus_musculus, default = homo_sapiens  
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
