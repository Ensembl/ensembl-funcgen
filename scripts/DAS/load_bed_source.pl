#!/usr/bin/env perl

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


=head1 NAME

load_bed_source.pl

=head1 SYNOPSIS

 LoadBedSources|load_bed_source.pl [options]

 e.g. Using the associated efg environment function

 LoadBedDASSources  -host $DB_HOST --user $DB_USER --pass WRITE_PASS --dbname efg_DAS --files ES_DNase_le1m_reads.bed.gz --names ES_DNase --profile --frag_length 150 --bin_size 25 -reads


=head1 DESCRIPTION

This script loads a DAS source from raw read alignments(bed) and/or profiles generated from them.
NOTE: Does not yet support paired end data.

=head1 OPTIONS

 DB Connection
  --host|h
  --port
  --user|u
  --pass|p
  --dbname|d

 Run Modes
  --reads         Flag to load raw read alignments, expects input files are bed
  --profile       Flag to generate and load read profile
  --no_load       Skips load step and just generates profile if specified
  --no_compress   Skips default compression of files after loading 
  --profile_input Skips profile generation, expects input(--files) is profile

 Input
  --files         Space separate list of input files
  --names         Optional space separated list of source names, must preferably end in valid 
                  feature type to enable automatic colour usage e.g. ES_DNase1
  --prefix        Optional prefix for source names
  --assembly      Assembly version, default is current default chromosome version.
  
 Profile Options, only required if --profile specified
  --bin_size      Bin size to compute profile scores over
  --frag_length   Length to extend single reads to
  Need to add paired end support here
		
 Other
  --mysql_sort_buffer Sets the minimum mysql_sort_buffer_size for this session(default=50331648 ~47MB)
                      Also modifies myisam_max_sort_file_size accordingly. Increasing this can help
                      speed up load times in 'Repair by KeyCache' is being used by mysql after the 
                      initial file load. Check the processlist in the mysql client if things are too slow.
  --help              Prints a short help message and exits
  --man|m             Prints the man page

=cut

#To do
#Integrate this with the collection code to get a set of bins
#No need to keep reads/profile file? Can always dump out from DB?
#Implement multiple windows sizes? Or matrix? MySQL COMPRESS?
#Change to use mysqlimport
#Does this work with paired end data???? 
#binary packing/windows? 
#How does this use maxbins?
#Merge this with bsub script
#Validate/warn last token of name wrt feature type colours in config ini.
#Integrate with BED parser


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file);
use DBI;

my ($pass,$host,$user,$dbname,$prefix, $input_file, @files, @names);
my ($no_load, $bin_size, $frag_length, $profile_input, $source_name);
my ($profile, $reads, %output_files, $assembly);
my $port = 3306;
my $no_compress = 0;
my $mysql_sort_buffer = 50331648;#6 * default 8388608  ~8GB >> ~47MB


my $params_msg = "Params are:\t@ARGV";

GetOptions (
    'host|h=s'         => \$host,
    'port:i'           => \$port,
    'user|u=s'         => \$user,
    'pass|p:s'         => \$pass,
    'dbname|d=s'       => \$dbname,
    'files=s{,}'       => \@files,
    'reads'            => \$reads,
    'profile'          => \$profile,
    'profile_input'    => \$profile_input,
    'prefix=s'         => \$prefix,
    'names=s{,}'       => \@names,
    'bin_size=i'       => \$bin_size,
    'no_load'          => \$no_load,
    'no_compress'      => \$no_compress,
	'assembly=s'       => \$assembly,
    'frag_length=i'    => \$frag_length,
    'mysql_sort_buffer=s' => \$mysql_sort_buffer,
    'help|?'           => sub { pos2usage(-exitval => 0, -message => $params_msg);},
    'man|m'            => sub { pos2usage(-exitval => 0, -message => $params_msg, verbose => 2);},
		   ) or pod2usage ( -exitval => 1,
							-message => $params_msg );

if (@ARGV){
  pod2usage( -exitval =>1,
			 -message => "You have specified incomplete options. $params_msg");
}


#Check params

if( ! ($host && $user && $pass && $dbname)){
  warn "$host && $user && $pass && $dbname";
  die("You must provide some DB connection paramters:\t --host --user --pass --dbname [ --port ]\n");
}

if(! ($reads || $profile || $profile_input)){
  die("Must provide at least one of the following run modes e.g. --reads or --profile or --profile_input\n");
}
elsif($no_load && ! $profile){
  die("You have selected options --no_load without specifying --profile, no action taken\n");
}
elsif($profile &&
	  ! ($bin_size && $frag_length)){
  die("You must provide a --bin_size and a --frag_length to generate a profile\n");
}
elsif(($bin_size || $frag_length) &&
	  (! $profile)){
  die("You have specified a --bin_size and/or --frag_length, did you want to load a --profile?\n")
}
elsif($profile_input && $reads){
  die("You have specified mutually exclusive options --profile_input and --reads\n");
}
elsif($profile_input && $profile){
  die("You have specified mutually exclusive options --profile_input and --profile, please select one or the other\n");
}

if(@names &&
   scalar(@names) != scalar(@files)){
  die('You have specified an unmatched number of source names, these must match the number and order of your -files');
}

#Get/Check file
if (exists $ENV{LSB_JOBINDEX}) {
  @files = ($files[$ENV{LSB_JOBINDEX}-1]);
  @names = ($names[$ENV{LSB_JOBINDEX}-1]);
}
else{

  if(scalar(@files) > 1){
	throw('You have specified more than one file, maybe you want to submit this to the farm using run_build_profile.sh|LoadBedDASSources');
  }

  @files = ($files[0]);
  @names = ($names[0]);
}

#We need to set this
$input_file = $files[0];
$source_name = $names[0];

#Set the input file correctly
if($profile_input){
  $output_files{profile} = $input_file;
}
elsif($reads){
  $output_files{reads} = $input_file;
}


if(! -f $input_file){
  die("File does not exist:\t$input_file\nMust provide at least one file path to build a profile");
}

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;port=$port",
					   "$user", "$pass",
					   {RaiseError => 1,
						mysql_auto_reconnect => 1});

#This should use the API
if($assembly){  #Validate assembly
 
  my ($chr_assm) = $dbh->selectrow_array("SELECT distinct(version) from coord_system where name='chromosome' and version='$assembly'");

  if($assembly ne $chr_assm){
	die("The assembly you have specified($assembly) is not present in the $dbname");
  }
}
else{#Get current default
  ($assembly) = $dbh->selectrow_array("SELECT distinct(version) from coord_system where name='chromosome' and is_current=1");
}


my (@bin, $start_bin, $start_bin_start, $end_bin, $end_bin_start, $curr_seq,
	$seq, $read_start, $read_end, $read_length, $ori, $read_extend);

if($profile && ! $profile_input){


  #Build profile
  print ":: Building profile for:\t$input_file\n";
  warn "WARNING:\tThis does not yet support paired end data\n";
  

  #gzip expects .z or .gz
  #This should be in Helper?
  
  my $compressed_data =  `file -L $input_file` or die "Can't execute 'file -L $input_file'";

  my $gzip = 1 if $compressed_data =~ /gzip/;

  my $fh;

  if($gzip){
	#open(FILE, "gzip -dc $input_file | sort -k 1,3 |") or die ("Can't open compressed file:\t$input_file");
	$fh = open_file($input_file, "gzip -dc %s | sort -k 1,3 |");
  }
  elsif($compressed_data =~ /compressed/){
      die("This script only handles gzip compressed files, please uncompress $input_file manually before rerunning");
  }
  else{
	#open(FILE, "sort -k 1,3 $input_file |") or die ("Cannot open file:$input_file");
	$fh = open_file($input_file, "sort -k 1,3 $input_file |");
  }


  my $binsize = sprintf("%03d", $bin_size);

  #We need to validate this input file name
  #Let's not depend on _reads.bed


  my $output_file = $input_file."_profile_${binsize}";
  #Doing this will prevent us from having : in the filename
  #Either we don't allow this or we always add the assembly?
  #This is also growing the table name...we need another registry table to manage this
  #This would be part of DAS subset of tables, why bother?  Let's just have them all and just use das table
  #maybe have meta entry is_das_db? or is_only_das_db
  #The only down side to this is that people won't be able to run script in isolation to create ensembl independant files
  #Will have to use API to display. Maybe we could just have an option to produce non seq_region flat files?
  #Would load direct into result_feature as result feature_set.
  #Then we can just reimplement result_feature table as matrix files.
  $output_file .= "_$assembly";# if $assembly;
  
  #Need to test if compressed here!


  #No need to compress file here as we are most likely just going to decompress it straight away to load
  #open(OUT, "| gzip -c > $output_file")
  #Actually should compress by default if no load?

  open(OUT, "> $output_file") or throw ("Can't open out file $output_file");


  #Can we test size of file here?
  #And print a progress counter?
  #This need to move to the Bed parser and use the collection code


  while (<$fh>) {
    chomp;#should be chump for windows safety?
	
	#undef is score
	#Not dealing with any other trailing fields here
	($seq, $read_start, $read_end, undef, $ori) = split/\t/, $_;
	$read_length = $read_end-$read_start+1;

	if (defined $curr_seq && ($curr_seq ne $seq)) {
	  &write_bins();
	  @bin = ();
	}
	
	$curr_seq = $seq;
    
    die("read is longer ($read_length) than specified fragment length ($frag_length)") if ($frag_length<$read_length);

	$read_extend = $frag_length-$read_length;
	
	# extend reads to given fragment length
	if ($ori eq '+') {
	  $read_end+=$read_extend;
    } else {
	  $read_start-=$read_extend;
        $read_start=1 if ( $read_start < 1 );
    }

    # update read length
    $read_length = $read_end-$read_start+1;
    
    # determine bins that are covered by the read and add 
    # coverage to bin score
	$start_bin = sprintf("%d", ($read_start-1)/$bin_size);
    #start pos of start bin
    $start_bin_start = ($start_bin*$bin_size)+1;

    $end_bin = sprintf("%d", ($read_end-1)/$bin_size);
    #start pos of end bin
    $end_bin_start = ($end_bin*$bin_size)+1;

    #printf "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", 
    #$seq, $read_start, $read_end, $read_length,	$ori,
    #$start_bin, $start_bin_start, 
    #$end_bin, $end_bin_start;

    if ($start_bin == $end_bin) { ## read start and end together in one bin
	  $bin[$start_bin] += $read_length/$bin_size;
	} 
	else {
	  $bin[$start_bin] += (($start_bin_start+$bin_size)-$read_start)/$bin_size;
        #print "(($start_bin_start+$bin_size)-$read_start)/$bin_size = "
        #    .(($start_bin_start+$bin_size)-$read_start)/$bin_size."\n";
        
        for (my $i=$start_bin+1; $i<$end_bin; $i++) {

            #print $i, "\n";
            $bin[$i]++;

        }

        $bin[$end_bin] += ($read_end-$end_bin_start+1)/$bin_size;
        #print "($read_end-$end_bin_start+1)/$bin_size = "
        #    .($read_end-$end_bin_start+1)/$bin_size."\n";
	  
    }

  }

  &write_bins;

  close($fh);
  close OUT;


  
  $output_files{profile} = $output_file;
  #push @files, $out;
  #Okay this is not working with the names array!
  
  


}


if( ! $no_load){

  #warn("No Hydra source name prefix specified!\n") if (! $prefix);

  #my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;port=$port",
  #						 "$user", "$pass",
  #						 {RaiseError => 1,
  #						  mysql_auto_reconnect => 1});
  



  my @sort_vars= $dbh->selectrow_array(q{select @@myisam_sort_buffer_size});#, @@myisam_max_sort_file_size});
  
  if($sort_vars[0] < $mysql_sort_buffer){
      
      #Could we generate this based on the actual file sizes and the suggested defaults ratios?

      #Generate the ratio of increase
      #my $mysql_sort_file = $sort_vars[1] * ( $mysql_sort_buffer / $sort_vars[0]);
      print ":: Setting MySQL session vars for optimised load time:\tmyisam_sort_buffer_size=$mysql_sort_buffer\n";#\tmyisam_max_sort_file_size=$mysql_sort_file\n";
      $dbh->do("SET SESSION myisam_sort_buffer_size=$mysql_sort_buffer");  #This is GLOBAL and SESSION
      #$dbh->do("SET GLOBAL myisam_max_sort_file_size=$mysql_sort_file");  #This is GLOBAL hence should not set unless we can reset it afterwards safely
      #Cannot do this as script may exit before we get a chance to reset

      #Need to check if actually reset

      @sort_vars= $dbh->selectrow_array(q{select @@myisam_sort_buffer_size});#, @@myisam_max_sort_file_size});

      if($sort_vars[0] != $mysql_sort_buffer){ 
	  warn ":: WARNING\tCould not reset misam_sort_buffer_size to $mysql_sort_buffer, using ".$sort_vars[0]."\n";
      }

      #As we cannot rely on reseting the max sort file size here we should warn that it will take a long time

      
      #if($sort_vars[1] != $mysql_sort_file){                                                                                                                           
      #    warn ":: WARNING\tCould not reset misam_max_sort_dile_size to $mysql_sort_file, using ".$sort_vars[1]."\n";                                                    
      #}   

      #real    real    197m51.489s  Using defaults
      #real    36m54.789s           With max size rest but buffer on default? No zip unzip between profile generation and load. Is this using file sort?


  }

  #Now check the @@myisam_max_sort_file_size is not smaller than the file size
  #This relates to the mysql file size, not the input file size, but we can use this as a rough guide
  @sort_vars= $dbh->selectrow_array(q{select @@myisam_max_sort_file_size}); 
  my $max_sort_size = $sort_vars[0]/1024;
  


  #Validate/identify file type, not by name!
  foreach my $type(keys %output_files){
	my $table_name = $source_name;
	$table_name = "${prefix}_${table_name}" if $prefix;
	$table_name = "bed_${type}_${table_name}";

	my $file = $output_files{$type};
	

	my $compressed_data =  `file -L $file` or die "Can't execute 'file -L $file'";  
	my $gzip = 1 if $compressed_data =~ /gzip/;

	my $link;
	$link = 1 if -l $file;
	
	if ($gzip) {
	  print ":: Decompressing:\t$file\n";

	  #Get the suffix?

	  #Backup link
	  if($link){
		system("cp -fP $file ${file}.backup") == 0 or die('Failed to backup link');
	  }
	  

	  system("gzip -df $file") == 0
	      or die "Can't decompress file $file";
	  $file =~ s/\.gz$//;
	  die("$file does not exit, expected suffix is .gz") if (! -f $file);
	}
	elsif($compressed_data =~ /compressed/){
	    die("This script only handles gzip compressed files, please uncompress $file manually before rerunning:\t$compressed_data\n");  
	}

	print ":: Loading $type file:\t$file\n";    
	my $file_size = -s $file;

	if ($file_size  < $max_sort_size){
	    
	    warn "WARNING: $file is larger than myisam_max_sort_file_size($file_size > $max_sort_size).\tWARNING: This may slow down import due to using 'Repair by KeyCache' rather than 'Repair by sorting'. If your import is taking a long time(>30 mins for ~3 million reads) check the repair type on the mysql process list.\nYou may need to get you DB admin to increase the GLOBAL session variable (N>".($file_size*1024).")e.g. mysql> SET GLOBAL myisam_max_sort_file_size=N\n";
	}

	
	#Can we split this into something more readable/useable
	#We need to be able to identify these table in an funcgen schema
	#Stanard prefix
	#bed_reads|profile_PREFIX_NAME

	#Maximal match/remove path and bed suffix
	
	#if(! $name){
	#  ($name=$file) =~ s,^(.*/)?(.+)\.bed$,$2,;
	#  #remove .'s for MySQL
	#  $name =~ s,\.,_,g;
	#}

	#$name = $prefix.'_'.$name if($prefix);
	#$name = 'bed_'.$format.'_'.$name;

	if(length($table_name) >64){
	  die("Table name exceeded MySQL maximum of 64 characters:\t$table_name\n".
		  'Please rename your file or choose a shorter --prefix or --names to rectify');
	}

	print ":: Table name:\t\t$table_name\n";

	my $sth = $dbh->do("DROP TABLE IF EXISTS `$table_name`;");
		
	$sth = $dbh->do("CREATE TABLE `$table_name` (
    `feature_id`    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
    `seq_region`    VARCHAR(20) NOT NULL,
    `start`         INT(10) UNSIGNED NOT NULL DEFAULT '0',
    `end`           INT(10) UNSIGNED NOT NULL DEFAULT '0',
    `name`          VARCHAR(40) NOT NULL,
    `score`         FLOAT NOT NULL DEFAULT '0',
    `strand`        ENUM('0','+','-') DEFAULT '0',
    `note`          VARCHAR(40) DEFAULT NULL,
    PRIMARY KEY     (`feature_id`),
    KEY `seq_region_idx` (`seq_region`, `start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;");

#latin1_bin is case sensitive collation! Why do we need this?

  
  #This needs to change to mysqlimport! In case mysql server does not have access to local file
  
  if ($type eq 'reads') {
	
	#what is @mm for?

      #Wouldn't it be faster to get unix to do this replace?
      #Also mysql automatically disable keys when loading data infile into an empty table
      #once loaded it will attempt to to Repair with keycache
      #However Repair with filesort would be 100-1000Xs faster
      #Although doesn't work with primary/unique keys? Is this an issue?
      #It chooses which method to use based on your myisam_sort_buffer_size, 
      #myisam_max_sort_file_size and myisam_max_extra_sort_file_size.  Have 
      #you increased the size of these?  Keep in mind these are SESSION 
      #variables, so they can be set on the connection right before you LOAD 
      #DATA INFILE.

      # Also keep in mind that Repair by sort doesn't work for UNIQUE or 
      # PRIMARY KEYs.

      #mymac - default setting
      #| myisam_max_sort_file_size | 2146435072 | # 2GB
      #| myisam_sort_buffer_size   | 8388608    | #Which is about 8MB..this is tiny!
      #Could probably up both of these default, but should we just set it temporarily for this connection at 5X?


      #ens-genomics1
      #| myisam_max_sort_file_size | 9223372036853727232 | #8589934590 GB??? wtf?
      #| myisam_sort_buffer_size   | 67108864            |  #This is only about 60 MB?
      #


	
	$sth = $dbh->do("LOAD DATA LOCAL INFILE '$file' INTO TABLE $table_name 
               (seq_region,start,end,name,\@mm,strand,score) 
               set seq_region=replace(seq_region, 'chr', ''), note=concat('mm=',\@mm);");
	
  } 
	elsif ($type eq 'profile') {
	  
	  $sth = $dbh->do("LOAD DATA LOCAL INFILE '$file' INTO TABLE $table_name 
               (seq_region,start,end,name,score,strand) 
               set seq_region=replace(seq_region, 'chr', '');");
	  
	}

	$dbh->disconnect();


	#Finally set IMPORTED_ASSEMBLY status

	my $sql = "INSERT IGNORE into status_name(name) values('IMPORTED_${assembly}')";
	$dbh->do($sql);
	$sql =  "INSERT IGNORE into status(table_name, table_id, status_name_id) select '$table_name', 1, sn.status_name_id from status_name sn where sn.name='IMPORTED_${assembly}'";
	$dbh->do($sql);

	print ":: Finished loading:\t$file\n";


	#We should default to gzipping all files here
	#Unless we pass -no_ziop?

	if ($gzip && $link) {
	  
	  #We need to remove if the file and regenrate link if file was link?
	  
	    #if($link){
	    print ":: Restoring link\n";
	    system("cp -f ${file}.gz.backup ${file}.gz") == 0 or die('Failed to restore link');
	}
	elsif(! $no_compress){
	    print ":: Compressing file $file...\n";
	    system("gzip $file") == 0
		or die "Can't compress file $file: $!";
	}	
	#}
  }
}




sub write_bins () {

    my ($bin_start, $bin_end);

    for (my $i=0; $i<=$#bin; $i++) {

        $bin_start = $i*$bin_size+1;
        $bin_end = $bin_start+$bin_size-1;
        
        if (defined $bin[$i]) {
            printf OUT "%s\t%d\t%d\t.\t%.1f\n", 
            $seq, $bin_start, $bin_end, $bin[$i];
        }

    }

    return 1;

}

1;
