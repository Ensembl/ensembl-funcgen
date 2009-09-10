#!/software/bin/perl


=head1 NAME

load_bed_source.pl

=head1 SYNOPSIS

 LoadBedSources|load_bed_source.pl [options]

 e.g. Using the associated efg environment function

 LoadBedDASSources  -host $DB_HOST --user $DB_USER --pass WRITE_PASS --dbname efg_DAS --files ES_DNase_le1m_reads.bed.gz --names ES_DNase --profile --frag_length 150 --bin_size 25 -reads


=head1 DESCRIPTION

This script loads a DAS source from raw read alignments(bed) and/or 
profiles generated from them.


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
  --profile_input Skips profile generation, expects input(--files) is profile

 Input
  --files         Space separate list of input files
  --names         Optional space separated list of source names, must preferably end in valid 
                  feature type to enable automatic colour usage e.g. ES_DNase1
  --prefix        Optional prefix for source names
  
 Profile Options, only required if --profile specified
  --bin_size      Bin size to compute profile scores over
  --frag_length   Length to extend single reads to
  Need to add paired end support here
		
 Other
  --help          Prints a short help message and exits
  --man|m         Prints the man page


=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see
mysqlro
    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


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



use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use DBI;

my ($pass,$host,$user,$dbname,$prefix, $input_file, @files, @names);
my ($no_load, $bin_size, $frag_length, $profile_input, $source_name);
my ($profile, $reads, %output_files);
my $port = 3306;



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
            'frag_length=i'    => \$frag_length,
            'help|?'           => sub { pos2usage(-exitval => 0, -message => $params_msg);},
            'man|m'            => sub { pos2usage(-exitval => 0, -message => $params_msg, verbose => 2);},
		   ) or pod2usage ( -exitval => 1,
							-message => $params_msg
						  );



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


my (@bin, $start_bin, $start_bin_start, $end_bin, $end_bin_start,
	$seq, $read_start, $read_end, $read_length, $ori, $read_extend);

if($profile && ! $profile_input){


  #Build profile
  print ":: Building profile for:\t$input_file\n";
  warn "WARNING:\tThis does not yet support paired end data\n";
  
  open(CMD, "file -L $input_file |")
	or die "Can't execute command: $!";
  my $gzip = grep {/gzip compressed data/} (<CMD>);
  close CMD;

  if($gzip){
	open(FILE, "gzip -dc $input_file |") or die ("Can't open compressed file:\t$input_file");
  }
  else{
	open(FILE, $input_file) or die ("Cannot open file:$input_file");
  }


  my $binsize = sprintf("%03d", $bin_size);

  #We need to validate this input file name
  #Let's not depend on _reads.bed

  my $output_file = $input_file."_profile_${binsize}";
  
  #Need to test if compressed here!
  

  open(OUT, "| gzip -c > $output_file")
    or throw ("Can't open out file $output_file");

  while (<FILE>) {
    chomp;
    my @col = split("\t");
  
    if (defined $seq && $seq ne $col[0]) {
	  &write_bins();
	  @bin = ();
	}

    $seq = $col[0];
    $read_start = $col[1];
    $read_end = $col[2];
    $read_length = $read_end-$read_start+1;
    $ori = $col[5];

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

  close FILE;
  close OUT;


  
  $output_files{profile} = $output_file;
  #push @files, $out;
  #Okay this is not working with the names array!
  
  


}


if( ! $no_load){

  #warn("No Hydra source name prefix specified!\n") if (! $prefix);

  my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;port=$port",
						 "$user", "$pass",
						 {RaiseError => 1,
						  mysql_auto_reconnect => 1});
  





  #Validate/identify file type, not by name!
  foreach my $type(keys %output_files){
	my $table_name = $source_name;
	$table_name = "${prefix}_${table_name}" if $prefix;
	$table_name = "bed_${type}_${table_name}";

	my $file = $output_files{$type};
	
	print ":: Loading $type file:\t$file\n";
	open(CMD, "file -L $file |")
	  or die "Can't execute command: $!";
	my $gzip = grep {/gzip compressed data/} (<CMD>);
	close CMD;

	my $link;
	$link = 1 if -l $file;
	
	if ($gzip) {
	  print ":: Decompressing $file\n";

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

	print ':: Table name: ', $table_name, "\n";

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
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;");#why do we need this COLLATE latin1_bin here?

  
  #This needs to change to mysqlimport!!!
  
  if ($type eq 'reads') {
	
	#what is @mm for?
	
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

	print ":: Finished loading $file\n";

	if ($gzip) {
	  
	  #We need to remove if the file and regenrate link if file was link?
	  
	  if($link){
		print ":: Restoring link\n";
		system("cp -f ${file}.gz.backup ${file}.gz") == 0 or die('Failed to restore link');
	  }
	  else{
		print ":: Compressing file $file...\n";
		system("gzip $file") == 0
		  or die "Can't compress file $file: $!";
	  }	
	}
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
