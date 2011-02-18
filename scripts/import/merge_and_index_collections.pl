#!/usr/bin/env perl

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

 merge_and_index_collections.pl
  
=head1 SYNOPSIS

 merge_and_index_collections.pl [options]

=head1 OPTIONS

 Mandatory params:

  --dbname          
  --dbuser          
  --dbpass          
  --dbport            
  --dbhost            
  --data_dir        The data dir of the input seq_region .col files

 Optional params:
  --result_set_name The name of the corresponding ResultSet. Default is --data_dir name.
  #--output_dir      
  --force           Forces over-writing of existing output files
  --packed_size     Default is 4bytes i.e. perl pack f(loat) encoding.


=head1 DESCRIPTION

B<This program> merges seq_region specific BLOB collection files into one indexed file.

=cut


# To do
# 1 Add target dir to copy to and updte dbfile_registry_dir?
# 2 Add some states to the DB?

use Getopt::Long;
use strict;
use warnings;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my ($dbuser, $dbpass, $dbname, $dbport, $dbhost);
my ($force, $data_dir, $rset_name, $packed_size);
my @tmp_args = @ARGV;

GetOptions (
			#DB connection
			"dbname=s" => \$dbname,
			"dbpass=s" => \$dbpass,
			"dbport=i" => \$dbport,
			"dbhost=s" => \$dbhost,
			"dbuser=s" => \$dbuser,

			
			'data_dir=s'        => \$data_dir, #with result_features.rset_name.window_size.sr_id.col files
			'result_set_name=s' => \$rset_name,
			'packed_size=i'     => \$packed_size,
			
			#Run modes
			'force'             => \$force,
		
		   ) or die("Specified invalid params: @tmp_args");

print "blob_test.pl @tmp_args\n";


if ( !( $dbuser && $dbname && $dbhost) ) {
   die('You must provide mysql -dbuser, -dbhost and -dbname arguments');
}


if(! -d $data_dir){
  die('You must provide a valid -data_dir argument which contains seq_region .col files');
}

#Defaults
$packed_size ||= 4;#perl pack f(loat) encoding

if(! $rset_name){
  ($rset_name = $data_dir) =~ s/\/$//;
  $rset_name               =~ s/.*\///g;
  print "Defaulting to ResultSet name:\t$rset_name\n";
}



#Set up DBAdaptor here
print "Setting up DBAdaptor:\t$dbuser\@$dbname:$dbhost:$dbport\n";
my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host   => $dbhost,
													  -dbname => $dbname,
													  -user   => $dbuser,
													  -pass   => $dbpass,
													  -port   => $dbport,
													  -dnadb_host => 'ensdb-archive',
													  -dnadb_user => 'ensro',
													  -dnadb_port => '5304',
);
#test connection
print "Testing DB connection...";
$db->dbc->db_handle;

my $slice_a = $db->get_SliceAdaptor;
my @slices    = @{$slice_a->fetch_all('toplevel', undef, undef, 1)};#inc dups


#Need to create registry of unrepresented slice names for each rset.
#And mandatory unrepresented slice names i.e. chroms!
my (%index, $fh, @rm_files);

#Need to get these from ResultFeatureAdaptor!
my $rf_adaptor   = $db->get_ResultFeatureAdaptor;
my $rset_adaptor = $db->get_ResultSetAdaptor;
my ($rset)       = @{$rset_adaptor->fetch_all_by_name($rset_name)};

if(! $rset){
  die('Could not find ResultSet $rset_name in the DB. Failed to acquire window_sizes');
}

if(! -d $data_dir){
  die('You need to specify a valid -data_dir');
}

$rf_adaptor->set_collection_defs_by_ResultSet($rset);
my @wsizes = @{$rf_adaptor->window_sizes};
my $first_wsize = 1;

foreach my $wsize(@wsizes){

  #Check for target file first
  my $merged_file =  "${data_dir}/".join('.', ('result_features',  $rset_name, $wsize, 'col'));

  if(-f $merged_file){
	if(! $force){
	  die("Please specify -force to over-write existing output file:\t$merged_file\n$!");
	}else{
	  print "Over-writing existing output file:\t$merged_file\n";
	}
  }

  my $total_length = 0;
  my $off_set      = 0;
  my $last_length  = 0;
  my $seen_y       = 0;
  my (@efg_sr_ids, @file_list);
   
  for my $slice(@slices){

	my $sr_name    = $slice->seq_region_name;
	my $slice_file = $data_dir.'/'.join('.', ('result_features',  $rset_name, $wsize, $sr_name, 'col'));

	if(! -f $slice_file){

	  if(($sr_name =~ /^[0-9]+$/) ||
		 ($sr_name eq 'X') ||
		 ($sr_name eq 'Y')) {
		die("Found absent col file for mandatory seq_region:\t$sr_name");
	  }
	  next;
	}


	my $efg_sr_id = $rf_adaptor->get_seq_region_id_by_Slice($slice);
	push @efg_sr_ids, $efg_sr_id;
	push @file_list, $slice_file; 

	#Set new byte off set for this slice
	$off_set           = $off_set + ($last_length);
	$index{$efg_sr_id} = $off_set;

	#Calculate last length
	$last_length     = $slice->length / $wsize;#rather than end, just in case we don't start at 1
	my $tmp_last_length = int($last_length);
	$last_length     = $tmp_last_length + 1 if ($last_length != $tmp_last_length);
	$last_length    *= $packed_size;

	$total_length    = $off_set + $last_length;

	#Now validate file length
	#This should have already been done in the Collector
	$rf_adaptor->validate_file_length($slice_file, $last_length, 1);#1 is binmode;
  }




  #Could probably do this via th dump by ordering by seq_region_name/id
  #But mysql sort is not reliably the same as unix/perl sort
  #So let's do it here to be sure.
  
  print "Total $wsize window_size data length:\t\t$total_length\n";

  #Encode releae version in the file name for visibility.
  #How big does the index need to be? I (32 bit unsigned) endian order depends on arch
  #v for size of index
  #v for index key i.e. sr_id
  #V for offset (4bytes?, max unsigned is 4,294,967,295)
  #Max offset is??
  #=>Nchrs * (max chr length) / (min window size) * packed_size
  #=> ~23(ignore non-chr for now) * 249250621 / 30 * 4 (also ignore index length)
  #=>764368568 < 4294967295
  #              404840556
  #Plenty of room for ~ 6* more data!
  
  #Encode index as key(sr_id - v 2 bytes) value(offset - V 4 bytes) pairs

  #V (long) is 4 bytes, which is actually an int? long is 8 bytes according to config!
  #
  
  
  
  
  my $num_keys      = scalar(keys %index);
  my $index_size    = scalar(keys %index)*(6);
  #Full index is index_size + key-values pairs
  #Could add file length here to enable validation
  #i.e. seek to full length and make sure we don't get and EOF until we seek one byte more
  
  #This is not right??
  #Seems to be 1/2 the size.
  #perl -V:{short,int,long{,long}}size
  #shortsize='2';
  #intsize='4';
  #longsize='8';
  #longlongsize='8';
  #Actual size can change due to compiler
  #	   use Config;
  #       print $Config{shortsize},    "\n";
  #       print $Config{intsize},      "\n";
  #       print $Config{longsize},     "\n";
  #       print $Config{longlongsize}, "\n";

  
  #Does this also mean that changing data between perl builds 32 -> 64 will screw the reads?
  #Or do we just have to change the unpack format?
  #Let's just test that they match the compiling platform?
  #This has to be done for each collection?

  
  #It seems like the size on disk are 1/2 of those listed by perl!
  #So the first read was offsetting the next read by 1/2 of the original read length!
  
  
  #No I just got the size wrong here as we're using f not v!!!!!
	
  
  
  my $total_index_size = 2 + ($num_keys * (2 + 4));
  
  #Adjust the offsets to account for the total index size
  foreach my $key(keys %index){
	#warn "Adjusting $key index value from $index{$key}\t";
	$index{$key} += $total_index_size;
	#warn "to ".$index{$key};
  }


  my $pack_template = 'v(vV)'.$num_keys;
  #warn "pack template is $pack_template";
  my @index = %index;
  my $packed_index  = pack($pack_template, ($index_size, @index));
  
    
  #Let's just validate the index here
  #my @unpacked_index = unpack($pack_template, $packed_index);
  #my $unpacked_index_size = shift(@unpacked_index);
  #warn "Index size($index_size) unpacked is $unpacked_index_size";
  #my %unpacked_index = @unpacked_index;
  #warn "Keys ".$num_keys.' unpacked is '.scalar(keys(%unpacked_index));
  #for my $key(keys(%unpacked_index)){
  #  warn "key $key value ".$index{$key}.' unpacked is '.$unpacked_index{$key};
  #}


  #Dumps the index to file and add it to the start of the cat list
  my $index_file = "${data_dir}/".join('.', ('result_features',  $rset_name, $wsize, 'idx'));  
  print "Writing index:\t\t\t$index_file\n";
  open($fh, '>', $index_file) or die("Cannot open $index_file");
  binmode($fh);
  print $fh $packed_index;
  close($fh);


  #Validate index file
  $fh        = $rf_adaptor->get_filehandle($index_file, {-binmode => 1});
  #Do this directly here as we have little need for grabbing the whole index normally
  my $index_ref = $rf_adaptor->{file_cache}{$index_file}{off_sets};

  #### Now validate written index

  #Test num keys is same
  if($num_keys != scalar(keys %{$index_ref})){
	die("Original index has $num_keys whilst read and unpacked index has ".scalar(keys %{$index_ref}));
  }
  
  for my $key(keys(%{$index_ref})){ #Test values match
	
	if($index{$key} != $index_ref->{$key}){
	  die("Original key $key value ".$index{$key}.
		  ' does not match read and unpacked value '.$index_ref->{$key});
	}
  }
  

  #Generate merged file
  my $cat = "cat $index_file @file_list > $merged_file";
  print "Generating indexed file:\t$merged_file\n";
  system($cat) == 0 or die ("Cannot cat files:\n$cat");

  #Validate length
  $total_length += $total_index_size;
  $rf_adaptor->validate_file_length($merged_file, $total_length, 1);#1 is binmode;

  push @rm_files, ($index_file, @file_list);
}

		
#Remove index and slice files
print "Removing index and sr_name col files\n";
system("rm -f @rm_files") == 0 or 
  warn("Failed to remove tmp files:\nrm -f @rm_files\n");


#Should we set some states here?

1;
