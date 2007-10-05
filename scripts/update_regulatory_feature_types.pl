#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file);
use Getopt::Long;


my $reg = "Bio::EnsEMBL::Registry";

my $species = 'homo_sapiens';
my ($pass, $test_slice, $schema_build, $type, $help, @slices);
my $port = $ENV{'PORT'};
my $user = $ENV{'WRITE_USER'};
my $host = $ENV{'HOST'};
my $out_dir = './';
my $log_file = $out_dir.'RegulatoryFeatures.type.log';

GetOptions( 'help'             => \$help,
            'host|h=s'         => \$host,
            'port=i'           => \$port,
            'user|u=s'         => \$user,
            'pass|p=s'         => \$pass,
			'slice=s'          => \$test_slice,
            'species|d=s'      => \$species,
			'schema_build|s=s' => \$schema_build,
		  );




my %regex_types =
  (
   #old definitions
   #'1...1.....' => 'Promoter associated',
   #DNase1','H3K27me3',
   
   #'1.0.001...' => 'Non-gene associated',
   #DNase1', 'H4K20me3', 'H3K36me3', 'H3K4me3', 'H3K79me3',
   
   #'11..01....' => 'Gene associated',
   #DNase1', 'CTCF', 'H3K36me3', 'H3K4me3', 
   
   
   #DNase1', 'CTCF', 'H4K20me3', 'H3K27me3', 
   #'H3K36me3', 'H3K4me3', 'H3K79me3', 'H3K9me3', 'TSS Proximal', 'TES Proximal'); 
   
   
   #new string
   #CD4_CTCF (31), CD4_DNASE_IMPORT (20), CD4_H2AZ (29), CD4_H2BK5me1 (30), CD4_H3K27me1 (32), CD4_H3K27me2 (33), CD4_H3K27me3 (34),CD4_H3K36me1 (35), CD4_H3K36me3 (51), CD4_H3K4me1 (36), CD4_H3K4me2 (37), CD4_H3K4me3 (49), CD4_H3K79me1 (38), CD4_H3K79me2 (39), CD4_H3K79me3 (40), CD4_H3K9me1 (41), CD4_H3K9me2 (42), CD4_H3K9me3 (43), CD4_H3R2me1 (44), CD4_H3R2me2 (45), CD4_H4K20me1 (46), CD4_H4K20me3 (47), CD4_H4R3me2 (50), CD4_PolII (48), GM06990_DNASE_IMPORT (21), Nessie_NG_STD_2_ctcf_ren_BR1 (22), Wiggle_H3K27me3 (13), Wiggle_H3K36me3 (14), Wiggle_H3K4me3 (15), Wiggle_H3K79me3 (16), Wiggle_H3K9me3 (17), Wiggle_H4K20me3 (12)
   
   #old defintions reconstituted on new string with cell specificity added
   
   #dnase indexs 1 & 24
   '.1..................0...1.00.1..' => 'Non-gene Associated',
   '.1..................0...0.00.1..' => 'Non-gene Associated - Cell type specific',
   '.0..................0...1.00.1..' => 'Non-gene Associated - Cell type specific',
   
   #gene
   '11......................1..01...' => 'Gene Associated',
   '11......................0..01...' => 'Gene Associated - Cell type specific',
   '10......................1..01...' => 'Gene Associated - Cell type specific',
   
   #promoter
   '.1......................1.1.....' => 'Promoter Associated',
   '.1......................0.1.....' => 'Promoter Associated - Cell type specific',
   '.0......................1.1.....' => 'Promoter Associated - Cell type specific',
);


my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $species.'_funcgen_'.$schema_build,
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														 );


die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_{$schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $efg_db;


my $ftype_adaptor = $efg_db->get_FeatureTypeAdaptor();
my %counts;
my %ftypes = (
			  'Unclassified' => undef,
			  'Unclassified - Cell type specific ' => undef,		
			  'Promoter Associated' => undef,
			  'Gene Associated' => undef,
			  'Non-gene Associated - Cell type specific' => undef,
			  'Promoter Associated - Cell type specific' => undef,
			  'Gene Associated - Cell type specific' => undef,
			  'Non-gene Associated - Cell type specific' => undef,
			 );


foreach $type(keys %ftypes){
  $counts{$type} = 0,

  my $ftype = $ftype_adaptor->fetch_by_name($type);

  throw('FeatureType $type not present in the DB, use import_type.pl before running update') if ! defined $ftype;

  $ftypes{$type} = $ftype;
}

my $handle = open_file($log_file, '>');

my $rf_fset = $efg_db->get_FeatureSetAdaptor()->fetch_by_name('RegulatoryFeatures');

#should we split this up into slices?
my $sql1 = 'UPDATE regulatory_feature set feature_type_id=';
my $sql2 = ' where regulatory_feature_id=';
my $cnt = 0;
my $total_cnt = 0;
my $warnings = 0;

if(defined $test_slice){
  @slices = ($efg_db->get_SliceAdaptor->fetch_by_name($test_slice));
}
else{
  @slices = @{$efg_db->get_SliceAdaptor->fetch_all('toplevel')};
}


print "Assigning new feature types:\n\t".join("\n\t", keys(%ftypes))."\n";

foreach my $slice(@slices){

  my @features = @{$rf_fset->get_Features_by_Slice($slice)};
  $total_cnt += scalar(@features);

  print "\nAssigning to ".scalar(@features)." RegulatoryFeatures on slice:\t".$slice->name."\n";

  #convert this to use sql directly?
  #or is this now redundant, do operation on flat files and reimport from file


  
 FEATURE: foreach my $rfeat(@features){
	my $assigned = 0;
	
	foreach my $regex(keys %regex_types){
	  
	  
	  #access directly just in case we put a method hack inplace to retain the binary string, but display something moire meaningful
	  if($rfeat->{'display_label'} =~ /$regex/){
		
		if($assigned && $rfeat->display_label ne $regex_types{$regex}){
		  print $handle 'WARNING: Skipping RegulatoryFeature with multiple type assignments:\tdbID:'.$rfeat->dbID."\t".$rfeat->feature_type->name().' & '.$regex_types{$regex}."\n";
		  $warnings ++;
		next FEATURE;
		}
		
		$rfeat->feature_type($ftypes{$regex_types{$regex}});
		$counts{$regex_types{$regex}} ++;
		$assigned = 1;
	  }
	}
	
	if(! $assigned){
	  
	  my @bits = split//, $rfeat->{'display_label'};
	  
	  $type = ($bits[1] != $bits[24]) ? 'Unclassified - Cell type specific' : 'Unclassified';
	  
	  $counts{$type} ++;
	  $rfeat->feature_type($ftypes{$type});
	}
		
	$cnt ++;
	$efg_db->dbc->db_handle->do($sql1.$rfeat->feature_type->dbID().$sql2.$rfeat->dbID());
  
	print "Updated $cnt RegulatoryFeatures\n" if(! ($cnt  % 10000));
  }
}

map {print 'Annotated '.$counts{$_}." $_ regulatory features\n"} keys %counts;
print "Finished updating a total of $cnt / $total_cnt RegulatoryFeature\n";
print "Found $warnings features which had mutliple assignments.  Please check $log_file\n" if $warnings;
