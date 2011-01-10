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

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($ftype_name, $ctype_name, $pnames, $pass, $status, $chr);
my ($fset, $exp_name, $db_name, $cdb_name, $help, $man, @features, $slices, $not_status);
my $port = $ENV{'EFG_PORT'};
my $host = $ENV{'EFG_HOST'};
my $user = $ENV{'EFG_READ_USER'};
my $cdb_port = $ENV{'CORE_PORT'};
my $cdb_host = $ENV{'CORE_HOST'};
my $cdb_user = $ENV{'CORE_USER'};

my $anal_name = 'VSN_GLOG';
my $probe_names = 0;
my $out_dir = "./";
my $format = 'ProbeFeature';





GetOptions (
			#"feature_type=s"   => \$ftype_name,
			#"file=s"           => \$file,
			"exp_name=s"       => \$exp_name,
			#"cell_type=s"      => \$ctype_name,
			#"feature_set_id=i" => \$fset_id,
			"pass=s"           => \$pass,
			"port=s"           => \$port,
            "dbname=s"         => \$db_name,
			"host=s"           => \$host,
			"format=s"         => \$format,
			"cdbname=s"       => \$cdb_name,
			"cdbhost=s"       => \$cdb_host,
			"cdbuser=s"       => \$cdb_user,
			'probe_names'      => \$probe_names,
			"analysis=s"       => \$anal_name,
			"outdir=s"         => \$out_dir,
			"help|?"           => \$help,
			"man|m"            => \$man,
			"slice=s"          => \$slices,
			"chr=s"	           => \$chr,
			"status=s"         => \$status,
			"not_status"       => \$not_status,
		   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

### Check conflicting params
$format = lc($format);

if($format eq 'regulatoryfeature' && defined $exp_name){
  #Does throw exit from a script?
  throw("Conflicting parameters:\n\t-format\t$format\n\t-exp_name\t$exp_name");
  die("Throw not exiting");
}

die("Conflicting parameters:\n\t-slice\t$slice\n\t-chr\t$chr") if($chr && $slices);

#Mandatory params

if(! $cdb_name  || ! $db_name ){
  throw("You must provide a funcgen(-dbname) and a core(-cdbname) dbname");
}
throw("Must define your funcgen dbhost -dbhost") if ! $host;
#throw("Must supply an input file with -file") if ! $file;
#throw("Must supply a password for your fungen db") if ! $pass;

#Keep this as local dbname not likely to work with auto dnadb

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host   => $cdb_host,
   -dbname => $cdb_name,
   -user   => $cdb_user,
  # -pass   => "",
   -group => 'core',
   -port => $cdb_port,
  );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host   => $host,
   -dbname => $db_name,
   -user   => $user,
   -group  => 'funcgen',
   -pass   => $pass,
   -dnadb  => $cdb,
   -port   => $port,
  );

#should check db's here
#Check db_handle|dbc?



#Validate formats and set adaptor here
my $feature_adaptor;


if($format eq 'probefeature'){
  $feature_adaptor = $db->get_ProbeFeatureAdaptor();
}
elsif($format eq 'resultfeature'){
  $feature_adaptor = $db->get_ResultSetAdaptor();
}
elsif($format eq 'regulatoryfeature'){
  $feature_adaptor = $db->get_RegulatoryFeatureAdaptor();
}



### Set up adaptors and FeatureSet 

my $rsa = $db->get_ResultSetAdaptor();
my $exa = $db->get_ExperimentAdaptor();
my $ex_obj = $exa->fetch_by_name($exp_name);#"Stunnenberg_all_OID_1963");
throw("Not a valid Experiment name in the DB:\t$exp_name") if ! $ex_obj;

my $analysis = $db->get_AnalysisAdaptor()->fetch_by_logic_name($anal_name) 
  || die ("Not a valid analysis logic_name:\t$anal_name");

if(defined $chr){
  $slices = [ $cdb->get_SliceAdaptor->fetch_by_region('chromosome', $chr) ];
}
elsif(defined $slices){
  $slices = [ $cdb->get_SliceAdaptor->fetch_by_name($slices) ];
}
else{
  print "Defaulting to use all toplevel slices\n";
  $slices = $cdb->get_SliceAdaptor->fetch_all('toplevel');
}


die("You have no valid slices.  Did you specify a valid chromosome?") if(! @$slices);


my @result_sets = @{$rsa->fetch_all_by_Experiment_Analysis($ex_obj, $analysis)};


foreach my $set (@result_sets) {
  
  if(! defined $status){
	&print_data($set);
  }
  elsif ($set->has_status($status)  && (! $not_status)){
	&print_data($set, $status);
  }
  elsif((! $set->has_status($status)) && ($not_status)){
	&print_data($set, 'not_'.$status);
  }
}




sub print_data{
  my ($set, $status) = @_;

  my $logic_name  = $set->analysis->logic_name;


  my $status_string = ($status) ? " with status '$status'" : '';
  $status = ($status) ? "_${status}" : '';

  foreach my $slice(@$slices){
	my ($result, @features);
	my $no_score = 0;
	my $probe_id_name = '';

	print "Getting ".$set->name." $format results${status_string} for slice:\t".$slice->name."\n";

	my $out_string = "";
	
	my $ofile = $out_dir.$set->name()."${status}_${format}.".$logic_name.'.'.$slice->name;
	open (OUT, ">$ofile") || throw("Cannot open output:\t$ofile");


	if($format eq 'ResultFeature'){
	  @features =  @{$set->get_ResultFeatures_by_Slice($slice)};#Change this to ResultSetAdaptor
	}
	elsif($format eq 'ProbeFeature'){
	  @features = @{$pfa->fetch_all_by_Slice_ExperimentalChips($slice, $set->get_ExperimentalChips())};
	}
	elsif($format eq 'RegulatoryFeature'){

	}

	print $set->name." has ". @features ." ${format}s for slice:\t".$slice->name."\n";
    
	foreach my $feat (@features) {
	  
	  if($format eq 'ResultFeature'){
		$result = $feat->score;
	  }
	  elsif($format eq 'ProbeFeature'){
		$result = $feat->get_result_by_ResultSet($set);
		$probe_id_name = (($probe_names) ? $feat->probe->get_probename() : $feat->probe_id())."\t";
	  }


	  if($result){
		$out_string .= $slice->seq_region_name()."\t".$probe_id_name;
		$out_string .= join("\t", ($feat->start(), $feat->end(), $result))."\n";
	  }
	  else{
		$no_score++;
	  }
	}

	if($no_score){
	  print "Found $no_score features with no result associated for slice:\t".$slice->name.
		"\nThis is probably due to importing a subset of an array and remapping the probes\n";
	}
	
	print OUT $out_string;
	close(OUT);
  }

  return;
}



1;
