#!/usr/local/ensembl/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($file, $ofile, $ftype_name, $ctype_name, $pass, $line, $fset_id, %slice_cache, $status);
my ($fset, $dbhost, $exp_name, $dbname, $cdbname, $help, $pids, $man, @features, $chr, $not_status);
my $port = 3306;
my $anal_name = 'Nessie';
my $out_dir = "./";


GetOptions (
	    "feature_type=s"   => \$ftype_name,
	    "file=s"           => \$file,
	    "exp_name=s"       => \$exp_name,
	    "cell_type=s"      => \$ctype_name,
	    "feature_set_id=i" => \$fset_id,
	    "pass=s"           => \$pass,
	    "port=s"           => \$port,
            "dbname=s"         => \$dbname,
	    "dbhost=s"         => \$dbhost,
	    "probe_ids"        => \$pids,
            "cdbname=s"        => \$cdbname,
	    "analysis_name=s"  => \$anal_name,
	    "outdir=s"         => \$out_dir,
	    "help|?"           => \$help,
	    "man|m"            => \$man,
	    "chr=s"	       => \$chr,
	    "status=s"         => \$status,
	    "not_status"       => \$not_status,
	   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


### Set up adaptors and FeatureSet 

if(! $cdbname  || ! $dbname ){
  throw("You must provide a funcgen(-dbname) and a core(-cdbname) dbname");
}
throw("Must define your funcgen dbhost -dbhost") if ! $dbhost;
#throw("Must supply an input file with -file") if ! $file;
throw("Must supply a password for your fungen db") if ! $pass;


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					      -host => "ens-staging",
					      -dbname => $cdbname,
					      #-species => "homo_sapiens",
					      -user => "ensro",
					      -pass => "",
					      #	-group => 'funcgen',
					      -port => '3306',
					     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
						      -host => $dbhost,
						      -dbname => $dbname,
						      #-species => "homo_sapiens",
						      -user => "ensadmin",
						      -pass => $pass,
						      -dnadb => $cdb,
						      -port => '3306',
						     );

#should check db's here


my $rsa = $db->get_ResultSetAdaptor();
my $exa = $db->get_ExperimentAdaptor();
my $ex_obj = $exa->fetch_by_name($exp_name);#"Stunnenberg_all_OID_1963");
throw("Not a valid Experiment name in the DB:\t$exp_name") if ! $ex_obj;

my $analy_obj = $db->get_AnalysisAdaptor()->fetch_by_logic_name("VSN_GLOG");
my $slice = $cdb->get_SliceAdaptor()->fetch_by_region('chromosome', $chr); 
my $pfa = $db->get_ProbeFeatureAdaptor();
my @result_sets = @{$rsa->fetch_all_by_Experiment_Analysis($ex_obj, $analy_obj)};

foreach my $set (@result_sets) {
  
  if ($set->has_status($status)  && (! $not_status)){
        
    if($pids){
      &print_Probe_results($set, $status);
    }else{
      &print_ResultFeatures($set, $status);
    }
  }elsif((! $set->has_status($status)) && ($not_status)){
    if($pids){
      &print_Probe_results($set, "not_${status}");
    }else{
      &print_ResultFeatures($set, "not_${status}");
    }
  }
}

sub print_Probe_results{
  my ($set, $status) = @_;

  warn "Getting Probe results with status: $status\n";

  my $out_string = "";
  my $ofile = $out_dir.$exp_name."_".$set->dbID()."_${status}_ProbeFeature_results.${chr}";
  open (OUT, ">$ofile") || throw("Cannot open output:\t$ofile");
      

  my @features = @{$pfa->fetch_all_by_Slice_ExperimentalChips($slice, $set->get_ExperimentalChips())};
  warn "Feature array size: ". @features ."\n";
      
  foreach my $pfeature(@features){
    my $result = $pfeature->get_result_by_ResultSet($set);
    
    #this is restrieving results in a feature centric manner, so some may not have results
    if($result){
      $out_string .= "chr${chr}\t".$pfeature->probe_id(),"\t".$pfeature->start().
	"\t".$pfeature->end()."\t${result}\n";
    }

    #warn "chr${chr}\t".$pfeature->probe_id()."\t".$pfeature->start()."\t".$pfeature->end()."\t${result}\n";
  }

  print OUT $out_string;
      
  close(OUT);
  return;
}

sub print_ResultFeatures{
  my ($set, $status) = @_;

  warn "Getting ResultFeatures with status: $status\n";
  my $out_string = "";
  my $ofile = $out_dir.$exp_name."_".$set->dbID()."_${status}_ResultFeatures.${chr}";
  open (OUT, ">$ofile") || throw("Cannot open output:\t$ofile");
  #    warn "Displayable set ".$set->dbID()." has analysis ".$set->analysis->logic_name()."\n";
  
  my @features =  @{$set->get_ResultFeatures_by_Slice($slice)};
  warn "ResultFeature array size: ". @features ."\n";
    
  foreach my $feat (@features) {
     $out_string .= "chr${chr}". "\t". $feat->start() ."\t". $feat->end() ."\t". $feat->score() ."\n";
  }


  print OUT $out_string;
  close(OUT);
  return;
}


exit;

__END__

my $pfa = $db->get_PredictedFeatureAdaptor();
my $fset_adaptor = $db->get_FeatureSetAdaptor();

if($fset_id){
  $fset = $fset_adaptor->fetch_by_dbID($fset_id);


  throw("Could not retrieve FeatureSet with dbID $fset_id") if ! $fset;

  warn("You are loading PredictedFeature using a previously stored FeatureSet:\n".
       "\tCellType:\t".$fset->cell_type->name()."\n".
       "\tFeatureType:\t".$fset->feature_type->name()."\n");

  #should also check types and anal if the have been set
  #should add more on which experiment/s this is associated with and feature_set/dat_set_name when we have implemented it.
  #also ask continute question or use force_import flag
  #we could also do a count on the PFs in the set to make sure we're know we're adding to a populated set.
  #hard to repair if we do load ontop of another feature set


}elsif(! ($ftype_name && $ctype_name)){
  throw("Must provide a FeatureType and a CellType name to load your PredictedFeatures");
}else{
  my $anal =  $db->get_AnalysisAdaptor->fetch_by_logic_name($anal_name);
  my $ftype = $db->get_FeatureTypeAdaptor->fetch_by_name($ftype_name);
  my $ctype = $db->get_CellTypeAdaptor->fetch_by_name($ctype_name);


  throw("No valid CellType available for $ctype_name") if ! $ctype; 
  throw("No valid FeatureType available for $ftype_name") if ! $ftype;
  throw("No valid Analysis available for $anal_name") if ! $anal;
  
  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
    (
     -CELL_TYPE => $ctype,
     -FEATURE_TYPE => $ftype,
     -ANALYSIS => $anal,
    );

  ($fset) = @{$fset_adaptor->store($fset)};

  warn("Generated FeatureSet\n");
}


open (FILE, $file) || die "Unable to open input\t$file";

warn("Loading PredictedFeatures from $file\n");

@features =();
while ($line = <FILE>){

  chomp $line;
  my @tmp = split /\s+/, $line;
  my $start = $tmp[1];
  my $end = $tmp[2];
  my $score = $tmp[4];
  my $text = "enriched_site";
  my $chr = $tmp[0];

  $chr =~ s/chr//;

  #	print STDERR "$cdb->get_SliceAdaptor()->fetch_by_region(\'chromosome\', $chr);\n";
  

  if (! exists  $slice_cache{$chr}){
    $slice_cache{$chr} = $cdb->get_SliceAdaptor()->fetch_by_region('chromosome', $chr);
    throw("Could not generate slice for chromosome $chr") if ! $slice_cache{$chr};
  }

  my $pfeature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
    (
     -SLICE         => $slice_cache{$chr},
     -START         => $start,
     -END           => $end,
     -STRAND        => 1,
     -DISPLAY_LABEL => $text,
     -SCORE         => $score,
     -FEATURE_SET   => $fset,
    );
  
  push @features, $pfeature;

}

$pfa->store(@features);

warn("Loaded ".($.)." PredictedFeatures onto chromosomes ".(keys %slice_cache)."\n");



__END__
