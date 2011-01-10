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

ensembl-efg create_transcript_result_features.pl
  
=head1 SYNOPSIS

create_transcript_result_features.pl [options]

Options:

Mandatory
  -pass|p            The password for the target DB, if not defined in GroupDefs.pm
  -dbname|d          Defines the eFG dbname
  -user|u            The user name
  -species|s         Species name e.g. homo_sapiens
  -port|l            efg DB port number
  -no_dump           No text file dump flag
  -no_load           No Feature load flag
  -slice_name        Runs only on this slice e.g. chromosome:NCBI36:18:1:10000
  -quick             This employs a optimised but minimal query which will not return any probe name information.
                     A coverage value will be annotated instead.
  -overlap           % overlap require of probe on given feature level i.e. transcript or exon.
  -rollback          Forces rollback of existing FeatureSets
  -pol_crop          Polymerase 5' crop base pairs e.g. 2000 default = 5000 NOT YET IMPLEMENTED!
  -probe_length      Used to calculate overlap. Default = 50. NOTE, this is done generically for all probes, not dynamic!
  -opass
  -odbname
  -oport
  -ouser
Optional
  -assembly_version  Ensembl sssembly version of the experimental data you want to use e.g. 36 for homo_sapiens_core_48_36j, which is NCBI36
  -log_file          Defines the log file
  -help              Brief help message
  -man               Full documentation

=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-species|s>

Species name for the array.


=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut


use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;

use Getopt::Long;
use Pod::Usage;
use strict;

$| = 1;


### TO DO
#Minimum coverage? as a function of feature slice length?
#Full transcript type/5' crop config
#Get all Echips and fetch features only once for all rsets
#Then simply only get the result for each rset
#implement subs - see bottom of this script.
#Separate out DB

#my $reg = "Bio::EnsEMBL::Registry";
#Should we use registry here to get species alias
#Would have to delete and reset the efg DB and core DB from the registry
#This needs to be a method in the Funcgen::DBAdaptor

#Param definitions and defaults
my ($dbname, $pass, $port, $host, $species, $user, $clobber, $help, $man, $assembly_version, $exon_set);
my ($odbname, $opass, $oport, $ohost, $ouser, $no_dump, $no_load, $slice_name, @slices, $quick);
my $logic_name = 'VSN_GLOG';
my $outdir = './';
my $overlap = 80;
my $pol_crop = 5000;
my $probe_length = 50;
my $with_probe = 1;
my $min_coverage = 0;#Should we calculate this as a function of length

my %full_transcript_types = (
							#FeatureType->name => 5' crop?
							'H3K36me3' => undef,
							'PolII'    => undef, #5000,
						   );


GetOptions (
			"dbname|d=s"   => \$dbname,
			"pass|p=s"     => \$pass,
			"port|l=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"species|s=s"  => \$species,
			"assembly=s"   => \$assembly_version,
			"exon_set"     => \$exon_set,
			"no_dump"      => \$no_dump,
			"no_load"      => \$no_load,
			"quick"        => \$quick,
			"overlap=s"    => \$overlap,
			"slice_name=s" => \$slice_name,
			"min_coverage=s" => \$min_coverage,
			"pol_crop=i"   => \$pol_crop,
			"probe_length=i" => \$probe_length,
			"analysis=s"   => \$logic_name,
			"outdir|o=s"   => \$outdir,
			"rollback"     => \$clobber,
			"help|?"       => \$help,
			"man|m"        => \$man,
			"tee"          => \$main::_tee,
			"log_file=s"   => \$main::_log_file,
		   );

my @result_set_names = @ARGV;

die('Maybe you want to specify some ResultSet names to process?') if ! @result_set_names;

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

$with_probe = 0 if $quick;

#Mandatory param checking here

#This provides logging and rollback functionality
#All $main vars will take effect in Helper
my $helper =  Bio::EnsEMBL::Funcgen::Utils::Helper->new();


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => $host,
   -dbname => $dbname, #'HDACC_mus_musculus_funcgen_47_37a',#homo_sapiens_funcgen_49_36k',
   -species => $species,
   ##-dnadb => $cdb,
   -user => $user,
   -pass => $pass,
   -group => 'funcgen',
   -port => $port,
  );


$db->set_dnadb_by_assembly_version($assembly_version);

#check dbs
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;


##TO DO## 
#Need to set output DB here?
die('Does not yet accomodate writing to the "out DB"') if defined $odbname;

#We are currently arbitrarily setting the overlap by expanding the slice
#This only works if all probes are guranteed to be of $probe_length
#Implement dynamic overlap calc


if($no_load && $no_dump){
  die('You have specified -no_load and -no_dump, maybe you want this script to do something?');
}

die('Cannot specify an overlap percentage of greater than 100') if $overlap >100;

### What are we doing ? ###
#We want to use a result set to calculate averaged exon/transcript values, 
#where a probe contributes to an exon when it overlaps by 80%
#can be fairly low as we're dealing with 50mers

#We want 2 feature sets
#1. Retaining an exon context to highlight potential alternative transcripts
#2. Averaging over a transcript
#Resolve transcript to gene relationship in post processing...take both into ccount or just take longest?

#Each of the features has to have a link back to the exon/transcript
#Providng easy access to gene name
#DBEntry xref or simply use display label= GeneStableID:Exon/TranscriptStableID:ProbedbIDs
#We also want to capture which probes overlapped

#The we want to be able to dump this information in a nice tabular format.


#Let's go...grab all the adaptors
my $af_adaptor = $db->get_AnnotatedFeatureAdaptor();
my $pf_a = $db->get_ProbeFeatureAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $anal_a = $db->get_AnalysisAdaptor();
my $rset_a = $db->get_ResultSetAdaptor();
my $trans_adaptor = $db->dnadb->get_TranscriptAdaptor;
my $gene_adaptor = $db->dnadb->get_GeneAdaptor;
my $rset_analysis = $anal_a->fetch_by_logic_name('VSN_GLOG');
my $fset_analysis = $anal_a->fetch_by_logic_name('ProbeFeatureWindow');

if(! defined $fset_analysis && ! $no_load){

  #This could do with added parameter values
  #However, parameters will not get stored if analysis already present
  #Need separate paramters table?
  #Or just non-unique for logic_name and change implementation to fetch_all
  #The use set parameters
  #or fetch_by_logic_name_paramters?

  $fset_analysis = Bio::EnsEMBL::Analysis->new(
											   -logic_name      => 'ProbeFeatureWindow',
											   -db              => 'NULL',
											   -db_version      => 'NULL',
											   -db_file         => 'NULL',
											   -program         => 'NULL',
											   -program_version => 'NULL',
											   -program_file    => 'NULL',
											   -gff_source      => 'NULL',
											   -gff_feature     => 'NULL',
											   -module          => 'NULL',
											   -module_version  => 'NULL',
											   -parameters      => 'NULL',
											   -created         => 'NULL',
											   -description     => 'Average score of probes within a give feature window',
											   -display_label   => 'ProbeFeatureWindow',
											   -displayable     => 1,
											  );

  $anal_a->store($fset_analysis);
  $fset_analysis = $anal_a->fetch_by_logic_name('ProbeFeatureWindow');
}


#Get slices
if($slice_name){
  @slices = ($slice_a->fetch_by_name($slice_name));
}
else{
  #toplevel is all top level assembled regions i.e. chromosomes or an supercontigs which still haven't been assembled yet
  #1 arg is to include non-reference regions e.g. haplotypes
  @slices = @{$slice_a->fetch_all('toplevel', 1)};
}
#test slices here?

#This can give more than one result set if there are differing cell/feature_types with the same result_set name
my ($full_transcript_types, $exon_transcript_types, @rsets);

foreach my $rset_name(@result_set_names){
  my @tmp_rsets = @{$rset_a->fetch_all_by_name_Analysis($rset_name, $rset_analysis)};

  if(scalar(@tmp_rsets) == 1){
	push @rsets, $tmp_rsets[0];
  }
  elsif(scalar(@tmp_rsets) == 0){
	die("Could not find ResultSet named $rset_name");
  }
  else{
	die('Found '.scalar(@tmp_rsets).' Feature/CellType variants for '.$rset_name.
		  "\nThis script does not yet accomodate non unique ResultSet names");
  }

  #Set flag to build exon_transcript features
  #Saves a lot of Exon/ProbeFeature querying if we only have full_transcript_types e.g. H3K36me3

  my $ftype = $tmp_rsets[0]->feature_type->name;

  if(! $exon_transcript_types){
	$exon_transcript_types = 1 if ! grep(/^${ftype}$/, keys %full_transcript_types);
  }

  if(! $full_transcript_types){
	$full_transcript_types = 1 if grep(/^${ftype}$/, keys %full_transcript_types);
  }
}




#Set up/clobber old Data/FeatureSets
my %feature_sets;
my @feature_levels = ('transcript');
push @feature_levels, 'exon' if $exon_set;

foreach my $rset(@rsets){
  
  foreach my $level(@feature_levels){

	my $set_name = $level.'_features_'.$rset->name;

	if(! $no_load){
	  
	  my %set_params = (
						'-name'         => $set_name,
						'-feature_type' => $rset->feature_type,
						'-cell_type'    => $rset->cell_type,
						'-analysis'     => $fset_analysis,
						'-supporting_sets' => [$rset],
					   );
	  
	  my $dset = $helper->define_and_validate_sets($db, \%set_params, $clobber);
	  $feature_sets{$set_name}{'feature_set'} = $dset->product_FeatureSet;
	}

	#Create dump file handles
	if(! $no_dump){
	  $feature_sets{$set_name}{'file'} = open_file($outdir."/${set_name}.tab", '>');
	  #print header here?
	}


  }
}

$helper->log('');
$helper->log("Processing @feature_levels level features for the following ResultSets:\t@result_set_names");


#Define some vars and defaults 
my ($score, $total_probes, $display_label);
my $total_exon_cnt = 0;
my $total_transcript_cnt = 0;
my $expand_length = $probe_length * (1 - ($overlap/100)); #Translates to 80% overlap
my ($gene, $gid, $gname, @probe_features);

#Lets go!
foreach my $slice(@slices){
  $helper->log("Processing slice:\t".$slice->name);
  my $exon_cnt = 0;
  my $transcript_cnt = 0;
  my $chr = $slice->seq_region_name;

  my @trans = @{$trans_adaptor->fetch_all_by_Slice($slice)};
  $helper->log('Found '.scalar(@trans).' transcripts');

  foreach my $transcript(@trans){
	my $trans_id = $transcript->stable_id;
	my %result_cache;
	my ($textend_slice, $coverage);

	#warn ("Processing $trans_id");
	#Things slow down dramatically when we hit ENST00000399121 
	#Is this where the data starts?
	#Doesn't seem to be dumping for transcripts following this id...flush issue?

	if($full_transcript_types){
	  $textend_slice=$transcript->feature_Slice->extend($expand_length, $expand_length);
	}


	if(! $no_dump){
	  $gene = $gene_adaptor->fetch_by_transcript_id($transcript->dbID);
	  $gid = $gene->stable_id;
	  $gname = $gene->external_name;
	}

	
	#### Exon level caculations
	#Not here as we'll still be bringing back
	#all the features for the exon even if we only have full_transcript_type sets
	#We need to check the resultsets first to make sure we don't have any and set a flag

	if($exon_transcript_types){

	  foreach my $exon(@{$transcript->get_all_Exons}){
		my $exon_id = $exon->stable_id;

		#Could use RangeRegistry here?
		#Just use simple overlap rules for now
		
		
		#my $estart = $exon->seq_region_start;
		#my $eend   = $exon->seq_region_end;
		#Don't need this as we assuming all probes are 50mers, therefore all returned from
		#an exon slice expand by 80% of the probe length should meet our requirements
		#Will this also bring back probes which overhang end of extended slice?
		#If so then need to check start >=0 and end <= $exon->length
		
		my $expand_slice = $exon->feature_Slice->expand($expand_length, $expand_length);
		
		#warn "Processing $exon_id ".$expand_slice->length."\n";


		#Cache all overlapping probe ID and and result information


		#Can we bring back the features here given all the rsets
		#Then just get the result for each rset individually.
		#This will only speed up if using multiple sets

		#if(! $quick){
		#  @probe_features = @{$pf_a->fetch_all_by_Slice_ExperimentalChips($expand_slice, \@all_echips)});
		#}


		
		foreach my $rset(@rsets){

		  next if exists $full_transcript_types{$rset->feature_type->name};

		  $score = 0;
		  

		  #warn "processing exon probe features for rset ".$rset->name."\n";

		  #This is bogging down even when there are no probes
		  #Is this because this is the first long chromosome?
		  #Need to implement ProbeFeatures for all rsets once using all Echips
		  #Also need to implement new ResultFeature with nested Probe object
		  #Can we do something clever here and break up the Probe _obj_from_sth method
		  #such that we can reuse some of it in other adaptors e.g. ResultSetAdaptor
	
		  my $probe_annotation;
		  $coverage = 0;
		  @probe_features = @{$rset->get_ResultFeatures_by_Slice($expand_slice, undef, $with_probe)};
		  
		  next if scalar(@probe_features) == 0;#No features, next rset

		  foreach my $feature(@probe_features){

			#Filter probes which are not entirely contained within slice
			next if($feature->start < 0 || $feature->end > $expand_slice->length);
					
			if($quick){
			  $score += $feature->score;
			}
			else{
			  #Do we want probe_id here or probe_name?
			  #Using array to account for multi-matching probes
			  push @{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probes'}}, $feature->probe->get_probename;
			  #We need the probes array to build the exon transcript features
			  
			  #ResultFeature->probe testing
			  #$score += $feature->get_result_by_ResultSet($rset);
			  $score += $feature->score;

			  $display_label .= $feature->probe->get_probename.';';
			}

			$coverage++;
		  }
		  

		  #Finished all probes for this exon&result_set so create the exon features for each rset
		  if($coverage > 0){ 
			$result_cache{$rset->name}{$trans_id}{$exon_id}{'probes'} = $coverage if $quick;
			$result_cache{$rset->name}{$trans_id}{$exon_id}{'score'} = $score;
			$score = $score/$coverage;
			
			


			if(! $no_load){
			  my $exon_feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
				(
				 -slice => $slice,#This may need changing to top level???????
				 -start => $expand_slice->start,
				 -end   => $expand_slice->end,
				 #No strand information
				 -strand => 0,
				 -display_label => $trans_id.':'.$exon_id.':Coverage='.$coverage,#.':'.
#				 join(':', @{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probes'}}),
				 -score => $score,
				 -feature_set => $feature_sets{'exon_features_'.$rset->name}{'feature_set'},
				);
			
			  $af_adaptor->store($exon_feature);
			}
			
	
			#Dump text
			if(! $no_dump){
			  #Annoyingly print want use a hash element as a file handle
			  my $file_handle = $feature_sets{'exon_features_'.$rset->name}{'file'};

			  #GeneSID display name TranscriptSID ExonSIDchrom start end score list of probe name
			  print $file_handle join("\t", ($gid, $gname, $trans_id, $exon_id, $chr, $expand_slice->start, $expand_slice->end, $score, $coverage, $display_label))."\n";
			}
			
			$exon_cnt++;
		  }
		}
	  }
  
	  #### Finished all exons for this transcript ####
	  #So let's create the transcript features
	  foreach my $rset(@rsets){#Build (exon|full)transcript level info
		next if ! exists $result_cache{$rset->name}{$trans_id};

		my $ftype = $rset->feature_type->name;
		$score = 0;
		$total_probes = 0;
		$display_label = '';
		$coverage = 0;
		my ($start, $end);
		

		### Get the transcript level info ###
		if(exists $full_transcript_types{$ftype}){#full transcript coverage
		  my $ftextend_slice = $textend_slice;
		  
		  #redefine extend slice here dependant on PolII?
		  #Could defined by feature_type class ew POLYMERASE?
		  #Or have default extend lengths in full_transcript_types hash
		  $ftextend_slice = $ftextend_slice->extend(-1*$pol_crop) if $ftype eq 'PolII';
			
		  # Get probe/result features
		  @probe_features = @{$rset->get_ResultFeatures_by_Slice($ftextend_slice, undef, $with_probe)};

		  next if scalar(@probe_features) == 0;#No features, next rset
			
			
		  # Calculate start end score etc... 
			foreach my $feature(@probe_features){
			  
			  #Filter probes which are not entirely contained within slice
			  next if($feature->start < 0 || $feature->end> $ftextend_slice->length);
			  
			  if($quick){
				$coverage++;
				$score += $feature->score;
				#display label????????????????????????????????????????????????????????????????????????????
			  }
			  else{

				$display_label .= $feature->probe->get_probename.';';

				#ResultFeature->probe testing
				#$score += $feature->get_result_by_ResultSet($rset);
				$score += $feature->score;
			  }
			  $start = $ftextend_slice->start;
			  $end   = $ftextend_slice->end;
			
			}
		  }
		  else{#exon transcript coverage
	
			$start = ($transcript->seq_region_start - $expand_length);
			$end =  ($transcript->seq_region_end + $expand_length);
			
			my @exon_ids = sort keys %{$result_cache{$rset->name}{$trans_id}};
			next if scalar(@exon_ids) == 0;#next rset
		  
			foreach my $exon_id(@exon_ids){

			  my $exon_hash = $result_cache{$rset->name}{$trans_id}{$exon_id};
			  $score += $exon_hash->{'score'};

			  if($quick){
				$coverage += $exon_hash->{'probes'};
				#display label????????????????????????????????????????????????????????????????????????????
				#Is transcript level coverage enough, or do we want to dump exon level coverage?
			  }
			  else{
				$coverage += scalar(@{$exon_hash->{'probes'}});
				$display_label .= $exon_id.':'.join(',', @{$exon_hash->{'probes'}}).';';
			  }
			}
		  }
		

		  ### Now build/store/dump the transcript feature ###
		  if($coverage > 0){
			$score = $score/$coverage;
			
			#Store feature
			if(! $no_load){
			  my $transcript_feature =  Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
				(
				 -slice => $slice,#This may need changing to top level???????
				 -start => $start,
				 -end   => $end,
				 #No strand information
				 -strand => 0,
				 -display_label => $trans_id.':Coverage='.$coverage,#.':'.$display_label,
				 -score => $score,
				 -feature_set => $feature_sets{'transcript_features_'.$rset->name}{'feature_set'},
				);
			  
			  $af_adaptor->store($transcript_feature);
			}
			
			#Dump text
			if(! $no_dump){
			  #Annoyingly print want use a hash element as a file handle
			  my $file_handle = $feature_sets{'transcript_features_'.$rset->name}{'file'};
			  
			  #GeneSID display name TranscriptSID ExonSIDchrom start end score list of probe name
			  print $file_handle join("\t", ($gid, $gname, $trans_id, $chr, ($transcript->seq_region_start - $expand_length), 
											 ($transcript->seq_region_end + $expand_length), $score, $coverage, $display_label))."\n";
			  
			}
			
			$transcript_cnt++;
		  }
	  }
	}
  }


  #Counting and logging
  $total_exon_cnt += $exon_cnt;
  $total_transcript_cnt += $transcript_cnt;

  $helper->log('Finished processing '.$slice->name);
  $helper->log('Generated '.$transcript_cnt.' transcript level features for '.scalar(@rsets).' ResultSets');
  $helper->log('Generated '.$exon_cnt.' exon level features for '.scalar(@rsets).' ResultSets') if $exon_set;

}

$helper->log("Finished processing @feature_levels ProbeWindowFeatures for ResultSets:\t@result_set_names");
$helper->log("Total transcript features:\t".$total_transcript_cnt);
$helper->log("Total exon features:\t".$total_exon_cnt) if $exon_set;




#subs
#get_feature_stats_by_ResultFeatures
#get_feature_stats_by_ProbeFeatures
#store_and_dump_Feature
1;
