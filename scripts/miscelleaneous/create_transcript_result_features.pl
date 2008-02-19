
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
  -overlap           % overlap require of probe on given feature level i.e. transcript or exon.
  -rollback          Forces rollback of existing FeatureSets
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

#my $reg = "Bio::EnsEMBL::Registry";
#Should we use registry here to get species alias
#Would have to delete and reset the efg DB and core DB from the registry
#This needs to be a method in the Funcgen::DBAdaptor

#Param definitions and defaults
my ($dbname, $pass, $port, $host, $species, $user, $clobber, $help, $man, $assembly_verion, $exon_set);
my ($odbname, $opass, $oport, $ohost, $ouser, $no_dump, $no_load);
my $logic_name = 'VSN_GLOG';
my $outdir = './';
my $overlap = 80;

GetOptions (
			"dbname|d=s"     => \$dbname,
			"pass|p=s"     => \$pass,
			"port|l=s"     => \$port,
			"host|h=s"     => \$host,
			"species|s=s"  => \$species,
			"assembly"     => $assembly_version,
			"exon_set"     => \$exon_set,
			"no_dump"      => \$no_dump,
			"no_load"      => \$no_load,
			"overlap=s"    => \$overlap,
			"analysis=s"   => \$logic_name,
			"outdir|o=s"   => \$outdir,
			"rollback"      => \$clobber,
			"help|?"       => \$help,
			"man|m"        => \$man,
			"tee"          => \$main::_tee,
			"log_file=s"   => \$main::_log_file,
		   );

my @rset_names = @ARGV;

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#Mandatory param checking here

#This provides logging and rollback functionality
#All $main vars will take effect in Helper
my $helper =  Bio::EnsEMBL::Funcgen::Utils::Helper->new();


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => $host,
   -dbname => $dbname, #'HDACC_mus_musculus_funcgen_47_37a',#homo_sapiens_funcgen_49_36k',
   #-species => "Homo_sapiens",
   -species => 'mus_musculus',
   #-dnadb => $cdb,
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
throw('Does not yet accomodate writing to the "out DB"') if defined $odbname;

warn("NEED TO CHECK OVERLAP/EXPAND RULE WORKS!");


if($no_load && $no_dump){
  throw('You have specified -no_load and -no_dump, maybe you want this script to do something?');
}

throw('Cannot specify an overlap percentage of greater than 100') if $overlap >100;

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
my $af_a = $db->get_AnnotatedFeatureAdaptor();
my $pf_a = $db->get_ProbeFeatureAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $anal_a = $db->get_AnalysisAdaptor();
my $rset_a = $db->get_ResultSetAdaptor();
my $trans_adaptor = $db->dnadb->get_TranscriptAdaptor;
my $gene_adaptor = $db->dnadb->get_GeneAdaptor;

my @result_set_names = ();
my $rset_analysis = $anal_a->fetch_by_logic_name('VSN_GLOG');
my $fset_analysis = $anal_a->fetch_by_logic_name('ProbeFeatureWindow');

if(! defined $fset_analysis){

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

  $anal_a->store($fset_anal);
  $fset_analysis = $anal_a->fetch_by_logic_name('ProbeFeatureWindow');
}





}
#toplevel is all top level assembled regions i.e. chromosomes or an supercontigs which still haven't been assembled yet
#1 arg is to include non-reference regions e.g. haplotypes
my @slices = @{$slice_a->fetch_all('toplevel', 1)};

#This can give more than one result set if there are differing cell/feature_types with the same result_set name
my @rsets = @{map $rset_a->fetch_all_by_name_analysis($_, $rset_analysis), @result_set_names)};


#Set up/clobber old Data/FeatureSets
my %feature_sets;
my @feature_level = ('transcript');
push @feature_level, 'exon' if $exon_set;

foreach my $rset(@rsets){

  foreach my $level(@feature_levels){

	my $set_name = $level.'_features_'.$rset->name;

	if(! $no_load){
	  
	  my %set_params = (
						'-name'         => $set_name,
						'-feature_type' => $rset->feature_type,
						'-cell_type'    => $rset->cell_type,
						'-analysis'     => $fset_anlaysis,
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

#Define some vars and defaults 
my ($score, $total_probes, $display_label);
my $total_exon_cnt = 0;
my $total_transcript_cnt = 0;
my $probe_length = 50;#Should we defien the dynamically for each probe? WOuld impact on overlap logic.
my $expand_length = $probe_length * (1 - ($overlap/100)); #Translates to 80% overlap
my ($gene, $gid, $gname);

#Lets go!
foreach my $slice(@slices){
  $helper->log("Processing slice:\t".$slice->name);
  my $exon_cnt = 0;
  my $transcript_cnt = 0;
  my $chr = $slice->seq_region_name;

  foreach my $transcript(@{$trans_adaptor->fetch_all_by_Slice($slice)}){
	my $trans_id = $transcript->stable_id;
	my %result_cache;

	if(! $no_dump){
	  $gene = $gene_adaptor->fetch_by_transcript_id($transcript->dbID);
	  $gid = $gene->stable_id;
	  $gname = $gene_display_id;
	}

	#This is only used for Xn sets
	#How are -ve strand transcripts caught by +ve strand design tiling array?
	#my $strand = $transcript->strand;
	#Don't need this as Xn protocol generates double stranded cDNA which is hybed

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

	  my $overlap_slice = $exon->feature_Slice->expand($expand_length, $expand_length);

	  #Cache all overlapping probe ID and and result information
	  
	  foreach my $rset(@rsets){
		$score = 0;

		foreach my $probe_feature(@{$pf_a->fetch_all_by_Slice_ExperimentalChips
									  ($overlap_slice, $rset->get_ExperimentalChips)}){
		  

		  #Do we want probe_id here or probe_name?
		  #Using array to account for multi-matching probes
		  push @{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probenames'}}, $probe_feature->probe->get_probename;
		  
		  $score += $probe_feature->get_result_by_ResultSet($rset);
		}

		#Finished all probes for this exon&result_set so create the exon features for each rset
		$result_cache{$rset->name}{$trans_id}{$exon_id}{'score'} = $score;
		$score = $score/scalar(@{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probenames'}});
		
		if(! $no_load){


		  my $exon_feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
			(
			 -slice => $slice,#This may need changing to top level???????
			 -start => $expand_slice->seq_region_start,
			 -end   => $expand_slice->seq_region_end,
			 #No strand information
			 -strand => 0,
			 -display_label => $trans_id.':'.$exon_id.':'.
			 join(':', @{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probenames'}}),
			 -score => $score,
			 -feature_set => $feature_sets{'exon_'.$rset->name}{'feature_set'},
			);
		  
		  $af_adaptor->store($exon_feature);
		}

	
		#Dump text
		if(! $no_dump){
		  #GeneSID display name TranscriptSID ExonSIDchrom start end score list of probe name
		  print $feature_sets{'exon_'.$rset->name}{'file'} join("\t", ($gid, $gname, $trans_id, $exon_id, $chr, $expand_slice->seq_region_start, $expand_slice->seq_region_end, $score, join(';', @{$result_cache{$rset->name}{$trans_id}{$exon_id}{'probenames'}})))."\n";
		}

		$exon_cnt++;
	  }
	}
  
	#Finished all exons for this transcript
	#So let's create the features
	$score = 0;
	$total_probes = 0;
	$display_label = '';

	#Build transcript level info
	foreach my $exon_id(sort keys %{$result_cache{$rset->name}{$trans_id}}){
	  my $exon_hash = $result_cache{$rset->name}{$trans_id}{$exon_id};
	  $score += $exon_hash->{'score'};
	  $total_probes += scalar(@{$exon_hash->{'probenames'}});
	  $display_label .= $exon_id.':'.join(',', @{$exon_hash->{'probenames'}}).';';
	}

	#Store feature
	if(! $no_load){
	  
	  my $transcript_feature =  Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
		(
		 -slice => $slice,#This may need changing to top level???????
		 -start => ($transcript->seq_region_start - $expand_length),
		 -end   => ($transcript->seq_region_end + $expand_length),
		 #No strand information
		 -strand => 0,
		 -display_label => $trans_id.':'.$display_label,
		 -score => $score,
		 -feature_set => $feature_sets{'transcript_'.$rset->name}{'feature_set'},
		);
	  
	  $af_adaptor->store($transcript_feature);
	}

	#Dump text
	if(! $no_dump){
	  #GeneSID display name TranscriptSID ExonSIDchrom start end score list of probe name
	  print $feature_sets{'transcript_'.$rset->name}{'file'} join("\t", ($gid, $gname, $trans_id, $chr, $expand_slice->seq_region_start, $expand_slice->seq_region_end, $score, $display_label))."\n";

	}

	$transcript_cnt++;
  }


  #Counting and logging
  $total_exon_cnt += $exon_cnt;
  $total_transcript_cnt += $transcript_cnt;

  $helper->log('Finished processing '.$slice->name);
  $helper->log('Generated '.$transcript_cnt.' transcript level features for '.scalar(@rsets).' ResultSets');
  $helper->log('Generated '.$exon_cnt.' exon level features for '.scalar(@rsets).' ResultSets') if $exon_Set;

}

$helper->log("Finished processing @feature_level ProbeWindowFeatures for ResultSets:\t@result_set_names");
$helper->log("Total transcript features:\t".$total_transcript_cnt);
$helper->log("Total exon features:\t".$total_exon_cnt) if $exon_set;




1;
