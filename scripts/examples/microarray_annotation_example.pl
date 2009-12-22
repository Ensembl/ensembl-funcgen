#!/software/bin/perl -w

#
# Ensembl example script for microarray annotation access
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

micro_array_annotation_example.pl

=head1 DESCRIPTION

This script extracts probes mapping and probeset annotations for the Mouse
Affy array MOE430A, on the chromosome 2:26771920:26789448


=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

# Remember to add the ensembl and ensembl-functgenomics packages to your PERL5LIB


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

####################################
# Use the Registry
# Load database details
# Then get the slice/array adaptors
####################################

my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(
								 -host => 'ens-staging2',#'ensembldb.ensembl.org',
								 -user => 'ensro',#'anonymous',
#								 -db_version => 56,
				 );

#print "Using api: ",$registry->software_version,"\n";

my $slice_adaptor = $registry->get_adaptor("mus musculus","core","Slice");
my $pfa           = $registry->get_adaptor("mus musculus","funcgen","ProbeFeature");
my $aa            = $registry->get_adaptor("mus musculus","funcgen","Array");
my $pba           = $registry->get_adaptor("mus musculus","funcgen","ProbeSet");
my $tx_adaptor    = $registry->get_adaptor("mus musculus","core","Transcript");
my $gene_adaptor  = $registry->get_adaptor("mus musculus","core","Gene");

#####################################################
# Get the slice from chr 20, start 400000, end 500000
# Get the oligo array object for MOE430A
#####################################################

my $slice = $slice_adaptor->fetch_by_region('chromosome','2',26771920, 26789448);
my $array = $aa->fetch_by_name_vendor('MOE430A', 'AFFY');


#########################################
# Get mapped probes for that slice
#########################################

my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n\n";

foreach my $pf(@$oligo_features){
    my $probeset   = $pf->probeset->name;
    my $probename  = $pf->probe->get_probename($array->name);
	my $slice_name = $pf->feature_Slice->name;
	#extended cigar_string as defined by SAMTools group
    print "\tProbe $probeset $probename is aligned on $slice_name with cigar string ".$pf->cigar_string."\n";
  }


############################################
# Get all ProbeFeatures annotated on the probesets annotated to Exons
# Note: ProbeFeature level DBEntries are used to build Probe/ProbeSet
# level transcript annotatations.
my $chro2 = $slice_adaptor->fetch_by_region('chromosome','2');
print "\nProbeFeatures xref'd to Transcripts on slice ".$chro2->name."\n";
my @genes = @{ $slice->get_all_Genes() };

foreach my $gene  (@genes) {
  print "\n",$gene->external_name," ",$gene->display_id(),"\n";
  
  foreach my $transcript(@{$gene->get_all_Transcripts}){ 
	print  "  ",$transcript->external_name," ",$transcript->display_id()," has the following probes matches:\n";
	
	#Grab all the ProbeFeatures which are xref'd to this transcript
	foreach my $pf(@{$pfa->fetch_all_by_external_name($transcript->stable_id)}){

	  #Now filter ProbeFeatures which are on our array
	  my $on_our_array = 0;

	  foreach my $pf_array(@{$pf->probe->get_all_Arrays}){
		
		if($pf_array->name eq $array->name){
		  $on_our_array = 1;
		  last;
		}
	  }

	  next if ! $on_our_array;#Go to next ProbeFeature

	  my $probeset = $pf->probeset->name;
	  
	  #Get the DBEntry information
	  #Filter for transcript by passing optional argument
	  foreach my $dbe(@{$pf->get_all_Transcript_DBEntries($transcript)}){
		#Probe can have different names on different arrays
		#So have to specify which array we are dealing with
		my $complete_name = $pf->probe->get_complete_name($array->name);
		my $slice_name = $pf->feature_Slice->name;
		my $dbe_info = $dbe->linkage_annotation;
		print "\t$complete_name $slice_name\t\tHits $dbe_info\n";
	  }
	}
  }
}


############################################
# Get All the probesets annotated 
# to a given transcripts
############################################
print "\nAll probesets annotated on ENSMUST00000015011:\n";
#Could also do this with the probe adaptor to get annotations for non probe set arrays
#e.g. Illumina etc
my $transcript = $tx_adaptor->fetch_by_stable_id('ENSMUST00000015011');
my @probesets = @{$pba->fetch_all_by_external_name($transcript->stable_id)};

foreach my $probeset (@probesets){

  my $arrays_string = join(', ', (map $_->name, @{$probeset->get_all_Arrays}));
  my $dbe_info;

  #Now get linkage_annotation
  foreach my $dbe(@{$probeset->get_all_Transcript_DBEntries($transcript)}){
	#This will return all ProbeSet DBEntries for this transcript
	#There should really be one max per transcript per probeset/probe
	$dbe_info = $dbe->linkage_annotation;
  }
  

  print "\t".$probeset->name." on arrays ".$arrays_string." with Probe hits $dbe_info\n";
  
}

############################################
# Probeset centric access
############################################
print "\nProbeset centric access. All transcripts annotated to Mouse430A_2 ProbeSet 1418625_s_at";
my $probeset = $pba->fetch_by_array_probeset_name('Mouse430A_2', '1418625_s_at');

my @trans_dbentries = @{$probeset->get_all_Transcript_DBEntries}; 

foreach my $dbe (@trans_dbentries){
  
  #Fetch transcript using primary ID of DBEntry object
  my $tx = $tx_adaptor->fetch_by_stable_id($dbe->primary_id);
  my $dbe_info = $dbe->linkage_annotation;

  #Grab the gene using the transcript stable ID stored in the DBEntry primary ID
  my $gene_sid = $gene_adaptor->fetch_by_transcript_stable_id($dbe->primary_id)->stable_id;

  print "\t".$dbe->db_display_name.' '.$dbe->primary_id."($gene_sid) at ".$tx->feature_Slice->name." with Probe hits $dbe_info\n";
 
}



##########################################
# Get details of unmapped objects
##########################################

print "\nUnmapped objects\n";

foreach my $ps_name('1437721_at', '1435628_x_at', '1418625_s_at'){
  my $ps = $pba->fetch_by_array_probeset_name($array->name, $ps_name);

  print "\n\tUnmappedObjects for ProbeSet $ps_name are:\n";
  &list_uos($ps->get_all_UnmappedObjects);

  foreach my $probe(@{$ps->get_all_Probes}){
	my $pname = $probe->get_complete_name($array->name);
	my @uos = @{$probe->get_all_UnmappedObjects};
	
	if(@uos){
	  	print "\tUnmappedObjects for Probe $pname are:\n";
		&list_uos(\@uos);
	  }

	foreach my $pf(@{$probe->get_all_ProbeFeatures}){
	  my $pname = $pf->probe->get_complete_name($array->name);
	  @uos = @{$pf->get_all_UnmappedObjects};
	
	  if(@uos){
		print "\tUnmappedObjects for ProbeFeature $pname ".$pf->feature_Slice->name.":\n";
		&list_uos(\@uos);
	  }
	}
  }
}

sub list_uos{
  my ($uos) = @_;

  foreach my $uo(@$uos){
	print "\t".$uo->identifier."\t".$uo->description."\n";

  }
}
