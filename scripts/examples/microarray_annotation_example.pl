#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

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

=head1 NAME

micro_array_annotation_example.pl

=head1 DESCRIPTION

This script extracts probes mapping and probeset annotations for the Mouse
Affy array MOE430A, on the chromosome 2:26771920:26789448

=cut

# Remember to add the ensembl and ensembl-funcgen packages to your PERL5LIB


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Config::Tiny;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use feature qw(say);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(dump_data);

####################################
# Use the Registry
# Load database details
# Then get the slice/array adaptors
####################################

my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',

my $slice_adaptor = $registry->get_adaptor("mus musculus","core","Slice");
my $pfa           = $registry->get_adaptor("mus musculus","funcgen","ProbeFeature");
my $aa            = $registry->get_adaptor("mus musculus","funcgen","Array");
my $pba           = $registry->get_adaptor("mus musculus","funcgen","ProbeSet");
my $tx_adaptor    = $registry->get_adaptor("mus musculus","core","Transcript");
my $gene_adaptor  = $registry->get_adaptor("mus musculus","core","Gene");

#my $cfg = Config::Tiny->new;
#   $cfg = Config::Tiny->read('../test_db.ini');
#
#my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new (
#  -user       => $cfg->{staging}->{user},
#  # -pass       => $cfg->{staging}->{pass},
#  -host       => $cfg->{staging}->{host},
#  -port       => $cfg->{staging}->{port},
#  -dbname     => $cfg->{staging}->{dbname},
#  -dnadb_host => $cfg->{staging}->{host},
#  -dnadb_user => $cfg->{staging}->{user},
#  );
#
#my $dna_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new (
#  -user       => $cfg->{dna_db}->{user},
#  # -pass       => $cfg->{staging}->{pass},
#  -host       => $cfg->{dna_db}->{host},
#  -port       => $cfg->{dna_db}->{port},
#  -dbname     => $cfg->{dna_db}->{dbname},
#  );


#print "Using api: ",$registry->software_version,"\n";

my $slice_adaptor = $db->get_SliceAdaptor;
my $pfa           = $db->get_ProbeFeatureAdaptor;
my $aa            = $db->get_ArrayAdaptor;
my $pba           = $db->get_ProbeSetAdaptor;
my $tx_adaptor    = $dna_db->get_TranscriptAdaptor;
my $gene_adaptor  = $dna_db->get_GeneAdaptor;

#####################################################
# Get the slice for gene Cntnap1 with 100000bp flanks
# Get the MOE430A array
#####################################################

my ($gene)  = @{$gene_adaptor->fetch_all_by_external_name('Cntnap1')};
my $slice = $gene->feature_Slice->expand(10000, 10000);
my $array = $aa->fetch_by_name_vendor('HuEx-1_0-st-v2', 'AFFY');

############################################
# Get mapped probes for that slice and array
############################################

my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);

print "\n".scalar(@$oligo_features)." ProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n\n";

foreach my $pf(@$oligo_features){
  my $probeset   = $pf->probeset->name;
  my $probename  = $pf->probe->get_probename($array->name);
	my $slice_name = $pf->feature_Slice->name;
	#extended cigar_string as defined by SAMTools group
  print "\tProbe $probeset $probename is aligned on $slice_name with cigar string ".$pf->cigar_string."\n";
}

############################################
# Get all ProbeFeatures annotated on the probesets annotated to Transcripts
# Note: ProbeFeature level DBEntries are used to build Probe/ProbeSet
# level transcript annotatations.
print "\nProbeFeatures xref'd to Transcripts on slice ".$slice->name."\n";
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
print "\nAll probesets annotated on ENST00000515683:\n";
#Could also do this with the probe adaptor to get annotations for non probe set arrays
#e.g. Illumina etc
my $transcript = $tx_adaptor->fetch_by_stable_id('ENST00000515683');
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

#############################################################
# Get all probesets annotated to transcripts for a given gene
#############################################################

print "\nAll probesets associated with transcript of Cntnap1:\n";


foreach my $pset(@{$pba->fetch_all_by_linked_transcript_Gene($gene)}){
	print $pset->name."\n";

	#Use the loop above to get all of the xref info for each probeset,
	#which may include mappings to other genes/transcript
}




############################################
# Probeset centric access
############################################
print "\nProbeset centric access. All transcripts annotated to  ProbeSet 1418625_s_at";
my $probeset = $pba->fetch_by_array_probeset_name('HuGene-2_0-st-v1', '16740948');

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
###########################################
#
#print "\nUnmapped objects\n";
#
#foreach my $ps_name('1960566' ){
#  my $ps = $pba->fetch_by_array_probeset_name($array->name, $ps_name);
#
#  print "\n\tUnmappedObjects for ProbeSet $ps_name are:\n";
#  &list_uos($ps->get_all_UnmappedObjects);
#
#  foreach my $probe(@{$ps->get_all_Probes}){
#	my $pname = $probe->get_complete_name($array->name);
#	my @uos = @{$probe->get_all_UnmappedObjects};
#
#	if(@uos){
#	  	print "\tUnmappedObjects for Probe $pname are:\n";
#		&list_uos(\@uos);
#	  }
#
#	foreach my $pf(@{$probe->get_all_ProbeFeatures}){
#	  my $pname = $pf->probe->get_complete_name($array->name);
#	  @uos = @{$pf->get_all_UnmappedObjects};
#
#	  if(@uos){
#		print "\tUnmappedObjects for ProbeFeature $pname ".$pf->feature_Slice->name.":\n";
#		&list_uos(\@uos);
#	  }
#	}
#  }
#}

sub list_uos{
  my ($uos) = @_;

  foreach my $uo(@$uos){
	print "\t".$uo->identifier."\t".$uo->description."\n";

  }
}
