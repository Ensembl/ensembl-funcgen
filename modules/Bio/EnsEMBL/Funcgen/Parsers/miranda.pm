=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Funcgen::Parsers::miranda;

use strict;

use Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::ExternalFeature;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser);



sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_, type => 'miRanda');

  #Set default feature_type and feature_set config
  $self->{static_config}{feature_types} = 
	{
	 'miRanda Target'   => {
							-name        => 'miRanda Target',
							-class       => 'RNA',
							-description => 'miRanda microRNA target',
						   },
	};

  $self->{static_config}{analyses} =
	{
	 miRanda => { 
				 -logic_name    => 'miRanda',
				 -description   => '<a href="http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/consset.html">miRanda microRNA target predictions</a>',
				 -display_label => 'miRanda Targets',
				 -displayable   => 1,
				},
	};

  $self->{static_config}{feature_sets}{'miRanda miRNA targets'} = 
	{
	 #analyses      => $self->{static_config}{analyses},
	 #feature_types => $self->{static_config}{feature_types},
	 feature_set => 
	 {
	  -feature_type      => 'miRanda Target',
	  -display_name      => 'miRanda Targets',
	  -description       => $self->{static_config}{analyses}{miRanda}{-description},
	  -analysis          => 'miRanda',
	 },
	 xrefs => 1,
   }; 


  
 
  #$self->validate_and_store_feature_types;
  $self->validate_and_store_config([keys %{$self->{static_config}{feature_sets}}]);
  $self->set_feature_sets;

  return $self;
}


#TO DO
# In loop logging should only be enable with verbose or change to debug?
# Use count methods?
# Sort input to enable cache clearing(caches max out at ~200MB so not essential)
# Optimise slice cache testing(load only takes 8 min sso not essential)
# Add verbose logging/debug to show reassigned xrefs?

sub parse_and_load{
  my ($self, $files, $old_assembly, $new_assembly) = @_;

  #Add num files to config and check this in BaseImporter(generically)
  if(scalar(@$files) != 1){
	 throw('You must currently define a single file to load miRanda features from:\t'.join(' ', @$files));
  }

  my $file = $files->[0];
  $self->log_header("Parsing miRanda data from:\t$file");

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
  my $ftype_adaptor    = $self->db->get_FeatureTypeAdaptor();
  my $extf_adaptor     = $self->db->get_ExternalFeatureAdaptor;
  my $trans_adaptor    = $self->db->dnadb->get_TranscriptAdaptor;
  my $dbentry_adaptor  = $self->db->get_DBEntryAdaptor; 
  my $set              = $self->{static_config}{feature_sets}{'miRanda miRNA targets'}{feature_set};
 
  my ($dbentry, $schema_build, %features_by_name, %slice_cache, $ens_display_name, %feature_cache);
  my (@txs, @sid_dnames, %xref_cache, %seen_mifeat_trans, %skipped_features, %failed_reassignment);
  my ($overlaps_3_utr);
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'miRandaProjection');

  #Some hairy counting to make sure the XREF reassignment
  #is doing something sensible
  #Can probably change all this to use the count methods
  my $slice_skipped    = 0;
  my $old_mid_skipped  = 0;
  my $cnt              = 0;
  my $proj_skipped     = 0;
  my $seq_skipped      = 0;
  my $row_cnt          = 0;
  my $new_xrefs        = 0;
  my $old_as_new_xrefs = 0;
  my $failed_new_sids  = 0;
  my $utr_failed_old_sids = 0;
  my $no_proj_failed_old_sids = 0;
  my $previously_skipped = 0;
  my $with_overlapping_tx = 0;
  my $retired_sids         = 0;
  my $existing_sids    = 0;
  my $total_xrefs = 0;
  my $species          = $self->db->species;
  my $analysis         = $set->analysis;


  #Need to allow this to be passed
  $schema_build ||= $self->db->_get_schema_build($self->db->dnadb);

  if(! $species){
	throw('Must define a species to define the external_db');
  }
  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($species = lc($species)) =~ s/ /_/;

  my $edb_name        = $species.'_core_Transcript',



  open (FILE, "<$file") || die "Can't open $file";

  #We used to have redundant target features wrt xrefs
	  #Similarity      mmu-miR-192     miRanda miRNA_target    5       127885168       127885188       -       .       15.2253 4.048320e-03    ENSMUST00000031367      Slc15a4
	  #Similarity      mmu-miR-192     miRanda miRNA_target    5       127885168       127885188       -       .       15.2397 3.945600e-03    ENSMUST00000075376      Slc15a4
	  
  #Have now changed this to nr features with multiple xrefs
  #Hence have had to remove transcript sid from diplay_label

  #Several 'seen' caches have be implemented to prevent redundant storing
  #Currently these span the whole data set, but could be cleaned after every feature
  #if the input is sorted correctly

  #Need to make sure we are counting correctly
  #skipped counts will reflect lines in input, not nr features

  
  #Could add UnmappedObjects in here

 LINE: while (<FILE>) {
	next LINE if ($_ =~ /^\s*\#/o || $_ =~ /^\s*$/o);
	$row_cnt ++;
	#Added next for old miRbase IDs.

	#Sanger
	##GROUP SEQ     METHOD  FEATURE CHR     START   END     STRAND  PHASE   SCORE   PVALUE_OG       TRANSCRIPT_ID   EXTERNAL_NAME
	#Similarity      mmu-miR-707     miRanda miRNA_target    2       120824620       120824640       +       .       15.3548 2.796540e-02    ENST00000295228 INHBB


    my ($group, $id, $method, $feature, $chr, $start, $end, $strand, undef, undef, undef, $ens_id, $display_name) = split;
	#We never use $diplay_name now, as it never match the transcript display name 
	#which as the transcript number suffic attached to the gene display name
	#e.g. BRCA2-001

	#%seen_mifeat_trans handles redundancy between existant(i.e. not retired) sids in input
	#file and those identified when trying to reannotated a retired sid
	
	#ALREADY ANNOTATED OR SKIPPED
	if(exists  $seen_mifeat_trans{$id.':'.$chr.':'.$start.':'.$end}{$ens_id} ){
	  $self->log("Skipping previously reannotated miRNA target:\t".$id.':'.$chr.':'.$start.':'.$end.' - '.$ens_id);
	  $old_as_new_xrefs ++;
	}
	elsif(exists $skipped_features{$id.':'.$chr.':'.$start.':'.$end}){
	  $self->log("Skipping previous failed feature:\t".$id.':'.$chr.':'.$start.':'.$end."\t".
				 $skipped_features{$id.':'.$chr.':'.$start.':'.$end});
	  $previously_skipped++;
	  #Could increment count in here to count old xrefs failed for each fail type
	  next;
	}


	#Added next for old miRbase IDs.
	if ( $id =~ /\*$/o ){
	  $self->log("Skipping old miRbase ID:\t$id");

	  $skipped_features{$id.':'.$chr.':'.$start.':'.$end} = 'Old invalid miRbase ID';
	  $old_mid_skipped ++;

	  #$old_xrefs_skipped ++;
	  #Just want to make sure we have comparable
	  #old skipped xrefs and re-annotated xrefs
	  #so this is not useful here
	  

	  next LINE;
	}
	
    $strand = ($strand =~ /\+/o) ? 1 : -1;

	#Now moved transript xref info exclusively to xrefs
    ##my $id = $ens_id =~ s/[\"\']//g;  # strip quotes
	#my $id = $ens_id.':'.$seq;


	#change this to only test once
	#if exists and not defined then skip

	if(! defined $slice_cache{$chr}){
	
	  #Was originally limiting to chromosome

	  if($old_assembly){
		$slice_cache{$chr} = $self->slice_adaptor->fetch_by_region(undef, 
																	$chr, 
																	undef, 
																	undef, 
																	undef, 
																	$old_assembly);
	  }else{
		$slice_cache{$chr} = $self->slice_adaptor->fetch_by_region(undef, $chr);
	  }

	  if(! defined $slice_cache{$chr}){
		warn "Can't get slice $chr for sequence $id\n";

		$slice_skipped ++;
		$skipped_features{$id.':'.$chr.':'.$start.':'.$end} = "Failed to fetch slice $chr";
		
		#Add UnmappedObject here?
		next LINE;
	  }
	}

	
	#We can add coding xref to feature type based on the miRbase name
	#.e.g hsa-mir-24-1
	#However, this isn't stored as an xref
	#It is stored in the gene.description
	#e.g. hsa-mir-24-1 [Source:miRBase;Acc:MI0000080]
	#Not easy to fetch as descriptions not indexed!
	#

	#Cache/store FeatureType

	if(! exists $features_by_name{$id}){
	  $features_by_name{$id} = $ftype_adaptor->fetch_by_name($id);
	  
	  if(! defined $features_by_name{$id}){
		($features_by_name{$id}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
															 (
															  -name  => $id,
															  -class => 'RNA',
															  -description => $method.' '.$feature,
															 ))};

		#Need to add source gene xref here to enable target conequences implied by source variation

	  }
	}


	#Make sure we have the target transcript before we store the feature
	#Have to do this as we can't always run with the correct core DB
	#as it may be too old. Hence we have to hard code the edb.release
	
	
	##This should enever happen, as the search regions are defined by ens transcript
	##i.e there is always a ensembl iD
	
	if (! $ens_id) {
      warn("No xref available for miRNA $id\n");
      $skipped_xref++;
      next;
    }
   


	if(exists $feature_cache{$id.':'.$chr.':'.$start.':'.$end}){
	  $feature = $feature_cache{$id.':'.$chr.':'.$start.':'.$end};
	}
	else{

	  $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
		(
		 -display_label => $id,
		 -start         => $start,
		 -end           => $end,
		 -strand        => $strand,
		 -feature_type  => $features_by_name{$id},
		 -feature_set   => $set,
		 -slice         => $slice_cache{$chr},
		);

	  # project if necessary
	  if ($new_assembly) {
		
		$feature = $self->project_feature($feature, $new_assembly);
		
		if (! defined $feature) {
		  $proj_skipped ++;
		  $skipped_features{$id.':'.$chr.':'.$start.':'.$end} = 'Failed projection';
		  next;
		}		
		
		#This was failing as old assembly seqs are not stored, hence are returned as Ns
		#if ($old_seq ne $feature->seq){
		#  $skipped_features{$id.':'.$chr.':'.$start.':'.$end} = 'Projected sequence mismatch';
		#  $seq_skipped++;
		#  next;
		#}

		#Assembly mapping tends to ignore single seq mismatches, but does handle gaps.
		#This is not generally a problemas a new assembly is normally a reshuffle of existing
		#clones/contigs which will have exactly the same seq(actually there is a small amount of change)
		#but this only affected a handful of transcripts
		#old super(contigs) which have be integrated into a chromosome are not mapped
		#so we will lose these
		#the rest (retired clones/contigs and new clones/contigs) are aligned
		#but again mismatches tend to be ignored.
		
		#There may be a tendency for this later case to contain more mismatches due to an entirely new
		#clone seqeunce, hence let's filter/count these?

		#Tested and import using the following hack using old 67 DB to test seq mismatches
		#only 1 failed.
		#my $old_slice = $v67_sa->fetch_by_region( undef, $chr, $start, $end, $strand);

		#if ($old_slice->seq ne $feature->seq){
		#  $skipped_features{$id.':'.$chr.':'.$start.':'.$end} = 'Projected sequence mismatch';
		#  $seq_skipped++;
		#  next LINE;
		#}
	  }
	}



	### DEFINE TRANSCRIPT XREF INFO

	#Check transcript exists
	$display_name = $self->get_core_display_name_by_stable_id($self->db->dnadb, $ens_id, 'transcript');
	

	# ATTEMPT TO REANNOTATE
	if (! defined $display_name){ # Transcript does not exist in current release

	  #Set to 0 as we are not storing this sid
	  $seen_mifeat_trans{$id.':'.$chr.':'.$start.':'.$end}{$ens_id} = 0;

	  $self->log("$id $ens_id stable ID has been retired");
	  $retired_sids ++;

	  #Try and re-annotate on newer overlapping transcripts
	  @txs = @{$trans_adaptor->fetch_all_by_Slice($feature->feature_Slice)};
	  @sid_dnames = ();
	  
	  #COUNT unique miRNA features which have failed reannotated
	  if(! exists $failed_reassignment{$id.':'.$chr.':'.$start.':'.$end}){
		$failed_reassignment{$id.':'.$chr.':'.$start.':'.$end} = 1;
	  }



	  if (@txs) { #OVERLAPPING TRANSCRIPTS
		$with_overlapping_tx ++;

		foreach my $tx(@txs){

		  #Check we have previously seen this xref
		  if(! exists $seen_mifeat_trans{$id.':'.$chr.':'.$start.':'.$end}{$tx->stable_id}){

			 #Do UTR checking here
			 $overlaps_3_utr = 0;

			 #Transcript will always return wrt to feature_Slice of miRNA feature.
			 #As miRNA are complimentary to the the mRNA, which is complimentary
			 #to the sense strand of the gene, these should always be on the same strand.

			 #Could add some more detailed strand/UTR fail counts in here
			 #but leave for now
			 #$same_strand = 1;

			 if($tx->seq_region_strand != $strand){
			   #$same_strand = 0;
			 }
			 elsif($strand == 1){
			   
			   if( ($end <= $tx->seq_region_end) &&
				   ($start > $tx->coding_region_end) ){
				 $overlaps_3_utr = 1;
			   }

			 }
			 else{#Must be -1
			   
			   if( ($end < $tx->coding_region_start) &&
				   ($start >= $tx->seq_region_start) ){
				 $overlaps_3_utr = 1;
			   }
			 }
			
			 #could count no utr match here if we set same_strand boolean

			 if($overlaps_3_utr){
			   $display_name = $self->get_core_display_name_by_stable_id($self->db->dnadb, 
																		 $tx->stable_id, 
																		 'transcript');
			   push @sid_dnames, [$tx->stable_id, $display_name];
			   $seen_mifeat_trans{$id.':'.$chr.':'.$start.':'.$end}{$tx->stable_id} = 1;
			   #want to count re-annotated xrefs here
			   #These may have been represented in the original file
			   $new_xrefs ++;
			 }
			 else{
			   # FAILED annotate new sid
			   # This currently includes +ve and -ve strand transcripts
			   $failed_new_sids ++;
			 }
		   }
		   
		   if(! @sid_dnames){
			 #FAILED TO REANNOTATE retired sid xref
			 $utr_failed_old_sids ++;

			 next LINE;
		   }
		 }
	  }
	  else{ #FAILED TO REANNOTATE XREF - NO OVERLAPPING TRANSCRIPTS
		$no_proj_failed_old_sids ++;
		next LINE;
	  }


	  #COUNT unique miRNA features which have failed reannotated
	  if(@sid_dnames){
		$failed_reassignment{$id.':'.$chr.':'.$start.':'.$end} = 0
	  }

	}
	else { #Add xref for existing transcript	  
	  $seen_mifeat_trans{$id.':'.$chr.':'.$start.':'.$end} = {$ens_id => 1};
	  #$display_name = $ens_display_name;
	  @sid_dnames = ([$ens_id, $display_name]);
	  $existing_sids++;
	}


	#shouldn't need this if as we call next for everycase above
	#if (@xref_details) {
	#Only if we have target transcripts

	# STORE FEATURE
	if (! exists $feature_cache{$id.':'.$chr.':'.$start.':'.$end}){
	  ($feature) = @{$extf_adaptor->store($feature)};
	  $feature_cache{$id.':'.$chr.':'.$start.':'.$end} = $feature;
	  $cnt++;
	}

 
	# STORE XREFS
	foreach my $xref_info(@sid_dnames) {
	  
	  #Handle release/version in xref version as stable_id version?

	  $dbentry = Bio::EnsEMBL::DBEntry->new
		(
		 -dbname                 => $edb_name,
		 -release                => $schema_build,
		 #-release                => '58_37k',#'46_36h', #Hard coded due to schema to old to use with API
		 -status                 => 'KNOWNXREF',
		 #-display_label_linkable => 1,
		 -db_display_name        => 'EnsemblTranscript',
		 -type                   => 'MISC',
		 -primary_id             => $xref_info->[0],
		 -display_id             => $xref_info->[1],
		 -info_type              => 'MISC',
		 -info_text              => 'TRANSCRIPT',
		 -linkage_annotation     => 'miRanda target - negative influence',
		 #could have version here if we use the correct dnadb to build the cache
		 -analysis               =>  $analysis,
		);
	
	  $dbentry_adaptor->store($dbentry, $feature->dbID, 'ExternalFeature', 1); #1 is ignore release flag  
	  $total_xrefs++;
	}
  }

  close FILE;
  
  #scalar context returns count
  my $miRNA_failed_reassignment = grep/1/, values %failed_reassignment;
  my $miRNA_skipped = $old_mid_skipped + $slice_skipped + 
	$proj_skipped + $seq_skipped + $miRNA_failed_reassignment;
 
  $self->log_header($set->name." Import Report");
  $self->log(sprintf("%-090s", "miRNA feature:target pairs seen(i.e. input rows):").$row_cnt);
  $self->log(sprintf("%-090s","Stored NR miRanda miRNA ExternalFeatures:").$cnt);
  $self->log(sprintf("%-090s","Total skipped miRanda miRNA target features(inc reassigned):").$miRNA_skipped);
  $self->log(sprintf("%-090s","Old miRbase IDs skipped").$old_mid_skipped);
  $self->log(sprintf("%-090s","Skipped features on unknown slice:").$slice_skipped."\n");

  if($new_assembly){
	$self->log(sprintf("%-090s","Skipped due to failed assembly projection:").$proj_skipped);
	$self->log(sprintf("%-090s","Skipped due to seq mismatch for assembly projection:")."$seq_skipped\n");
  }

  $self->log("The following numbers are counted from the mappable/valid miRNA target features");
  $self->log(sprintf("%-090s", "Total stored Transcript xrefs(current and retired re-assigned):"). $total_xrefs);
  $self->log(sprintf("%-090s", "Total current Transcript xrefs:"). $existing_sids);
  $self->log(sprintf("%-090s", "Total skipped miRanda miRNA target xrefs:").($previously_skipped	+ $miRNA_skipped));
  $self->log(sprintf("%-090s", "Retired Transcript xrefs:").$retired_sids);
  $self->log(sprintf("%-090s", "Unique miRNA features which completely failed reassignment:").$miRNA_failed_reassignment);
  $self->log(sprintf("%-090s", "Total new Xrefs assigned due to retired Transcript:").$new_xrefs); 
  $self->log(sprintf("%-090s", "Previously assigned Transcript Xrefs skipped:").$old_as_new_xrefs);  
  $self->log(sprintf("%-090s", "Retired Transcript Xrefs with no new overlapping Transcript:").$no_proj_failed_old_sids);
  $self->log(sprintf("%-090s", "Retired Transcript Xrefs with new overlapping Transcript(s):").$with_overlapping_tx);
  $self->log(sprintf("%-090s", "Retired Transcript Xrefs with new overlapping Transcript(s), all fail 3'UTR/strand test:").
			 $utr_failed_old_sids);	
  #This figure may include transcripts which would have failed in the original set!
  $self->log(sprintf("%-090s","Total new Transcript xrefs considered which fail 3' UTR/strand test:"). $failed_new_sids);
  $self->log(sprintf("%-090s","True new Xref assignments due to retired Transcripts:").($new_xrefs - $old_as_new_xrefs));

  #We also want unique miRNA feature which we re-assigned
  #Hard to calculate as the re-assignment might be to a current transcript in the input

  #Now set states
  foreach my $status(qw(DISPLAYABLE MART_DISPLAYABLE)){
	$set->adaptor->store_status($status, $set);
  }

  return;
}

1;
