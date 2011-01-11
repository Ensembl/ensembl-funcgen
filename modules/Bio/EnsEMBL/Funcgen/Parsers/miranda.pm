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

package Bio::EnsEMBL::Funcgen::Parsers::miranda;

use strict;

#This was in the description
#For MSKCC miRanda targets
#Contact is Doron Betel (or maybe Mark Levenstien levenstm@mskcc.org)
#http://www.microrna.org/microrna/getDownloads.do
#Current Human/Mouse/Rat



#This is what we have actually hosted
#For Sanger miRanda targets
#contact svd@sanger.ac.uk stijn van dongen
#From our download page
#http://microrna.sanger.ac.uk/cgi-bin/targets/v5/download.pl
#and then
#ftp://ftp.sanger.ac.uk/pub/mirbase/targets/v5/arch.v5.gff.homo_sapiens.zip



#  #<GROUP>	<SEQ>	<METHOD>	<FEATURE>	<CHR>	<START>	<END>	<STRAND>	<PHASE>	<SCORE>	
#  Similarity	hsa-miR-23b	miRanda	miRNA_target	1	919788	919807	+	.	69	transcript id "ENST00000310998"
#  Similarity	hsa-miR-23a	miRanda	miRNA_target	1	919787	919807	+	.	71	transcript id "ENST00000310998"

##GROUP SEQ     METHOD  FEATURE CHR     START   END     STRAND  PHASE   SCORE   PVALUE_OG       TRANSCRIPT_ID   EXTERNAL_NAME
#Similarity      mmu-miR-707     miRanda miRNA_target    2       120824620       120824640       +       .       15.3548 2.796540e-02    ENST00000295228 INHBB

#But we now want to host a reduced set from Anton

use Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::ExternalFeature;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser);

# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors




sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #Set default feature_type and feature_set config
  $self->{'feature_types'} = {
							  'miRanda Target'   => {
											  name        => 'miRanda Target',
											  class       => 'RNA',
											  description => 'miRanda microRNA target',
											 },
							 };

  $self->{feature_sets} = {
						   'miRanda miRNA targets' => 
						   {
							feature_type      => \$self->{'feature_types'}{'miRanda Target'},
							display_name      => 'miRanda Targets',
							analysis          => 
							{ 
							 -logic_name    => 'miRanda',
							 #-description   => 'miRanda microRNA target prediction (http://www.microrna.org)',
							 -description => 'miRanda microRNA target prediction (http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/)',
							 -display_label => 'miRanda Target',
							 -displayable   => 1,
							},
							xrefs => 1,
						   },						   
						  };

  
 
  $self->validate_and_store_feature_types;
  $self->set_feature_sets;

  return $self;
}




sub parse_and_load{
  my ($self, $files, $old_assembly, $new_assembly) = @_;

  my $file = $files->[0];
  $self->log_header("Parsing miRanda data from:\t$file");  
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
  my $ftype_adaptor    = $self->db->get_FeatureTypeAdaptor();
  my $extf_adaptor     = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_adaptor  = $self->db->get_DBEntryAdaptor; 
  my $set              = $self->{'feature_sets'}{'miRanda miRNA targets'};
  my %features_by_name; # name -> feature_type
  my %slice_cache;
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'miRandaProjection');
  my $skipped = 0;
  my $cnt = 0;
  my $skipped_xref = 0;
  my $species = $self->db->species;
  if(! $species){
	throw('Must define a species to define the external_db');
  }
  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($species = lc($species)) =~ s/ /_/;



  open (FILE, "<$file") || die "Can't open $file";

 LINE: while (<FILE>) {
	next LINE if ($_ =~ /^\s*\#/o || $_ =~ /^\s*$/o);

	#Sanger
	##GROUP SEQ     METHOD  FEATURE CHR     START   END     STRAND  PHASE   SCORE   PVALUE_OG       TRANSCRIPT_ID   EXTERNAL_NAME
	#Similarity      mmu-miR-707     miRanda miRNA_target    2       120824620       120824640       +       .       15.3548 2.796540e-02    ENST00000295228 INHBB


	#MSKCC
	#UCSC ID	mRNA	Gene ID	miRNA acc	miRNA	miRNA align	alignment	gene align	align score	conservation	miRNA start	miRNA end	gene start	gene end	%ID	%Similar	energy	organism	prediction_date
	#uc001abs.1	AK091100	LOC643837	MIMAT0000062	hsa-let-7a	uuGAUAUGUUGGAUGAUGGAGu	||:  || |:| ||:|||| 	guCUGCUCACCUUCCUGCCUCa	144	0.615759	2	21	396	417	63	78	0	9606	2008-05-16
	#uc001afh.1	NM_001039577	CCNL2	MIMAT0000062	hsa-let-7a	uugauAUGUUGGAUGAUGGAGu	|::||||| ||:|||| 	cucacUGUAACCU-CUGCCUCc	147	0.767611	2	18	2290	2310	75	93	0	9606	2008-05-16



    my ($group, $seq, $method, $feature, $chr, $start, $end, $strand, undef, undef, undef, $ens_id, $display_name) = split;
    $strand = ($strand =~ /\+/o) ? 1 : -1;
    ##my $id = $ens_id =~ s/[\"\']//g;  # strip quotes
	my $id = $ens_id.':'.$seq;

	if(! defined $slice_cache{$chr}){
	
	  if($old_assembly){
		$slice_cache{$chr} = $self->slice_adaptor->fetch_by_region('chromosome', 
																	$chr, 
																	undef, 
																	undef, 
																	undef, 
																	$old_assembly);
	  }else{
		$slice_cache{$chr} = $self->slice_adaptor->fetch_by_region('chromosome', $chr);
	  }

	  if(! defined $slice_cache{$chr}){
		warn "Can't get slice $chr for sequence $id\n";
		$skipped++;
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

	if(! exists $features_by_name{$seq}){
	  $features_by_name{$seq} = $ftype_adaptor->fetch_by_name($seq);
	  
	  if(! defined $features_by_name{$seq}){
		($features_by_name{$seq}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
															 (
															  -name  => $seq,
															  -class => 'RNA',
															  -description => $method.' '.$feature,
															 ))};

		#Need to add source gene xref here to enable target conequences implied by source variation

	  }
	}

   
	$feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -display_label => $id,
	   -start         => $start,
	   -end           => $end,
	   -strand        => $strand,
	   -feature_type  => $features_by_name{$seq},
	   -feature_set   => $set,
	   -slice         => $slice_cache{$chr},
	  );

	

    # project if necessary
    if ($new_assembly) {
      $feature = $self->project_feature($feature, $new_assembly);

	  if(! defined $feature){
		$skipped ++;
		next;
	  }
    }

	($feature) = @{$extf_adaptor->store($feature)};
	$cnt++;

 
	#Build xref

	#This should enever happen, as the search regions are defined by ens transcript
	if (! $ens_id) {
      warn("No xref available for miRNA $id\n");
      $skipped_xref++;
      next;
    }


	#use external_name first, else try and get it from the core DB
	#should we just get it from the core DB anyway?

	if(! defined $display_name){
	  $display_name = $self->get_core_display_name_by_stable_id($self->db->dnadb, $ens_id, 'transcript');
	}


	#Handle release/version in xref version as stable_id version?

	my $dbentry = Bio::EnsEMBL::DBEntry->new(
											 -dbname                 => $species.'_core_Transcript',
											 #-release                => $self->db->dnadb->dbc->dbname,
											 -release => '46_36h', #Hard coded due to schema to old to use with API
											 #-release => '46_36g',
											 -status                 => 'KNOWNXREF',
											 #-display_label_linkable => 1,
											 -db_display_name        => 'EnsemblTranscript',
											 -type                   => 'MISC',
											 -primary_id             => $ens_id,
											 -display_id             => $display_name,
											 -info_type              => 'MISC',
											 -info_text              => 'TRANSCRIPT',
											 #-linkage_annotation     => 'miRanda miRNA negative influence',
											 -linkage_annotation     => 'miRanda target - negative influence',
											 #could have version here if we use the correct dnadb to build the cache
											 -analysis  =>  $self->{'feature_sets'}{'miRanda miRNA targets'}{analysis},
											);

	$dbentry_adaptor->store($dbentry, $feature->dbID, 'ExternalFeature', 1);#1 is ignore release flag  
  }

  close FILE;

  $self->log("Stored $cnt miRanda miRNA ExternalFeatures");
  $self->log("Skipped $skipped miRanda miRNA imports");
  $self->log("Skipped an additional $skipped_xref DBEntry imports");


  #Now set states
  foreach my $status(qw(DISPLAYABLE MART_DISPLAYABLE)){
	$set->adaptor->store_status($status, $set);
  }

  return;
}

1;
