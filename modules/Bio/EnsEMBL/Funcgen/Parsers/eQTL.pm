#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::eQTL
#

#Could this be based on a Generic Flat file parser?

=head1 NAME

Bio::EnsEMBL::Funcgen::Parsers::eQTL

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::eQTL";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  eQTL contains meta 
data and methods specific to data in bed format, to aid 
parsing and importing of experimental data.

=head1 AUTHOR

This module was created by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

#To do
# 1 Remove the RegulatoryFeature xrefs? No we just need a script to calculate these independantly.


package Bio::EnsEMBL::Funcgen::Parsers::eQTL;

use Bio::EnsEMBL::Funcgen::Parsers::ExperimentalSet;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use strict;


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::ExperimentalSet);


my $reg = 'Bio::EnsEMBL::Registry';

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for eQTL class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::eQTL
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;
  
  my $class = ref($caller) || $caller;

  #define default fields here and pass
  #We also need to be able to take custom attrs mappings

  #keys are array index of field, key are Feature paramter names
  #reverse this?
  #Unless we have a compound field which we name accordingly
  #And then call e.g. parse_attrs
  #Which will return a hash with the relevant Feature attributes

  #Is splitting this up simply going to make the parse slower due to acessor methods?

  #Pass or just set directly here?
  #SNP_label       GENE_Label      GENE_ID SNP_chromosome  GENE_chromosome SNP_location    GENE_location   Distance        pvalue  -log10(pvalue)  Adjusted_R^2    Gradient



  my $self  = $class->SUPER::new(@_);
  
  my ($store_regf_xrefs) = rearrange(['STORE_REG_FEAT_XREFS'], @_);

  #my ($display_label_field) = rearrange(['DISPLAY_LABEL_FIELD'], @_);

  #We need to define meta header method, starting with '##'
  #Also need to skip comments '#' at begining or end of line
  #Do we also need to skip field header? No methinks not.

  #Define result method
 # $self->{'file_ext'} => 'gff';#Could use vendor here?
  
  #define this if we want to override the generic method in Simple
  #$self->{'config'}{'results_data'} => ["and_import_gff"];  

  #$self->display_label_field($display_label_field);
  

  $self->{'store_reg_feat_xrefs'} = 1; #$store_regf_xrefs;

  return $self;
}

#we surely use thes only once in the code
#is it faster to cache like this over the reg method
#or should we just use reg directly?

sub regulatory_feature_adaptor{
  my $self = shift;
  return $self->{'regulatory_feature_adaptor'};
}

sub gene_adaptor{
  my $self = shift;
  return $self->{'gene_adaptor'};
}

sub variation_adaptor{
  my $self = shift;
  return $self->{'variation_adaptor'};
}

sub variation_feature_adaptor{
  my $self = shift;
  return $self->{'variation_feature_adaptor'};
}


=head2 set_config

  Example    : my $self->set_config;
  Description: Sets attribute dependent config
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_config{
  my $self = shift;

  $self->SUPER::set_config;

  #some convenience methods
  $self->{'regulatory_feature_adaptor'} = $self->db->get_RegulatoryFeatureAdaptor;
  $self->{'variation_adaptor'} = $reg->get_adaptor($self->species, 'variation', 'Variation');
  $self->{'variation_feature_adaptor'} = $reg->get_adaptor($self->species, 'variation', 'VariationFeature');
  $self->{'gene_adaptor'} = $self->db->dnadb->get_GeneAdaptor;

   
  my $species = $self->db->species;
  if(! $species){
	throw('Must define a species to define the external_db');
  }
  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($self->{'species'} = lc($species)) =~ s/ /_/;
  return;
}

#The flow of this module requires a gene sorted list
#default is position sorted.

sub pre_process_file{
  my ($self, $filepath) = @_;

  #$self->backup_file($file_name)
  #Let's just work with a copy rather than backing up.



  #We could just do the sort on the fly as with gzip
  my $new_path = $filepath.'.sorted';
  my $sort = 'sort -k 2 '.$filepath.' > '.$new_path;
  run_system_cmd($sort);

  return $new_path;
}

sub parse_line{
  my ($self, $line) = @_;

  return if $line =~ /^SNP_label/;
  #can we test for $. here instead, as header may change.
  #we'll have to recode for the changed header anyway?
  #is a $. check quick than a ~= /^SNP_label/ test?
  #we should validate the header

  

  #my %fields = (
  #				0 => 'fetch_slice',
  #				1 => 'get_source',
  #				2 => 'get_feature_type',
  #				3 => '-start',
  #			4 => '-end',
  #				5 => '-strand',#Will most likely be , need to convert to -.+ > -1 0 1
  #6 => 'frame',#will most likely be .
  #				7 => 'get_attributes',
  #			   );
  
 #SNP_label       GENE_Label      GENE_ID SNP_chromosome  GENE_chromosome SNP_location    GENE_location   Distance        pvalue  -log10(pvalue)  Adjusted_R^2    Gradient

  my ($gene);
    
  my ($rsid, $gene_name, undef, $gene_chr, $snp_chr, $snp_start, $gene_start) = split/\t/o, $line;			
  #currently ignoring GENE_ID, distance, pvalue, -log10(pvalue), Adjusted_R^2 and Gradient.
  

  #This appears to slow things down a lot
  #Hack to handle scientific notation
  #$gene_chr = sprintf('%d', $gene_chr);
  #$snp_chr = sprintf('%d', $snp_chr);
  #$gene_start = sprintf('%d', $gene_start);
  #$snp_start = sprintf('%d', $snp_start);


  if($self->ucsc_coords){
	$gene_start += 1;
	$snp_start  += 1;
  }

  #
  #we need to return feature_params, DBEntry params( and seq if defined?)
  #return seq in feature_params and delete from hash before storing

  #Store as experimental sets with annotated features, with DBEntries to Variation(or DB SNP?) and core.
  #How do we want to represent this as a feature? Simply as the SNP? loci?
  #Or as the region between all of the SNPs and the associated gene.
  #This then becomes a multiline format
  #We need to test that SNPs are on sme chromosome, and that e-QTLs between genes do not overlap
  #Can we assume overlap of whole region with reg feat is enough to assign to gene, or do we need SNP overlap?
  
  #So we need to maintain internal params cache for multiline format which we will return explicitly one discovery of a new gene name
  #Then we just need to empty the cache in Simple once we have finished the file.

  #We are going to be calling SNP and Gene adaptor for each record
  #Will internal cache be faster than registry
  #Also need to load registry from ensembldb to access variation DB.

  #Need internal attr last_gene, or generic 'feature_separator'

  
  if( ! defined $self->feature_separator || ($self->feature_separator ne $gene_name) ){
	#Copy params data

	#We need to have an internal method here to do the final processing
	#This will call the store methods in the parent, but also allow and 
	#format specific stats for the last record, which would other wise be done in the caller?
	#or can 

	
	#Only process if we have seen a gene
	$self->process_params if defined $self->feature_separator;
  
	
	#Now set default params for new feature
	$self->set_feature_separator($gene_name);




	
	#Get gene and validate
	my @genes = @{$self->gene_adaptor->fetch_all_by_external_name($gene_name)};
	$self->log('Found '.scalar(@genes).' genes for external name '.$gene_name) 	if(scalar(@genes) > 1);	

	my $warn_txt = '';
	my $disp_gene_cnt = 0;
	my $nondisp_gene_cnt = 0;
	my ($tmp_gene, $display_name);

	foreach my $ext_gene(@genes){
	  
	  #validate gene loci
	  if($ext_gene->seq_region_name ne $gene_chr && $ext_gene->start != $gene_start){
		$warn_txt .= "eQTL $gene_name gene loci does not match Ensembl(".$ext_gene->stable_id.") loci:\t".
		  $gene_chr.':'.$gene_start.' vs '.$ext_gene->seq_region_name.':'.$ext_gene->start."\n";
	  }
	  else{
	
		if($ext_gene->display_xref && ($ext_gene->display_xref->display_id eq $gene_name)){
		  $disp_gene_cnt++;
		  $gene = $ext_gene;
		}
		else{
		  $warn_txt .= "eQTL $gene_name is not display_xref for loci validated Ensembl gene ".$ext_gene->stable_id."\n";
		  $nondisp_gene_cnt++;
		  $tmp_gene = $ext_gene;
		}
	  }
	}


	if($gene){
	  throw('Found more than one loci/display_xref validated gene for '.$gene_name) if $disp_gene_cnt > 1;
	}
	else{
	  $self->log($warn_txt) if $warn_txt;
	  
	  if(! $tmp_gene){
		#log missed gene.
		#Still build the eQTL
		#have counts hash attr, but allow dynamic naming of count elements, to allow global counting and reset in Simple
		$self->count('skipped_genes');
		$self->log('Failed to find comparable Ensembl gene for '.$gene_name.' Building sparse eQTL');
	  }
	  else{
		$gene = $tmp_gene;
				
		if($nondisp_gene_cnt == 1){
		  #Only reset the dispay name if it maps to one which is not the display_xref
		  $display_name = $gene->display_xref->display_id if $gene->display_xref;
		}
	  }
	}

	$display_name ||= $gene_name;


	#Set generic feature params
	$self->feature_params->{'-feature_type'}  = $self->data_set->product_FeatureSet->feature_type;
	$self->feature_params->{'-display_label'} = 'e-QTL:'.$display_name;	
	#These are generic to the whole import so move to load?
	$self->feature_params->{'-feature_set'}   = $self->data_set->product_FeatureSet;
	$self->feature_params->{'-analysis'}      = $self->data_set->product_FeatureSet->analysis;
	#No score attributable to feature?
	#we could simply count the number of SNPs as the score?
	$self->feature_params->{'-score'} = 0;
	#setting strand to 0 for now
	$self->feature_params->{'-strand'} = 0;


	if($gene){
	  #???????????????????????????????????????????????????????????????
	  #Use ensembl's label to avoid propogation of xref records?
	  #This is masking original data?
	  #Altho' it is loci based.
	  
	  #This shouldn't need testing as is should have a display xref
	  #my $display_label = $gene->display_xref->display_label;
	  
	  #Now set up first DBEntry params
	  my %gene_xref_params = 
		(
		 -dbname                 => $self->{'species'}.'_core_Gene',
		 #-release                => $self->db->dnadb->dbc->dbname,
		 -status                 => 'KNOWNXREF',
		 #-display_label_linkable => 1,
		 #-db_display_name        => $self->db->dnadb->dbc->dbname,
		 -db_display_name        => 'EnsemblGene',
		 -type                   => 'MISC',#edb.type???
		 -primary_id             => $gene->stable_id,
		 #This may cause duplications in xref, as we are assigning the gene_name
		 #to display_id which is part of the unique key
		 -display_id             => ($gene->display_xref) ? $gene->display_xref->display_id : $gene_name,
		 -info_type              => 'MISC',
		 #Always need linkage annotation here
		 #Otherwise ther is no way to trace how these object_xrefs were geneerated
		 -linkage_annotation     => 'eQTL Target',#better text here?
		 -info_text              => 'GENE',
		 #store method param
		 feature_type            => 'AnnotatedFeature',
		);
	  
	  
	  #need sub cache_dbentry_params?
	  push @{$self->{'_dbentry_params'}}, \%gene_xref_params;
	  
	}
	  
	#Now deal with SNP
	  
	#we need to test whether we are on the same slice
	#Need to log if SNP and chrs are on different slices
	#Need to build new features if SNP are on different slice
	if($gene_chr ne $snp_chr){
	  $self->log("Identified Trans-Chromosomal eQTL $gene_name($gene_chr):$rsid($snp_chr)");
	  $self->count('trans-chromosomal');
	}

#Hs.390503

	#Validate SNP first
	my $snp = $self->variation_adaptor->fetch_by_name($rsid);

	if(! defined $snp){
	  $self->log("Failed to retrieve $rsid SNP associated with $gene_name");
	  $self->count('unknown_snps');
	}
	else{

	  #This returns a feat on a core slice, not a genotyped slice
	  my @snp_feats = @{$self->variation_feature_adaptor->fetch_all_by_Variation($snp)};

	  #if (scalar(@snp_feats > 1)){
	  #These are most likely from different haplotype slices
	  #we just need to make sure one of these features loci validates, and then use the 
	  #loci from the file

	  
	  if(scalar(@snp_feats) == 0){
		$self->log("Failed to retrieve $rsid SNP features");
		#As we don't want a region with no xrefs and just a gene name we don't have in the db.
		#Can we store these as unmapped objects?
		$self->count('absent_snps');
	  }
	  else{
		my $snp_feat;

		foreach my $sfeat(@snp_feats){
		  
		  if($sfeat->seq_region_name eq $snp_chr && $sfeat->start == $snp_start){
			$snp_feat = $sfeat;
			last;
		  }
		  #else{
		#	warn "Ensembl $rsid ".$sfeat->feature_Slice->name;
		#  }
		}
	  
		  
	  	  
		if(! defined $snp_feat){
		  $self->log("Failed to retrieve loci validated $rsid SNP feature.");
		  #As we don't want a region with no xrefs and just a gene name we don't have in the db.
		  #Can we store these as unmapped objects?
		  $self->count('invalid_loci_snps');
		}
		else{
	 			
		  #Now test the new and last SNP chrs match
		  if($self->last_snp_chr){
			
			if($self->last_snp_chr ne $snp_chr){
			  throw("Found Multi-Chromosomal eQTL $gene_name:".$self->last_snp_chr.
					":$snp_chr\nNeed to implement multiple features")
			}
		  }
		  else{
			#Set intial feature loci params
			$self->feature_params->{'-slice'} = $snp_feat->slice;
			$self->feature_params->{'-start'} = $snp_start;
			$self->feature_params->{'-end'}   = $snp_start;
		
			#Now we need to extend the feature to the gene loci(SNP may be internal)
			#Either from file or from Ensembl gene
			#This approach is debatable
			#There may be some trans effects going on
			#Hence LD might not be a string factor?

			my $start = $gene_start;
			my $end   = $gene_start;

			if($gene){
			  $start = $gene->start;
			  $end   = $gene->end;
			}

			$self->feature_params->{'-start'} = $start if $start < $self->feature_params->{'-start'};
			$self->feature_params->{'-end'}   = $end   if $end   > $self->feature_params->{'-end'};
		  }
		  
		  $self->last_snp_chr($snp_chr);
		  

		  #Modify the feature params
		  #Do we want it to start directly on the SNP, or include some flank?
		  $self->feature_params->{'-start'} = $snp_start if($snp_start < $self->feature_params->{'-start'});
		  $self->feature_params->{'-end'}   = $snp_start if($snp_start > $self->feature_params->{'-end'});
		  $self->feature_params->{'-score'}++;
		  
		  #Create the SNP xref
		  
		  my %snp_xref_params = (
								 -dbname                 => $self->{'species'}.'_variation_Variation',
								 #-release                => $self->db->dnadb->dbc->dbname,
								 -status                 => 'KNOWNXREF',
								 #-display_label_linkable => 1,
								 #-db_display_name        => $self->db->dnadb->dbc->dbname,
								 -db_display_name        => 'EnsemblVariation',
								 -type                   => 'MISC',#???#this is external_db.type
								 -primary_id             => $snp_feat->variation_name,##?????????This should be the ENSVAR?
								 #Is this a one to one mapping? Do we have ENSVARS for SNP we don't call?
								 -display_id             => $snp_feat->display_id,#This is the same as variation name.
								 -info_type              => 'MISC',
								 -linkage_annotation     => 'eQTL SNP',#DIRECT?	
								 -info_text              => 'VARIATION',
								 #store method param
								 feature_type            => 'AnnotatedFeature',
								);
		
		  #need sub cache_dbentry_params?
		  push @{$self->{'_dbentry_params'}}, \%snp_xref_params;
		}
	  }
	}
  }
  return;
}
 

#This is to allow post-processing of features
#in particular last feature
#And any other format specific end of record operations

sub process_params{
  my $self = shift;
  
  $self->count('Total features');

  #Now we have the possiblity of a sparse feature(no gene identified)
  #with no SNPs retrieved
  if(! defined $self->feature_params->{'-slice'}){
	#Could really do with access to gene_name and rsids here
	#display
	$self->count('failed_eQTLs');
	$self->log('WARNING: Could not build eQTL for '. $self->feature_params->{'-display_label'});

	#Still need to clean params
	#Clean data cache
	$self->{'_feature_params'} = {};
	$self->{'_dbentry_params'} = [];
  }
  else{

	#count eQTL which only have one SNP and no valid gene i.e. 1bp eQTL
	$self->count('1 SNP no gene eQTLs') if scalar(@{$self->dbentry_params}) == 1;


	#Can we genericise this and pass to the load method?
	#We would need to account for xreffing to other eFG feats apart from the one being generated
	#i.e. to reg feats
	#just do here for now

	my $slice = $self->slice_adaptor->fetch_by_region
	  (
	   $self->feature_params->{-slice}->coord_system_name,
	   $self->feature_params->{-slice}->seq_region_name,
	   $self->feature_params->{-start},
	   $self->feature_params->{-end}
	  );

	my @reg_feats = @{$self->regulatory_feature_adaptor->fetch_all_by_Slice($slice)};
	

	$self->count('Total features with RegulatoryFeature xrefs') if (scalar(@reg_feats)>0);

	foreach my $reg_feat(@reg_feats){
	  
	  #$self->count($reg_feat->feature_type->name);
	  #$self->count('Total RegulatoryFeature xrefs');

	  if($self->{'store_reg_feat_xrefs'}){
		my $first_dbentry = $self->dbentry_params->[0];

		#test whether it is the gene xref
		if($first_dbentry->{-dbname} eq 'ensembl_core_Gene'){
		  #defined when assigning

		  $self->count($reg_feat->feature_type->name);
		  $self->count('Total RegulatoryFeature xrefs');
		  
		  my $dbentry = Bio::EnsEMBL::DBEntry->new(%{$first_dbentry});
		  $self->dbentry_adaptor->store($dbentry, $reg_feat->dbID, 'RegulatoryFeature', 1);#1 is ignore release flag
		  #We maybe actually want to use release flag here?
		}
	  }
	}

	my $feature = $self->load_feature_and_xrefs;
	#Not actually using $feature anymore
  }

  $self->{'last_snp_chr'} = undef;
  

}

sub last_snp_chr{
  my ($self, $chr) = @_;

  $self->{'last_snp_chr'} = $chr if defined $chr;

  return $self->{'last_snp_chr'};
}

1;
  
