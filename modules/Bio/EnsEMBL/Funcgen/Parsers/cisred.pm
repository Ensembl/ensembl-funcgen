package Bio::EnsEMBL::Funcgen::Parsers::cisred;

use strict;

use File::Basename;

# To get files for CisRed data, download the following 2 files (e.g. via wget):
#
# http://www.cisred.org/content/databases_methods/human_2/data_files/motifs.txt
#
# http://www.cisred.org/content/databases_methods/human_2/data_files/search_regions.txt


#No longer valid urls, now use the following for ensembl formats for all species:
#http://www.bcgsc.ca/downloads/cisred/temp/cisRED4Ensembl/
#naminf may not be obvious, may have to cross reference with above previous urls to get build info

# Format of motifs.txt (note group_name often blank)

#name    chromosome      start   end     strand  group_name      ensembl_gene
#craHsap1        1       168129978       168129997       -1      1       ENSG00000000457
#craHsap2        1       168129772       168129781       -1      2       ENSG00000000457
#craHsap3        1       168129745       168129756       -1      3       ENSG00000000457
#craHsap4        1       168129746       168129753       -1      4       ENSG00000000457
#craHsap5        1       168129745       168129752       -1      5       ENSG00000000457
#craHsap6        1       168129741       168129757       -1      6       ENSG00000000457


# Format of search_regions.txt
# name	chromosome	start	end	strand	ensembl_gene_id
# 1	17	39822200	39824467	-1	ENSG00000005961
# 8	17	23151483	23153621	-1	ENSG00000007171
# 14	1	166434638	166437230	-1	ENSG00000007908
# 19	1	23602820	23605631	-1	ENSG00000007968


use Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::ExternalFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw );


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser);





sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #Set default feature_type and feature_set config
  $self->{'feature_types'} = {
							   'cisRED Search Region'   => {
															class       => 'Region',
															description => 'cisRED search region',
												   },
							   'cisRED Motif' => {
												  class       => 'Regulatory Motif',
												  description => 'cisRED group motif',
												 },
							 };
  
  $self->{feature_sets} = {
						   'cisRED search regions' => {
													   feature_type      => \$self->{'feature_types'}{'cisRED Search Region'},
													   analysis          => 
													   { 
														-logic_name    => 'cisRED',
														-description   => 'cisRED motif search (www.cisred.org)',
														-display_label => 'cisRED',
														-displayable   => 1,
													   },
													  },
						   'cisRED group motifs' => {
													 feature_type      => \$self->{'feature_types'}{'cisRED Motif'},
													 analysis          => 
													 { 
													  -logic_name    => 'cisRED',
													  -description   => 'cisRED motif search (www.cisred.org)',
													  -display_label => 'cisRED',
													  -displayable   => 1,
													 },
													},
						  };
 
  $self->validate_and_store_feature_types;
  $self->set_feature_sets;

  return $self;
}





# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors

sub parse_and_load {
  my ($self, $file, $old_assembly, $new_assembly) = @_;
  print ":: Parsing cisRED data\n";

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
  my %features_by_group; # name -> factor_id
  my %slice_cache;
  my $extf_adaptor  = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_adaptor = $self->db->get_DBEntryAdaptor;
  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
  #my $display_name_cache = $self->build_display_name_cache('gene');
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'CisRedProjection');

  # ----------------------------------------
  # We need a "blank" factor for those features which aren't assigned factors
  # Done this way to maintain referential integrity
  #my $blank_factor_id = $self->get_blank_factor_id($db_adaptor);

  # Parse motifs.txt file
  print ":: Parsing cisRED motifs from $file\n";
  my $skipped = 0;
  my $skipped_xref = 0;
  #my $coords_changed = 0;
  my $cnt = 0;
	
  open (FILE, "<$file") || die "Can't open $file";
  <FILE>; # skip header

  while (<FILE>) {
    next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);
	chomp;
    # name	chromosome	start	end	strand	group_name   ensembl_gene_id
    my ($motif_name, $chromosome, $start, $end, $strand, $group_name, $gene_id) = split (/\t/);
    #($gene_id) = $gene_id =~ /(ENS.*G\d{11})/;

	if(! exists $slice_cache{$chromosome}){
	
	  if($old_assembly){
		$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', 
																	$chromosome, 
																	undef, 
																	undef, 
																	undef, 
																	$old_assembly);
	  }else{
		$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', $chromosome);
	  }

	  if(! defined 	$slice_cache{$chromosome}){
		warn "Can't get slice $chromosome for motif $motif_name\n";
		$skipped++;
		next;
	  }
	}

	#get feature_type first

	#we are not maintaining this link in the DB!
	#Do we need another xref for this or a different table?

	
	if ($group_name && $group_name ne '' && $group_name !~ /\s/) {

	  if(! exists $features_by_group{$group_name}){
		$features_by_group{$group_name} = $ftype_adaptor->fetch_by_name('crtHsap'.$group_name);

		if(! defined $features_by_group{$group_name}){
		  ($features_by_group{$group_name}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
																	  (
																	   -name  => 'crtHsap'.$group_name,
																	   -class => 'Regulatory Motif',
																	   -description => 'cisRED group motif',
																	  ))};
		}
	  }
	}
	#}else{
	#  throw("Found cisRED feature $motif_name with no group_name, unable to defined feature_type");
	#}


	my $ftype = (defined $features_by_group{$group_name}) ? $features_by_group{$group_name} : $self->{'feature_sets'}{'cisRED group motifs'}->feature_type;

	my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -display_label => $motif_name,
	   -start         => $start,
	   -end           => $end,
	   -strand        => $strand,
	   -feature_type  => $ftype,
	   -feature_set   => $self->{'feature_sets'}{'cisRED group motifs'},
	   -slice         => $slice_cache{$chromosome},
	  );


	
    # project if necessary
    if ($new_assembly) {
      $feature = $self->project_feature($feature, $new_assembly);

	  if(! defined $feature){
		$skipped ++;
		next;
	  }


      #$coords_changed++ if ($projected_feature->start() != $start || $projected_feature->end() != $end);
    }

	($feature) = @{$extf_adaptor->store($feature)};
	$cnt++;


	#need to validate gene_id here
	#retrieve from db?
	#could have version here if we use the correct dnadb to build the cache
	
	#Build Xref here
	if (! $gene_id) {
      warn("No xref available for motif $motif_name\n");
      $skipped_xref++;
      next;
    }

	my $display_name = $self->get_display_name_by_stable_id($gene_id);

	#Handle release/version in xref version as stable_id version?

	my $dbentry = Bio::EnsEMBL::DBEntry->new(
											 -dbname                 => 'core_gene',
											 #-release                => $self->db->dnadb->dbc->dbname,
											 -status                 => 'KNOWNXREF',
											 #-display_label_linkable => 1,
											 -#db_display_name        => $self->db->dnadb->dbc->dbname,
											 -db_display_name        => 'ensembl_core_gene',
											 -type                   => 'MISC',
											 -primary_id             => $gene_id,
											 -display_id             => $display_name,
											 -info_type              => 'MISC',
											 -description            => 'cisRED motif gene xref',
											 #could have version here if we use the correct dnadb to build the cache
											);
	$dbentry_adaptor->store($dbentry, $feature->dbID, 'ExternalFeature', 1);#1 is ignore release flag
  }
  
	
  close FILE;

  print ":: Stored $cnt cisRED ExternalFeature motif\n";
  print ":: Skipped $skipped cisRED ExternalFeature motif imports\n";
  print ":: Skipped an additional $skipped_xref DBEntry imports\n";

  # ----------------------------------------
  # Search regions 
  # read search_regions.txt from same location as $file

  #my $search_regions_file = dirname($file) . "/search_regions.txt";
  my $search_regions_file;
  ($search_regions_file = $file) =~ s/motifs/searchregions/;

  $skipped = 0;
  $cnt = 0;
  $skipped_xref = 0;
	
  print ":: Parsing cisRED search regions from $search_regions_file\n";
  open (SEARCH_REGIONS, "<$search_regions_file") || die "Can't open $search_regions_file";
  <SEARCH_REGIONS>; # skip header

  while (<SEARCH_REGIONS>) {
    chomp;
    my ($id, $chromosome, $start, $end, $strand, $gene_id) = split;
    my $display_id = $self->get_display_name_by_stable_id($gene_id);
	my $name = "CisRed_Search_$id";

	if(! exists $slice_cache{$chromosome}){
	  
	  if($old_assembly){
		$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', 
																		  $chromosome, 
																		  undef, 
																		  undef, 
																		  undef, 
																		  $old_assembly);
	  }else{
		$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', $chromosome);
	  }

	  if(! defined 	$slice_cache{$chromosome}){
		warn "Can't get slice $chromosome for for search region $name\n";
		next;
	  }
	}


	

	my $search_feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -display_label => $name,
	   -start         => $start,
	   -end           => $end,
	   -strand        => $strand,
	   -feature_type  => $self->{'feature_sets'}{'cisRED search regions'}->feature_type,
	   -feature_set   => $self->{'feature_sets'}{'cisRED search regions'},
	   -slice         => $slice_cache{$chromosome},
	  );
																


    # project if necessary
    if ($new_assembly) {
      $search_feature = $self->project_feature($search_feature);
	  
	  if(! defined $search_feature){
		$skipped ++;
		next;
	  }
	}

	$extf_adaptor->store($search_feature);
	$cnt++;

	#Build Xref here
	#need to validate gene_id here!!

	if (! $gene_id) {
      warn("Can't get internal ID for $gene_id\n");
      $skipped_xref++;
      next;
    }

	my $display_name = $self->get_display_name_by_stable_id($gene_id);
	
	my $dbentry = Bio::EnsEMBL::DBEntry->new(
											 -dbname                 => 'core',
											 -release                => $self->db->dnadb->dbc->dbname,
											 -status                 => 'KNOWNXREF',
											 #-display_label_linkable => 1,
											 -db_display_name        => $self->db->dnadb->dbc->dbname,
											 -type                   => 'MISC',
											 -primary_id             => $gene_id,
											 -display_id             => $display_name,
											 -info_type              => 'MISC',
											 -description            => 'cisRED search region gene xref',
											 #could have version here if we use the correct dnadb to build the cache
											  );
	$dbentry_adaptor->store($dbentry, $search_feature->dbID, 'ExternalFeature');  
  }

  close(SEARCH_REGIONS);

  
  print ":: Stored $cnt cisRED search region ExternalFeatures\n";
  print ":: Skipped $skipped cisRED search region ExternalFeatures\n";
  print ":: Skipped an additional $skipped_xref cisRED search region DBEntry imports\n";

  #print "$coords_changed features had their co-ordinates changed as a result of assembly mapping.\n" if ($new_assembly);

  return;

}



1;
