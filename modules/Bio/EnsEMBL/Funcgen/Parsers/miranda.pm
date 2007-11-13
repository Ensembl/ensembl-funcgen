package Bio::EnsEMBL::Funcgen::Parsers::miranda;

use strict;

# Parse data from miRanda analyses; format:
#  #<GROUP>	<SEQ>	<METHOD>	<FEATURE>	<CHR>	<START>	<END>	<STRAND>	<PHASE>	<SCORE>	
#  Similarity	hsa-miR-23b	miRanda	miRNA_target	1	919788	919807	+	.	69	transcript id "ENST00000310998"
#  Similarity	hsa-miR-23a	miRanda	miRNA_target	1	919787	919807	+	.	71	transcript id "ENST00000310998"

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
							  'miRanda'   => {
														   class       => 'RNA',
														   description => 'miRanda microRNA',
														  },
							 };
  $self->{feature_sets} = {
						   'miRanda miRNA' => {
											   feature_type      => \$self->{'feature_types'}{'cisRED Search Region'},
											   analysis          => 
											   { 
												-logic_name    => 'miRanda',
												-description   => 'miRanda microRNA target prediction (http://cbio.mskcc.org/mirnaviewer/)',
												-display_label => 'miRanda',
												-displayable   => 1,
											   },
											  },						   
						  };

 
 
  $self->validate_and_store_feature_types;
  $self->set_feature_sets;

  return $self;
}




sub parse_and_load{

  my ($self, $file, $old_assembly, $new_assembly) = @_;

  print ":: Parsing miRanda data from:\t$file\n";

  my $analysis_adaptor = $db_adaptor->get_AnalysisAdaptor();
  my %features_by_name; # name -> feature_type
  my %slice_cache;
  my $display_name_cache = $self->build_display_name_cache('gene');
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'miRandaProjection');
  my $skipped = 0;


  open (FILE, "<$file") || die "Can't open $file";

  while (<FILE>) {
	next if ($_ =~ /^\s*\#/ || $_ =~ /^\s*$/);

    my ($group, $seq, $method, $feature, $chr, $start, $end, $str, $phase, $score, $pvalue, $type, $id_ignore, $id) = split;
    my $strand = ($str =~ /\+/ ? 1 : -1);
    $id =~ s/[\"\']//g;  # strip quotes
	$id .= ':'.$seq;

	if(! exists $slice_cache{$chromosome}){
	
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

	  if(! defined 	$slice_cache{$chr}){
		warn "Can't get slice $chr for sequence $id\n";
		$skipped++;
		next;
	  }
	}





    # ----------------------------------------
    # Feature name

    # For miRNA_target, individual features don't have a unique name, so create
    # a composite one. Also set influence.

  
    #$feature{INFLUENCE} = "negative";#???????????????????????????

    # ----------------------------------------
    # Factor

	if(! exists $features_by_name{$seq}){
	  $features_by_name{$seq} = $ftype_adaptor->fetch_by_name($seq);
	  
	  if(! defined $features_by_name{$seq}){
		($features_by_name{$seq}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
															 (
															  -name  => $seq,
															  -class => 'RNA',
															  -description => 'miRanda RNA',
															 ))};
	  }
	}

   
	my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -display_label => $seq,
	   -start         => $start,
	   -end           => $end,
	   -strand        => $str,
	   -feature_type  => $features_by_name{$seq}
	   -feature_set   => $self->{'feature_sets'}{'miRanda miRNA'},
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

 
    # ----------------------------------------
    # Ensembl object

    my ($ensembl_type) = $type =~ /(gene|transcript|translation)/;
	
    if (!$ensembl_type) {
      print STDERR "Can't get ensembl type from $type, skipping\n";
      next;
    }

    my $ensembl_id = $self->get_display_name_by_stable_id($id, $ensembl_type);
    $ensembl_type = ucfirst(lc($ensembl_type));

    if (!$ensembl_id) {
      print STDERR "Can't get ensembl internal ID for $id, skipping\n";
      next;
    }

    $feature{ENSEMBL_TYPE} = $ensembl_type;
    $feature{ENSEMBL_ID} = $ensembl_id;

    # ----------------------------------------
    # Feature internal ID
    # note this is not the "id" referred to above

    $feature{INTERNAL_ID} = $feature_internal_id++;

    # ----------------------------------------
    # Evidence

    $feature{EVIDENCE} = "";

    # ----------------------------------------
    # Add to object to be returned

    push @features, \%feature;

    #print "Feature: ";
    #foreach my $field (keys %feature) {
    #  print $field . ": " . $feature{$field} . " ";
    #}
    #print "\n";

  }

  close FILE;

  $result{FEATURES} = \@features;
  $result{FACTORS} = \@factors;

  print "Parsed " . scalar(@{$result{FEATURES}}) . " features and " . scalar(@{$result{FACTORS}}) . " factors\n";

  return \%result;

}



sub new {

  my $self = {};
  bless $self, "RegulatoryFeatureParser::miranda";
  return $self;

}

1;
