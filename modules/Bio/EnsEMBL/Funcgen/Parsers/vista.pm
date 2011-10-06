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

package Bio::EnsEMBL::Funcgen::Parsers::vista;

use strict;

# Parse data from LBL enhancers, see http://enhancer.lbl.gov/cgi-bin/imagedb.pl?show=1;search.result=yes;form=search;search.form=no;action=search;search.sequence=1
# e.g.
#
# >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
# AACTGAAGGGACCCCGTTAGCATAtaaacaaaaggtggggggtagccccgagcctcttct
# ctgacagccagtggcggcagtgatgaatttgtgaagttatctaattttccactgttttaa
# ttagagacttgggctctgaggcctcgcagctggcttctttgtgctgtattctgttgcctg
# acagag

use Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser);



sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #Set default feature_type and feature_set config
  $self->{'feature_types'} = {(
							   'VISTA Target'   => {
													name        => 'VISTA Target', 
													class       => 'Search Region',
													description => 'VISTA target region',
												   },
							   'VISTA Enhancer' => {
													name        => 'VISTA Enhancer', 
													class       => 'Enhancer',
													description => 'Enhancer identified by positive VISTA assay',
												   },
							   'VISTA Target - Negative' => {
															 name        => 'VISTA Target - Negative', 
															 class => 'Search Region',
															 description => 'Enhancer negative region identified by VISTA assay',
															},
							  )};
  
  $self->{feature_sets} = {
						   'VISTA enhancer set' => {
													feature_type      => \$self->{'feature_types'}{'VISTA Target'},
													display_label     => 'VISTA Enhancers',
													analysis          => 
													{ 
													 -logic_name => 'VISTA',
													 -description   => 'VISTA Enhancer Assay (http://enhancer.lbl.gov/)',
													 -display_label => 'VISTA',
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



sub parse_and_load{
  my ($self, $file, $old_assembly, $new_assembly) = @_;

  my %result;


  #we want to set up the new external_feature set here or in the caller?

  $self->log_header("Parsing and loading LBNL VISTA enhancer data from:\t$file");

  #my $feature_internal_id = ($self->find_max_id("external_feature")) + 1;
  
  #This is now feature_type
  #but we don't want to just import a flat file for this, as we'd have to remove all the previous
  #entries.  Either check each one, either via the API or by dumping to file and checking the sorted file?
  #my $highest_factor_id = ($self->find_max_id("feature_type")) + 1;

  my $extfeat_adaptor = $self->db->get_ExternalFeatureAdaptor;

  my @features;
  my %feature_objects;
  my %anal;

  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'EnhancerProjection');

  my $feature_positive = $self->{'feature_types'}{'VISTA Enhancer'};
  my $feature_negative = $self->{'feature_types'}{'VISTA Target - Negative'};
  my $set              = $self->{'feature_sets'}{'VISTA enhancer set'}; 


  #we should separate this into rna_feature_type for cisRED?
  #The adaptor should then conditionally look up the rna_feature_type table instead of the feature_type table

  # read file

  open (FILE, "<$file") || die "Can't open $file";
  my $cnt = 0;
  my $skipped = 0;


  while (<FILE>) {

    next if ($_ !~ /^>/o); # only read headers

    my %feature;

    # >chr16:84987588-84988227 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
    my ($coords, $element, $posneg, @stuff) = split /\s+\|\s+/o;

    # parse co-ordinates
    my ($chr, $start, $end) = $coords =~ /chr([^:]+):(\d+)-(\d+)/o;

    # parse element name
    my ($element_number) = $element =~ /\s*element\s*(\d+)/o;

    # ----------------------------------------
    # Feature name
    #$feature{NAME} = "LBNL-$element_number";
    # ----------------------------------------
    # Analysis
    #$feature{ANALYSIS_ID} = $posneg eq 'positive' ? $analysis_positive : $analysis_negative;
	#Feature type
	#$feature{FEATURE_TYPE_ID} = $posneg eq 'positive' ? $feature_positive_id : $feature_negative_id;


    # ----------------------------------------
    # Seq_region ID and co-ordinates

    my $chr_slice;

    if ($old_assembly) {
      $chr_slice = $self->slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, undef, $old_assembly);
    } else {
      $chr_slice = $self->slice_adaptor->fetch_by_region('chromosome', $chr);
    }

    if (!$chr_slice) {
      warn "Can't get slice for chromosme $chr\n";
      next;
    }

    my $seq_region_id = $chr_slice->get_seq_region_id;
	throw("Can't get seq_region_id for chromosome $chr") if (!$seq_region_id);


    #$feature{SEQ_REGION_ID} = $seq_region_id;

    # Assume these are all on the positive strand? Is this correct?
    my $strand = 1;

	my $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
	  (
	   -start         => $start,#is this in UCSC coords?
	   -end           => $end,  #is this in UCSC coords?
	   -strand        => $strand,
	   -feature_type  => $posneg eq 'positive' ? $feature_positive : $feature_negative,
	   -slice         => $self->slice_adaptor->fetch_by_region('chromosome', $chr, undef, undef, $strand, $old_assembly),
	   -display_label => "LBNL-$element_number",
	   -feature_set   => $set,
	  );
	

    # project if necessary
    if ($new_assembly) {

      $feature = $self->project_feature($feature, $dummy_analysis, $new_assembly);

	  if(! defined $feature){
		$skipped ++;
		next;
	  }
    }

	$cnt ++;
	$extfeat_adaptor->store($feature);
  }

  close FILE;

  #$result{FEATURES} = \@features;

  $self->log('Parsed '.($cnt+$skipped).' features');
  $self->log("Loaded $cnt features");
  $self->log("Skipped $skipped features");

  #Now set states
  foreach my $status(qw(DISPLAYABLE MART_DISPLAYABLE)){
	$set->adaptor->store_status($status, $set);
  }
  

  return;

}


1;
