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

package Bio::EnsEMBL::Funcgen::Parsers::biotiffin;

use strict;

use File::Basename;

# To get files for bioTIFFIN, download the following GFF file (e.g. via wget):
#
#  http://td-blade.gurdon.cam.ac.uk/tad26/fly-tiffinScan-tiffin12.dm3.gff.gz

# Thomas Down <thomas.down@gurdon.cam.ac.uk>

# 
# 3R      MotifScanner    TIFDMEM0000001  936391  936401  0.0     +       0
# 3R      MotifScanner    TIFDMEM0000001  13455911        13455921        0.0     -       0
# 3R      MotifScanner    TIFDMEM0000001  17062830        17062840        0.0     +       0

use Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::ExternalFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Parsers::BaseExternalParser);


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  #Set default feature_type and feature_set config

  #We need to capture version/release/data of external feature sets.
  #This can be nested in the description?  Need to add description to feature_set?

  $self->{'feature_types'} = {
			'BioTIFFIN Motif'   => {
					name        => 'BioTIFFIN Motif',
					class       => 'Regulatory Motif',
					description => 'BioTIFFIN motif',
			}
	};
  
  $self->{feature_sets} = {
			'BioTIFFIN Motif' => {
					feature_type      => \$self->{'feature_types'}{'BioTIFFIN Motif'},
					analysis          => 
					{ 
							-logic_name    => 'BioTIFFIN Motif',
							-description   => 'BioTIFFIN regulatory motif database',
							-display_label => 'BioTIFFIN motifs',
							-displayable   => 1,
					},
					xrefs => 0,
			}
	};
 

  #Move xref flag here?
  $self->{config} =  {
					  'BioTIFFIN Motif' => {
										file  => $ENV{'EFG_DATA'}.'/input/BioTIFFIN/fly-tiffinScan-tiffin12.dm3.gff',
										gff_attrs => {
													  'ID' => 1,
													 },
									   },					  
					 };
  
  
  #Default feature_set names
  if(! defined $self->import_sets){
	@{$self->{'import_sets'}} = keys %{$self->{'feature_sets'}};
  }
  else{#validate

	foreach my $import_fset(@{$self->import_sets}){
	  
	  if(! exists $self->{'feature_sets'}{$import_fset}){
		throw("$import_fset is not a valid import feature set. Maybe you need to add this to the config in:\t".ref($self));
	  }
	}
  }


  #Need to change this so we can just (re)load the one set.
  
  #Change this so we only call it from parse_and_load?
  #Should we validate all first, so we fail at the earliest possible moment?
  #Or serially?

  $self->validate_and_store_feature_types;
  $self->set_feature_sets;

  return $self;
}



# Parse file and return hashref containing:
#
# - arrayref of features
# - arrayref of factors




sub parse_and_load {
  my $self = shift;

  my ($fset_name, $old_assembly, $new_assembly, $file) = rearrange(['FEATURE_SET', 'OLD_ASSEMBLY', 'NEW_ASSEMBLY', 'FILE'], @_);

  warn "file arg not yet fully supported, loading defaults import sets";

  my %slice_cache;
  my $extf_adaptor  = $self->db->get_ExternalFeatureAdaptor;
  my $dbentry_adaptor = $self->db->get_DBEntryAdaptor;
  my $ftype_adaptor = $self->db->get_FeatureTypeAdaptor;
  # this object is only used for projection
  my $dummy_analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'BioTIFFINProjection');#do we need this?
  my $species = $self->db->species;
  if(! $species){
	throw('Must define a species to define the external_db');
  }
  #Just to make sure we hav homo_sapiens and not Homo Sapiens
  ($species = lc($species)) =~ s/ /_/;



  foreach my $import_set(@{$self->import_sets}){
	$self->log_header("Parsing $import_set data");

	my %motif_cache; # name -> factor_id
	my $config = $self->{'config'}{$import_set};
	my $fset =  $self->{'feature_sets'}{$import_set};
	my %gff_attrs =  %{$config->{'gff_attrs'}};
	
	
	# Parse motifs.txt file
	my $file =  $config->{'file'};
	my $skipped = 0;
	my $motif_cnt = 0;
	my $factor_xref_cnt = 0;
	my $feature_cnt = 0;
	my $feature_target_cnt = 0;
	
	open (FILE, "<$file") || die("Can't open $file\n$!\n");

	LINE: while (my $line = <FILE>) {
	  chomp $line;

	  #GFF3
		#3R      MotifScanner    TIFDMEM0000001  936391  936401  0.0     +       0
		#3R      MotifScanner    TIFDMEM0000001  13455911        13455921        0.0     -       0
    #3R      MotifScanner    TIFDMEM0000001  17062830        17062840        0.0     +       0
    #3R      MotifScanner    TIFDMEM0000001  17973965        17973975        0.0     +       0

	  #seq_name, source, feature, start, end, score, strand, frame, [attrs]
	  my ($chromosome, $program, $feature, $start, $end, $score, $strand, undef) = split /\t/o, $line;

	  if(! exists $slice_cache{$chromosome}){
				
				if($old_assembly){
						$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', 
																																							$chromosome, 
																																							undef, 
																																							undef, 
																																							undef, 
																																							$old_assembly);
				} else {
						$slice_cache{$chromosome} = $self->slice_adaptor->fetch_by_region('chromosome', $chromosome);
				}
	  }
		
	  if(! defined  $slice_cache{$chromosome}){
				warn "Can't get slice $chromosome for motif $feature;\n";
				$skipped++;
				next;
	  }
				
		if(! exists $motif_cache{$feature}){
				
				$motif_cache{$feature} = $ftype_adaptor->fetch_by_name($feature);
				
				if(! defined $motif_cache{$feature}){
						
						($motif_cache{$feature}) = @{$ftype_adaptor->store(Bio::EnsEMBL::Funcgen::FeatureType->new
																															 (
																																-name  => $feature,
																																-class => $fset->feature_type->class,
																																-description => $fset->feature_type->description,
																															 ))};						
						
						$motif_cnt ++;
				}
	  }
		
		my $feature_type = $motif_cache{$feature};
		
	  #Now build actual feature

	  $feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new
				(
				 -display_label => $feature,
				 -start         => $start,
				 -end           => $end,
				 -strand        => (($strand eq '+') ? 1 : -1),
				 -feature_type  => $feature_type,
				 -feature_set   => $fset,
				 -slice         => $slice_cache{$chromosome},
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
	  $feature_cnt++;


	}
	
	close FILE;
	
	$self->log("Loaded ".$fset->name);
	$self->log("$motif_cnt feature types");
	$self->log("$feature_cnt features");
	$self->log("Skipped $skipped features");

  }

  return;
}

1;
