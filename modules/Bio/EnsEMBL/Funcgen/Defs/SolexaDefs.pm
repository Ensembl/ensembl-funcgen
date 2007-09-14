#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Defs::SolexaDefs
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Defs::SolexaDefs

=head1 SYNOPSIS

  my $defs_type = "Bio::EnsEMBL::Funcgen::Defs::SolexaDefs";
  push @INC, $defs_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally set by the Importer as the parent class.  SolexaDefs contains meta 
data and methods specific to Solexa clustered sequencing data, to aid 
parsing and importing of experimental data.

=head1 AUTHOR

This module was written by Nathan Johnson.

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::Defs::SolexaDefs;

use Bio::EnsEMBL::Funcgen::ExperimentalSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Helper;
use strict;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Helper);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for SolexaDefs class
  Returntype : Bio::EnsEMBL::Funcgen::Defs::SolexaDefs
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self  = $class->SUPER::new();

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	
  ($self->{'name'}) = rearrange(['NAME'], @_);

  #throw('Must provide an Experiment name for a Solexa import') if ! defined $name;

  #should we provide args override for all of these?

  $self->{'defs'} =  
    {(
      #order of these data arrays is important!
      array_data   => [],#['experiment'],
      probe_data   => [],#["probe"],
      results_data => ["and_import_bed"],
      #sample_key_fields => ['DESIGN_ID', 'CHIP_ID', 'DYE', 'PROMOT_SAMPLE_TYPE'],# 'SAMPLE_LABEL'],label now optional
      # 'SAMPLE_DESCRIPTION removed due to naming disparities
      #ndf_fields      => ['CONTAINER', 'PROBE_SEQUENCE', 'MISMATCH', 'FEATURE_ID', 'PROBE_ID'],
      #pos_fields      => ['CHROMOSOME', 'PROBE_ID', 'POSITION', 'COUNT'],
      #result_fields   => ['PROBE_ID', 'PM', 'X', 'Y'],
      #notes_fields   => ['DESIGN_ID', 'DESIGN_NAME', 'DESCRIPTION'],
      norm_method => undef,
      

	  #Need to make these definable?
	  #have protocolfile arg and just parse tab2mage protocol section format
#	  protocols => {(
#					 grow          => {(
#										accession => 'GROW_NIMB',
#										name      => 'GROW NIMBLEGEN CULTURE CONDITIONS',
#										text      => 'Nimblegen culture conditions description here. Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										paramters => undef,
#									   )},
#					 treatment     => {(
#										accession => 'CROSSLINK_NIMB',
#										name      => 'NIMBLEGEN CHROMATIN PREPARATION',
#										text      => 'Nimblegen X-linking and DNA extraction protocol.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										paramters => undef,
#									   )},
#					 extraction    => {(
#										accession => 'CHROMATIN_IP_NIMB',
#										name      => 'NIMBLEGEN CHROMATIN IMMUNOPRECIPITATION and DNA RECOVERY',
#										text      => 'Nimblegen chromatin immunoprecipitation and DNA extraction protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										paramters => undef,
#									   )},
#					 labeling      => {(
#										accession => 'LABELLING_NIMB',
#										name      => 'NIMBLEGEN LABELLING',
#										text      => 'Nimblegen labelling protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										paramteres => undef,
#									   )},
#					 hybridization => {(
#										accession => 'HYBRIDISATION_NIMB',
#										name      => 'NIMBLEGEN HYBRIDISATION',
#										text      => 'Nimblegen chip hybridisation protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										parameters => undef,
#									   )},
#					 scanning      => {(
#										accession => 'SCANNING_NIMB',
#										name      => 'NIMBLESCAN',
#										text      => 'Nimblegen Nimblescan protocol here.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.Padding text here to avoid description too short warnings.',
#										paramters => undef,
#									   )},
#					)},

     )};
	
  return $self;
}


=head2 set_defs

  Example    : my $self->set_defs;
  Description: Sets attribute dependent defs
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub set_defs{
  my $self = shift;

#  #dir are not set in defs to enable generic get_dir method access

#  $self->{'design_dir'} = $self->get_dir('data').'/input/'.
  #    $self->vendor().'/'.$self->name().'/DesignFiles';
    
  #$self->{'defs'}{'bed_file'} = $self->get_dir('data').'/input/'.
  #  $self->vendor().'/'.$self->name().'/SampleKey.txt';
    
#  $self->{'defs'}{'notes_file'} = $self->get_dir('data').'/input/'.
#    $self->vendor().'/'.$self->name().'/DesignNotes.txt';
  
#  $self->{'defs'}{'tab2mage_file'} = $self->get_dir('data').'/output/'.
#    $self->vendor().'/'.$self->name().'/E-TABM-'.$self->name().'.txt';

#  $self->{'defs'}{'mage_xml_file'} = $self->get_dir('data').'/output/'.
#    $self->vendor().'/'.$self->name().'/{UNASSIGNED}.xml';

#  $self->{'results_dir'} = $self->get_dir('data').'/input/'.
#    $self->vendor().'/'.$self->name().'/PairData';

  return;
}



 
sub read_and_import_bed_data{
  my $self = shift;
  
  $self->log("Reading and importing ".$self->vendor()." data");
  my (@header, @data, @design_ids, @lines);
  my ($fh, $file);
 
  my $eset_adaptor = $self->db->get_ExperimentalSetAdaptor();
  my $af_adaptor = $self->db->get_AnnotatedFeatureAdaptor();
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("Parzen");
  my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("VendorMap");
  my $new_data = 0;
 
  my $eset = $eset_adaptor->fetch_by_name($self->experimental_set_name());
  
  if(! defined $eset){
	$eset = Bio::EnsEMBL::Funcgen::ExperimentalSet->new(
														-name         => $self->experimental_set_name(),
														-experiment   => $self->experiment(),
														-feature_type => $self->feature_type(),
														-cell_type    => $self->cell_type(),
														-vendor       => $self->vendor(),
														-format       => $self->format(),
													   );
	($eset)  = @{$eset_adaptor->store($eset)};
  }


  #we need a way to define replicates on a file basis when we have no meta file!
  #can we make this generic for application to array imports?
  #currently we have to do a separate import for each replicate, specifying the result files each time
  #we need to add a experimental_set_name option


  #Get file
  if (! @{$self->result_files()}) {
    my $list = "ls ".$self->input_dir().'/'.$self->name().'*.bed';
    my @rfiles = `$list`;
	throw("Found more than one cluster file:\n@rfiles\nNeed to implement ExperimentalSubset rollback before removing this!") if (scalar(@rfiles) >1);
	
    $self->result_files(\@rfiles);
  }
  
  if (scalar(@{$self->result_files()}) >1) {
	throw("Found more than one cluster file:\n".
		  join("\n", @{$self->result_files()})."\nSolexaDefs does not yet handle replicates\n".
		  "we need to resolve how we are going handle replicates with random cluster IDs");
	#do we even need to?
  }


  #how are we going to track import of files if they are being directly imported into annotated_feature?
  #Current solution is to create dummy chips, but we want something neater
  #do we need to track import as closely as with chips?
  
  foreach my $filepath(@{$self->result_files()}) {
	chomp $filepath;
	my $filename;
	my $roll_back = 0;
    ($filename = $filepath) =~ s/.*\///;
	my $sub_set;

	$self->log("Found SOLEXA results file\t$filename");

	if($sub_set = $eset->get_subset_by_name($filename)){
	  $roll_back = 1;
	}else{
	  $sub_set = $eset->add_new_subset($filename);
	}
	
	#store if not already, skips if stored
	$eset_adaptor->store_ExperimentalSubsets([$sub_set]);

	warn "Got stored eset with dbID ".$eset->dbID();
  
	if ($sub_set->adaptor->has_status('IMPORTED', $sub_set)){
	  $self->log("ExperimentalSubset(${filename}) has already been imported");
	} 
	else {
	  $new_data = 1;

	  if ($self->recovery() && $roll_back) {
		$self->log("Rolling back results for ExperimentalSubset:\t".$filename);

		warn "Cannot yet rollback for just an ExperimentalSubset, rolling back entire set\n";
		throw("Need to implement annotated_feature rollback!\n");
		#$self->db->rollback_results($cc_id);
	  }
		  
	  $self->log("Reading SOLEXA cluster file:\t".$filename);
	  my $fh = open_file($filepath);
	  my @lines = <$fh>;
	  close($fh);
		  
	  #my $rfile_path = $self->get_dir("norm")."/result.Parzen.".$echip->unique_id().".txt";
	  #my $rfile = open_file($rfile_path, '>');
	  #my $r_string = "";
	  my ($line, $f_out);
	  my $fasta = '';
	  
	  #warn "we need to either dump the pid rather than the dbID or dump the fasta in the DB dir";
	  my $fasta_file = $ENV{'EFG_DATA'}."/fastas/".$self->experiment->name().'.'.$filename.'.fasta';

	  if($self->dump_fasta()){
		$self->backup_file($fasta_file);
		$f_out = open_file($fasta_file, '>');
	  }
		  
	  $self->log("Parsing file:\t$filename");

	  foreach my $line (@lines) {
		$line =~ s/\r*\n//o;
		next if $line =~ /^#/;
		
		my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;				  
		#change from UCSC to EnsEMBL coords
		$start +=1;
		$end +=1;
		
		if(!  $self->cache_slice($chr)){
		  warn "Skipping AnnotatedFeature import, cound non standard chromosome: $chr";
		}else{
		  
		  #this is throwing away the encode region which could be used for the probeset/family?	
		  my $feature = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new
			(
			 -START         => $start,
			 -END           => $end,
			 -STRAND        => 1,
			 -SLICE         => $self->cache_slice($chr),
			 -ANALYSIS      => $fanal,
			 -DISPLAY_LABEL => $pid,
			);
		  
		  $af_adaptor->store($feature);
		  
		  #dump fasta here
		  if ($self->dump_fasta()){
			$fasta .= '>'.$pid."\n".$self->cache_slice($chr)->sub_Slice($start, $end, 1)->seq()."\n";
		  }
		}
	  }


	  if ($self->dump_fasta()){
		print $f_out $fasta;
		close($f_out);
	  }


	  $self->log("Finished importing:\t$filepath");
	  $sub_set->adaptor->set_status('IMPORTED', $sub_set);
	}
  }

  $self->log("No new data, skipping result parse") if ! $new_data;
  
  $self->log("Finished parsing and importing results");
  
  return;
}
  


1;
