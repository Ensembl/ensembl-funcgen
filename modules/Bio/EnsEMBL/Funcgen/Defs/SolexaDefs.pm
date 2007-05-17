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

use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Helper;
#use Devel::Size::Report qw(report_size);
#use Devel::Size qw( size total_size);
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
  my ($fh, $file, %roll_back, %chip_files);
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $a_adaptor = $self->db->get_ArrayAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();
  my $probe_a = $self->db->get_ProbeAdaptor();
  my $probef_a = $self->db->get_ProbeFeatureAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  my $anal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("Parzen");
  my $fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("VendorMap");
  
  
  #Get file
  if (! @{$self->result_files()}) {
    my $list = "ls ".$self->input_dir().'/'.$self->name().'*.bed';
    my @rfiles = `$list`;
	throw("Found more than one cluster file:\n@rfiles") if (scalar(@rfiles) >1);
	
    $self->result_files(\@rfiles);
  }
  
  if (scalar(@{$self->result_files()}) >1) {
	throw("Found more than one cluster file:\n".
		  join("\n", @{$self->result_files()})."\nSolexaDefs does not yet handle replicates\n".
		  "we need to resolve how we are going handle replicates with random cluster IDs");
	#do we even need to?
  }


  

  foreach $file(@{$self->result_files()}) {
	my ($chip_uid, $array, $array_chip);

	chomp $file;
    ($chip_uid = $file) =~ s/.*\///;
    $chip_uid =~ s/\..*//;

	$self->log("Found SOLEXA results file for $chip_uid:\t$file");
    $chip_files{$chip_uid} = $file;

	#Store/retrieve dummy chips
	warn "SOLEXA chip_uids are not guaranteed to be unique";

	#do we need to force this to use Experiment to?
	my $echip =  $ec_adaptor->fetch_by_unique_id_vendor($chip_uid, 'SOLEXA');
  
	if ($echip) {
	  #this assumes a successful dummy Array import has already been succesful
	  #or should we just check the experiment_id here?

	  $array_chip = $ac_adaptor->fetch_by_dbID($echip->dbID());
	  $array = $a_adaptor->fetch_by_dbID($array_chip->dbID());

	  
      if (! $self->recovery()) {
		throw("ExperimentalChip(".$echip->unqiue_id().") already exists in the database\nMaybe you want to recover?");
      } else {#log pre-reg'd chips for rollback
		$roll_back{$echip->dbID()} = 1;
	  }
    } else {
	  #we need to make this handle replicates
	  #this needs to be done for each replicate, otherwise we're going to run into problems with resolving the probe cache

	  #create dummy Array and ArrayChip
	  $array = Bio::EnsEMBL::Funcgen::Array->new
		(
		 -NAME        => $self->name().":${chip_uid}",
		 -FORMAT      => 'DUMMY',
		 -VENDOR      => uc($self->vendor()),
		 -TYPE        => 'SEQUENCING',
		 -DESCRIPTION => 'Dummy array for sequencing data',
		);
	  
	  ($array) = @{$a_adaptor->store($array)};  
	  
	  $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
															 -NAME      => $array->name(),
															 -DESIGN_ID => $array->name(),
															 -ARRAY_ID  => $array->dbID(),
															);
	  
	  ($array_chip) = @{$ac_adaptor->store($array_chip)};
	  $array->add_ArrayChip($array_chip);
	  $self->add_Array($array);
	  
	  $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
		(
		 -EXPERIMENT_ID  => $self->experiment->dbID(),
		 -ARRAY_CHIP_ID  => $array_chip->dbID(),
		 -UNIQUE_ID      => $chip_uid,
		);
	
	  ($echip) = @{$ec_adaptor->store($echip)};	
	  $self->experiment->add_ExperimentalChip($echip); #if we need a contains method in  here , always add!!
	
	  
	  
	  foreach my $type ('DUMMY_TOTAL', 'DUMMY_EXPERIMENTAL') {
		
		my $channel = $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());
		
		if ($channel) {
		  if (! $self->recovery()) {
			throw("Channel(".$echip->unique_id().":$type) already exists in the database\nMaybe you want to recover?");
		  }
		} else {
		  
		  $channel =  Bio::EnsEMBL::Funcgen::Channel->new
			(
			 -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
			 -TYPE                 => $type,
			);
		  
		  ($channel) = @{$chan_adaptor->store($channel)};
		}
	  }
	}
  
	my $rset = $self->get_import_ResultSet($anal, 'experimental_chip');


	if ($rset) {				#we have some new data

	  foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {

		if ($echip->has_status('IMPORTED_Parzen', $echip)) {
		  $self->log("ExperimentalChip(".$echip->unique_id().") has already been imported");
		} else {
		  
		  my $cc_id = $rset->get_chip_channel_id($echip->dbID());
		  
		  if ($self->recovery() && $roll_back{$echip->dbID()}) {
			$self->log("Rolling back results for ExperimentalChip:\t".$echip->unique_id());
			$self->db->rollback_results($cc_id);
		  }
		  
		  $self->log("Reading SOLEXA cluster file for ".$echip->unique_id().":\t".$chip_files{$echip->unique_id()});
		  my $fh = open_file($chip_files{$echip->unique_id()});
		  my @lines = <$fh>;
		  close($fh);
		  
		  my $rfile_path = $self->get_dir("norm")."/result.Parzen.".$echip->unique_id().".txt";
		  my $rfile = open_file($rfile_path, '>');
		  my $r_string = "";
		  my ($line, $f_out);
		  my $fasta = '';

		  #warn "we need to either dump the pid rather than the dbID or dump the fasta in the DB dir";
		  my $fasta_file = $ENV{'EFG_DATA'}."/fastas/probe.".$array_chip->name().".fasta";

		  if($self->dump_fasta()){
			$self->backup_file($fasta_file);
			$f_out = open_file($fasta_file, '>');
		  }
		  
		  $self->log("Parsing file:\t $rfile_path");


		  foreach my $line (@lines) {
			$line =~ s/\r*\n//o;

			next if $line =~ /^#/;
		  
			my ($chr, $start, $end, $pid, $score) = split/\t/o, $line;
				  
			#change from UCSC to EnsEMBL coords
			$start +=1;
			$end +=1;
			  

			#$ratio = '\N' if $ratio eq 'NA'; #NULL is still useful info to store in result
			#my ($x, $y) = @{$self->get_probe_x_y_by_name($pid)};
		  
			#this is throwing away the encode region which could be used for the probeset/family?	


			my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
													-NAME          => $pid,
													-LENGTH        => ($end - $start),
													-ARRAY         => $array,
													-ARRAY_CHIP_ID => $array_chip->dbID(),
													-CLASS         => 'CLUSTER',
												   );


			($probe) = @{$probe_a->store($probe)};

			if(!  $self->cache_slice($chr)){
			  warn "Skipping ProbeFeature import, cound non standard chromosome: $chr";
			}else{
			
			  my $feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
				(
				 -START         => $start,
				 -END           => $end,
				 -STRAND        => 1,
				 -SLICE         => $self->cache_slice($chr),
				 -ANALYSIS      => $fanal,
				 -MISMATCHCOUNT => 0,
				 -PROBE         => $probe,
				);

			  $probef_a->store($feature);
			  

			  #dump fasta here

			  if ($self->dump_fasta()){
				$fasta .= '>'.$pid."\n".$self->cache_slice($chr)->sub_Slice($start, $end, 1)->seq()."\n";
			  }

			}




			$r_string .= '\N'."\t".$probe->dbID()."\t${score}\t${cc_id}\t".'\N'."\t".'\N'."\n";#${x}\t${y}\n";
		  }
				
		  print $rfile $r_string;
		  close($rfile);

		  if ($self->dump_fasta()){
			print $f_out $fasta;
			close($f_out);
		  }


		  $self->log("Importing:\t$rfile_path");
		  $self->db->load_table_data("result",  $rfile_path);
		  $self->log("Finished importing:\t$rfile_path");
		  $echip->adaptor->set_status('IMPORTED_Parzen', $echip);


		  #get probe cache from DB for completness
		  #used in probe mapping
		  $self->get_probe_cache_by_Array($array, 1) || throw('Failed to build probe cache');
		}
	  }
	} else {
	  $self->log("No new data, skipping result parse");
	}
  }

  $self->log("Finished parsing and importing results");
  
  return;
}
  


1;
