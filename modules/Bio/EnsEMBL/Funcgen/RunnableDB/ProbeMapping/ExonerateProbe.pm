=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe - 

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateProbe->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, to
align probes to a target genome. (The resulting alignments will be stored in an 
ansembl Funcgen db as Bio::EnsEMBL::ProbeFeature objects.)

NOTE: the ProbeFeature objects refer to Probe id's, and they in turn 
refer to ArrayChip and Array id's. Hence, Arrays, ArrayChips and Probes
should be pre-loaded into the ensembl db: there are separate RunnableDB
/RunnableDB's to do this from the Affymetrix data sets r use the EFG Importer 
to load other arrays e.g. Nimblegen or Sanger. 
This runnable just creates fake Probes in order to create reasonable-looking
affy features???????

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ExonerateProbe;

use warnings ;
use strict;

# use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::BaseExonerate;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::BaseExonerate);

sub new {
  my ( $class, @args ) = @_;
  #my $self = $class->SUPER::new(@args);

  #slightly rearranged order as we want to pass some different defaults to BaseExonerate

  my ($max_mismatches, $mapping_type, $basic_options, $probe_align_helper) = rearrange(['max_mismatches', 'mapping_type', 'basic_options', 'probe_align_helper'], @args);
  
  my %basic_opts;

  #We could change mapping_type to same_strand hits filter
  #As this is the only thing we're using it for here


  if(! defined $basic_options){
	#parse result depends on the output format options
	#only override if you intend overload or rewrite the parse_results method.
	#Now let's reset the default BaseExonerate options to remove vulgar and add scores
	#RESULT: 3020922 0 50 + ENSMUST00000111559 964 1014 + 250 100.00 50 3184 0 scores:0:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0:

	$basic_opts{'-basic_options'} = "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n\" ";
  }

  my $self = $class->SUPER::new(@args, %basic_opts);

  if(! defined $max_mismatches){
	throw("Must provide a -max_mismatches parameter e.g. 0, 1 or 2");
  }

  if(! (defined $mapping_type && $mapping_type =~ /genomic|transcript/)){
	throw("Must provide a valid -mapping_type parameter e.g. genomic or transcript");
  }

  #reset to the new default exonerate
  if (not $self->program) {
	#This is an architecture specific build path
    $self->program('/software/ensembl/compara/exonerate/exonerate');
  }
 


  #if(! (defined $filter_method && (ref($filter_method) ne 'CODE' || $self->can($filter_method)))){
#	throw('You must pass a -filter_method name or CODEREF to filter the ProbeFeatures');
#  }

  ##Set code ref or pointer to internal method
  #We will always have explicitly pass self to the coderef
  

  #$self->{'filter_method'} = (ref($filter_method) eq 'CODE') ? $filter_method : $self->can($filter_method);
  $self->{'max_mismatches'} = $max_mismatches;
  $self->{'mapping_type'}   = $mapping_type;
  
  $self->probe_align_helper($probe_align_helper);

  return $self;
}


sub max_mismatches{
  return $_[0]->{max_mismatches};
}

sub mapping_type{
  return $_[0]->{mapping_type};
}

sub probe_align_helper {
  my ($self, $probe_align_helper) = @_;
  if ($probe_align_helper) {
    confess('Type error') unless ($probe_align_helper->isa('Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlignHelper'));
    $self->{'probe_align_helper'} = $probe_align_helper;
  }
  return $self->{'probe_align_helper'};
}

#
# Implementation of method in abstract superclass
#
sub parse_results {
  my ( $self, $fh ) = @_;
  
  my @features;
  my @promiscuous_probes;
  #my $filter_method = $self->{'match_rules'};
  #my $filter_method = $self->filter_method;
  #No, now uses code ref to allow definition in config

  #print "Parsing results from fh ".Data::Dumper::Dumper($fh)."\n";

  
  my ($tag, $probe_id, $q_start, $q_end, $q_strand, 
	 $t_id, $t_start, $t_end, $t_strand, $score, $tscore,
	 $perc_id, $q_length, $t_length, $mismatch_count, $scores,
	$match_length, $align_mismatch, $total_mismatches, $tmp);

  my $max_mismatches = $self->max_mismatches;
  my $mapping_type   = $self->mapping_type;
  $self->{counts}{total_probe_alignments} =0;
  $self->{counts}{total_probe_mismatches} =0;
  
  #$self->probe_align_helper->dbc->disconnect_when_inactive(0);

  # Memory saving strategy, only put into global @features, if it is not a 
  # promiscuous match.
  #
  my @current_features;
  
  my $current_query;
  my $current_query_num_hits = 0;
  my $previous_query;
  my $query_has_changed_this_iteration;
  
  while (<$fh>){
    #print STDERR $_ if $self->_verbose;
	#warn "\n".$_;
	
    next unless /^RESULT:/;
    chomp;
    
	#Vulgar blocks are also report in in-between coords! So need to add 1(only to start???)
	#Shows  the  alignments  in "vulgar" format.  Vulgar is Verbose Useful Labelled Gapped Alignment Report, This format also starts with the same 9 fields as sugar output (see above), and is followed by a series of <label, query_length, target_length> triplets.  The label may be one of the following:
	
#              M      Match
#              G      Gap
#              N      Non-equivalenced region
#              5      5' splice site
#              3      3' splice site
#              I      Intron
#              S      Split codon
#              F      Frameshift

	#Should always be M probe_length target_length(which would be equal to probe length for ungapped alignment))

	#Can we have a way of dynamically linking the ryo string to a particular method?
	#We shouldn't need different formats!!!
	#RESULT: %S %pi %ql %tl %em scores:{%Ps:}\n
	#RESULT: 3020922 0 50 + ENSMUST00000111559 964 1014 + 250 100.00 50 3184 0 scores:0:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0:
	
	#So we use %em to figure out if we don't have a perfect match, then we use %em the info string to figure out where the mismatches are
	#< 25M effectively means mismatches at the start of end of a query sequence
	
	
    (undef, $probe_id, $q_start, $q_end, $q_strand, 
	 $t_id, $t_start, $t_end, $t_strand, $score, 
	 $perc_id, $q_length, $t_length, $mismatch_count, $scores) = split;

    #print "Current probe seq id: $probe_id current_query = $current_query\n";
	 
    if (! defined $current_query) {
      $current_query = $probe_id;
    }
    if ($current_query ne $probe_id) {
    
      $previous_query = $current_query;
      $previous_query = $probe_id;
      $query_has_changed_this_iteration = 1;
      $current_query = $probe_id;
      
      #print "query_has_changed_this_iteration!\n";
    
    } else {
    
      $query_has_changed_this_iteration = undef;
      $current_query_num_hits++;
    
    }
    use Data::Dumper;    
    
    if ($query_has_changed_this_iteration) {
    
      my $sequence_matches_promiscuously = $current_query_num_hits > 100;
      
      if ($sequence_matches_promiscuously) {
      
	my @probe_analysis_pair = $self->probe_align_helper->fetch_probe_and_analysis_for_probe_feature_linking(
	  $probe_id,
	  $self->mapping_type
	);
	
# 	print "Promiscuous\n";
# 	$Data::Dumper::Maxdepth = 2;
# 	print Dumper(\@probe_analysis_pair);
	
	push @promiscuous_probes, @probe_analysis_pair;
      
      } else {
      
	#print "Not promiscuous\n";
      
 	push @features, @current_features;
# 	$Data::Dumper::Maxdepth = 1;
# 	print Dumper(\@current_features);
      
      }
      # Clear for features of the next probe
      @current_features = ();
      $current_query_num_hits = 0;
    }
    
	
	#Account for half open coords	

	if($t_strand eq '+'){
	  $t_start += 1;	
	}
	else{
	  #Is actually ensembl start
	  $t_end +=1;
	}
	
	if(!($probe_id =~ /\d+/)){
	  throw "Probe headers MUST be the internal db ids of the Probes for this parser to work!\n";
	}
  
	if($mapping_type eq 'transcript' && $t_strand eq '-'){

	  if($q_strand eq '-'){
		#Should never happen
		throw("Found -ve to -ve mapping for probe $probe_id, need to account for this");
	  }
	  else{
		#count here?
		next;
	  }
	}

	$self->{counts}{total_probe_alignments}++;
	$self->{counts}{total_probes}{$probe_id} ||=0;
	$self->{counts}{total_probes}{$probe_id}++;
	
	$match_length = $q_end - $q_start;
	$align_mismatch = $q_length - $match_length;

	$total_mismatches = $mismatch_count + $align_mismatch;
	my @soft_cigar_line;
	
	#5' unaligned
	#1 25
	if($align_mismatch && 
	  ($q_length == $q_end)){
	  push @soft_cigar_line, $q_start.'X' if $q_start;#set this to the value of start if not 0
	  #We want to subtract from start if +ve hit
	  #else we want to add to end if -ve strand, end is actually start in ensembl terms

	  #warn "$t_strand  + ($t_start - $q_start) : - ($t_start + $q_start)";
	  $t_start = ($t_strand eq '+') ? ($t_start - $q_start) : ($t_start + $q_start);
	}

	#mismatches
	#if($mismatch_count){
	if($total_mismatches){
	  my @scores = split/:/, $scores;
	  #remove scores:0
	  shift @scores;
	  shift @scores;
	  #remove last :0
	  pop @scores;
	  #Is this quicker that a splice?

	  my $cnt        = 0;
	  my $prev_score = 5;#always starts with a match otherwise it wouldn't have been reported.

	  while($tscore = shift @scores){
		
		if($tscore == $prev_score){
		  $cnt += 1;
		}else{
		  $tmp = ($prev_score == 5) ?  $cnt.'=' :  $cnt.'X';
		  push @soft_cigar_line, $tmp;
		  $cnt = 1;
		  $prev_score       = $tscore;
		}
	  }

	  #handle last align length
	  $tmp = ($prev_score == 5) ? $cnt.'=' : $cnt.'X';
	  push @soft_cigar_line, $tmp;
	}

	#if(! $t_end || $t_end < $t_start){
	#  throw("GRRR");
	#}


	#3' unaligned
	#0 24
	if($align_mismatch != $q_start){
	  #Add end mismatch if
	  #not accounted for by 5' mismatch
	  my $three_mismatch = ($q_length - $q_end);
	  push @soft_cigar_line, $three_mismatch.'X';


	  #warn "+ ($t_end + $three_mismatch) - ($t_end - $three_mismatch)";

	  $t_end = ($t_strand eq '+') ? ($t_end + $three_mismatch) : ($t_end - $three_mismatch);
	}

	
	#warn "after 3' unaligned start end $t_start $t_end";

	if($total_mismatches <= $max_mismatches){
	#if(1){

	  if($q_strand eq '+'){
		if($t_strand eq '+'){
		  $t_strand = 1;
		  #Alter start/end if we have flanking mismatches here or above?
		}
		elsif($t_strand eq '-'){
		  $t_strand = -1;

		  #Do we really want to do this?
		  #Are we interested in -ve strand exon hits?
		  #Isn't -ve strand transcription normally intronic?

		  #Yes but we may have an ST type array so don't discard these away.

		  ($t_start, $t_end) = reverse($t_start, $t_end);

		  #Reverse the cigar line
		  #As match scores are reported wrt to strand of hit
		  #And we want cigar line wrt to +ve strand of hit
		  @soft_cigar_line = reverse(@soft_cigar_line);
		}
		else{
		  throw "Unrecognised target strand symbol: $t_strand\n";
		}
	  }
	  elsif($q_strand eq '-'){

		throw('We have found a -ve query match');
		#Exonerate only reports +ve query strand matches??? Really?
		#Are we even interested in -ve query strand matches?

		#if($t_strand eq '-'){
		#  $t_strand = 1;
		#}elsif($t_strand eq '+'){
		#  $t_strand = -1;
		#}else{
		#  throw "unrecognised target strand symbol: $t_strand\n";
		#}
	  }else{    
		throw "unrecognised query strand symbol: $q_strand\n";
	  }

	  
	  #warn "final $t_start $t_end ".join(':', @soft_cigar_line) || $match_length.'M';
#--------------------------------------------------------------------------------------------------------------

     my @probe_analysis_pair = $self->probe_align_helper->fetch_probe_and_analysis_for_probe_feature_linking(
      $probe_id,
      $self->mapping_type
    );
    foreach my $current_probe_analysis_pair (@probe_analysis_pair) {
    
      my $probe    = $current_probe_analysis_pair->{probe};
      my $analysis = $current_probe_analysis_pair->{analysis};
      
      foreach my $current_analysis (@$analysis) {
      
	use Bio::EnsEMBL::Funcgen::ProbeFeature;
	
	#push @features, Bio::EnsEMBL::Funcgen::ProbeFeature->new(
	push @current_features, Bio::EnsEMBL::Funcgen::ProbeFeature->new(	
		-PROBE         => $probe,
		-MISMATCHCOUNT => $total_mismatches,
		-START         => $t_start,
		-END           => $t_end,
		-STRAND        => $t_strand,
		-ANALYSIS      => $current_analysis,
		-CIGAR_STRING  => join(':', @soft_cigar_line) || $match_length.'=',
		-seqname => $t_id,
	);
      }
    }
#--------------------------------------------------------------------------------------------------------------
	  
# 	  push @features, new Bio::EnsEMBL::Funcgen::ProbeFeature
# 		(
# 		 #-probe => $probe,
# 		 -probe_id => $probe_id,
# 		 -start => $t_start,
# 		 -end => $t_end,
# 		 -strand => $t_strand,
# 		 -mismatchcount => $total_mismatches,
# 		 -cigar_string => join(':', @soft_cigar_line) || $match_length.'=',
# 		 -seqname => $t_id,
# 		);

	  #warn "After strand $t_start $t_end ";
	  # attach the slice name onto the feature: let the runnabledb
	  # sort out whether it's valid.
	  #$feature->seqname($t_id);
	  
	  $self->{counts}{probe_matches}{$probe_id} ||=0;
	  $self->{counts}{probe_matches}{$probe_id}++;

	}
	else{
      #print "Feature from probe :$probe_id does not match well enough $total_mismatches $max_mismatches \n";
	  $self->{counts}{probe_mismatches}{$probe_id} ||=0;
	  $self->{counts}{probe_mismatches}{$probe_id}++;
	  $self->{counts}{total_probe_mismatches}++;
	}

	#Need to replace this with some count we can report at the end of each run.

  }

  #$self->probe_align_helper->dbc->disconnect_when_inactive(1);
  
  print "Total alignments Probe/ProbeFeature:\t".
	(keys (%{$self->{counts}{total_probes}})).'/'.$self->{counts}{total_probe_alignments}."\n";
  
  print "Failed alignments Probe/ProbeFeatures:\t".
	(keys (%{$self->{counts}{probe_mismatches}})).'/'.$self->{counts}{total_probe_mismatches}."\n";
  
  print "Parsed unfiltered Probe/ProbeFeatures:\t".
	(keys (%{$self->{counts}{probe_matches}})).'/'.scalar(@features)."\n";

  return {
    features => \@features, 
    promiscuous_probes => \@promiscuous_probes
  }
  ;
}

#we need the full unmatched query length to calculat the mismatches!
#Can't assume will always be 25 or 50
#This should be calculated based on % id?

#This can now coderef'd in the analysis config hash
#Altho this obfuscates this module slighlty, it gather all element which need to be 
#configured/added into one place, making this truly generic/extendable

sub filter_mismatches{
  
  my ($self, $query_match_length, $q_length, $score) = @_;

  my $mismatch;
  
  #score is + 5 for a match, -4 for mismatch: hence -9 from total score for one mismatch.
  my $full_score     = $q_length * 5;
  # Should add this as a config var?
  
  #my $mismatch_length = $q_length - 1;
  #my $mismatch_score = $mismatch_length * 5;
  
  #Can we not rewrite this to use the allowed mismatches conf var

  #We could just use the mismatch output instead of calculating all this nonsense!
  #Then we just build the soft cigar line from the info string and the match length query start!


  if($query_match_length == $q_length){
	
	if($score == $full_score){
	  $mismatch = 0;
	}
  }

  if(! defined $mismatch){
	
	my $max_mismatch = $self->max_mismatches;
	
	for my $i(1..$max_mismatch){
	  
	  my $mismatch_length = $q_length - $i;
	  my $mismatch_score = $mismatch_length * 5;
	  
	  if($query_match_length == $q_length){
		
		if ($score == ($full_score - ($i*9))) {
		  $mismatch = $i;
		}
	  }
	  elsif($query_match_length == $mismatch_length){
		$mismatch = $i if ($score == $mismatch_score);
	  }
	}
  }

  return $mismatch;
}



#This simply contains/returns the coderef which has been set to an internal method or sub specified in the RunnableDB config
#And passed to ExonerateProbe via the -filter_method param.

sub filter_method{
  my $self = shift;
  return $self->{'filter_method'};
}



1;

