#!/usr/bin/env perl

# Copyright [1999-2018] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use strict;
use REST::Client;
use JSON;
use Data::Dumper;
use feature qw(switch);

#----------------
#Check parameter
#----------------
my $url = $ARGV[0];

if (!$url) {
	usage();
	exit;
}


#--------------------------------------------
#Add parameters to get the JSON object format
#--------------------------------------------
$url.='&frame=object&format=json';

#-----------------------------------------------
#Init REST client and array to store hashes of data
#-----------------------------------------------
my $clientSum = REST::Client->new();
my @arrEpi; #array to store the rows


#----------------------
#create columns headers
#----------------------
my %epiTitle = assign_values('project', 'species', 'accession', 'ihec', 'epigenome', 'reference stage', 'ATAC-Seq', 'ChIP-Seq Input', 'DNase-Seq', 
	'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3', 'CTCF', 'Available Experiments', 'other Histone Mod', 'TFs');

#----------------------------
#add the hash to the array
#----------------------------
#push @arrClsEpi, $clsEpi;
push @arrEpi, \%epiTitle;

#-------------
#Call the Api
#-------------
my $headers = {Accept => 'application/json'};
$clientSum->GET($url, $headers);
my $responseSum = decode_json($clientSum->responseContent());

#----------------
#loop over series
#----------------	

foreach my $serie (@{$responseSum->{'@graph'}}){
	my $project = get_project($serie->{'award'});
	
	#species
	my $species = $serie->{'organism'}[0];
	$species = (split '/', $species)[-1];
	my $scSpecie;
	if (uc $species eq 'MOUSE'){
		$scSpecie = 'Mus+musculus'
	}
	if (uc $species eq 'HUMAN'){
		$scSpecie = 'Homo+sapiens'
	}

	my $accession = $serie->{'accession'};
	my $ihec = $serie->{'dbxrefs'}[0];
	$ihec = (split ':', $ihec)[-1];
	my $epigenome = $serie->{'biosample_term_name'}[0];
	my $refStage = '';
	my $atacSeq = 'NO';
	my $chipSeqInput = 'NO';
	my $DNaseSeq='NO';
	my $H3K27ac = 'NO';
	my $H3K27me3 = 'NO';
	my $H3K36me3 = 'NO';
	my $H3K4me1 = 'NO';
	my $H3K4me3 = 'NO';
	my $H3K9me3 = 'NO';
	my $others = '';
	my $availableExp = 0;
	
	#-------------------
	#look for DNase-Seq
	#-------------------
	foreach my $assayTerm (@{$serie->{'assay_term_name'}}){
		if (uc $assayTerm eq 'DNASE-SEQ') {
			$DNaseSeq = 'YES';
		}
		if (uc $assayTerm eq 'ATAC-SEQ') {
			$atacSeq = 'YES';
		}
	}
	
	#---------------------------------------------------------
	#find Histone modification, Chip-Seq Input and count them
	#Also checks for reference satge consistance
	#---------------------------------------------------------
	my $numH3K27ac = 0;
	my $numH3K27me3 = 0;
	my $numH3K36me3 = 0;
	my $numH3K4me1 = 0;
	my $numH3K4me3 = 0;
	my $numH3K9me3 = 0;
	my $numControl = 0;
	
	my ($refStage, $lstTargets) = get_avilable_targets($serie->{'related_datasets'});
	foreach my $target (@{$lstTargets}){
		given (uc $target){
			when ('H3K27AC') {
				$H3K27ac = 'YES';
				if ($numH3K27ac == 0){
					$numH3K27ac=1;
				}
			}
			when ('H3K27ME3') {
				$H3K27me3 = 'YES';
				if ($numH3K27me3 == 0){
					$numH3K27me3=1;
				}
			}
			when ('H3K36ME3') {
				$H3K36me3 = 'YES';
				if ($numH3K36me3 == 0){
					$numH3K36me3=1;
				}
			}
			when ('H3K4ME1') {
				$H3K4me1 = 'YES';
				if ($numH3K4me1 == 0){
					$numH3K4me1=1;
				}
			}
			when ('H3K4ME3') {
				$H3K4me3 = 'YES';
				if ($numH3K4me3 == 0){
					$numH3K4me3=1;
				}
			}
			when ('H3K9ME3') {
				$H3K9me3 = 'YES';
				if ($numH3K9me3 == 0){
					$numH3K9me3=1;
				}
			}	
			when ('CONTROL') {
				$chipSeqInput = 'YES';
				if ($numControl == 0){
					$numControl=1;
				}
			}
			default {
				if ($others ne '') {
					$others .= ', ';
				}
				$others .= $target;
			}
		}
		
	}
	$availableExp = $numH3K27ac + $numH3K27me3 + $numH3K36me3 + $numH3K4me1 + $numH3K4me3 + $numH3K9me3 + $numControl;
	my @stage = (split ',', $refStage);
	my $life_stage = @stage[0];
	my $age_stage = @stage[1];
	$age_stage =~ s/^\s+|\s+$//g;
	my ($otherTFs, $CTCF) = get_TFs($epigenome, $scSpecie, $life_stage, $age_stage);
	
	
	#-------------------------------------
	#load the data in the hash
	#-------------------------------------		
	
	my %epi = assign_values($project, $species, $accession, $ihec, $epigenome, $refStage, $atacSeq, $chipSeqInput, $DNaseSeq, 
	$H3K27ac, $H3K27me3, $H3K36me3, $H3K4me1, $H3K4me3, $H3K9me3, $CTCF, $availableExp, $others, $otherTFs);

	
	#----------------------------
	#add the hash to the array
	#----------------------------	
	push @arrEpi, \%epi;

}

#----------------------------
#print results
#----------------------------	

print_result(\@arrEpi);





sub assign_values{
	my ($project, $species, $accession, $ihec, $epigenome, $refStage, $ATAC_Seq, $ChIP_Seq_Input, 
		$DNase_Seq, $H3K27ac, $H3K27me3, $H3K36me3, $H3K4me1, $H3K4me3, $H3K9me3, $CTCF, $available_Exp, $other_Histone_Mod, $TFs) = @_;

	my %epi;
	$epi{project} = $project;
	$epi{species} = $species;
	$epi{accession} = $accession;
	$epi{ihec} = $ihec;
	$epi{epigenome} = $epigenome;
	$epi{reference_stage} = $refStage;
	$epi{ATAC_Seq} = $ATAC_Seq;
	$epi{ChIP_Seq_Input} = $ChIP_Seq_Input;
	$epi{DNase_Seq} = $DNase_Seq;
	$epi{H3K27ac} = $H3K27ac;
	$epi{H3K27me3} = $H3K27me3;
	$epi{H3K36me3} = $H3K36me3;
	$epi{H3K4me1} = $H3K4me1;
	$epi{H3K4me3} = $H3K4me3;
	$epi{H3K9me3} = $H3K9me3;
	$epi{CTCF} = $CTCF;
	$epi{Available_Experiments} = $available_Exp;
	$epi{other_Histone_Mod} = $other_Histone_Mod;
	$epi{TFs} = $TFs;

	return %epi;

}
sub print_result {
	my $addEpi = shift;
	my @arrEpi = @{$addEpi};

	foreach my $epi (@arrEpi){
		my %REpi = %{$epi};
		my $csv_row = join ("\t", $REpi{project}, $REpi{species}, $REpi{accession}, $REpi{ihec}, $REpi{epigenome}, $REpi{reference_stage}, 
			$REpi{ATAC_Seq}, $REpi{ChIP_Seq_Input}, $REpi{DNase_Seq}, $REpi{H3K27ac}, $REpi{H3K27me3}, $REpi{H3K36me3}, $REpi{H3K4me1}, 
			$REpi{H3K4me3}, $REpi{H3K9me3}, $REpi{CTCF}, $REpi{Available_Experiments}, $REpi{other_Histone_Mod}, $REpi{TFs});
		
		print $csv_row."\n";	
	}

}

sub get_TFs {
	my $epigenome = shift;
	my $specie = shift;
	my $life_stage = shift;
	my $age = shift;

	my $TFurl = 'https://www.encodeproject.org/search/?searchTerm='.$epigenome;
	$TFurl .= '&type=Experiment&assay_slims=DNA+binding';
	$TFurl .= '&replicates.library.biosample.donor.organism.scientific_name='.$specie;
	$TFurl .= '&limit=all';
	$TFurl .= '&target.investigated_as=transcription+factor';
	$TFurl .= '&assay_title=ChIP-seq';
	$TFurl .= '&replicates.library.biosample.life_stage='.$life_stage;
	$TFurl .= '&replicates.library.biosample.age='.$age;
	$TFurl .= '&frame=object&format=json';

	
	#my $TFurl = 'https://www.encodeproject.org/search/?searchTerm='.$epigenome;
	#$TFurl .= '&type=Experiment&assay_slims=DNA+binding&assay_title=ChIP-seq&replicates.library.biosample.donor.';
	#$TFurl .= 'organism.scientific_name='.$specie.'&limit=all&target.investigated_as=transcription+factor&frame=object&format=json';
	my $clientAwd = REST::Client->new();
	my $headersAwd = {Accept => 'application/json'};
	$clientAwd->GET($TFurl, $headersAwd);
	my $responseAwd = decode_json($clientAwd->responseContent());
	my $otherTFs='';
	my $CTCF='NO';
	my @Tfs;
	foreach my $results (@{$responseAwd->{'@graph'}}){
		if (uc $results->{'status'} eq 'RELEASED'){
			my $rawTarget = (split '/', $results->{'target'})[-1];
			my $target = (split '-', $rawTarget)[0];
			if ($target){
				if (not grep(/^$target$/, @Tfs)){
					push @Tfs, $target;
					if ($otherTFs ne '') {
						$otherTFs .= ', ';
					}
					$otherTFs .= $target;

					if ($target eq 'CTCF'){
						$CTCF='YES';
					}
				}
			}
		}
	}

	return $otherTFs, $CTCF;

}

sub get_project {
	my $url = shift;
	my $clientAwd = REST::Client->new();
	my $headersAwd = {Accept => 'application/json'};
	my $baseUrl = 'https://www.encodeproject.org';
	my $awdUrl = $baseUrl.$url.'?frame=object&format=json';
	$clientAwd->GET($awdUrl, $headersAwd);
	my $responseAwd = decode_json($clientAwd->responseContent());
	my $project = "";
	if ($responseAwd->{'project'}){
		$project = $responseAwd->{'project'};
	}
	
	return $project;
}

sub get_avilable_targets {
	my $experiments = shift;
	my @lstExperiments = @{$experiments};
	my $clientExp = REST::Client->new();
	my $headersExp = {Accept => 'application/json'};
	my $baseUrl = 'https://www.encodeproject.org';
	my @lstReleasedExps;
	my $ret='';
	my @lifeStage;
	my @age;
	foreach my $exp (@lstExperiments){
		my $expUrl = $baseUrl.$exp.'?frame=object&format=json';
		$clientExp->GET($expUrl, $headersExp);
		my $responseExp = decode_json($clientExp->responseContent());
		#-----------------------------------------------------
		#filter for ChIP-Seq, DNase-Seq and status = released
		#-----------------------------------------------------		
		if ((uc $responseExp->{'status'} eq 'RELEASED') && ((uc $responseExp->{'assay_title'} eq 'CHIP-SEQ') || (uc $responseExp->{'assay_title'} eq 'DNASE-SEQ'))){
			#check for target of chip-seq assays
			if (uc $responseExp->{'assay_title'} eq 'CHIP-SEQ'){
				my $rawTarget = (split '/', $responseExp->{'target'})[-1];
				my $target = (split '-', $rawTarget)[0];
				if ($target){
					push @lstReleasedExps, $target;
				}
			}
			#-------------------------------------
			#consistency of reference_stage 
			#-------------------------------------			

			my @replicates = @{$responseExp->{'replicates'}};

			foreach my $rep (@replicates) {
				my $urlReplicate = $baseUrl.$rep.'?frame=embedded&format=json';
				my $clientRep = REST::Client->new();
				$clientRep->GET($urlReplicate, $headersExp);
				my $responseRep = decode_json($clientRep->responseContent());
				
				#add LifeStage if is not been already added
				if ( !grep { $_ eq $responseRep->{'library'}->{'biosample'}->{'life_stage'}} @lifeStage )
				{
				  push @lifeStage, $responseRep->{'library'}->{'biosample'}->{'life_stage'};
				}

				#add Age if is not been already added
				if ( !grep { $_ eq $responseRep->{'library'}->{'biosample'}->{'age'}} @age )
				{
				  push @age, $responseRep->{'library'}->{'biosample'}->{'age'};
				}					
			}				

		}		
	}
	
	#-------------------------------------
	#check consistency of reference_stage 
	#-------------------------------------		
	my $numLifeStage = @lifeStage;
	my $numAge = @age;
	if ($numLifeStage > 1 || $numAge > 1){
		$ret='WARNING -- not consistent';
		if ($numLifeStage > 1){
			$ret .= ' -- Life Stages: ';
			foreach my $lStage (@lifeStage){
				$ret .= $lStage.' - '
			}
		}
		if ($numAge > 1){
			$ret .= ' -- Ages: ';
			foreach my $age (@age){
				$ret .= $age.' - '
			}
		}
	}else{
		$ret = @lifeStage[0].', '.@age[0];
	}
	
	return $ret, \@lstReleasedExps;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: encode_reference_epigenomes.pl url [look_for_reference_stage_consistency]

Options:
url: The encode search url

 
END_USAGE

    say $usage;

    return 1;
}


