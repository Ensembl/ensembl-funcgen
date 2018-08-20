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
use Bio::EnsEMBL::Funcgen::clsRegisterMetadata;
use REST::Client;
use JSON;
use Data::Dumper;
use Getopt::Long qw(GetOptions);
use feature qw(say);
use Config::Tiny;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBEntry;

##################### MAIN FUNCTION ##############################

#Get Parameters
my $pathFile;
my $localFilePath;
my @lstTargets;
my @lstAssays;
my $help;
my $cfgFile;
my @lstErrors;
my $skipPath;
my %skipRefEpi;
my %featureNotFound;
my %features_done;

my %targetsNotFound;

GetOptions(
			'f=s' => \$pathFile, 
			'p=s' => \$localFilePath, 
			'a=s' => \@lstAssays,
			'c=s' => \$cfgFile,
			't=s' => \@lstTargets,
			's=s' => \$skipPath, 
			'h=s' => \$help,);

#check parameters			
if ( $help || !$pathFile || !$localFilePath || !@lstAssays || !$cfgFile) {
    usage();
    exit;
}
if (not(-f $pathFile)){
	#file not found
	die "$pathFile file can not be found\n";
}

if (not(-f $cfgFile)){
	#file not found
	die "$cfgFile file can not be found\n";
}

if ($skipPath){
	if (not(-f $skipPath)){
		#file not found
		die "$skipPath file can not be found\n";
	}
	#load skip list
	open(my $sfh, '<:', $skipPath) or die "Could not open file $skipPath\n";
	while (my $fullRow = <$sfh>) {
		chomp $fullRow;
		$skipRefEpi{$fullRow} = $fullRow;
	}
	close($sfh);

}




#read config file
my $cfg = Config::Tiny->read($cfgFile);

#get adaptors
my $adaptors = fetch_adaptors($cfg);

#create hashes to compare and normalize the data
my ($hAnalysis, $hFeatureType, $hGender) = get_compare_hashes($adaptors);
my %hshAnalysis = %{$hAnalysis};
my %hshFeatureType = %{$hFeatureType};
my %hshGender = %{$hGender};

#delete 'Features_not_found.txt' if exists
if(-e 'Features_not_found.txt') 
{
   unlink 'Features_not_found.txt';
}
#delete 'Features_not_imported.txt' if exists
if(-e 'Features_not_imported.txt') 
{
   unlink 'Features_not_imported.txt';
}
#delete 'new_Features_to_be_imported.txt' if exists
if(-e 'new_Features_to_be_imported.txt') 
{
   unlink 'new_Features_to_be_imported.txt';
}


#compare parameter values with DB values (assays and targets) to check if they exists and to normalize them
my @lstDbAssays;
foreach my $assayVal (@lstAssays){
	my $dbAssay = check_db_value($assayVal, \%hshAnalysis);
	if (!$dbAssay){
		push @lstErrors, "parameter: assay\tvalue: $assayVal\terror: Not found in analysis table";
	}else {
		push @lstDbAssays, $dbAssay;
	}
}

my @lstDbTargets;
foreach my $targetVal (@lstTargets){
	my $dbTarget = check_db_value($targetVal, \%hshFeatureType);
	if (!$dbTarget){
		push @lstErrors, "parameter: target\tvalue: $targetVal\terror: Not found in feature_type table";
	}else {
		push @lstDbTargets, $dbTarget;
	}
}

#if there are any errors in the parameters exits and show the errors
if (@lstErrors){
	foreach my $errorVal (@lstErrors){
		print $errorVal."\n";
	}
	exit;
}

#add the final bar if necessary
if ((substr $localFilePath, -1) ne '/'){
	$localFilePath.='/';
}

my @lstRegMeta;

#Create a clsRegisterMetadata object with the column headers and add it to the array
my $clsRegMetaHead = store_row('epi_accession', 'accession', 'experiment_accession', 'epigenome', 'feature_type', 'biological_replicate', 'new_biological_replicate', 'technical_replicate', 'new_technical_replicate', 'gender', 'md5_checksum', 'local_url', 'analysis', 'experimental_group', 'assay_xrefs', 'ontology_xrefs', 'xrefs', 'epigenome_description', 'control_id', 'paired', 'paired_end_tag', 'read_length', 'multiple', 'paired_with','download_url', 'info');
push @lstRegMeta, $clsRegMetaHead;


#open file
open(my $fh, '<:', $pathFile) or die "Could not open file $pathFile\n";

#read file 
while (my $fullRow = <$fh>) {
	chomp $fullRow;
  	my @row = (split "\t", $fullRow);
  	my $epiAccess = @row[0];
  	if ($skipRefEpi{$epiAccess}){
  		next;
  	}
	
	#init varaiables
	my $epi_accession= '-';
	my $accession= '-';
	my $expAccession= '-';
	my $epigenome = '-';
	my $featureType = '-';
	my $bioReplicate = '-';
	my $newBioReplicate = '-';
	my $techReplicate = '-';
	my $newTechReplicate = '-';
	my $gender = '-';
	my $md5Check = '-';
	my $localUrl = '-';
	my $analysis = '-';
	my $expGroup = '-'; 
	my $assayXrefs ='-';
	my $ontXrefs = '-';
	my $xrefs = '-';
	my $epiDesc = '-';
	my $controlId = '-';
	
	my $downUrl = '-';
	my $info = '-';
	my $fileName ='';
	
	my $epiFeature='';
	my $swBioRep = 1;
	my $swTecRep = 1;
	
	#get epigenome data
	my $rEData = get_reference_epigenome_data($epiAccess);
	if ($rEData){
		#xref
		foreach my $ref (@{$rEData->{'dbxrefs'}}){
			my @splref = split ":", $ref;
			my $dbname = @splref[0].'-';
			my $idxref;
			my $lenSplRef = scalar (@splref);
			for (my $i=1; $i < $lenSplRef; $i++) {
				if ($i > 1){
					$idxref += ':';
				}
				$idxref .= @splref[$i];
			} 
			my $tmpXref = $dbname.$idxref;
			
			if ($xrefs ne '-'){
				$xrefs .= ';'.$tmpXref;
			}else{
				$xrefs = $tmpXref;
			}	
			
		}

		$epi_accession = $rEData->{'accession'};
		
		#experimental group
		$expGroup = $rEData->{'award'}->{'project'};


	
		#Experiments
		my @experiments = @{$rEData->{'related_datasets'}};


		#add extra experiments (ctcf experiments)
		my $extraEpigenome = $rEData->{'biosample_term_name'}[0];
		my $extraSpecie = $rEData->{'organism'}[0]->{'scientific_name'};
		my $specieName;
		if (uc $extraSpecie eq 'MUS MUSCULUS'){
			$specieName = 'Mouse';
		}else{
			if (uc $extraSpecie eq 'HOMO SAPIENS'){
				$specieName = 'Human';
			}
		}

		my @extraDone;
		my @extraExperiments;

		foreach my $exp (@experiments){
			my $life_stage ='';
			my $age = '';
			if (uc $exp->{'status'} eq 'RELEASED'){
				my @replicates = @{$exp->{'replicates'}};
				my $rep = @replicates[0]; #get the information from the first replicate
				$life_stage = $rep->{'library'}->{'biosample'}->{'life_stage'};
				$age = $rep->{'library'}->{'biosample'}->{'age'};
				$age =~ s/^\s+|\s+$//g;
				my $searchParameters = $extraEpigenome.$life_stage.$age;
				if ( not grep( /^$searchParameters$/, @extraDone) ){
					
					@extraExperiments = get_extra_experiments($extraEpigenome, $extraSpecie, $life_stage, $age);
					
					push @extraDone, $searchParameters;
				}
			} 

		}

		
		
		#add extra Experiments

		if (@extraExperiments){
			push @experiments, @extraExperiments;
		}

		#only experiments with known targets or GO annotation like "DNA binding transcription factor activity"
		#also exclude experiments which target is FLAG- *or eGFP-*
		@experiments = only_known_and_filtered_tagets($adaptors, \@experiments, \%hshFeatureType, \%features_done);


		#only experiments with known targets
		#@experiments = only_known_tagets(\@experiments, \%hshFeatureType, \%featureNotFound);

		#get extra controls not contained in the reference epigenome set
		my @extraCtls = get_extra_controls(\@experiments);

		if (@extraCtls){
			push @experiments, @extraCtls;
		}


		#get control experiments info
		my $infoControl = get_control_info(\@experiments);

		my @lstControls;

		#loop experiments
		foreach my $exp (@experiments) {

			my $target;
			$info = '-';
			#my %lstMultiple;
			if ($exp->{'target'}->{'label'}){
				$target = $exp->{'target'}->{'label'};
			}

			my $swTarget =0;
			if (@lstDbTargets){
				#Target parameter has been informed
				if ($target){
					if (not( (grep( /^$target$/, @lstDbTargets)) || (index(uc($target), 'CONTROL') != -1))) {

=pod
					if ( not ((grep( /^$target$/, @lstDbTargets) ) || (uc $target eq 'CONTROL') || (uc $target eq 'RABBIT-IGG-CONTROL') || (uc $target eq '
MOUSE-IGG-CONTROL')) ) {
=cut
						#it is not a control or one of our targets
						$swTarget =1;
					} 
				}				
			}
			
			#change Control target to WCE
			if (index(uc($target), 'CONTROL') != -1) {
    			$target = 'WCE';
			}
			#if ((uc $target eq 'CONTROL') || (uc $target eq 'RABBIT-IGG-CONTROL')){
			#	$target = 'WCE';
			#}
			
			#filter experiments by assay
			my $assay = $exp->{'assay_title'};

			my $expReleased=0;
			
			#check status of the experiment
			if (uc $exp->{'status'} eq 'RELEASED'){
				$expReleased=1;
			}

			if ((grep( /^$assay$/i, @lstDbAssays)) && ($swTarget==0) && ($expReleased==1)){
				
				#experiment accession
				$expAccession = $exp->{'accession'};

				
				#normalize the assay term
				$assay = check_db_value($assay, \%hshAnalysis);

				#check and normalize target
				if ($target){
					my $DBtarget = check_db_value($target, \%hshFeatureType);
					if (!$DBtarget){
						if ($info eq '-'){
							$info = "Target: $target Not found in feature_type table";
						}else{
							$info .= "; Target: $target Not found in feature_type table";
						}
=pod
						if (!$targetsNotFound{$target}){
							$targetsNotFound{$target} = $target;
						}
=cut						
						#push @lstErrors, "experiment: ".$exp->{'accession'}."\ttarget: $target\terror: Not found in feature_type table";
					}else{
						$target = $DBtarget;
					}
				}else{
					#no target, check if assay is DNase-seq
					if (uc $assay eq 'DNASE-SEQ'){
						$target = 'DNase1';
					}else{
						#Error no target.
						push @lstErrors, "experiment: ".$exp->{'accession'}."\ttarget: Empty\terror: No target in ENCODE data";
					}
					
				}
				if ($target){
					$featureType = $target;
					$analysis = $exp->{'assay_title'};
					#check and normalize analysis
					my $DBAnalysis = check_db_value($analysis, \%hshAnalysis);
					if (!$DBAnalysis){
						if ($info eq '-'){
								$info = "Analysis: $analysis error: Not found in analysis table";
							}else{
								$info .= "; Analysis: $analysis error: Not found in analysis table";
							}
						
						#push @lstErrors, "experiment: ".$exp->{'accession'}."\tanlysis: $analysis\terror: Not found in analysis table";
					}
					$analysis = $DBAnalysis;
					#$controlId = $controlAccession;
					#$info ='-';
					if (uc $target eq 'WCE'){
						$controlId = '-';
						$info = $infoControl;
					}
				}
				
				$ontXrefs = $exp->{'biosample_term_id'};
				#replicates (Biosamples)
				
				my $epilabel;
				my @replicates = @{$exp->{'replicates'}};
				my $rep = @replicates[0]; #get the information from the first replicate
				my $termName = '-';
				if ($specieName){
					$epiDesc = $specieName;
				}
				if (uc $exp->{'biosample_term_name'} ne 'UNKNOWN'){
					$termName = $exp->{'biosample_term_name'};
					if ($epiDesc ne '-'){
						$epiDesc .= ' '.$termName;
					}else{
						$epiDesc = $termName;
					}

					$epilabel = $termName;

				}
				my $abr;
				my $lifeStage = '-';
				if (uc $rep->{'library'}->{'biosample'}->{'life_stage'} ne 'UNKNOWN'){
					$lifeStage = $rep->{'library'}->{'biosample'}->{'life_stage'};
					$epiDesc .= ' from '.$lifeStage;
					
					if (lc $lifeStage eq 'embryonic'){
						$abr = 'E';
					}else{
						if (lc $lifeStage eq 'postnatal'){
							$abr = 'P0';
						}else{
							if (lc $lifeStage eq 'adult'){
								$abr = 'adult';
							}
						}
					}

					#$epilabel .= ' '.$abr;

				}

				my $age = '-';
				if (uc $rep->{'library'}->{'biosample'}->{'age'} ne 'UNKNOWN'){
					$age = $rep->{'library'}->{'biosample'}->{'age'};
					$epiDesc .= ' '.$age;
					if ($abr ne 'P0' && $abr ne 'adult'){
						#$epilabel .= $age;
					}
					
				}
				my $units='-';
				if (uc $rep->{'library'}->{'biosample'}->{'age_units'} ne 'UNKNOWN' && $age ne '-' ){
					$units = $rep->{'library'}->{'biosample'}->{'age_units'};
					if ($age != 1){
						$units .= 's';
					}
					#$epiDesc .= ' '.$units;

				}



				#$epigenome = join (':', $termName, $lifeStage, $age, $units);
				$epigenome = $termName."_".$epi_accession;
				#$epiDesc .=';'.$epilabel;
				$epiDesc = $termName;
				$assayXrefs = $exp->{'assay_term_id'};
				
				#Normalize gender
				$gender = check_db_value($rep->{'library'}->{'biosample'}->{'sex'}, \%hshGender);
					
				#files
				my @files = @{$exp->{'files'}};
				foreach my $file (@files){
					my $paired = 'No',
					my $paired_end_tag = '-';
					my $read_length = '-';
					my $multiple = '-';
					my $paired_with = '-';
					my $fileStatus=0;
					#check file status
					if (uc $file->{'status'} eq 'RELEASED'){
						$fileStatus=1;
					}
					#get only fastq files
					if ((uc $file->{'file_format'} eq 'FASTQ') && ($fileStatus == 1)){
						
						#accession
						$accession = $file->{'accession'};
						#md5 checksum
						$md5Check = $file->{'md5sum'};
						#replicates
						$bioReplicate = $file->{'biological_replicates'}[0];
						$newBioReplicate = $bioReplicate;
						$techReplicate = (split '_', $file->{'technical_replicates'}[0])[-1];
						$newTechReplicate = $techReplicate;
						#paired and paired end tag
						if (uc $file->{'run_type'} eq 'PAIRED-ENDED' ){
							$paired ='YES';
							$paired_end_tag = $file->{'paired_end'};
							$paired_with = (split '/', $file->{'paired_with'})[-1];
						}
						
						#read length
						$read_length = $file->{'read_length'};
						
						#get controled by accessions (only in ChIP-Seq)
						$controlId='-';
						if (uc $assay eq 'CHIP-SEQ'){
							my @controlledId;
							if ($file->{'controlled_by'}){
								@controlledId = @{$file->{'controlled_by'}};
								$controlId='';
								foreach my $controledBy (@controlledId){
									if ($controlId ne ''){
										$controlId.=', ';
									}
									$controlId.= (split '/', $controledBy)[-1];
								}
							}else{
								#no controled_by field
								#check if there is a control at experiment level


								$controlId = get_controls_from_experiment($exp->{'possible_controls'});

								#push @lstErrors, "file: ".$accession."\tControlled_by: Missing\terror: Controlled_by field missing";
							}
						}
					


						#dowload url
						if ($file->{'href'}){
							$downUrl = 'https://www.encodeproject.org'.$file->{'href'};
							$localUrl = $localFilePath.(split '/', $downUrl)[-1];
						}
						
						if (uc $target eq 'WCE'){
							#assign values to object Register Metadata
							$controlId = '-';
							my $clsCtlRegMeta = store_row($epi_accession, $accession, $expAccession, $epigenome, $featureType, $bioReplicate, $newBioReplicate, $techReplicate, $newTechReplicate, $gender, $md5Check, $localUrl, $analysis, $expGroup, $assayXrefs, $ontXrefs, $xrefs, $epiDesc, $controlId, $paired, $paired_end_tag, $read_length, $multiple, $paired_with, $downUrl, $info);
							#Add the object to the list of control rows
							push @lstControls, $clsCtlRegMeta;
							
						}else{
							#assign values to object Register Metadata
							my $clsRegMeta = store_row($epi_accession, $accession, $expAccession, $epigenome, $featureType, $bioReplicate, $newBioReplicate, $techReplicate, $newTechReplicate, $gender, $md5Check, $localUrl, $analysis, $expGroup, $assayXrefs, $ontXrefs, $xrefs, $epiDesc, $controlId, $paired, $paired_end_tag, $read_length, $multiple, $paired_with, $downUrl, $info);
							#Add the object to the list of rows
							push @lstRegMeta, $clsRegMeta;
							
						}
						
					}
				}
					
			}

		}
		#add control rows
		push (@lstRegMeta, @lstControls);
	}
}

#close file
close($fh);



#print data or Errors
if (@lstErrors){
	foreach my $errorVal (@lstErrors){
		print $errorVal."\n";
	}
}else{
	#final check and updates


	@lstRegMeta = check_duplicated_entries(\@lstRegMeta);

	@lstRegMeta = join_controls_controlling_same_experiment(\@lstRegMeta);

	@lstRegMeta = update_paired_end_controlled_by_column(\@lstRegMeta);

	@lstRegMeta = remove_paired_signal_and_controls_with_missing_parts(\@lstRegMeta);

	@lstRegMeta = remove_unused_controls(\@lstRegMeta);

	



	###@lstRegMeta = check_replicate_gaps(\@lstRegMeta);

	###@lstRegMeta = check_duplicated_replicates(\@lstRegMeta);

	###@lstRegMeta = update_column_multiple(\@lstRegMeta);




	@lstRegMeta = warns_paired_and_single_end_mixed(\@lstRegMeta);

	@lstRegMeta = check_control_analysis($adaptors, \@lstRegMeta);




	###@lstRegMeta = check_paired_end_replicates($adaptors, \@lstRegMeta);
	###@lstRegMeta = change_multiple_and_tech_replicat_for_DNSA1($adaptors, \@lstRegMeta);



	

	@lstRegMeta = update_analysis_for_paired_controls(\@lstRegMeta);

	

	@lstRegMeta = remove_signal_files_without_control(\@lstRegMeta);

	@lstRegMeta = check_gender(\@lstRegMeta);


	@lstRegMeta = check_replicate_numbers_and_multiple_column(\@lstRegMeta);



	foreach my $reg (@lstRegMeta) {
		print $reg->csv_row."\n";
	}

}

##################### END MAIN FUNCTION ##############################

sub check_gender {
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my $current_epigenome;
	my $current_gender;
	my %hshGender;
	foreach my $regMeta (@lstRegMeta){
		if ($current_epigenome ne $regMeta->get('epi_accession')){
			$current_epigenome = $regMeta->get('epi_accession');
			$current_gender = $regMeta->get('gender');
		}else{
			if ($current_gender ne $regMeta->get('gender')){
				if (!exists($hshGender{$current_epigenome})){
					$hshGender{$current_epigenome} = 'mixed';
				}
			}
		}
	}

	foreach my $regMeta (@lstRegMeta){
		if (exists $hshGender{$regMeta->get('epi_accession')}){
			$regMeta->set_gender($hshGender{$regMeta->get('epi_accession')});
		}

	}

	return @lstRegMeta;
}

sub remove_paired_signal_and_controls_with_missing_parts {
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my $index = 0;
	my @to_remove;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('paired') eq 'YES' && $index > 0){
			my $paired_found = get_row_by_accession($regMeta->get('paired_with'), \@lstRegMeta);
			if (!$paired_found){
				push @to_remove, $index;
			}
		}
		$index +=1;

	}

	my @sorted_to_remove = sort{$b <=> $a} @to_remove;
	foreach my $del (@sorted_to_remove) {
		splice @lstRegMeta, $del, 1;
	}

	return @lstRegMeta;
}

sub remove_signal_files_without_control {
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my $index = 0;
	my @to_remove;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('feature_type') ne 'WCE' && $index > 0 && uc $regMeta->get('feature_type') ne 'DNASE1'){
			my $fullCtl = $regMeta->get('control_id');
			if ($fullCtl && $fullCtl ne '-'){
				my @lstCtl = (split ",", $fullCtl);
				foreach my $ctl (@lstCtl){
					$ctl =~ s/^\s+|\s+$//g;
					my $found = search_control($ctl, \@lstRegMeta);
					if ($found == 0){
						push @to_remove, $index;
						#splice @lstRegMeta, $index, 1;
					}
				}				
			}else{
				push @to_remove, $index;
				#splice @lstRegMeta, $index, 1;
			}

		}
		$index += 1;
	}

	my @sorted_to_remove = sort{$b <=> $a} @to_remove;
	foreach my $del (@sorted_to_remove) {
		splice @lstRegMeta, $del, 1;
	}

	return @lstRegMeta;
}

sub search_control{
	my $ctlAccession = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my $found = 0;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('feature_type') eq 'WCE'){
			if (uc $regMeta->get('accession') eq uc $ctlAccession){
				$found = 1;
				last;
			}
		}
	}
	return $found;
}

sub get_controls_from_experiment {
	my $possible_ctls = shift;
	my $controlId='-';
	foreach my $ctl (@{$possible_ctls}) {
		$controlId='';
		foreach my $file (@{$ctl->{'files'}}) {
			my $ctl_file = (split '/', $file)[-1];
			my $file_format = get_file_format($ctl_file);
			if (uc $file_format ne 'FASTQ'){
				next;
			}
			if ($controlId ne ''){
				$controlId.=', ';
			}
			$controlId.= $ctl_file;
		}
	}

	return $controlId;
}

sub get_file_format {
	my $file_accession = shift;
	my $f_format = '-';
	my $url = 'https://www.encodeproject.org/files/'.$file_accession.'/?frame=embedded&format=json';
	my $client = REST::Client->new();
	my $headers = {Accept => 'application/json'};
	$client->GET($url, $headers);
	my $response = decode_json($client->responseContent());
	if ($response){
		if (uc $response->{'status'} eq 'RELEASED'){
			$f_format = $response->{'file_format'};
		}
	}

	return $f_format;
}

sub only_known_and_filtered_tagets{

	my $adaptors = shift;
	my $lstExp=shift;
	my $features = shift;
	my $done_features = shift;
	
	my %db_features = %{$features};
	my @experiments = @{$lstExp};
	my @filetred_exps;
	
	open(my $f_not_imported, '>>', 'Features_not_imported.txt') or die "Could not open file 'Features_not_imported.txt' $!";
	open(my $f_to_be_imported, '>>', 'new_Features_to_be_imported.txt') or die "Could not open file 'new_Features_to_be_imported.txt' $!";

	foreach my $exp (@experiments){

		my $target = $exp->{'target'}->{'label'};

		if (uc $exp->{'assay_title'} eq 'DNASE-SEQ'){
			$target = 'DNase1';
		}

		if (!$target){
			next;
		}

		if (index(uc($target), 'CONTROL') != -1) {
			push @filetred_exps, $exp;
			next;
		}

		if (exists ($done_features->{uc $target})){
			if ($done_features->{uc $target} == 1){
				push @filetred_exps, $exp;
			}
			next;
		}

		if (exists ($db_features{uc $target})){
			$done_features->{uc $target} = 1;
			push @filetred_exps, $exp;
			next;
		}

		$done_features->{uc $target} = must_be_imported($target, $adaptors);

		if ($done_features->{uc $target} == 1){
			push @filetred_exps, $exp;
			say $f_to_be_imported $target;
		}else{
			say $f_not_imported $target;
		}

	}

	close $f_to_be_imported;
	close $f_not_imported;
	return @filetred_exps;

}

sub must_be_imported {
	my $target = shift;
	my $adaptors = shift;
	my $add = 0;
	
	my $lstgene = $adaptors->{Gene}->fetch_all_by_external_name($target);

	LBL_GENE: {
		foreach my $gene (@{$lstgene}){
			my $lstDbEntries = $gene->get_all_DBLinks('GO');
			foreach my $dbentry (@{$lstDbEntries}){
				if (index(uc $dbentry->description(), 'DNA BINDING TRANSCRIPTION FACTOR ACTIVITY') != -1) {
	    			$add = 1;
	    			last LBL_GENE;
				} 
			}
	
		}
	}
	
	return $add;

}

sub check_replicate_numbers_and_multiple_column{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @epigenomeDone;
	my $firstRow=1;

	foreach my $regMeta (@lstRegMeta){
		if ($firstRow == 1){
			$firstRow=0;
			next;
		}
		my $epigenome = $regMeta->get('epigenome');
		if (not grep( /^$epigenome$/, @epigenomeDone)){
			push  (@epigenomeDone, $epigenome);
			my @epiObjects = get_filtered_objects('epigenome', $epigenome, \@lstRegMeta);
			my @featureDone;

			foreach my $feature (@epiObjects){
				my $feature_type = $feature->get('feature_type');
				if (not grep( /^$feature_type$/, @featureDone)){
					push (@featureDone, $feature_type);
					my @featureEpigenome = get_filtered_objects('feature_type', $feature_type, \@epiObjects);
					
					#create the tree
					my %replicate_tree = create_replicate_tree(\@featureEpigenome);

					#print Dumper %replicate_tree;
					#exit;
					
					#numerate column "multiple"
					update_column_multiple_in_tree(\@featureEpigenome, \%replicate_tree);

					#numerate replicates
					numerate_replicates(\%replicate_tree);

				}
			}
		}

	}
	return @lstRegMeta;
}

sub numerate_replicates{
	my $experiment = shift;
	my $bio_count =1;
	foreach my $exp (keys %{$experiment}){
		foreach my $bio (keys %{$experiment->{$exp}}){
			my $tec_count = 1;
			foreach my $tec (keys %{$experiment->{$exp}{$bio}}){
				my @objs = @{$experiment->{$exp}{$bio}{$tec}};
				foreach my $obj (@objs){
					$obj->set_new_bio_replicate($bio_count);
					$obj->set_new_tech_replicate($tec_count);
				}
				$tec_count +=1;
			}
			$bio_count +=1;
		}
	}

	return;
}

sub update_column_multiple_in_tree{
	my $objs = shift;
	my $experiment = shift;
	my @featureEpigenome = @{$objs};
	my %obj_done;
	
	
	foreach my $fe (@featureEpigenome){
		if (lc $fe->get('acession') eq 'accession'){
			next;
		}
		if (not(exists $obj_done{$fe->get('experiment_accession').'_'.$fe->get('biological_replicate').'_'.$fe->get('technical_replicate')})){
			$obj_done{$fe->get('experiment_accession').'_'.$fe->get('biological_replicate').'_'.$fe->get('technical_replicate')}=1;

			my @objs = @{$experiment->{$fe->get('experiment_accession')}{$fe->get('biological_replicate')}{$fe->get('technical_replicate')}};
			my $counter=1;
			my %paired_num;
			my $paired_count = 1;
			foreach my $obj (@objs){
				if (uc $obj->get('paired') eq 'NO'){
					$obj->set_multiple($counter);
					$counter +=1;
				}else{
					if (uc $obj->get('paired') eq 'YES'){
						if (not(exists $paired_num{$obj->get('accession')})){
							$paired_num{$obj->get('paired_with')} = $paired_count;
							$obj->set_multiple ($paired_count);
							$paired_count += 1;
						}else{
							$obj->set_multiple ($paired_num{$obj->get('accession')});
						}
						
					}
				}
			}
		}
	}
	return;
}


sub create_replicate_tree{
	my $objs = shift;
	my @featureEpigenome = @{$objs};
	my %experiment;
	foreach my $fe (@featureEpigenome){
		#if (not (exists $experiment{$fe->get('experiment_accession')}{$fe->get('biological_replicate')}{$fe->get('technical_replicate')})){
		#	$experiment{$fe->get('experiment_accession')}{$fe->get('biological_replicate')}{$fe->get('technical_replicate')}=[];
		
		#}
		push @{$experiment{$fe->get('experiment_accession')}{$fe->get('biological_replicate')}{$fe->get('technical_replicate')}}, $fe;

	} 
	return %experiment;
}

sub get_filtered_objects{
    my $filterName = shift;
    my $filterValue = shift;
    my $lstObjs = shift;
    my @lstRegMeta = @{$lstObjs};
    my @filterObjs;
    foreach my $obj (@lstRegMeta){
        if ($filterValue eq $obj->get($filterName)) {
            push @filterObjs, $obj;
        }
    }
    return @filterObjs;
}



sub update_analysis_for_paired_controls{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('paired') eq 'YES' && $regMeta->get('feature_type') eq 'WCE' && $regMeta->get('paired_end_tag') eq '2'){
			my $paired1ctl = get_row_by_accession($regMeta->get('paired_with'), \@lstRegMeta);
			$regMeta->set_analysis($paired1ctl->get('analysis'));
		}
	}
	return @lstRegMeta;
}


sub update_paired_end_controlled_by_column{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('paired') eq 'YES' && $regMeta->get('control_id') eq '-'){
			my $fileFound = get_row_by_accession($regMeta->get('paired_with'), \@lstRegMeta);
			if ($fileFound){
				$regMeta->set_control_id($fileFound->{'control_id'});
			}
		}

	}
	return @lstRegMeta;
}

sub get_row_by_accession{
	my $accession = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my $fileFound;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('accession') eq uc $accession){
			$fileFound = $regMeta;
			last;
		}
	}

	return $fileFound;
}

sub get_extra_controls{
	my $exps = shift;
	my @experiments = @{$exps};
	my %expDone;
	my @ctlsDone;
	my @extraCtls;

	foreach my $exp (@experiments){
		if (not (exists $expDone{$exp->{'accession'}})){
			$expDone{$exp->{'accession'}}=$exp->{'accession'};
			my $target = $exp->{'target'}->{'label'};
			if (not (index(uc($target), 'CONTROL') != -1)) {

				foreach my $ctl (@{$exp->{'possible_controls'}}){
					my $ctlAccession = $ctl->{'accession'};
					if ( ($ctlAccession) && (not grep( /^$ctlAccession$/, @ctlsDone)) ){
						push @ctlsDone, $ctl->{'accession'};
						if (not (ctl_alredy_included($ctlAccession, \@experiments))){

							my $ctlExp = get_control($ctlAccession);
							
							push @extraCtls, $ctlExp;
						}

					}
				}
			}
		}
		
	}

	return @extraCtls;

}

sub ctl_alredy_included{
	my $ctl = shift;
	my $exps = shift;
	my $found;
	foreach my $exp (@{$exps}){
		my $target = $exp->{'target'}->{'label'};
		if (index(uc($target), 'CONTROL') != -1) {
			if (uc ($ctl) eq uc ($exp->{'accession'}) ){
				$found=1;
				last;
			}
		}
	}
	return $found;
}

sub only_known_tagets{

	my $lstExp=shift;
	my $features = shift;
	my $notFound = shift;
	my @expsKnownTarget = @{$lstExp};
	my %hshCompare = %{$features};
	my @onlyFoundTargets;
	
	open(my $f_not_found, '>>', 'Features_not_found.txt') or die "Could not open file 'Features_not_found.txt' $!";


	foreach my $extraExp (@expsKnownTarget){
		my $extraTarget;
		if ($extraExp->{'target'}->{'label'}){
			$extraTarget = $extraExp->{'target'}->{'label'};
		}
		if (not (index(uc($extraTarget), 'CONTROL') != -1)) {
=pod
		if ( not (uc $extraTarget eq 'CONTROL' || uc $extraTarget eq 'RABBIT-IGG-CONTROL') || (uc $extraTarget eq '
MOUSE-IGG-CONTROL')) {
=cut
			#it is not a control

			if (uc $extraExp->{'assay_title'} eq 'DNASE-SEQ'){
				$extraTarget = 'DNase1';
			}
			if ($hshCompare{uc $extraTarget}){
				push @onlyFoundTargets, $extraExp;
			}else{
				#target not found in the DB
				if ($extraTarget){
					if (not exists($notFound->{$extraTarget})){
						say $f_not_found $extraTarget;
						$notFound->{$extraTarget}=$extraTarget;
					}

				}
				

			}
		}else{
			push @onlyFoundTargets, $extraExp;
		}

	}

	close $f_not_found;
	return @onlyFoundTargets;

}

sub change_multiple_and_tech_replicat_for_DNSA1{
	my %adaptors = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstExp;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('paired') eq 'NO' && uc $regMeta->get('feature_type') eq 'DNASE1'){
			my $ext_exp = $regMeta->get('experiment_accession');
			if (not grep( /^$ext_exp$/, @lstExp)){
				push @lstExp, $ext_exp;
				my $exps = get_experiment_objects($ext_exp, \@lstRegMeta);
				my %reps;
				foreach my $exp (@{$exps}){
				
					if (!$reps{$exp->get('new_bio_replicate')}){
						$reps{$exp->get('new_bio_replicate')}=1;
					}
					if(uc $exp->get('paired') eq 'NO'){
						$exp->set_new_tech_replicate($reps{$exp->get('new_bio_replicate')});
						$exp->set_multiple(1);
					}


					$reps{$exp->get('new_bio_replicate')} = $reps{$exp->get('new_bio_replicate')} + 1;


				}

			}
		}
	}
	return @lstRegMeta;
}




sub check_paired_end_replicates{
	my %adaptors = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstPaired;
	foreach my $regMeta (@lstRegMeta){
		#get only paired end files
		if (uc $regMeta->get('paired') eq 'YES'){
			push @lstPaired, $regMeta;
		}

	}
	my @lstExps;
	foreach my $regMeta (@lstPaired){
		my $ext_exp = $regMeta->get('experiment_accession');
		if (not grep( /^$ext_exp$/, @lstExps)){
			push @lstExps, $ext_exp;
			#get the ones with the same encode experiment
			my $exps = get_experiment_objects($ext_exp, \@lstPaired);
			my @paired;
			my @rep_num;
			my %tech_rep;
			foreach my $exp (@{$exps}){
				if(uc $exp->get('paired') eq 'YES'){
					my $accession = $exp->get('accession');
					my $acc_paired = $exp->get('paired_with');
					if (not grep( /^$accession$/, @paired)){
						my $bio_rep = $exp->get('new_bio_replicate');
						if (!$tech_rep{$bio_rep}){
							$tech_rep{$bio_rep}=1;
						}

						$exp->set_new_tech_replicate($tech_rep{$bio_rep});
						$exp->set_multiple(1);
						foreach my $paired_row (@{$exps}){
							if ($paired_row->get('accession') eq $acc_paired){
								$paired_row->set_new_tech_replicate($tech_rep{$bio_rep});
								$paired_row->set_multiple(1);
							}
						}
				


						#add them to done list
						push @paired, $accession;
						push @paired, $acc_paired;

						$tech_rep{$bio_rep} = $tech_rep{$bio_rep} + 1;

					}

				}



			}

		}
		
	}

	return @lstRegMeta;


}


sub check_control_analysis{
	my %adaptors = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};

	foreach my $regMeta (@lstRegMeta){
		if ($regMeta->get('feature_type') eq 'WCE'){
			$regMeta->set_analysis('');
			my @lstFt = get_histone_TF($regMeta->get('accession'), \%adaptors, \@lstRegMeta);
			my $first =0;
			my $feature='-';
			foreach my $hTF (@lstFt){
				if ($first != 0){
					$feature .= ', ';
				}else{
					$feature ='';
				}

				$feature .= $hTF;
				$first = 1;
			}
			$regMeta->set_analysis($feature);
		}

	}
	return @lstRegMeta;

}

sub get_histone_TF{
	my $ctlAccession = shift;
	my $adap = shift;
	my $lstObjs = shift;
	my %adaptors = %{$adap};
	my @lstRegMeta = @{$lstObjs};
	my $lstFeatureTypes = $adaptors->{FeatureType}->fetch_all();
	my %hshFeatureType;
	my @hTF;
	foreach my $objFeatureType (@{$lstFeatureTypes}){
		$hshFeatureType{uc $objFeatureType->name}=$objFeatureType->class;
	} 


	foreach my $regMeta (@lstRegMeta){
		if ($regMeta->get('feature_type') ne 'WCE'){

			my @ctl_by = split (',', $regMeta->get('control_id'));

			foreach my $controlled_by (@ctl_by){
				$controlled_by =~ s/\s//g; #remove blank spaces
				if ($controlled_by eq $ctlAccession){

					my $ftClass = $hshFeatureType{uc $regMeta->get('feature_type')};
					if (uc $ftClass eq 'HISTONE'){
							my $ctlHist = 'ChIP-Seq_control_hist';
							if (not grep( /^$ctlHist$/, @hTF)){
								push @hTF, $ctlHist;
							}
					}else{
						if (uc $ftClass eq 'TRANSCRIPTION FACTOR'){
							my $ctlTF = 'ChIP-Seq_control_TF';
							if (not grep( /^$ctlTF$/, @hTF)){
								push @hTF, $ctlTF;
							}
							
						}else{
							#************ warning TF harcoded *************
							#In this case, All new Feature types are TFs 
							my $ctlTF = 'ChIP-Seq_control_TF';
							if (not grep( /^$ctlTF$/, @hTF)){
								push @hTF, $ctlTF;
							}
							#**********************************************

							my $info = $regMeta->get('info');
							if ($info ne '-'){
								$info .= '; Feature Class not Histone or TF';
							}else{
								$info = 'Feature Class not Histone or TF';
							}
							$regMeta->set_info($info);
						}
					}
				}

			}

		}
	}

	return @hTF;



}

sub get_extra_experiments {
	my $epigenome = shift;
	my $specie = shift;
	my $life_stage = shift;
	my $age = shift;
	

	#replace blank spaces
	my $find = " ";
	my $replace = "+";
	$find = quotemeta $find; # escape regex metachars if present
	$epigenome =~ s/$find/$replace/g;
	$specie =~ s/$find/$replace/g;

	my $TFurl = 'https://www.encodeproject.org/search/?searchTerm='.$epigenome;
	$TFurl .= '&type=Experiment&assay_slims=DNA+binding';
	$TFurl .= '&replicates.library.biosample.donor.organism.scientific_name='.$specie;
	$TFurl .= '&limit=all';
	$TFurl .= '&target.investigated_as=transcription+factor';
	$TFurl .= '&assay_title=ChIP-seq';
	$TFurl .= '&replicates.library.biosample.life_stage='.$life_stage;
	$TFurl .= '&replicates.library.biosample.age='.$age;
	$TFurl .= '&status=released';
	$TFurl .= '&frame=embedded&format=json';

	my $clientAwd = REST::Client->new();
	my $headersAwd = {Accept => 'application/json'};
	$clientAwd->GET($TFurl, $headersAwd);
	my $responseAwd = decode_json($clientAwd->responseContent());

	my @extraExperiments;
	my @ctlsDone;
	foreach my $results (@{$responseAwd->{'@graph'}}){
		my $treatment = $results->{'replicates'}[0]->{'library'}->{'biosample'}->{'treatments'}[0]->{'treatment_term_name'};
		
		if (!$treatment){ #not treated
			push @extraExperiments, $results;

			#get controls
			foreach my $ctl (@{$results->{'possible_controls'}}){
				if ( ($ctl->{'accession'}) && (not grep( /^$ctl$/, @ctlsDone)) ){
					my $ctlExp = get_control($ctl->{'accession'});
					push @ctlsDone, $ctl->{'accession'};
					push @extraExperiments, $ctlExp;
				}

			}

		}
	}

	return @extraExperiments;
}


sub join_controls_controlling_same_experiment {
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @signalDone;
	my $count=0;
	my $count1=0;
	foreach my $regMeta (@lstRegMeta){
		
		if ($regMeta->get('feature_type') ne 'WCE'){
			
			my $expsig = $regMeta->get('experiment_accession');
			if (not grep( /^$expsig$/, @signalDone)){
				push @signalDone, $expsig;
				
				my @expSignal = get_objects_by_experiments($expsig, \@lstRegMeta); #get signal rows controlled by the same experiment
				
				my @lstctlExp;
				foreach my $exp (@expSignal){
					my @multi_ctl = split (',', $exp->get('control_id'));
					foreach my $controlled_by (@multi_ctl){
						$controlled_by =~ s/\s//g; #remove blank spaces
						if ($controlled_by ne '-'){
							my $ctlExp = get_ctl_experiment($controlled_by, \@lstRegMeta);
			                if (not grep( /^$ctlExp$/, @lstctlExp)){
			                    push @lstctlExp, $ctlExp;
			                }
						}

					}
				}
				my $numExp = @lstctlExp;
				
				if ($numExp > 1){
					
					my $newExpAccession;
		            foreach my $expAcc (@lstctlExp){
		                if ($newExpAccession){
							$newExpAccession .= '_'.$expAcc;
	                	}else{
	                		$newExpAccession = $expAcc;
	                	}
	 	            }
	 	            foreach my $expAcc (@lstctlExp){
						my @lstExp = get_objects_by_experiments($expAcc, \@lstRegMeta);
						foreach my $exp (@lstExp){
							$exp->set_experiment_accession($newExpAccession);
						}

	 	            }
				}

			}
		}


	}

	
	return @lstRegMeta;
}


sub get_objects_by_experiments{
    my $expAcc = shift;
    my $lstObjs = shift;
    my @lstRegMeta = @{$lstObjs};
    my @lstSig;
    foreach my $sig (@lstRegMeta){
        if ($expAcc eq $sig->get('experiment_accession')) {
            push @lstSig, $sig;
        }
    }
    return @lstSig;
}

sub get_ctl_experiment{
    my $ctlExp = shift;
    my $lstObjs = shift;
    my @lstRegMeta = @{$lstObjs};
    my $ctlExp2;
    foreach my $ctl (@lstRegMeta){
		if ($ctl->get('feature_type') eq 'WCE'){
	        if ($ctl->{accession} eq $ctlExp){
	            $ctlExp2= $ctl->get('experiment_accession');
	        }
		}


    }
    return $ctlExp2;
}


sub warns_paired_and_single_end_mixed{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @expDone;
	foreach my $regMeta (@lstRegMeta){
		my $exp = $regMeta->{'experiment_accession'};
		if (not grep( /^$exp$/, @expDone)){
			push @expDone, $exp;
			my $lstExperiment = get_experiment_objects($exp, \@lstRegMeta);
			my $firstExp=-1;
			my $paired;
			my $multiple = 0;
			foreach my $objSelected (@{$lstExperiment}){
				if ($firstExp == -1){
					$paired = $objSelected->{'paired'};
					$firstExp = 0;
				}else{
					if (uc $paired ne uc $objSelected->{'paired'} && $multiple == 0){
						$multiple = -1;
					}
				}
			}
			if ($multiple == -1){
				foreach my $objSelected (@{$lstExperiment}){
					my $info = $objSelected->{'info'};
					$info .= 'WARNING: Multiple paired';
					$objSelected->set_info($info);
				}
			} 

		}
	}
	return @lstRegMeta;
}

sub get_control{
	my $ctlExpAccession = shift;
	my $ctlUrl = 'https://www.encodeproject.org/experiments/'.$ctlExpAccession.'/?frame=embedded&format=json';
	my $client = REST::Client->new();
	my $headers = {Accept => 'application/json'};
	$client->GET($ctlUrl, $headers);
	my $response = eval{decode_json($client->responseContent());};

	unless($response) {
		print $ctlUrl."\n";
    	print $@."\n";
    	exit;
	}
	return $response;
}

sub check_duplicated_entries{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @storedAccessions;
	
	#load controls in a hash
	my @to_remove;
	my $ind_obj = 0;
	foreach my $regMeta (@lstRegMeta){
		my $accession = $regMeta->get('accession');
		if (not grep( /^$accession$/, @storedAccessions)){
			push @storedAccessions, $accession;
		}else{
			push @to_remove, $ind_obj;
		}
		$ind_obj += 1;
	}

	my @sorted_to_remove = sort{$b <=> $a} @to_remove;
	foreach my $del (@sorted_to_remove) {
		splice @lstRegMeta, $del, 1;
	}


	return @lstRegMeta;
}

sub remove_unused_controls{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @ctls;
	#get list of used controls
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('feature_type') ne 'WCE'){
			my $fullCtl = $regMeta->get('control_id');
			my @lstCtl = (split ",", $fullCtl);
			foreach my $ctl (@lstCtl){
				$ctl =~ s/^\s+|\s+$//g;
				if (not (grep(/^$ctl$/, @ctls))){
					push @ctls, $ctl;
				}
			}
		}
	}

	#add to the ctls list the second control file of paired controls
	foreach my $regMeta (@lstRegMeta){
		my $accession = $regMeta->get('accession');
		if (grep(/^$accession$/, @ctls)){
			my $pairedCtl = $regMeta->get('paired_with');
			if ($pairedCtl ne '-'){
				if (not (grep(/^$pairedCtl$/, @ctls))){
					push @ctls, $pairedCtl;
				}
			}
		}
	}

	my @to_remove;
	my $ind_count=0;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('feature_type') eq 'WCE'){
			my $accession = $regMeta->get('accession');
			if (not (grep(/^$accession$/, @ctls))){
				#the control is not used
				push @to_remove, $ind_count;
				#splice @lstRegMeta, $index, 1;
			}
		}
		$ind_count += 1;
	}

	my @sorted_to_remove = sort{$b <=> $a} @to_remove;
	foreach my $del (@sorted_to_remove) {
		splice @lstRegMeta, $del, 1;
	}


	return @lstRegMeta;
}


sub check_duplicated_replicates{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	
	my @lstDone;
	my @lstRepeated;

	foreach my $regMeta (@lstRegMeta){
		if ($regMeta->get('accession') ne 'accession'){ #scape title row
			if ($regMeta->get('feature_type') ne 'WCE'){ #exclude controls
				my $epiFeature = $regMeta->get('epigenome').$regMeta->get('feature_type'); #create a key combining epigenome and feature type
				if ( not (grep(/^$epiFeature$/, @lstRepeated))){ 
					#combination epigenome/feature still not checked
					push @lstRepeated, $epiFeature;
					my @lstExp = get_experiments_epigenome_feature($regMeta->get('epigenome'), $regMeta->get('feature_type'), \@lstRegMeta);
					my $numberExperiments = @lstExp; #count experiments for the same epigenome and feature type
					if ($numberExperiments > 1){ #if there is more than one experiment, renumerate replicates if necessary
						renumerate_replicates(\@lstExp, \@lstRegMeta);
					}

					#check experiments controlled for more than one control experiment
					my $ExpCtl =  get_control_experiments_for_experiment(\@lstExp, \@lstRegMeta);
					my %hshExpCtl = %{$ExpCtl};
					foreach my $exp (@lstExp){
						my @lstCtls = $hshExpCtl{$exp};
						my $numberExpCtls = @lstCtls;
						if ($numberExpCtls > 1){
							renumerate_replicates(\@lstCtls, \@lstRegMeta);
						}
					}
				}
			}
		}
	}

	return @lstRegMeta;
}

sub renumerate_replicates{
	my $exps = shift;
	my $regs = shift; 
	my @lstExp = @{$exps};
	my @lstRegMeta = @{$regs};
	my %renumerated;
	my $lstObjExperiments = get_multi_experiment_objects(\@lstExp, \@lstRegMeta);
	foreach my $objSelected (@{$lstObjExperiments}){
		my $combinatedKey = $objSelected->get('epigenome').$objSelected->get('feature_type').$objSelected->get('new_bio_replicate').$objSelected->get('new_tech_replicate');
		if ($renumerated{$combinatedKey}){
			#already exists
			my $exit=0;
			my $index=$objSelected->get('new_bio_replicate');
			my $new_combinatedKey;

			while ($exit==0){
				$index ++;
				$new_combinatedKey = $objSelected->get('epigenome').$objSelected->get('feature_type').$index.$objSelected->get('new_tech_replicate');
				if (!$renumerated{$new_combinatedKey}){
					$renumerated{$new_combinatedKey}=$index;
					$exit=1;
				}
			}
			$objSelected->set_new_bio_replicate($renumerated{$new_combinatedKey});
			#print $objSelected->get('accession').'_'.$objSelected->get('experiment_accession').'_'.$objSelected->get('epigenome').'_'.$objSelected->get('feature_type').'_new_'.$objSelected->get('new_bio_replicate')."\n";
		}else{
			#does not exists
			$renumerated{$combinatedKey}=$objSelected->get('new_bio_replicate');
			#print $objSelected->get('accession').'_'.$objSelected->get('experiment_accession').'_'.$objSelected->get('epigenome').'_'.$objSelected->get('feature_type').'_old_'.$objSelected->get('new_bio_replicate')."\n";
		}

	}

	return;
}

sub check_replicate_gaps{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstDone;
	foreach my $regMeta (@lstRegMeta){
		my $epigenome = $regMeta->get('epigenome');
		my $feature = $regMeta->get('feature_type');
		
		my $epi_feature = $epigenome.$feature;

		if ( not (grep(/^$epi_feature$/, @lstDone)) && ($regMeta->get('accession') ne 'accession')) {
			#not checked yet		
			#Add the experiment to the cheked list
			push @lstDone, $epi_feature;
			my @lstExpEigenome = get_experiments_epigenome_feature($epigenome, $feature, \@lstRegMeta);
			foreach my $expEpigenomeFeature (@lstExpEigenome){
				my $lstExperiment = get_experiment_objects($expEpigenomeFeature, \@lstRegMeta);
				my @bioRep;
				my @techRep;
				my $i=0;
				#load replicate numbers
				foreach my $objSelected (@{$lstExperiment}){
					$bioRep[$i][0] = $objSelected->get('biological_replicate');
					$bioRep[$i][1] = $objSelected;
					
					$techRep[$i][0] = $objSelected->get('technical_replicate');
					$techRep[$i][1] = $objSelected;
					$i++;
				}
				#sort replicates
				my @sortedBio = sort { $a->[0] <=> $b->[0] } @bioRep;
				my @sortedTech = sort { $a->[0] <=> $b->[0] } @techRep;
				#find gaps in Bio_Replicates and Fix them
				my $numRep=0;
				my $oldVal;
				for my $rowBioIndex ( 0..$#sortedBio ){
					my $val = $sortedBio[$rowBioIndex][0];
					if ($oldVal != $val){
						$oldVal = $val;
						$numRep ++;
					}
					if ($val > $numRep){
						#wrong enumeration, need to be fixed
						my $regObj = $sortedBio[$rowBioIndex][1];
						$regObj->set_new_bio_replicate($numRep);
					}				
					
				}

				#find gaps in Tech_Replicates and Fix them
				$numRep=0;
				$oldVal=-1;
				for my $rowTechIndex ( 0..$#sortedTech ){
					my $val = $sortedTech[$rowTechIndex][0];
					if ($oldVal != $val){
						$oldVal = $val;
						$numRep ++;
					}
					if ($val > $numRep){
						#wrong enumeration, need to be fixed
						my $regObj = $sortedTech[$rowTechIndex][1];
						$regObj->set_new_tech_replicate($numRep);
					}
				}
			}
		
		}
	} 
	
	return @lstRegMeta;
}
=pod
sub update_column_multiple{
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstDone;
	foreach my $regMeta (@lstRegMeta){
		my $epigenome = $regMeta->get('epigenome');
		my $feature = $regMeta->get('feature_type');
		
		my $epi_feature = $epigenome.$feature;
		if ( not (grep(/^$epi_feature$/, @lstDone)) && ($regMeta->get('accession') ne 'accession')) {
			push @lstDone, $epi_feature;
			my @lstExpEigenome = get_experiments_epigenome_feature($epigenome, $feature, \@lstRegMeta);
			foreach my $expEpigenomeFeature (@lstExpEigenome){
				
				my $lstExperiment = get_experiment_objects($expEpigenomeFeature, \@lstRegMeta);

				#update multiple column based on the new replicate enumeration columns
				my %hshMultiple;
				my $valMultiple = 0;
				foreach my $objSelected (@{$lstExperiment}){
					my $nBioRep = $objSelected->get('new_bio_replicate');
					my $nTechRep = $objSelected->get('new_tech_replicate');
					my $nPairedEnd = $objSelected->get('paired_end_tag');
					my $myKey = $nBioRep.'_'.$nTechRep.'_'.$nPairedEnd;
					
					if (exists $hshMultiple{$myKey}) {
						$valMultiple = $hshMultiple{$myKey};
						$valMultiple ++;
						$hshMultiple{$myKey} = $valMultiple;
					}else{
						$hshMultiple{$myKey} = 1;
						$valMultiple = 1;
					}
					$objSelected->set_multiple($valMultiple);
				}
			}

		}
	}
	return @lstRegMeta;
}
=cut
sub get_experiments_epigenome_feature{
	my $epigenome = shift;
	my $feature = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstExperiments;
	my @lstDone;
	foreach my $regMeta (@lstRegMeta){
		if ($regMeta->get('epigenome') eq $epigenome && $regMeta->get('feature_type') eq $feature && $regMeta->get('accession') ne 'accession') {
			my $exp = $regMeta->get('experiment_accession');
			if ( not (grep(/^$exp$/, @lstDone))){
				push @lstDone, $exp;
				push @lstExperiments, $exp;
			}
		}
	}
	return @lstExperiments;
}

sub get_control_experiments_for_experiment{
	my $exps = shift;
	my $regs = shift;
	my @lstExp = @{$exps};
	my @lstRegMeta = @{$regs};
	my %expCtl;
	foreach my $exp (@lstExp){
		my @lstCtlExps;
		foreach my $regMeta (@lstRegMeta){

			if (uc $regMeta->get('experiment_accession') eq uc $exp){
				my @ctlIds = (split ",", $regMeta->get('control_id'));
				my $index = 0;
				my $numIds = @ctlIds;
				while ($index < $numIds){
					my $controlId = $ctlIds[$index];
					$controlId =~ s/^\s+|\s+$//g;
					my @lstExpCtls = get_expeiment_accession($controlId, \@lstRegMeta);
					push @lstCtlExps, @lstExpCtls;
					$index++;
				}
			}
		}
		#remove duplicates
		my @filtered = uniq(@lstCtlExps);
		$expCtl{$exp}=@filtered;
	}

	return \%expCtl;
	

}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}



sub get_expeiment_accession{
	my $control_id = shift;
	my $regs = shift;
	my @lstRegMeta = @{$regs};
	my @lstExp;
	foreach my $regMeta (@lstRegMeta){
		if (uc $regMeta->get('accession') eq uc $control_id){
			my $exp = $regMeta->get('experiment_accession');
			if (not grep(/^$exp$/, @lstExp)){
				push @lstExp, $exp;
			}

		}

	}
	return @lstExp;
}


sub get_experiment_objects{
	my $expAccession = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @lstSelected;
	foreach my $regMeta (@lstRegMeta){
		if ($regMeta->get('experiment_accession') eq $expAccession){
			push @lstSelected, $regMeta;
		}
	}
	return \@lstSelected;
}

sub get_multi_experiment_objects{
	my $expAccession = shift;
	my $lstObjs = shift;
	my @lstExpAccessions = @{$expAccession};
	my @lstRegMeta = @{$lstObjs};
	my @lstSelected;
	foreach my $regMeta (@lstRegMeta){
		my $exp = $regMeta->get('experiment_accession');
		if (grep(/^$exp$/, @lstExpAccessions)){
			push @lstSelected, $regMeta;
		}
	}
	return \@lstSelected;
}



=pod
sub get_num_multiple{
	my $experiment = shift;
	my $bioRep = shift;
	my $techRep = shift;
	my $paired_end = shift;
	my $multiple = shift;
	my %lstMultiple = %{$multiple};
	my $ret=0;
	
	my $myKey = $experiment.'_'.$bioRep.'_'.$techRep.'_'.$paired_end;
	#print "$myKey\n";
	print Dumper(%lstMultiple);
	if (exists $lstMultiple{$myKey}) {
		$ret = $lstMultiple{$myKey};
		$ret ++;
	}else{
		$ret = 1;
	}
	$lstMultiple{$myKey} = $ret;
	
	return $ret, \%lstMultiple;
}
=cut


sub get_control_info{
	my $experiments = shift;
	my $info;
	my $numExp = 0;
	foreach my $exp (@{$experiments}){
		if (index($exp->{'target'}->{'label'}, 'control') != -1) {
=pod
		if ((uc $exp->{'target'}->{'label'} eq 'CONTROL') || (uc $exp->{'target'}->{'label'} eq 'RABBIT-IGG-CONTROL') || (uc $exp->{'target'}->{'label'} eq '
MOUSE-IGG-CONTROL')){
=cut
			$numExp ++;
			$info .='Exp:'.$exp->{'accession'};
			my @files = @{$exp->{'files'}};
			foreach my $file (@files){
				if (uc $file->{'file_format'} eq 'FASTQ'){
					$info .=' File:'.$file->{'accession'}.'; ';
				}
			} 
			
		}
		
	}
	if ($numExp < 2){
		$info = '-';
	}
	return $info;
}

sub get_reference_epigenome_data{
	my $refEpiAccession = shift;
	my $url = 'https://www.encodeproject.org/reference-epigenomes/'.$refEpiAccession.'/?frame=embedded&format=json';
	my $headers = {Accept => 'application/json'};
	my $rEClient = REST::Client->new();
	$rEClient->GET($url, $headers);
	my $rEResponse = decode_json($rEClient->responseContent());
	return $rEResponse;
}


sub store_row{
	my $epi_accession = shift;
	my $accession= shift;
	my $expAccession= shift;
	my $epigenome = shift;
	my $featureType = shift;
	my $bioReplicate = shift;
	my $newBioReplicate = shift;
	my $techReplicate = shift;
	my $newTechReplicate = shift;
	my $gender = shift;
	my $md5Check = shift;
	my $localUrl = shift;
	my $analysis = shift;
	my $expGroup = shift; 
	my $assayXrefs = shift;
	my $ontXrefs = shift;
	my $xref = shift;
	my $epiDesc = shift;
	my $controlId = shift;
	my $paired = shift;
	my $paired_end_tag = shift;
	my $read_length = shift;
	my $multiple= shift;
	my $paired_with= shift;
	my $downUrl = shift; 
	my $info = shift;
	my $clsRegMeta = clsRegisterMetadata->new();
	$clsRegMeta->set_epi_accession($epi_accession);
	$clsRegMeta->set_accession($accession);
	$clsRegMeta->set_experiment_accession($expAccession);
	$clsRegMeta->set_epigenome($epigenome);
	$clsRegMeta->set_feature_type($featureType);
	$clsRegMeta->set_biological_replicate($bioReplicate);
	$clsRegMeta->set_new_bio_replicate($newBioReplicate);
	$clsRegMeta->set_technical_replicate($techReplicate);
	$clsRegMeta->set_new_tech_replicate($newTechReplicate);
	$clsRegMeta->set_gender($gender);
	$clsRegMeta->set_md5_checksum ($md5Check);
	$clsRegMeta->set_local_url ($localUrl);
	$clsRegMeta->set_analysis ($analysis);
	$clsRegMeta->set_experimental_group ($expGroup);
	$clsRegMeta->set_assay_xrefs ($assayXrefs);
	$clsRegMeta->set_ontology_xrefs ($ontXrefs);
	$clsRegMeta->set_xrefs ($xref);
	$clsRegMeta->set_epigenome_description ($epiDesc);
	$clsRegMeta->set_control_id ($controlId);
	$clsRegMeta->set_paired ($paired);
	$clsRegMeta->set_paired_end_tag ($paired_end_tag);
	$clsRegMeta->set_read_length ($read_length);
	$clsRegMeta->set_multiple ($multiple);
	$clsRegMeta->set_paired_with ($paired_with);
	$clsRegMeta->set_download_url ($downUrl);
	$clsRegMeta->set_info ($info);
	return $clsRegMeta;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: createRegMetaInputFile.pl -f <source_file> -p <local_store_path> -a <assay> -c <config_file> [-t <target>]

Options:
-f source_file:			file that contains the accessions of the Epigenome Summaries
-p local_store_path:	path where the files will be allocated
-a assay:				analysis description to filter experiments
-c config_file:			configuration file that contains the database connection details
-h help:				help message
[-t] target:			feature_type to filter by
 
END_USAGE

    say $usage;

    return 1;
}

sub fetch_adaptors {
    my ($cfg) = @_;
    my %adaptors;

    my $dbAdaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user    => $cfg->{efg_db}->{user},
        -host    => $cfg->{efg_db}->{host},
        -port    => $cfg->{efg_db}->{port},
        -dbname  => $cfg->{efg_db}->{dbname}
    );

    my $dbCoreAdaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user    => $cfg->{dna_db}->{user},
        -host    => $cfg->{dna_db}->{host},
        -port    => $cfg->{dna_db}->{port},
        -dbname  => $cfg->{dna_db}->{dbname}
    );

	
	$adaptors{analysis} = $dbAdaptor->get_adaptor("Analysis");
	$adaptors{FeatureType} = $dbAdaptor->get_adaptor("FeatureType");

	$adaptors{Gene} = $dbCoreAdaptor->get_adaptor("Gene");

	
	
    #my $dba = $dbAdaptor->db();

    #$adaptors{epigenome}    = $dba->get_EpigenomeAdaptor();
    #$adaptors{feature_type} = $dba->get_FeatureTypeAdaptor();
    #$adaptors{analysis}     = $dba->get_AnalysisAdaptor();

    return \%adaptors;
}

sub get_compare_hashes{
	my %adaptors = @_;
	my $lstAnalysis = $adaptors->{analysis}->fetch_all();
	my $lstFeatureTypes = $adaptors->{FeatureType}->fetch_all();
	
	my %hshAnalysis;
	foreach my $objAnalysis (@{$lstAnalysis}){
		my $logName = $objAnalysis->logic_name();
		$hshAnalysis{uc $logName}=$logName;
	}
	
	my %hshFeatureType;
	foreach my $objFeatureType (@{$lstFeatureTypes}){
		my $ftNamne = $objFeatureType->name;
		$hshFeatureType{uc $ftNamne}=$ftNamne;
	} 
	
	my %hshGender;
	$hshGender{'MALE'}='male';
	$hshGender{'FEMALE'}='female';
	$hshGender{'HERMAPHRODITE'}='hermaphrodite';
	$hshGender{'MIXED'}='mixed';
	$hshGender{'UNKNOWN'}='unknown';
	
	return \%hshAnalysis, \%hshFeatureType, \%hshGender;
}

sub check_db_value(){
	my $value= shift;
	my $hCompare = shift;
	my $ret;
	my %hshCompare = %{$hCompare};
	
	if ($hshCompare{uc $value}){
		#exists
		$ret = $hshCompare{uc $value};
	}
	return $ret;
	
}