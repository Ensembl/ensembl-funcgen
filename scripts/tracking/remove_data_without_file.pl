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
use List::MoreUtils qw(uniq);

my $filename = $ARGV[0];

if (!$filename) {
	usage();
	exit;
}


#=====================================
#Load data
#=====================================
my @data_set = load_data($filename);

#==============================================
# Remove rows with missing files
#==============================================
@data_set = remove_data_without_files(\@data_set);

#=================================
#delete controls without signal
#=================================
@data_set = delete_controls_without_signal(\@data_set);

#==================================
#delete signals without controls
#==================================
@data_set = delete_signal_without_controls(\@data_set);

#==================================
#check paired integrity
#==================================
@data_set = remove_paired_signal_and_controls_with_missing_parts(\@data_set);

#===============================
#numerate replicates
#===============================
@data_set = check_replicate_numbers_and_multiple_column(\@data_set);

#===============================
#Print
#===============================
foreach my $final_row (@data_set){
	print $final_row->csv_row."\n";
}

#==============================
#END
#==============================
#print "FIN\n";

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

sub delete_signal_without_controls{
	my $data_set = shift;
	my @lstRegMeta = @{$data_set};
	my @accessions_to_remove;
	foreach my $row (@lstRegMeta){
		if (uc $row->get('epi_accession') ne 'EPI_ACCESSION'){
			if ($row->get('feature_type') ne 'WCE' && $row->get('control_id') ne '-'){
				my $is_controlled = check_signal($row->get('control_id'), \@lstRegMeta);

				if ($is_controlled == 0){
					if (uc $row->get('paired') eq "YES"){
						$is_controlled = check_signal($row->get('paired_with'), \@lstRegMeta);
					}
					if ($is_controlled == 0){
						my @multiple = get_multiple($row, \@lstRegMeta);
						if (@multiple){
							foreach my $multiple_file (@multiple){
								$is_controlled = check_signal($multiple_file, \@lstRegMeta);
								if ($is_controlled == 1){
									last;
								}
							}
						}
					}

					if ($is_controlled == 0){
						push @accessions_to_remove, $row->get('accession');
					}

				}
 
			}
		}

	}

	my @uniq_data_to_remove = uniq @accessions_to_remove;

	my @final_data;
	foreach my $row (@lstRegMeta){
		my $accession = $row->get('accession');
		if (not (grep(/^$accession$/, @uniq_data_to_remove))){
			push  @final_data, $row;
		}
	}

	return @final_data;
}


sub delete_controls_without_signal{
	my $data_set = shift;
	my @lstRegMeta = @{$data_set};
	my @accessions_to_remove;
	foreach my $row (@lstRegMeta){
		if ($row->get('feature_type') eq 'WCE'){
			my $is_controlling = check_control($row->get('accession'), \@lstRegMeta);

			if ($is_controlling == 0){

				if (uc $row->get('paired') eq "YES"){
					#check if the paired controls a signal
					$is_controlling = check_control($row->get('paired_with'), \@lstRegMeta);
				}

				if ($is_controlling == 0){
					#check if any of the multiple files controls a signal
					my @multiple = get_multiple($row, \@lstRegMeta);
					if (@multiple){
						foreach my $multiple_file (@multiple){
							$is_controlling = check_control($multiple_file, \@lstRegMeta);
							if ($is_controlling == 1){
								last;
							}
						}
					}

				}

				if ($is_controlling == 0){
					push @accessions_to_remove, $row->get('accession');
				}				
			}

		}
	}

	my @uniq_data_to_remove = uniq @accessions_to_remove;

	my @final_data;
	foreach my $row (@lstRegMeta){
		my $accession = $row->get('accession');
		if (not (grep(/^$accession$/, @uniq_data_to_remove))){
			push  @final_data, $row;
		}
	}

	return @final_data;

}

sub remove_data_without_files{
	my $data_set = shift;
	my @lstRegMeta = @{$data_set};
	my @data_to_remove;
	foreach my $row (@lstRegMeta){
		if (uc $row->get('epi_accession') eq "EPI_ACCESSION"){
			next;
		}
		
		my $file_path =$row->get('local_url');
		if (!(-e $file_path)) {
			push @data_to_remove, $row->get('accession');

			#check paired files
			if (uc $row->get('paired') eq "YES"){
				push @data_to_remove, $row->get('paired_with');
			}

			#check multiple files
			my @multiple = get_multiple($row, \@lstRegMeta);
			if (@multiple){
				push @data_to_remove, @multiple;
			}

		}
	}

	my @uniq_data_to_remove = uniq @data_to_remove;

	my @final_data;
	foreach my $row (@lstRegMeta){
		my $accession = $row->get('accession');
		if (not (grep(/^$accession$/, @uniq_data_to_remove))){
			push @final_data, $row;
		}
		
	}
	return @final_data;
}

sub get_multiple {
	my $source_row = shift;
	my $lstObjs = shift;
	my @lstRegMeta = @{$lstObjs};
	my @multiple_rows;

	foreach my $obj (@lstRegMeta){

		if ($obj->get('accession') ne $source_row->get('accession')){
			if ($obj->get('experiment_accession') eq $source_row->get('experiment_accession')){
				if ($obj->get('new_bio_replicate') eq $source_row->get('new_bio_replicate')){
					if ($obj->get('new_tech_replicate') eq $source_row->get('new_tech_replicate')){
						if ($obj->get('paired_end_tag') eq $source_row->get('paired_end_tag')){
							push @multiple_rows, $obj->get('accession');
						}
					}
				}
			}
		}
	}

	return @multiple_rows;
}

sub load_data{
	my $filename = shift;
	my @data_set;
	open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
	while (my $row = <$fh>) {
	    chomp $row;

	    #split each line into array
	    my @line = split(/\t/,$row);

	    my $obj_data = store_row(\@line);

	    push @data_set, $obj_data;

	 }

	 close($fh);

	 return @data_set;

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


sub check_control{
	my $ctl_accession = shift;
	my $data = shift;
	my $found = 0;

	FOO: {
		foreach my $ctl (@{$data}){
			if ($ctl->get('feature_type') ne 'WCE' && ($ctl->get('epi_accession') ne 'epi_accession')){
				my @ctl_ids = split(/\,/,$ctl->get('control_id'));
				foreach my $control_by (@ctl_ids){
					if ($ctl_accession eq $control_by){
						$found = 1;
						last FOO;
					}
				} 
			}
		}
	}

	return $found;
}

sub check_signal{
	my $ctl_accession = shift;
	my $data = shift;
	my $found = 0;

	my @ctl_ids = split(/\,/,$ctl_accession);

	foreach my $control_by (@ctl_ids){
		foreach my $ctl (@{$data}){
			if ($ctl->get('feature_type') eq 'WCE'){
				if ($control_by eq $ctl->{'accession'}){
					$found = 1;
					last;
				}
			}
		}
		if ($found == 0){
			#print "Not found: ".$control_by."\n";
			last;
		}
	}

	return $found;

}




sub store_row{
	my $line = shift;

	my $epi_accession = @{$line}[0];
	my $accession= @{$line}[1];
	my $expAccession= @{$line}[2];
	my $epigenome = @{$line}[3];
	my $featureType = @{$line}[4];
	my $bioReplicate = @{$line}[5];
	my $newBioReplicate = @{$line}[6];
	my $techReplicate = @{$line}[7];
	my $newTechReplicate = @{$line}[8];
	my $gender = @{$line}[9];
	my $md5Check = @{$line}[10];
	my $localUrl = @{$line}[11];
	my $analysis = @{$line}[12];
	my $expGroup = @{$line}[13]; 
	my $assayXrefs = @{$line}[14];
	my $ontXrefs = @{$line}[15];
	my $xref = @{$line}[16];
	my $epiDesc = @{$line}[17];
	my $controlId = @{$line}[18];
	my $paired = @{$line}[19];
	my $paired_end_tag = @{$line}[20];
	my $read_length = @{$line}[21];
	my $multiple= @{$line}[22];
	my $paired_with= @{$line}[23];
	my $downUrl = @{$line}[24]; 
	my $info = @{$line}[25];
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

Usage: remove_data_without_file.pl input_file 

Options:
input_file: Input file with the data to register 

 
END_USAGE

    say $usage;

    return 1;
}



