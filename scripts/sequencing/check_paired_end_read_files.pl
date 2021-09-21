#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 405ress or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 

=head1 SYNOPSIS

=head1 DESCRIPTION

paired_end_execution_plan_ids="155 156 161 170 187 194 201 207 212 263 279 295 296 328 329 340 341 358 401 488 621"

for id in $paired_end_execution_plan_ids
do
    bsub perl scripts/sequencing/check_paired_end_read_files.pl \
        --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
        --species mus_musculus \
        --execution_plan_id $id
done

perl scripts/sequencing/check_paired_end_read_files.pl \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
  --species mus_musculus \
  --execution_plan_id 155

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);
use Carp;
use Bio::EnsEMBL::Utils::Logger;

my $registry;
my $species;
my $execution_plan_id;

my %config_hash = (
  'registry'      => \$registry,
  'species'       => \$species,
  'execution_plan_id'  => \$execution_plan_id,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
  'execution_plan_id=s',
);

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $execution_plan_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'executionplan');

$logger->info("Fetching execution plan\n");

my $execution_plan_obj = $execution_plan_adaptor->fetch_by_dbID($execution_plan_id);

if (! defined $execution_plan_obj) {
  $logger->error("Can't find peak execution plan with id $execution_plan_id!");
  $logger->error(Dumper($execution_plan_adaptor));
  
  die;
}

my $execution_plan = $execution_plan_obj->execution_plan_deserialised;

lock_execution_plan($execution_plan);

print Dumper($execution_plan);

my $complete_alignment_names = get_all_complete_alignment_names($execution_plan);

print Dumper($complete_alignment_names);

my $all_paired_end_read_files = get_all_paired_end_read_files($execution_plan, $complete_alignment_names);

print (Dumper($all_paired_end_read_files));

my $read_file_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'readfile');

READ_FILE:
foreach my $all_paired_end_read_file (@$all_paired_end_read_files) {

    if ($all_paired_end_read_file->{type} eq SINGLE_END) {
        next READ_FILE;
    }

    my $read_file_name_1 = $all_paired_end_read_file->{1};
    my $read_file_name_2 = $all_paired_end_read_file->{2};
    
    my $read_file_1 = $read_file_adaptor->fetch_by_name($read_file_name_1);
    my $read_file_2 = $read_file_adaptor->fetch_by_name($read_file_name_2);
    
    $logger->info($read_file_1->file . " " . $read_file_2->file . "\n" );
    check_number_of_reads_identical ($read_file_1, $read_file_2);
}

READ_FILE:
foreach my $all_paired_end_read_file (@$all_paired_end_read_files) {

    if ($all_paired_end_read_file->{type} eq SINGLE_END) {
        next READ_FILE;
    }

    my $read_file_name_1 = $all_paired_end_read_file->{1};
    my $read_file_name_2 = $all_paired_end_read_file->{2};
    
    my $read_file_1 = $read_file_adaptor->fetch_by_name($read_file_name_1);
    my $read_file_2 = $read_file_adaptor->fetch_by_name($read_file_name_2);
    
    $logger->info($read_file_1->file . " " . $read_file_2->file . "\n" );
    my $error_message = check_read_names_match_up       ($read_file_1, $read_file_2);
    
    if (@$error_message) {
        
        my $read_file_notes = join "\n", @$error_message;
        
        $read_file_1->notes($read_file_notes);
        $read_file_2->notes($read_file_notes);
        
        $read_file_adaptor->update($read_file_1);
        $read_file_adaptor->update($read_file_2);
    }
}

$logger->finish_log;

sub check_read_names_match_up {

  my @read_file_objects = @_;

  use Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser;
  my $parser = Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser->new;
  
  my @fastq_files   = map { $_->file                                        } @read_file_objects;
  my @cmds          = map { "zcat $_ |"                                     } @fastq_files;
  my @file_handles  = map { open my $fh, $_ or die("Can't execute $_"); $fh } @cmds;

  use List::Util qw( none any uniq );
  my $all_records_read = undef;
  
  my @error_message;
  my $errors_found         =  0;
  my $max_errors_to_report = 10;
  
  PAIR_OF_READS:
  while (
      (! $all_records_read)
    ) {
    
    my @fastq_records = map  { $parser->parse_next_record($_) } @file_handles;
    $all_records_read = none { defined $_ } @fastq_records;
    
    if ($all_records_read) {
        last PAIR_OF_READS;
    }
    my @read_names = map {
        
        # @SOLEXA2_0007:7:1:0:17#0/1
        #
        my $name_found = $_->[0] =~ /^(\@.+?)[ |\/]/;
        
        if (!$name_found) {
            die("Can't parse name in $_->[0]");
        }
        my $read_name = $1;
        $read_name;
    
    } @fastq_records;

    my $reads_match = uniq(@read_names) == 1;
    
    if (! $reads_match) {
        
        my $error_message = "Names don't match up: " . Dumper(\@read_names);
        
        warn $error_message;
        
        push @error_message, $error_message;
        $errors_found++;
        
        if ($errors_found >= $max_errors_to_report) {
            last PAIR_OF_READS;
        }
    }
    next PAIR_OF_READS;
  }
  map { $_->close or die "Can't close file!" } @file_handles;
  return \@error_message;
}

sub check_number_of_reads_identical {

  my @read_file_objects = @_;

  use List::Util qw( uniq );
  my $same_number_of_reads_in_every_file 
    = 1 == uniq map { $_->number_of_reads } @read_file_objects;
  
  if (! $same_number_of_reads_in_every_file) {
    confess(
      "These read files are paired, but don't have the same number of reads:\n"
      . "\n"
      . ( join 
            "\n    and\n", 
            map { 
                "    read file name: "  . $_->name            . ", "
                  . "number of reads: " . $_->number_of_reads . ", "
                  . "file size: "       . $_->file_size       . ", "
                  . "file: "            . $_->file
            } @read_file_objects 
        )
      . "\n"
      . "\n"
    );
  }
}

sub get_all_paired_end_read_files {

    my $execution_plan  = shift;
    my $alignment_names = shift;

    my @all_paired_end_read_files;

    foreach my $complete_alignment_name (@$alignment_names) {

        my $alignment_plan = $execution_plan->{alignment}->{$complete_alignment_name};
        my $read_files = $alignment_plan->{input}->{read_files};
        
        my @paired_end_read_files = map {
            $_->{type} eq PAIRED_END; $_
        } @$read_files;
        
        push @all_paired_end_read_files, @paired_end_read_files;
    }
    return \@all_paired_end_read_files;
}

sub get_all_complete_alignment_names {

    my $execution_plan = shift;
    
    my $alignment_plans = $execution_plan->{alignment};
    my $alignment_names = [ keys %$alignment_plans ];

    my @complete_alignment_names;

    PLANNED_ALIGNMENT:
    foreach my $alignment_name (@$alignment_names) {

        my $alignment_plan = $alignment_plans->{$alignment_name};
        
        if ($alignment_plan->{name} eq NO_CONTROL_FLAG) {
            next PLANNED_ALIGNMENT;
        }
        if ($alignment_plan->{analysis} ne ALIGNMENT_ANALYSIS) {
            next PLANNED_ALIGNMENT;
        }
        if ($alignment_plan->{is_complete} ne TRUE) {
            next PLANNED_ALIGNMENT;
        }
        push @complete_alignment_names, $alignment_name;
    }
    return \@complete_alignment_names

}
