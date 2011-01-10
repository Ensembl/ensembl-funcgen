#!/usr/bin/env perl

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


=head1 NAME

write_pipeline_config.pl -- generates analysis pipeline config files

=head1 SYNOPSIS

write_pipeline_config.pl -module <name> [-module <module2>] (-file|-slice)

Mandatory:

   -module name      Name of Runnable, like SWEmbl or Nessie

   -slice|file       Select SubmitType, i.e. SubmitSlice or SubmitFile 
                     respectively

Optional:
   
   -regexp <string>  Regular expression to select a subset of config keys
                     (paramter sets)

   -overwrite        Overwrite outfiles if they already exist

   -queue <string>   Set name of queue (default: normal)
   -batch_size <int> Set size of batch (default: 10)

   -help|?           Print this 

=head1 DESCRIPTION

This script generates three config files in the $SCRIPTSDIR/conf directory 
from the Runnable analysis config modules specified on the commandline. 

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($host,$port,$user,$pass,$dbname,$species,$assembly,$help,
    $module,$modules,$file,$slice,$array,$regexp,$overwrite);

### sone defaults for BatchQueue conf
my $queue = 'normal';
my $batch_size = 10;

GetOptions 
	(
	 'module=s@'    => \$modules, 
	 'file'         => \$file,
	 'slice'        => \$slice,
	 'array'        => \$array,
	 'regexp=s'     => \$regexp, 
	 'overwrite'    => \$overwrite,
	 
	 'queue=s'      => \$queue,
	 'batch_size=i' => \$batch_size,
	 
	 'help|?'       => \$help,
	);

system("perldoc $0") if (defined $help);


$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

throw("Need to pass at least one runnable/module name (-module)")
    unless (defined $modules);

my $Submit;
if ($slice) {
    $Submit = 'Slice';
} elsif ($file) {
    $Submit = 'File';
} elsif ($array) {
    $Submit = 'Array';
} else {
    throw("Need to pass a Submit type (either -file or -slice).");
}
   
#if ( ! -e $analysis_conf ) {
#    $write_Submit = 1;
#}

my $dbhost;
if ($ENV{EFG_HOST} eq "127.0.0.1") {
	$dbhost = 'localhost';
} else {
	($dbhost = $ENV{EFG_HOST}) =~ s,-,_,g;
}

# create outfiles 
my $analysis_conf = $ENV{ANALYSIS_CONFIG};
throw("File '$analysis_conf' already exists. (Re)move this file ".
      "or specify option -overwrite on the command line.")
    if (-e $analysis_conf && !defined $overwrite);

open(CONF, "> $analysis_conf")
    or throw ("Can't open file $analysis_conf");

# Write SubmitType
print CONF 
    "[Submit$Submit]\n".
    "input_id_type=$Submit\n\n";


my $rules_conf = $ENV{RULES_CONFIG};
throw("File '$rules_conf' already exists. (Re)move this file ".
      "or specify option -overwrite on the command line.")
    if (-e $rules_conf && !defined $overwrite);

open(RULES, "> $rules_conf")
    or throw ("Can't open file $rules_conf");


my $batchq_conf = $ENV{BATCHQ_CONFIG};
throw("File '$batchq_conf' already exists. (Re)move this file ".
      "or specify option -overwrite on the command line.")
    if (-e $batchq_conf && !defined $overwrite);

open(BATCHQ, "> $batchq_conf")
    or throw ("Can't open file $batchq_conf");


# write analysis config
foreach my $module (@$modules) {

    # Read at run time Runnable config for specified module from 
    # Bio::EnsEMBL::Analysis::Config::Funcgen::${module}.pm

    my $config_module = "Bio::EnsEMBL::Analysis::Config::Funcgen::${module}";
    warn("Reading $config_module");

    eval "require $config_module";
    throw("Couldn't require $config_module") if ($@);
    $config_module->import(qw(CONFIG));

    my $CONFIG = $main::CONFIG;
    #warn("\nConfig read from $config_module\n".Dumper $CONFIG);

    foreach (sort keys %{$CONFIG}) {
        
        next if m/^DEFAULT$/;
        next if (defined $regexp && ! m/$regexp/);
        #print Dumper $_;

        my $program = $CONFIG->{$_}->{'PROGRAM'} || $CONFIG->{'DEFAULT'}->{'PROGRAM'};
        my $module = $CONFIG->{$_}->{'MODULE'} || $CONFIG->{'DEFAULT'}->{'MODULE'};
        
        print CONF 
            "[$_]\n".
            "program=$program\n".
            "module=$module\n".
            "input_id_type=$Submit\n\n";
        
        print RULES
            "[$_]\n".
            "condition=Submit$Submit\n\n";
        
        
        print BATCHQ "
            {
             logic_name => '$_',
             queue => '$queue', 
             batch_size => $batch_size, 
             resources => 'select[type==X86_64 && my$dbhost<80] rusage[my$dbhost=10:duration=10]',
             retries => 3,
             runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Funcgen',
             cleanup => 'no',
            },\n";
   
    }
}

close CONF;
close RULES;
close BATCHQ;

warn("\nDon't forget to copy and paste the config in $batchq_conf into ".
     "your Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm. Here add these ".
     "lines to the QUEUE_CONFIG list in the Config hash.\n");

1;
