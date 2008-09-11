#!/usr/bin/env perl

=head1 NAME

write_pipeline_config.pl -- generates analysis pipeline config files

=head1 SYNOPSIS

write_pipeline_config.pl -module <name> [-module <module2>] (-file|-slice)

Mandatory:

   -module name    Name of Runnable, like SWEmbl or Nessie

   -slice|file     Select SubmitType, i.e. SubmitSlice or SubmitFile 
                   respectively

Optional:
   
   -regexp string  regular expression to select a subset of config keys
                   (paramter sets)

   -overwrite      Overwrite outfiles if they already exist

   -help|?         Print this 

=head1 DESCRIPTION

This script generates three config files in the $SCRIPTSDIR/conf directory 
from the Runnable analysis config modules specified on the commandline. 

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($host,$port,$user,$pass,$dbname,$species,$assembly,$help,
    $module,$modules,$file,$slice,$regexp,$overwrite);

GetOptions (
            'module=s@'    => \$modules, 
            'file'         => \$file,
            'slice'        => \$slice,
            'regexp=s'     => \$regexp, 
            'overwrite'    => \$overwrite,
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
} else {
    throw("Need to pass a Submit type (either -file or -slice).");
}
   
#if ( ! -e $analysis_conf ) {
#    $write_Submit = 1;
#}


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
             queue => 'normal', 
             batch_size => 10, 
             resources => 'select[type==X86_64 && myens_genomics1<80] rusage[myens_genomics1=10:duration=10]',
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
