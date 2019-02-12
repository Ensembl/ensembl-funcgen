#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;
use Carp;

my $cmd;

# run_ChromHMM_skip_nonsense_errors.pl --cmd "ls -lah"
# run_ChromHMM_skip_nonsense_errors.pl --cmd "java -Xmx30000m -jar /nfs/production/panda/ensembl/funcgen/picard/picard.jar"

my %config_hash = (
  'cmd' => \$cmd,
);

my $result = GetOptions(
  \%config_hash,
  'cmd=s',
);

my $non_issue_error_message = 'Exception in thread "main" java.lang.NumberFormatException: For input string:';

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;
$logger->info("Running: $cmd\n");

require Capture::Tiny;

my $return_value;

# Capture:Tiny has weird behavior if 'require'd instead of 'use'd
# see, for example,http://www.perlmonks.org/?node_id=870439 
my $stderr = Capture::Tiny::tee_stderr(sub {
    $return_value = system($cmd);
});

my $error_indicated = $return_value != 0;
my $error_is_non_issue = $stderr =~ /$non_issue_error_message/;

if ($error_indicated) {
    $logger->info("Not ok: The exit code was: $return_value\n");
} else {
    $logger->info("Ok: The exit code was: $return_value\n");
}

if ($error_indicated && $error_is_non_issue) {
    $logger->info("Ok: The error message indicates that this is not an issue.\n");
}

if ($error_indicated && !$error_is_non_issue) {
  $logger->info("Not ok: The error message indicates that this might be an issue.\n");
}

my $ChromHMM_exited_ok = 
        !$error_indicated
    || ( $error_indicated && $error_is_non_issue )
;

if (! $ChromHMM_exited_ok) {
    confess("Error running command\n" . $cmd . "\n");
}


$logger->info("Done.\n");
$logger->finish_log;
