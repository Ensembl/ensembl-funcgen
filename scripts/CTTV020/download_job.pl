#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Storable;
use Expect;

my $data_dir = $ENV{'DATA_DIR'};
my $job_index = $ENV{'LSB_JOBINDEX'};
my $user = $ENV{'USER'};
my $password = $ENV{'SANG_PASS'};

my @filenames = @{retrieve($data_dir . '/filenames_variable')};

my $filename = $filenames[$job_index - 1];
my ( $parent_dir, $junk ) = split /_/, $filename;


# deal with iRODS initialisation
my $session = Expect->spawn('kinit') or die $!;

if ($session->expect(5, "Password for $user\@INTERNAL.SANGER.AC.UK:")){
	print $session "$password\r";
}
else {
	say STDERR 'kinit failed';
}

$session->soft_close();


# run actual download command
my $iget_cmd   = 'iget -K /seq/' . $parent_dir . '/' . $filename;
my $exit_code = system($iget_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $iget_cmd;
}


