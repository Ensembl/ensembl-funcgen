#!/usr/bin/env perl

# ----------------------------------------------------------------------------
# This script retrieves and stores the meta data for each .cram file that is
# related to the study of interest, ie. CTTV020 epigenomes of cell lines PILOT.
# ----------------------------------------------------------------------------

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

# -------------------------------------------------------------------
# Read the list that contains the files related to our study.
# The list is the ouput of command
# imeta qu -z seq -d study = 'CTTV020 epigenomes of cell lines PILOT'
# and target = 1 and manual_qc = 1 >filteredFileList'
# -------------------------------------------------------------------
open my $file_list, '<', $ARGV[0];

my ( $dir, $filename, $extension, $command );

$/ = "----\n";    # read entry by entry, not line by line

while (<$file_list>) {
    if (/^collection: (\/seq\/\S+)\ndataObj: (\S+)\.(\S+)\n/) {
        $dir       = $1;
        $filename  = $2;
        $extension = $3;
    }
    else {
        next;
    }

    my $cram_file_path = $dir . '/' . $filename . '.' . $extension;
    my $meta_file_path = $ARGV[1] . '/' . $filename . '.meta';

    if ( !-e $meta_file_path ) {

        # retrieve and store meta data for each cram file
        $command = 'imeta ls -d ' . $cram_file_path . ' >' . $meta_file_path;
        my $exit_code = system($command);
        if ( $exit_code != 0 ) {
            say 'This system call has failed: ' . $command;
        }
    }

}

close($file_list);