#!/usr/bin/env perl
=head1

=head2 Description

  Checks that a bam file has a valid end of file marker.
  
  This script it required by the ERSA pipeline.
  
  Based on advice pasted in some internet forum:
  
  http://seqanswers.com/forums/showthread.php?t=15363

=head2 Example

perl check_bam_file_has_EOF_marker.pl \
  --bam_file /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/dbfiles//homo_sapiens/GRCh38/funcgen/alignment/091/ersa_signal/bam/ENCODE/GM12878/DNase1/GM12878_DNase1_ENCODE_TR2_no_duplicates.bam

=cut
use strict;
use Getopt::Long;

my %options;
GetOptions (
    \%options,
    "bam_file|s=s",
);

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $bam_file = $options{'bam_file'};

if (! -e $bam_file) {
  die("Can't find $bam_file!");
}

use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( run_cmd );

my $cmd = qq(tail -c 28 $bam_file | hexdump -x);
my $stdout = run_cmd($cmd);

my $expected =<<EOF_MARKER
0000000    8b1f    0408    0000    0000    ff00    0006    4342    0002
0000010    001b    0003    0000    0000    0000    0000                
000001c
EOF_MARKER
;

if ($stdout ne $expected) {
  die (
    "$bam_file doesn't match!\n"
    . "Got:\n\n"
    . $stdout . "\n"
  );
}

print "ok\n";
