#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

sam2bed.pl

=head1 SYNOPSIS

sam2bed.pl [ -on_farm -one_based | -man | -help ] -files (file.sam[.gz])+

=head1 DESCRIPTION

This script converts sam (inc gzipped) format files to bed format.
Unmapped reads are filtered as appropriate and the MAPQ score is used as the bed score.

=head1 OPTIONS

=over

=item B<-on_farm>

Flag to submit each file as a separate job to the farm via bsub

=item B<-one_based>

Flag to over-ride default 0-based (half open) bed format start loci i.e. forces 1 based, ensembl style start loci

=item B<-files>

Takes a list of space separated sam files. These can be gzipped i.e. have a suffix of '.gz'

=item B<-man>

Prints the full manual page

=item B<-help>

Gives a short help message

=back

=cut


use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_file_format is_gzipped open_file);


# To do
# Utilise Bio::DB::Sam
# Add filters on quality etc. See IDs code
# Add GetOpt::long support:
#     -No zip output?
#     -submit to farm

#It works using bwa output. Not sure if it is fully SAM compliant!
#This does not use the Binary format!

my ($on_farm, @files, $one_based);

my @tmp_args = @ARGV;

GetOptions 
  (
   'on_farm'    => \$on_farm,
   'files=s{,}' => \@files,
   '1_based'    => \$one_based,         
   'man'        => sub { pod2usage(-exitval => 0, -verbose => 2); },
   'help'       => sub { pod2usage(-exitval => 0, -verbose => 1, -message => "Params are:\t@tmp_args"); }
  ) or pod2usage( -exitval => 1,
                  -message => "Params are:\t@tmp_args"
                );


print "sam2bed.pl @tmp_args\n" if $on_farm;

#Not implemented option functions yet

if (! @files){
  pod2usage( -exitval => 1,
			 -message => "You must supply some sam filenames to convert");
}


if((scalar(@files) > 1) &&
   ! $on_farm){
  die("You must supply the -on_farm flag if you want to run on more than one file");
  
}

my (@non_sams, $file_operator, $gz, $infile, $outfile, $regex);

foreach my $file(@files){
  
  #Let's die before we submit to farm

  if(&get_file_format($file) ne 'sam'){
    push @non_sams, $file;
  }

  if($on_farm){
    my $bsub_cmd = "bsub -q normal -J sam2bed:$file -o ${file}.sam2bed.out -e ${file}.sam2bed.err ".$ENV{'EFG_SRC'}."/scripts/miscellaneous/sam2bed.pl -files $file";
    system($bsub_cmd);
    next;
  }


  $gz = (&is_gzipped($file)) ? '.gz' : '';
  $file_operator = ($gz) ? "gzip -dc %s |" : '';
  
  my $infile = open_file($file, $file_operator);
  my $outfile = $file;
  my $regex = ($gz) ? '\.sam\.gz' : '\.sam$';
  $outfile =~ s/${regex}/\.bed/;
  $outfile .= '.gz';
  $outfile = open_file($outfile, '| gzip -c > %s');

  print "Converting file to bed format:\t$file\n";
  my ($strand, @cache, $name, $flag, $slice_name, $pos, $mapq);
  my ($read, $start, $seq_region_name);
   

  while(<$infile>){
    next if /^@/; #Skip header
    chomp;

    #(Query,flag,ref,pos,mapq,cigar,mrnm,mpos,isize,seq,quak,opt)
    ($name, $flag, $slice_name, $pos, $mapq, undef, undef, undef, undef, $read) = split("\t");
    next if $flag & 4; #Unmapped read.
	
    #Query strand is in the 5th(index == 4) bit...  assuming the reference strand never changes
    #i.e. 2*4 = 16 bitwise 16 & 16 == 16 else 0
    $strand = ($flag & 16) ? '-' : '+';
    (undef, undef , $seq_region_name) = split(':', $slice_name);
	
    #Can we put something better than 100 in score position?
    #Do we really want the names in the bed file?
    #SEQ_REGION_NAME, START, END, FEATURE_NAME, SCORE STRAND
    $start = ($one_based) ? $pos : ($pos - 1);
    push @cache, join("\t", ($seq_region_name, $start, ($pos + length($read) -1), $name, $mapq, $strand));


    if(scalar(@cache) == 1000){
      print $outfile join("\n", @cache)."\n";
      @cache = ();
    }
  }

  print $outfile join("\n", @cache)."\n";
  
  close $infile;
  close $outfile;
}


if(@non_sams){
  my $all_or_some = 'Some';

  if(scalar(@non_sams) == scalar(@files)){
    $all_or_some = 'All';
  }
  
  pod2usage(-exitval => 1,
            -message => "$all_or_some of your specified files are not sam format:\n\t".
            join("\n\t", @non_sams));
}

1;
