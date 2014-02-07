#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

Simple overlap utility which takes two files of tab separated values,
each line corresponding to a feature. Each file must have as its first
three columns <region identifier> <start-coord> <end_coord>. The
region identifiers must be in the same format in both files. Any
number of subsequent columns can be present in either file.

The files are sorted into temporary files and then an overlap analysis
is done. If two features overlap one another then both lines are
output on a single tab-separated line.

By default output appears on standard output. The -o option allows an
output file to be specified.

The standard output can be piped into other unix utilities for instance

  cut or awk to reduce or re-order the fields in the output. 

  awk to see if the two features are on the same strand.eg if the fourth
  field in each file is the strand specifier, awk '{if($4 == $8) print}'
  only prints the features which are on the same strand.

  uniq -u to get the non-overlapping cases eg
  file a.tab has 4 columns
  cooccur.pl a.tab b.tab |cut -f1-4 > a_overlapping_b.tab
  cat a.tab  a_overlapping_b.tab |sort|uniq -u > a_NOT_overlapping_b.tab


For very big files the sort may need more temporary space than is
available in the /tmp directory. cooccur.pl can use a different
temporary directory specified with the -T option.

If you are doing a lot of work on the same very big files it is worth creating sorted versions of them to speed things up. It is vital that both files are sorted in the same way... typically using 

  sort -k1,1 -k2,2g -k3,3g < unsorted_file > sorted_file

Then use the -m flag with cooccur.pl to indicate that the sort utility can simply merge the files.

If you are dealing with .bed files or other files which use the UCSC half-open coordinate convention (ie sequence start = 0) then specify this with the -u flag. Both files need to use the same coordinate convention.

=head1 USAGE

=head1 EXAMPLES

mysql -uensro -h ensdb-1-12 -P3304 -BN -e"select sr.name,g.seq_region_start,g.seq_region_end,x.display_label,g.seq_region_strand from  gene g,seq_region sr, xref x where sr.seq_region_id = g.seq_region_id and x.xref_id = g.display_xref_id and g.biotype in ('protein_coding','V_segment','J_segment','C_segment','D_segment') " homo_sapiens_core_47_36i > protein_coding_genes.tab

mysql -uensro -h ensdb-1-12 -P3304 -BN -e"select sr.name,g.seq_region_start,g.seq_region_end,x.display_label from  gene g,seq_region sr, xref x where sr.seq_region_id = g.seq_region_id and x.xref_id = g.display_xref_id and g.biotype in ('miRNA') " homo_sapiens_core_47_36i > micro_RNA_genes.tab

 cooccur.pl micro_RNA_genes.tab protein_coding_genes.tab |cut -f4,8

gives a two column list showing which micro_RNAs occur within the boundaries of which protein coding gene.


mysql -uensro -h ensdb-1-12 -P3304 -BN -e"select sr.name,if(g.seq_region_strand = 1,g.seq_region_start-500,g.seq_region_end+1) as promoter_start, if(g.seq_region_strand = 1,g.seq_region_start-1,g.seq_region_end+500) as promoter_end,x.display_label,g.seq_region_strand from  gene g,seq_region sr, xref x where sr.seq_region_id = g.seq_region_id and x.xref_id = g.display_xref_id and g.biotype in ('protein_coding','V_segment','J_segment','C_segment','D_segment') " homo_sapiens_core_47_36i >protein_coding_promoters.tab

alternatively

cat protein_coding_genes.tab| awk '{if($5 == 1){ print $1"\t"$2-500"\t"$2-1"\t"$4"\t"$5}else{print $1"\t"$3+1"\t"$3+500"\t"$4"\t"$5} }' > protein_coding_promoters.tab


 cooccur.pl micro_RNA_genes.tab protein_coding_promoters.tab  | wc

Only 13 micro_RNA genes lie inside the 500bp promoter regions of protein coding genes.
 
  mysql -uensro -h ensdb-archive -P5304 -BN -e"select sr.name,r.seq_region_start,r.seq_region_end,rc.repeat_name,rc.repeat_class,rc.repeat_type,a.logic_name from  repeat_feature r,seq_region sr, repeat_consensus rc, analysis a where sr.seq_region_id = r.seq_region_id and r.repeat_consensus_id = rc.repeat_consensus_id and r.analysis_id = a.analysis_id  " homo_sapiens_core_50_36l > all_repeats.tab
  grep RepeatMask all_repeats.tab > RepeatMask.tab
  grep TRF all_repeats.tab > TRF.tab
  grep Dust all_repeats.tab > Dust.tab

  cooccur.pl Dust.tab TRF.tab > potentially_redundant_annotation.tab

The above may take a minute or so to run.

=head1 SEE ALSO

=cut


use strict;
use Env;
use Getopt::Std;
use IO::Handle;
use IO::File;

my $outfile='';
my $verbose = 0; # for debugging
my $fh;
my $sort_mem;
my $sort_str='';
my $sort_m = '';
my $sort_dir = '';
my $sort_merge = '';
#my $flip=0;
my $half_open=0;

my %opt;

# if this script is called from a script which uses the DBI module it breaks.
# The DBI module somehow causes $SIG{'PIPE'} to be set to 'IGNORE' in the
# child shell. This causes grep or awk to ignore the interrupt when piped
# into head -1
$SIG{'PIPE'}='DEFAULT';


if ($ARGV[0]){
&Getopt::Std::getopts('umh:o:S:T:', \%opt) || &help_text("Invalid argument");
}else{
&help_text; 
}

if($verbose){ $| = 1; } #no output buffer

&process_arguments;
unless(@ARGV == 2){&help_text("You must give two filenames")}
my($file1,$file2)= @ARGV;

#precalculate array index for print out
#my $fields0 = &backtick("grep -v '^#' $file1 |head -n 1");
my $fields0 = &backtick("awk /^[^#]/  $file1 |head -n 1");
$fields0 = $fields0 =~ tr/\t/\t/; # count the tabs
#my $fields1 = &backtick("grep -v '^#' $file2 |head -n 1");
my $fields1 = &backtick("awk /^[^#]/ $file2 |head -n 1");
$fields1 = $fields1 =~ tr/\t/\t/; 

my $ext = time."_$$";
if($sort_dir){$sort_str = '-T '.$sort_dir}
if($sort_mem){$sort_m = '-S '.$sort_mem}

if($sort_merge){
    &backtick("awk '{print \$0\"\t\"0}' < $file1 > $sort_dir"."together0_$ext");
    &backtick("awk '{print \$0\"\t\"1}' < $file2 > $sort_dir"."together1_$ext");
    &backtick("sort $sort_merge $sort_m $sort_str -k1,1 -k2,2g -k3,3g $sort_dir"."together0_$ext $sort_dir"."together1_$ext |grep -v '^#'  > $sort_dir"."sorted_$ext");
    unlink("$sort_dir"."together0_$ext");
    unlink("$sort_dir"."together1_$ext");
}else{
# add a file identifier to each row then sort 
    &backtick("awk '{print \$0\"\t\"0}' < $file1 > $sort_dir"."together_$ext");
    &backtick("awk '{print \$0\"\t\"1}' < $file2 >> $sort_dir"."together_$ext");
    &backtick("sort $sort_m $sort_str -k1,1 -k2,2g -k3,3g $sort_dir"."together_$ext |grep -v '^#'  > $sort_dir"."sorted_$ext");
    unlink("$sort_dir"."together_$ext");
}

# we add one here and subtract one on printout
if($half_open){
    # add one to the start coord of all lines
    &backtick("awk '{print \$1\"\t\"\$2+1\"\t\"\$0}' < $sort_dir"."sorted_$ext |cut -f 1,2,5- > $sort_dir"."together_$ext");
    &backtick("rm -f $sort_dir"."sorted_$ext");
    &backtick("mv $sort_dir"."together_$ext $sort_dir"."sorted_$ext");
}


open(IN,"$sort_dir"."sorted_$ext") or die "failed to open temporary file $sort_dir"."sorted_$ext";

if($outfile){
     open($fh,">$outfile") or die "couldn't open file $outfile\n";
}else{
        #or write to STDOUT
     $fh = new IO::File;
     $fh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}

my @feats;
my $last_chr = '';
my $last_end= -1;
my $mixed = 0; 
my $last_set = undef;
while(<IN>){
    print $_ if $verbose > 1;
    chop;
    my @row = split("\t",$_);


    if($row[0] ne $last_chr){
	print "chrom change \n" if $verbose > 1;
        if($mixed){&process_block($fh,\@feats,$fields0,$fields1,$half_open)} 
	$last_chr = $row[0];
	$last_end = $row[2];
        $mixed = 0;
        $last_set = $row[-1];
	@feats = ();
        print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
	next;
    }

    # if start of this feat is greater than last end its end of a block
    if($row[1] > $last_end){
	print "gap\n" if $verbose > 1;
        if($mixed){&process_block($fh,\@feats,$fields0,$fields1,$half_open)} 
	$last_chr = $row[0];
	$last_end = $row[2];
        $mixed = 0;
        $last_set = $row[-1];
	@feats = ();
        print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
    }else{
        print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
	$last_end = ($row[2] > $last_end)? $row[2]:$last_end;
        if(defined $last_set && $row[-1] != $last_set){$mixed = 1}
    }
    print "last_set = $last_set mixed = $mixed\n" if $verbose > 1;
}
if($mixed){&process_block($fh,\@feats,$fields0,$fields1,$half_open)}


unlink("$sort_dir"."sorted_$ext");

exit;


 
###################################################################
sub process_block{
    my($fh,$feats_ref,$fields0,$fields1,$half_open ) = @_;


    my @set1;
    my @set0;

    print "block\n" if $verbose > 1;
    foreach my $aref (@$feats_ref){
	if($aref->[-1] == 0){
	    push @set0,$aref;
	}else{
	    push @set1,$aref;
	}
	print "\t". join("\t",@$aref)."\n" if $verbose > 1;
    }

		

    foreach my $ref0 (@set0){
	foreach my $ref1 (@set1){

	    if($ref0->[1] <= $ref1->[2] && $ref0->[2] >= $ref1->[1]){

		if($half_open){
                    print $fh $ref0->[0]."\t".($ref0->[1]-1)."\t".
			join("\t",@$ref0[2..$fields0])."\t";
                    print $fh $ref1->[0]."\t".($ref1->[1]-1)."\t".
			join("\t",@$ref1[2..$fields1])."\n";
		}else{

		    print $fh join("\t",@$ref0[0..$fields0],
                                    @$ref1[0..$fields1])."\n" or die 
                                    "printing out failed";
		}
	    }
	}
    }
}


sub backtick{
    my $command = shift;
#    my $verbose=2;
    warn "executing $command \n" if ($verbose ==2);

    my $res = `$command`;

    # grep gives an error if it finds nothing
    # this is not always what we want it to do
    if($command =~ 'grep' && $? == 256){
#	return '';
    }


    if($?){
        warn "failed to execute $command \n";
        warn "output :-\n $res \n";
	die "exit code $?";
    }

    return $res;
}


sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }

    if (exists $opt{u}){
        $half_open = $opt{u};
    }   

    if (exists $opt{o}){
        $outfile = $opt{o};
    }  

    if (exists $opt{T}){
        $sort_dir = $opt{T}.'/';
        $sort_dir =~ s /\/\//\//;# change // to /	
    }

    if (exists $opt{S}){
	$sort_mem = $opt{S};
    }

    if  (exists $opt{m}){
        $sort_merge = '-m';
    } 

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    cooccur.pl [options]  <file1> <file2>

                  [-h] for help

                  [-o] <file> - name of a file for output. default=STDOUT
                  [-u] flag indicating coordinates are UCSC style half-open
                  [-m] <flag> - indicating both files are sorted in the same way
                  [-S] <string> - memory use option for sort command eg '10%'
                  [-T] <string> - directory for writing temporary files 
                                  eg '/nfs/misc/encode/tmp'

END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
