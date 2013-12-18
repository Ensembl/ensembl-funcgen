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
  <helpdesk@ensembl.org>.

=head1 DESCRIPTION

Takes the elements from an arbitrary number of .bed files, as specified on the command line or as a list in a file, and clusters them into larger elements which are output in three column bed format. Optionally splits long elements into shorter ones. 

Clustering is by single linkage using a single base pair overlap. 

By default output is to STDOUT but an output file can be specified on the command line.

On request a distribution of the element lengths is printed to a separate file.

By default expects 0-based half-open UCSC convention coordinates but can work with 1 based coords if the -u flag is used.

=head1 USAGE

=head1 EXAMPLES

slink.pl -1 -o junk -d distrib -M4000 -m max all_dnase
slink.pl -1 -o junk -d distrib -M4000 -m max dnase_1 dnase_2

slink.pl -f Cmyc_files >> Cmyc.slc

slink.pl -o all_dnase.slc -d distrib all_dnase

slink.pl -1 -o junk -d distrib -M4000 -m max all_dnase

=head1 SEE ALSO


=head1 TO DO



=cut


use strict;
use Env;
use Getopt::Std;
use IO::Handle;
use IO::File;
#use Statistics::Descriptive;

my $outfile='';
my $list_file;
my $temp_dir = './';
my $verbose = 0;
my $distrib = '';
my $sort_mem = '';
my $sort_dir='';
my $one_based;
my $meth = 'none';
my $half_open = 1;
my $max_feature_length=100000000;
my @in_files;

my %opt;

if ($ARGV[0]){
    &Getopt::Std::getopts('d:H:M:m:o:f:T:U', \%opt) || die ;
}else{
    &help_text; 
}
&process_arguments;

if($list_file){
    @in_files = &get_lines_from_file($list_file);
}else{
    if(scalar @ARGV > 0){
	@in_files = @ARGV;
    }else{
	&help_text("you must supply at least one filename for input");
    }
}

if($distrib){
    require Statistics::Descriptive or die "Sorry you have requested a a distribution. Distributions need the perl module  Statistics::Descriptive but this does not appear to be present on your system";
}


my $ofh;
if($outfile){
     open($ofh,">$outfile") or die "couldn't open file $outfile\n";
}else{
        #or write to STDOUT
#     $ofh = new IO::File;
#     $ofh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}


#my $ifh;
#open($ifh,"$infile") or die "couldn't open file $infile\n";


#goto PROCESS; #debugging

my $together = $temp_dir."together_$$";
my $sorted =  $temp_dir."sorted_$$";
&backtick("rm -rf $together");
my $decr_clause = '';
unless($half_open){$decr_clause=" | awk '{print \$1\"\t\"\$2-1\"\t\"\$3}' "}

foreach my  $infile (@in_files){

    # copy the first three cols for each element file into the together 
    # file
    &backtick("cut -f1-3 ".$infile." |grep -v chrom >> $together");
}

&backtick("sort $sort_mem $sort_dir -k1,1 -k2,2g -k3,3g $together |grep -v '^#' ".$decr_clause." > $sorted");
unlink("$together");




PROCESS:


open(IN,"$sorted") or die "failed to open temporary file $sorted";
my $fh;
if($outfile){
    open($fh,">$outfile") or die "couldn't open file $outfile\n";
}else{
    $fh = new IO::File;
    $fh->fdopen(fileno(STDOUT), "w") || die "couldn't write to STDOUT\n";
}    

my @feats;
my $last_chr;
my $last_end=-1;
while(<IN>){
    #print $_ if $verbose > 1;
    chop;
    my @row = split("\t",$_);
    my $decimal_format = sprintf("%s\t%d\t%d",@row);
    @row = split("\t",$decimal_format);


    if(defined $last_chr && $row[0] ne $last_chr){
	print STDERR "chrom change \n" if $verbose > 1;
        &process_block($fh,\@feats,$max_feature_length,$meth,$last_end);
	$last_chr = $row[0];
	$last_end = $row[2];
	@feats = ();
        #print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
	next;
    }


    # if start of this feat is greater than last end its end of a block
    if($last_end > -1  && $row[1] > $last_end){
	print STDERR "gap\n" if $verbose > 1;
        &process_block($fh,\@feats,$max_feature_length,$meth,$last_end);
	$last_chr = $row[0];
	$last_end = $row[2];
	@feats = ();
        #print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
    }else{
        #print "adding to feats\n" if $verbose > 1;
	push @feats,\@row;
	$last_end = ($row[2] > $last_end)? $row[2]:$last_end;
    }
    
}
&process_block($fh,\@feats,$max_feature_length,$meth,$last_end);
close($fh);


unless($half_open){

    &backtick("cat $outfile | awk '{print \$1\"\t\"\$2+1\"\t\"\$3}' > $sorted ");
    &backtick("rm -f $outfile");
    &backtick("mv $sorted $outfile");
}



#unlink("$sorted");

if($distrib){
    &feature_distrib($outfile,$distrib, $max_feature_length);
}

###################################################################
sub feature_distrib{
    my($infile,$outfile, $max_feature_length)=@_;

    my $ifh;
    open($ifh,"$infile") or die "couldn't open file $infile\n";
    my $ofh;
    open($ofh,"> $outfile") or die "couldn't open file $outfile\n";
   

    &feature_stats($ifh,$ofh, $max_feature_length);

    close($ifh);
    close($ofh);
}

sub feature_stats{
    my($ifh,$ofh, $max_feature_length)=@_;

    my %further_checks;
    my $feature_count;
    my $max_len = -1;
    my $min_len = $max_feature_length;
    my %names;
    my @lengths;
    my $no_missing_vals = 1;
    my $features_without_val=0;
    my $overlapping_features;
    my $neg_length=0;
    my $last_end=0;
    my $last_chr='none';
    while(<$ifh>){
	
	$feature_count++;
        my ($name,$start,$end,@field) = split("\t", $_); 
	my $len = $end - $start;
	push @lengths,$len;
	$max_len = ($len > $max_len)? $len:$max_len;
	$min_len = ($len < $min_len)? $len:$min_len;
	if($len < 0){$neg_length += 1}
	$names{$name} += 1;

    	
    }

    

    print $ofh "INFO: feature_count $feature_count \n";

    print $ofh "INFO: min feature length $min_len\n";

    print $ofh "INFO: max feature length $max_len\n";

    if($distrib){
	print $ofh "INFO: Distribution of feature lengths \n";
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@lengths);
	my @bins = (1); # get count of 0 length features
	my $bl = int($max_len/50);
	for(my $i=100;$i<=4000;$i+=100){

	    push @bins,$i;
	}
	push @bins,$max_len;
	my %f = $stat->frequency_distribution(\@bins);
	for (sort {$a <=> $b} keys %f) {
	    print $ofh "up to $_ bp\t$f{$_}\t".
                  ( (int(100*$f{$_}/$feature_count))."%\n" );
	    #printf("<%6s %s\n",$_,('*' x int(100*$f{$_}/$feature_count) ));
	}
        # for pasting into excel
	for (sort {$a <=> $b} keys %f) {
	    print $ofh "$_\t$f{$_}\n";
	}

    }

    return \%further_checks;
}



sub process_block{
    my($fh,$feats_aaref,$max_len,$meth,$end) = @_;

    my $start = $feats_aaref->[0]->[1];
    my $len=$end -$start;

    if($len <1){return}  # filter out rubbish

    if($len <= $max_len ){
        print $fh join("\t",($feats_aaref->[0]->[0],
                             $feats_aaref->[0]->[1],$end))."\n";
    }
    elsif($meth =~ 'none'){
        print $fh join("\t",($feats_aaref->[0]->[0],
                             $feats_aaref->[0]->[1],$end))."\n";	
    }
    elsif($meth eq 'max'){
        &default_split($fh,$feats_aaref,$max_len,$end);
    }
    elsif($meth eq 'depth1_middle'){
        &depth1_middle($fh,$feats_aaref,$max_len,$end,$len);
    }
    elsif($meth =~ 'depth1_ends'){
	&depth1_ends($fh,$feats_aaref,$max_len,$end,$meth);
    }
    elsif($meth =~ 'depth_change'){
	&depth_change($fh,$feats_aaref,$max_len,$end);
    }
    elsif($meth =~ 'valley_percent'){
	&valley_percent($fh,$feats_aaref,$max_len,$end,$meth);
    }    
    else{
	die "unrecognised block splitting method $meth\n";
    }

}

sub depth_change{
    my($fh,$feats_aaref,$max_len,$end)=@_;


    my $start = $feats_aaref->[0]->[1];
    my $chr = $feats_aaref->[0]->[0];
 
    #print "length $len start $start\n";
    # NB &cover_depth changes the values in $feats_aaref
    my $depths_aaref = &cover_depth($feats_aaref);
    #my $n_segs = scalar(@$depths_aaref);

    foreach my $aref (@$depths_aaref){
	print $fh join("\t",$chr,$aref->[0],$aref->[1])."\n";
    }

}

sub valley_percent{
    my($fh,$feats_aaref,$max_len,$end,$meth)=@_;

    my($perc) = $meth =~ /valley_percent_([0-9]+)/;
    $perc = $perc/100;


    my $start = $feats_aaref->[0]->[1];
    my $chr = $feats_aaref->[0]->[0];
 
    #print "length $len start $start\n";
    # NB &cover_depth changes the values in $feats_aaref
    my $depths_aaref = &cover_depth($feats_aaref);
    my $n_segs = scalar(@$depths_aaref);


    my $peak_index = 0;
    my $peak_height=$depths_aaref->[0]->[3];
    my $peak_start = $depths_aaref->[0]->[0];
    my $found_val =0;
    for(my $i=0;$i<$n_segs;$i++){
	print "i = $i n_segs = $n_segs ".join("\t",@{$depths_aaref->[$i]})."\n";

	if($depths_aaref->[$i]->[3] >= $peak_height){
	    $peak_height = $depths_aaref->[$i]->[3];
	    $peak_index = $i;
	}

        if($depths_aaref->[$i]->[3] <= $peak_height * $perc ){
            # we've hit a valley
	    $found_val = 1;
            # find the lowest seg in this valley
	    my $val_depth = $depths_aaref->[$i]->[3];
	    my $val_start_index = $i;
	    my $val_end_index = $i;
	    while($depths_aaref->[$i]->[3] <= $peak_height * $perc){
		if($depths_aaref->[$i]->[3] < $val_depth){
		    $val_depth = $depths_aaref->[$i]->[3];
		    $val_start_index = $i;
		    $val_end_index = $i;
		}
		if($depths_aaref->[$i]->[3] == $val_depth){
		    $val_end_index = $i;
		}

	    $i++;


		if($i == $n_segs){last}
	print "i = $i n_segs = $n_segs ".join("\t",@{$depths_aaref->[$i]})."\n";
	    }
            # we are out of the valley or off the end of the block ie no valley
            if($i == $n_segs){
                #off the end of the block
		print $fh join("\t", $chr,$peak_start,
                                      $depths_aaref->[-1]->[1])."\n";

		print join("\t", $chr,$peak_start,
                                      $depths_aaref->[-1]->[1])."\n";
	    }else{
                my $val_mid = int(($depths_aaref->[$val_start_index]->[0] + 
                                   $depths_aaref->[$val_end_index]->[1])/2);
                print $fh join("\t", $chr,$peak_start,
                                      $val_mid)."\n";
                print  join("\t", $chr,$peak_start,
                                      $val_mid)."\n";


                # start a new peak
		$peak_start = $val_mid;
		$peak_height = $depths_aaref->[$i]->[3];
		$found_val = 0;
	    }
	}

    }
    # if there were no valleys or the last peak was not high enough to make
    # the end of the block look like a valley
    unless($found_val){
		print $fh join("\t", $chr,$peak_start,
                                      $depths_aaref->[-1]->[1])."\n";

		print join("\t", $chr,$peak_start,
                                      $depths_aaref->[-1]->[1])."\n";
    }

}



sub depth1_ends{
    my($fh,$feats_aaref,$max_len,$end,$len)=@_;

    my $start = $feats_aaref->[0]->[1];
    my $chr = $feats_aaref->[0]->[0];
    # get this now cos &cover_depth changes the values in $feats_aaref
    my $last = $feats_aaref->[-1]->[2];

    #print "length $len start $start\n";

    my $depths_aaref = &cover_depth($feats_aaref);
    my $n_segs = scalar(@$depths_aaref);

    if($n_segs ==1){
        print $fh $chr."\t".$start."\t".$depths_aaref->[0]->[1]."\n";
	return;
    }
    
    if($depths_aaref->[0]->[3] == 1){
        # if the first seg is depth 1 print it out
        print $fh $chr."\t".$start."\t".$depths_aaref->[0]->[1]."\n";
        # and print the start of the depth2+ seg
        print $fh $chr."\t".$depths_aaref->[0]->[1]."\t";
	
    }else{
        # only print the chr and start
        print $fh $chr."\t".$start."\t";
    }

    for(my $i = 1;$i<@$depths_aaref;$i++){ 
	#print join("\t", @$feat_ref)."\n";

        if($depths_aaref->[$i]->[3] == 1){

	    # this is end of depth2+ seg
	    my $new_end= $depths_aaref->[$i]->[0] ;
	    print $fh $new_end."\n";
	    # print entire depth1 seg
	    print $fh $chr."\t".$depths_aaref->[$i]->[0]."\t".
				$depths_aaref->[$i]->[1]."\n";

	    #if there are more segs
	    if($i < $n_segs-1){
		#print start of next depth2+ seg
		print $fh $chr."\t".$depths_aaref->[$i]->[1]."\t";
	    }else{
		#we've finished with a depth1 seg
		return;

	    }
	    
	}

    }
    # we finished in a depth2+ seg
    print $fh $depths_aaref->[-1]->[1]."\n"

}



sub depth1_middle{
    my($fh,$feats_aaref,$max_len,$end,$len)=@_;

    my $start = $feats_aaref->[0]->[1];
    my $chr = $feats_aaref->[0]->[0];
    # get this now cos &cover_depth changes the values in $feats_aaref
    #my $last = $feats_aaref->[-1]->[2]; # incorrect

    print "length $len start $start\n";

    my $depths_aaref = &cover_depth($feats_aaref);


    print $fh $chr."\t".$start."\t";

    foreach my $aref (@$depths_aaref){
	print join("\t", $chr,@$aref)."\n";
    }

    for(my $i = 1;$i<@$depths_aaref;$i++){ 
	#print join("\t", @{$depths_aaref->[$i]})."\n";

        # dont split last depth region
	if($i+1 == scalar(@$depths_aaref)){last}

        if($depths_aaref->[$i]->[3] == 1){
	    my $new_end= $depths_aaref->[$i]->[0] +
                         int($depths_aaref->[$i]->[2]/2);
	    print $fh $new_end."\n".$chr."\t".$new_end."\t";
	}

    }
    print $fh $depths_aaref->[-1]->[1]."\n";


}


# returns an array of arrays with each array containing
# start, end ,length and depth of each region with a given depth
sub cover_depth{
    my($feats_aaref) = @_;
    my @depth;

    my $start = $feats_aaref->[0]->[1];
    #print $start."\n";
    foreach my $feat_ref (@$feats_aaref){ 
	$feat_ref->[1] -= $start;
	$feat_ref->[2] -= $start;
	for(my $i=$feat_ref->[1];$i<$feat_ref->[2];$i++){
	    $depth[$i]++;
	}
    }

    my @segs;
    my $f_start = 0;
    my $f_end = -1;
    my $last_depth = $depth[0];
    my $i;
    for($i=0;$i<@depth;$i++){
	
        if($depth[$i] != $last_depth){
            $f_end = $i;
            my @arr = ($f_start+$start,$f_end+$start,
                       $f_end - $f_start,$last_depth);
	    push @segs,\@arr;
	    $last_depth = $depth[$i];
	    $f_start = $i; # half open
	}
    }
    my @arr = ($f_start+$start,$i+$start,$i - $f_start,$last_depth);
    push @segs,\@arr;


    return \@segs;

}






# not used
sub depth1{
    my($fh,$feats_aaref,$max_len,$end,$len)=@_;

    my $start = $feats_aaref->[0]->[1];
    my @depth;

    
    foreach my $feat_ref (@$feats_aaref){ 
	$feat_ref->[1] -= $start;
	$feat_ref->[2] -= $start;
	for(my $i=$feat_ref->[1];$i<$feat_ref->[2];$i++){
	    $depth[$i]++;
	}
    }

    print join("\n",@depth)."\n";
    my @segs;
    my $f_start = -1;
    my $f_end = -1;
    for(my $i=0;$i<$len;$i++){

	if($f_start == -1 && $depth[$i] == 1){
	    $f_start = $i;
	}

        if($f_start > -1 && $depth[$i] > 1){
            #end of run of 1s
	    $f_end = $i;
            print "$f_start, $f_end\n";
	    my @seg = ($f_start,$f_end);
	    push @segs, \@seg;
	    $f_start= -1;
	}

    }
    if($f_start > -1){
	$f_end = scalar(@depth);
        print "$f_start, $f_end\n";
	my @seg = ($f_start,$f_end);
        push @segs, \@seg;
    }
    
}


sub percent_split{
    my($fh,$feats_aaref,$max_len,$end,$meth)= @_;

    my($perc) = $meth =~ /([0-9]+)_percent/;
    my $start = $feats_aaref->[0]->[1];

# not finished
}


sub default_split{
    my($fh,$feats_aaref,$max_len,$end) = @_;

    my $start = $feats_aaref->[0]->[1];
    my $new_end = $start + $max_len;
    while($new_end < $end){
        print $fh join("\t",($feats_aaref->[0]->[0],
                             $start,$new_end))."\n";
	$start += $max_len;
        $new_end = $start + $max_len;
    }
    print $fh join("\t",($feats_aaref->[0]->[0],
                             $start,$end))."\n";

}



sub too_short_inputs{
    my($infile,$ofh,$need_method)=@_;

    my $ifh;
    open($ifh,"$infile") or die "couldn't open file $infile\n";
    my $feature_href = &feature_header($ifh,$ofh);
    my $col_attrib_aref = &data_header($ifh,$ofh,$need_method);
    close($ifh);

    #print join(' ',@$col_attrib_aref)." col\n";

    foreach my $col_href (@$col_attrib_aref){
	unless(defined $col_href){next} # first 4 elements are empty

	my($n_zero,$n_neg)=&get_too_short_counts($col_href,$feature_href);
	if($n_zero){
            print $ofh "INPUT DATA ERROR: FILE ".$col_href->{'#data'}.
                       " contains $n_zero zero length features\n";
	}

	if($n_neg){
            print $ofh "INPUT DATA ERROR: FILE ".$col_href->{'#data'}.
                       " contains $n_neg features with negative lengths\n";
	}
    }
}


sub get_too_short_counts{
    my($col_href,$feature_href)=@_;

    open(IN,$col_href->{'#data'}) or die "failed to open file ".$col_href->{'#data'};

    my $max = $feature_href->{max_feature_length};

    my $zero_count=0;
    my $neg_count = 0;
    my @examples;
    my $header_line = <IN>;
    while(<IN>){
	chop;
	my($chr,$start,$end,@field)=split("\t",$_);
	my $len = $end-$start;
	#if($len > $max){
	if($len == 0){
	    #push @examples, join('_',($chr,$start,$end));
	    $zero_count++;
	}
	if($len < 0){
	    #push @examples, join('_',($chr,$start,$end));
	    $neg_count++;
	}

    }
    close(IN);


    return($zero_count,$neg_count);
}




# returns a hashref indicating the type of problems encountered eg
# hash keys :-
# min_len_0 - there are features with length 0
# neg_lengths -  there are features with length < 0
sub feature_matrix{
    my($ifh,$ofh,$feature_href)=@_;

    my %further_checks;
    my $feature_count;
    my $max_len = -1;
    my $min_len = $feature_href->{max_feature_length};
    my %names;
    my @lengths;
    my $no_missing_vals = 1;
    my $features_without_val=0;
    my $overlapping_features;
    my $neg_length=0;
    my $last_end=0;
    my $last_chr='none';
    while(<$ifh>){
	
	$feature_count++;
        my ($name,$start,$end,@field) = split("\t", $_); 
	my $len = $end - $start;
	push @lengths,$len;
	$max_len = ($len > $max_len)? $len:$max_len;
	$min_len = ($len < $min_len)? $len:$min_len;
	if($len < 0){$neg_length += 1}
	$names{$name} += 1;

	if($_ =~ '.'){$no_missing_vals = 0}
	unless(join('',@field) =~ /[1-9]/){
	    $features_without_val++;
	    #print $_;
	}

	if($name eq $last_chr){
	    if($start < $last_end){ # half open
	        $overlapping_features = "feature $name $start $end overlaps ".
                    "preceeding feature which ends at $last_end\n";
	    }
	}else{
	    $last_chr = $name;
	}

	$last_end = $end;
      	
    }

    

    print $ofh "INFO: feature_count $feature_count \n";

    if($min_len >0){
        print $ofh "INFO: min feature length $min_len\n";
    }else{
	print $ofh "ERROR: min feature length $min_len\n";
	$further_checks{min_len_0} = 1;
    }

    if($max_len <= $feature_href->{max_feature_length}){
        print $ofh "INFO: max feature length $max_len\n";
    }else{
	print $ofh "ERROR: max feature length $max_len ".
        "should be at most ".$feature_href->{max_feature_length}."\n";
    }

    if($neg_length){
	print $ofh "ERROR: $neg_length features have negative lengths\n";
	$further_checks{neg_lengths} = 1;
    }

    print $ofh "INFO: Sequence Names and feature count for that sequence:-\n";
    print $ofh join(' ',%names)."\n";

    unless($no_missing_vals){
        print $ofh "WARNING: No features contain missing values ie '.' \n";
    }

    if($features_without_val){
        print $ofh "WARNING: $features_without_val features have no values above or below 0\n";
    }

    if($overlapping_features){
        print $ofh "WARNING: Some features overlap others eg.\n".
                 $overlapping_features;   
    }

    #print $ofh "INFO: \n";

    if($distrib){
	print $ofh "INFO: Distribution of feature lengths \n";
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@lengths);
	my @bins = (1); # get count of 0 length features
	my $bl = int($feature_href->{max_feature_length}/50);
	for(my $i=100;$i<=$feature_href->{max_feature_length};$i+=$bl){
	#for(my $i=100;$i<=5000;$i+=100){ #hack due to error in latest file
	    push @bins,$i;
	}
	my %f = $stat->frequency_distribution(\@bins);
	for (sort {$a <=> $b} keys %f) {
	    print $ofh "up to $_ bp\t$f{$_}\t".
                  ( (int(100*$f{$_}/$feature_count))."%\n" );
	    #printf("<%6s %s\n",$_,('*' x int(100*$f{$_}/$feature_count) ));
	}
    }
    return \%further_checks;
}


sub get_lines_from_file{
    my $file=shift;

    open(IN, "< $file") or die "couldn't open list file $file";

    my @ret;
    while( <IN> ){
        chop;
        push  @ret, $_ ;
    }

    return @ret;

}



sub backtick{
    my $command = shift;
#    $verbose = 2;
    warn "executing $command \n" if ($verbose ==2);

    my $res = `$command`;
    if($?){
        warn "failed to execute $command \n";
        warn "output :-\n $res \n";
        die "exit code $?";
    }

    return $res;
}



sub commentary{
    print STDERR "$_[0]";
}

   
sub err{
    print STDERR "$_[0]\n";
}
  


sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }


    if (exists $opt{o}){
        $outfile = $opt{o};
    }else{
        &help_text("Please supply the name of a file for output");
    } 

    if (exists $opt{U}){
        $half_open = 0;
    }

    if (exists $opt{f}){
        $list_file = $opt{f};
    }

    if  (exists $opt{m}){
        $meth = $opt{m};
    } 


    if  (exists $opt{M}){
        $max_feature_length = $opt{M};
	unless($opt{m}){
	    &help_text("If you specify a max feature length you must also specify a method for splitting features which are longer than this");
	}
    } 


    if  (exists $opt{d}){
        $distrib = $opt{d};
	unless($outfile){
            &help_text("use option -o if you want a distribution");
        }
    }

     if (exists $opt{T}){
        $temp_dir = $opt{T}.'/';
        $temp_dir =~ s/\/\//\//g; # change // to /
	$sort_dir = '-T '.$temp_dir;
    }


} 



sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    slink.pl -o <output_filename> [options] filename [filenames]

                  [-h] for help
                   -o <output file> - name of a file for output
                  [-f] <filename> - file containing a list of input filenames
                  [-M] <integer> - max length of output elements default:no max
                  [-d] <filename> for the distribution of feature lengths
                  [-m] <method> for splitting long elements eg 'depth1_middle',
                       'depth1_ends', 'max'
                       default = 'none'
                  [-T] <dir_path> temp directory. default= './'
                  [-U] flag - indicates coordsystem is 1-based (ie not UCSC)

END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
