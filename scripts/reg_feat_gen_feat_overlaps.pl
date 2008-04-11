#!/usr/local/ensembl/bin/perl -w

=head1 DESCRIPTION

Classifies regulatory features on the basis of their pattern of attributes and their overlap with various classes of genomic feature.

Assumes the genomic feature table name and the feature_type in the table are the same.

=head1 AUTHOR(S)

dkeefe@ebi.ac.uk

=head1 USAGE

set up the database for overlap analysis using scripts
reg_feats_4_classification.pl
and
gen_feats_4_classification.pl

mock_reg_feat_gen.pl -e dk_genomic_features_36k -i regulatory_features_filtered -o mockreg_features_filtered

reg_feat_gen_feat_overlaps.pl -e dk_genomic_features_36k -v1 -c /lustre/work1/ensembl/dkeefe/func_gen/reg_feat_gen_feat_conf.1    

=head1 EXAMPLES

 reg_feat_gen_feat_overlaps.pl -e dk_genomic_features_36k -v2  
        
=head1 SEE ALSO

=head1 TO DO

 Tidy Up
 POD
 Finish off.

=head1 CVS

 $Log: not supported by cvs2svn $







=cut


use strict;
use DBI;
use Env;
use Getopt::Std;
use IO::Handle;
use IO::File;
use lib '/nfs/acari/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
use DBSQL::Utils;
use DBSQL::DBS; # to use databases listed in ~/dbs.

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port);
my @temp_tables;
#my $scratch_dir = "/lustre/scratch1/ensembl/dkeefe/overlap_$$"."/";
#my $scratch_dir = "/lustre/scratch1/ensembl/dkeefe/overlap_14273"."/";
my $scratch_dir = "./junk"."/";

my $sp;
my $infile;
my $verbose = 0;
my $jump_label = '';
my $reg_feat_table = 'regulatory_features_filtered';
my $mock_reg_table = 'mockreg_features_filtered';
my $gen_feat_table = 'genomic_features';
my $overlaps_table = 'reg_feat_gen_feat_overlaps';
my $mock_olaps_table = 'mock_feat_gen_feat_overlaps';
my $overlap_results_file = $scratch_dir.'overlap_results.tab';
my $flags_table = 'regulatory_feature_association_flags';
my $types_table = 'regulatory_features_classified';
my $pat_bits = 32;
my $patt_count_thresh = 100;
my $assoc_thresh = 75;
my $second_thresh = 50;
my $rerun = 0;
my %opt;

if ($ARGV[0]){
&Getopt::Std::getopts('v:u:p:s:H:h:e:j:P:c:r', \%opt) || die ;
}else{
&help_text; 
}

# get configuration from environment variables
&config; # this may fail but config can be on command line


my $enc_db = 'dk_genomic_features_36k';# default, can be overridden by args
&process_arguments;

# hook up with the server
my $dbh = &make_contact($enc_db);
my $dbu = DBSQL::Utils->new($dbh);


#my $count =$dbu->get_count("select count(*) from  reg_feat_gen_feat_overlaps where pattern like '_______1________1_______________' ");
#print "$count<\n";



# we need access to the latest ensembl funcgen
#my $source_dbs = DBSQL::DBS->new('current_funcgen');
#my $source_dbh = $source_dbs->connect();
#my $source_dbu = DBSQL::Utils->new($source_dbh);

&source_checks($dbh,$dbu);

# get the list of gen feats in the analysis
# each item is the name of a table
my @gen_feats;
if($infile){
    @gen_feats = &get_lines_from_file($infile);
}else{
    @gen_feats=@ARGV;
}

if($jump_label){
    goto $jump_label;
}

# create a scratch directory
&backtick("rm -rf $scratch_dir");
&backtick("mkdir $scratch_dir");


# put all gen feats in one table and assign feature type ids
# ids begin at 0 and are in same order as  @gen_feats array
&create_gen_feats_table($dbh,$dbu,\@gen_feats,$gen_feat_table);

# get the list of seq_region_names
my $q = "select distinct seq_region_name from $gen_feat_table";
my $chrs_ref = $dbh->selectcol_arrayref($q);
unless(defined $chrs_ref && @$chrs_ref > 0){
    die "No seq_region_names extracted from $gen_feat_table using \n $q\n".$dbh->errstr;
}

# open a file for writing the results
my $reg_results_file = $scratch_dir.'reg_results.tab';
my $mock_results_file = $scratch_dir.'mock_results.tab';

my $reg_ofh;
my $mock_ofh;
open($reg_ofh,"> $reg_results_file") or die "couldn't open $reg_results_file";
open($mock_ofh,"> $mock_results_file") or die "couldn't open $mock_results_file";
# filename for temp storage of genomic features
my $gen_file = $scratch_dir.'genomic_features';

foreach my $chr (@$chrs_ref){

    &commentary("processing chromosome $chr\n") if $verbose;
    &write_gen_feats_file($chr,$gen_file);
    &overlap_analysis($dbh,$gen_file,$chr,$reg_feat_table,$reg_ofh,scalar(@gen_feats));
    &overlap_analysis($dbh,$gen_file,$chr,$mock_reg_table,$mock_ofh,scalar(@gen_feats));

}
close($reg_ofh);
close($mock_ofh);

# create a table for the overlap results and fill it from the file
# also adds the reg_feat attribute patterns to the table
# and at the moment cuts off the two gene bits
&commentary("loading overlap results to DB\n") if $verbose;
&load_results($dbh,$dbu,$reg_results_file,$overlaps_table,$reg_feat_table,\@gen_feats);
&load_results($dbh,$dbu,$mock_results_file,$mock_olaps_table,$mock_reg_table,\@gen_feats);






# get the patterns we wish to investigate
my @pats;
if($rerun){
    &commentary("reading bit patterns from overlap_summary_table\n") if $verbose;

    my $q = "select distinct pattern from pattern_overlap_summary";
    my $aref = $dbh->selectcol_arrayref($q);
    unless(defined $aref && @$aref > 0){
        die "failed to get any data from :\n$q\n".$dbh->errstr;
    }

    @pats = @$aref;

}else{
    &commentary("generating bit patterns\n") if $verbose;
    @pats = &combos_of_n_from_m($dbh,$dbu,7,$pat_bits,$overlaps_table);
    my $n_pats = scalar(@pats);
    # the only valid regfeats with just one bit set are the focus feature only ones
    # we know these are not associated with any particular genomic feature
    @pats = @pats[$pat_bits..($n_pats-1)];
}



#@pats = @pats[20000..($n_pats-1)];

#print join("\n",@pats[0..32]);

# open a file for the pattern overlap results

my $overlap_ofh;
open($overlap_ofh,"> $overlap_results_file") or die "couldn't open $overlap_results_file";

#for each pattern we want the number of reg feats with that pattern  
# and then for each genomic feature we want the overlap counts for the
# real and mock reg feats with that pattern
&commentary("calculating overlap counts\n") if $verbose;
my($reg_count,$real_count,$mock_count);
foreach my $pat (@pats){
    $reg_count = $dbu->get_count("select count(1) from  $overlaps_table where pattern like '$pat'");
 

    if($reg_count >= $patt_count_thresh){
        print $overlap_ofh "$pat\t$reg_count";
        foreach my $gen (@gen_feats){

            $real_count = $dbu->get_count("select count(1) from  $overlaps_table where pattern like '$pat' and $gen");
            $mock_count = $dbu->get_count("select count(1) from  $mock_olaps_table where pattern like '$pat' and $gen");

	    print $overlap_ofh "\t$real_count\t$mock_count";

        }
	print $overlap_ofh "\n";
    }
    
}
close($overlap_ofh); 


FINAL_LOAD:
#create a table for the overlap results
&commentary("importing pattern overlap counts to DB\n") if $verbose;
&create_pattern_results_tables($dbh,$dbu, $overlap_results_file,\@gen_feats);
&commentary("creating association tables\n") if $verbose;
&association_tables($dbh,$dbu,\@gen_feats,'pattern_overlap_summary','pattern_overlap_summary_chi',$assoc_thresh,$second_thresh);

FLAGS:
&create_flags_table($dbh,$dbu,$flags_table,\@gen_feats);

TYPES:
&create_types_table($dbh,$dbu,$flags_table,$types_table);



&clean_temp();
#&backtick("rm -rf $scratch_dir");
$dbh->disconnect;
exit;


 
###################################################################
sub create_types_table{
    my($dbh,$dbu,$flags_table,$types_table)=@_;

    my @sql;
    push @sql, "drop table if exists $types_table";
    push @sql, "create table $types_table select regulatory_feature_id,display_label,if(gm06990+cd4+imr90 = 1,1,0) as cell_type_specific,protein_coding_exon1_plus_enhancer as promoter_associated,protein_coding_gene_body as gene_associated, intergenic_2500 as non_gene_associated,0 as unclassified from $flags_table";

    # apply arbitrary rules to resolve conflicts
    push @sql, "update $types_table set promoter_associated = 0 where promoter_associated and gene_associated";
    # use unclassified col to flag conflicts
    push @sql, "update $types_table set unclassified = 1 where promoter_associated and non_gene_associated";
    push @sql, "update $types_table set unclassified = 1 where gene_associated and non_gene_associated";
    # set both of the conflicting cols to 0 where there is a conflict
    push @sql, "update $types_table set promoter_associated = 0 where unclassified";
    push @sql, "update $types_table set gene_associated = 0 where unclassified";    push @sql, "update $types_table set non_gene_associated = 0 where unclassified";
    # set unclassified for rows with no flags
   push @sql, "update $types_table set unclassified = 1 where gene_associated + non_gene_associated + promoter_associated = 0";


    &execute($dbh,@sql) or die;



    my $res = $dbu->get_count("select count(*) from  regulatory_features_classified where promoter_associated and cell_type_specific");
    &commentary("promoter_associated and cell_type_specific         $res\n");
    $res = $dbu->get_count(" select count(*) from  regulatory_features_classified where promoter_associated and not cell_type_specific");
    &commentary("promoter_associated and not cell_type_specific     $res\n");
    $res = $dbu->get_count(" select count(*) from  regulatory_features_classified where gene_associated and cell_type_specific");
    &commentary("gene_associated and cell_type_specific             $res\n");
    $res = $dbu->get_count(" select count(*) from  regulatory_features_classified where gene_associated and not cell_type_specific");
    &commentary("gene_associated and not cell_type_specific         $res\n");
    $res = $dbu->get_count(" select count(*) from  regulatory_features_classified where non_gene_associated and cell_type_specific");
    &commentary(" non_gene_associated and cell_type_specific        $res\n");
    $res = $dbu->get_count(" select count(*) from  regulatory_features_classified where non_gene_associated and not cell_type_specific");
    &commentary("non_gene_associated and not cell_type_specific     $res\n");


}


sub create_flags_table{
    my($dbh,$dbu,$flags_table,$gen_feat_ref)=@_;

    my @sql;
    push @sql,"drop table if exists $flags_table";
    push @sql, "create table $flags_table select regulatory_feature_id,display_label, 1 as cell_type_specific from regulatory_feature";
    push @sql, "alter table $flags_table add index(display_label)";
    &execute($dbh,@sql) or die;
    @sql = ();
    foreach my $feat (@$gen_feat_ref){
	&commentary("setting flags for $feat\n") if $verbose > 1;
        my $assoc_table = $feat."_assoc";
        my $not_assoc_table =  $feat."_not_assoc";

	unless($dbu->get_count("select count(*) from $assoc_table") > 0){
	    next;
	}

	push @sql, "alter table $flags_table add column $feat int(1) default 0";
	my $q = "select pattern from $assoc_table";
        my $pat_ref = $dbh->selectcol_arrayref($q);
	unless(defined $pat_ref){die "failed on \n$q\n".$dbh->errstr}
        foreach my $pat (@$pat_ref){
	    push @sql, "update $flags_table set $feat = 1 where display_label like '$pat'";
            if( ! &cell_type_specific_pattern($pat)){
                print "not specific\n";
                push @sql, "update $flags_table set cell_type_specific = 0 where display_label like '$pat'";
	    }

	}


	$q = "select pattern from $not_assoc_table";
        $pat_ref = $dbh->selectcol_arrayref($q);
	unless(defined $pat_ref){die "failed on \n$q\n".$dbh->errstr}
        foreach my $pat (@$pat_ref){
	    push @sql, "update $flags_table set $feat = 0 where display_label like '$pat'";
            push @sql, "update $flags_table set cell_type_specific = 1 where display_label like '$pat'";

	}


    }
    &execute($dbh,@sql) or die;

    # add cell type flags 
    @sql = ();
    push @sql, "alter table $flags_table add column cd4 int(1) default 0";
    push @sql, "update $flags_table set cd4 = 1 where display_label like '1%'";
    push @sql, "update $flags_table set cd4 = 1 where display_label like '_1%'";    push @sql, "alter table $flags_table add column gm06990 int(1) default 0";      push @sql, "update $flags_table set gm06990 = 1 where display_label like '________________________1%'";  
    push @sql, "alter table $flags_table add column imr90 int(1) default 0"; 
    push @sql, "update $flags_table set imr90 = 1 where display_label like   '_________________________1%'";
    push @sql, "alter table $flags_table add index(regulatory_feature_id)";
    &execute($dbh,@sql) or die;


}




sub cell_type_specific_pattern{
    my($pattern) = @_;

my @attrib = ('CD4_CTCF', 'CD4_DNASE_IMPORT', 'CD4_H2AZ', 'CD4_H2BK5me1', 'CD4_H3K27me1', 'CD4_H3K27me2', 'CD4_H3K27me3', 'CD4_H3K36me1', 'CD4_H3K36me3', 'CD4_H3K4me1', 'CD4_H3K4me2', 'CD4_H3K4me3', 'CD4_H3K79me1', 'CD4_H3K79me2', 'CD4_H3K79me3', 'CD4_H3K9me1', 'CD4_H3K9me2', 'CD4_H3K9me3', 'CD4_H3R2me1', 'CD4_H3R2me2', 'CD4_H4K20me1', 'CD4_H4K20me3', 'CD4_H4R3me2', 'CD4_PolII', 'GM06990_DNASE_IMPORT', 'Nessie_NG_STD_2_ctcf_ren_BR1', 'Wiggle_H3K27me3', 'Wiggle_H3K36me3', 'Wiggle_H3K4me3', 'Wiggle_H3K79me3', 'Wiggle_H3K9me3', 'Wiggle_H4K20me3');

# not cell type specific
#         1         2         3
#12345678901234567890123456789012
#1........................1......
#1.........................1.....
#.1.......................1......
#.1........................1.....
#.........................11.....

my %not_cell_spec = ('1........................1......', 1,
                     '1.........................1.....', 1,
                     '.1.......................1......', 1,
                     '.1........................1.....', 1,
                     '.........................11.....', 1);

    my $specific = 1;
    foreach my $pat (keys(%not_cell_spec)){
        # add the two extra 
	$pat .= '..';

	#print $pattern."\n$pat\n\n";
        if( $pattern =~ /$pat/ ){$specific = 0}


    }
    return $specific;
}



sub association_tables{
    my($dbh,$dbu,$gen_feat_ref,$summary_table,$chi_table,$assoc_thresh,$second_thresh)=@_;
    
    
    foreach my $feat (@$gen_feat_ref){
	# scan for initial associations
	my $q = "select perc.pattern, perc.$feat from $summary_table perc, $chi_table chi where perc.$feat >= $assoc_thresh and chi.$feat > 8 and perc.pattern = chi.pattern";
	&commentary($q) if $verbose > 1;
	my $aaref = $dbh->selectall_arrayref($q);
	unless(defined $aaref){die "failed on \n$q\n".$dbh->errstr}
	foreach my $aref (@$aaref){
	    push @$aref,$feat;
	    print join("\t",@$aref)."\n" if $verbose >1;
	}

	my @sql;
	my $assoc_table = $feat."_assoc";
	my $not_assoc_table =  $feat."_not_assoc";
	my $zero_one_table = $feat."_0_1";
	push @sql,"drop table if exists $assoc_table";
	push @sql,"create table $assoc_table select * from $summary_table where 1 = 0";
	push @sql,"drop table if exists $not_assoc_table";
	push @sql,"create table $not_assoc_table select * from $summary_table where 1 = 0";
     push @sql,"drop table if exists $zero_one_table";
	push @sql,"create table $zero_one_table select pattern from $summary_table where 1 = 0";

	&execute($dbh,@sql) or die;
	foreach my $aref (@$aaref){
	    my $pattern = $aref->[0];
	    print join("\t",@$aref)."\n" if $verbose >1;
	    &execute($dbh, "insert into $assoc_table select * from $summary_table where pattern like '$pattern'") or die;
	    # if none of the patterns have percent < $second_thresh then the original pattern
	    # is OK and it stays in the table.  we remove all the
	    # other patterns which have percent < $assoc_thresh
	    $q = "select pattern from $assoc_table where  $feat < $second_thresh";
	    my $not_ref = $dbh->selectcol_arrayref($q);
	    die $q."\n".$dbh->errstr unless(defined $not_ref);

	    if(@$not_ref > 0){
		foreach my $pat (@$not_ref){
		    print "not pattern $pat\n" if $verbose >1;
		    &execute($dbh,"insert into $not_assoc_table select * from $assoc_table where  pattern = '$pat'") or die;
		    &execute($dbh,"delete from $assoc_table where pattern = '$pat'") or die;
		    my $zero_one = &add_nots($pattern,$pat);
		    print "$zero_one\n" if $verbose >1;
		    &execute($dbh,"insert into $zero_one_table values('$zero_one')") or die;
		}
	    }
	}

	@sql = ();
        # the assoc table currently contains patterns which lie between
        # the assoc_thresh and the second_thresh
	push @sql,"delete from $assoc_table where $feat < $assoc_thresh";
	# distinct assoc table
	# **** and add the two gene bits to the patterns ****
	push @sql,"drop table if exists temp_$$";
	push @sql,"create table  temp_$$ select distinct * from $assoc_table";
	push @sql,"drop table $assoc_table";
	push @sql,"alter table temp_$$ rename as $assoc_table";
	push @sql,"update $assoc_table set pattern = concat(pattern,'__')";
	push @sql,"alter table $assoc_table add index(pattern)";

	# distinct not_assoc table
	# **** and add the two gene bits to the patterns ****
	# distinct not_assoc table
	# **** and add the two gene bits to the patterns ****
	push @sql,"drop table if exists temp_$$";
	push @sql,"create table  temp_$$ select distinct * from $not_assoc_table";
	push @sql,"drop table $not_assoc_table";
	push @sql,"alter table temp_$$ rename as $not_assoc_table";
	push @sql,"update $not_assoc_table set pattern = concat(pattern,'__')";
	push @sql,"alter table $assoc_table add index(pattern)";



	&execute($dbh,@sql) or die;
    }


}


sub add_nots{
    my($pat,$not) = @_;

    my @pat = split('',$pat);
    my @not = split('',$not);

    for(my $i = 0;$i<@pat;$i++){
        if($pat[$i] ne '1' && $not[$i] eq '1'){
            $pat[$i] = 0;
        }
    }
    return join('',@pat);
}


sub create_pattern_results_tables{
    my($dbh,$dbu, $overlap_results_file,$gen_feats_ref)=@_;

    my @sql;
    push @sql,"drop table if exists raw_overlap_results";

    my $q = "create table raw_overlap_results (pattern varchar(45),total_reg_feats int(10)";
    foreach my $gen (@$gen_feats_ref){
	$q .= ",$gen"."_real int(10), $gen"."_mock int(10)";
    }
    $q .= ")";
    push @sql,$q;    
    push @sql, "load data local infile '$overlap_results_file' into table raw_overlap_results";


    push @sql,"drop table if exists pattern_overlap_summary";
    push @sql,"drop table if exists pattern_overlap_summary_chi";
    $q = "create table pattern_overlap_summary (pattern varchar(45)";
    my $q1 = "create table pattern_overlap_summary_chi (pattern varchar(45)"; 
    my $q2 = "insert into pattern_overlap_summary select pattern";
    my $q3 = "insert into pattern_overlap_summary_chi select pattern";
    foreach my $gen (@$gen_feats_ref){
        # using int(10) to remove non integer part of %age 
	$q .= ",$gen int(10)";
	$q1 .= ",$gen float";
        $q2 .= ",100* $gen".'_real/total_reg_feats '."as $gen";
$q3 .= ",($gen"."_real- $gen"."_mock)*($gen"."_real- $gen"."_mock)/ $gen"."_mock as $gen"
    }
    $q .= ")";
    $q1 .= ")";
    $q2 .= " from raw_overlap_results";
    $q3 .= " from raw_overlap_results";
    push @sql,$q;   
    push @sql,$q1;   
    push @sql,$q2;
    push @sql,$q3;
    push @sql,"alter table pattern_overlap_summary_chi add index(pattern)";
    push @sql,"alter table pattern_overlap_summary add index(pattern)";

    &execute($dbh,@sql) or die;

}


sub combos_of_n_from_m{
    my($dbh,$dbu,$n,$len,$table) = @_;

    my $str = '_' x $len;
    my @blank = split('',$str);

    my @pat_arrs;
    push @pat_arrs, \@blank;

    my @pat_store;

    for(my $b = 1;$b <= $n;$b++){
    &commentary("creating patterns with $b bits set\n") if $verbose;
        my %new_pats;
        # take each of the patterns
        foreach my $aref (@pat_arrs){

            # and add an extra bit in each of the unnoccupied places
            for(my $i=0;$i<$len;$i++){
                unless($aref->[$i] eq '1'){
                    my @new = @$aref;
                    $new[$i] = 1;
                    $new_pats{ join('',@new) } =1;
                }
            }
        }

        # check that there are at least $patt_count_thresh reg feats with 
        # the new patterns
	if($b > 2 && $b != $n){
	    foreach my $pat (keys(%new_pats)){

		my $q =  "select count(*) from  $table where pattern like '$pat'";


		#my $reg_count = $dbu->get_count("select count(*) from  $table where pattern like '$pat'");
    #print "$q\n"; 

		my @reg_count = $dbh->selectrow_array("$q");
		#print $reg_count[0]."<\n";
		unless(defined $reg_count[0]){
		    die "failed on : $q\n".$dbh->errstr;
		}

		unless($reg_count[0] >= $patt_count_thresh){
		    #print "deleting\n";
		    delete $new_pats{$pat};
		}
	    }
        }

        push @pat_store,keys(%new_pats);
        @pat_arrs = ();
        foreach my $pat (keys(%new_pats)){
            my @arr = split('',$pat);
            push @pat_arrs , \@arr;
        }
	print $b.' '.scalar(@pat_store)."\n";
    }


    &commentary(join(" pat \n",@pat_store)." pat \n") if $verbose > 1;
    return @pat_store;
}


sub combos_of_n_from_m_simple{
    my($n,$len) = @_;

    my $str = '_' x $len;
    my @blank = split('',$str);

    my @pat_arrs;
    push @pat_arrs, \@blank;

    my @pat_store;

    for(my $b = 1;$b <= $n;$b++){
        my %new_pats;
        # take each of the patterns
        foreach my $aref (@pat_arrs){

            # and add an extra bit in each of the unnoccupied places
            for(my $i=0;$i<$len;$i++){
                unless($aref->[$i] eq '1'){
                    my @new = @$aref;
                    $new[$i] = 1;
                    $new_pats{ join('',@new) } =1;
                }
            }
        }

        push @pat_store,keys(%new_pats);
        @pat_arrs = ();
        foreach my $pat (keys(%new_pats)){
            my @arr = split('',$pat);
            push @pat_arrs , \@arr;
        }

    }


    print join("\n",@pat_store)."\n" if $verbose > 1;
    return @pat_store;
}


sub all_bitwise_pairs{
    my @bin_strings = @_;

    if(length( $bin_strings[0]) > 32){die "string too long for this algorithm"}

    my @dec_patterns;
    foreach my $str (@bin_strings){
       push @dec_patterns,&bin2dec($str);
    }

    my $npatterns = scalar(@dec_patterns);

    my %common_patterns_dec;
    for(my $i = 0;$i<$npatterns - 2;$i++){

       for(my $j=$i+1;$j<$npatterns-1;$j++){
	   my $common = $dec_patterns[$i]  & $dec_patterns[$j];
	    if($common){ # if the pattern is not all 0
		$common_patterns_dec{$common} ++;
	    }
       }
    }

    my @reg_exps;
    while(my ($dec,$count) = each(%common_patterns_dec)){
       push @reg_exps, &dec2bin($dec);
       #print  &dec2bin($dec)."\n";
    }

    return(\@reg_exps);
}



#convert a scalar integer to binary string
sub dec2bin {
   my $str = unpack("B32", pack("N", shift));
   #$str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
   return $str;
}


sub bin2dec {
   return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}





sub load_results{
    my($dbh,$dbu,$results_file,$overlaps_table,$reg_feat_table,$gen_feat_ref)=@_;

    my @sql;

    #create the summary table
    push @sql,"drop table if exists $overlaps_table";
    $q = "create table $overlaps_table (regulatory_feature_id int(10) unsigned, ";
    $q .= join(' int(1) unsigned, ',@$gen_feat_ref);
    $q .= " int(1) unsigned)";
    push @sql, $q;

    push @sql, "load data local infile '$results_file' into table $overlaps_table";
    push @sql, "alter table $overlaps_table add index(regulatory_feature_id)";
    my $temp1 = &new_temp();
    push @sql,"drop table if exists $temp1";
    #******************* Cuts end off pattern *******************************
    push @sql,"create table $temp1 select o.*,substring(r.display_label,1,32) as pattern from $reg_feat_table r, $overlaps_table o where o.regulatory_feature_id = r.regulatory_feature_id";
    push @sql,"drop table $overlaps_table";
    push @sql,"alter table $temp1 rename as $overlaps_table ";
    
    &execute($dbh,@sql) or die;
    

}

sub overlap_analysis{
    my($dbh,$gen_file,$chr,$reg_feat_table,$ofh,$columns)=@_;

    my $zeros = '0' x $columns;
    my @zero = split('',$zeros);

    my $q = "select regulatory_feature_id,feature_start,feature_end from $reg_feat_table where seq_region_name = '$chr' order by feature_start,feature_end";
    my $reg_aaref = $dbh->selectall_arrayref($q);
    unless(defined $reg_aaref){die "failed on \n $q\n".$dbh->errstr}

    # regulatory features may not be present on all chrs
    if($reg_aaref eq NO_ROWS){return}

    # create a hash to store the overlap flags
    my %flags;
    foreach my $reg_ref (@$reg_aaref){
        my @zero = split('',$zeros);
	$flags{$reg_ref->[0]} = \@zero;
    }


    open(IN,$gen_file) or die "failed to open $gen_file";

    while(<IN>){
	chop;
	my($type,$start,$end) = split("\t",$_);

        foreach  my $reg_ref (@$reg_aaref){

            # if the reg and gen feats overlap set the flag for that gen feat
            if($reg_ref->[1] <= $end && $reg_ref->[2] >= $start){
                $flags{$reg_ref->[0]}->[$type] += 1;
	    }



            # trim the regulatory feature array
            # remove any reg feat whose end < this gen feat start
	    if($reg_ref->[2] < $start){shift(@$reg_aaref)}

            # if the start of the reg feat is > end of the gen feat
            # we can move onto the next gen feat
            if($reg_ref->[1] > $end){last}
	}

	#print "length of reg_feat_array = ".scalar(@$reg_aaref)."\n";

    }

    while(my($id,$aref)=each(%flags)){
	print $ofh $id."\t".join("\t",@$aref)."\n";
    }


    
}





sub write_gen_feats_file{
    my($chr,$gen_file)=@_;

my $command = "mysql -h ".$host.
                 " -u ".$user.
                 " -P ".$port.
                 " -p".$password.
                 " ".$enc_db.
                 " -BN".
                 ' -e"'.
                 "select feature_type_id,feature_start,feature_end from  genomic_features where seq_region_name = '$chr' order by feature_start,feature_end".'"'.
                 " > $gen_file";

    &backtick($command);

}

sub create_gen_feats_table{
    my($dbh,$dbu,$gen_feats_aref,$table)=@_;

    my @sql;
    # create a temporary table of gen_feat,id
    my $temp1 = &new_temp();
    push @sql,"drop table if exists $temp1";
    push @sql,"create table $temp1 (feature_type varchar(45),feature_type_id int(10) unsigned)";
    for(my $i=0;$i < @$gen_feats_aref;$i++){
	push @sql, "insert into $temp1 values( '".$gen_feats_aref->[$i]."','".$i."')";
    }
    push @sql,"alter table $temp1 add index(feature_type)";

    # concatenate all the gen feats into one table
    my $temp2 = &new_temp();
    push @sql,"drop table if exists $temp2";
    push @sql, "create table $temp2 (feature_type varchar(45), seq_region_name varchar(40), feature_strand tinyint(2), feature_start int(10) unsigned, feature_end int(10) unsigned)";
    foreach my $feat (@$gen_feats_aref){
	push @sql, "insert into $temp2 select distinct feature_type, seq_region_name, feature_strand, feature_start, feature_end from $feat";
    }


    push @sql,"alter table $temp2 add index(feature_type)";


    push @sql,"drop table if exists $table";
    push @sql, "create table $table select t2.*,t1.feature_type_id from $temp1 t1, $temp2 t2 where t1.feature_type = t2.feature_type";
    push @sql,"alter table $table add index(seq_region_name)";

    #push @sql,"alter table $temp1 rename as lookup"; #debugging

    push @sql,"drop table if exists $temp1";
    push @sql,"drop table if exists $temp2";

    &execute($dbh,@sql) or die;


}


sub source_checks{
    my($dbh,$dbu) = @_;

    unless($dbu->table_exists('sources')){
	die "you must run the scripts gen_feats_4_classification.pl and reg_feats_4_classification.pl before this one";
    }

    my $q = "select db_name from sources where db_type = 'core'";
    my $core_db = $dbh->selectrow_array($q);

    $q = "select db_name from sources where db_type = 'funcgen'";
    my $func_db = $dbh->selectrow_array($q);


    my ($core_release) = $core_db =~ /.*_([0-9][0-9]_.*)/;
    my ($func_release) = $func_db =~ /.*_([0-9][0-9]_.*)/;

    unless($core_release eq $func_release){
	die "database versions don't match $core_db $func_db";
    }

        
}


sub get_lines_from_file{
    my $file=shift;

    open(IN, "< $file") or die "couldn't open config file $file";

    my @ret;
    while( <IN> ){
        if($_ =~ /^#/){next}
        chop;
        push  @ret, $_ ;
    }

    return @ret;

}



# uses global array @temp_tables
sub new_temp{

    my $nom = 'temp_'.$$.'_'.scalar(@temp_tables);
    push @temp_tables,$nom;
    return $nom;

}

# uses global array @temp_tables
sub clean_temp{

    my @sql;
    foreach my $table (@temp_tables){
	push @sql,"drop table if exists $table";
    }

    &execute($dbh,@sql) or die;
}

sub backtick{
    my $command = shift;

    &commentary( "executing $command \n") if $verbose >1;

    my $res = `$command`;
    if($?){
        warn "failed to execute $command \n";
        warn "output :-\n $res \n";
	die "exit code $?";
    }

    return $res;
}


# execute
#
#  Arg [1]   : scalar database handle from DBI module
#  Arg [2]   : array of scalars containing lines of text in the format of SQL statements
#  Function  : submits each element of Arg2 to the $dbh->prepare and 
#              $dbh->execute methods and checks for successful execution. 
#              Returns 0 on failure, 1 on success. Emits error messages 
#              via &err.
#  Returntype: int
#  Exceptions: none
#  Example   : execute($dbh, @array);

sub execute{
    my $dbh = shift;
    my (@array)=@_;

    my $sth;

    &commentary("processing SQL.") if $verbose;

    foreach my $query(@array){
	
    	&commentary("Executing  -  $query\n") if $verbose;

	    unless($sth=$dbh->prepare($query)){
               &err("Error: preparation of statement failed on: $query\n");
               &err("database error message: ".$dbh->errstr."\n");
               return(0);
	    }

        unless($sth->execute){ # returns true on success
            &err("Error: statement execution failed on: $query\n");
            &err("statement handle error message:".$sth->errstr."\n");
            return(0);
	    }
    }

    &commentary("\n");

    return(1);
}


# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print  STDERR "$_[0]";
}



sub config{
 
($user =     $ENV{'ENSFGUSER'}) or return(0); # ecs1dadmin
($password = $ENV{'ENSFGPWD'}) or return(0); #
($host   =   $ENV{'ENSFGHOST'}) or return(0); #localhost
($port =     $ENV{'ENSFGPORT'}) or return(0); #3360
($driver  =  $ENV{'ENSFGDRIVER'}) or return(0); #mysql
 
}
   
sub err{
    print STDERR "$_[0]\n";
}
  


sub make_contact{
    my $db = shift;

    unless($driver && $db && $host && $port && $user){
	&err("DB connection parameters not set");
        exit(1);
    }

    # hook up with the server
    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port;mysql_local_infile=1";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 0,PrintError=>0});
    if ($dbh){
        print STDERR ("connected to $db\n");
    }else{
        &err("failed to connect to database $db");
        exit(1);
 
    }   
            
    return $dbh;
}



sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }

    if (exists $opt{H}){
        $host = $opt{H}; 
    }

    if (exists $opt{p}){
        $password = $opt{p}; 
    }

    if (exists $opt{P}){
        $port = $opt{P}; 
    }

    if (exists $opt{u}){
        $user = $opt{u}; 
    }

    $sp = 'homo_sapiens';
    if (exists $opt{s}){
        $sp = lc($opt{s});
    }

    if (exists $opt{v}){
        $verbose = $opt{v};
    }


    if  (exists $opt{c}){
        $infile = $opt{c};
    } 


    if (exists $opt{j}){
        $jump_label = $opt{j} ;
    }

    if  (exists $opt{e}){
        $enc_db = $opt{e};
    } 

    if  (exists $opt{r}){
        $rerun = 1;
    } 

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    .pl [-h] for help
                  [-e] use a particular local encode catalog
                  [-H] <host machine> eg ecs2
                  [-u] <database user>
                  [-j] <label> jump to label then start execution
                  [-p] <mysql password> 
                  [-P] <mysql port> 
                  [-s] <species> eg -smus_musculus, default = homo_sapiens
                  [-v] <integer> verbosity level 0,1 or 2 
                  [-c] <filename> list of genomic feature tables
                  [-r] flag - indicates a re-run on same set of reg feats
                       as last time. Patterns can be taken  from the existing
                       pattern_overlap_summary table
                  [-] 
                  [-] <> 
                  [-] <> 


END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
