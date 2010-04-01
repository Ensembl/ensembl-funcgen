#!/usr/local/ensembl/bin/perl -w

=head1 DESCRIPTION

Copies, via dumps to a scratch directory, and then denormalises data from an ensembl core database and creates tables of genomic features in a second database for use in classifying regulatory features.

Tables all have a standard set of columns.

=head1 AUTHOR(S)

dkeefe@ebi.ac.uk

=head1 USAGE


At the mysql command prompt, create a database for the classification run.

  create database dk_reg_feat_classify_49

edit the file ~dkeefe/dbs/current_core or  ~dkeefe/dbs/current_mouse_core to point at the latest ensembl core database

  gen_feats_4_classification.pl -e dk_reg_feat_classify_49 -s mus_musculus


=head1 EXAMPLES

=head1 SEE ALSO

=head1 TO DO

 Tidy up.
 POD

=head1 CVS



 $Log: not supported by cvs2svn $
 Revision 1.6  2010-02-12 10:29:53  dkeefe
 check in test

 Revision 1.5  2009-08-04 08:26:59  dkeefe
 mods for farm2

 Revision 1.4  2009/06/03 10:18:20  dkeefe
 created separate tables for each RNA type

 Revision 1.3  2009/05/26 10:46:48  dkeefe
 lib directory name change

 Revision 1.2  2008/08/20 08:54:18  dkeefe
 tidied and more pod

 Revision 1.1  2008/07/28 14:20:33  dkeefe
 moved from scripts directory.
 IG biotypes updated.

 Revision 1.1  2008/04/11 10:29:18  dkeefe
 Copies and then denormalises data from an ensembl core database and creates tables of genomic features in a second database for use in classifying regulatory features.





=cut


use strict;
use DBI;
use Env;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Std;
use IO::Handle;
use IO::File;
use lib '/nfs/users/nfs_d/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
use DBSQL::Utils;
use DBSQL::DBS; # to use databases listed in ~/dbs.

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port);
my @temp_tables;
my $dump_dir = "/lustre/scratch103/ensembl/dkeefe/gen_feat_junk_$$/";
my $slim_table = 'goslim_goa_acc_list';
my $id_list;
my $sp = 'homo_sapiens';
my $jump_label = '';
my $verbose = 2;
my $do_go = 0;
my $do_repeats = 0;

my %opt;

if ($ARGV[0]){
&Getopt::Std::getopts('v:j:u:p:s:H:h:e:P:gr', \%opt) || die "option problem" ;
}else{
&help_text; 
}


# get configuration from environment variables
&config; # this may fail but config can be on command line


my $enc_db = 'dk_funcgen_classify_55_1';# default, can be overridden by args
&process_arguments;

# hook up with the server
my $dbh = &make_contact($enc_db);
my $dbu = DBSQL::Utils->new($dbh);

#we need access to the latest ensembl core
#my $core_spec = 'ens_staging'; 
my $core_spec = 'current_core';
if($sp eq 'mus_musculus'){$core_spec = 'current_mouse_core'}
my $core_dbs = DBSQL::DBS->new($core_spec);



my $core_dbh = $core_dbs->connect();
my $core_dbu = DBSQL::Utils->new($core_dbh);
# we also want data from the corresponding GO database which is on same server
# so need the release number
my ($release) = $core_dbs->name =~ /.*_([0-9][0-9])_/;
my $go_db = 'ensembl_go_'.$release;

if($jump_label){
    goto $jump_label;
}
# we need a subset of the core tables so we dump them and read them in to
# the classify db
&commentary("clearing and creating dump directory $dump_dir") if $verbose;
&backtick("rm -rf $dump_dir");
&backtick("mkdir $dump_dir");

&make_sources_table($dbh,$core_dbs->name);


my $dump_templ = "mysqldump --opt --skip-lock-tables -h ".$core_dbs->host.
                                 " -u ".$core_dbs->user.
                                 " -P ".$core_dbs->port.
                                 " ".$core_dbs->name.
                                 ' %s '.
                                 " > $dump_dir".'%s'.".dump";

my $load_templ = "mysql -h ".$host.
                 " -u ".$user.
                 " -P ".$port.
                 " -p".$password.
                 " ".$enc_db.
                 " < $dump_dir".'%s'.".dump";

my @go_tables = ('term','term2term','graph_path');

my @core_tables = ('transcript',
                   'transcript_stable_id',
                   'gene',
                   'gene_attrib',
                   'attrib_type',
                   'exon',
                   'exon_transcript',
                   'coord_system',
                   'meta_coord',
                   'seq_region',
                   'seq_region_attrib',
                   'assembly',
#                   'ditag_feature',
#                   'ditag',
                   'seq_region_attrib',
                   'simple_feature',
                   'analysis',
                   'karyotype',
                   'xref',
                   'object_xref',
                   'external_db',
                   'translation',
                   'analysis_description'
		   );

# currently we always  need the repeats in order to filter the reg features
#if($do_repeats){
    push @core_tables, ('repeat_feature','repeat_consensus');
#}



foreach my $table (@core_tables){
                   
    unless( $core_dbu->table_exists($table)){die ("ERROR: The core database does not contain table $table")}

    my $command = sprintf($dump_templ,$table,$table);
    &backtick($command);
    $command = sprintf($load_templ,$table);
    #print $command."\n";
    &backtick($command);
}



$dump_templ = "mysqldump --opt --skip-lock-tables -h ".$core_dbs->host.
                                 " -u ".$core_dbs->user.
                                 " -P ".$core_dbs->port.
                                 " ".$go_db.
                                 ' %s '.
                                 " > $dump_dir".'%s'.".dump";

if($do_go){
    foreach my $table (@go_tables){

	unless( $core_dbu->table_exists($go_db.'.'.$table)){die ("ERROR: The $go_db database does not contain table $table")}

	my $command = sprintf($dump_templ,$table,$table);
	&backtick($command);
	$command = sprintf($load_templ,$table);
	#print $command."\n";
	&backtick($command);
    }
}



POST_IMPORT:

# create transcript-based feature tables
&transcript_features($dbh);


#POST_IMPORT:
&repeat_features($dbh,$dbu,$do_repeats);
#&repeat_features_old_version($dbh,$dbu,$do_repeats); # for contig coords


if($dbu->table_exists('de_ferrari_gene_classification')){
    &housekeeping_tissue_specific($dbh,$dbu);
}

#POST_IMPORT:
&gene_features($dbh);
&exon_features($dbh);

POST_EXON:

&intergenic_features($dbh,$dbu);
#&cage_ditag_transcript_tss($dbh,$dbu);
POST_DITAG:
&karyotype_features($dbh,$dbu);
POST_KARYOTYPE:

# simple features contains FirstEF CpG, Eponine and tRNAscan
# logic_name for CpG is 'CpG'

#POST_IMPORT:

&cpg_features($dbh,$dbu);
#&cpg_features_old_version($dbh,$dbu); # for contig coords

if($dbu->table_exists($slim_table) && $do_go){
    &go_term_features($dbh,$dbu,$slim_table);
}else{
    &commentary("Can't create go features cos $slim_table is not in database\n This is only needed for research purposes\n");
}



&clean_temp();
&backtick("rm -rf $dump_dir");
$dbh->disconnect;
$core_dbh->disconnect;
exit;


 
###################################################################
sub make_sources_table{
    my($dbh,$core_db) = @_;

    &commentary("creating sources table\n") if $verbose > 1;
    my @sql;
    push @sql,"drop table if exists sources";
    push @sql,"create table sources (db_type varchar(40),db_name varchar(40))";
    push @sql, "insert into sources values('core','$core_db')";
    &execute($dbh,@sql) or die;

}

sub in_list{

    return "('".join("','",@_)."')";

}


sub go_term_features{
    my($dbh,$dbu,$slim_table)=@_;

    &commentary("creating GO classified transcript features\n") if $verbose > 1;
    die "The GO tables have changed this method needs updating";

    my @sql;

    my $temp1 = &new_temp();
    push @sql,"drop table if exists $temp1";

    push @sql,"create table $temp1 select x.dbprimary_acc as acc,t.transcript_id,t.translation_id,tm.id as term_id,tm.name as go_term,tm.term_type from  object_xref ox, xref x ,external_db ed,translation t,term tm where ox.xref_id = x.xref_id and ox.ensembl_object_type = 'Translation' and x.external_db_id = ed.external_db_id and ed.db_name = 'GO' and ox.ensembl_id = t.translation_id and tm.acc = x.dbprimary_acc";
    push @sql,"alter table $temp1 add index(transcript_id)";
    push @sql,"alter table $temp1 add index(acc)";
    &execute($dbh,@sql);

    my $q = "select acc from $slim_table";
    my $slims_ref = $dbh->selectcol_arrayref($q);
    unless(defined  $slims_ref){
	die "query failed on :\n$q\n"->$dbh->errstr;
    }

    foreach my $slim (@$slims_ref){


	$q="SELECT DISTINCT descendant.acc FROM term, graph_path gp, term AS descendant where term.id=gp.term1_id and descendant.id=gp.term2_id and term.acc = '$slim'";
        my $aref = $dbh->selectcol_arrayref($q);
        unless(defined  $aref){
	    die "query failed on :\n$q\n"->$dbh->errstr;
        }
	my $in_list = &in_list(@$aref);

        $q = "select distinct name from term where acc = '$slim'";
	my $term = $dbh->selectrow_array($q);
        unless(defined  $term){
	    die "query failed on :\n$q\n"->$dbh->errstr;
        }

	$term =~ tr/ /_/;
	$term =~ s/\,//g;
        # this term is too long for a table name
	if($term eq 'nucleobase_nucleoside_nucleotide_and_nucleic_acid_metabolic_process'){
            $term = 'nucleo_etc_metabolic_process';
	}
        my $table = 'go_'.$term;
	print "$table\t$table\t$table\t0\t0\t0\n";

	@sql = ();
	push @sql,"drop table if exists $table";
        #push @sql,"create table $table select distinct sr.seq_region_id,sr.name as seq_region_name,t.transcript_id as feature_id,t.biotype as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr, $temp1 t1 where t.seq_region_id = sr.seq_region_id and t.transcript_id = t1.transcript_id and t1.acc in $in_list";
        # this gets the transcript plus its upstream 2500
        push @sql,"create table $table select distinct sr.seq_region_id,sr.name as seq_region_name,t.transcript_id as feature_id,t.biotype as feature_type, if(t.seq_region_strand = 1,t.seq_region_start-2500,t.seq_region_start) as feature_start, if(t.seq_region_strand = 1,t.seq_region_end,t.seq_region_end + 2500) as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr, $temp1 t1 where t.seq_region_id = sr.seq_region_id and t.transcript_id = t1.transcript_id and t1.acc in $in_list";
        &execute($dbh,@sql) or die;


    }

    @sql = ();
    push @sql,"drop table if exists $temp1";
    &execute($dbh,@sql) or die;

}


sub housekeeping_tissue_specific{
    my($dbh,$dbu) = @_;

    &commentary("creating housekeeping and tissue_specific gene features\n") if $verbose > 1;

    my @sql;

    my $temp1 = &new_temp();
    push @sql,"create table $temp1 select t.*,ts.stable_id as transcript_stable_id from protein_coding_transcript t, transcript_stable_id ts where t.feature_id = ts.transcript_id";
    push @sql,"drop table if exists protein_coding_transcript_hk";
    push @sql,"create table protein_coding_transcript_hk select distinct t1.* from $temp1 t1, de_ferrari_gene_classification d where t1.transcript_stable_id = d.transcript_stable_id and d.is_hk";
    push @sql, &col_types_and_indices('protein_coding_transcript_hk');
    push @sql,"update protein_coding_transcript_hk set feature_type = 'protein_coding_transcript_hk'"; 
    
    push @sql,"drop table if exists protein_coding_transcript_ts";
    push @sql,"create table protein_coding_transcript_ts select distinct t1.* from $temp1 t1, de_ferrari_gene_classification d where t1.transcript_stable_id = d.transcript_stable_id and d.is_hk =0"; 
    push @sql, &col_types_and_indices('protein_coding_transcript_ts');
push @sql,"update protein_coding_transcript_ts set feature_type = 'protein_coding_transcript_ts'";

    push @sql, "drop table $temp1";
    &execute($dbh,@sql) or die;

}

sub karyotype_features{
    my($dbh,$dbu)=@_;

    &commentary("creating karyotype features\n") if $verbose > 1;

    my $q = "select distinct stain from karyotype";
    my $stain_aref = $dbh->selectcol_arrayref($q) ;
    unless(defined $stain_aref && @$stain_aref > 0){
	die "got no data from query\n$q\n".$dbh->errstr;
    }

    my @sql;
    my @pos;
    foreach my $stain (@$stain_aref){
	my $table = "karyo_".$stain;
	if($stain =~ 'pos'){push @pos,$table}

	push @sql, "drop table if exists $table";
	push @sql, "create table $table select k.seq_region_id,sr.name as seq_region_name,karyotype_id as feature_id, '$table' as feature_type,seq_region_start as feature_start,seq_region_end as feature_end,0 as feature_strand from karyotype k, seq_region sr where k.seq_region_id = sr.seq_region_id and stain = '$stain'";
        push @sql,&col_types_and_indices($table);

    }
    push @sql, "drop table if exists karyo_gpos";
    push @sql, "create table karyo_gpos select * from karyo_gneg where 1=0";

    foreach my $table (@pos){
	push @sql, "insert into karyo_gpos select * from $table";
	push @sql, "drop table $table";

    }
    push @sql,&col_types_and_indices('karyo_gpos');

    &execute($dbh, @sql);



}


#for repeats on top level sequence
sub repeat_features{
    my($dbh,$dbu,$all_repeats)=@_;

    &commentary("creating repeat features\n") if $verbose > 1;


    my @sql;
    my $temp1 = &new_temp();
    # we lose a few mappings because some mappings are on parts of the contigs
    # which are not included in the assembly
    push @sql,"drop table if exists $temp1";
    
        my $q= "create table $temp1
        select f.*,
               f.seq_region_start as chr_start,
               f.seq_region_end as chr_end,
               f.seq_region_id as chromosome_id,
               asm_sr.name as chr_name 
from repeat_feature f,
     analysis an,
     seq_region asm_sr
where f.analysis_id = an.analysis_id
and an.logic_name = 'RepeatMask'
and f.seq_region_id = asm_sr.seq_region_id ";

	$q =~ tr/\n/ /;
#	print $q."\n";
	push @sql, $q;

	push @sql,"alter table $temp1 add index(repeat_consensus_id)";
        # now add the repeat type 
        push @sql,"drop table if exists all_repeats";
        push @sql,"create table all_repeats select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,repeat_feature_id as feature_id,concat(rc.repeat_class,'_repeat') as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1, repeat_consensus rc where rc.repeat_consensus_id = t1.repeat_consensus_id ";
	push @sql, "alter table all_repeats add index(feature_type)";
        #push @sql,&col_types_and_indices("all_repeats");
        push @sql,"drop table $temp1";

    
    &execute($dbh,@sql) or die;

    # now create separate tables for each type 

    $q = "select distinct feature_type from all_repeats";
    my $aref = $dbh->selectcol_arrayref($q);
    unless(defined $aref && @$aref > 0){
	die "failed to get any repeat types :\n$q\n".$dbh->errstr;
    }


    unless($all_repeats){$aref = ['Satellite_repeat','Satellite/centr_repeat'] }

    @sql = ();
    foreach my $type (@$aref){
        @sql = ();
	my $table = $type;
	$table =~ tr/\//_/;
	push @sql,"drop table if exists $table";
	push @sql,"create table $table select * from all_repeats where feature_type = '$type'";
	push @sql,"update $table set feature_type = '$table'";
        push @sql,&col_types_and_indices("$table");
        &execute($dbh,@sql) or die;
	print "$table\t$table\t$table\t0\t0\t0\n";
    }


}


# for repeats on contigs
sub repeat_features_old_version{
    my($dbh,$dbu)=@_;

    &commentary("creating repeat features\n") if $verbose > 1;
    # repeats should be mapped to contigs...this checks that they are
    my($coord_level,$coord_system_id)=&coord_system_check($dbh,'repeat_feature','contig');

    my @sql;
    my $temp1 = &new_temp();
    # we lose a few mappings because some mappings are on parts of the contigs
    # which are not included in the assembly
    push @sql,"drop table if exists $temp1";
    if($coord_level eq 'sequence_level'){
        my $q= "create table $temp1
        select f.*,
               if(a.ori=1,(a.asm_start+f.seq_region_start-1),(a.asm_start+a.cmp_end-f.seq_region_end)) as chr_start,
               if(a.ori=1,(a.asm_start+f.seq_region_end-1),(a.asm_start+a.cmp_end-f.seq_region_start)) as chr_end,
               a.asm_seq_region_id as chromosome_id,
               asm_sr.name as chr_name 
from repeat_feature f,
     analysis an,
     assembly a,
     seq_region asm_sr,
     seq_region_attrib sra , 
     attrib_type  at
where f.analysis_id = an.analysis_id
and an.logic_name = 'RepeatMask'
and f.seq_region_id = a.cmp_seq_region_id
and f.seq_region_start >= a.cmp_start
and f.seq_region_end <= a.cmp_end
and   a.asm_seq_region_id = asm_sr.seq_region_id
and sra.attrib_type_id = at.attrib_type_id 
and sra.seq_region_id = a.asm_seq_region_id
and at.code = 'toplevel' 
and sra.value = 1 ";

	$q =~ tr/\n/ /;
#	print $q."\n";
	push @sql, $q;

	push @sql,"alter table $temp1 add index(repeat_consensus_id)";
        # now add the repeat type 
        push @sql,"drop table if exists all_repeats";
        push @sql,"create table all_repeats select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,repeat_feature_id as feature_id,concat(rc.repeat_class,'_repeat') as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1, repeat_consensus rc where rc.repeat_consensus_id = t1.repeat_consensus_id ";
	push @sql, "alter table all_repeats add index(feature_type)";
        #push @sql,&col_types_and_indices("all_repeats");
        push @sql,"drop table $temp1";

    }
    &execute($dbh,@sql) or die;

    # now create separate tables for each type 

    my $q = "select distinct feature_type from all_repeats";
    my $aref = $dbh->selectcol_arrayref($q);
    unless(defined $aref && @$aref > 0){
	die "failed to get any repeat types :\n$q\n".$dbh->errstr;
    }

    @sql = ();
    foreach my $type (@$aref){
        @sql = ();
	my $table = $type;
	$table =~ tr/\//_/;
	$table =~ tr/-/_/;
	push @sql,"drop table if exists $table";
	push @sql,"create table $table select * from all_repeats where feature_type = '$type'";
	push @sql,"update $table set feature_type = '$table'";
        push @sql,&col_types_and_indices("$table");
        &execute($dbh,@sql) or die;
	print "$table\t$table\t$table\t0\t0\t0\n";
    }


}

sub cpg_features{
    my($dbh,$dbu)=@_;

    &commentary("creating CpG features\n") if $verbose > 1;

    my @sql;
    my $temp1 = &new_temp();
    # we lose a few mappings because some mappings are on parts of the contigs
    # which are not included in the assembly
    push @sql,"drop table if exists $temp1";
    
        my $q= "create table $temp1
        select f.*,
               f.seq_region_start as chr_start,
               f.seq_region_end as chr_end,
               f.seq_region_id as chromosome_id,
               asm_sr.name as chr_name 
from simple_feature f,
     analysis an,
     seq_region asm_sr 
where f.analysis_id = an.analysis_id
and an.logic_name = 'CpG'
and asm_sr.seq_region_id = f.seq_region_id  ";

	$q =~ tr/\n/ /;
#	print $q."\n";
	push @sql, $q;
        push @sql,"drop table if exists cpg_island";
        push @sql,"create table cpg_island select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,simple_feature_id as feature_id,'cpg_island' as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1 ";
        push @sql,&col_types_and_indices("cpg_island");
        push @sql,"drop table $temp1";

    
    &execute($dbh,@sql) or die;

    # sub classes of CpGs
    unless($dbu->table_exists('RNA_transcript')){
        warn "Can't generate CpG sub class... no RNA_transcript table present";
    }
    @sql = ();
    push @sql,"drop table if exists RNA_cpg";
    push @sql,"drop table if exists $temp1";
    # create temp table with RNA_transcript's  upstream  500bp
    push @sql,"create table $temp1 select seq_region_name,seq_region_id,if(feature_strand = 1,feature_start-500,feature_end) as feature_start, if(feature_strand = 1, feature_start,feature_end +500) as feature_end from RNA_transcript";
    
    push @sql,"alter table $temp1 add index(seq_region_name)";
    push @sql,"create table RNA_cpg select c.* from cpg_island c, $temp1  t where c.seq_region_name = t.seq_region_name and c.feature_end >= t.feature_start and c.feature_start <= t.feature_end ";
    push @sql,&col_types_and_indices("RNA_cpg");
    &execute($dbh,@sql);


    unless($dbu->table_exists('protein_coding_transcript')){
        warn "Can't generate CpG sub class... no protein_coding_transcript table present";
    }
    @sql = ();
    push @sql,"drop table if exists protein_coding_cpg";
    push @sql,"drop table if exists $temp1";
    # create temp table with protein_coding_transcript's upstream  500bp
    push @sql,"create table $temp1 select seq_region_name,seq_region_id,if(feature_strand = 1,feature_start-500,feature_end) as feature_start, if(feature_strand = 1, feature_start,feature_end +500) as feature_end from protein_coding_transcript";
    push @sql,"alter table $temp1 add index(seq_region_name)";
    push @sql,"create table protein_coding_cpg select c.* from cpg_island c, $temp1  t where c.seq_region_name = t.seq_region_name and c.feature_end >= t.feature_start and c.feature_start <= t.feature_end ";
    push @sql,&col_types_and_indices("protein_coding_cpg");
    &execute($dbh,@sql);

    @sql = ();
    push @sql,"drop table if exists non_promoter_cpg";
    push @sql,"create table non_promoter_cpg select * from cpg_island";
    push @sql, "alter table non_promoter_cpg add index(feature_id)";
    push @sql,"delete n from non_promoter_cpg n, RNA_cpg r where r.feature_id = n.feature_id";
    push @sql,"delete n from non_promoter_cpg n, protein_coding_cpg r where r.feature_id = n.feature_id";
    push @sql,&col_types_and_indices("non_promoter_cpg");
    &execute($dbh,@sql) or die;

}



sub cpg_features_old_version{
    my($dbh,$dbu)=@_;

    &commentary("creating CpG features\n") if $verbose > 1;

    # cpg islands should be mapped to contigs...this checks that they are
    my($coord_level,$coord_system_id)=&coord_system_check($dbh,'simple_feature','contig');

    my @sql;
    my $temp1 = &new_temp();
    # we lose a few mappings because some mappings are on parts of the contigs
    # which are not included in the assembly
    push @sql,"drop table if exists $temp1";
    if($coord_level eq 'sequence_level'){
        my $q= "create table $temp1
        select f.*,
               if(a.ori=1,(a.asm_start+f.seq_region_start-1),(a.asm_start+a.cmp_end-f.seq_region_end)) as chr_start,
               if(a.ori=1,(a.asm_start+f.seq_region_end-1),(a.asm_start+a.cmp_end-f.seq_region_start)) as chr_end,
               a.asm_seq_region_id as chromosome_id,
               asm_sr.name as chr_name 
from simple_feature f,
     analysis an,
     assembly a,
     seq_region asm_sr,
     seq_region_attrib sra , 
     attrib_type  at
where f.analysis_id = an.analysis_id
and an.logic_name = 'CpG'
and f.seq_region_id = a.cmp_seq_region_id
and f.seq_region_start >= a.cmp_start
and f.seq_region_end <= a.cmp_end
and   a.asm_seq_region_id = asm_sr.seq_region_id
and sra.attrib_type_id = at.attrib_type_id 
and sra.seq_region_id = a.asm_seq_region_id
and at.code = 'toplevel' 
and sra.value = 1 ";

	$q =~ tr/\n/ /;
#	print $q."\n";
	push @sql, $q;
        push @sql,"drop table if exists cpg_island";
        push @sql,"create table cpg_island select t1.chromosome_id as seq_region_id,t1.chr_name as seq_region_name,simple_feature_id as feature_id,'cpg_island' as feature_type,chr_start as feature_start,chr_end as feature_end , 0 as feature_strand from $temp1 t1 ";
        push @sql,&col_types_and_indices("cpg_island");
        push @sql,"drop table $temp1";

    }
    &execute($dbh,@sql) or die;

    # sub classes of CpGs
    unless($dbu->table_exists('RNA_transcript')){
        warn "Can't generate CpG sub class... no RNA_transcript table present";
    }
    @sql = ();
    push @sql,"drop table if exists RNA_cpg";
    push @sql,"drop table if exists $temp1";
    # create temp table with RNA_transcript's  upstream  500bp
    push @sql,"create table $temp1 select seq_region_name,seq_region_id,if(feature_strand = 1,feature_start-500,feature_end) as feature_start, if(feature_strand = 1, feature_start,feature_end +500) as feature_end from RNA_transcript";
    
    push @sql,"alter table $temp1 add index(seq_region_name)";
    push @sql,"create table RNA_cpg select distinct c.* from cpg_island c, $temp1  t where c.seq_region_name = t.seq_region_name and c.feature_end >= t.feature_start and c.feature_start <= t.feature_end ";
    push @sql,&col_types_and_indices("RNA_cpg");
    &execute($dbh,@sql);


    unless($dbu->table_exists('protein_coding_transcript')){
        warn "Can't generate CpG sub class... no protein_coding_transcript table present";
    }
    @sql = ();
    push @sql,"drop table if exists protein_coding_cpg";
    push @sql,"drop table if exists $temp1";
    # create temp table with protein_coding_transcript's upstream  500bp
    push @sql,"create table $temp1 select seq_region_name,seq_region_id,if(feature_strand = 1,feature_start-500,feature_end) as feature_start, if(feature_strand = 1, feature_start,feature_end +500) as feature_end from protein_coding_transcript";
    push @sql,"alter table $temp1 add index(seq_region_name)";
    push @sql,"create table protein_coding_cpg select distinct c.* from cpg_island c, $temp1  t where c.seq_region_name = t.seq_region_name and c.feature_end >= t.feature_start and c.feature_start <= t.feature_end ";
    push @sql,&col_types_and_indices("protein_coding_cpg");
    &execute($dbh,@sql);

    @sql = ();
    push @sql,"drop table if exists non_promoter_cpg";
    push @sql,"create table non_promoter_cpg select * from cpg_island";
    push @sql, "alter table non_promoter_cpg add index(feature_id)";
    push @sql,"delete n from non_promoter_cpg n, RNA_cpg r where r.feature_id = n.feature_id";
    push @sql,"delete n from non_promoter_cpg n, protein_coding_cpg r where r.feature_id = n.feature_id";
    push @sql,&col_types_and_indices("non_promoter_cpg");
    &execute($dbh,@sql) or die;

}


sub coord_system_check{
    my $dbh = shift;
    my $table = shift;
    my $expected = shift;

    my $q = "select mc.coord_system_id,cs.attrib,cs.name from meta_coord mc, coord_system cs where mc.table_name = 'simple_feature' and mc.coord_system_id = cs.coord_system_id";


    my $level;
    my @arr =  $dbh->selectrow_array($q) or die($dbh->errstr);
     
    my $problem = "$table coords not on expected coord_system ie $expected";
    if($arr[1] =~ 'default_version,sequence_level' &&
       $arr[2] =~ $expected ){
        $level = 'sequence_level';
        $problem = '';

    }

    if($problem){
	die("$problem");
    }

    my $coord_system_id= $arr[0];
    return ($level,$coord_system_id)
}



sub exon_features{
    my($dbh)=@_;
    &commentary("creating exon features\n") if $verbose > 1;
    my @sql;
    my $temp1 = &new_temp();

    # denormalise exon, transcript and exon_transcript
    push @sql, "drop table if exists $temp1";
    push @sql, "create table $temp1 select e.*,et.rank,et.transcript_id,t.biotype from exon e, exon_transcript et, transcript t where e.exon_id = et.exon_id and et.transcript_id = t.transcript_id order by et.transcript_id,e.seq_region_start";
    push @sql, "alter table $temp1 add index(transcript_id)";

    # get rank of last exon for each transcript
    my $temp2 = &new_temp();
    push @sql, "drop table if exists $temp2";
    push @sql, "create table $temp2 select transcript_id,max(rank) as max_rank from $temp1 group by transcript_id";
    push @sql, "alter table $temp2 add index(transcript_id)";
    my $temp3 = &new_temp();
    push @sql, "drop table if exists $temp3";
    push @sql, "create table $temp3 select t1.*,t2.max_rank from $temp1 t1, $temp2 t2 where t1.transcript_id = t2.transcript_id";
    push @sql, "drop table if exists $temp2";
    push @sql, "drop table if exists $temp1";
    push @sql, "alter table $temp3 add index(transcript_id)";
    push @sql, "alter table $temp3 add index(max_rank)";
    push @sql, "alter table $temp3 add index(rank)";
    push @sql, "alter table $temp3 add index(seq_region_id)";
    

    # get first intron for each transcript
    my $temp4 = &new_temp();
    push @sql, "create table $temp4 select t3.seq_region_id,sr.name as seq_region_name,t3.transcript_id as feature_id, concat(t3.biotype,'_intron1') as feature_type,if(t3.seq_region_strand = 1,t3.seq_region_end+1,b.seq_region_end+1) as feature_start, if(t3.seq_region_strand = 1,b.seq_region_start-1,t3.seq_region_start -1) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, $temp3 b, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.transcript_id = b.transcript_id and t3.rank = 1 and b.rank = 2"; 
    push @sql,"drop table if exists intron1";
    push @sql,"alter table $temp4 rename as intron1";
    push @sql,&col_types_and_indices("intron1");
    push @sql,&split_by_biotype("intron1");
    &execute($dbh,@sql);



    # promoter defined as 500 bp upstream
    push @sql, "drop table if exists exon1_plus_promoter";
    push @sql, "create table exon1_plus_promoter select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_exon1_plus_promoter') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.rank = 1";
    push @sql,&col_types_and_indices("exon1_plus_promoter");
    push @sql,&split_by_biotype("exon1_plus_promoter");

    
    # enhancer defined as 2500 bp upstream
    push @sql, "drop table if exists exon1_plus_enhancer";
    push @sql, "create table exon1_plus_enhancer select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_exon1_plus_enhancer') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-2500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +2500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.rank = 1";
    push @sql,&col_types_and_indices("exon1_plus_enhancer");
    push @sql,&split_by_biotype("exon1_plus_enhancer");


 
    push @sql, "drop table if exists single_exon_gene";
    push @sql, "create table single_exon_gene select t3.seq_region_id,sr.name as seq_region_name,exon_id as feature_id, concat(t3.biotype, '_single_exon_gene') as feature_type,t3.seq_region_start as feature_start, t3.seq_region_end as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank = 1";
    push @sql,&col_types_and_indices("single_exon_gene");    
    push @sql,&split_by_biotype("single_exon_gene");

    push @sql, "drop table if exists single_exon_gene_plus_enhancer";
    push @sql, "create table single_exon_gene_plus_enhancer select  t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, concat(t3.biotype,'_single_exon_gene_plus_enhancer') as feature_type,if(seq_region_strand = 1,t3.seq_region_start-2500,seq_region_start) as feature_start, if(seq_region_strand = 1,seq_region_end,seq_region_end +2500) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank = 1";
    push @sql,&col_types_and_indices("single_exon_gene_plus_enhancer");    
    push @sql,&split_by_biotype("single_exon_gene_plus_enhancer");

    # gene body is all exons and introns except the first of each
    push @sql, "drop table if exists gene_body";
    push @sql, "create table gene_body select t3.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,  concat(t3.biotype,'_gene_body') as feature_type,min(t3.seq_region_start) as feature_start, max(t3.seq_region_end) as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id and t3.max_rank != 1 and t3.rank != 1 group by transcript_id";
    push @sql,&col_types_and_indices("gene_body");  
    push @sql,&split_by_biotype("gene_body");


    push @sql, "drop table if exists exon_plus_flanks_500";
    push @sql, "create table exon_plus_flanks_500 select distinct t3.seq_region_id,sr.name as seq_region_name,exon_id as feature_id, concat(t3.biotype,'_exon_plus_flanks_500') as feature_type,t3.seq_region_start-500 as feature_start, seq_region_end + 500 as feature_end,t3.seq_region_strand as feature_strand from $temp3 t3, seq_region sr where t3.seq_region_id = sr.seq_region_id ";
    push @sql,&col_types_and_indices("exon_plus_flanks_500");    
    push @sql,&split_by_biotype("exon_plus_flanks_500");


#    push @sql, "drop table if exists $temp3";
    &execute($dbh,@sql) or die;
    
}

sub split_by_biotype{
    my($feat_name) = @_;

    my @sql;

    my $prots = "'protein_coding_$feat_name', 'IG_V_gene_$feat_name', 'IG_C_gene_$feat_name', 'IG_J_gene_$feat_name', 'IG_D_gene_$feat_name'";
    push @sql,"drop table if exists protein_coding_$feat_name";
    push @sql,"create table protein_coding_$feat_name select * from $feat_name where feature_type in ($prots)";
    push @sql,"update protein_coding_$feat_name set feature_type = 'protein_coding_$feat_name'"; 
    push @sql,&col_types_and_indices("protein_coding_$feat_name");

    push @sql,"drop table if exists RNA_gene_$feat_name";
    push @sql,"create table RNA_gene_$feat_name select * from $feat_name where feature_type like '%RNA_".$feat_name."'";
    #************* lump all RNA sub types together *********************
    #push @sql,"update RNA_gene_$feat_name set feature_type = 'RNA_gene_$feat_name'"; 
    push @sql,&col_types_and_indices("RNA_gene_$feat_name");

    if($feat_name eq 'exon1_plus_enhancer'){
	push @sql,"drop table if exists snRNA_gene_$feat_name";
	push @sql,"create table snRNA_gene_$feat_name select * from $feat_name where feature_type like '%snRNA_".$feat_name."'";
	push @sql,"update snRNA_gene_$feat_name set feature_type = 'snRNA_gene_$feat_name'";
	push @sql,&col_types_and_indices("snRNA_gene_$feat_name");

	push @sql,"drop table if exists snoRNA_gene_$feat_name";
	push @sql,"create table snoRNA_gene_$feat_name select * from $feat_name where feature_type like '%snoRNA_".$feat_name."'";
	push @sql,"update snoRNA_gene_$feat_name set feature_type = 'snoRNA_gene_$feat_name'";
	push @sql,&col_types_and_indices("snoRNA_gene_$feat_name");

	push @sql,"drop table if exists miRNA_gene_$feat_name";
	push @sql,"create table miRNA_gene_$feat_name select * from $feat_name where feature_type like '%miRNA_".$feat_name."'";
	push @sql,"update miRNA_gene_$feat_name set feature_type = 'miRNA_gene_$feat_name'";
	push @sql,&col_types_and_indices("miRNA_gene_$feat_name");

	push @sql,"drop table if exists miscRNA_gene_$feat_name";
	push @sql,"create table miscRNA_gene_$feat_name select * from $feat_name where feature_type like '%miscRNA_".$feat_name."'";
	push @sql,"update miscRNA_gene_$feat_name set feature_type = 'miscRNA_gene_$feat_name'";
	push @sql,&col_types_and_indices("miscRNA_gene_$feat_name");

	push @sql,"drop table if exists rRNA_gene_$feat_name";
	push @sql,"create table rRNA_gene_$feat_name select * from $feat_name where feature_type like '%rRNA_".$feat_name."'";
	push @sql,"update rRNA_gene_$feat_name set feature_type = 'rRNA_gene_$feat_name'";
	push @sql,&col_types_and_indices("rRNA_gene_$feat_name");
    }

    my $pseuds = "'pseudogene_$feat_name','repeat_$feat_name','retrotransposed_$feat_name'";
    push @sql,"drop table if exists pseudogene_$feat_name";
    push @sql,"create table pseudogene_$feat_name select * from $feat_name where feature_type in ($pseuds)";
    push @sql,&col_types_and_indices("pseudogene_$feat_name");




    return @sql;
}




sub cage_ditag_transcript_tss{
    my($dbh,$dbu)=@_;

    &commentary("creating cage_ditag regions \n") if $verbose > 1;

    my @sql;
    my $temp1 = &new_temp();
    push @sql, "drop table if exists $temp1";
    push @sql, "create table $temp1 select f.seq_region_id,sr.name as seq_region_name,if(f.seq_region_strand = 1,f.seq_region_start,f.seq_region_end) as feature_start,f.seq_region_strand as feature_strand from transcript f , seq_region sr where sr.seq_region_id = f.seq_region_id and f.biotype not like '%pseudogene%' and f.biotype not in('repeat','retrotransposed')";

    push @sql, "insert into $temp1  select f.seq_region_id,sr.name as seq_region_name,if(f.seq_region_strand = 1,f.seq_region_start,f.seq_region_end) as feature_start,f.seq_region_strand from ditag_feature f , seq_region sr where sr.seq_region_id = f.seq_region_id";

    my $temp2 = &new_temp();
    push @sql, "drop table if exists $temp2";
    push @sql, "create table $temp2 select distinct * from $temp1";

    # we may want to do clustering here to reduce the total number of features



    push @sql, "drop table if exists tss_upstream_500";
    push @sql, "create table tss_upstream_500 select seq_region_id,seq_region_name,if(feature_strand = -1,feature_start,feature_start-500) as feature_start,if(feature_strand = -1,feature_start+500,feature_start) as feature_end,feature_strand,'tss_upstream_500' as feature_type from $temp2";
push @sql,"alter table tss_upstream_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,&col_types_and_indices("tss_upstream_500");


    push @sql, "drop table if exists tss_downstream_500";
    push @sql, "create table tss_downstream_500 select seq_region_id,seq_region_name,if(feature_strand = -1,feature_start-500,feature_start) as feature_start,if(feature_strand = -1,feature_start,feature_start+500) as feature_end,feature_strand,'tss_downstream_500' as feature_type from $temp2";
push @sql,"alter table tss_downstream_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,&col_types_and_indices("tss_downstream_500");


    push @sql, "drop table if exists tss_centred_500";
    push @sql, "create table tss_centred_500 select seq_region_id,seq_region_name,feature_start-250 as feature_start,feature_start+250 as feature_end,feature_strand,'tss_centred_500' as feature_type from $temp2";
push @sql,"alter table tss_centred_500 add column feature_id int(10) not null auto_increment primary key";
    push @sql,&col_types_and_indices("tss_centred_500");


    # 
    &execute($dbh,@sql) or die;

    @sql = ();
    push @sql,"drop table if exists $temp1";
    push @sql,"drop table if exists $temp2";
    &execute($dbh,@sql) or die;

}




sub intergenic_features{
    my ($dbh,$dbu) = @_;
    &commentary("creating intergenic regions\n") if $verbose > 1;
    my @sql;

    # we process one seq_region ata a time
    # if we want to include more biotypes of gene then we should select from 
    # the gene table
    my $q = "select distinct seq_region_id from protein_coding_gene";
    my $regions_aref = $dbh->selectcol_arrayref($q);
    unless(@$regions_aref > 0){ die "no data returned by :\n$q\n" }


    @sql = ();
    my $temp1 = &new_temp();
    push @sql, "drop table if exists $temp1";
    # start and end are signed so we can use -ve numbers as flags
    push @sql, "create table $temp1 (seq_region_id int(10) unsigned,feature_start int(10), feature_end int(10))";
    &execute($dbh,@sql) or die;

    # first we get the regions which lie around/between the genes
    my @intergenics;
    foreach my $id (@$regions_aref){

	@intergenics = ();

	my $region_length = $dbu->get_count("select length from seq_region where seq_region_id = $id");

	$q = "select feature_start,feature_end from protein_coding_gene where seq_region_id = $id order by feature_start";
	my $genes_aaref = $dbh->selectall_arrayref($q);
        unless(defined $genes_aaref){
            die "query failed on:\n$q\n".$dbh->errstr;
	}

	my $last_end;
	#print "gene 0 ".join("\t",@{$genes_aaref->[0]})."\n";
        if($genes_aaref->[0]->[0] == 1){
	    $last_end = $genes_aaref->[0]->[1];
	}else{
	    my @intergenic;
            $intergenic[0] = -1; # use -ve to flag chromosome end
            $intergenic[1] = $genes_aaref->[0]->[0] -1;
	    push @intergenics,\@intergenic;
            $last_end = $genes_aaref->[0]->[1];
	}

        
        for(my $i=1;$i< @$genes_aaref;$i++){

	    #print "gene $i ".join("\t",@{$genes_aaref->[$i]})."\n";
	    my @intergenic;
	    if($genes_aaref->[$i]->[0] > $last_end+1){
                $intergenic[0] = $last_end+1;
                $intergenic[1] = $genes_aaref->[$i]->[0] -1;
 	        push @intergenics,\@intergenic;
	    }

            if($genes_aaref->[$i]->[1] > $last_end){
                $last_end = $genes_aaref->[$i]->[1];
	    }

	}

        if($last_end < $region_length){
	    my @intergenic;
            $intergenic[0] = $last_end+1;
            $intergenic[1] = - $region_length; # use -ve to flag chromosome end
 	    push @intergenics,\@intergenic;
	}

	@sql = ();
        foreach my $aref (@intergenics){
	    my $q = "insert into $temp1 values($id,";
	    $q .= join(',',@$aref);
	    $q .= ")";
	    push @sql,$q;
	    #print "inter ".join("\t",@$aref)."\n" ;
	}
	&execute($dbh,@sql) or die;
	
    }

    # now we create gene-distal regions by shortening the regions in $temp1
    &intergenic_variant($dbh,$temp1,2500);
    &intergenic_variant($dbh,$temp1,5000);
    &intergenic_variant($dbh,$temp1,10000);

    &execute($dbh,"drop table $temp1") or die $dbh->errstr;

}

sub intergenic_variant{
    my($dbh,$temp1,$dist) = @_;

    &commentary("creating intergenic regions $dist\n") if $verbose > 1;

    my $final_table = 'intergenic_'.$dist;
    my @sql;
    push @sql,"drop table if exists $final_table";
    my $temp2 = &new_temp();
    push @sql,"drop table if exists $temp2";
    push @sql,"create table $temp2 select seq_region_id,feature_start,feature_end from $temp1 where 1=0";
    # only add/subtract $dist if we are not at the chromosome end
    push @sql,"insert into $temp2 select seq_region_id,if(feature_start> 0,feature_start +$dist,feature_start * -1) as feature_start,if(feature_end > 0,feature_end - $dist, feature_end * -1) as feature_end from $temp1";
    
    push @sql,"delete from $temp2 where feature_end <= feature_start";
    push @sql,"create table $final_table select t2.seq_region_id,sr.name as seq_region_name,'$final_table' as feature_type,t2.feature_start,t2.feature_end,'0' as feature_strand from $temp2 t2, seq_region sr where t2.seq_region_id = sr.seq_region_id";
    push @sql,"alter table $final_table add column feature_id int(10) not null auto_increment primary key";

    push @sql,&col_types_and_indices($final_table);

    &execute($dbh,@sql) or die;
    &execute($dbh,"drop table $temp2") or die;

}

sub col_types_and_indices{
    my $table = shift;

    my @sql;

    push @sql,"alter table $table modify column feature_start int(10) unsigned";
    push @sql,"alter table $table modify column feature_end int(10) unsigned";
    push @sql,"alter table $table modify column feature_strand tinyint(2)";
    push @sql,"alter table $table modify column feature_type varchar(45)";
    push @sql,"alter table $table add index(seq_region_id)";
    push @sql,"alter table $table add index(seq_region_name)";
    push @sql,"alter table $table add index(feature_strand)";
    push @sql,"alter table $table add index(feature_start)";
    push @sql,"alter table $table add index(feature_end)";
    return @sql;
}

sub gene_features{
    my ($dbh) = @_;
# not all transcripts are on chromosomes but if we only want those that are
#select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, from transcript t, seq_region sr, coord_system cs where t.seq_region_id = sr.seq_region_id and cs.coord_system_id = sr.coord_system_id and cs.name = 'chromosome' ;

#if we want UTRs the data needed is in the translation and exon tables
# select transcript_id,if(e1.seq_region_strand = 1,e1.seq_region_start,e2.seq_region_start+seq_end) as 5utr_start,if(e1.seq_region_strand = 1,e1.seq_region_start+tl.seq_start-2,e2.seq_region_end) as 5utr_end from translation tl,exon e1,exon e2 where e1.exon_id = tl.start_exon_id and e2.exon_id = tl.end_exon_id ; 


    my @sql;
    push @sql,"drop table if exists protein_coding_gene";
    push @sql,"create table protein_coding_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,'protein_coding_gene' as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_gene add index(seq_region_name)";
    push @sql,"alter table protein_coding_gene add index(seq_region_id)";


    push @sql,"drop table if exists pseudogene";
    push @sql,"create table pseudogene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('pseudogene','repeat','retrotransposed')";
    push @sql,"alter table pseudogene add index(seq_region_name)";
    push @sql,"alter table pseudogene add index(seq_region_id)";


    push @sql,"drop table if exists RNA_gene";
    push @sql,"create table RNA_gene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype not like '%pseudogene%'";
    push @sql,"alter table RNA_gene add index(seq_region_name)";
    push @sql,"alter table RNA_gene add index(seq_region_id)";


    push @sql,"drop table if exists RNA_pseudogene";
    push @sql,"create table RNA_pseudogene select sr.seq_region_id,sr.name as seq_region_name,gene_id as feature_id,concat(biotype,'_gene') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from gene t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype like '%pseudogene%'";
    push @sql,"alter table RNA_pseudogene add index(seq_region_name)";
    push @sql,"alter table RNA_pseudogene add index(seq_region_id)";



    &execute($dbh,@sql) or die;


}

sub transcript_features{
    my ($dbh) = @_;
# not all transcripts are on chromosomes but if we only want those that are
#select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id, from transcript t, seq_region sr, coord_system cs where t.seq_region_id = sr.seq_region_id and cs.coord_system_id = sr.coord_system_id and cs.name = 'chromosome' ;

    my @sql;
    push @sql,"drop table if exists protein_coding_transcript";
    push @sql,"create table protein_coding_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,'protein_coding_transcript' as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_transcript add index(seq_region_name)";
    push @sql,"alter table protein_coding_transcript add index(seq_region_id)";



    push @sql,"drop table if exists protein_coding_transcript_downstream_2500";
    push @sql,"create table protein_coding_transcript_downstream_2500 select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,'protein_coding_transcript_downstream_2500' as feature_type,if(seq_region_strand = 1, t.seq_region_end+1,t.seq_region_start -2500) as feature_start,if(seq_region_strand = 1,t.seq_region_end+2500,t.seq_region_start - 1) as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('protein_coding')";
    push @sql,"alter table protein_coding_transcript_downstream_2500 add index(seq_region_name)";
    push @sql,"alter table protein_coding_transcript_downstream_2500 add index(seq_region_id)";





    push @sql,"drop table if exists pseudogene_transcript";
    push @sql,"create table pseudogene_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype in ('pseudogene','repeat','retrotransposed')";
    push @sql,"alter table pseudogene_transcript add index(seq_region_name)";
    push @sql,"alter table pseudogene_transcript add index(seq_region_id)";


    push @sql,"drop table if exists RNA_transcript";
    push @sql,"create table RNA_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype not like '%pseudogene%'";
    push @sql,"alter table RNA_transcript add index(seq_region_name)";
    push @sql,"alter table RNA_transcript add index(seq_region_id)";


    push @sql,"drop table if exists RNA_pseudogene_transcript";
    push @sql,"create table RNA_pseudogene_transcript select sr.seq_region_id,sr.name as seq_region_name,transcript_id as feature_id,concat(biotype,'_transcript') as feature_type, t.seq_region_start as feature_start,t.seq_region_end as feature_end, t.seq_region_strand as feature_strand from transcript t, seq_region sr where t.seq_region_id = sr.seq_region_id and t.biotype like '%RNA%' and t.biotype like '%pseudogene%'";
    push @sql,"alter table RNA_pseudogene_transcript add index(seq_region_name)";
    push @sql,"alter table RNA_pseudogene_transcript add index(seq_region_id)";



    &execute($dbh,@sql) or die;
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

    warn "executing $command \n" if ($verbose > 1);

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

    &commentary("processing SQL.");

    foreach my $query(@array){
	
    	&commentary(".");

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
        &err("failed to connect to database $db using\n$dsn $user $password");
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

    if (exists $opt{v}){
        $verbose = $opt{v}; 
    }

    $sp = 'homo_sapiens';
    if (exists $opt{s}){
        $sp = lc($opt{s});
    }


    if (exists $opt{r}){
        $do_repeats = 1;
    }

    if (exists $opt{g}){
	$do_go = 1;
    }


    if (exists $opt{j}){
	$jump_label = $opt{j} ;
    }


    if  (exists $opt{e}){
        $enc_db = $opt{e};
    } 

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
  }

    print STDERR <<"END_OF_TEXT";

    gen_feats_4_classification.pl [-h] for help
                  [-e] <db_name> use this database for the end product
                  [-g] flag - create the GO classified gene tables
                  [-r] flag - create the repeat_feature tables
                  [-j] <label> jump to label then start execution
                       POST_IMPORT
                  [-H] <host machine> eg ens-genomics2
                  [-u] <database user> 
                  [-p] <mysql password> 
                  [-P] <mysql port> 
                  [-s] <species> eg -smus_musculus, default = homo_sapiens 
                  [-v] <integer> verbosity level 0,1 or 2; default=0 
                  [-] 
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
