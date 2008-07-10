#!/software/bin/perl 

=head1 DESCRIPTION



=head1 AUTHOR(S)

dkeefe@ebi.ac.uk

=head1 USAGE


=head1 EXAMPLES

   not_masked_not_probed_in_trfdust.pl -e dk_5prime_features

=head1 SEE ALSO


=head1 CVS

 $Log: not supported by cvs2svn $



=cut


use strict;
use DBI;
use Env;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Std;
use IO::Handle;
use IO::File;
use lib '/nfs/acari/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
use DBSQL::Utils;
use DBSQL::DBS; # to use databases listed in ~/dbs.
use Encode::TBAExtractor;

use constant  NO_ROWS => '0E0';
$| = 1; #no output buffer

my($user, $password, $driver, $host, $port);
my $outfile='';
my $infile='';
my $id_list;
my $sp;
my $verbose = 0;



my %opt;

if ($ARGV[0]){
&Getopt::Std::getopts('u:p:s:e:H:h:e:n:N:o:i:P:', \%opt) || die ;
}else{
&help_text; 
}


# get configuration from environment variables
&config; # this may fail but config can be on command line

my $enc_db = 'dk_barsky_acetylation';# default, can be overridden by args
&process_arguments;

# hook up with the server
my $dbh = &make_contact($enc_db);
my $dbu = DBSQL::Utils->new($dbh);


#we need access to a seq_region table
#my $msa_dbs = DBSQL::DBS->new('human_core_36');
#my $msa_dbh = $msa_dbs->connect();
#my $msa_dbu = DBSQL::Utils->new($msa_dbh);
#unless( $msa_dbu->table_exists('seq_region')){die ("ERROR: seq_region table is missing from core db ".$msa_dbs->name)}


my @files = ('H2AK5ac','H2AK9ac','H2BK5ac','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K9ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K36ac','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K91ac');

my $stem = 'CD4-%s.bed';

#get the mappings into the DB
foreach my $mod (@files){
#$mod =  'H2AK9ac';   

    $infile = sprintf($stem,$mod);

    my $table = 'cd4_'.lc($mod);
    print "creating table $table from file $infile\n";
    my @sql;
    push @sql, "drop table if exists $table";
    push @sql, "create table $table (chromosome varchar(40) , feature_start int(11), feature_end int(11),map_type varchar(40),score int(10), strand char(1))";
    push @sql, "load data local infile '$infile' into table $table";

    push @sql, "update $table set feature_start = feature_start+1";
    push @sql, "alter table $table add column chr_name varchar(40)";
    push @sql, "update  $table set  chr_name = substring_index(chromosome,'chr',-1)";
    push @sql, "alter table $table add index(chr_name)";

    &execute($dbh,@sql) or die;
    
#exit;    
}


$dbh->disconnect;
exit;


 
###################################################################



sub create_regions{
    my($dbh,$dbu,$chr,$len,$source_table,$table) = @_;

 
 
    my @sql;
    my @g;
    print "initialising\n";
    for(my $i=1;$i<=$len;$i++){
        $g[$i] = 0;
    }

    # add the TRF regions - table already contains trimmed encode coords
    print "adding\n";
    &add_score_to_array($dbh,$chr,$len,\@g,$source_table);

    # add the DUSTed regions - table already contains trimmed encode coords
    #&populate_array($dbh,$chr,$len,\@g,'dusted_regions');


=head1
    # remove repeats
    @sql;
    push @sql,"drop table if exists temp";
    # rmsk already has encode coords just trim ends to be within encode regions
    push @sql,"create table temp select t.encode_region_name,
                if(t.feature_enc_start < 1,1,t.feature_enc_start ) as feature_enc_start,
                if(t.feature_enc_end > e.chrom_end,e.chrom_end,t.feature_enc_end ) as feature_enc_end from rmsk t, homo_sapiens_encode_regions e where t.encode_region_name = e.encode_region_name and e.encode_region_name = '$chr' ";
    &execute($dbh,@sql) or die;


    &de_populate_array($dbh,$chr,$len,\@g,'temp');

    # remove probes
    @sql = ();
    push @sql,"drop table if exists temp";
    push @sql,"create table temp select t.*, e.encode_region_name,
                if(t.feature_chrom_start > e.chrom_start, t.feature_chrom_start-e.chrom_start +1,1) as feature_enc_start,
                if(t.feature_chrom_end < e.chrom_end,t.feature_chrom_end-e.chrom_start+1,e.chrom_end-e.chrom_start+1 ) as feature_enc_end from affy_probes t, homo_sapiens_encode_regions e where t.feature_chrom_start <= e.chrom_end and t.feature_chrom_end >= e.chrom_start and t.chromosome = e.chromosome and e.encode_region_name = '$chr' ";
    &execute($dbh,@sql) or die;


    &de_populate_array($dbh,$chr,$len,\@g,'temp');

=cut

    print "reading\n";
    my $thresh=5;
    my $inside = 0;
    my $start;
    my $end;
    
    @sql = ();
    my $tot = 0;
    for(my $i=1;$i<=$len;$i++){

        if($g[$i] == 1){ $tot++ }

        # we are starting a new region
        if( $inside == 0 && $g[$i] >= $thresh){

	    $start = $i;
	    $inside = 1;
	    

	}

	if( $inside == 1 && $g[$i] < $thresh){
        # we have just left a region
	    $inside = 0;
	    $end = $i-1;
            push @sql,"insert into $table values('$chr','$start','$end')";
        }
    }
    # if we are still in a region at the end record it
    if( $inside == 1){
        push @sql,"insert into $table values('$chr','$start','$len')";
    }

    print "inserting\n";
    &execute($dbh,@sql) or die;

    return $tot;
}


sub add_score_to_array{
    my($dbh,$chr_name,$chr_len,$aref,$table)=@_;

    #print "\n $chr_name";

    my %weight = (333 => 1,
                  800 => 3,
                  1000 => 4);


    my $q = "select feature_start,feature_end,score from $table where chr_name = '$chr_name' order by feature_start,feature_end ";
    my $coords_ref =$dbh->selectall_arrayref($q);

    foreach my $ref (@$coords_ref){
	my $start = $ref->[0];
	my $end =  $ref->[1];
	my $score =  $ref->[2];

	#print $table." $start $end\n";
        for(my $i=$start;$i<=$end;$i++){
    	    $aref->[$i] += $weight{$score};
	    if($end > $chr_len){print "end greater than chr_len\n"}
        }
    }
}



sub add_to_array{
    my($dbh,$chr_name,$chr_len,$aref,$table)=@_;

    #print "\n $chr_name";


    my $q = "select feature_start,feature_end from $table where chr_name = '$chr_name' order by feature_start,feature_end ";
    my $coords_ref =$dbh->selectall_arrayref($q);

    foreach my $ref (@$coords_ref){
	my $start = $ref->[0];
	my $end =  $ref->[1];

	#print $table." $start $end\n";
        for(my $i=$start;$i<=$end;$i++){
    	    $aref->[$i] += 1;
	    if($end > $chr_len){print "end greater than chr_len\n"}
        }
    }
}




sub populate_array{
    my($dbh,$chr_name,$chr_len,$aref,$table)=@_;

    #print "\n $chr_name";


    my $q = "select feature_enc_start,feature_enc_end from $table where encode_region_name = '$chr_name' order by feature_enc_start,feature_enc_end ";
    my $coords_ref =$dbh->selectall_arrayref($q);

    foreach my $ref (@$coords_ref){
	my $start = $ref->[0];
	my $end =  $ref->[1];

	#print $table." $start $end\n";
        for(my $i=$start;$i<=$end;$i++){
    	    $aref->[$i] = 1;
	    if($end > $chr_len){print "end greater than enc_len\n"}
        }
    }
}

sub and_with_array{
    my($dbh,$chr_name,$chr_len,$aref,$table)=@_;

    #print "\n $chr_name";


    my $q = "select feature_enc_start,feature_enc_end from $table where encode_region_name = '$chr_name' order by feature_enc_start,feature_enc_end ";
    my $coords_ref =$dbh->selectall_arrayref($q);

    foreach my $ref (@$coords_ref){
	my $start = $ref->[0];
	my $end =  $ref->[1];

	#print $table." $start $end\n";
        for(my $i=$start;$i<=$end;$i++){
    	    $aref->[$i] = ($aref->[$i] && 1);
	    if($end > $chr_len){die "end greater than enc_len\n"}
        }
    }
}


sub de_populate_array{
    my($dbh,$chr_name,$chr_len,$aref,$table)=@_;

    #print "\n $chr_name";


    my $q = "select feature_enc_start,feature_enc_end from $table where encode_region_name = '$chr_name' order by feature_enc_start,feature_enc_end ";
    my $coords_ref =$dbh->selectall_arrayref($q);

    foreach my $ref (@$coords_ref){
	my $start = $ref->[0];
	my $end =  $ref->[1];

	#print $table." $start $end\n";
        for(my $i=$start;$i<=$end;$i++){
    	    $aref->[$i] = 0;
	    if($end > $chr_len){print "end greater than enc_len\n"}
        }
    }
}


sub perc_gc{
    my($msa,$enc,$len) = @_;

    my $ref_slice = uc($msa->reference_slice($enc,1,$len,'degap'));
   
    unless(length($ref_slice) == $len){die "sequence wrong length ".length($ref_slice)." should be $len" }

    my $g = $ref_slice =~ tr/G/G/;
    my $c =  $ref_slice =~ tr/C/C/;
my $n =  $ref_slice =~ tr/N/N/;

    print "$enc $g,$c,$n\n";

    return( sprintf("%2.2f",(100*($g+$c))/$len));


}





sub backtick{
    my $command = shift;

    warn "executing $command \n" if ($verbose ==2);

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
    print  "$_[0]";
}



sub get_names_from_file{
    my $file=shift;

    open(IN, "< $file") or die "couldn't open list file $file";

    my @ret;
    while( <IN> ){
        chop;
	push  @ret, $_ ;
    }

    my $text_list = join(",",@ret);
    #print $text_list;
    return $text_list;

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
  



sub print_mart_results{
    my $res_aref = shift;
    my $ma; #dummy to avoid undef err

    if( $ma->result_type() eq 'sequence'){
	foreach my $aref (@$res_aref){
	    print $aref->[0]."\n";
	    print $aref->[1]."\n";
	}
    }else{
	foreach my $aref (@$res_aref){
	    foreach my $field (@$aref){
	           ($field)? print $field." ":print "unavailable ";
	    }
	    print "\n";
	}
    }
}


sub make_contact{
    my $db = shift;

    unless($driver && $db && $host && $port && $user){
	&err("DB connection parameters not set");
        exit(1);
    }

    # hook up with the server
    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port;mysql_local_infile=1";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 0});
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


    if (exists $opt{o}){
        $outfile = $opt{o};
    }  


    if (exists $opt{i}){
        $infile = $opt{i};
    }  


    if (exists $opt{n}){
        $id_list = $opt{n};
    }

    if (exists $opt{N}){
	$infile = $opt{N};
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

    non_rmsk_regions.pl [-h] for help
                  [-e] use a particular local encode catalog
                  [-H] <host machine> eg ecs2
                  [-u] <database user> 
                  [-o] <output file> - name of a file for output
                  [-p] <mysql password> 
                  [-P] <mysql port> 
            


END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
