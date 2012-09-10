#!usr/bin/env perl
use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Test::More;

#change database parameters if needed

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => 'ensadmin',
    -pass       => 'ensembl',
    -driver     => 'mysql',
    -DNADB_USER => 'ensro',
    -DNADB_PORT => 3306,
   
    #human
    -species    => 'homo_sapiens',
    -dbname     => 'dev_homo_sapiens_funcgen_69_37',
    -host       => 'ens-genomics2',
    -DNADB_HOST => 'ens-livemirror',
    -DNADB_NAME => 'homo_sapiens_core_67_37',
    
    #mouse
    #-species    => 'mus_musculus',
    #-host       => 'ens-genomics1',
    #-dbname     => 'dev_mus_musculus_funcgen_69_38',
    #-DNADB_HOST => 'ens-genomics1',
    #-DNADB_NAME => 'ia3_mus_musculus_core_68_38',

);


##define data root
my $data_root =
  "/lustre/scratch103/ensembl/funcgen/output/dev_homo_sapiens_funcgen_69_37/";

#my $data_root="/lustre/scratch103/ensembl/funcgen/output/dev_mus_musculus_funcgen_69_38/";



my $rsa = $efgdba->get_adaptor("resultset");
$rsa->dbfile_data_root($data_root);

#my $ftype_adaptor = $efgdba->get_adaptor('featuretype');
#my $ftype         = $ftype_adaptor->fetch_by_name('5mC');

#my @result_sets = @{ $rsa->fetch_all_by_FeatureType($ftype) };

my @result_sets = @{ $rsa->fetch_all_by_feature_class('dna_methylation', {status => 'DISPLAYABLE'}) };


plan tests => scalar @result_sets * 3;

my $dnaa          = $efgdba->get_adaptor("DNAMethylationFeature");
my $slice_adaptor = $efgdba->get_adaptor("slice");
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 1, 1, 3011550 );

foreach my $resultset (@result_sets) {
    print "\ntesting resultset " . $resultset->name . "...\n";


    
    #This is not a a ResultSetAdaptor test
    #change this to test the DNAMethylationFeatire::result_set method
    #ok(
    #    $resultset->isa('Bio::EnsEMBL::Funcgen::ResultSet'),
    #    "the object belongs to the right class - ResultSet"
    #);

    my $dna_meth_features = $dnaa->fetch_all_by_Slice_ResultSet($slice, $resultset);

    #generate random index to test one DNAMethylationFeature object
    my $index = int( rand( scalar @{$dna_meth_features} ) );
    my $df    = $dna_meth_features->[$index];

    ok(
        $df->isa('Bio::EnsEMBL::Funcgen::DNAMethylationFeature'),
        "the object belongs to the right class - DNAMethylationFeature"
    );


    #This is not quite true as there maybe many other DNA meth analyses

    ok(
        (
                 $df->analysis->logic_name eq 'RRBS_merged_filtered_10'
              or $df->analysis->logic_name eq 'Whole_Genome_methylC_seq'
        ),
        "the analysis method is correct: " . $df->analysis->logic_name
    );


    #This also needs to match seq from DB!
    #only need to match context to DB seq
    #testing the actual seq is a data test, not a code test
    ok( $df->context =~ m/CG|CHG|CHH/,
        "the sequence context is correct - " . $df->context );

    #Should test feature_type class is of DNAModification?
    #Altho this is a data test, not a code test


    #add tests for fetch method constraint methods
    #$params = {min_read_depth => 5, context => 'CT', percent_methylation => 50};

}

print "\n\nTested following result_sets...\n\n";

map { print $_->name, "\t..ok\n" } @result_sets;

print "\n\n\nTested "
  . scalar @result_sets
  . " datasets for DNAMethylationFeature\n\ndatabase name: "
  . $efgdba->dbc->dbname . "\n"
  . "database host: "
  . $efgdba->dbc->host . "\n";
 
















######
=i
my $efgdba_mouse = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => 'ensadmin',
    -pass       => 'ensembl',
    -dbname     => 'dev_mus_musculus_funcgen_69_38',
    -host       => 'ens-genomics1',
    -species    => 'mus_musculus',
    -driver     => 'mysql',
    -DNADB_USER => 'ensro',
    -DNADB_HOST => 'ens-genomics1',
    -DNADB_NAME => 'ia3_mus_musculus_core_68_38',
    -DNADB_PORT => 3306
);

my @resultset_names_human =
  qw /Ag04449_5mC_ENCODE_Uw AG04450_5mC_ENCODE_Uw CMK_5mC_ENCODE_Uw
  NHBE_5mC_ENCODE_Uw
  AG09309_5mC_ENCODE_Uw AG09319_5mC_ENCODE_Uw AG10803_5mC_ENCODE_Uw
  GM12878_5mC_ENCODE_Uw
  GM12891_5mC_ENCODE_Uw GM12892_5mC_ENCODE_Uw GM19239_5mC_ENCODE_Uw
  GM19240_5mC_ENCODE_Uw
  h1ESC_5mC_ENCODE_Uw HCF_5mC_ENCODE_Uw HCM_5mC_ENCODE_Uw
  HCT116_5mC_ENCODE_Uw HEK293_5mC_ENCODE_Uw
  HEPG2_5mC_ENCODE_Uw HMEC_5mC_ENCODE_Uw HRE_5mC_ENCODE_Uw
  HSMM_5mC_ENCODE_Uw HSMMtube_5mC_ENCODE_Uw
  IMR90_5mC_ENCODE_Uw JURKAT_5mC_ENCODE_Uw K562_5mC_ENCODE_Uw
  MCF7_5mC_ENCODE_Uw NB4_5mC_ENCODE_Uw
  SAEC_5mC_ENCODE_Uw SKMC_5mC_ENCODE_Uw SKNSHRA_5mC_ENCODE_Uw
  HCPEpiC_5mC_ENCODE_Uw HEEpiC_5mC_ENCODE_Uw
  HIPEpiC_5mC_ENCODE_Uw HNPCEpiC_5mC_ENCODE_Uw HPAEpiC_5mC_ENCODE_Uw
  HRCEpiC_5mC_ENCODE_Uw HRPEpiC_5mC_ENCODE_Uw
  Melano_5mC_ENCODE_Uw NH-A_5mC_ENCODE_Uw Osteobl_5mC_ENCODE_Uw
  PanIslets_5mC_ENCODE_Uw Fibrobl_5mC_ENCODE_Uw
  HAEpiC_5mC_ENCODE_Uw NHDF-neo_5mC_ENCODE_Uw NT2-D1_5mC_ENCODE_Uw
  H1ESC_5mC_Publication_Lister2009_PMID19829295
  IMR90_5mC_Publication_Lister2009_PMID19829295/;

my @resultset_names_mouse =
  qw /ES_5mC_Publication_Stadler2011_PMID22170606	NPC_5mC_Publication_Stadler2011_PMID22170606/;

my $rsa_mouse = $efgdba_mouse->get_adaptor("resultset");
$rsa_mouse->dbfile_data_root($data_root_mouse);
my $dnaa_mouse = $efgdba_mouse->get_adaptor("DNAMethylationFeature");



map { push( @result_sets, $rsa->fetch_all_by_name($_)->[0] ) }
  (@resultset_names_human);
map { push( @result_sets, $rsa_mouse->fetch_all_by_name($_)->[0] ) }
  (@resultset_names_mouse);
  
  
      my $dnamadaptor =
        $resultset->adaptor->dbc->dbname eq 'dnamet_homo_sapiens_funcgen_69_37'
      ? $dnaa
      : $dnaa_mouse;


=cut






#####################################################################################################################







####works fine#####
=i
my $rsa =$efgdba->get_adaptor("resultset");
$rsa->dbfile_data_root("/lustre/scratch103/ensembl/funcgen/output/dnamet_mus_musculus_funcgen_69_38/");
my @a =@{$rsa->fetch_all_by_name('ES_5mC_Publication_Stadler2011_PMID22170606')};

my $dnaa =$efgdba->get_adaptor("DNAMethylationFeature"); 
$dnaa->load_resultset($a[0]);


my $slice_adaptor = $efgdba->get_adaptor("slice");
my $slice = $slice_adaptor->fetch_by_region('chromosome',1,1,3011550); 

my $dna_meth_features= $dnaa -> get_DNAMethylationFeatures (-SLICE => $slice);

foreach my $df (@{$dna_meth_features})
{
print "defined keys.....\n\n";
foreach my $key (keys (%$df))
{
print $key . "=>" . $df->{$key} ."\t" if defined $df->{$key};
#print "cell type " 
print $df->methylated_reads, "\t", $df->unmethylated_reads if $key eq 'methylated_reads';
print $df->total_reads if $key eq 'total_reads';
print $df->percent_methylation if $key eq 'percent_methylation';
print $df->context if $key eq 'context';
print $df->get_bigbed_adaptor if $key eq 'bigbed_adaptor';print "\n";

}
print "display_label  " . $df->display_label ."\n";
print "Analysis logic name  " . $df->analysis->logic_name ."\n";
print "chr",$df->seq_region_name,":",$df->seq_region_start,"-",$df->seq_region_end,"\n";
print "__________________________________\n\n\n";
}


=cut


###############
=i

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
$efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor -> new (-user   => 'ensadmin', -pass => 'ensembl', -dbname => 'dnamet_mus_musculus_funcgen_69_38',-host   => 'ens-genomics1', -species => 'mus_musculus',-driver => 'mysql', -DNADB_USER   => 'ensro', -DNADB_HOST   => 'ens-genomics1', -DNADB_NAME   => 'ia3_mus_musculus_core_68_38', -DNADB_PORT => 3306 );
$rsa =$efgdba->get_adaptor("resultset");                                                                                                                                                   
@a =@{$rsa->fetch_all_by_name('ES_5mC_Stadler2011_PMID22170606')};                                                                                                                         
$dnaa =$efgdba->get_adaptor("DNAMethylationFeature");                                                                                                                                      
$dnaa->load_resultset($a[0]);                                                                                                                                                              
$slice_adaptor = $efgdba->get_adaptor("slice");                                                                                                                                            
$slice = $slice_adaptor->fetch_by_region('chromosome',1,3010820,61265651);                                                                                                                  
@dna_meth_faetures= @{$dnaa -> get_DNAMethylationFeatures (-SLICE => $slice )};
                                                                                         
$df=$dna_meth_faetures[0];                                     
print $df->cell_type->name ;
print $df->feature_type->name ;
print $df->context;
print $df->methylated_reads;
print $df->unmethylated_reads;
print $df->total_reads;
print $df->percent_methylation;
print $df->methylation_level;
print $df->strand;
print $df->display_label;
print $df->set->adaptor->dbc->dbname
  

#use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor;
#use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;

#my $dnamethfeatadaptor='Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor';


my $dnamethfeatadaptor= Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor->new (-BIG_BED => '/lustre/scratch101/ensembl/funcgen/ia3/test/fastq/mus_musculus/stadler_downloaded/BiSeq/BigBed/GSE30202_BisSeq_ES_CpGmeth_sorted1.bb');


my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
    #The Registry automatically picks port 5306 for ensembl db
    #-verbose => 1, #specificy verbose to see exactly what it loaded
#);

my $slice_adaptor = $registry->get_adaptor("Mouse","funcgen", "slice");
my $fset_adaptor = $registry->get_adaptor("Mouse","funcgen","featureset");

my $slice = $slice_adaptor->fetch_by_region('chromosome',1,3000574,3002599);

my $fset = $fset_adaptor -> fetch_by_name ("ES_H3K4me3_Mikkelsen2007_PMID17603471_SWEmbl_R015_D150");
#my $display_label = "ES_Biseq_

#DNAMethylationFeatureAdaptor
my $list_head =  $dnamethfeatadaptor -> get_DNAMethylationFeatures (-SLICE => $slice , -NAME => "ES_Biseq" , -STRAND => 0);#, -FEATURE_SET   => $fset);

print "dnamethfeatadaptor is\t" , ref $dnamethfeatadaptor ," and the adaptor is\t", ref $dnamethfeatadaptor->bigBed_File_Adaptor () ,"\n";


#$list_head = $test->bigBed_File_Adaptor -> bigBedIntervalQuery(1,3000574,3002599); # the start is same as slice start and end is slice end +1. the return cordinates do not need to be changed. Note bigBed cordinates should be 0 based

print ref $list_head,"\t", scalar @$list_head ,"\n";



for (my $i=$list_head->head;$i;$i=$i->next) {
         print join("\t",$i->start,$i->end,$i->rest),"\n";
     }



foreach my $df (@{$list_head})
{
print "defined keys.....\n\n";
foreach my $key (keys (%$df))
{
print $key . "=>" . $df->{$key} ."\t" if defined $df->{$key};
print $df->methylated_reads, "\t", $df->unmethylated_reads if $key eq 'methylated_reads';
print $df->total_reads if $key eq 'total_reads';
print $df->methylation_level if $key eq 'methylation_level';
print $df->context if $key eq 'context';
print $df->get_bigbed_adaptor if $key eq 'bigbed_adaptor';print "\n";
}
print "chr",$df->seq_region_name,":",$df->seq_region_start,"-",$df->seq_region_end,"\n";
print "__________________________________\n\n\n";
}


foreach my $df (@{$list_head})
{
print "undefined keys.....\n\n";
foreach my $key (keys (%$df))
{
print $key ."=>\n" if !defined $df->{$key};
}
print "/////////////////////////////////\n\n\n";
}
=cut
