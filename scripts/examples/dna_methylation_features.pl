=head1 DESCRIPTION

This is an example script to retrieve and use DNAMethylationFeature objects from UCSC BigBed file format. 
The Ensembl methylation data can be downloaded via FTP from the Ensembl site (http://www.ensembl.org/info/data/ftp/index.html ). 
The BigBed files (.bb) in this directory represent processed DNA methylation information from whole genome (WGBS) and reduced 
representation (RRBS) bisulphite sequencing experiments. 
For information about BigBed files please see http://genome.ucsc.edu/FAQ/FAQformat.html. Each line of the source bed file 
represents one genomic cytosine (DNAMethylationFeature). The fourth field of each line contains information about cytosine 
context and total coverage, and the fifth field gives percentage of methylated reads multiplied by a factor of 10 that generates 
a methylation score value for the DNAMethylationFeature and ranges from 0-1000. Strand orientation of the feature is present in the 
sixth field. An example line that represents a cytosine in CG context with 19 methylated reads out of a total of 20 reads would be as under:

chr1    3000573 3000574 CG/20   950     +

Access to these data requires installation of EnsEMBL::ExternalData modules from  Ensembl CVS and Bio::DB::BigFile from CPAN. 
DNAMethylationFeatureAdaptor::fetch_all_by_Slice_ResultSet API method can be used to fetch the data. This requires that the 
ResultSetAdaptor:dbfile_data_root is set explicitly before fetching any ResultSets.

It is important that the directory structure provided here for the files
is maintained wherever you wish to store them, i.e they must be kept in

species_name/assembly_name/feature_class/

e.g. if you download the files to the following directory:

/home/ensembl/data_files/homo_sapiens/GRCh37/dna_methylation_feature 

Your dbfile_data_root is /home/ensembl/data_files/

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry         = 'Bio::EnsEMBL::Registry';
my $dbfile_data_root = '/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $dmf_adaptor =
  $registry->get_adaptor( 'Human', 'funcgen', 'DNAMethylationFeature' );
my $slice_adaptor = $registry->get_adaptor( 'Human', 'core', "slice" );
my $result_set_adaptor =
  $registry->get_adaptor( 'Human', 'funcgen', "resultset" );
$result_set_adaptor->dbfile_data_root($dbfile_data_root);

#fetch all result_sets for DNAMethylation data
my @result_sets = @{
    $result_set_adaptor->fetch_all_by_feature_class( 'dna_methylation',
        { status => 'DISPLAYABLE' } )
};

#fetch result_set by name
my ($rset) =
  @{ $result_set_adaptor->fetch_all_by_name('AG04449_5mC_ENCODE_Uw') };

my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 1, 1, 3011550 );
my $constraints = {
    min_read_depth  => 15,
    context         => 'CG',
    min_methylation => '25',
    max_methylation => '75.5',
};
my $dna_meth_features =
  $dmf_adaptor->fetch_all_by_Slice_ResultSet( $slice, $rset, $constraints );

foreach my $dmf ( @{$dna_meth_features} ) {
    print "\n-------------------------------------\n";
    print "Display label:\t" . $dmf->display_label . "\n";
    print "Location:\t" . $dmf->feature_Slice->name . "\n";
    print "Methylated reads:\t" . $dmf->methylated_reads . "\n";
    print "Total reads:\t" . $dmf->total_reads . "\n";
    print "Percent Methylation:\t" . $dmf->percent_methylation . "\n";
    print "Context:\t" . $dmf->context . "\n";
}