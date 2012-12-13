=head1 DESCRIPTION

This is an example script for the DNAMethylationFeature API.  This class is supported by the database and also by
UCSC BigBed format files, and hence has some additional requirements:

    ensembl-external       http://www.ensembl.org/info/docs/api/api_cvs.html
    ensembl-webcode        http://www.ensembl.org/info/docs/api/api_cvs.html
    Bio-BigFile-1.07       http://search.cpan.org/~lds/Bio-BigFile-1.07/lib/Bio/DB/BigFile.pm
    UCSC Kent source tree  http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
      or more docs here    https://sites.google.com/site/xzhou82/Home/bioinformatics-computational-biology/ucsc-genome-browser

The Ensembl methylation data can be downloaded from the Ensembl FTP site via the 'data_files' link available here:
    www.ensembl.org/info/data/ftp/index.html
Or more preferably direct access e.g.:
    ftp.ensembl.org/pub/release-69/data_files/homo_sapiens/GRCh37/dna_methylation_feature/

It is important that the sub directory structure provided here for the files is maintained i.e. dna_methylation_feature and 
all of it's contents. This is because these paths are referenced from the database table 'dbfile_registry'. More information 
regarding the file format is available in the README file on the FTP site.

For API access, the location (or dbfile_data_root) of the dna_methylation_feature directory will require configuration in the 
ResultSetAdaptor.

e.g. if you download the files to the following directory:

    /home/ensembl/data_files/homo_sapiens/GRCh37/dna_methylation_feature

<<<<<<< dna_methylation_features.pl
Then

   my $dbfile_data_root='/home/ensembl/data_files/homo_sapiens/GRCh37';

See below for how this is set and how features are retrieved.
=======
Your dbfile_data_root is /home/ensembl/data_files/homo_sapiens/GRCh37/
>>>>>>> 1.2

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry         = 'Bio::EnsEMBL::Registry';
<<<<<<< dna_methylation_features.pl
my $dbfile_data_root = '/home/ensembl/data_files/homo_sapiens/GRCh37/';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous'
  );

my $dmf_adaptor        = $registry->get_adaptor( 'Human', 'funcgen', 'DNAMethylationFeature' );
my $slice_adaptor      = $registry->get_adaptor( 'Human', 'core', "slice" );
my $result_set_adaptor = $registry->get_adaptor( 'Human', 'funcgen', "resultset" );
=======

#change for db_file_data_root
my $dbfile_data_root = '/nfs/ensnfs-dev/staging/homo_sapiens/GRCh37/';
>>>>>>> 1.2

#Set the dbfile_data_root!
$result_set_adaptor->dbfile_data_root($dbfile_data_root);


#fetch ResultSet by name
my ($rset) = @{ $result_set_adaptor->fetch_all_by_name('AG04449_5mC_ENCODE_Uw') };

#pick a genomic region of interest
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 1, 1, 3011550 );

#Optionally define some contraints
my $constraints = 
  {
   min_read_depth  => 15,
   context         => 'CG',
   min_methylation => 25,
   max_methylation => '75.5',
  };

#Fetch the features
my $dna_meth_features = $dmf_adaptor->fetch_all_by_Slice_ResultSet( $slice, $rset, $constraints );

#Print some information
foreach my $dmf ( @{$dna_meth_features} ) {
    print "\n-------------------------------------\n";
    print "Display label:\t" . $dmf->display_label . "\n";
    print "Location:\t" . $dmf->feature_Slice->name . "\n";
    print "Methylated reads:\t" . $dmf->methylated_reads . "\n";
    print "Total reads:\t" . $dmf->total_reads . "\n";
    print "Percent Methylation:\t" . $dmf->percent_methylation . "\n";
    print "Context:\t" . $dmf->context . "\n";
}
