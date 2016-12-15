package Bio::EnsEMBL::Funcgen::Hive::Ftp::WriteReadmeFile;

use strict;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self    = shift;
  my $destination = $self->param('destination');
  
  use File::Basename;
  my $destination_dir = dirname($destination);

  use File::Path qw(make_path);
  make_path($destination_dir);

  use IO::File;
  my $output_fh = IO::File->new(">$destination");
  
  my $readme_content = $self->generate_readme_content();

  $output_fh->print($readme_content);
}

sub generate_readme_content {

  my $self = shift;
  
  return <<README_FILE_CONTENT
The data files on the regulation FTP site follow a naming convention, which is described 
in greater detail here:

https://github.com/FAANG/faang-metadata/blob/master/docs/faang_analysis_metadata.md#file-naming

The file names include the following information separated with a dot ('.'):

  - species,
  - assembly version,
  - cell type (if applicable),
  - feature type (if applicable),
  - analysis name,
  - results type,
  - data freeze date and
  - file format.

E.g.: homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20161111.gff.gz

The data available on our FTP site are:

Peaks
-----

The set of peaks for transcription factors, histone modifications and 
variants that are part of our regulatory resources. In previous releases 
these used to be collated in one file, called ‘AnnotatedFeatures.gff.gz’, but 
with our recent expansion to 88 human cell types with ChIP-seq data, the file 
became too big. Therefore, we split it into separate files by cell and 
feature type in the ‘Peaks’ subdirectory. The peaks are now available in gff, 
bed and bigBed format.

Quality scores
--------------

The outcome of our quality checks from processing the ChiP-seq data that 
yielded the peaks. They are in JSON format in the 'QualityChecks'
subdirectory:

  - the number of mapped reads, 
  - the estimated fragment length, the NSC and RSC values using 
    phantompeakqualtools (https://www.encodeproject.org/software/phantompeakqualtools),
  - the proportion of reads in peaks and 
  - the enrichment of the ChIP over the Input using CHANCE 
    (https://www.ncbi.nlm.nih.gov/pubmed/23068444).

Regulatory build
----------------

The current set of regulatory features along with their predicted activity in 
every cell type. We provide one gff file per cell type in the 
'regulatory_features' subdirectory.

Transcription factor motifs
---------------------------

The transcription factor motifs identified using position weight matrices 
from JASPAR (http://jaspar.genereg.net) in enriched regions identified by our ChIP-seq analysis pipeline 
in gff format.


README_FILE_CONTENT

}

1;
