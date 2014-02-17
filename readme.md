# ensembl-funcgen
The Ensembl Funcgen API (Application Programme Interface) serves as a middle layer between the underlying MySQL database and aims to encapsulate the database layout by providing high level access to the database. The team providing the data is Ensembl Regulation, but for historic reasons, the API and databases are called Funcgen.

## Content
The Regulation database contains currently 4 different types of data which can be accessed through the API
#### 1. Regulatory Features
Using our Regulatory Build, which is described here: http://www.ensembl.org/info/genome/funcgen/regulatory_build.html, we identify regions which can be involved in gene transcript regulation. The build process uses publicly available data sets based on ChIP-Seq, DNAse1-Seq and Transcription factor binding motifs.

#### 2. Segmentation
Ensembl also provides genome segmentation analyses for six human cells lines from the ENCODE project, an analogous summary of a region's regulatory function.

#### 3. Microarray Probe Mappings
Microarray probe mappings for several species, for example Affymetrix, Codelink, Agilent and Illumina.

#### 4. External Regulatory data
The Regulation database also contains data loaded from external sources like Tarbase (human, mouse), human enhancers from the VISTA Browser, human RRBS and WGBS assays, and regulatory element predictions from cisRED (human, mouse) to name some.

More detailed information about the Regulatory Build can be found at: 
http://www.ensembl.org/info/genome/funcgen/regulatory_build.html

Please also note that we are redesigning our build process. More information about this can be found in our blog:

http://www.ensembl.info/blog/2013/12/26/the-new-ensembl-regulatory-annotation/
and

http://www.ensembl.info/blog/2014/01/01/computing-ensembls-new-regulatory-annotation/

## Download
To clone the Ensembl Regulation API, use the following command:

```
git clone https://github.com/Ensembl/ensembl-funcgen.git
```

## API requirements
In order to use the Regulation API, an installation of the same version of the Ensembl-core API is mandatory. A guide for installing all Ensembl APIs and their respective prerequisites  is available here:
http://www.ensembl.org/info/docs/api/api_installation.html

## API tutorial
An extensive tutorial for the Regulation API can be found at:
http://www.ensembl.org/info/docs/api/funcgen/regulation_tutorial.html


## Contact us
Please email comments or questions to the public Ensembl developers list at 
<http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.
