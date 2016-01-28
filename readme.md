# ensembl-funcgen
The Ensembl Funcgen API (Application Programme Interface) serves as a middle layer between the underlying MySQL database and aims to encapsulate the database layout by providing high level access to the database. The team providing the data is Ensembl Regulation, but for historic reasons, the API and databases are called Funcgen.

## Content
The Funcgen database contains currently 4 different types of data which can be accessed through the API.

#### 1. Regulatory Features
Using our Regulatory Build, which is described here: http://www.ensembl.org/info/genome/funcgen/regulatory_build.html, we identify regions which can be involved in gene transcript regulation. The build process uses publicly available data sets based on ChIP-Seq, DNAse1-Seq and Transcription factor binding motifs.

#### 2. Segmentation
Ensembl also provides the complementary genome segmentation analyses which are used as inputs to the Ensembl Regulatory Build: http://www.ensembl.org/info/genome/funcgen/regulatory_segmentation.html

#### 3. Microarray Probe Mappings
Microarray probe mappings for several species, for example Affymetrix, Codelink, Agilent, Illumina and Phalanx.

#### 4. External Regulatory Data
The Funcgen database also contains data loaded from external sources such as: Tarbase (human, mouse); human enhancers from the VISTA Browser; human RRBS and WGBS assays and the Fantom 5 human regulatory element predictions.


More information can be found on our blog, including details of the recent regulatory build redesign process:
http://www.ensembl.info/blog/2013/12/26/the-new-ensembl-regulatory-annotation/
http://www.ensembl.info/blog/2014/01/01/computing-ensembls-new-regulatory-annotation/

## Download
To clone the Ensembl Regulation API, use the following command:

```
git clone https://github.com/Ensembl/ensembl-funcgen.git
```

## API requirements
In order to use the Regulation API, an installation of the same version of the Ensembl-core API is mandatory. A guide for installing all Ensembl APIs and their respective prerequisites is available here:
http://www.ensembl.org/info/docs/api/api_installation.html

## API tutorial
An extensive tutorial for the Regulation API can be found at:

http://www.ensembl.org/info/docs/api/funcgen/regulation_tutorial.html

## Contributions

If you wish to contribute to this repository or any Ensembl repository, please refer to [our contribution guide](https://github.com/Ensembl/ensembl/blob/master/CONTRIBUTING.md).


## Contact us
Please email comments or questions to the public Ensembl developers list at 

<http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at

<http://www.ensembl.org/Help/Contact>.
