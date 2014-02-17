# ensembl-funcgen
The Ensembl Funcgen API (Application Programme Interface) serves as a middle layer between the underlying MySQL database and aims to encapsulate the database layout by providing high level access to the database. The team providing the data is Ensembl Regulation, but for historic reasons, the API and databases are called Funcgen.

## Content
The Regulation database contains 4 different types of data:
### 1. Regulatory Features
  Using our Regulatory Build, which is described here: http://www.ensembl.org/info/genome/funcgen/regulatory_build.html, we identify regions which can be involved in gene transcript regulation. The build process uses publicly available data sets based on ChIP-Seq, DNAse1-Seq and Transcription factor binding motifs.

2. Seqmentation

3. Microarray Probe Mappings

4. Other / External Regulatory data


## Download
To clone the Ensembl Regulation API, use the following command:

```
git clone https://github.com/Ensembl/ensembl-funcgen.git
```

## API requirements
In order to use the Regulation API, an installation of the same version of the Ensembl-core API is mandatory. A guide for installing all Ensembl APIs and their prerequisits is available here:
http://www.ensembl.org/info/docs/api/api_installation.html

## API tutorial
An extensive tutorial for the Regulation API can be found at:
http://www.ensembl.org/info/docs/api/funcgen/regulation_tutorial.html


## Contact us
Please email comments or questions to the public Ensembl developers list at 
<http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.
