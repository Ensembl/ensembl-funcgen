Overview
########

The Ensembl Regulatory Build (`Zerbino et al. 2015 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0621-5>`_) contains a genome-wide set of regions that are likely to be involved in gene regulation. 

The process is divided into

* Registration of metadata
* Download and QC of sequence data
* Alignment to the genome
* Peak calling to identify enriched regions
* Segmentation


* The script was written in a time when epigenomes where cell types.
* Variables assigned to the options hash should not be changed once assigned.
