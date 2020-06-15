Regulatory Features
###################

The Ensembl Regulatory Build (Zerbino et al. 2015) contains a genome-wide set of regions that are likely to be involved in gene regulation. 

The process is divided into

* Registration of metadata
* Download and QC of sequence data
* Alignment to the genome
* Peak calling to identify enriched regions
* Segmentation


* The script was written in a time when epigenomes where cell types.
* Variables assigned to the options hash should not be changed once assigned.

.. toctree::
  :hidden:

  rf_compute_celltype_tf_sites
  rf_shared_methods
  rf_directories
  rf_input
  rf_tools


