Compute celltype-specific TF binding sites
###########################################

Identify transcription factor (TF) binding sites in the genome by running WiggleTools on TF peaks.

.. figure:: images/rf_compute_celltype_tf_sites.png
   :scale: 50 %
   :alt: workflow compute_celltype_tf_sites

   Workflow showing how transcription factor bidning sites are calculated.



.. _rf_m_compute_celltype_tf_sites:

compute_celltype_tf_sites
-------------------------
Input
  * $options->{cell_type_tfs}
    Contains cell types TFs in some format

Creates two directories, *celltype_tf* and *celltype_dnase*.
Iterates through each cell type defined in $options->{cell_type_tfs} (see get_metadata) and calls *compute_celltype_tf_sites_2*.


.. _rf_m_compute_celltype_tf_sites_2:

compute_celltype_tf_sites_2
---------------------------

Input
  * Transcription Factor peaks
  * Open Chromantin peaks [optional]

Output
  * working_dir/celltype_tf/$celltype.bed
  * working_dir/celltype_dnase/$celltype.bed [condition: open chromatin peaks]

If open chromatin data exists for this epigenome, create bedgraph containing all open areas. Then overlap with TF data and find all areas where open chromatin and TF signal overlaps (sum > 1). Could also be where either exist.

If not open chromatin data exist, run only on those.

All files are trimmed to chromosome length after being created. 