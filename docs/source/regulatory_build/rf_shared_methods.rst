Shared Methods
##############

Methods used multiple times in the script which are mostly self-explanatory

.. _rf_sm_must_compute:

must_compute
""""""""""""
Decides whether a computation is needed
    * As long as all checked files exist, they are not recomputed
    * As soon as a checked file needs to be recomputed, everything after that is recomputed


.. _rf_sm_convert_to_bigBed:

convert_to_bigBed
"""""""""""""""""


.. _rf_sm_run:

run
"""
 Wrapper function for system calls. Uses elements of EnsEMBL::Hive.