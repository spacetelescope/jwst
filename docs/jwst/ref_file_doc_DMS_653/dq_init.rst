===========================
Data Quality Initialization
===========================

The DQ initialization step propagates pixel-dependent flags from a static pixel mask reference 
file into the 2-D DQ array of the science data. The 2-D pixel mask is first translated from the 
8-bit DQ array of the reference file into the 32-bit DQ array specified by the master DQ list, 
then propagated into the 2-D DQ array of the science data file using a bit-wise OR operation. 


.. toctree::
   :maxdepth: 4

   dq_init_reference_files
