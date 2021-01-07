Description
============

``exp_to_source`` is a data reorganization tool that is used to convert
Stage 2 exposure-based data products to Stage 3 source-based data products.
It is only used when there is a known source list for the exposure data,
which is required in order to reorganize the data by source. Hence it is
only useable for NIRSpec MOS, NIRSpec fixed-slit, NIRCam WFSS, and NIRISS
WFSS data. Details on the operation for each mode are given below.

The tool is run near the beginning of the :ref:`calwebb_spec3 <calwebb_spec3>`
pipeline, immediately after the :ref:`master background <master_background_step>`
step.

In general, the tool takes as input multiple exposure-based "cal" products
created during Stage 2 spectroscopic (:ref:`calwebb_spec2 <calwebb_spec2>`)
processing and reorganizes the data in them to create a set of output
source-based "cal" products, which are then processed through the remaining
steps of the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline. For example,
if the input consists of a set of 3 exposure-based "cal" files (from a
3-point nod dither pattern, for example), each one of which contains data
for 5 defined sources, then the output consists of a set of 5
source-based "cal" products (one per source), each of which contains the
data from the 3 exposures for each source. This is done as a convenience,
in order to have all the data for a given source contained in a single
product. All data arrays associated with a given source, e.g. SCI, ERR, DQ,
WAVELENGTH, VAR_POISSON, etc., are copied from each input product into
the corresponding output product.

NIRSpec MOS
^^^^^^^^^^^
NIRSpec MOS observations are created at the APT level by defining a
configuration of MSA slitlets with a source assigned to each slitlet.
The source-to-slitlet linkage is carried along in the information contained
in the MSA metadata file used during :ref:`calwebb_spec2 <calwebb_spec2>`
processing. Each slitlet instance created by the :ref:`extract_2d <extract_2d_step>`
step stores the source ID (a simple integer number) in the SOURCEID keyword of
the SCI extension header for the slitlet. The ``exp_to_source`` tool uses
the SOURCEID values to sort the data from each input product into an
appropriate source-based output product.

NIRSpec Fixed-Slit
^^^^^^^^^^^^^^^^^^
NIRSpec fixed-slit observations do not have sources identified with each
slit, so the slit names, e.g. S200A1, S1600A1, etc., are mapped to predefined
source ID values, as follows:

=========  =========
Slit Name  Source ID
=========  =========
S200A1         1
S200A2         2
S400A1         3
S1600A1        4
S200B1         5
=========  =========

The assigned source ID values are used by ``exp_to_source`` to sort the data
from each slit into the appropriate source-based output product.

NIRCam and NIRISS WFSS
^^^^^^^^^^^^^^^^^^^^^^
Wide-Field Slitless Spectroscopy (WFSS) modes do not have a predefined
source list, but a source list is created by the
:ref:`calwebb_image3 <calwebb_image3>` pipeline when it processes the
direct image exposures that are included in a WFSS observation. That
source catalog is then used during :ref:`calwebb_spec2 <calwebb_spec2>`
processing of the grism exposures to define and create cutouts for each
identified source. Like NIRSpec MOS mode, each "slit" instance is identified
by the source ID number from the catalog and is used by the ``exp_to_source``
tool to reorganize exposures into source-based products.

File Name Syntax
^^^^^^^^^^^^^^^^
The input exposure-based "cal" products have file names that follow the
standard Stage 2 exposure syntax, such as::

 jw93065002001_02101_00001_nrs1_cal.fits

See :ref:`exposure-based file names <exp_file_names>` for more details
on the meaning of each field in the file names.

The FITS file naming scheme for the source-based "cal" products follows
the standard Stage 3 syntax, such as::

 jw10006-o010_s00061_nirspec_f170lp-g235m_cal.fits

where "s00061" in this example is the source ID.
See :ref:`source-based file names <src_file_names>` for more details
on the meaning of each field in this type of file name.
