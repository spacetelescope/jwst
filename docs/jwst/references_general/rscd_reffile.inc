.. _rscd_reffile:

RSCD Reference File
-------------------

:REFTYPE: RSCD
:Data model: `~jwst.datamodels.RSCDModel`

The RSCD reference file contains coefficients used to compute the
correction.

.. include:: ../references_general/rscd_selection.inc

.. include:: ../includes/standard_keywords.inc

Type Specific Keywords for RSCD
+++++++++++++++++++++++++++++++
In addition to the standard reference file keywords listed above,
the following keywords are *required* in RSCD reference files,
because they are used as CRDS selectors
(see :ref:`rscd_selectors`):

=========  ==============================
Keyword    Data Model Name
=========  ==============================
DETECTOR   model.meta.instrument.detector
=========  ==============================

Reference File Format
+++++++++++++++++++++
RSCD reference files are FITS format, with 1 BINTABLE extension.
The FITS primary HDU does not contain a data array.
The BINTABLE extension uses the identifier EXTNAME = "RSCD" and
the characteristics of the table columns are as follows:

===========  ==========  =============================================
Column name  Data type   Notes
===========  ==========  =============================================
subarray     char\*13    FULL or subarray name
readpatt     char\*4     SLOW or FAST
rows         char\*4     EVEN or ODD
tau          float32     e-folding time scale, in units of frames
ascale       float32     :math:`b{1}` in equation 2
pow          float32     :math:`b{2}` in equation 2
illum_zp     float32     :math:`illum_{zpt}` in equation 2.1
illum_slope  float32     :math:`illum_{slope}` in equation 2.1
illum2       float32     :math:`illum2` in equation 2.1
param3       float32     :math:`b{3}` in equation 2
crossopt     float32     :math:`Crossover Point` in equation 2.2
sat_zp       float32     :math:`sat_\text{zp}` in equation 3.1
sat_slope    float32     :math:`sat_\text{slope}` in equation 3.1
sat_2        float32     :math:`sat_2` in equation 3.1
sat_mzp      float32     :math:`sat_\text{mzp}` in equation 3
sat_rowterm  float32     :math:`evenrow_{corrections}` in equation 3.1
sat_scale    float32     :math:`scale_\text{sat}` in equation 3
===========  ==========  =============================================

The entries in the first 3 columns of the table are used as row selection criteria, matching
the exposure properties and row type of the data. The remaining 14 columns contain the
parameter values for the double-exponential correction function.

The general form of  the correction to be added to the input data is::

   corrected data = input data + dn_accumulated * scale * exp(-T / tau)  (Equation 1)

where:

 - T is the time since the last group in the previous integration
 - tau is the exponential time constant found in the RSCD table
 - dn_accumulated is the DN level that was accumulated for the pixel from the previous integration.

In cases where the last integration does not saturate, the :math:`scale` term in equation 1
is determined according to:

  :math:`scale = b{1}* [Counts{2}^{b{2}} * [1/exp(Counts{2}/b{3}) -1 ]\; \; (Equation \; 2)`

The following two additional equations are used in Equation 2:

  :math:`b{1} = ascale * (illum_{zpt} + illum_{slope}*N + illum2* N^2) \; \; (Equation \; 2.1)`

  :math:`Counts{2} = Final \, DN \, in \, the \,  last \, group \, in \; the \, previous \, integration \, - Crossover \, Point \; \; (Equation \; 2.2)`
      
where:

 - N in equation 2.1 is the number of groups per integration
 - Crossover Point in equation 2.2 is column CROSSOPT in the RSCD table.

If the previous integration saturates, :math:`scale` is no longer calculated using equations
2 - 2.2. Instead it is calculated using equations 3 and 3.1:

   :math:`scale_\text{sat} = slope * Counts{3} + sat_\text{mzp} \; \; (Equation \; 3)`
 
   :math:`slope = sat_{zp} + sat_{slope} * N + sat_2*N^2 + evenrow_{corrections} \; \; (Equation \; 3.1)`

where:

    - :math:`Counts{3}` is an estimate of what the last group in the
      previous integration would have been if saturation did not exist
    - :math:`scale_\text{sat}` is sat_scale in the RSCD table
    - :math:`sat_\text{mzp}` is sat_mzp in the RSCD table
    - :math:`sat_\text{zp}` is sat_zp in the RSCD table
    - :math:`sat_\text{slope}` is sat_slope in the RSCD table
    - :math:`sat_2` is sat2 in the RSCD table
    - :math:`evenrow_{corrections}` is sat_rowterm in the RSCD table
    - N is the number of groups per integration
 
