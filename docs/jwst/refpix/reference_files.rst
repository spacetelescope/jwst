Reference File Types
====================

The refpix step only uses the refpix reference file when processing
NIRSpec exposures that have been acquired using an IRS2 readout
pattern. No other instruments or exposure modes require a reference
file for this step.


CRDS Selection Criteria
-----------------------
Refpix reference files are selected by DETECTOR and READPATT.


Reference File Format
---------------------

A single extension, with a EXTNAME keyword of 'IRS2', 
provides the complex coefficients for the correction,
and contains 8 columns ALPHA_0, ALPHA_1, ALPHA_2, ALPHA_3, BETA_0, BETA_1,
BETA_2, and BETA_3.  The ALPHA arrays contains correction multipliers to the
data, and the BETA arrays contains correction multiplier to the reference
output. Both arrays have 4 components - one for each sector.
