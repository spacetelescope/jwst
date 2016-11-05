Introduction
============

Purpose and Scope
-----------------
This document specifies the format of each calibration reference file used by the 
JWST Calibration pipeline for DMS Build 7, satisfying requirement DMS-653
("The format of each calibration reference file shall be specified in the JWST 
Calibration Reference File Specification Document.")  Many calibration steps in 
the DMS Build 7 Calibration Pipeline require reference files retrieved from CRDS.  
This document is intended to be a reference guide to the formats of reference 
files for steps requiring them, and is not intended to be a detailed description of each of 
those pipeline steps.

Data Quality Flags
------------------
Within science data files, the PIXELDQ flags are stored as 32-bit integers; 
the GROUPDQ flags are 8-bit integers.  The meaning of each bit is specified 
in a separate binary table extension called DQ_DEF.  The binary table has the 
format presented in Table 1, which represents the master list of DQ flags.  
Only the first eight entries in the table below are relevant to the 
GROUPDQ array. All calibrated data from a particular instrument and observing mode 
have the same set of DQ flags in the same (bit) order. For Build 7, this master 
list will be used to impose this uniformity.  We may eventually use different master 
lists for different instruments or observing modes.


Within reference files for some steps, the Data Quality arrays for some steps are 
stored as 8-bit integers to conserve memory.  Only the flags actually used by a reference 
file are included in its DQ array.  The meaning of each bit in the DQ array is stored in 
the DQ_DEF extension, which is a binary table having the following fields: Bit, Value, 
Name, and Description.


Table 1. Flags for the PIXELDQ and GROUPDQ Arrays (Format of DQ_DEF Extension)

===  ==========    ================  ===========================================
Bit  Value         Name              Description
===  ==========    ================  ===========================================
0    1	           DO_NOT_USE        Bad pixel. Do not use.
1    2             SATURATED         Pixel saturated during exposure
2    4             JUMP_DET          Jump detected during exposure
3    8             DROPOUT           Data lost in transmission
4    16            RESERVED	 
5    32	           RESERVED	 
6    64            RESERVED	 
7    128           RESERVED	 
8    256           UNRELIABLE_ERROR  Uncertainty exceeds quoted error
9    512           NON_SCIENCE       Pixel not on science portion of detector
10   1024          DEAD              Dead pixel
11   2048          HOT               Hot pixel
12   4096          WARM              Warm pixel
13   8192          LOW_QE            Low quantum efficiency
14   16384         RC                RC pixel
15   32768         TELEGRAPH         Telegraph pixel
16   65536         NONLINEAR         Pixel highly nonlinear
17   131072        BAD_REF_PIXEL     Reference pixel cannot be used
18   262144        NO_FLAT_FIELD     Flat field cannot be measured
19   524288        NO_GAIN_VALUE     Gain cannot be measured
20   1048576       NO_LIN_CORR       Linearity correction not available
21   2097152       NO_SAT_CHECK      Saturation check not available
22   4194304       UNRELIABLE_BIAS   Bias variance large
23   8388608       UNRELIABLE_DARK   Dark variance large
24   16777216      UNRELIABLE_SLOPE  Slope variance large (i.e., noisy pixel)
25   33554432      UNRELIABLE_FLAT   Flat variance large
26   67108864      OPEN              Open pixel (counts move to adjacent pixels)
27   134217728     ADJ_OPEN          Adjacent to open pixel
28   268435456     UNRELIABLE_RESET  Sensitive to reset anomaly
29   536870912     MSA_FAILED_OPEN   Pixel sees light from failed-open shutter
30   1073741824    OTHER_BAD_PIXEL   A catch-all flag
===  ==========    ================  ===========================================


