Introduction
============

This document is intended to be a core reference guide to the formats, naming convention and 
data quality flags used by the reference files for pipeline steps requiring them, and is not
intended to be a detailed description of each of those pipeline steps. It also does not give 
details on pipeline steps that do not use reference files. 
The present manual is referred to by several other documentation pages, 
such as the JWST pipeline and JDocs.

Reference Files Naming Convention
=================================

Before the reference files are ingested into CRDS, they are renamed following a 
convention used by the pipeline. As any other changes undergone by the reference files 
the information of the previous names is kept in the keader keywords, so the Instrument Teams
can easily track which delivered file is being used by the pipeline on each step. 

The HISTORY header keyword of each reference file includes details on specific processing 
undergone by the files before being ingested in CRDS.

The reference files names used by the different pipeline steps are define in the table below.
Details on the files formats used and the CRDS selection criteria are given in the 
correspondent pipeline sections defined in this document.

==================  ======================================== 
Pipeline Step        Name (*)         
==================  ========================================  
**ami_analyze**      jwst_niriss_throughput_xxxx.fits
**assign_wcs**       jwst_nirspec_camera_xxx.rmap
..          	     jwst_nirspec_camera_xxx.rmap
..                   jwst_nirspec_collimator_xxxx.rmap
..                   jwst_nirspec_disperser_xxxx.rmap
..                   jwst_niriss_distortion_xxxx.rmap
..                   jwst_nircam_distortion_xxxx.rmap
..                   jwst_miri_distortion_xxxx.rmap
..                   jwst_miri_filteroffset_xxxx.rmap
..                   jwst_nirspec_fore_xxxx.rmap
..                   jwst_nirspec_fpa_xxxx.rmap
..                   jwst_nirspec_ifufore_xxxx.rmap
..                   jwst_nirspec_ifupost_xxxx.rmap
..                   jwst_nirspec_ifuslicer_xxxx.rmap
..                   jwst_nirspec_msa_xxxx.rmap
..                   jwst_nirspec_ote_xxxx.rmap
..                   jwst_niriss_specwcs_xxxx.asdf
..                   jwst_miri_specwcs_xxxx.asdf
..                   jwst_miri_regions_xxxx.asdf
..                   jwst_miri_v2v3_xxxx.asdf
..                   jwst_miri_wavelengthrange_xxxx.asdf
..                   jwst_nirspec_wavelengthrange_xxxx.asdf
**dark_current**     miri_dark_xxxx.fits
..                   nircam_dark_xxxx.fits 
..                   niriss_dark_xxxx.fits 
..                   nirspec_dark_xxxx.fits
**dq_init**          jwst_miri_mask_xxxx.fits
..                   jwst_nircam_mask_xxxx.fits
..                   jwst_niriss_mask_xxxx.fits
..                   jwst_nirspec_mask_xxxx.fits
**extract_1d**       jwst_miri_extract1d_xxxx.json
..                   jwst_niriss_extract1d_xxxx.json
..                   jwst_nirspec_extract1d_xxxx.json
**flatfield**        jwst_miri_flat_xxxx.fits
..                   jwst_nircam_flat_xxxx.fits
..                   jwst_niriss_flat_xxxx.fits
..                   jwst_nirspec_dflat_xxxx.fits
..                   jwst_nirspec_fflat_xxxx.fits
..                   jwst_nirspec_sflat_xxxx.fits
**fringe**		     jwst_miri_fringe_xxxx.fits
**jump**             jwst_miri_gain_xxxx.fits
..                   jwst_miri_readnoise_xxxx.fits
..                   jwst_nircam_gain_xxxx.fits
..                   jwst_nircam_readnoise_xxxx.fits
..                   jwst_niriss_gain_xxxx.fits
..                   jwst_niriss_readnoise_xxxx.fits
..                   jwst_nirspec_gain_xxxx.fits
..                   jwst_nirspec_readnoise_xxxx.fits
**linearity**        jwst_miri_linearity_xxxx.fits 
..                   jwst_nircam_linearity_xxxx.fits
..                   jwst_niriss_linearity_xxxx.fits
..                   jwst_nirspec_linearity_xxxx.fits
**pathloss**         jwst_nirspec_pathloss_xxxx.fits
**photom**           jwst_miri_photom_xxxx.fits
..                   jwst_miri_area_xxxx.fits
..                   jwst_nircam_photom_xxxx.fits
..                   jwst_nircam_area_xxxx.fits
..                   jwst_niriss_photom_xxxx.fits
..                   jwst_niriss_area_xxxx.fits
..                   jwst_nirspec_photom_xxxx.fits
..                   jwst_nirspec_area_xxxx.fits
**ramp_fitting**     jwst_miri_gain_xxxx.fits
..                   jwst_miri_readnoise_xxxx.fits
..                   jwst_nircam_gain_xxxx.fits
..                   jwst_nircam_readnoise_xxxx.fits
..                   jwst_niriss_gain_xxxx.fits
..                   jwst_niriss_readnoise_xxxx.fits
..                   jwst_nirspec_gain_xxxx.fits
..                   jwst_nirspec_readnoise_xxxx.fits
**refpix**           jwst_nirspec_refpix_xxxx.fits
**rcsd**             jwst_miri_rscd_xxxx.fits
**saturation**       jwst_miri_saturation_xxxx.fits
..                   jwst_nircam_saturation_xxxx.fits
..                   jwst_niriss_saturation_xxxx.fits
..                   jwst_nirspec_saturation_xxxx.fits
**straylight**       jwst_miri_straymask_xxxx.fits
**superbias**        jwst_nircam_superbias_xxxx.fits
..                   jwst_niriss_superbias_xxxx.fits
..                   jwst_nirspec_superbias_xxxx.fits
==================  ========================================

(*) xxx indicates different version numbers for the reference files.

Common Keywords to All Reference Files
--------------------------------------

At present, most JWST science and reference files are FITS files with image or table extensions. 
The FITS primary data unit is always empty. The primary header contains all keywords not specific to individual extensions. Keywords specific to a particular extension are contained in the header of that extension.

The required Keywords Documenting Contents of Reference Files are:

========  =========================================================================
Keyword   Comment
========  =========================================================================
REFTYPE   Required values are listed in the discussion of each pipeline step.
DESCRIP   Summary of file content and/or reason for delivery
AUTHOR    Person(s) who created the file
USEAFTER  YYYY-MM-DDThh:mm:ss Date and time after the reference files will be used. 
          The T is required. Time string may NOT be omitted; use T00:00:00 if no 
          meaningful value is available.
PEDIGREE  Options are
          'SIMULATION'
          'GROUND'
          'DUMMY'
          'INFLIGHT YYYY-MM-DD YYYY-MM-DD'
HISTORY   'Description of Reference File Creation'
HISTORY   DOCUMENT: Name of document describing the strategy and algorithms used to 
          create file
HISTORY   SOFTWARE: Description, version number, location of software used to create 
          file
HISTORY   DATA USED: Data used to create file
HISTORY   DIFFERENCES: How is this version different from the one that it replaces?
HISTORY   If your text spills over to the next line,
HISTORY   begin it with another HISTORY keyword, as in this example.
========  =========================================================================

A pipeline module may require separate reference files for each instrument, detector, 
filter, observation date, etc.  The values of these parameters must be included in the 
reference file header.  The observing-mode keyword values are vital to the process of 
ingesting reference files into CRDS, as they are used to establish the mapping between 
observing modes and specific reference files. Some observing-mode keywords are also 
used in the pipeline processing steps.  If an observing-mode keyword is irrelevant to a 
particular observing mode (such as GRATING for the MIRI imager mode or the NIRCam and NIRISS 
instruments), then it may be omitted from the file header. The Keywords Documenting the Observing 
Mode are:

========  ==================  =============================================================================================
Keyword   Sample Value        Comment
========  ==================  =============================================================================================
TELESCOP  JWST     
INSTRUME  MIRI                Instrument name. Allowed values: FGS, NIRCAM, NIRISS, NIRSPEC, MIRI
PUPIL     NRM                 Pupil wheel element. Required only for NIRCam and NIRISS.
                              NIRCam allowed values: CLEAR, F162M, F164N, F323N, F405N, F466N, F470N, GRISMV2, GRISMV3
                              NIRISS allowed values: CLEARP, F090W, F115W, F140M, F150W, F158M, F200W, GR700XD, NRM
FILTER    F2100W              Filter wheel element. Allowed values: too many to list here
GRATING   G395M               Required only for NIRSpec.

                              NIRSpec allowed values: G140M, G235M, G395M, G140H, G235H, G395H, PRISM, MIRROR
EXP_TYPE  MIR_MRS             Exposure type.

                              FGS allowed values: FGS_IMAGE, FGS_FOCUS, FGS_SKYFLAT, FGS_INTFLAT, FGS_DARK

                              MIRI allowed values: MIR_IMAGE, MIR_TACQ, MIR_LYOT, MIR_4QPM, MIR_LRS-FIXEDSLIT, 
                              MIR_LRS-SLITLESS, MIR_MRS, MIR_DARK, MIR_FLATIMAGE, MIR_FLATMRS, MIR_CORONCAL

                              NIRCam allowed values: NRC_IMAGE, NRC_GRISM, NRC_TACQ, NRC_TACONFIRM, NRC_CORON, 
                              NRC_TSIMAGE, NRC_TSGRISM, NRC_FOCUS, NRC_DARK, NRC_FLAT, NRC_LED

                              NIRISS allowed values: NIS_IMAGE, NIS_TACQ, NIS_TACONFIRM, NIS_WFSS, NIS_SOSS, NIS_AMI, 
                              NIS_FOCUS, NIS_DARK, NIS_LAMP

                              NIRSpec allowed values: NRS_TASLIT, NRS_TACQ, NRS_TACONFIRM, NRS_CONFIRM, NRS_FIXEDSLIT, 
                              NRS_AUTOWAVE, NRS_IFU, NRS_MSASPEC, NRS_AUTOFLAT, NRS_IMAGE, NRS_FOCUS, NRS_DARK, NRS_LAMP, 
                              NRS_BOTA, NRS_BRIGHTOBJ, NRS_MIMF
DETECTOR  MIRIFULONG          Allowed values:
                              GUIDER1, GUIDER2

                              NIS

                              NRCA1, NRCA2, NRCA3, NRCA4, NRCB1, NRCB2, NRCB3, NRCB4, NRCALONG, NRCBLONG

                              NRS1, NRS2

                              MIRIFULONG, MIRIFUSHORT, MIRIMAGE

CHANNEL   12                  MIRI MRS (IFU) channel. Allowed values: 1, 2, 3, 4, 12, 34
                              SHORT   NIRCam channel. Allowed values: SHORT, LONG
BAND      MEDIUM              IFU band. Required only for MIRI. Allowed values are SHORT, MEDIUM, LONG, and N/A, as well 
                              as any allowable combination of two values (SHORT-MEDIUM, LONG-SHORT, etc.). (Also used as 
                              a header keyword for selection of all MIRI Flat files, Imager included.)
READPATT  FAST                Name of the readout pattern used for the exposure. Each pattern represents a particular 
                              combination of parameters like nframes and groups. For MIRI, FAST and SLOW refer to the rate 
                              at which the detector is read.

                              MIRI allowed values: SLOW, FAST, FASTGRPAVG, FASTINTAVG

                              NIRCam allowed values: DEEP8, DEEP2, MEDIUM8, MEDIUM2, SHALLOW4, SHALLOW2, BRIGHT2, BRIGHT1, 
                              RAPID

                              NIRSpec allowed values: NRSRAPID, NRS, NRSN16R4, NRSIRS2RAPID

                              NIRISS allowed values: NIS, NISRAPID

                              FGS allowed values: ID, ACQ1, ACQ2, TRACK, FINEGUIDE, FGS60, FGS840, FGS7850, FGSRAPID, FGS
NRS_NORM  16                  Required only for NIRSpec.
NRS_REF   4                   Required only for NIRSpec.
SUBARRAY  FULL                MIRI allowed values: FULL, GENERIC, MASK1140, MASK1550, MASK1065, MASKLYOT, BRIGHTSKY, SUB256, 
                              SUB128, SUB64, SLITLESSPRISM
P_XXXXXX  P_READPA            pattern keywords used by CRDS for JWST to describe the intended uses of a reference file 
                              using or'ed combinations of values. Only a subset of :ref:`p-patterns` 
                              are supported.
SUBSTRT1  1                   Starting pixel index along axis 1 (1-indexed)
SUBSIZE1  2048                Size of subarray along axis 1
SUBSTRT2  1                   Starting pixel index along axis 2 (1-indexed)
SUBSIZE2  2048                Size of subarray along axis 2
FASTAXIS  1                   Fast readout direction relative to image axes for Amplifier #1 (1 = +x axis, 2 = +y axis,
                              -1 = -x axis, -2 = -y axis) SEE NOTE BELOW.
SLOWAXIS  2                   Slow readout direction relative to image axes for all amplifiers (1 = +x axis, 2 = +y axis,
                               -1 = -x axis, -2 = -y axis)
========  ==================  =============================================================================================

Note: For the NIR detectors, the fast readout direction changes sign from one amplifier to the next.  It is +1, -1, +1, and -1, for amps 1, 2, 3, and 4, respectively.  The keyword FASTAXIS refers specifically to amp 1.  That way, it is entirely correct for single-amp readouts and correct at the origin for 4-amp readouts.  For MIRI, FASTAXIS is always +1.


Tracking Pipeline Progress
++++++++++++++++++++++++++

As each pipeline step is applied to a science data product, it will record a status indicator in a header keyword of the science data product. The current list of step status keyword names is given in the following table. These status keywords may be included in the primary header of reference files, in order to maintain a history of the data that went into creating the reference file. Allowed values for the status keywords are 'COMPLETE' and 'SKIPPED'. Absence of a particular keyword is understood to mean that step was not even attempted.

Table 3.  Keywords Documenting Which Pipeline Steps Have Been Performed

=========   ========================================
S_IPC       IPC correction  
S_RESET     MIRI reset correction
S_SUPERB    Superbias subtraction   
S_IMPRNT    NIRSpec MSA imprint subtraction
S_MSAFLG    NIRSpec MSA failed shutter flagging 
S_EXTR1D    1-D spectral extraction
S_LASTFR    MIRI last frame correction  
S_DQINIT    DQ initialization
S_REFPIX    Reference pixel correction  
S_ERRINI    ERR initialization
S_DARK      Dark subtraction    
S_SATURA    Saturation check
S_LINEAR    Linearity correction    
S_JUMP      Jump detection
S_RAMP      Ramp fitting    
S_WCS       WCS assignment
S_FLAT      Flat-fielding   
S_FRINGE    Fringe correction
S_PERSIS    Persistence correction  
S_STRAY     Straylight correction
S_TELEMI    Telescope emission  
S_PHOTOM    Photometric (absolute flux) calibration
S_EXTR1D    1-D extraction  
S_EXTR2D    2-D spectral extraction
S_RESAMP    Image resampling    
S_BKDSUB    Background subtraction
S_SLOSS     Slit-loss correction         
=========   ========================================



Orientation of Detector Image
+++++++++++++++++++++++++++++

All steps in the pipeline assume the data are in the DMS (science) orientation, not the native readout orientation. The pipeline does NOT check or correct for the orientation of the reference data. It assumes that all files ingested into CRDS have been put into the science orientation.  All header keywords documenting the observing mode (Table 2) should likewise be transformed into the DMS orientation.   For square data array dimensions it's not possible to infer the actual orientation directly so reference file authors must manage orientation carefully.   

    Correct values for FASTAXIS and SLOWAXIS for each detector are:
=========== ======== ========
DETECTOR    FASTAXIS SLOWAXIS
=========== ======== ========
MIRIMAGE      1       2
MIRIFULONG    1       2
MIRIFUSHORT   1       2
NRCA1        -1       2
NRCA2         1      -2
NRCA3        -1       2
NRCA4         1      -2
NRCALONG     -1       2
NRCB1         1      -2
NRCB2        -1       2
NRCB3         1      -2
NRCB4        -1       2
NRCBLONG      1      -2
NRS1          2       1
NRS2         -2      -1
NIS          -2      -1
GUIDER1      -2      -1
GUIDER2       2      -1
=========== ======== ========

Differing values for these keywords will be taken as an indicator that neither the keyword value nor the array orientation are correct.

.. _p-patterns:

P_pattern keywords
------------------
P_ pattern keywords used by CRDS for JWST to describe the intended uses of a reference file using or’ed combinations

For example, if the same NIRISS SUPERBIAS should be used for

    READPATT=NIS

or

    READPATT=NISRAPID

the definition of READPATT in the calibration s/w datamodels schema does not allow it. READPATT can specify one or the other but not both.

To support expressing combinations of values, CRDS and the CAL s/w have added “pattern keywords” which nominally begin with P_ followed by the ordinary keyword, truncated as needed to 8 characters. In this case, P_READPA corresponds to READPATT.

Pattern keywords override the corresponding ordinary keyword for the purposes of automatically updating CRDS rmaps. Pattern keywords describe intended use.

In this example, the pattern keyword:

    P_READPA = ‘NIS | NISRAPID |‘

can be used to specify the intent “use for NIS or for NISRAPID”.

Only or-ed combinations of the values used in ordinary keywords are valid for pattern keywords.

Patterns appear in a slightly different form in rmaps than they do in P_ keywords. The value of a P_ keyword always ends with a trailing or-bar. In rmaps, no trailing or-bar is used so the equivalient of the above in an rmap is:

    ‘NIS|NISRAPID’
    
    From a CRDS perspective, the ``P_ pattern`` keywords and their corresponding datamodels paths currently supported can be found in the 
    `JWST Pattern Keywords section of the CRDS documentation. <https://jwst-crds.stsci.edu/static/users_guide/reference_conventions.html#id2>`_ 

Currently all ``P_`` keywords correspond to basic keywords found only in the primary headers of reference files and are typically only valid for FITS format..

The traslation from these ``P_pattern`` keywords are completely generic in CRDS and can apply to any reference file type so they should be assumed to 
be reserved whether a particular type uses them or not. Defining non-pattern keywords with the prefix ``P_`` is strongly discouraged.

Data Quality Flags
==================

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

Note: Words like "highly" and "large" will be defined by each instrument team.  They are likely to vary from one detector to another – or even from one observing mode to another.  
