Introduction
============

This document is intended to be a core reference guide to the formats, naming convention and
data quality flags used by the reference files for pipeline steps requiring them, and is not
intended to be a detailed description of each of those pipeline steps. It also does not give
details on pipeline steps that do not use reference files.
The present manual is referred to by several other documentation pages,
such as the JWST pipeline and JDocs.

Reference File Naming Convention
================================

Before reference files are ingested into CRDS, they are renamed following a
convention used by the pipeline. As with any other changes undergone by the reference files,
the previous names are kept in header keywords, so the Instrument Teams
can easily track which delivered file is being used by the pipeline in each step.

The naming of reference files uses the following syntax::

 jwst_<instrument>_<reftype>_<version>.<extension>

where

- ``instrument`` is one of "fgs", "miri", "nircam", "niriss", and "nirspec"
- ``reftype`` is one of the type names listed in the table below
- ``version`` is a 4-digit version number (e.g. 0042)
- ``extension`` gives the file format, such as "fits" or "asdf"

An example NIRCam GAIN reference file name would be "jwst_nircam_gain_0042.fits".

The HISTORY header keyword of each reference file includes details on specific processing
undergone by the files before being ingested in CRDS.

.. _reference_file_types:

Reference File Types
====================

Most reference files have a one-to-one relationship with calibration steps, e.g.
there is one step that uses one type of reference file. Some steps, however, use
several types of reference files and some reference file types are used by more
than one step. The tables below show the correspondence between pipeline steps and
reference file types. The first table is ordered by pipeline step, while the second
is ordered by reference file type. Links to the reference file types provide detailed
documentation on each reference file.

+-----------------------------------------------+--------------------------------------------------+
| Pipeline Step                                 | Reference File Type (REFTYPE)                    |
+===============================================+==================================================+
| :ref:`align_refs <align_refs_step>`           | :ref:`PSFMASK <psfmask_reffile>`                 |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`ami_analyze <ami_analyze_step>`         | :ref:`THROUGHPUT <throughput_reffile>`           |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`assign_wcs <assign_wcs_step>`           | :ref:`CAMERA <camera_reffile>`                   |
+                                               +--------------------------------------------------+
|                                               | :ref:`COLLIMATOR <collimator_reffile>`           |
+                                               +--------------------------------------------------+
|                                               | :ref:`DISPERSER <disperser_reffile>`             |
+                                               +--------------------------------------------------+
|                                               | :ref:`DISTORTION <distortion_reffile>`           |
+                                               +--------------------------------------------------+
|                                               | :ref:`FILTEROFFSET <filteroffset_reffile>`       |
+                                               +--------------------------------------------------+
|                                               | :ref:`FORE <fore_reffile>`                       |
+                                               +--------------------------------------------------+
|                                               | :ref:`FPA <fpa_reffile>`                         |
+                                               +--------------------------------------------------+
|                                               | :ref:`IFUFORE <ifufore_reffile>`                 |
+                                               +--------------------------------------------------+
|                                               | :ref:`IFUPOST <ifupost_reffile>`                 |
+                                               +--------------------------------------------------+
|                                               | :ref:`IFUSLICER <ifuslicer_reffile>`             |
+                                               +--------------------------------------------------+
|                                               | :ref:`MSA <msa_reffile>`                         |
+                                               +--------------------------------------------------+
|                                               | :ref:`OTE <ote_reffile>`                         |
+                                               +--------------------------------------------------+
|                                               | :ref:`SPECWCS <specwcs_reffile>`                 |
+                                               +--------------------------------------------------+
|                                               | :ref:`REGIONS <regions_reffile>`                 |
+                                               +--------------------------------------------------+
|                                               | :ref:`WAVELENGTHRANGE <wavelengthrange_reffile>` |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`background <background_step>`           | :ref:`WFSSBKG <wfssbkg_reffile>`                 |
+                                               +--------------------------------------------------+
|                                               | :ref:`WAVELENGTHRANGE <wavelengthrange_reffile>` |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`barshadow <barshadow_step>`             | :ref:`BARSHADOW <barshadow_reffile>`             |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`cube_build <cube_build_step>`           | :ref:`CUBEPAR <cubepar_reffile>`                 |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`dark_current <dark_current_step>`       | :ref:`DARK <dark_reffile>`                       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`dq_init <dq_init_step>`                 | :ref:`MASK <mask_reffile>`                       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`extract_1d <extract_1d_step>`           | :ref:`EXTRACT1D <extract1d_reffile>`             |
+                                               +--------------------------------------------------+
|                                               | :ref:`APCORR <apcorr_reffile>`                   |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`extract_2d <extract_2d_step>`           | :ref:`WAVECORR <wavecorr_reffile>`               |
+                                               +--------------------------------------------------+
|                                               | :ref:`WAVELENGTHRANGE <wavelengthrange_reffile>` |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`flatfield <flatfield_step>`             | :ref:`FLAT <flat_reffile>`                       |
+                                               +--------------------------------------------------+
|                                               | :ref:`DFLAT <dflat_reffile>`                     |
+                                               +--------------------------------------------------+
|                                               | :ref:`FFLAT <fflat_reffile>`                     |
+                                               +--------------------------------------------------+
|                                               | :ref:`SFLAT <sflat_reffile>`                     |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`fringe <fringe_step>`                   | :ref:`FRINGE <fringe_reffile>`                   |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`gain_scale <gain_scale_step>`           | :ref:`GAIN <gain_reffile>`                       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`ipc <ipc_step>`                         | :ref:`IPC <ipc_reffile>`                         |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`jump <jump_step>`                       | :ref:`GAIN <gain_reffile>`                       |
+                                               +--------------------------------------------------+
|                                               | :ref:`READNOISE <readnoise_reffile>`             |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`linearity <linearity_step>`             | :ref:`LINEARITY <linearity_reffile>`             |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`msaflagopen <msaflagopen_step>`         | :ref:`MSAOPER <msaoper_reffile>`                 |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`pathloss <pathloss_step>`               | :ref:`PATHLOSS <pathloss_reffile>`               |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`persistence <persistence_step>`         | :ref:`PERSAT <persat_reffile>`                   |
+                                               +--------------------------------------------------+
|                                               | :ref:`TRAPDENSITY <trapdensity_reffile>`         |
+                                               +--------------------------------------------------+
|                                               | :ref:`TRAPPARS <trappars_reffile>`               |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`photom <photom_step>`                   | :ref:`PHOTOM <photom_reffile>`                   |
+                                               +--------------------------------------------------+
|                                               | :ref:`AREA <area_reffile>`                       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`ramp_fitting <ramp_fitting_step>`       | :ref:`GAIN <gain_reffile>`                       |
+                                               +--------------------------------------------------+
|                                               | :ref:`READNOISE <readnoise_reffile>`             |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`refpix <refpix_step>`                   | :ref:`REFPIX <refpix_reffile>`                   |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`resample <resample_step>`               | :ref:`DRIZPARS <drizpars_reffile>`               |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`reset <reset_step>`                     | :ref:`RESET <reset_reffile>`                     |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`residual_fringe <residual_fringe_step>` | :ref:`FRINGEFREQ <fringefreq_reffile>`           |
+                                               +--------------------------------------------------+
|                                               | :ref:`REGIONS <regions_reffile>`                 |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`rscd <rscd_step>`                       | :ref:`RSCD <rscd_reffile>`                       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`saturation <saturation_step>`           | :ref:`SATURATION <saturation_reffile>`           |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`source_catalog <source_catalog_step>`   | :ref:`APCORR <apcorr_reffile>`                   |
+                                               +--------------------------------------------------+
|                                               | :ref:`ABVEGAOFFSET <abvegaoffset_reffile>`       |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`straylight <straylight_step>`           | :ref:`MRSXARTCORR <mrsxartcorr_reffile>`         |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`superbias <superbias_step>`             | :ref:`SUPERBIAS <superbias_reffile>`             |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`tso_photometry <tso_photometry_step>`   | :ref:`TSOPHOT <tsophot_reffile>`                 |
+-----------------------------------------------+--------------------------------------------------+
| :ref:`wavecorr <wavecorr_step>`               | :ref:`WAVECORR <wavecorr_reffile>`               |
+-----------------------------------------------+--------------------------------------------------+

+--------------------------------------------------+-----------------------------------------------+
| Reference File Type (REFTYPE)                    | Pipeline Step                                 |
+==================================================+===============================================+
| :ref:`ABVEGAOFFSET <abvegaoffset_reffile>`       | :ref:`source_catalog <source_catalog_step>`   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`APCORR <apcorr_reffile>`                   | :ref:`extract_1d <extract_1d_step>`           |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`source_catalog <source_catalog_step>`   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`AREA <area_reffile>`                       | :ref:`photom <photom_step>`                   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`BARSHADOW <barshadow_reffile>`             | :ref:`barshadow <barshadow_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`CAMERA <camera_reffile>`                   | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`COLLIMATOR <collimator_reffile>`           | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`CUBEPAR <cubepar_reffile>`                 | :ref:`cube_build <cube_build_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`DARK <dark_reffile>`                       | :ref:`dark_current <dark_current_step>`       |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`DFLAT <dflat_reffile>`                     | :ref:`flatfield <flatfield_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`DISPERSER <disperser_reffile>`             | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`DISTORTION <distortion_reffile>`           | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`DRIZPARS <drizpars_reffile>`               | :ref:`resample <resample_step>`               |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`EXTRACT1D <extract1d_reffile>`             | :ref:`extract_1d <extract_1d_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FFLAT <fflat_reffile>`                     | :ref:`flatfield <flatfield_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FILTEROFFSET <filteroffset_reffile>`       | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FLAT <flat_reffile>`                       | :ref:`flatfield <flatfield_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FORE <fore_reffile>`                       | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FPA <fpa_reffile>`                         | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FRINGE <fringe_reffile>`                   | :ref:`fringe <fringe_step>`                   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`FRINGEFREQ <fringefreq_reffile>`           | :ref:`residual_fringe <residual_fringe_step>` |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`GAIN <gain_reffile>`                       | :ref:`gain_scale <gain_scale_step>`           |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`jump <jump_step>`                       |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`ramp_fitting <ramp_fitting_step>`       |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`IFUFORE <ifufore_reffile>`                 | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`IFUPOST <ifupost_reffile>`                 | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`IFUSLICER <ifuslicer_reffile>`             | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`IPC <ipc_reffile>`                         | :ref:`ipc <ipc_step>`                         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`LINEARITY <linearity_reffile>`             | :ref:`linearity <linearity_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`MASK <mask_reffile>`                       | :ref:`dq_init <dq_init_step>`                 |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`MRSXARTCORR <mrsxartcorr_reffile>`         | :ref:`straylight <straylight_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`MSA <msa_reffile>`                         | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`MSAOPER <msaoper_reffile>`                 | :ref:`msaflagopen <msaflagopen_step>`         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`OTE <ote_reffile>`                         | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`PATHLOSS <pathloss_reffile>`               | :ref:`pathloss <pathloss_step>`               |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`PERSAT <persat_reffile>`                   | :ref:`persistence <persistence_step>`         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`PHOTOM <photom_reffile>`                   | :ref:`photom <photom_step>`                   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`PSFMASK <psfmask_reffile>`                 | :ref:`align_refs <align_refs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`READNOISE <readnoise_reffile>`             | :ref:`jump <jump_step>`                       |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`ramp_fitting <ramp_fitting_step>`       |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`REFPIX <refpix_reffile>`                   | :ref:`refpix <refpix_step>`                   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`REGIONS <regions_reffile>`                 | :ref:`assign_wcs <assign_wcs_step>`           |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`residual_fringe <residual_fringe_step>` |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`RESET <reset_reffile>`                     | :ref:`reset <reset_step>`                     |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`RSCD <rscd_reffile>`                       | :ref:`rscd <rscd_step>`                       |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`SATURATION <saturation_reffile>`           | :ref:`saturation <saturation_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`SFLAT <sflat_reffile>`                     | :ref:`flatfield <flatfield_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`SPECWCS <specwcs_reffile>`                 | :ref:`assign_wcs <assign_wcs_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`SUPERBIAS <superbias_reffile>`             | :ref:`superbias <superbias_step>`             |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`THROUGHPUT <throughput_reffile>`           | :ref:`ami_analyze <ami_analyze_step>`         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`TRAPDENSITY <trapdensity_reffile>`         | :ref:`persistence <persistence_step>`         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`TRAPPARS <trappars_reffile>`               | :ref:`persistence <persistence_step>`         |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`TSOPHOT <tsophot_reffile>`                 | :ref:`tso_photometry <tso_photometry_step>`   |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`WAVELENGTHRANGE <wavelengthrange_reffile>` | :ref:`assign_wcs <assign_wcs_step>`           |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`background <background_step>`           |
+                                                  +-----------------------------------------------+
|                                                  | :ref:`extract_2d <extract_2d_step>`           |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`WAVECORR <wavecorr_reffile>`               | :ref:`wavecorr <wavecorr_step>`               |
+--------------------------------------------------+-----------------------------------------------+
| :ref:`WFSSBKG <wfssbkg_reffile>`                 | :ref:`background <background_step>`           |
+--------------------------------------------------+-----------------------------------------------+

Step Parameters Reference Types
+++++++++++++++++++++++++++++++

When each ``Step`` is instantiated, a CRDS look-up, based on the ``Step`` class
name and input data, is made to retrieve a parameter file. The ``reftype``
for such parameter files is ``pars-<class name>``. For example, for the step
``jwst.persistence.PersistenceStep``, the ``reftype`` would be
``pars-persistencestep``.

For more information, see :ref:`parameter_files`.

.. _`Standard Required Keywords`:

Standard Required Keywords
==========================

At present, most JWST science and reference files are FITS files with image or table extensions.
The FITS primary data unit is always empty. The primary header contains all keywords not specific to individual extensions. Keywords specific to a particular extension are contained in the header of that extension.

The required Keywords Documenting Contents of Reference Files are:

========  ==================================================================================
Keyword   Comment
========  ==================================================================================
REFTYPE   `WFSSBKG    Required values are listed in the discussion of each pipeline step.`
DESCRIP   `Summary of file content and/or reason for delivery`
AUTHOR    `Fred Jones     Person(s) who created the file`
USEAFTER  `YYYY-MM-DDThh:mm:ss Date and time after the reference files will
          be used. The T is required. Time string may NOT be omitted;
          use T00:00:00 if no meaningful value is available.`
PEDIGREE  `Options are
          'SIMULATION'
          'GROUND'
          'DUMMY'
          'INFLIGHT YYYY-MM-DD YYYY-MM-DD'`
HISTORY   `Description of Reference File Creation`
HISTORY   `DOCUMENT: Name of document describing the strategy and algorithms
          used to create file.`
HISTORY   `SOFTWARE: Description, version number, location of software used
          to create file.`
HISTORY   `DATA USED: Data used to create file`
HISTORY   `DIFFERENCES: How is this version different from the one that
          it replaces?`
HISTORY   `If your text spills over to the next line,
          begin it with another HISTORY keyword, as in this example.`
TELESCOP  `JWST   Name of the telescope/project.`
INSTRUME  `FGS   Instrument name. Allowed values: FGS, NIRCAM, NIRISS,
          NIRSPEC, MIRI`
SUBARRAY  `FULL, GENERIC, SUBS200A1, ...   (XXX abstract technical description
          of SUBARRAY)`
SUBSTRT1  `1        Starting pixel index along axis 1 (1-indexed)`
SUBSIZE1  `2048     Size of subarray along axis 1`
SUBSTRT2  `1        Starting pixel index along axis 2 (1-indexed)`
SUBSIZE2  `2048     Size of subarray along axis 2`
FASTAXIS  `1        Fast readout direction relative to image axes for
          Amplifier #1 (1 = +x axis, 2 = +y axis, -1 = -x axis, -2 = -y axis)
          SEE NOTE BELOW.`
SLOWAXIS  `2        Slow readout direction relative to image axes for
          all amplifiers (1 = +x axis, 2 = +y axis, -1 = -x axis, -2 = -y axis)`
========  ==================================================================================


Observing Mode Keywords
=======================

A pipeline module may require separate reference files for each instrument, detector,
filter, observation date, etc.  The values of these parameters must be included in the
reference file header.  The observing-mode keyword values are vital to the process of
ingesting reference files into CRDS, as they are used to establish the mapping between
observing modes and specific reference files. Some observing-mode keywords are also
used in the pipeline processing steps.  If an observing-mode keyword is irrelevant to a
particular observing mode (such as GRATING for the MIRI imager mode or the NIRCam and NIRISS
instruments), then it may be omitted from the file header.

The Keywords Documenting the Observing Mode are:

========  ==================  =============================================================================================
Keyword   Sample Value        Comment
========  ==================  =============================================================================================
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
P_XXXXXX  P_READPA            pattern keywords used by CRDS for JWST to describe the intended uses of a reference file
                              using or'ed combinations of values. Only a subset of :ref:`p-patterns`
                              are supported.
========  ==================  =============================================================================================

Note: For the NIR detectors, the fast readout direction changes sign from one amplifier to the next.  It is +1, -1, +1, and -1, for amps 1, 2, 3, and 4, respectively.  The keyword FASTAXIS refers specifically to amp 1.  That way, it is entirely correct for single-amp readouts and correct at the origin for 4-amp readouts.  For MIRI, FASTAXIS is always +1.


Tracking Pipeline Progress
++++++++++++++++++++++++++

As each pipeline step is applied to a science data product, it will record a status indicator in a
header keyword of the science data product. The current list of step status keyword names is given
in the following table. These status keywords may be included in the primary header of reference
files, in order to maintain a history of the data that went into creating the reference file.
Allowed values for the status keywords are 'COMPLETE' and 'SKIPPED'. Absence of a particular keyword
is understood to mean that step was not even attempted.

Table 1.  Keywords Documenting Which Pipeline Steps Have Been Performed.

=========   ========================================
S_AMIANA    AMI fringe analysis
S_AMIAVG    AMI fringe averaging
S_AMINOR    AMI fringe normalization
S_BARSHA    Bar shadow correction
S_BKDSUB    Background subtraction
S_COMB1D    1-D spectral combination
S_DARK      Dark subtraction
S_DQINIT    DQ initialization
S_ERRINI    ERR initialization
S_EXTR1D    1-D spectral extraction
S_EXTR2D    2-D spectral extraction
S_FLAT      Flat field correction
S_FRINGE    Fringe correction
S_FRSTFR    MIRI first frame correction
S_GANSCL    Gain scale correction
S_GRPSCL    Group scale correction
S_GUICDS    Guide mode CDS computation
S_IFUCUB    IFU cube creation
S_IMPRNT    NIRSpec MSA imprint subtraction
S_IPC       IPC correction
S_JUMP      Jump detection
S_KLIP      Coronagraphic PSF subtraction
S_LASTFR    MIRI last frame correction
S_LINEAR    Linearity correction
S_MRSMAT    MIRI MRS background matching
S_MSAFLG    NIRSpec MSA failed shutter flagging
S_OUTLIR    Outlier detection
S_PERSIS    Persistence correction
S_PHOTOM    Photometric (absolute flux) calibration
S_PSFALI    Coronagraphic PSF alignment
S_PSFSTK    Coronagraphic PSF stacking
S_PTHLOS    Pathloss correction
S_RAMP      Ramp fitting
S_REFPIX    Reference pixel correction
S_RESAMP    Resampling (drizzling)
S_RESET     MIRI reset correction
S_RSCD      MIRI RSCD correction
S_SATURA    Saturation check
S_SKYMAT    Sky matching
S_SRCCAT    Source catalog creation
S_SRCTYP    Source type determination
S_STRAY     Straylight correction
S_SUPERB    Superbias subtraction
S_TELEMI    Telescope emission correction
S_TSPHOT    TSO imaging photometry
S_TWKREG    Tweakreg image alignment
S_WCS       WCS assignment
S_WFSCOM    Wavefront sensing image combination
S_WHTLIT    TSO white-light curve generation
=========   ========================================


Orientation of Detector Image
+++++++++++++++++++++++++++++

All steps in the pipeline assume the data are in the DMS (science) orientation, not the native readout orientation. The pipeline does NOT check or correct for the orientation of the reference data. It assumes that all files ingested into CRDS have been put into the science orientation.  All header keywords documenting the observing mode (Table 2) should likewise be transformed into the DMS orientation.   For square data array dimensions it's not possible to infer the actual orientation directly so reference file authors must manage orientation carefully.

Table 2.  Correct values for FASTAXIS and SLOWAXIS for each detector.

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

``P_`` pattern keywords used by CRDS for JWST to describe the intended uses of a reference file using or’ed combinations

For example, if the same NIRISS SUPERBIAS should be used for

    READPATT=NIS

or

    READPATT=NISRAPID

the definition of READPATT in the calibration s/w datamodels schema does not allow it. READPATT can specify one or the other but not both.

To support expressing combinations of values, CRDS and the CAL s/w have added “pattern keywords” which nominally begin with ``P_`` followed by the ordinary keyword, truncated as needed to 8 characters. In this case, P_READPA corresponds to READPATT.

Pattern keywords override the corresponding ordinary keyword for the purposes of automatically updating CRDS rmaps. Pattern keywords describe intended use.

In this example, the pattern keyword:

  P_READPA = NIS | NISRAPID |

can be used to specify the intent “use for NIS or for NISRAPID”.

Only or-ed combinations of the values used in ordinary keywords are valid for pattern keywords.

Patterns appear in a slightly different form in rmaps than they do in ``P_`` keywords. The value of a ``P_ keyword`` always ends with a trailing or-bar. In rmaps, no trailing or-bar is used so the equivalent of the above in an rmap is:

    ‘NIS|NISRAPID’

    From a CRDS perspective, the ``P_ pattern`` keywords and their corresponding datamodels paths currently supported can be found in the
    `JWST Pattern Keywords section of the CRDS documentation. <https://jwst-crds.stsci.edu/static/users_guide/reference_conventions.html#id2>`_

Currently all ``P_`` keywords correspond to basic keywords found only in the primary headers of reference files and are typically only valid for FITS format..

The translation from these ``P_`` pattern keywords are completely generic in CRDS and can apply to any reference file type so they should be assumed to
be reserved whether a particular type uses them or not. Defining non-pattern keywords with the prefix ``P_`` is strongly discouraged.

.. _`Data Quality Flags`:

Data Quality Flags
==================

Within science data files, the PIXELDQ flags are stored as 32-bit integers;
the GROUPDQ flags are 8-bit integers.  The meaning of each bit is specified
in a separate binary table extension called DQ_DEF.  The binary table has the
format presented in Table 3, which represents the master list of DQ flags.
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


Table 3. Flags for the DQ, PIXELDQ, and GROUPDQ Arrays (Format of DQ_DEF Extension).

===  ==========    ================  ===========================================
Bit  Value         Name              Description
===  ==========    ================  ===========================================
0    1             DO_NOT_USE        Bad pixel. Do not use.
1    2             SATURATED         Pixel saturated during exposure
2    4             JUMP_DET          Jump detected during exposure
3    8             DROPOUT           Data lost in transmission
4    16            OUTLIER           Flagged by outlier detection
5    32            PERSISTENCE       High persistence
6    64            AD_FLOOR          Below A/D floor
7    128           CHARGELOSS        Charge Migration
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
28   268435456     FLUX_ESTIMATED    Pixel flux estimated due to missing/bad data
29   536870912     MSA_FAILED_OPEN   Pixel sees light from failed-open shutter
30   1073741824    OTHER_BAD_PIXEL   A catch-all flag
31   2147483648    REFERENCE_PIXEL   Pixel is a reference pixel
===  ==========    ================  ===========================================

Note: Words like "highly" and "large" will be defined by each instrument team.  They are likely to vary from one detector to another – or even from one observing mode to another.

.. _`dq_parameter_specification`:

Parameter Specification
=======================

There are a number of steps, such as :ref:`OutlierDetectionStep
<outlier_detection_step>` or :ref:`SkyMatchStep <skymatch_step>`, that define
what data quality flags a pixel is allowed to have to be considered in
calculations. Such parameters can be set in a number of ways.

First, the flag can be defined as the integer sum of all the DQ bit values from
the input images DQ arrays that should be considered "good". For example, if
pixels in the DQ array can have combinations of 1, 2, 4, and 8 and one wants to
consider DQ flags 2 and 4 as being acceptable for computations, then the
parameter value should be set to "6" (2+4). In this case a pixel having DQ values
2, 4, or 6 will be considered a good pixel, while a pixel with a DQ value, e.g.,
1+2=3, 4+8="12", etc. will be flagged as a "bad" pixel.

Alternatively, one can enter a comma-separated or '+' separated list of integer
bit flags that should be summed to obtain the final "good" bits. For example,
both "4,8" and "4+8" are equivalent to a setting of "12".

Finally, instead of integers, the JWST mnemonics, as defined above, may be used.
For example, all the following specifications are equivalent:

`"12" == "4+8" == "4, 8" == "JUMP_DET, DROPOUT"`

.. note::
   - The default value (0) will make *all* non-zero
     pixels in the DQ mask be considered "bad" pixels and the
     corresponding pixels will not be used in computations.

   - Setting to `None` will turn off the use of the DQ array
     for computations.

   - In order to reverse the meaning of the flags
     from indicating values of the "good" DQ flags
     to indicating the "bad" DQ flags, prepend '~' to the string
     value. For example, in order to exclude pixels with
     DQ flags 4 and 8 for computations and to consider
     as "good" all other pixels (regardless of their DQ flag),
     use a value of ``~4+8``, or ``~4,8``. A string value of
     ``~0`` would be equivalent to a setting of ``None``.
