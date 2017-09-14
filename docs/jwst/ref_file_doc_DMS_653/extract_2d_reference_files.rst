Reference File Types
--------------------

The extract_2d step uses the following reference files:

===================    ==========================================================   ============================
reftype                                     description                              Instrument
===================    ==========================================================   ============================
**wavecorr**           NIRSPEC wavelength zero-point correction                      NIRSPEC
===================    ==========================================================   ============================


CRDS Selection Criteria For Each Reference File Type 
----------------------------------------------------

WAVECORR
::::::::
The WAVECORR reference file is selected based on EXP_TYPE of the science data.
The reference file s relevant only for Nirspec observations with EXP_TYPE of
NRS_FIXEDSLIT, NRS_MSASPEC, NRS_BRIGHTOBJ.


Reference File Formats For Each Reference File Type 
---------------------------------------------------

WAVECORR
::::::::
The WAVECORR file contains reference data about the NIRSPEC wavelength zero-point correction.
The reference data is described in the NIRSPEC Technical Note ESA-JWST–SCI-NRS-TN-2016-018.

:apertures:
   :aperture_name: Aperture name.
    :variance: Estimated variance on the zero-point offset.
    :width: Aperture width [SLIT] or pitch [MOS].
    :zero_point_offset: Zero-point offset as a function of wavelength (in m)
                        and source offset within the aperture (in units of fraction of the aperture width
                        [SLIT] or pitch [MOS]).
    

