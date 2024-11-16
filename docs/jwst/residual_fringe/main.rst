Description
===========

:Class: `jwst.residual_fringe.ResidualFringeStep`
:Alias: residual_fringe

The JWST pipeline contains multiple steps to mitigate the impact of fringing on science spectra, which
generally suffice to reduce the fringe signal to below a few percent of the target flux.

The first correction is applied by default in the :ref:`fringe <fringe_step>` step in the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline
and  consists of dividing the uncalibrated "rate" image by a static fringe flat constructed from observations of a
bright source that fills the entire MRS field of view. For more details see the :ref:`fringe <fringe_step>` step.
This  step generally does a good job of removing the strongest fringes from an astronomical scene, particularly
for nearly-uniform extended sources. Since the fringe signal is different for point sources, however, and varies
as a function of the location of a point source within the FOV, the static fringe flat cannot fully correct
such objects and the default high level data products will therefore still
show appreciable fringes.

The pipeline also includes two optional residual fringe correction steps whose purpose is to find and remove signals
whose periodicity is consistent with known fringe frequencies (set by the optical thickness of the detectors
and dichroics) using a Lomb-Scargle periodogram. The number of fringe components to be removed is governed
by a Bayesian evidence calculation.
The first of these residual fringe correction steps is a 2-D correction that can be applied to the flux-calibrated detector data
in the :ref:`residual_fringe <residual_fringe_step>` step. This step is part of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, but currently
it is skipped by default. To apply this step set the step parameter,  ``--skip = False``. This step is applied after
:ref:`photom <photom_step>`, but before :ref:`cube_build <cube_build_step>`.


The ``residual_fringe`` step can accept several different forms of input data, including:

#. a single file containing a 2-D IFU image

#. a data model (`~jwst.datamodels.IFUImageModel`) containing a 2-D IFU image

#. an association table (in json format) containing a single input file
   
The second of the residual fringe correction steps is a 1-D correction  that can be applied to one-dimensional
spectra extracted from MRS data cubes by setting the optional parameter ``extract_1d.ifu_rfcorr = True``
in the :ref:`extract_1d <extract_1d_step>` step.  Empirically, the 1-D correction step has been found to work
better than the 2-D correction step if it is applied to per-band spectra.
For more details on this step see :ref:`extract_1d <extract_1d_step>` step. 


Assumptions
-----------
This step only works on MIRI MRS data.


Fringe Background Information
-----------------------------
As is typical for spectrometers, the MIRI MRS detectors are affected by fringes.  These are periodic gain modulations caused by
standing waves between parallel surfaces in the optical path, acting as a slow-finesse Fabry-Pérot etalons. In the MRS,
the principal fringe sources are the detector layers. A detailed  detailed discussion on these fringe components
can be found in Argyriou, I., Wells, M., Glasse, A., et al. 2020, A&A, 641, A150 and
Wells, M., Pel, J.-W., Glasse, A., et al. 2015, PASP, 127, 646.


The primary MRS fringe, observed in all MRS bands, is caused by the etalons between the anti-reflection coating
and lower layers, encompassing the detector substrate and the infrared-active layer. Since the thickness of the substrate
is not the same in the SW and LW detectors, the fringe frequency will differ in the two detectors. Up to 16 microns, this
fringe is produced by the anti-reflection coating and  pixel metalization etalons, whereas above 16 microns it is
produced by the anti-reflection coating and  bottom contact etalon, resulting in a different fringe frequency.
The information in the fringe frequency
reference file  is used to determine, for each MRS band, the frequencies to fit to this main fringe component.
The residual fringes are corrected for by fitting and removing sinusoidal gain to the detector level data.
While the fringe frequencies are well known, amplitudes can vary due to beating between the different fringe components
and additionally are sensitive to the detailed location and intensity of objects within a given astronomical scene.
Fringing thus cannot be corrected in its entirety for an arbitrary astronomical scene without forward modeling.



