Description
===========

:Class: `jwst.residual_fringe.ResidualFringeStep`
:Alias: residual_fringe

The JWST pipeline contains two steps devoted to the removal of fringes on MIRI MRS images.
The first correction is applied in the ``fringe_step`` in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline
and  consists in dividing
detector-level data by a fringe-flat and is described in the :ref:`fringe <fringe_step>` step.
Applying the fringe flat should eliminate fringes from spectra of spatially extended sources, however
residual fringes can remain. For spatially unresolved (point) sources or extended sources with structure,
applying the fringe flat will undoubtedly leave residual fringes since these produce different fringe patterns
on the detector than accounted for by the fringe flat. The second step for fringe removal is the
``residual_fringe_step``. This step is part of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, but currently
it is skipped by default. To apply this step set the step parameter,  ``--skip = False``. This step is applied after
:ref:`photom <photom_step>`, but before :ref:`cube_build <cube_build_step>`.



The ``residual_fringe`` step can accept several different forms of input data, including:

  - a single file containing a 2-D IFU image

  - a data model (IFUImageModel) containing a 2-D IFU image

  - an association table (in json format) containing a single input file


Assumptions
-----------
This step only works on MIRI MRS data.
It is assumed that the calwebb_spec2 pipeline has been run on the data. In addition, the detection of residual fringes
are  better determined if the ``mrs_imatch``  step has also been applied to the data.



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

