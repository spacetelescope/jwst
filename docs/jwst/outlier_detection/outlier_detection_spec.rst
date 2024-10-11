.. _outlier-detection-spec:

Outlier Detection for Slit-like Spectroscopic Data
==================================================

This module serves as the interface for applying ``outlier_detection`` to slit-like
spectroscopic observations. The algorithm shares many similarities with the
:ref:`imaging algorithm <outlier-detection-imaging>`, and much of the same code is used.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a :py:class:`~jwst.datamodels.ModelContainer` object,
serves as the input and output to this step, and the :py:class:`~jwst.datamodels.ModelContainer`
is converted into a :py:class:`~jwst.datamodels.ModelLibrary` object where processing is
shared with the imaging mode.

This routine performs identical operations to the imaging mode, with the following exceptions:

#. Resampling is handled by a different class, :py:class:`~jwst.resample.resample_spec.ResampleSpecData`
   instead of :py:class:`~jwst.resample.resample.ResampleData`.

#. The resampled images are written out to disk with suffix "outlier_s2d" instead of 
   "outlier_i2d" if the ``save_intermediate_results`` parameter is set to `True`.

#. The ``in_memory`` parameter has no effect, and all operations are performed in memory.

.. automodapi:: jwst.outlier_detection.spec
