.. _outlier-detection-spec:

Slit-like Spectroscopic Data
============================

This module serves as the interface for applying ``outlier_detection`` to slit-like
spectroscopic observations. The algorithm is very similar to the
:ref:`imaging algorithm <outlier-detection-imaging>`, and much of the same code is used.
Please refer to those docs for more information.
A :ref:`Stage 3 association <asn-level3-techspecs>`,
which is loaded into a :py:class:`~jwst.datamodels.ModelContainer` object,
serves as the input and output to this step, and the :py:class:`~jwst.datamodels.ModelContainer`
is converted into a :py:class:`~jwst.datamodels.ModelLibrary` object to allow sharing code
with the imaging mode.

This routine performs identical operations to the imaging mode, with the following exceptions:

#. Error thresholding is handled differently: the error arrays are resampled and median-combined
   along with the data arrays, and the median error image is used to identify outliers
   instead of the input error images for each exposure. This median error image is included
   alongside the median datamodel (in the ``err`` extension) if ``save_intermediate_results``
   is `True`.

#. Resampling is handled by a different class, :py:class:`~jwst.resample.resample_spec.ResampleSpec`
   instead of :py:class:`~jwst.resample.resample.ResampleImage`.

#. The resampled images are written out to disk with suffix "outlier_s2d" instead of
   "outlier_i2d" if the ``save_intermediate_results`` parameter is set to `True`.

#. The ``in_memory`` parameter has no effect, and all operations are performed in memory.

.. automodapi:: jwst.outlier_detection.spec
