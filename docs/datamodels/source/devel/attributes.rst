.. py:module:: jwst_lib.models

Data model attributes
`````````````````````

The purpose of the data model is to abstract away the peculiarities of
the underlying file format.  The same data model may be used for data
created from scratch in memory, or loaded from FITS files or some
future FITS-replacement format.

Calling sequences of models
===========================

List of current models
----------------------

The current models are as follows:

    `AmiLgModel`, `AsnModel`, `ContrastModel`,
    `CubeModel`, `DarkModel`, `DrizParsModel`, `NircamDrizParsModel`,
    `MiriImgDrizParsModel`, `DrizProductModel`, `FilterModel`,
    `FlatModel`, `FringeModel`, `GainModel`, `GLS_RampFitModel`,
    `ImageModel`, `IPCModel`, `LastFrameModel`, `LinearityModel`,
    `MaskModel`, `MIRIRampModel`, `MultiSlitModel`, `MultiSpecModel`,
    `PhotomModel`, `NircamPhotomModel`, `NirissPhotomModel`,
    `NirspecPhotomModel`, `MiriImgPhotomModel`, `MiriMrsPhotomModel`,
    `RampModel`, `RampFitOutputModel`, `ReadnoiseModel`, `ResetModel`,
    `SaturationModel`, `SpecModel`, `StrayLightModel`

Commonly used attributes
------------------------

Here are a few model attributes that are used by many of the pipeline
steps.  Note that, following FORTRAN and FITS conventions, the
starting pixel numbers in X and Y are one-indexed.  Getting the number
of integrations and the number of groups from the first and second
axes assumes that the input data array is 4-D data.  Much of the
jwst_pipeline step code assumes that the data array is 4-D.

    - number of integrations = input_model.data.shape[0]
    - number of groups = input_model.data.shape[1]
    - number of frames = input_model.meta.exposure.nframes
    - group gap = input_model.meta.exposure.groupgap
    - starting pixel in X = input_model.meta.subarray.xstart
    - starting pixel in Y = input_model.meta.subarray.ystart
    - number of columns = input_model.meta.subarray.xsize
    - number of rows = input_model.meta.subarray.ysize

The `data`, `err`, `dq`, etc., attributes of most models are assumed to be
numpy.ndarray arrays, or at least objects that have some of the attributes
of these arrays.  numpy is used explicitly to create these arrays in some
cases (e.g. when a default value is needed).  The `data` and `err` arrays
are a floating point type, and the data quality arrays are an integer type.

Some of the step code makes assumptions about image array sizes.  For
example, full-frame MIRI data have 1032 columns and 1024 rows, and all
other detectors have 2048 columns and rows; anything smaller must be a
subarray.  Also, full-frame MIRI data are assumed to have four columns of
reference pixels on the left and right sides (the reference output array
is stored in a separate image extension).  Full-frame data for all other
instruments have four columns or rows of reference pixels on each edge
of the image.

Model classes
-------------

Base class
''''''''''

.. autoclass:: jwst_lib.models.DataModel
   :members:

Concrete model classes
''''''''''''''''''''''

.. automodule:: jwst_lib.models
   :members: AmiLgModel, AsnModel, ContrastModel,
    CubeModel, DarkModel, DrizParsModel, NircamDrizParsModel,
    MiriImgDrizParsModel, DrizProductModel, FilterModel,
    FlatModel, FringeModel, GainModel, GLS_RampFitModel,
    ImageModel, IPCModel, LastFrameModel, LinearityModel,
    MaskModel, MIRIRampModel, MultiSlitModel, MultiSpecModel,
    NircamPhotomModel, NirissPhotomModel,
    NirspecPhotomModel, MiriImgPhotomModel, MiriMrsPhotomModel,
    RampModel, RampFitOutputModel, ReadnoiseModel, ResetModel,
    SaturationModel, SpecModel, StrayLightModel
