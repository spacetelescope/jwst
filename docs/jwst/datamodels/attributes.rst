Data model attributes
=====================
The purpose of the data model is to abstract away the peculiarities of
the underlying file format.  The same data model may be used for data
created from scratch in memory, loaded from FITS or ASDF files, or from
some other future format.

Calling sequences of models
===========================

List of current models
----------------------

The current models are as follows:

    'ABVegaOffsetModel',
    'AmiLgModel',
    'FgsImgApcorrModel', 'MirImgApcorrModel', 'NrcImgApcorrModel', 'NisImgApcorrModel',
    'MirLrsApcorrModel', 'MirMrsApcorrModel', 'NrcWfssApcorrModel', 'NisWfssApcorrModel',
    'NrsMosApcorrModel', 'NrsFsApcorrModel', 'NrsIfuApcorrModel',
    'AsnModel',
    'BarshadowModel', 'CameraModel', 'CollimatorModel',
    'CombinedSpecModel', 'ContrastModel', 'CubeModel',
    'DarkModel', 'DarkMIRIModel',
    'DisperserModel', 'DistortionModel', 'DistortionMRSModel',
    'DrizParsModel',
    'Extract1dImageModel',
    'Extract1dIFUModel',
    'FilteroffsetModel',
    'FlatModel', 'NirspecFlatModel', 'NirspecQuadFlatModel',
    'FOREModel', 'FPAModel',
    'FringeModel', 'GainModel', 'GLS_RampFitModel',
    'GuiderRawModel', 'GuiderCalModel',
    'IFUCubeModel',
    'NirspecIFUCubeParsModel', 'MiriIFUCubeParsModel',
    'IFUFOREModel', 'IFUImageModel', 'IFUPostModel', 'IFUSlicerModel',
    'ImageModel', 'IPCModel', 'IRS2Model', 'LastFrameModel', 'Level1bModel',
    'LinearityModel', 'MaskModel', 'ModelContainer', 'MSAModel',
    'MultiCombinedSpecModel', 'MultiExposureModel',
    'MultiExtract1dImageModel', 'MultiSlitModel',
    'MultiSpecModel',
    'NIRCAMGrismModel', 'NIRISSGrismModel',
    'OTEModel',
    'OutlierParsModel',
    'PathlossModel',
    'PersistenceSatModel',
    'PixelAreaModel', 'NirspecSlitAreaModel', 'NirspecMosAreaModel', 'NirspecIfuAreaModel',
    'FgsImgPhotomModel',
    'MirImgPhotomModel', 'MirLrsPhotomModel', 'MirMrsPhotomModel',
    'NrcImgPhotomModel', 'NrcWfssPhotomModel',
    'NisImgPhotomModel', 'NisSossPhotomModel', 'NisWfssPhotomModel',
    'NrsFsPhotomModel', 'NrsMosPhotomModel',
    'PsfMaskModel',
    'QuadModel', 'RampModel', 'MIRIRampModel',
    'RampFitOutputModel', 'ReadnoiseModel',
    'ReferenceFileModel', 'ReferenceCubeModel', 'ReferenceImageModel', 'ReferenceQuadModel',
    'RegionsModel', 'ResetModel',
    'ResolutionModel', 'MiriResolutionModel',
    'RSCDModel', 'SaturationModel', 'SlitDataModel', 'SlitModel', 'SpecModel',
    'SegmentationMapModel', 'SourceModelContainer',
    'SpecKernelModel',
    'SpecProfileModel', 'SpecProfileSingleModel',
    'SpecTraceModel', 'SpecTraceSingleModel',
    'SpecwcsModel',
    'StrayLightModel', 'SuperBiasModel',
    'ThroughputModel',
    'TrapDensityModel', 'TrapParsModel', 'TrapsFilledModel',
    'TsoPhotModel',
    'WavelengthrangeModel', 'WaveCorrModel',
    'WaveMapModel', 'WaveMapSingleModel',
    'WfssBkgModel'

Commonly used attributes
------------------------
Here are a few model attributes that are used by some of the pipeline
steps.

For uncalibrated data `_uncal.fits`.  Getting the number of integrations
and the number of groups from the first and second axes assumes that the
input data array is 4-D data.  Pixel coordinates in the data extensions are
1-indexed as in FORTRAN and FITS headers, not 0-indexed as in Python.

    - ``input_model.data.shape[0]``: number of integrations
    - ``input_model.data.shape[1]``: number of groups
    - ``input_model.meta.exposure.nframes``: number of frames per group
    - ``input_model.meta.exposure.groupgap``: number of frames dropped
        between groups
    - ``input_model.meta.subarray.xstart``: starting pixel in X (1-based)
    - ``input_model.meta.subarray.ystart``: starting pixel in Y (1-based)
    - ``input_model.meta.subarray.xsize``: number of columns
    - ``input_model.meta.subarray.ysize``: number of rows

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

DataModel Base Class
--------------------

.. autoclass:: jwst.datamodels.DataModel
   :members:
   :noindex:
