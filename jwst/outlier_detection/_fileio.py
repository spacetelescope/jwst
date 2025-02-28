import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(median_model, make_output_path):
    """
    Save a median model.

    Output suffix is 'median'.

    Parameters
    ----------
    median_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save

    make_output_path : function
        A function to be used to create an output path from
        an input path and an optional suffix.
    """
    _save_intermediate_output(median_model, "median", make_output_path)


def save_drizzled(drizzled_model, make_output_path):
    """
    Save a drizzled model.

    Input model is expected to have a filename that ends with
    either 'outlier_i2d.fits' or 'outlier_s2d.fits'.

    Parameters
    ----------
    drizzled_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save.

    make_output_path : function
        A function to be used to create an output path from
        an input path and an optional suffix.
    """
    expected_tail = "outlier_?2d.fits"
    suffix = drizzled_model.meta.filename[-len(expected_tail) : -5]
    _save_intermediate_output(drizzled_model, suffix, make_output_path)


def save_blot(input_model, blot, blot_err, make_output_path):
    """
    Save a blotted model with output suffix 'blot'.

    Parameters
    ----------
    input_model : ~jwst.datamodels.ImageModel
        An input model corresponding to the blotted data,
        containing metadata to copy.

    blot : array-like
        The blotted science data.

    blot_err : array-like or None
        The blotted error data, if available.

    make_output_path : function
        A function to be used to create an output path from
        an input path and an optional suffix.
    """
    blot_model = _make_blot_model(input_model, blot, blot_err)
    _save_intermediate_output(blot_model, "blot", make_output_path)


def _make_blot_model(input_model, blot, blot_err):
    """
    Assemble a blot model.

    Parameters
    ----------
    input_model : ~jwst.datamodels.ImageModel
        An input model corresponding to the blotted data,
        containing metadata to copy.

    blot : array-like
        The blotted science data.

    blot_err : array-like or None
        The blotted error data, if available.

    Returns
    -------
    blot_model : ~jwst.datamodels.ImageModel
        An image model containing the blotted data.
    """
    blot_model = type(input_model)()
    blot_model.data = blot
    if blot_err is not None:
        blot_model.err = blot_err
    blot_model.update(input_model)
    return blot_model


def _save_intermediate_output(model, suffix, make_output_path):
    """
    Save an intermediate output from outlier detection.

    Ensure all intermediate outputs from OutlierDetectionStep have
    consistent file naming conventions.

    Parameters
    ----------
    model : ~jwst.datamodels.ImageModel
        The intermediate datamodel to save.

    suffix : str
        A suffix to add to the output file name.

    make_output_path : function
        A function to be used to create an output path from
        an input path and an optional suffix.

    Notes
    -----
    self.make_output_path() is updated globally for the step in the main pipeline
    to include the asn_id in the output path, so no need to handle it here.
    """
    # outlier_?2d is not a known suffix, and make_output_path cannot handle an
    # underscore in an unknown suffix, so do a manual string replacement
    input_path = model.meta.filename.replace("_outlier_", "_")

    # Add a slit name to the output path for MultiSlitModel data if not present
    if hasattr(model, "name") and model.name is not None:
        if "_" + model.name.lower() not in input_path:
            suffix = f"{model.name.lower()}_{suffix}"

    output_path = make_output_path(input_path, suffix=suffix)
    model.save(output_path)
    log.info(f"Saved {suffix} model in {output_path}")
