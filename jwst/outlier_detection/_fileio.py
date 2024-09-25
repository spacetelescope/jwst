import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(median_model, make_output_path):
    '''
    Save median if requested by user

    Parameters
    ----------
    median_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save
    '''
    _save_intermediate_output(median_model, "median", make_output_path)


def save_drizzled(drizzled_model, make_output_path):
    expected_tail = "outlier_?2d.fits"
    suffix = drizzled_model.meta.filename[-len(expected_tail):-5]
    _save_intermediate_output(drizzled_model, suffix, make_output_path)


def save_blot(input_model, blot, make_output_path):
    blot_model = _make_blot_model(input_model, blot)
    _save_intermediate_output(blot_model, "blot", make_output_path)


def _make_blot_model(input_model, blot):
    blot_model = type(input_model)()
    blot_model.data = blot
    blot_model.update(input_model)
    return blot_model


def _save_intermediate_output(model, suffix, make_output_path):
    """
    Ensure all intermediate outputs from OutlierDetectionStep have consistent file naming conventions
    
    Notes
    -----
    self.make_output_path() is updated globally for the step in the main pipeline
    to include the asn_id in the output path, so no need to handle it here.
    """

    # make_output_path cannot handle suffix with an underscore inside it,
    # e.g. _outlier_s2d.fits, so do a manual string replacement
    input_path = model.meta.filename.replace("_outlier_", "_")

    # Add a slit name to the output path for MultiSlitModel data if not present
    if hasattr(model, "name") and model.name is not None:
        if "_"+model.name.lower() not in input_path:
            suffix = f"{model.name.lower()}_{suffix}"

    output_path = make_output_path(input_path, suffix=suffix)
    model.meta.filename = output_path
    model.save(output_path)
    log.info(f"Saved {suffix} model in {output_path}")
