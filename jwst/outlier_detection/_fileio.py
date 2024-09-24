import re
import logging
from jwst.datamodels import ImageModel
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(median_data, median_wcs, example_model, make_output_path=None):
    '''
    Save median if requested by user

    Parameters
    ----------
    median_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save
    '''
    median_model = ImageModel(median_data)
    median_model.update(example_model)
    median_model.meta.wcs = median_wcs
    output_path = intermediate_output_path(example_model, "median", make_output_path)
    _save_to_path(median_model, output_path, suffix="median")


def save_drizzled(drizzled_model, make_output_path):
    expected_tail = "outlier_?2d.fits"
    suffix = drizzled_model.meta.filename[-len(expected_tail):-5]
    output_path = intermediate_output_path(drizzled_model, suffix, make_output_path)
    _save_to_path(drizzled_model, output_path, suffix=suffix)


def save_blot(input_model, blot, make_output_path):

    blot_model = _make_blot_model(input_model, blot)
    output_path = intermediate_output_path(input_model, "blot", make_output_path, suffix_to_remove=r"_cal.*")
    _save_to_path(blot_model, output_path, suffix="blot")


def _make_blot_model(input_model, blot):
    blot_model = type(input_model)()
    blot_model.data = blot
    blot_model.update(input_model)
    return blot_model


def intermediate_output_path(input_model, suffix, make_output_path, suffix_to_remove=r"_outlier.*"):
    """Ensure all intermediate outputs from OutlierDetectionStep have consistent file naming conventions"""
    if hasattr(input_model, "name") and input_model.name is not None:
        if "_"+input_model.name.lower()+"_" not in input_model.meta.filename:
            replacement = f"_{input_model.name.lower()}"
        else:
            replacement = ""
    else:
        replacement = ""
    basepath = re.sub(suffix_to_remove, replacement, input_model.meta.filename)

    if make_output_path is None:
        def make_output_path(basepath, suffix):
            return basepath.replace(".fits", f"_{suffix}.fits")

    return make_output_path(basepath, suffix=suffix)


def _save_to_path(model, output_path, suffix=""):
    model.meta.filename = output_path
    model.save(output_path)
    log.info(f"Saved {suffix} model in {output_path}")