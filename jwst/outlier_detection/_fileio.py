import re
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def save_median(median_model, make_output_path, asn_id=None):
    '''
    Save median if requested by user

    Parameters
    ----------
    median_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save
    '''
    # need unique slit ID for MultiSlitModel inputs to step
    slit_id =  getattr(median_model, "name", "").lower()
    if slit_id:
        slit_id = "_"+slit_id
    # regex to handle i2d and s2d as though they are the same string
    default_suffix = r"_outlier_.2d\.fits"
    suffix_to_remove = default_suffix if asn_id is None else fr"_{asn_id}{default_suffix}"
    basepath = re.sub(suffix_to_remove, ".fits", median_model.meta.filename)
    median_model_output_path = make_output_path(
        basepath=basepath,
        suffix='median')
    median_model.save(median_model_output_path)
    log.info(f"Saved model in {median_model_output_path}")
