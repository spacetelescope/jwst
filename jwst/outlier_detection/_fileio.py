import os

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def remove_file(fn):
    if isinstance(fn, str) and os.path.isfile(fn):
        os.remove(fn)
        log.info(f"Removing file {fn}")


def save_median(median_model, make_output_path, asn_id=None):
    '''
    Save median if requested by user

    Parameters
    ----------
    median_model : ~jwst.datamodels.ImageModel
        The median ImageModel or CubeModel to save
    '''
    default_suffix = "_outlier_i2d.fits"
    if asn_id is None:
        suffix_to_remove = default_suffix
    else:
        suffix_to_remove = f"_{asn_id}{default_suffix}"
    median_model_output_path = make_output_path(
        basepath=median_model.meta.filename.replace(suffix_to_remove, '.fits'),
        suffix='median')
    median_model.save(median_model_output_path)
    log.info(f"Saved model in {median_model_output_path}")
