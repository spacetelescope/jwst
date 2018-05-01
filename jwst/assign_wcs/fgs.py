"""
FGS WCS pipeline - depends on EXP_TYPE.
"""
import logging

from astropy import units as u
from astropy import coordinates as coord
from gwcs import coordinate_frames as cf

from .util import not_implemented_mode, subarray_transform
from . import pointing
from ..datamodels import DistortionModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_pipeline(input_model, reference_files):
    """
    Create a pipeline list based on EXP_TYPE.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Either an ImageModel or a CubeModel
    reference_files : dict
        {reftype: file_name} mapping
        In the pipeline it's returned by CRDS.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    log.info("Creating a FGS {0} pipeline with references {1}".format(
        exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    The FGS imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky.

    reference_files={'distortion': 'jwst_fgs_distortioon_0001.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(name='world', reference_frame=coord.ICRS())
    # V2, V3 to sky
    tel2sky = pointing.v23tosky(input_model)

    subarray2full = subarray_transform(input_model)
    if reference_files:
        imdistortion = imaging_distortion(input_model, reference_files)
        distortion = subarray2full | imdistortion
        distortion.bounding_box = imdistortion.bounding_box
        del imdistortion.bounding_box
    else:
        distortion = subarray2full

    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    dist = DistortionModel(reference_files['distortion'])
    transform = dist.model
    try:
        bb = transform.bounding_box
    except NotImplementedError:
        shape = input_model.data.shape
        # Note: Since bounding_box is attached to the model here
        # it's in reverse order.
        """
        A CubeModel is always treated as a stack (in dimension 1)
        of 2D images, as opposed to actual 3D data. In this case
        the bounding box is set to the 2nd and 3rd dimension.
        """
        if isinstance(input_model, CubeModel):
            bb = ((-0.5, shape[1] - 0.5),
                  (-0.5, shape[2] - 0.5))
        elif isinstance(input_model, ImageModel):
            bb = ((-0.5, shape[0] - 0.5),
                  (-0.5, shape[1] - 0.5))
        else:
            raise TypeError("Input is not an ImageModel or CubeModel")

        transform.bounding_box = bb
    dist.close()
    return transform


exp_type2transform = {'fgs_image': imaging,
                      'fgs_focus': imaging,
                      'fgs_skyflat': not_implemented_mode,
                      'fgs_intflat': not_implemented_mode,
                      'fgs_dark': not_implemented_mode
                      }
