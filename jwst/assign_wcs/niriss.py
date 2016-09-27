from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging

from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Identity, Const1D

import gwcs.coordinate_frames as cf
from . import pointing
from .util import not_implemented_mode


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def create_pipeline(input_model, reference_files):
    '''
    get reference files from crds

    '''

    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def soss_distortion_b7(reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree[1]
    return distortion


def nis_soss_b7(input_model, reference_files):
    """
    The NIRISS SOSS pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'soss_wavelengths_configuration.asdf'}
    """

    # Get the target RA and DEC, they will be used for setting the WCS RA and DEC based on a conversation
    # with Kevin Volk.
    target_ra, target_dec = 0.0, 0.0
    try:
        target_ra = float(input_model['meta.target.ra'])
        target_dec = float(input_model['meta.target.dec'])
    except:
        # There was an error getting the target RA and DEC, so we are not going to continue.
        raise ValueError('Problem getting the TARG_RA or TARG_DEC from input model {}'.format(input_model))

    # Define the frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS())

    # Define the transforms
    distortion = soss_distortion_b7(reference_files)

    # Do a quick check to make sure the orientation of the data is the same
    # in each file (DMS etc etc)
    #ys, xs = mgrid[:2048, :2048]
    #ww = distortion(xs, ys)
    #if nonzero(input_model.data>0):


    # Now set the RA and DEC based on the target.
    distortion_rdl = (Const1D(target_ra) & Const1D(target_dec) & distortion)

    pipeline = [(detector, distortion_rdl),
                (sky, None)
                ]

    return pipeline


def nis_soss_b7p1(input_model, reference_files):
    """
    The NIRISS SOSS pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'soss_wavelengths_configuration.asdf'}
    """

    log.info('Input is {}'.format(input_model))
    log.info('\tinput.data.shape {}'.format(input_model.data.shape))
    log.info('reference_files is {}'.format(reference_files))

    # We are going to assign the target RA, DEC as the output RA, DEC, per Kevin Volk.
    target_ra, target_dec = 0.0, 0.0

    # We need the PWCPOS (angle in degrees) from the SCI data
    pwcpos_sci = 0

    # Need the PWCPOS (angle in degrees) from the wavelength correction file
    pwcpos_wave_corr = 0

    # Need the PWCPOS (angle in degrees) from the wavelength standard file
    pwcpos_wave_standard = 0

    # Define the frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    focal = cf.Frame2D(name='focal', axes_order=(0, 1), unit=(u.arcmin, u.arcmin))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS())

    # Define the transforms
    distortion = imaging_distortion(input_model, reference_files)


    pipeline = [(detector, distortion),
                (focal, fitswcs_transform),
                (sky, None)
                ]

    return pipeline


def imaging(input_model, reference_files):
    """
    The NIRISS imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'jwst_niriss_distortioon_0001.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    focal = cf.Frame2D(name='focal', axes_order=(0, 1), unit=(u.arcmin, u.arcmin))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS())
    distortion = imaging_distortion(input_model, reference_files)
    fitswcs_transform = pointing.create_fitswcs_transform(input_model)
    pipeline = [(detector, distortion),
                (focal, None)]
                #(sky, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    return distortion


exp_type2transform = {'nis_image': imaging,
                      'nis_wfss': not_implemented_mode,#        ?? WFSS spec
                      'nis_soss': nis_soss_b7,
                      'nis_ami': not_implemented_mode#        ?? imaging
                      }
