from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os.path
import logging
import numpy as np
from asdf import AsdfFile
from astropy.modeling import models
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits

import gwcs.coordinate_frames as cf
from gwcs import selector
from . import pointing
from ..transforms import models as jwmodels
from .util import not_implemented_mode, subarray_transform

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_pipeline(input_model, reference_files):
    '''
    Create the WCS pipeline for MIRI modes.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    '''
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """
    Create MIRI Imagng WCS.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    The MIRI imaging pipeline includes 3 coordinate frames - detector,
    focal plane and sky

    reference_files={'distortion': 'test.asdf', 'filter_offsets': 'filter_offsets.asdf'}
    """

    # Create the Frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')

    # Create the transforms
    subarray2full = subarray_transform(input_model)
    imdistortion = imaging_distortion(input_model, reference_files)
    distortion = subarray2full | imdistortion
    tel2sky = pointing.v23tosky(input_model)

    # TODO: remove setting the bounding box when it is set in the new ref file.
    try:
        bb = distortion.bounding_box
    except NotImplementedError:
        shape = 1032,1024
        # Note: Since bounding_box is attached to the model here it's in reverse order.
        distortion.bounding_box = ((-0.5, shape[1] - 0.5), (3.5, shape[0] - 4.5))

    # Create the pipeline
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)
                ]

    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create pixe2sky and sky2pixel transformation for the MIRI imager.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    1. Filter dependent shift in (x,y) (!with an oposite sign to that delivered by the IT)
    2. Apply MI
    3. Apply Ai and BI matrices
    4. Apply the TI matrix (this gives Xan/Yan coordinates)
    5. Aply the XanYan --> V2V3 transform
    6. Apply V2V3 --> sky transform

    ref_file: filter_offset.asdf - (1)
    ref_file: distortion.asdf -(2,3,4)
    """
    # Read in the distortion.
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    obsfilter = input_model.meta.instrument.filter

    # Add an offset for the filter
    with AsdfFile.open(reference_files['filteroffset']) as filter_offset:
        if obsfilter in filter_offset.tree:
            filter_corr = filter_offset.tree[obsfilter]
            distortion = models.Shift(filter_corr['column_offset']) & models.Shift(
                filter_corr['row_offset']) | distortion

    # scale to degrees (remove this once pipeline can handle v2,v3 in arcsec)
    distortion = distortion | models.Scale(1 / 3600) & models.Scale(1 / 3600)
    return distortion


def lrs(input_model, reference_files):
    """
    Create the WCS pipeline for a MIRI fixed slit observation.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    reference_files = {
        "specwcs": 'MIRI_FM_MIRIMAGE_P750L_DISTORTION_04.02.00.fits'
    }
    """

    # Setup the frames.
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    spec = cf.SpectralFrame(name='wavelength', axes_order=(2,), unit=(u.micron,),
                            axes_names=('lambda',))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(), name='sky')
    world = cf.CompositeFrame(name="world", frames=[sky, spec])


    # Determine the distortion model.
    subarray2full = subarray_transform(input_model)
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    # Distortion is in arcsec.  Convert to degrees
    full_distortion = subarray2full | distortion | models.Scale(1 / 3600.) & models.Scale(1 / 3600.)

    # Load and process the reference data.
    with fits.open(reference_files['specwcs']) as ref:
        lrsdata = np.array([l for l in ref[1].data])

        # Get the zero point from the reference data.
        # The zero_point is X, Y  (which should be COLUMN, ROW)
        # TODO: Are imx, imy 0- or 1-indexed?  We are treating them here as
        # 0-indexed.  Since they are FITS, they are probably 1-indexed.
        if input_model.meta.exposure.type.lower() == 'mir_lrs-fixedslit':
            zero_point = ref[1].header['imx'], ref[1].header['imy']
        elif input_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
            #zero_point = ref[1].header['imxsltl'], ref[1].header['imysltl']
            zero_point = [35, 442]  # [35, 763] # account for subarray

    # Create the bounding_box
    x0 = lrsdata[:, 3]
    y0 = lrsdata[:, 4]
    x1 = lrsdata[:, 5]

    bb = ((x0.min() - 0.5 + zero_point[0], x1.max() + 0.5 + zero_point[0]),
          (y0.min() - 0.5 + zero_point[1], y0.max() + 0.5 + zero_point[1]))
    # Find the ROW of the zero point which should be the [1] of zero_point
    row_zero_point = zero_point[1]

    # Compute the v2v3 to sky.
    tel2sky = pointing.v23tosky(input_model)

    # Compute the V2/V3 for each pixel in this row
    # x.shape will be something like (1, 388)
    y, x = np.mgrid[row_zero_point:row_zero_point + 1, 0:input_model.data.shape[1]]

    spatial_transform = full_distortion | tel2sky
    radec = np.array(spatial_transform(x, y))[:, 0, :]

    ra_full = np.matlib.repmat(radec[0], bb[1][1] + 1 - bb[1][0], 1)
    dec_full = np.matlib.repmat(radec[1], bb[1][1] + 1 - bb[1][0], 1)

    ra_t2d = models.Tabular2D(lookup_table=ra_full, name='xtable',
        bounds_error=False, fill_value=np.nan)
    dec_t2d = models.Tabular2D(lookup_table=dec_full, name='ytable',
        bounds_error=False, fill_value=np.nan)

    # Create the model transforms.
    lrs_wav_model = jwmodels.LRSWavelength(lrsdata, zero_point)

    # Incorporate the small rotation
    angle = np.arctan(0.00421924)
    rot = models.Rotation2D(angle)
    radec_t2d = ra_t2d & dec_t2d | rot

    # Account for the subarray when computing spatial coordinates.
    xshift = -bb[0][0]
    yshift = -bb[1][0]
    det2world = models.Mapping((1, 0, 1, 0, 0, 1)) | models.Shift(yshift, name='yshift1') & \
              models.Shift(xshift, name='xshift1') & \
              models.Shift(yshift, name='yshift2') & models.Shift(xshift, name='xshift2') & \
              models.Identity(2) | radec_t2d & lrs_wav_model
    det2world.bounding_box = bb[::-1]
    # Now the actual pipeline.
    pipeline = [(detector, det2world),
                (world, None)
                ]

    return pipeline


def ifu(input_model, reference_files):
    """
    Create the WCS pipeline for a MIRI IFU observation.
    Goes from 0-indexed detector pixels (0,0) middle of lower left reference pixel
    to V2,V3 in arcsec.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.
    """

    #reference_files = {'distortion': 'jwst_miri_distortion_00001.asdf', #files must hold 2 channels each
                        #'specwcs': 'jwst_miri_specwcs_00001.asdf',
                        #'regions': 'jwst_miri_regions_00001.asdf',
                        #'wavelengthrange': 'jwst_miri_wavelengthrange_0001.asdf'}
    # Define reference frames
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    alpha_beta = cf.Frame2D(name='alpha_beta_spatial', axes_order=(0, 1), unit=(u.arcsec, u.arcsec), axes_names=('alpha', 'beta'))
    spec_local = cf.SpectralFrame(name='alpha_beta_spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    miri_focal = cf.CompositeFrame([alpha_beta, spec_local], name='alpha_beta')
    v23_spatial = cf.Frame2D(name='V2_V3_spatial', axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('v2', 'v3'))
    spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,), axes_names=('lambda',))
    v2v3 = cf.CompositeFrame([v23_spatial, spec], name='v2v3')
    icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
                             axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
    world = cf.CompositeFrame([icrs, spec], name='world')

    # Define the actual transforms
    det2abl = (detector_to_abl(input_model, reference_files)).rename(
        "detector_to_abl")
    abl2v2v3l = (abl_to_v2v3l(input_model, reference_files)).rename("abl_to_v2v3l")
    tel2sky = pointing.v23tosky(input_model) & models.Identity(1)

    # Put the transforms together into a single transform
    shape = input_model.data.shape
    det2abl.bounding_box = ((-0.5, shape[0] - 0.5), (-0.5, shape[1] - 0.5))
    pipeline = [(detector, det2abl),
                (miri_focal, abl2v2v3l),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline

def detector_to_abl(input_model, reference_files):
    """
    Create the transform from detector to alpha, beta frame.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    forward transform:
      RegionsSelector
        label_mapper is the regions array
        selector is {slice_number: alphs_model & beta_model & lambda_model}
    backward transform:
      RegionsSelector
        label_mapper is LabelMapperDict
           {channel_wave_range (): LabelMapperDict}
                                   {beta: slice_number}
        selector is {slice_number: x_transform & y_transform}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    f = AsdfFile.open(reference_files['distortion'])
    # The name of the model indicates the output coordinate
    alpha_model = f.tree['alpha_model']
    beta_model = f.tree['beta_model']
    x_model = f.tree['x_model']
    y_model = f.tree['y_model']
    bzero = f.tree['bzero']
    bdel = f.tree['bdel']
    f.close()
    f = AsdfFile.open(reference_files['specwcs'])
    lambda_model = f.tree['model']
    f.close()
    f = AsdfFile.open(reference_files['regions'])
    regions = f.tree['regions'].copy()
    f.close()
    label_mapper = selector.LabelMapperArray(regions)
    transforms = {}

    for sl in alpha_model:
        forward = models.Mapping([1, 0, 0, 1, 0]) | \
                alpha_model[sl] & beta_model[sl] & lambda_model[sl]
        inv = models.Mapping([2, 0, 2, 0]) | x_model[sl] & y_model[sl]
        forward.inverse = inv
        transforms[sl] = forward

    f = AsdfFile.open(reference_files['wavelengthrange'])
    # the following should go in the asdf reader
    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    wr = {}
    for ch, r in zip(wave_channels, wave_range):
        wr[ch] = r
    f.close()
    ch_dict = {}
    for c in channel:
        cb = c + band
        mapper = jwmodels.MIRI_AB2Slice(bzero[cb], bdel[cb], c)
        lm = selector.LabelMapper(inputs=('alpha', 'beta', 'lam'),
                                  mapper=mapper, inputs_mapping=models.Mapping((1,), n_inputs=3))
        ch_dict[tuple(wr[cb])] = lm


    alpha_beta_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), ch_dict,
                                                  models.Mapping((2,)))
    label_mapper.inverse = alpha_beta_mapper


    det2alpha_beta = selector.RegionsSelector(('x', 'y'), ('alpha', 'beta', 'lam'),
                                              label_mapper=label_mapper, selector=transforms)
    return det2alpha_beta

def abl_to_v2v3l(input_model, reference_files):
    """
    Create the transform from (alpha,beta,lambda) to (V2,V3,lambda) frame.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model.
    reference_files : dict
        Dictionary {reftype: reference file name}.

    forward transform:
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: ab2v2 & ab2v3}
    bacward_transform
      RegionsSelector
        label_mapper is LabelMapperDict()
        {channel_wave_range (): channel_number}
        selector is {channel_number: v22ab & v32ab}
    """
    band = input_model.meta.instrument.band
    channel = input_model.meta.instrument.channel
    # used to read the wavelength range
    channels = [c + band for c in channel]

    f = AsdfFile.open(reference_files['distortion'])
    v23 = f.tree['abv2v3_model']
    f.close()
    f = AsdfFile.open(reference_files['wavelengthrange'])
    # the following should go in the asdf reader
    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    wr = dict(zip(wave_channels, wave_range))
    f.close()

    dict_mapper = {}
    sel = {}
    # Since there are two channels in each reference file we need to loop over them
    for c in channels:
        ch = int(c[0])
        dict_mapper[tuple(wr[c])] = models.Mapping((2,), name="mapping_lam") | \
                   models.Const1D(ch, name="channel #")
        ident1 = models.Identity(1, name='identity_lam')
        ident1._inputs = ('lam',)
        chan_v23 = v23[c]
        v23chan_backward = chan_v23.inverse
        del chan_v23.inverse
        # This is the spatial part of the transform; tack on additional conversion to degrees
        # Remove this degrees conversion once pipeline can handle v2,v3 in arcsec
        v23_spatial=chan_v23 | models.Scale(1 / 3600) & models.Scale(1 / 3600)
        v23_spatial.inverse = models.Scale(3600) & models.Scale(3600) | v23chan_backward
        # Tack on passing the third wavelength component
        v23c = v23_spatial & ident1
        sel[ch] = v23c

    wave_range_mapper = selector.LabelMapperRange(('alpha', 'beta', 'lam'), dict_mapper,
                                                  inputs_mapping=models.Mapping([2, ]))
    wave_range_mapper.inverse = wave_range_mapper.copy()
    abl2v2v3l = selector.RegionsSelector(('alpha', 'beta', 'lam'), ('v2', 'v3', 'lam'),
                                      label_mapper=wave_range_mapper,
                                      selector=sel)

    return abl2v2v3l


exp_type2transform = {'mir_image': imaging,
                      'mir_tacq': imaging,
                      'mir_lyot': imaging,
                      'mir_4qpm': imaging,
                      'mir_coroncal': imaging,
                      'mir_lrs-fixedslit': lrs,
                      'mir_lrs-slitless': lrs,
                      'mir_mrs': ifu,
                      'mir_flatmrs': not_implemented_mode,
                      'mir_flatimage': not_implemented_mode,
                      'mir_flat-mrs': not_implemented_mode,
                      'mir_flat-image': not_implemented_mode,
                      'mir_dark': not_implemented_mode,
                      }


def get_wavelength_range(input_model, path=None):
    """
    Return the wavelength range used for computing the WCS.

    Needs access to the reference file used to construct the WCS object.

    Parameters
    ----------
    input_model : `jwst.datamodels.ImagingModel`
        Data model after assign_wcs has been run.
    path : str
        Directory where the reference file is. (optional)
    """
    fname = input_model.meta.ref_file.wavelengthrange.name.split('/')[-1]
    if path is None and not os.path.exists(fname):
        raise IOError("Reference file {0} not found. Please specify a path.".format(fname))
    else:
        fname = os.path.join(path, fname)
        f = AsdfFile.open(fname)

    wave_range = f.tree['wavelengthrange'].copy()
    wave_channels = f.tree['channels']
    f.close()

    wr = dict(zip(wave_channels, wave_range))
    channel = input_model.meta.instrument.channel
    band = input_model.meta.instrument.band

    return dict([(ch+band, wr[ch+band]) for ch in channel ])
