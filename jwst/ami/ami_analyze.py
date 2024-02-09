#  Module for applying the LG-PLUS algorithm to an AMI exposure
import logging
import numpy as np
import copy
import synphot

from .find_affine2d_parameters import find_rotation
from . import instrument_data
from . import nrm_core
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_LG_plus(input_model, 
                oversample, rotation,
                psf_offset, rotsearch_parameters, 
                src, bandpass, usebp, firstfew, 
                chooseholes, affine2d, run_bpfix
                ):
    """
    Short Summary
    -------------
    Applies the image plane algorithm to an AMI image

    Parameters
    ----------
    input_model : data model object
        AMI science image to be analyzed
    oversample : integer
        Oversampling factor
    rotation : float (degrees)
        Initial guess at rotation of science image relative to model
    psf_offset : string (two floats)
        PSF offset values to use to create the model array\
    rotsearch_parameters : string ('start stop step')
        Rotation search parameters
    src : string
        Source spectral type for model
    bandpass : synphot spectrum or array
        Synphot spectrum or array to override filter/source
    usebp : boolean
        If True, exclude pixels marked DO_NOT_USE from fringe fitting
    firstfew : integer
        If not None, process only the first few integrations
    chooseholes : string
        If not None, fit only certain fringes e.g. ['B4','B5','B6','C2']
    affine2d : user-defined Affine2D object
        None or user-defined Affine2d object
    run_bpfix : boolean
        Run Fourier bad pixel fix on cropped data

    Returns
    -------
    oifitsmodel: AmiOIModel object
        AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
    oifitsmodel_multi: AmiOIModel object
        AMI tables of observables for each integration from LG algorithm fringe fitting in OIFITS format
    amilgmodel: AmiLGFitModel object
        AMI cropped data, model, and residual data from LG algorithm fringe fitting

    """
    # Create copy of input_model to avoid overwriting input
    input_copy = copy.deepcopy(input_model)

    # If the input image is 2D, expand all relevant extensions to be 3D
    # Incl. those not currently used?
    if len(input_model.data.shape) == 2:
        input_copy.data = np.expand_dims(input_copy.data, axis=0)
        input_copy.dq = np.expand_dims(input_copy.dq, axis=0)
        # input_copy.err = np.expand_dims(input_copy.err, axis=0)
        # input_copy.var_poisson = np.expand_dims(input_copy.var_poisson, axis=0)
        # input_copy.var_rnoise = np.expand_dims(input_copy.var_rnoise, axis=0)
        # input_copy.var_flat = np.expand_dims(input_copy.var_flat, axis=0)



    # If the input data were taken in full-frame mode, extract a region
    # equivalent to the SUB80 subarray mode to make execution time acceptable.
    if input_model.meta.subarray.name.upper() == 'FULL':
        log.info("Extracting 80x80 subarray from full-frame data")
        xstart = 1045
        ystart = 1
        xsize = 80
        ysize = 80
        xstop = xstart + xsize - 1
        ystop = ystart + ysize - 1
        input_copy.data = input_copy.data[:, ystart - 1:ystop, xstart - 1:xstop]
        input_copy.dq = input_copy.dq[:, ystart - 1:ystop, xstart - 1:xstop]
        input_copy.err = input_copy.err[:, ystart - 1:ystop, xstart - 1:xstop]

    data = input_copy.data
    dim = data.shape[-1] # 80 px 

    # Initialize transformation parameters:
    #   mx, my: dimensionless magnifications
    #   sx, sy: dimensionless shears
    #   x0, y0: offsets in pupil space
    mx = 1.0
    my = 1.0
    sx = 0.0
    sy = 0.0
    xo = 0.0
    yo = 0.0

    psf_offset_ff = None
    # get filter, pixel scale from input_model,
    # make bandpass array for find_rotation, instrument_data calls
    filt = input_copy.meta.instrument.filter
    pscaledegx, pscaledegy = utils.degrees_per_pixel(input_copy)
    # model requires single pixel scale, so average X and Y scales
    # (this is done again in instrument_data?)
    pscale_deg = np.mean([pscaledegx, pscaledegy])
    PIXELSCALE_r = np.deg2rad(pscale_deg)
    holeshape = 'hex'
    # Throughput (combined filter and source spectrum) calculated here
    if bandpass is not None:
        log.info('User-defined bandpass provided: OVERWRITING ALL NIRISS-SPECIFIC FILTER/BANDPASS VARIABLES')
        # bandpass can be user-defined synphot object or appropriate array
        if isinstance(bandpass, synphot.spectrum.SpectralElement):
            log.info('User-defined synphot spectrum provided')
            wl, wt = bandpass._get_arrays(bandpass.waveset)
            bandpass = np.array((wt,wl)).T
        else:
            log.info('User-defined bandpass array provided')
            bandpass = np.array(bandpass)

    else:
        # get the filter and source spectrum
        log.info(f'Getting WebbPSF throughput data for {filt}.')
        filt_spec = utils.get_filt_spec(filt)
        log.info(f'Getting source spectrum for spectral type {src}.')
        src_spec = utils.get_src_spec(src) # always going to be A0V currently
        nspecbin = 19 # how many wavelngth bins used across bandpass -- affects runtime
        # **NOTE**: As of WebbPSF version 1.0.0 filter is trimmed to where throughput is 10% of peak
        # For consistency with WebbPSF simultions, use trim=0.1
        bandpass = utils.combine_src_filt(filt_spec, 
                                      src_spec, 
                                      trim=0.01, 
                                      nlambda=nspecbin,
                                      verbose=False, 
                                      plot=False) 
            

    rotsearch_d = np.append(np.arange(rotsearch_parameters[0], rotsearch_parameters[1], rotsearch_parameters[2]),
                            rotsearch_parameters[1])

    log.info(f'Initial values to use for rotation search: {rotsearch_d}')
    if affine2d is None:
        # affine2d object, can be overridden by user input affine.
        # do rotation search on uncropped median image (assuming rotation constant over exposure)
        # replace remaining NaNs in median image with median of surrounding 8 (non-NaN) pixels
        # (only used for centroiding in rotation search)
        meddata = np.median(data,axis=0)
        nan_locations = np.where(np.isnan(meddata))
        log.info(f'Replacing {len(nan_locations[0])} NaNs in median image with median of surrounding pixels')
        box_size = 3
        hbox = int(box_size / 2)
        for i_pos in range(len(nan_locations[0])):
            y_box = nan_locations[0][i_pos]
            x_box = nan_locations[1][i_pos]
            box = meddata[y_box - hbox:y_box + hbox + 1, x_box - hbox: x_box + hbox + 1]
            median_fill = np.nanmedian(box)
            if np.isnan(median_fill):
                median_fill = 0 # not ideal
            meddata[y_box, x_box] = median_fill

        affine2d = find_rotation(meddata, psf_offset, rotsearch_d,
                                 mx, my, sx, sy, xo, yo,
                                 PIXELSCALE_r, dim, bandpass, oversample, holeshape)

    niriss = instrument_data.NIRISS(filt, 
                                    bandpass=bandpass,
                                    affine2d=affine2d,
                                    src=src,
                                    firstfew=firstfew,
                                    usebp=usebp,
                                    chooseholes=chooseholes,
                                    run_bpfix=run_bpfix)

    ff_t = nrm_core.FringeFitter(niriss, 
                                psf_offset_ff=psf_offset_ff,
                                oversample=oversample)

    oifitsmodel, oifitsmodel_multi, amilgmodel = ff_t.fit_fringes_all(input_copy)


    # Copy header keywords from input to outputs
    oifitsmodel.update(input_model, only="PRIMARY")
    oifitsmodel_multi.update(input_model, only="PRIMARY")
    amilgmodel.update(input_model, only="PRIMARY")
    
    return oifitsmodel, oifitsmodel_multi, amilgmodel
