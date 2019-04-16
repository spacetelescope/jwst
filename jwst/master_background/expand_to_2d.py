import logging

import numpy as np

from gwcs.wcstools import grid_from_bounding_box
from .. import datamodels
from .. assign_wcs import nirspec               # for NIRSpec IFU data
from .. datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def expand_to_2d(input, m_bkg_spec):
    """Expand a 1-D background to 2-D.

    Parameters
    ----------
    input : `~jwst.datamodels.DataModel`
        The input science data.

    m_bkg_spec : str or `~jwst.datamodels.DataModel`
        Either the name of a file containing a 1-D background spectrum,
        or a data model containing such a spectrum.

    Returns
    -------
    background : `~jwst.datamodels.DataModel`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D.
    """

    with datamodels.open(m_bkg_spec) as bkg:
        if hasattr(bkg, 'spec'):                # MultiSpecModel
            if len(bkg.spec) > 1:
                log.warning("The input 1-D spectrum contains multiple spectra")
            spec_table = bkg.spec[0].spec_table
        else:                                   # CombinedSpecModel
            spec_table = bkg.spec_table
        tab_wavelength = spec_table['wavelength'].copy()
        try:
            tab_npixels = spec_table['npixels'].copy()
        except KeyError:
            tab_npixels = np.ones_like(tab_wavelength)
        tab_background = spec_table['flux'] / tab_npixels

    # We're going to use np.interp, so tab_wavelength must be strictly
    # increasing.
    if not np.all(np.diff(tab_wavelength) > 0):
        index = np.argsort(tab_wavelength)
        tab_wavelength = tab_wavelength[index].copy()
        tab_background = tab_background[index].copy()

    # Handle associations, or input ModelContainers
    if isinstance(input, datamodels.ModelContainer):
        background = bkg_for_container(input, tab_wavelength, tab_background)

    else:
        background = create_bkg(input, tab_wavelength, tab_background)

    return background


def bkg_for_container(input, tab_wavelength, tab_background):
    """Create a 2-D background for a container object.

    Parameters
    ----------
    input : JWST association or `~jwst.datamodels.ModelContainer`
        The input science data.

    tab_wavelength : 1-D ndarray
        The wavelength column read from the 1-D background table.

    tab_background : 1-D ndarray
        The flux column read from the 1-D background table, divided by
        the npixels column.

    Returns
    -------
    background : `~jwst.datamodels.ModelContainer`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D.
    """

    background = datamodels.ModelContainer()
    for input_model in input:
        temp = create_bkg(input_model, tab_wavelength, tab_background)
        background.append(temp)

    return background


def create_bkg(input, tab_wavelength, tab_background):
    """Create a 2-D background.

    Parameters
    ----------
    input : `~jwst.datamodels.DataModel`
        The input science data.

    tab_wavelength : 1-D ndarray
        The wavelength column read from the 1-D background table.

    tab_background : 1-D ndarray
        The flux column read from the 1-D background table, divided by
        the npixels column.

    Returns
    -------
    background : `~jwst.datamodels.DataModel`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D.
    """

    # Handle individual NIRSpec FS, NIRSpec MOS
    if isinstance(input, datamodels.MultiSlitModel):
        background = bkg_for_multislit(input, tab_wavelength, tab_background)

    # Handle MIRI LRS
    elif isinstance(input, datamodels.ImageModel):
        background = bkg_for_image(input, tab_wavelength, tab_background)

    # Handle MIRI MRS and NIRSpec IFU
    elif isinstance(input, datamodels.IFUImageModel):
        background = bkg_for_ifu_image(input, tab_wavelength, tab_background)

    else:
        # Shouldn't get here.
        raise RuntimeError("Input type {} is not supported."
                           .format(type(input)))

    return background


def bkg_for_multislit(input, tab_wavelength, tab_background):
    """Create a 2-D background for a MultiSlitModel.

    Parameters
    ----------
    input : `~jwst.datamodels.MultiSlitModel`
        The input science data.

    tab_wavelength : 1-D ndarray
        The wavelength column read from the 1-D background table.

    tab_background : 1-D ndarray
        The flux column read from the 1-D background table, divided by
        the npixels column.

    Returns
    -------
    background : `~jwst.datamodels.MultiSlitModel`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D.
    """

    background = input.copy()
    min_wave = np.amin(tab_wavelength)
    max_wave = np.amax(tab_wavelength)

    for (k, slit) in enumerate(input.slits):
        wl_array = get_wavelengths(slit)
        if wl_array is None:
            raise RuntimeError("Can't determine wavelengths for {}"
                               .format(type(slit)))

        # Wherever the wavelength is NaN, the background flux should to be set
        # to 0.  We replace NaN elements in wl_array with -1, so that np.interp
        # will detect that those values are out of range (note the `left`
        # argument to np.interp) and set the output to 0.
        wl_array[np.isnan(wl_array)] = -1.
        # flag values outside of background wavelength table
        mask_limit = (wl_array > max_wave) | (wl_array < min_wave)
        wl_array[mask_limit] = -1
        # bkg_flux will be a 2-D array, because wl_array is 2-D.
        bkg_flux = np.interp(wl_array, tab_wavelength, tab_background,
                             left=0., right=0.)

        background.slits[k].data[:] = bkg_flux.copy()
        background.slits[k].dq[mask_limit] = np.bitwise_or(background.slits[k].dq[mask_limit],
                                                           dqflags.pixel['DO_NOT_USE'])

    return background


def bkg_for_image(input, tab_wavelength, tab_background):
    """Create a 2-D background for an ImageModel.

    Parameters
    ----------
    input : `~jwst.datamodels.ImageModel`
        The input science data.

    tab_wavelength : 1-D ndarray
        The wavelength column read from the 1-D background table.

    tab_background : 1-D ndarray
        The flux column read from the 1-D background table, divided by
        the npixels column.

    Returns
    -------
    background : `~jwst.datamodels.ImageModel`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D.
    """

    background = input.copy()
    min_wave = np.amin(tab_wavelength)
    max_wave = np.amax(tab_wavelength)
    wl_array = get_wavelengths(input)
    if wl_array is None:
        raise RuntimeError("Can't determine wavelengths for {}"
                           .format(type(input)))

    wl_array[np.isnan(wl_array)] = -1.
    # flag values outside of background wavelength table
    mask_limit = (wl_array > max_wave) | (wl_array < min_wave)
    wl_array[mask_limit] = -1
    # bkg_flux will be a 2-D array, because wl_array is 2-D.
    bkg_flux = np.interp(wl_array, tab_wavelength, tab_background,
                         left=0., right=0.)

    background.data[:] = bkg_flux.copy()
    background.dq[mask_limit] = np.bitwise_or(background.dq[mask_limit],
                                              dqflags.pixel['DO_NOT_USE'])

    return background


def bkg_for_ifu_image(input, tab_wavelength, tab_background):
    """Create a 2-D background for an IFUImageModel

    Parameters
    ----------
    input : `~jwst.datamodels.IFUImageModel`
        The input science data.

    tab_wavelength : 1-D ndarray
        The wavelength column read from the 1-D background table.

    tab_background : 1-D ndarray
        The flux column read from the 1-D background table, divided by
        the npixels column.

    Returns
    -------
    background : `~jwst.datamodels.IFUImageModel`
        A copy of `input` but with the data replaced by the background,
        "expanded" from 1-D to 2-D. The dq flags are set to DO_NOT_USE
        for the pixels outside the region provided in the X1D background
        wavelength table.

    """

    background = input.copy()
    background.data[:, :] = 0.
    min_wave = np.amin(tab_wavelength)
    max_wave = np.amax(tab_wavelength)

    if input.meta.instrument.name.upper() == "NIRSPEC":
        list_of_wcs = nirspec.nrs_ifu_wcs(input)
        for ifu_wcs in list_of_wcs:
            x, y = grid_from_bounding_box(ifu_wcs.bounding_box)
            wl_array = ifu_wcs(x, y)[2]
            wl_array[np.isnan(wl_array)] = -1.

            # mask wavelengths not covered by the master background
            mask_limit = (wl_array > max_wave) | (wl_array < min_wave)
            wl_array[mask_limit] = -1

            # mask_limit is indices into each WCS slice grid.  Need them in
            # full frame coordinates, so we use x and y, which are full-frame
            full_frame_ind = y[mask_limit].astype(int), x[mask_limit].astype(int)
            # TODO - add another DQ Flag something like NO_BACKGROUND when we have space in dqflags
            background.dq[full_frame_ind] = np.bitwise_or(background.dq[full_frame_ind],
                                                          dqflags.pixel['DO_NOT_USE'])

            bkg_flux = np.interp(wl_array, tab_wavelength, tab_background,
                                 left=0., right=0.)
            background.data[y.astype(int), x.astype(int)] = bkg_flux.copy()

    elif input.meta.instrument.name.upper() == "MIRI":
        shape = input.data.shape
        grid = np.indices(shape, dtype=np.float64)
        wl_array = input.meta.wcs(grid[1], grid[0])[2]
        # first remove the nans from wl_array and replace with -1
        mask = np.isnan(wl_array)
        wl_array[mask] = -1.
        # next look at the limits of the wavelength table
        mask_limit = (wl_array > max_wave) | (wl_array < min_wave)
        wl_array[mask_limit] = -1

        # TODO - add another DQ Flag something like NO_BACKGROUND when we have space in dqflags
        background.dq[mask_limit] = np.bitwise_or(background.dq[mask_limit],
                                                  dqflags.pixel['DO_NOT_USE'])
        bkg_flux = np.interp(wl_array, tab_wavelength, tab_background,
                             left=0., right=0.)
        background.data[:, :] = bkg_flux.copy()
    else:
        raise RuntimeError("Exposure type {} is not supported."
                           .format(input.meta.exposure.type))

    return background


def get_wavelengths(model):
    """Read or compute wavelengths.

    Parameters
    ----------
    model : `~jwst.datamodels.DataModel`
        The input science data.

    Returns
    -------
    wl_array : 2-D ndarray
        An array of wavelengths corresponding to the data in `model`.
    """

    if hasattr(model, "wavelength"):
        wl_array = model.wavelength.copy()
        got_wavelength = True                   # may be reset below
    else:
        wl_array = None
    if (wl_array is None or len(wl_array) == 0 or
        np.nanmin(wl_array) == 0. and np.nanmax(wl_array) == 0.):
            got_wavelength = False
            wl_array = None

    if hasattr(model.meta, "wcs") and not got_wavelength:
        wcs = model.meta.wcs
        shape = model.data.shape
        grid = np.indices(shape[-2:], dtype=np.float64)
        wl_array = wcs(grid[1], grid[0])[2]

    return wl_array
