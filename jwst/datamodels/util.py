"""
Various utility functions and data types
"""

import sys
import warnings
import os
from os.path import basename
from platform import system as platform_system
import psutil
import traceback
import logging

import asdf

import numpy as np
from astropy.io import fits
from stdatamodels import filetype
from stdatamodels import s3_utils

from ..lib.basic_utils import bytes2human


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())



class NoTypeWarning(Warning):
    pass

def open(init=None, memmap=False, **kwargs):
    """
    Creates a DataModel from a number of different types

    Parameters
    ----------
    init : shape tuple, file path, file object, astropy.io.fits.HDUList,
           numpy array, dict, None

        - None: A default data model with no shape

        - shape tuple: Initialize with empty data of the given shape

        - file path: Initialize from the given file (FITS , JSON or ASDF)

        - readable file object: Initialize from the given file object

        - astropy.io.fits.HDUList: Initialize from the given
          `~astropy.io.fits.HDUList`

        - A numpy array: A new model with the data array initialized
          to what was passed in.

        - dict: The object model tree for the data model

    memmap : bool
        Turn memmap of FITS file on or off.  (default: False).  Ignored for
        ASDF files.

    kwargs : dict
        Additional keyword arguments passed to lower level functions. These arguments
        are generally file format-specific. Arguments of note are:

        - FITS

           skip_fits_update - bool or None
              `True` to skip updating the ASDF tree from the FITS headers, if possible.
              If `None`, value will be taken from the environmental SKIP_FITS_UPDATE.
              Otherwise, the default value is `True`.

    Returns
    -------
    model : DataModel instance
    """

    from . import model_base

    # Initialize variables used to select model class

    hdulist = {}
    shape = ()
    file_name = None
    file_to_close = None

    # Get special cases for opening a model out of the way
    # all special cases return a model if they match

    if init is None:
        return model_base.JwstDataModel(None)

    elif isinstance(init, model_base.JwstDataModel):
        # Copy the object so it knows not to close here
        return init.__class__(init)

    elif isinstance(init, (str, bytes)) or hasattr(init, "read"):
        # If given a string, presume its a file path.
        # if it has a read method, assume a file descriptor

        if isinstance(init, bytes):
            init = init.decode(sys.getfilesystemencoding())

        file_name = basename(init)
        file_type = filetype.check(init)

        if file_type == "fits":
            if s3_utils.is_s3_uri(init):
                hdulist = fits.open(s3_utils.get_object(init))
            else:
                hdulist = fits.open(init, memmap=memmap)
            file_to_close = hdulist

        elif file_type == "asn":
            # Read the file as an association / model container
            from . import container
            return container.ModelContainer(init, **kwargs)

        elif file_type == "asdf":
            if s3_utils.is_s3_uri(init):
                asdffile = asdf.open(s3_utils.get_object(init), **kwargs)
            else:
                asdffile = asdf.open(init, **kwargs)

            # Detect model type, then get defined model, and call it.
            new_class = _class_from_model_type(asdffile)
            if new_class is None:
                # No model class found, so return generic DataModel.
                return model_base.JwstDataModel(asdffile, **kwargs)

            return new_class(asdffile)

    elif isinstance(init, tuple):
        for item in init:
            if not isinstance(item, int):
                raise ValueError("shape must be a tuple of ints")
        shape = init

    elif isinstance(init, np.ndarray):
        shape = init.shape

    elif isinstance(init, fits.HDUList):
        hdulist = init

    elif is_association(init) or isinstance(init, list):
        from . import container
        return container.ModelContainer(init, **kwargs)

    # If we have it, determine the shape from the science hdu
    if hdulist:
        # So we don't need to open the image twice
        init = hdulist
        info = init.fileinfo(0)
        if info is not None:
            file_name = info.get('filename')

        try:
            hdu = hdulist[('SCI', 1)]
        except (KeyError, NameError):
            shape = ()
        else:
            if hasattr(hdu, 'shape'):
                shape = hdu.shape
            else:
                shape = ()

    # First try to get the class name from the primary header
    new_class = _class_from_model_type(hdulist)
    has_model_type = new_class is not None

    # Special handling for ramp files for backwards compatibility
    if new_class is None:
        new_class = _class_from_ramp_type(hdulist, shape)

    # Or get the class from the reference file type and other header keywords
    if new_class is None:
        new_class = _class_from_reftype(hdulist, shape)

    # Or Get the class from the shape
    if new_class is None:
        new_class = _class_from_shape(hdulist, shape)

    # Throw an error if these attempts were unsuccessful
    if new_class is None:
        raise TypeError("Can't determine datamodel class from argument to open")

    # Log a message about how the model was opened
    if file_name:
        log.debug(f'Opening {file_name} as {new_class}')
    else:
        log.debug(f'Opening as {new_class}')

    # Actually open the model
    model = new_class(init, **kwargs)

    # Close the hdulist if we opened it
    if file_to_close is not None:
        model._files_to_close.append(file_to_close)

    if not has_model_type:
        class_name = new_class.__name__.split('.')[-1]
        if file_name:
            warnings.warn(f"model_type not found. Opening {file_name} as a {class_name}",
                NoTypeWarning)
        try:
            delattr(model.meta, 'model_type')
        except AttributeError:
            pass

    return model


def _class_from_model_type(init):
    """
    Get the model type from the primary header, lookup to get class

    Parameter
    ---------
    init: AsdfFile or HDUList

    Return
    ------
    new_class: str or None
    """
    from . import _defined_models as defined_models

    if init:
        if isinstance(init, fits.hdu.hdulist.HDUList):
            primary = init[0]
            model_type = primary.header.get('DATAMODL')
        elif isinstance(init, asdf.AsdfFile):
            try:
                model_type = init.tree['meta']['model_type']
            except KeyError:
                model_type = None

        if model_type is None:
            new_class = None
        else:
            new_class = defined_models.get(model_type)
    else:
        new_class = None

    return new_class


def _class_from_ramp_type(hdulist, shape):
    """
    Special check to see if file is ramp file
    """
    if not hdulist:
        new_class = None
    else:
        if len(shape) == 4:
            try:
                hdulist['DQ']
            except KeyError:
                from . import ramp
                new_class = ramp.RampModel
            else:
                new_class = None
        else:
            new_class = None

    return new_class


def _class_from_reftype(hdulist, shape):
    """
    Get the class name from the reftype and other header keywords
    """
    if not hdulist:
        new_class = None

    else:
        primary = hdulist[0]
        reftype = primary.header.get('REFTYPE')
        if reftype is None:
            new_class = None

        else:
            from . import reference
            if len(shape) == 0:
                new_class = reference.ReferenceFileModel
            elif len(shape) == 2:
                new_class = reference.ReferenceImageModel
            elif len(shape) == 3:
                new_class = reference.ReferenceCubeModel
            elif len(shape) == 4:
                new_class = reference.ReferenceQuadModel
            else:
                new_class = None

    return new_class


def _class_from_shape(hdulist, shape):
    """
    Get the class name from the shape
    """
    if len(shape) == 0:
        from . import model_base
        new_class = model_base.JwstDataModel
    elif len(shape) == 4:
        from . import quad
        new_class = quad.QuadModel
    elif len(shape) == 3:
        from . import cube
        new_class = cube.CubeModel
    elif len(shape) == 2:
        try:
            hdulist[('SCI', 2)]
        except (KeyError, NameError):
            # It's an ImageModel
            from . import image
            new_class = image.ImageModel
        else:
            # It's a MultiSlitModel
            from . import multislit
            new_class = multislit.MultiSlitModel
    else:
        new_class = None

    return new_class


def can_broadcast(a, b):
    """
    Given two shapes, returns True if they are broadcastable.
    """
    for i in range(1, min(len(a), len(b)) + 1):
        adim = a[-i]
        bdim = b[-i]

        if not (adim == 1 or bdim == 1 or adim == bdim):
            return False

    return True


def to_camelcase(token):
    return ''.join(x.capitalize() for x in token.split('_-'))


def is_association(asn_data):
    """
    Test if an object is an association by checking for required fields
    """
    if isinstance(asn_data, dict):
        if 'asn_id' in asn_data and 'asn_pool' in asn_data:
            return True
    return False


def check_memory_allocation(shape, allowed=None, model_type=None, include_swap=True):
    """Check if a DataModel can be instantiated

    Parameters
    ----------
    shape : tuple
        The desired shape of the model.

    allowed : number or None
        Fraction of memory allowed to be allocated.
        If None, the environmental variable `DMODEL_ALLOWED_MEMORY`
        is retrieved. If undefined, then no check is performed.
        `1.0` would be all available memory. `0.5` would be half available memory.

    model_type : DataModel or None
        The desired model to instantiate.
        If None, `open` will be used to guess at a model type depending on shape.

    include_swap : bool
        Include available swap in the calculation.

    Returns
    -------
    can_instantiate, required_memory : bool, number
        True if the model can be instantiated and the predicted memory footprint.
    """
    # Determine desired allowed amount.
    if allowed is None:
        allowed = os.environ.get('DMODEL_ALLOWED_MEMORY', None)
        if allowed is not None:
            allowed = float(allowed)

    # Create the unit shape
    unit_shape = (1,) * len(shape)

    # Create the unit model.
    if model_type:
        unit_model = model_type(unit_shape)
    else:
        unit_model = open(unit_shape)

    # Size of the primary array.
    primary_array_name = unit_model.get_primary_array_name()
    primary_array = getattr(unit_model, primary_array_name)
    size = primary_array.nbytes
    for dimension in shape:
        size *= dimension

    # Get available memory
    available = get_available_memory(include_swap=include_swap)
    log.debug(f'Model size {bytes2human(size)} available system memory {bytes2human(available)}')

    if size > available:
        log.warning(
            f'Model {model_type} shape {shape} requires {bytes2human(size)} which is more than'
            f' system available {bytes2human(available)}'
        )

    if allowed and size > (allowed * available):
        log.debug(
            f'Model size greater than allowed memory {bytes2human(allowed * available)}'
        )
        return False, size

    return True, size


def get_available_memory(include_swap=True):
    """Retrieve available memory

    Parameters
    ----------
    include_swap : bool
        Include available swap in the calculation.

    Returns
    -------
    available : number
        The amount available.
    """
    system = platform_system()

    # Apple MacOS
    log.debug(f'Running OS is "{system}"')
    if system in ['Darwin']:
        return get_available_memory_darwin(include_swap=include_swap)

    # Default to Linux-like:
    return get_available_memory_linux(include_swap=include_swap)


def get_available_memory_linux(include_swap=True):
    """Get memory for a Linux system

    Presume that the swap space as reported is accurate at the time of
    the query and that any subsequent allocation will be held the value.

    Parameters
    ----------
    include_swap : bool
        Include available swap in the calculation.

    Returns
    -------
    available : number
        The amount available.
    """
    vm_stats = psutil.virtual_memory()
    available = vm_stats.available
    if include_swap:
        swap = psutil.swap_memory()
        available += swap.total
    return available


def get_available_memory_darwin(include_swap=True):
    """Get the available memory on an Apple MacOS-like system

    For Darwin, swap space is dynamic and will attempt to use the whole of the
    boot partition.

    If the system has been configured to use swap from other sources besides
    the boot partition, that available space will not be included.

    Parameters
    ----------
    include_swap : bool
        Include available swap in the calculation.

    Returns
    -------
    available : number
        The amount available.
    """
    vm_stats = psutil.virtual_memory()
    available = vm_stats.available
    if include_swap:

        # Attempt to determine amount of free disk space on the boot partition.
        try:
            swap = psutil.disk_usage('/private/var/vm').free
        except FileNotFoundError as exception:
            log.warn('Cannot determine available swap space.'
                     f'Reason:\n'
                     f'{"".join(traceback.format_exception(exception))}')
            swap = 0
        available += swap

    return available
