"""
Various utility functions and data types
"""

import sys
from os.path import basename

import numpy as np
from astropy.io import fits

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())

def open(init=None, extensions=None, **kwargs):
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

    extensions : list of AsdfExtension
        A list of extensions to the ASDF to support when reading
        and writing ASDF files.

    Returns
    -------
    model : DataModel instance
    """

    from . import model_base
    from . import filetype

    # Initialize variables used to select model class

    hdulist = {}
    shape = ()
    file_to_close = None

    # Get special cases for opening a model out of the way
    # all special cases return a model if they match

    if init is None:
        return model_base.DataModel(None)

    elif isinstance(init, model_base.DataModel):
        # Copy the object so it knows not to close here
        return init.__class__(init)

    elif isinstance(init, (str, bytes)) or hasattr(init, "read"):
        # If given a string, presume its a file path.
        # if it has a read method, assume a file descriptor

        if isinstance(init, bytes):
            init = init.decode(sys.getfilesystemencoding())

        file_type = filetype.check(init)

        if file_type == "fits":
            hdulist = fits.open(init)
            file_to_close = hdulist

        elif file_type == "asn":
            # Read the file as an association / model container
            from . import container
            return container.ModelContainer(init, extensions=extensions,
                                            **kwargs)

        elif file_type == "asdf":
            # Read the file as asdf, no need for a special class
            return model_base.DataModel(init, extensions=extensions,
                                        **kwargs)

    elif isinstance(init, tuple):
        for item in init:
            if not isinstance(item, int):
                raise ValueError("shape must be a tuple of ints")
        shape = init

    elif isinstance(init, np.ndarray):
        shape = init.shape

    elif isinstance(init, fits.HDUList):
        hdulist = init

    # If we have it, determine the shape from the science hdu
    if hdulist:
        # So we don't need to open the image twice
        init = hdulist

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
    if isinstance(init, str):
        log.debug('Opening {0} as {1}'.format(basename(init), new_class))
    else:
        log.debug('Opening as {0}'.format(new_class))

    # Actually open the model
    model = new_class(init, extensions=extensions, **kwargs)
    if not has_model_type:
        try:
            delattr(model.meta, 'model_type')
        except AttributeError:
            pass

    # Close the hdulist if we opened it
    if file_to_close is not None:
        model._files_to_close.append(file_to_close)

    return model


def _class_from_model_type(hdulist):
    """
    Get the model type from the primary header, lookup to get class
    """
    from . import _defined_models as defined_models

    if hdulist:
        primary = hdulist[0]
        model_type = primary.header.get('DATAMODL')

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
                # It's a RampModel or MIRIRampModel
                try:
                    hdulist['REFOUT']
                except KeyError:
                    # It's a RampModel
                    from . import ramp
                    new_class = ramp.RampModel
                else:
                    # It's a MIRIRampModel
                    from . import miri_ramp
                    new_class = miri_ramp.MIRIRampModel
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
        new_class = model_base.DataModel
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


def gentle_asarray(a, dtype):
    """
    Performs an asarray that doesn't cause a copy if the byteorder is
    different.  It also ignores column name differences -- the
    resulting array will have the column names from the given dtype.
    """
    out_dtype = np.dtype(dtype)
    if isinstance(a, np.ndarray):
        in_dtype = a.dtype
        # Non-table array
        if in_dtype.fields is None and out_dtype.fields is None:
            if np.can_cast(in_dtype, out_dtype, 'equiv'):
                return a
            else:
                return np.asanyarray(a, dtype=out_dtype)
        elif in_dtype.fields is not None and out_dtype.fields is not None:
            if in_dtype == out_dtype:
                return a
            if len(in_dtype) != len(out_dtype):
                raise ValueError(
                    "Wrong number of columns.  Expected {0}, got {1}".format(
                        len(out_dtype), len(in_dtype)))
            new_dtype = []
            # Change the dtype name to match the fits record names
            # as the mismatch causes case insensitive access to fail
            if hasattr(in_dtype, 'names') and hasattr(out_dtype, 'names'):
                out_dtype.names = in_dtype.names
            for i in range(len(out_dtype.fields)):
                in_type = in_dtype[i]
                out_type = out_dtype[i]
                if in_type.subdtype is None:
                    type_str = in_type.str
                else:
                    type_str = in_type.subdtype[0].str
                if np.can_cast(in_type, out_type, 'equiv'):
                    new_dtype.append(
                        (out_dtype.names[i],
                         type_str,
                         in_type.shape))
                else:
                    return np.asanyarray(a, dtype=out_dtype)
            return a.view(dtype=np.dtype(new_dtype))
        else:
            return np.asanyarray(a, dtype=out_dtype)
    else:
        try:
            a = np.asarray(a, dtype=out_dtype)
        except:
            raise ValueError("Can't convert {0!s} to ndarray".format(type(a)))
        return a

def get_short_doc(schema):
    title = schema.get('title', None)
    description = schema.get('description', None)
    if description is None:
        description = title or ''
    else:
        if title is not None:
            description = title + '\n\n' + description
    return description.partition('\n')[0]


def ensure_ascii(s):
    if isinstance(s, bytes):
        s = s.decode('ascii')
    return s


def create_history_entry(description, software=None):
    """
    Create a HistoryEntry object.

    Parameters
    ----------
    description : str
        Description of the change.
    software : dict or list of dict
        A description of the software used.  It should not include
        asdf itself, as that is automatically notated in the
        `asdf_library` entry.

        Each dict must have the following keys:

        ``name``: The name of the software
        ``author``: The author or institution that produced the software
        ``homepage``: A URI to the homepage of the software
        ``version``: The version of the software

    Examples
    --------
    >>> soft = {'name': 'jwreftools', 'author': 'STSCI',
                'homepage': 'https://github.com/spacetelescope/jwreftools', 'version': "0.7"}
    >>> entry = create_history_entry(description="HISTORY of this file", software=soft)

    """
    from asdf.tags.core import Software, HistoryEntry
    import datetime

    if isinstance(software, list):
            software = [Software(x) for x in software]
    elif software is not None:
        software = Software(software)

    entry = HistoryEntry({
        'description': description,
        'time': datetime.datetime.utcnow()
    })

    if software is not None:
        entry['software'] = software
    return entry
