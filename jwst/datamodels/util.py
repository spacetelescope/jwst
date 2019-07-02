"""
Various utility functions and data types
"""

import sys
import warnings
from os.path import basename

import numpy as np
from astropy.io import fits

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())

class NoTypeWarning(Warning):
    pass

def open(init=None, **kwargs):
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

    Returns
    -------
    model : DataModel instance
    """

    from . import model_base
    from . import filetype

    # Initialize variables used to select model class

    hdulist = {}
    shape = ()
    file_name = None
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

        file_name = basename(init)
        file_type = filetype.check(init)

        if file_type == "fits":
            hdulist = fits.open(init)
            file_to_close = hdulist

        elif file_type == "asn":
            # Read the file as an association / model container
            from . import container
            return container.ModelContainer(init, **kwargs)

        elif file_type == "asdf":
            # Read the file as asdf, no need for a special class
            return model_base.DataModel(init, **kwargs)

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
        log.debug('Opening {0} as {1}'.format(file_name, new_class))
    else:
        log.debug('Opening as {0}'.format(new_class))

    # Actually open the model
    model = new_class(init, **kwargs)

    # Close the hdulist if we opened it
    if file_to_close is not None:
        model._files_to_close.append(file_to_close)

    if not has_model_type:
        class_name = new_class.__name__.split('.')[-1]
        if file_name:
            errmsg = \
                "model_type not found. Opening {} as a {}".format(file_name, class_name)
            warnings.warn(errmsg, NoTypeWarning)

        try:
            delattr(model.meta, 'model_type')
        except AttributeError:
            pass

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


def is_association(asn_data):
    """
    Test if an object is an association by checking for required fields
    """
    if isinstance(asn_data, dict):
        if 'asn_id' in asn_data and 'asn_pool' in asn_data:
            return True
    return False


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
            # When a FITS file includes a pseudo-unsigned-int column, astropy will return
            # a FITS_rec with an incorrect table dtype.  The following code rebuilds
            # in_dtype from the individual fields, which are correctly labeled with an
            # unsigned int dtype.
            # We can remove this once the issue is resolved in astropy:
            # https://github.com/astropy/astropy/issues/8862
            if isinstance(a, fits.fitsrec.FITS_rec):
                new_in_dtype = []
                updated = False
                for field_name in in_dtype.fields:
                    table_dtype = in_dtype[field_name]
                    field_dtype = a.field(field_name).dtype
                    if np.issubdtype(table_dtype, np.signedinteger) and np.issubdtype(field_dtype, np.unsignedinteger):
                        new_in_dtype.append((field_name, field_dtype))
                        updated = True
                    else:
                        new_in_dtype.append((field_name, table_dtype))
                if updated:
                    in_dtype = np.dtype(new_in_dtype)

            if in_dtype == out_dtype:
                return a
            in_names = {n.lower() for n in in_dtype.names}
            out_names = {n.lower() for n in out_dtype.names}
            if in_names == out_names:
                # Change the dtype name to match the fits record names
                # as the mismatch causes case insensitive access to fail
                out_dtype.names = in_dtype.names
            else:
                raise ValueError(
                    "Column names don't match schema. "
                    "Schema has {0}. Data has {1}".format(
                        str(out_names.difference(in_names)),
                        str(in_names.difference(out_names))))

            new_dtype = []
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
        except Exception:
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
    >>> soft = {'name': 'jwreftools', 'author': 'STSCI', \
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
