# Copyright (C) 2012 Association of Universities for Research in Astronomy(AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
"""
A client library for CRDS
"""
import re
import gc

# ----------------------------------------------------------------------

import six

# ----------------------------------------------------------------------

import crds
from crds.core import log, config, exceptions, heavy_client
from crds.core import crds_cache_locking

# ----------------------------------------------------------------------

def _flatten_dict(nested):
    """Takes a hierarchical arbitrarily nested dictionary of dictionaries, much
    like a data model, and flattens it.  The end result has only non-dictionary
    values referred to by the dotted paths that describe the traversal down
    through the nested dictionaries to their values.
    """
    def flatten(root, path, output):
        """Private worker function for _flatten_dict()."""
        for key, val in root.items():
            if isinstance(key, six.string_types):
                if isinstance(val, dict):
                    flatten(val, path + [key], output)
                else:
                    output['.'.join(path + [key])] = val
    output = {}
    flatten(nested, [], output)
    return output

# ----------------------------------------------------------------------

def get_multiple_reference_paths(input_file, reference_file_types):
    """Aligns JWST pipeline requirements with CRDS library top
    level interfaces.

    get_multiple_reference_paths() layers these additional tasks onto
    crds.getreferences():

    It converts an input file into a flat dictionary of JWST data
    model dotted parameters for defining CRDS best references.

    Returns { filetype : filepath or "N/A", ... }
    """
    from .. import datamodels

    gc.collect()

    if not reference_file_types:   # [] interpreted in CRDS as *all types*, intercept to handle *no types*.
        return {}

    if isinstance(input_file, (six.string_types, datamodels.DataModel)):
        with datamodels.open(input_file) as dm:
            data_dict = dm.to_flat_dict(include_arrays=False)
    else:  # XXX not sure what this does... seems unneeded.
        data_dict = _flatten_dict(input_file)

    gc.collect()

    exc = None
    bestrefs = {}
    with crds_cache_locking.get_cache_lock():
        for reftype in reference_file_types:
            try:
                ref = crds.getreferences(data_dict, reftypes=[reftype], observatory="jwst")
                bestrefs.update(ref)
            except Exception as exc:
                log.error(str(exc))

    if exc is not None:
        raise exceptions.CrdsError("One or more reference file fetches failed,  review CRDS ERROR messages.")

    refpaths = {filetype: filepath if "N/A" not in filepath.upper() else "N/A"
                for (filetype, filepath) in bestrefs.items()}

    return refpaths


def check_reference_open(refpath):
    """Verify that `refpath` exists and is readable for the current user.

    Ignore reference path values of "N/A" or "" for checking.
    """
    if refpath != "N/A" and refpath.strip() != "":
        fd = open(refpath, "rb")
        fd.close()
    return refpath

def get_reference_file(input_file, reference_file_type):
    """
    Gets a reference file from CRDS as a readable file-like object.
    The actual file may be optionally overridden.

    Parameters
    ----------
    input_file : jwst.datamodels.ModelBase instance
        A model of the input file.  Metadata on this input file will
        be used by the CRDS "bestref" algorithm to obtain a reference
        file.

    reference_file_type : string
        The type of reference file to retrieve.  For example, to
        retrieve a flat field reference file, this would be 'flat'.

    Returns
    -------
    reference_filepath : string
        The path of the reference in the CRDS file cache.
    """
    return get_multiple_reference_paths(input_file, [reference_file_type])[reference_file_type]

def get_override_name(reference_file_type):
    """
    Returns the name of the override configuration parameter for the
    given reference file type.

    Parameters
    ----------
    reference_file_type : string
        A reference file type name

    Returns
    -------
    config_parameter_name : string
        The configuration parameter name to use to override the given
        reference file type.
    """
    if not re.match('^[_A-Za-z][_A-Za-z0-9]*$', reference_file_type):
        raise ValueError(
            "{0!r} is not a valid reference file type name. "
            "It must be an identifier".format(reference_file_type))
    return "override_{0}".format(reference_file_type)


def get_svn_version():
    """Return the CRDS s/w version used for determining best references."""
    return crds.__version__


def get_context_used():
    """Return the context (.pmap) used for determining best references."""
    _connected, final_context = heavy_client.get_processing_mode("jwst")
    return final_context

def reference_uri_to_cache_path(reference_uri):
    """Convert an abstract CRDS reference file URI into an absolute path for the
    file as located in the CRDS cache.

    e.g. 'crds://jwst_miri_flat_0177.fits'  -->  
            '/grp/crds/cache/references/jwst/jwst_miri_flat_0177.fits'

    Standard CRDS cache paths are typically defined relative to the CRDS_PATH
    environment variable.  See https://jwst-crds.stsci.edu guide and top level
    page for more info on configuring CRDS.

    The default CRDS_PATH value is /grp/crds/cache, currently on the Central Store.
    """
    if not reference_uri.startswith("crds://"):
        raise exceptions.CrdsError("CRDS reference URI's should start with 'crds://' but got", repr(reference_uri))
    basename = config.pop_crds_uri(reference_uri)
    return crds.locate_file(basename, "jwst")
