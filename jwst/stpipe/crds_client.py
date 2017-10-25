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

# import memory_profiler
# @memory_profiler.profile
def get_multiple_reference_paths(input_file, reference_file_types):
    """Aligns JWST pipeline requirements with CRDS library top
    level interfaces.

    `input_file` can be a simple data model, a filename, or a nested meta data dict.
    `reference_file_types` is a list of reference type names to fetch bestrefs for.

    Returns best references dict { filetype : filepath or "N/A", ... }
    """
    # gc.collect()
    
    data_dict = _get_data_dict(input_file)

    refpaths = _get_refpaths(data_dict, tuple(reference_file_types))

    # gc.collect()
    
    return refpaths

# ......................
    
def _input_name(input_file):
    """
    Return input_file.meta.filename IFF input_file is a DataModel and it works.
    Return input_file IF input_file is a string,  assumed to be a filename.
    Return input_file otherwise.
    """
    from .. import datamodels
    if isinstance(input_file, six.string_types):
        return input_file
    elif isinstance(input_file, datamodels.DataModel):
        try:
            return input_file.meta.filename
        except AttributeError:
            return input_file
    else:
        return input_file

# ......................
    
_HEADER_CACHE = {}

def _get_data_dict(input_file):
    """Return the data models header dictionary based on input_file.

    `input_file` can be a filepath, an open data model, or a nested
    data-model-like metadata dictionary.

    Returns a flat parameter dictionary used for CRDS bestrefs matching.
    """
    from .. import datamodels
    name = _input_name(input_file)
    if isinstance(name, six.string_types):
        if name in _HEADER_CACHE:
            log.verbose("Using cached CRDS matching header for", repr(name))
        else:
            log.verbose("Caching CRDS matching header for", repr(name))
            with datamodels.open(input_file) as model:
                _HEADER_CACHE[name] = model.to_flat_dict(include_arrays=False)
    else:  # XXX not sure what this does... seems unneeded.
        _HEADER_CACHE[name] = _flatten_dict(input_file)
    return _HEADER_CACHE[name]

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

# ......................

def _get_refpaths(data_dict, reference_file_types):
    """Tailor the CRDS core library getreferences() call to the JWST CAL code by
    adding locking and truncating expected exceptions.
    """
    if not reference_file_types:   # [] interpreted as *all types*.
        return {}
    try:
        with crds_cache_locking.get_cache_lock():
            bestrefs = crds.getreferences(data_dict, reftypes=reference_file_types, observatory="jwst")
    except crds.CrdsBadRulesError as exc:
        raise crds.CrdsBadRulesError(str(exc))
    except crds.CrdsBadReferenceError as exc:
        raise crds.CrdsBadReferenceError(str(exc))    
    refpaths = {filetype: filepath if "N/A" not in filepath.upper() else "N/A"
                for (filetype, filepath) in bestrefs.items()}
    return refpaths

# ----------------------------------------------------------------------

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

def get_context_used():
    """Return the context (.pmap) used for determining best references."""
    _connected, final_context = heavy_client.get_processing_mode("jwst")
    return final_context

def init_multiprocessing(common_dataset_filepaths=()):
    """Perform CRDS cached operations first in a multipocessing root process so that every 
    subprocess inherits cached CRDS information as a consequence of UNIX forking semantics
    and doesn't repeat expensive operations.

    `common_dataset_filepaths` is a list of filepaths used to define reference matching
    parameters used by more than one process.  It can be empty.
    """
    # Determine the context to be used based on CRDS cache, CRDS server, and CRDS_CONTEXT.
    with log.warn_on_exception("Failed determining context name"):
        final_context = get_context_used()

    # Perform the load of `final_context`, in this case for all instruments, questionable
    # since a header-based load will only mappings for load the required instrument.
    # Should be around 2-12 seconds depending on file system performance.
    with log.warn_on_exception("Failed loading context:", repr(final_context)):
        mapping = crds.get_symbolic_mapping(final_context)
        mapping.force_load()
    
    # Perform header reads for all dataset_filepaths shared between processes
    for filepath in common_dataset_filepaths:
        with log.warn_on_exception("Failed reading file header for", repr(filepath)):
            _get_data_dict(filepath)

    # probably good to prevent inherited garbage,  do last.
    gc.collect()
