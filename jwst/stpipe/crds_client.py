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

"""This module defines functions that connect the core CRDS package to the
JWST CAL code and STPIPE, tailoring it to provide results in the forms
required by STPIPE.

WARNING:  JWST CAL and CRDS have circular dependencies.  Do not use CRDS imports
directly in modules other than this crds_client so that dependency order and
general integration can be managed here.
"""

import re

import crds
from crds.core import config, heavy_client, log
from crds.core.exceptions import CrdsError
from crds.core import crds_cache_locking
from stdatamodels import s3_utils


__all__ = [
    "check_reference_open",
    "CrdsError",
    "get_multiple_reference_paths",
    "get_override_name",
    "get_reference_file",
    "get_svn_version",
    "reference_uri_to_cache_path",
    "get_context_used",
]


def get_multiple_reference_paths(parameters, reference_file_types, observatory):
    """Aligns JWST pipeline requirements with CRDS library top level interfaces.

    Parameters
    ----------
    parameters : dict
        Parameters used by CRDS to compute best references.  Can
        be obtained from DataModel.get_crds_parameters().

    reference_file_types : list of str
        List of reference file types.

    observatory : str
        CRDS observatory code.

    Returns
    -------
    dict
        Map of reference file type to path (or N/A)
    """
    if not isinstance(parameters, dict):
        raise TypeError("First argument must be a dict of parameters")

    log.set_log_time(True)
    refpaths = _get_refpaths(parameters, tuple(reference_file_types), observatory)
    return refpaths


def _get_refpaths(data_dict, reference_file_types, observatory):
    """Tailor the CRDS core library getreferences() call to the JWST CAL code by
    adding locking and truncating expected exceptions.   Also simplify 'NOT FOUND n/a' to
    'N/A'.  Re-interpret empty reference_file_types as "no types" instead of core
    library default of "all types."
    """
    if not reference_file_types:   # [] interpreted as *all types*.
        return {}
    with crds_cache_locking.get_cache_lock():
        bestrefs = crds.getreferences(
            data_dict, reftypes=reference_file_types, observatory=observatory)
    refpaths = {filetype: filepath if "N/A" not in filepath.upper() else "N/A"
                for (filetype, filepath) in bestrefs.items()}
    return refpaths


def check_reference_open(refpath):
    """Verify that `refpath` exists and is readable for the current user.

    Ignore reference path values of "N/A" or "" for checking.
    """
    if refpath != "N/A" and refpath.strip() != "":
        if s3_utils.is_s3_uri(refpath):
            if not s3_utils.object_exists(refpath):
                raise RuntimeError("S3 object does not exist: " + refpath)
        else:
            with open(refpath, "rb"):
                pass
    return refpath


def get_reference_file(parameters, reference_file_type, observatory):
    """
    Gets a reference file from CRDS as a readable file-like object.
    The actual file may be optionally overridden.

    Parameters
    ----------
    parameters : dict
        Parameters used by CRDS to compute best references.  Can
        be obtained from DataModel.get_crds_parameters().

    reference_file_type : string
        The type of reference file to retrieve.  For example, to
        retrieve a flat field reference file, this would be 'flat'.

    observatory: string
        telescope name used with CRDS,  e.g. 'jwst'.

    Returns
    -------
    reference_filepath : string
        The path of the reference in the CRDS file cache.


    See also get_multiple_reference_paths().
    """
    return get_multiple_reference_paths(
        parameters, [reference_file_type], observatory)[reference_file_type]


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


def reference_uri_to_cache_path(reference_uri, observatory):
    """Convert an abstract CRDS reference file URI into an absolute path for the
    file as located in the CRDS cache.

    Parameters
    ----------

    reference_uri :  string identifying an abstract CRDS reference file,
                      .e.g. 'crds://jwst_miri_flat_0177.fits'

    observatory : string naming the telescope/project for CRDS, e.g. 'jwst'
           `observatory` should be omitted, None,  or defined as a literal string value.
           Omission and the value None are interpreted as 'jwst'.

    Returns
    -------

    reference_path : absolute file path to reference_uri in the CRDS cache,
           .e.g. '/grp/crds/cache/references/jwst/jwst_miri_flat_0177.fits'

    Standard CRDS cache paths are typically defined relative to the CRDS_PATH
    environment variable.  See https://jwst-crds.stsci.edu guide and top level
    page for more info on configuring CRDS.

    The default CRDS_PATH value is /grp/crds/cache, currently on the Central Store.
    """
    if not reference_uri.startswith("crds://"):
        raise CrdsError(
            "CRDS reference URI's should start with 'crds://' but got", repr(reference_uri))
    observatory = (observatory or 'jwst').lower()
    basename = config.pop_crds_uri(reference_uri)
    return crds.locate_file(basename, observatory)


def get_context_used(observatory):
    """Return the context (.pmap) used for determining best references.

    observatory : string naming the telescope/project for CRDS, e.g. 'jwst'
           `observatory` should be omitted, None,  or defined as a literal string value.
           Omission and the value None are interpreted as 'jwst'.
    """
    _connected, final_context = heavy_client.get_processing_mode(observatory)
    return final_context
