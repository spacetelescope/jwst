#!/usr/bin/env python

# Copyright (C) 2018 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

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

import argparse

import re
import inspect
import numpy as np
from os import path as os_path
from asdf import ValidationError
from asdf.tags.core import ndarray

from stdatamodels import validate
from stdatamodels import properties
from stdatamodels import schema as mschema
from stdatamodels.jwst.datamodels import (ImageModel, IFUImageModel)

from ..pipeline import (GuiderPipeline, DarkPipeline, Ami3Pipeline,
                        Tso3Pipeline, Coron3Pipeline,
                        Detector1Pipeline, Image2Pipeline,
                        Image3Pipeline, Spec2Pipeline, Spec3Pipeline)

DEBUG = 0

comment = re.compile(r'#.*')
continued = re.compile(r'\\\s*$')
assignment = re.compile(r'=(?!=)')

single_quote = re.compile(r"'[^']+'")
double_quote = re.compile(r'"[^"]+"')
triple_quote = re.compile(r'""".*?"""', flags=re.DOTALL)

relative_import = re.compile(r'\s*from\s+\.+[\.\w]*\s+import')


def _main(instrument, mode, level):
    """
    Main procedure for make_header
    """
    cls = get_class(instrument, mode)
    im = cls()

    defaults = set_default_values(im, instrument, mode, level)
    pipeline = choose_pipeline(defaults)

    add_metadata(im, defaults, pipeline)
    fill_extensions(im, defaults)

    return im


def add_metadata(im, defaults, pipeline):
    """
    Add metadata values to match values used in pipeline
    """
    keywords = undefined_metadata(pipeline, 'meta')
    for keyword in keywords:
        if keyword == 'meta':
            pass
        elif keyword == 'meta.model_type':
            pass
        elif keyword == 'meta.wcs':
            set_value(im, defaults, 'meta.wcsinfo')
        elif not another_instrument(defaults, im, keyword):
            set_value(im, defaults, keyword)
            additional_metadata(im, defaults, keyword)


def additional_metadata(im, defaults, keyword):
    """
    Add additional metadata items from the same group
    """
    groups = []

    groups.append(set(('meta.subarray.name',
                       'meta.subarray.xstart',
                       'meta.subarray.ystart',
                       'meta.subarray.xsize',
                       'meta.subarray.ysize')))

    groups.append(set(('meta.instrument.name',
                       'meta.instrument.band',
                       'meta.instrument.channel',
                       'meta.instrument.module',
                       'meta.instrument.detector',
                       'meta.instrument.filter',
                       'meta.instrument.grating',
                       'meta.instrument.pupil',
                       'meta.exposure.readpatt')))

    groups.append(set(('meta.asn.pool_name',
                       'meta.asn.table_name')))

    groups.append(set(('meta.target.catalog_name',
                       'meta.target.ra',
                       'meta.target.dec')))

    groups.append(set(('meta.exposure.start_time',
                       'meta.exposure.mid_time',
                       'meta.exposure.end_time',
                       'meta.exposure.exposure_time')))

    groups.append(set(('meta.exposure.frame_divisor',
                       'meta.exposure.frame_time',
                       'meta.exposure.groupgap',
                       'meta.exposure.group_time',
                       'meta.exposure.integration_time',
                       'meta.exposure.nframes',
                       'meta.exposure.ngroups',
                       'meta.exposure.nints',
                       'meta.exposure.nresets_at_start',
                       'meta.exposure.nresets_between_ints')))

    for group in groups:
        if keyword in group:
            for additional_key in group:
                if not another_instrument(defaults, im, additional_key):
                    if additional_key != keyword:
                        set_value(im, defaults, additional_key)


def another_instrument(defaults, im, keyword):
    """
    Return True if metadata item is for
    another instrument or exposure type
    """
    modes = ((('NIRCAM',), 'NIRCAM'),
             (('NIRSPEC',), 'NIRSPEC'),
             (('NIRSPEC',), 'IRS2'),
             (('MIRI',), 'MIRI'),
             (('FGS',), 'FGS'),
             (('NIRISS',), 'NIRISS'),
             (('NRC_FOCUS',), 'WFSC'),
             (('NRS_IFU',), 'IFU'),
             (('NRS_MSASPEC',), 'MSA'),
             (('MIR_MRS',), 'MRS'),
             (('MIR_LRS-FIXEDSLIT', 'MIR_LRS-SLITLESS'), 'LRS'),
             )

    instrument = defaults['meta.instrument.name']
    exposure_type = defaults['meta.exposure.type']

    schema = im._schema
    for attr in keyword.split('.'):
        schema = properties._get_schema_for_property(schema, attr)
        if schema:
            title = schema.get('title', '').upper()
            for mode in modes:
                (mode_name, mode_value) = mode
                if not (instrument in mode_name or
                        exposure_type in mode_name):
                    if title.find(mode_value) >= 0:
                        return True
        else:
            break

    return False


def build_exposure_type(instrument, mode):
    """
    Build the exposure type from the instrument name and mode
    """
    abbreviations = {'FGS': 'FGS', 'MIRI': 'MIR', 'NIRISS': 'NIS',
                     'NIRCAM': 'NRC', 'NIRSPEC': 'NRS'}

    spec_mode = {'FGS': '', 'MIRI': 'MRS', 'NIRISS': 'WFSS',
                 'NIRCAM': 'WFSS', 'NIRSPEC': 'FIXEDSLIT'}

    abbrev = abbreviations[instrument]
    if mode[0:4] == 'IMAG':
        mode = 'IMAGE'
        mode = '_'.join((abbrev, mode))
    elif mode[0:4] == 'SPEC':
        mode = spec_mode[instrument]
        mode = '_'.join((abbrev, mode))

    return mode


def choose_pipeline(defaults):
    """
    Select the pipeline class from the exposure type and level
    """
    mode = {'FGS_ACQ1': 1, 'FGS_ACQ2': 1, 'FGS_DARK': 2,
            'FGS_FINEGUIDE': 1, 'FGS_FOCUS': 6, 'FGS_ID-IMAGE': 1,
            'FGS_ID-STACK': 1, 'FGS_IMAGE': 6, 'FGS_INTFLAT': 0,
            'FGS_SKYFLAT': 6, 'FGS_TRACK': 1, 'MIR_IMAGE': 6,
            'MIR_TACQ': 0, 'MIR_LYOT': 5, 'MIR_4QPM': 5,
            'MIR_LRS-FIXEDSLIT': 9, 'MIR_LRS-SLITLESS': 4, 'MIR_MRS': 9,
            'MIR_DARK': 2, 'MIR_FLAT-IMAGE': 6, 'MIR_FLATIMAGE': 6,
            'MIR_FLAT-MRS': 6, 'MIR_FLATMRS': 6, 'MIR_CORONCAL': 7,
            'NIS_AMI': 3, 'NIS_DARK': 2, 'NIS_FOCUS': 6, 'NIS_IMAGE': 6,
            'NIS_LAMP': 0, 'NIS_SOSS': 4, 'NIS_TACQ': 6,
            'NIS_TACONFIRM': 6, 'NIS_WFSS': 9, 'N/A': 6, 'ANY': 6,
            'NRC_IMAGE': 6, 'NRC_GRISM': 9, 'NRC_TACQ': 6,
            'NRC_CORON': 5, 'NRC_FOCUS': 6, 'NRC_DARK': 2, 'NRC_FLAT': 6,
            'NRC_LED': 0, 'NRC_TACONFIRM': 6, 'NRC_TSIMAGE': 4,
            'NRC_WFSS': 9, 'NRC_TSGRISM': 4, 'NRS_AUTOFLAT': 0,
            'NRS_AUTOWAVE': 0, 'NRS_BRIGHTOBJ': 4, 'NRS_CONFIRM': 0,
            'NRS_DARK': 2, 'NRS_FIXEDSLIT': 9, 'NRS_FOCUS': 0, 'NRS_IFU': 9,
            'NRS_IMAGE': 6, 'NRS_LAMP': 0, 'NRS_MIMF': 0, 'NRS_MSASPEC': 9,
            'NRS_MSATA': 0, 'NRS_TACONFIRM': 0, 'NRS_TACQ': 0, 'NRS_TASLIT': 9,
            'NRS_WATA': 0}

    pipelines = [None, GuiderPipeline, DarkPipeline,
                 Ami3Pipeline, Tso3Pipeline, Coron3Pipeline,
                 Detector1Pipeline, Image2Pipeline, Image3Pipeline,
                 Detector1Pipeline, Spec2Pipeline, Spec3Pipeline]

    exposure_type = defaults['meta.exposure.type']
    index = mode[exposure_type]
    if index == 0:
        raise ValueError("Pipeline not defined for mode " +
                         exposure_type)
    if index > 5:
        level = defaults['meta.pipeline_level']
        index += level - 1

    return pipelines[index]


def drop_prefix(string, name):
    """
    Drop the model name preceding the metadata keyword
    """
    parts = name.split('.')
    for i, part in enumerate(parts):
        if part == string:
            return '.'.join(parts[i:])
    return ''


def drop_suffix(name):
    """
    Drop a method name from the keyword name
    """
    parts = name.split('.')
    while parts:
        if parts[-1][0] != '_' and parts[-1][-1] != '(':
            break
        parts.pop()

    return '.'.join(parts)


def fill_extensions(im, defaults):
    """
    Set image extensions to their default values
    """
    def fill_hdu(subschema, path, combiner, shape, recurse):
        if 'fits_hdu' not in subschema:
            return

        dtype = subschema.get('datatype')
        if dtype is None:
            return

        ndim = subschema.get('ndim')
        if ndim and ndim != len(shape):
            return

        dtype = ndarray.asdf_datatype_to_numpy_dtype(dtype)

        keyword = '.'.join(path)
        default = subschema.get('default', 0.0)
        im[keyword] = np.full(shape, default, dtype=dtype)

    shape = get_shape(im, defaults)
    mschema.walk_schema(im._schema, fill_hdu, shape)


def get_class(instrument, mode):
    """
    Get the model class from the exposure type
    """
    instrument = instrument.upper()
    mode = mode.upper()
    exposure_type = build_exposure_type(instrument, mode)

    if exposure_type in ('NRS_IFU', 'MIR_MRS'):
        cls = IFUImageModel
    else:
        cls = ImageModel
    return cls


def get_package(obj):
    """
    Get the package name associated with an object
    """
    if inspect.isclass(obj):
        path = obj.__module__
    elif inspect.isfunction(obj):
        path = obj.__module__
    elif inspect.ismodule(obj):
        path = obj.__name__
    else:
        path = ''
    path = path.split('.')

    # Do not pop the module name if the source
    # comes from an `_init_.py` file.
    source_file = inspect.getfile(obj)
    if os_path.basename(source_file) != '__init__.py':
        path.pop()

    return '.'.join(path)


def get_shape(im, defaults):
    """
    Get the shape of the image
    """
    shape = im.shape
    if not shape or shape[0] == 0:
        xsize = defaults['meta.subarray.xsize']
        ysize = defaults['meta.subarray.ysize']
        shape = (xsize, ysize)

    return shape


def get_subschema(im, keyword):
    """
    Get the subschema for a keyword or None if not found
    """
    schema = im._schema
    for attr in keyword.split('.'):
        schema = properties._get_schema_for_property(schema, attr)
        if schema is None:
            break
    return schema


def import_objects(package, line):
    """
    Run an import command and save any objects generated by it
    """
    line = line.lstrip()
    local_variables = {}
    global_variables = {'__package__': package}
    try:
        exec(line, global_variables, local_variables)  # nosec
    except ImportError:
        pass

    objects = []
    for value in local_variables.values():
        if (inspect.isclass(value) or inspect.isfunction(value) or
                inspect.ismodule(value)):
            objects.append(value)
    return objects


def log_file(filename):
    """
    Log files checked to a log file
    """
    if DEBUG:
        path = filename.split('/')
        try:
            jpos = path.index('jwst')
        except ValueError:
            jpos = 0
        filename = '/'.join(path[jpos:])
        print("Searching " + '/'.join(path[jpos:]))


def log_results(title, strings):
    """
    Write intermediate results to a log file
    """
    if DEBUG:
        print("\n=== %s ===\n" % title)
        for s in strings:
            print(s)


def next_line(fd):
    """
    Read the next line from a file, including any continuations
    """
    qcount = 0
    pcount = 0
    more = True
    long_line = ''

    while more:
        more = False
        line = fd.readline()
        if line:
            qcount += line.count('"""')
            qcount = qcount % 2
            if qcount > 0:
                more = True

            else:
                line = unquote(line)
                if continued.search(line):
                    more = True

                pcount += line.count('(')
                pcount -= line.count(')')
                if pcount > 0:
                    more = True

        long_line += line

    return unquote(long_line)


def search_file(source_file, package, string, seen):
    """
    Search a file for lines containing a string
    """
    matches = []
    with open(source_file, 'r') as fd:
        while True:
            line = next_line(fd)
            if not line:
                break

            import_match = relative_import.match(line)
            if import_match:
                for obj in import_objects(package, line):
                    new_matches = search_source_lines(obj, string, seen)
                    matches.extend(new_matches)
            elif line.find(string) >= 0:
                matches.append(line)
    return matches


def search_source_lines(obj, string, seen):
    """
    Search all files contained in a pipeline for a string
    """
    matches = []
    package = get_package(obj)
    if package:
        source_file = inspect.getfile(obj)
        if source_file not in seen:
            seen.add(source_file)
            if source_file.endswith('.py'):
                log_file(source_file)
                new_matches = search_file(source_file, package,
                                          string, seen)
                matches.extend(new_matches)
    return matches


def set_default_values(im, instrument, mode, level):
    """
    Set the default values used to set header keywords
    """
    defaults = {}
    set_standard_defaults(im, defaults, instrument, mode, level)
    set_pipeline_defaults(im, defaults)
    set_wcs_defaults(im, defaults)
    return defaults


def set_pipeline_defaults(im, defaults):
    """
    Set header keyword defaults associated with a pipeline
    """
    generic_defaults = {'meta.telescope': 'JWST',
                        'meta.subarray.name': 'FULL',
                        'meta.subarray.xstart': 1,
                        'meta.subarray.ystart': 1,
                        'meta.target.catalog_name': 'SMC',
                        'meta.target.ra': 5.3196,
                        'meta.target.dec': -72.98605,
                        'meta.target.type': 'FIXED',
                        'meta.target.source_type': 'EXTENDED'
                        }

    nircam_defaults = {'meta.subarray.xsize': 2048,
                       'meta.subarray.ysize': 2048,
                       'meta.subarray.fastaxis': 1,
                       'meta.subarray.slowaxis': -2,
                       'meta.instrument.channel': 'SHORT',
                       'meta.instrument.module': 'B',
                       'meta.instrument.detector': 'NRCB1',
                       'meta.instrument.filter': 'F115W',
                       'meta.instrument.pupil': 'CLEAR'
                       }

    nirspec_defaults = {'meta.subarray.xsize': 2048,
                        'meta.subarray.ysize': 2048,
                        'meta.subarray.fastaxis': 1,
                        'meta.subarray.slowaxis': 2,
                        'meta.instrument.detector': 'NRS1',
                        'meta.instrument.filter': 'F100LP',
                        'meta.instrument.grating': 'MIRROR',
                        'meta.exposure.readpatt': 'NRSRAPID'
                        }

    miri_defaults = {'meta.subarray.xsize': 1032,
                     'meta.subarray.ysize': 1024,
                     'meta.subarray.fastaxis': 1,
                     'meta.subarray.slowaxis': 2,
                     'meta.subarray.name': 'FULL',
                     'meta.instrument.channel': '34',
                     'meta.instrument.band': 'SHORT',
                     'meta.instrument.detector': 'MIRIMAGE',
                     'meta.instrument.filter': 'N/A',
                     }

    fgs_defaults = {'meta.subarray.xsize': 2048,
                    'meta.subarray.ysize': 2048,
                    'meta.subarray.fastaxis': 2,
                    'meta.subarray.slowaxis': -1,
                    'meta.instrument.detector': 'GUIDER2'
                    }

    niriss_defaults = {'meta.subarray.xsize': 2048,
                       'meta.subarray.ysize': 2048,
                       'meta.subarray.fastaxis': -2,
                       'meta.subarray.slowaxis': -1,
                       'meta.instrument.detector': 'NIS',
                       'meta.instrument.filter': 'F380M'
                       }

    any_defaults = {}

    instrument_defaults = {'NIRCAM': nircam_defaults,
                           'NIRSPEC': nirspec_defaults,
                           'MIRI': miri_defaults,
                           'FGS': fgs_defaults,
                           'NIRISS': niriss_defaults,
                           'ANY': any_defaults
                           }

    mir_slitless_defaults = {'meta.subarray.name': 'SLITLESSPRISM',
                             'meta.subarray.xstart': 1,
                             'meta.subarray.ystart': 321,
                             'meta.subarray.xsize': 68,
                             'meta.subarray.ysize': 1024
                             }

    mir_fixedslit_defaults = {'meta.instrument.detector': 'MIRIFULONG'}

    mir_mrs_defaults = {'meta.instrument.detector': 'MIRIFULONG'}

    exposure_type_defaults = {'MIR_LRS-SLITLESS': mir_slitless_defaults,
                              'MIR_LRS-FIXEDSLIT': mir_fixedslit_defaults,
                              'MIR_MRS': mir_mrs_defaults
                              }

    defaults.update(generic_defaults)
    instrument = defaults['meta.instrument.name']
    defaults.update(instrument_defaults[instrument])
    exposure_type = defaults['meta.exposure.type']
    exposure_deltas = exposure_type_defaults.get(exposure_type)
    if exposure_deltas is not None:
        defaults.update(exposure_deltas)

    shape = im.shape
    if shape:
        if shape[0] != 0:
            defaults['meta.subarray.xsize'] = shape[0]
        if shape[1] != 0:
            defaults['meta.subarray.ysize'] = shape[1]


def set_standard_defaults(im, defaults, instrument, mode, level):
    """
    Set keyword defaults from instrument, mode, and level
    """
    instrument = instrument.upper()
    if validate_value(im, 'meta.instrument.name', instrument):
        defaults['meta.instrument.name'] = instrument
    else:
        raise ValueError('Unrecognized instrument name: ' + instrument)

    if defaults['meta.instrument.name'] == 'ANY':
        mode = 'ANY'
    else:
        mode = mode.upper()
    if validate_value(im, 'meta.exposure.type', mode):
        defaults['meta.exposure.type'] = mode
    else:
        value = build_exposure_type(instrument, mode)
        if not validate_value(im, 'meta.exposure.type', value):
            value = None

        if value is None:
            raise ValueError('Unrecognized mode: ' + mode)
        else:
            defaults['meta.exposure.type'] = value

    try:
        value = int(level)
    except ValueError:
        match = re.search(r'(\d+)', level)
        if match:
            value = int(match.group(1))
        else:
            value = 0
    if value < 1 or value > 3:
        raise ValueError('Unrecognized level: ' + level)
    else:
        defaults['meta.pipeline_level'] = value


def set_value(im, defaults, keyword):
    """
    Set a value for a header keyword from defaults and data type
    """
    schema = get_subschema(im, keyword)
    val = None
    if schema:
        if keyword in defaults:
            val = defaults[keyword]
        elif 'default' in schema:
            val = schema['default']
        elif 'enum' in schema:
            for choice in ('N/A', 'NONE', 'ANY'):
                try:
                    i = schema['enum'].index(choice)
                    break
                except ValueError:
                    i = 0
            val = schema['enum'][i]
        elif 'type' in schema:
            typ = schema['type']
            if typ == 'string':
                if keyword.startswith('meta.cal_step'):
                    val = 'COMPLETE'
                else:
                    val = ''
            elif typ == 'number':
                val = 0.0
            elif typ == 'integer':
                val = 1
            elif typ == 'boolean':
                val = True
            elif typ == 'object':
                for field in schema['properties']:
                    subkey = '.'.join((keyword, field))
                    set_value(im, defaults, subkey)
            else:
                val = None
        else:
            val = None

    if val is not None:
        im[keyword] = val


def set_wcs_defaults(im, defaults):
    """
    Set default values for the wcs header keywords
    """
    wcs_defaults = {
        'crval1': 5.3196, 'crval2': -72.98605,
        'ctype1': 'RA---TAN', 'ctype2': 'DEC--TAN',
        'cunit1': 'deg', 'cunit2': 'deg',
        'cdelt1': 5.5E-06, 'cdelt2': 5.5E-06}

    shape = get_shape(im, defaults)
    for i in (1, 2):
        name = "crpix%d" % i
        wcs_defaults[name] = 0.5 * shape[i - 1]

        for j in (1, 2):
            name = "pc%d_%d" % (i, j)
            if i == j:
                wcs_defaults[name] = 1.0
            else:
                wcs_defaults[name] = 0.0

    for name in wcs_defaults:
        wcs_name = 'meta.wcsinfo.' + name
        defaults[wcs_name] = wcs_defaults[name]


def undefined_metadata(module, string):
    """
    Find the metadata keywords that are undefined in the pipeline
    """
    name_set = []
    for i in range(2):
        name_set.append(set())

    pat = re.compile(r'([\w\.]+' + string + r'[\w\.]+\(?)')
    seen = set()
    matched = search_source_lines(module,
                                  '.' + string + '.',
                                  seen)
    log_results("matched lines", matched)

    for line in matched:
        sides = assignment.split(line)
        if len(sides) == 1:
            sides.append(sides[0][:])
            sides[0] = ''

        j = 1
        for i in range(len(sides), 0, -1):
            pat_match = pat.search(sides[i - 1])
            if pat_match is not None:
                for name in pat_match.groups():
                    if name is not None:
                        name = drop_suffix(name)
                        name_set[j].add(name)
            j = 0

    log_results("lhs metadata", name_set[0])
    log_results("rhs metadata", name_set[1])

    undefined = set()
    for name in name_set[1]:
        name = drop_prefix(string, name)
        if name:
            undefined.add(name)

    log_results("undefined metadata", undefined)
    return undefined


def unquote(line):
    """
    Remove quoted strings and comments from a line
    """
    line = triple_quote.sub('', line)
    line = double_quote.sub('""', line)
    line = single_quote.sub("''", line)
    line = comment.sub('', line)
    return line


def validate_value(im, name, value):
    """
    Test a value against the schema to see if it is valid
    """
    valid = True
    schema = get_subschema(im, name)
    if schema:
        try:
            validate._check_value(value, schema, im)
        except ValidationError:
            valid = False
    return valid


def main():
    long_description = """
    Create an image the contains the header keywords needed to run
    the JWST pipeline for a specified instrument, observing mode, and
    pipeline level. These three are specified by the first three arguments
    to the command. The fourth and final argument is the filename that the
    image will be saved as.
    """

    parser = argparse.ArgumentParser(description=long_description)
    parser.add_argument('instrument', help='The instrument name')
    parser.add_argument('mode', help='\'image\' or \'spectrum\' or an exposure type')
    parser.add_argument('level', help='The pipeline level: 1, 2, or 3')
    parser.add_argument('filename', help='The output image with the headers')
    args = parser.parse_args()

    im = _main(args.instrument, args.mode, args.level)
    im.save(args.filename)
    im.close()


if __name__ == '__main__':
    main()
