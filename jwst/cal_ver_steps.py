""" Module for generating CAL_VER related version reports/reference file for all
JWST calibration processing steps.

"""
from __future__ import print_function

import os, sys, inspect
import imp
import json

from datetime import datetime as dtime

from verhawk.scanner import Scanner

from . import version
from . import steps  # REQUIRED: in order to load all packages for inspect

from . import __name__ as jwst_pkg_name
from . import __file__ as jwst_pkg_file

__version__ = "0.9.3"

STDOUT = sys.stdout

class StepVersions(object):
    def __init__(self, author, **pars):
        """
        Parameters
        ==========
        author : str
            Name of user generating this reference file [Required]

        descrip : str, optional
            Basic description of this generated report. Default string will be
            generated if no user-supplied string is provided.

        history : str, optional
            Single line to provide info on when this file was created. Default
            string based on current date task was run will be used if no
            user-supplied string is provided.

        verbose : bool, optional
            Specify whether or not to generate additional diagnostic
            messages during operation. [Default: False}

        """
        # Initialize output reference file information based on user input
        descrip = pars.get('pars', "JWST calibration processing step version reference file")
        history = pars.get('history', "Created by cal_ver_steps version {}".format(__version__))
        self.verbose = pars.get('verbose', False)

#        useafter = pars.get('useafter', None)
#        if useafter is None:
#        useafter = dtime.isoformat(dtime.today())#"%Y-%m-%dT%H:%M:%S"

        self.output = {'reftype': "CALVER",
                      'instrument': "SYSTEM",
                      'useafter': "1900-01-01T00:00:00",
                      'telescope': "jwst",
                      'pedigree': "dummy",
                      'descrip': descrip,
                      'author': author,
                      'history': history,
                      'CAL_VER': version.__version__,
                      'versions': {}
                  }
        self.versions = {}

        # Use verhawk to extract all the version information from the package
        self.pkg_name = jwst_pkg_name

        try:
            with open(os.devnull, 'w') as devnull:
                sys.stdout = devnull
            jwst_pkg = imp.load_source(jwst_pkg_name, jwst_pkg_file)
            sys.stdout = STDOUT

            pvers = Scanner(jwst_pkg, recursive=True)
            pvers.scan()
            self.pkg_versions = pvers.versions

        except ImportError as e:
            print(e, file=sys.stderr)
            exit(1)

        # Identify which sub-packages are defined in this environment as
        # included in the 'steps' module
        self.modules = inspect.getmembers(jwst_pkg, inspect.ismodule)

        # Look for all modules within each sub-package which MAY contain a
        # processing step
        self.steps = {}
        for m in self.modules:
            mod_name = "{}.{}".format(self.pkg_name, m[0])
            self.steps[mod_name] = inspect.getmembers(m[1], inspect.ismodule)

    def scan(self):
        for s in self.steps:
            for mod in self.steps[s]:
                step_name = mod[0]
                step_module = mod[1]
                classname = None
                step_version = self.pkg_versions[s]
                if '_step' in step_name:
                    # Return all classes defined by step
                    classnames = inspect.getmembers(step_module, inspect.isclass)

                    # Extract actual processing step associated with this module
                    for c in classnames:
                        if c[0].endswith('Step') and c[0] != 'Step':
                            # store results
                            self.versions.update({c[0]: step_version})
                            if self.verbose:
                                print("Found Step class {} in {}".format(c[0], s))

        # Use this new information to update output versions
        self.output['versions'].update(self.versions)

    def as_json(self):
        """ Return json string for output version information.
        """
        return json.dumps(self.output, indent=4, sort_keys=True)

    def write_json(self, filename, clobber=False):
        """ Writes results out to json file.

        Parameters
        ==========
        filename : str
            Filename for output file.  If it does not end in .json, method will
            append '.json' to end of filename automatically.

        clobber : bool
            Specify whether or not this method should over-write a previous
            output file. If True, it will delete previous file.
            If False, it will raise an IOError exception when trying to overwrite
            a file.

        """

        if not filename.endswith('.json'):
            filename += '.json'

        if clobber:
            try:
                os.remove(filename)
            except:
                pass
            if self.verbose:
                print('Removing previous version of {}'.format(filename))
        else:
            if os.path.exists(filename):
                raise IOError("previous version of {} found, clobber not set".format(filename))

        with open(filename, 'w') as outfile:
            outfile.write(self.as_json())

        if self.verbose:
            print("Versions written to {}".format(filename))
