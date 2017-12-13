""" blendmeta - Merge metadata from multiple models to create
                   a new metadata instance and table

"""
from __future__ import absolute_import, print_function
import os
import glob
import copy

import numpy as np
from astropy.io import fits

from . import blender
from . import textutil

# Version of rules file format supported by this version of the code
# All changes should be backwards compatible to older rules versions
# so any rules file with Version >= __rules_version__ should work
# with this code
__rules_version__ = 2.1


# Custom blending functions
def multi(vals):
    """
    This will either return the common value from a list of identical values
    or 'MULTIPLE'
    """
    uniq_vals = list(set(vals))
    num_vals = len(uniq_vals)
    if num_vals == 0:
        return None
    if num_vals == 1:
        return uniq_vals[0]
    if num_vals > 1:
        return "MULTIPLE"


def multi1(vals):
    """
    This will either return the common value from a list of identical values
    or the single character '?'
    """
    uniq_vals = list(set(vals))
    num_vals = len(uniq_vals)
    if num_vals == 0:
        return None
    if num_vals == 1:
        return uniq_vals[0]
    if num_vals > 1:
        return "?"


def float_one(vals):
    """ Return a constant floating point value of 1.0
    """
    return 1.0


def int_one(vals):
    """ Return an integer value of 1
    """
    return int(1)


def zero(vals):
    """ Return a value of 0
    """
    return 0


def first(items):
    """ Return first item from list of values"""
    if len(items):
        return items[0]
    return None


def last(items):
    """ Return last item from list of values"""
    if len(items):
        return items[-1]
    return None


# translation dictionary for function entries from rules files
blender_funcs = {'first': first,
                 'last': last,
                 'float_one': float_one,
                 'int_one': int_one,
                 'zero': zero,
                 'multi': multi,
                 'multi?': multi1,
                 'mean': np.mean,
                 'sum': np.sum,
                 'max': np.max,
                 'min': np.min,
                 'stddev': np.std}

delete_command = '<delete>'


# Classes for managing keyword rules
class KeywordRules(object):

    rules_name_suffix = '_header.rules'

    def __init__(self, instrument, telescope='jwst', rules_file=None):
        """ Read in the rules used to interpret the keywords from the specified
            instrument image header.
        """
        self.instrument = instrument.lower()
        self.telescope = telescope.lower()
        self.new_header = None
        self.rules_version = None

        # Add support for user-specified rules file...
        self.rules_file = rules_file
        if self.rules_file is None:
            self.get_filename()  # define rules file

        self.rules_version, i = self.get_rules_header(self.rules_file)
        rfile = open(self.rules_file)
        self.rule_specs = rfile.readlines()
        rfile.close()

        self.rule_objects = []
        self.rules = []
        self.section_names = []
        self.delete_kws = []

    def get_filename(self):
        """ Return name of rules file to be used
        It will use a local copy if present, and use the installed version
        by default.
        Any local copy will take precendence over the default rules.

        This function will return the alphabetically first file that applies
        to the instrument and meets the version requirements
        """
        rules_file = None
        # get all potential local rules
        rfiles = glob.glob('*.rules')
        rfiles.sort()

        # Sort through list and find only applicable rules files
        # This would include picking up any rules files using the default
        # naming convention; namely, <instrument>_header.rules
        for r in rfiles:
            v, i = self.get_rules_header(r)
            if v is None or i is None:
                continue
            if v <= __rules_version__ and i == self.instrument.lower():
                rules_file = r
                break

        # Pipeline use should involve passing in rules_file read in
        # from CRDS.  If not, look for local user-defined file.
        if rules_file is None:
            # define default rules name installed with the software
            rules_name = self.instrument.lower()+self.rules_name_suffix
            rules_file = os.path.join(os.path.dirname(__file__), rules_name)
            if not os.path.exists(rules_file):
                rules_name = self.telescope+self.rules_name_suffix
                rules_file = os.path.join(os.path.dirname(__file__),
                                          rules_name)
                if not os.path.exists(rules_file):
                    rules_file = None

        if rules_file is None:
            errmsg = 'ERROR:\n'+'    No valid rules file found for:\n'
            errmsg += '    INSTRUMENT = %s\n' % (self.instrument)
            errmsg += '    RULES Version <= %s\n' % (__rules_version__)
            print(textutil.textbox(errmsg))
            raise ValueError

        self.rules_file = rules_file
        return rules_file

    def get_rules_header(self, filename):
        """
        Open a potential rules file and return the recognized
        version and instrument types provided in the file's first 2 lines
        """
        version = None
        instrument = None
        f = open(filename)  # open file in read-only mode
        for line in f.readlines():
            if line[0] == '!':
                if 'version' in line.lower():
                    version = float(line.strip('\n').split('=')[-1])
                if 'instrument' in line.lower():
                    instrument = line.strip('\n').split('=')[-1]
                    instrument = instrument.lower().strip()
        f.close()

        if version is not None and instrument is None:
            # try to extract instrument name from rules filename, if it
            # follows the default naming convention
            if 'header.rules' in filename.lower():
                inst = filename.split('_header.rules')[0].lower()
                if inst == self.instrument:
                    instrument = inst
        return version, instrument

    def interpret_rules(self, hdrs):
        """ Convert specifications for rules from rules file
            into specific rules for this header(instrument/detector)

            This allows for expansion rules to be applied to rules
            from the rules files (such as any wildcards or section titles).
        """
        if isinstance(hdrs, tuple):
            hdrs = list(hdrs)
        if not isinstance(hdrs, list):
            hdrs = [hdrs]

        # apply rules to headers
        for rule in self.rule_specs:
            if rule[0] in ['#', ' ', None, "None", "INDEF"]:
                continue
            kwr = KwRule(rule)
            duplicate_rule = False
            for robj in self.rule_objects:
                if kwr.rule_spec == robj.rule_spec:
                    duplicate_rule = True
                    break
            if not duplicate_rule:
                for hdr in hdrs:
                    kwr.interpret(hdr)
                self.rule_objects.append(kwr)

        for kwr in self.rule_objects:
            self.rules.extend(kwr.rules)
            self.delete_kws.extend(kwr.delete_kws)
            self.section_names.extend(kwr.section_name)

    def merge(self, kwrules):
        """
        Merge a new set of interpreted rules into the current set
        The new rules, kwrules, can either be a new class or a whole new
        set of rules (like those obtained from using self.interpret_rules with
        a new header).
        """
        if isinstance(kwrules, KeywordRules):
            kwrules = kwrules.rules

        # Determine what rules are specified in kwrules that
        #    are NOT in self.rules
        k = []
        # Delete these extraneous rules from input kwrules
        for r in kwrules:
            if r not in self.rules:
                k.append(r)

        # extend self.rules with additional rules
        self.rules.extend(k)

    def apply(self, headers, tabhdu=False):
        """ For a full list of metadata objects, apply the specified rules to
            generate a dictionary of new values and a table using
            blender.

            This method returns the new metadata object and summary table
            as `datamodels.model.ndmodel` and fits.binTableHDU objects.
        """
        # Apply rules to headers
        fbdict, fbtab = blender.fitsblender(headers, self.rules)

        # Determine which keywords are included in the table but not
        # the new dict(header). These will be removed from the output
        # header altogether
        tabcols = fbtab.dtype.names
        hdrkws = list(fbdict.keys())
        del_kws = list(set(tabcols) - set(hdrkws))

        # Start with a copy of the template as the new header
        # This will define what keywords need to be updated, as the rules
        # and input headers often include headers for multiple extensions in
        # order to build the complete table for all the keywords in the file
        # in one run
        new_header = copy.deepcopy(headers[0])

        # Delete all keywords from copy that are being moved into the table
        # However, this should only be done for those keywords which do are not
        # being kept in the header through fbdict (additional rules)
        for kw in del_kws:
            if (kw in new_header):
                try:
                    del new_header[kw]
                except KeyError:
                    pass

        # Remove section names from output header(s)
        for name in self.section_names:
            for indx, kw in zip(list(range(len(new_header), 0, -1)),
                                new_header[-1::-1]):
                if name in str(kw.value):
                    del new_header[indx-1]
                continue

        # Delete keywords marked in rules file
        for kw in self.delete_kws:
            if kw in new_header:
                try:
                    del new_header[kw]
                except KeyError:
                    pass

        # Apply updated/blended values into new header, but only those
        # keywords which are already present in the 'template' new header
        # this allows the rules to be used on all extensions at once yet
        # update each extension separately without making copies of kws from
        # one extension to another.
        for kw in fbdict:
            new_header[kw] = fbdict[kw]
        # Create summary table
        if len(tabcols) > 0:
            if tabhdu:
                new_table = fits.BinTableHDU.from_columns(fbtab)
                new_table.header['EXTNAME'] = 'HDRTAB'
            else:
                new_table = fbtab
        else:
            new_table = None
        return new_header, new_table

    def add_rules_kws(self, hdr):
        """
        Update PRIMARY header with HISTORY cards that report the exact
        rules used to create this header. Only non-comment lines from the
        rules file will be reported.
        """
        hdr['RULESVER'] = (self.rules_version,
                           'Version ID for header kw rules file')
        hdr['BLENDVER'] = (__version__, 'Version of blendheader software used')
        hdr['RULEFILE'] = (self.rules_file, 'Name of header kw rules file')
        hdr.add_history('='*60)
        hdr.add_history('Header Generation rules:')
        hdr.add_history('    Rules used to combine headers of input files')
        hdr.add_history('    Start of rules...')
        hdr.add_history('-'*60)
        for rule in self.rule_specs:
            if rule[0] in ['#', ' ', None, "None", "INDEF"]:
                continue
            hdr.add_history(rule.strip('\n'))

        hdr.add_history('-'*60)
        hdr.add_history('    End of rules...')
        hdr.add_history('='*60)

    def index_of(self, kw):
        """ Reports the index of the specified kw
        """
        indx = []
        for r, i in zip(self.rules, list(range(len(self.rules)))):
            if r[0] == kw:
                indx.append(i)
        return indx


class KwRule(object):
    """
    This class encapsulates the logic needed for interpreting a single keyword
    rule from a text file.

    The .rules attribute contains the interpreted set of rules that corresponds
    to this line.

    Example:
    Interpreting rule from:
        FILENAME    FILENAME    first

    into rule: [('FILENAME', 'FILENAME',
                <function first at 0x7fe505db7668>, 'ignore')]
    and sname: None
    and delkws: []

    """
    def __init__(self, line):
        self.rule_spec = line  # line read in from rules file
        self.rules = []
        self.delete_kws = []
        self.section_name = []

    def interpret(self, hdr):
        if self.rules:
            # If self.rules has already been defined for this rule, do not try
            # to interpret it any further with additional headers
            return
        irules, sname, delkws = interpret_line(self.rule_spec, hdr)

        # keep track of any section name identified for this rule
        if sname:
            self.section_name.append(sname)

        # also keep track of what keywords should be deleted based on this rule
        if delkws:
            self.delete_kws = delkws

        # Now, interpret rule based on presence of kw in hdr
        if irules:
            self.rules = irules


# Utility functions
def interpret_line(line, hdr):
    """ Generate the rule(s) specified by the input line from the rules file
    """
    # Initialize output values
    rules = []
    section_name = None
    delete_kws = []
    # Ignore comment lines in rules file
    if line[0] == '#' or len(line.strip()) == 0 or line[0] == '!':
        return rules, section_name, delete_kws
    # clean up input lines
    line = line.strip('\n')

    # strip off any comment from the line before parsing the line
    if '#' in line:
        line = line[:line.rfind('#')].strip()

    # Parse the line
    use_kws = True
    if delete_command in line:
        use_kws = False
        line = line.replace(delete_command, '').lstrip()

    if '/' in line:
        section_name = line.split('/')[1].strip()
        kws = find_keywords_in_section(hdr, section_name)
        if kws is not None:
            for kw in kws:
                if use_kws:
                    rules.append((kw, kw))
                else:
                    delete_kws.append(kw)
    else:
        kwnames = line.split()
        if '*' in line:
            kws = list(hdr[kwnames[0]].keys())
            kws2 = kws
        else:
            kws = [kwnames[0]]
            if len(kwnames) > 1:
                kws2 = [kwnames[1]]
            else:
                kws2 = kws

        lrule = None
        if len(kwnames) > 2:
            indx = line.index(kwnames[2])
            lrule = line[indx:].strip()  # rule used for new value in header

        # Interpret short-hand rules using dict
        if lrule is not None and len(lrule) > 0:
            if lrule in blender_funcs:
                lrule = blender_funcs[lrule]
            else:
                lrule = None
            # build separate rule for each kw
            for kw, kw2 in zip(kws, kws2):
                if use_kws:
                    new_rule = (kw, kw2, lrule, "ignore")
                    if new_rule not in rules:
                        rules.append(new_rule)
                else:
                    delete_kws.append(kw)
        else:
            for kw, kw2 in zip(kws, kws2):
                if use_kws:
                    new_rule = (kw, kw2)
                    if new_rule not in rules:
                        rules.append(new_rule)
                else:
                    delete_kws.append(kw)

    return rules, section_name, delete_kws


def find_keywords_in_section(hdr, title):
    """ Return a list of keyword names from hdr identified in the section
        with the specified section title.
    """
    # Indentify card indices of start and end of specified section
    sect_start = None
    sect_end = None
    for i, kw in enumerate(hdr.cards):
        if sect_start is None:
            if title in str(hdr[i]):
                sect_start = i
        else:
            if '/' in str(hdr[i]) and hdr[i] not in ['N/A', ' ', '']:
                sect_end = i
                break
    if sect_end is None:
        sect_end = len(hdr)
    if sect_start is None:
        return None

    # Now, extract the keyword names from this section
    # section_keys = hdr.ascard[sect_start+1:sect_end-1].keys()
    section_keys = list(hdr[sect_start+1:sect_end-1].keys())
    # remove any blank keywords
    while section_keys.count('') > 0:
        section_keys.remove('')

    return section_keys
