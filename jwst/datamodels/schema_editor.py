# The schema editor is desgned to be run as a command line script. It can be
# run either interactively or not. To run it non-inteactively, you must set
# the options when creating  a new editor object, as in the following
# example. The options are listed in the __doc__ comment for the Schema_editor
# class.
#
# from jwst.datamodels import schema_editor
# editor = schema_editor.Schema_editor(add=True, delete=True,
#                                      edit=True, rename=True)
# editor.change()
#
# To list the changes that would be made instead of making the changes,
# add the the and log options. The changes will be written to the file
# passed as the log option.
#
# from jwst.datamodels import schema_editor
# editor = schema_editor.Schema_editor(add=True, delete=True,
#                                      edit=True, rename=True,
#                                      list=True, log="changes.txt")
# editor.change()
#
# To run it interactively, first create an options object and pass it as an
# argument to the change method:
#
# from jwst.datamodels import schema_editor
# options = schema_editor.Options()
# editor = schema_editor.Schema_editor()
# editor.change(options)
#
# You can also save the options selected in a previous run to a file.
# Previous values are provided as defaults and schema changes that were
# previously again are not queried for again. The options file is written
# to your home directory unless the filename is absolute. The filename is
# passed to the options object when initialized:
#
# from jwst.datamodels import schema_editor
# options = schema_editor.Options(filename=".schema_editor_options")
# editor = schema_editor.Schema_editor()
# editor.change(options)
#
# It is possible to set all the options to the schema editor on the command
# line. Options that are not set are assumed to be false. To list the
# options, invoke the script with --help on the command line.

import argparse
from collections import OrderedDict
import datetime
import inspect
import os
import os.path
import re
import sys
from urllib.parse import urlparse
import yaml

from asdf import schema as aschema
from asdf import generic_io
from asdf import reference
from asdf import treeutil
from asdf.extension import get_default_resolver

from . import JwstDataModel


def enquote(value):
    """
    Put quote marks around value if it contains special characters
    """
    if re.search(r'[^a-zA-Z0-9_ \t\v\(\)]', value):
        value = re.sub("'", "''", value)
        value = "'" + value + "'"
    return value

def is_long_line(prefix, value, sep, max_length=80):
    """
    Test if a line is so long it needs to be folded onto the next,
    """
    value = str(value)
    if value.find(sep) >= 0:
        long_line = prefix + ": " + value
        return len(long_line) > max_length
    else:
        return False


def leading_length(prefix):
    return len(prefix) - len(prefix.lstrip())

def new_path(path, model_name):
    """
    Add  the model name to the path if it is not None
    """
    new_path = path[:]
    if model_name is not None:
        new_path.append(model_name)
    return new_path

def save_complex_list(fd, name, values, leading):
    fd.write(leading + name + ":\n")
    for value in values:
        save_dictionary(fd, value, leading=leading+'- ')

def save_dictionary(fd, schema, leading=''):
    delayed = []
    for (name, value) in schema.items():
        if name == "$schema":
            delayed.append(name)

        elif name == "properties" or isinstance(value, dict):
            fd.write(leading + name + ":\n")
            save_dictionary(fd, value, leading=leading+'  ')

        elif isinstance(value, list):
            if len(value) > 0 and isinstance(value[0], dict):
                save_complex_list(fd, name, value, leading)
            else:
                prefix = leading + name + ": "
                fd.write(save_simple_list(prefix, value))

        else:
            prefix = leading + name + ": "
            if name in ("title", "description"):
                value = enquote(value)
            fd.write(save_scalar(prefix, value))

        leading = leading.replace('-', ' ')

    for name in delayed:
        value = schema[name]
        if isinstance(value, list):
            fd.write(save_simple_list(leading, name, value))
        else:
            prefix = leading + name + ": "
            fd.write(save_scalar(prefix, value))

def save_long_line(prefix, value, sep1, sep2, max_length=80):
    """
    Converta long line into a folded string
    """
    long_line = prefix
    values = value.split(sep1)

    start = True
    line_start = 0
    leading = ' ' * (leading_length(prefix) + 2)

    for value in values:
        nl = value.find('\n')
        if nl >= 0:
            long_line += value[0:nl+1]
            value = value[nl+1:]
            line_start = len(long_line)

        if start:
            long_line += value
            start = False
        else:
            line_length = (len(long_line) + len(value) + len(sep1)
                            - line_start)

            if line_length < max_length:
                long_line += sep1 + value
            else:
                line_start = len(long_line) + len(sep2)
                long_line +=  sep2 + leading + value

    long_line += '\n'
    return long_line

def save_scalar(prefix, value):
    """
    Convert a scalar to a possibly folded string
    """
    sep1 = '|'
    sep2 = sep1 + '\\\n'
    if is_long_line(prefix, value, sep1):
        value = '"' + str(value) + '"'
        line = save_long_line(prefix, value, sep1, sep2)
    else:
        line = save_short_line(prefix, value)
    return line


def save_short_line(prefix, value):
    short_line = prefix + str(value) + '\n'
    return short_line


def save_simple_list(prefix, values, max_length=80):
    """
    Convert a list, possibly folding at the separator
    """
    sep1 = ', '
    sep2 = sep1 + '\n'
    vstr = [str(v) for v in values]
    value = '[' + sep1.join(vstr) + ']'
    if is_long_line(prefix, value, sep1):
        line = save_long_line(prefix, value, sep1, sep2)
    else:
        line = save_short_line(prefix, value)
    return line

class Keyword_db:
    def __init__(self, directory=""):
        """
        Combine all the top dbs in a directory into a single db

        Parameters
        ----------
        directory: The directory containing the downloaded schema files.
                   If blank, use the current directory
        """

        directory = self.find_directory(directory)
        if directory is None:
            raise ValueError("Cannot locate keyword database directory")

        self.schema = None
        for filename in os.listdir(directory):
            if filename.startswith("top."):
                keyword_db = os.path.abspath(os.path.join(directory, filename))
                schema = aschema.load_schema(keyword_db, resolve_references=False)
                if self.schema is None:
                    self.schema = schema
                else:
                    try:
                        self.combine_schemas(schema)
                    except ValueError as err:
                        raise ValueError(filename + ": " + str(err))

        self.resolve_references(directory)

    def combine_schemas(self, other_schema):
        """
        Combine another schema into the keyword databse schema
        """
        def combine_dictionaries(this_schema, other_schema):
            error_msg = "Unrecognized field in schema: "

            for (name, other_subschema) in other_schema.items():
                if "properties" not in other_subschema:
                    raise ValueError(error_msg + name)

                if name in this_schema:
                    this_subschema = this_schema[name]

                    if "properties" not in this_subschema:
                        raise ValueError(error_msg + name)

                    if ("allOf" in this_subschema["properties"] or
                        "allOf" in other_subschema["properties"]):
                        try:
                            combine_lists(this_subschema, other_subschema)
                        except ValueError:
                            raise ValueError(error_msg + name)

                    elif ("$ref" in this_subschema["properties"] and
                          "$ref" in other_subschema["properties"]):
                        if (this_subschema["properties"]["$ref"] !=
                            other_subschema["properties"]["$ref"]):
                            try:
                                combine_lists(this_subschema, other_subschema)
                            except ValueError:
                                raise ValueError(error_msg + name)

                    else:
                        raise ValueError(error_msg + name)

                else:
                    this_schema[name] = other_subschema

        def combine_lists(this_schema, other_schema):
            combination = {}
            for schema in (this_schema, other_schema):
                if "allOf" in schema["properties"]:
                    for item in schema["properties"]["allOf"]:
                        if "$ref" in item:
                            ref = item["$ref"]
                            combination[ref] = item
                        else:
                            raise ValueError

                elif "$ref" in schema["properties"]:
                    ref = schema["properties"]["$ref"]
                    combination[ref] = schema["properties"]

                else:
                    raise ValueError
            this_schema["properties"] = {"allOf" : list(combination.values())}

        this_schema = self.schema["properties"]["meta"]["properties"]
        other_schema = other_schema["properties"]["meta"]["properties"]
        combine_dictionaries(this_schema, other_schema)


    def create_dict(self):
        """
        Create a dictionary mapping fits keyword names to keyword db paths
        """
        def recurse(keyword_schema, keyword_dict, path):
            if "properties" in keyword_schema:
                keyword_schema = keyword_schema["properties"]

            for keyword_name, keyword_subschema in keyword_schema.items():
                if isinstance(keyword_subschema, dict):
                    if "fits_hdu" in keyword_subschema:
                        if keyword_subschema["fits_hdu"] in ('PRIMARY', 'SCI'):
                            if "fits_keyword" in keyword_subschema:
                                # Save the path to any dictionary
                                # with a fits_keyword field
                                keyword_fits_name = keyword_subschema["fits_keyword"]
                                keyword_fits_name = keyword_fits_name.strip()
                                keyword_dict[keyword_fits_name] = \
                                    '.'.join(new_path(path, keyword_name))

                    else:
                        # Skip "standard" and "misc",
                        # they are not in the model schema
                        if keyword_name not in ("standard", "misc"):
                            recurse(keyword_subschema,
                                    keyword_dict,
                                    new_path(path, keyword_name))

                elif (isinstance(keyword_subschema, list) and
                      keyword_name == "allOf"):
                    # Need to combine keyword db schemas if they are under
                    # an "allOf"
                    merged_subschema = self.merge_schemas(keyword_subschema)
                    recurse(merged_subschema, keyword_dict, path)

        path = []
        keyword_dict = {}
        keyword_schema = self.schema
        recurse(keyword_schema, keyword_dict, path)

        return keyword_dict


    def find_directory(self, directory):
        """
        Return full path to requested directory

        Notes
        -----
        There was originally some functionality involving a prefix.
        This functionality has been removed.
        """
        if len(directory) == 0:
            directory = os.getcwd()

        chosen_dir = os.path.join(os.path.abspath(directory), "")

        return chosen_dir

    def merge_schemas(self, schema):
        """
        Merge the contents of two schemas, including enums in both
        """
        def merge_dictionaries(merged_subschema, dictionary):
            for name, subdictionary in dictionary.items():
                if name in merged_subschema:
                    if ("enum" in merged_subschema[name] and
                        "enum" in subdictionary):
                            merge_enums(merged_subschema[name],
                                          subdictionary)
                    else:
                        merged_subschema[name] = subdictionary
                else:
                    merged_subschema[name] = subdictionary


        def merge_enums(merged_subschema, dictionary):
            merged_subschema["enum"] = list(set(merged_subschema["enum"]) |
                                            set(dictionary["enum"]))

        merged_subschema = OrderedDict()
        for dictionary in schema:
            if merged_subschema:
                merge_dictionaries(merged_subschema, dictionary)
            else:
                merged_subschema = dictionary

        return merged_subschema

    def resolve_references(self, url):
        """
        Resolve urls in the schema
        """
        resolver = get_default_resolver()

        def resolve_refs(node, json_id):
            if json_id is None:
                json_id = url
            if isinstance(node, dict) and '$ref' in node:
                suburl = generic_io.resolve_uri(json_id, node['$ref'])
                parts = urlparse(suburl)
                fragment = parts.fragment
                if len(fragment):
                    suburl_path = suburl[:-(len(fragment) + 1)]
                else:
                    suburl_path = suburl
                suburl_path = resolver(suburl_path)

                try:
                    subschema = aschema.load_schema(suburl_path,
                                                    resolver,
                                                    True)
                except IOError:
                    print("Could not read " + suburl_path)
                    subschema = OrderedDict()

                subschema_fragment = reference.resolve_fragment(
                    subschema, fragment)

                return subschema_fragment
            return node

        self.schema = treeutil.walk_and_modify(self.schema, resolve_refs)


class Model_db:
    """
    Load the list of datamodels schema files from the schema directory

    Parameters
    ----------
    exclude: None or [str[,...]]
        List of `DataModel` schemas to ignore.

    Attributes
    ----------
    base_url: str
        The location of the `DataModel` schemas

    schema_files: [str[...]]
        The schema files loaded.
    """

    def __init__(self, exclude=None):
        if exclude is None:
            exclude = []

        source_file = os.path.abspath(inspect.getfile(JwstDataModel))
        self.base_url = os.path.join(os.path.dirname(source_file),
                                     'schemas', '')
        self.schema_files = []

        for filename in os.listdir(self.base_url):
            if filename.endswith(".yaml") and filename not in exclude:
                self.schema_files.append(filename)


    def __len__(self):
        return len(self.schema_files)


    def __iter__(self):
        """
        Return an iterator to the list of schema files
        """
        return iter(self.schema_files)


    def read(self, schema_file):
        """
        Read a schema file into memory
        """
        fname = os.path.join(self.base_url, schema_file)
        schema = aschema.load_schema(fname, resolve_references=False)

        return schema


    def save(self, schema_file, schema):
        """
        Save the modified schema back to disk
        """
        with open(schema_file, mode='w') as fd:
            save_dictionary(fd, schema)


class Options:
    """
    Get options for running the schema editor from file, command line arguments,
    and interactive user input.
    """

    preamble = \
        """
         Welcome to the schema editor. The editor changes datamodel schema
         files to match information in the keyword database. You have control
         over which changes are made when the editor is run.

         Before starting, refer to the process document found at:
         https://innerspace.stsci.edu/display/SCSB/Keyword+Dictionary+JSON+vs.+YAML+Comparison+Process


         Retrieve the JWST Keyword Database as described. Then run this script.
         First it will ask you for the name of the directory containing the
         keyword database and the output directory you want the changed model
         schemas written to. Then it will ask you what kind of changes you wish
         to make to the schema files. (Additions to the schema, deletions, and
         so on.) Then it will determine the differences, display them one at a
         time, and ask if you want to make the change. If you say yes, the
         change will be made to the schema. Finally, it will create a new
         subdirectory in the output directory whose name starts with "schemas"
         and write all modified schemas to that directory.

         If schema are modified, and keywords are omitted, a file
         `aaa_omitted_keywords.yaml` will be created. This file can be used as
         the input to the `--omit_file` option

         The first two inputs are the names of the input and output directories.
         If you leave them blank, this script will use the current directory.
         All the following input to the script should either be y (yes) or n (no).
         Each question has a default answer which will be displayed as a
         capital letter (Y or N). If you hit return, the script will use
         the default value.

        """

    run_script = "do you want to continue running the editor"

    # Command line options. 2 or 3-tuple where:
    #    - name: Parameter name and full option name
    #    - help: The help text
    #    - Create a short option. If not present or True, a single character option,
    #      based on the first character of the name, is created.
    prompts = (
        ("input", "directory name containing keyword database"),
        ("output", "directory name model schemas will be written to"),
        ("add", "fields that exist in the keyword db but are not in the model schemas"),
        ("delete", "fields in the model but not in the keyword db from the model"),
        ("edit", "fields found in both to match the values in the keyword db"),
        ("rename", "fields in the model to match the names in the keyword db"),
        ("list", "changes without making them"),
        ("query", "for approval of changes before they are made"),
        ("omit_file", "containing the list of parameters to omit", False),
        ("exclude_file", "containing the list of DataModel schema to exclude from the comparison", False)
    )


    def __init__(self, filename=None):
        """
        Create option file reader / writer

        Parameters
        ----------
        filename: The name of the file stroing options from the previous run
                  If set to None, no file will be read or written.
        """
        self.batch = len(sys.argv) > 1

        if filename is None:
            self.filename = filename
            self.first_time = True

        else:
            if os.path.isabs(filename):
                self.filename = filename
            else:
                home = os.path.expanduser("~")
                self.filename = os.path.join(home, filename)
            self.first_time = not os.path.isfile(self.filename)

    def coerce_type(self, value, current_value):
        """
        Coerce value to type of the current value of a parameter
        """
        if isinstance(current_value, bool):
            if isinstance(value, str):
                value = value[0]
                value = value.lower()
                if value == "y":
                    value = True
                elif value == "n":
                    value = False
                else:
                    value = None

            elif not isinstance(value, bool):
                value = None

        elif isinstance(current_value, set):
            if hasattr(value, "__iter__"):
                value = set(value)
            else:
                value = set((value, ))

        return value

    def get(self, editor):
        """
        Initialize fields from defaults, command line, and options file

        Parameters
        ----------
        editor: The schema editor whose options are being set
        """
        parameters = self.get_parameters(editor)
        if self.filename is not None and not self.first_time:
            parameters = self.read(parameters)

        if self.batch:
            parameters = self.parse_command_line(parameters)
        else:
            if self.first_time:
                if not self.query_run(self.preamble, self.run_script):
                    return False

            parameters = self.query_options(parameters)

        self.set_parameters(editor, parameters)
        return True

    def get_parameters(self, editor):
        """
        Save the object's fields into a dictionary
        """
        parameters = {}
        for name, prompt, *abbrev in self.prompts:
            parameters[name] = getattr(editor, name)

        return parameters


    def parse_command_line(self, parameters):
        """
        Parse the parameters from the command line
        """
        parser = argparse.ArgumentParser()

        for prompt in self.prompts:
            if prompt[1]:
                name = prompt[0]
                help_text = " ".join(prompt[:2])
                full = "--" + name
                if len(prompt) > 2 and not prompt[2]:
                    option_names = (full,)
                else:
                    abbrev = full[1:3]
                    option_names = (abbrev, full)
                if isinstance(parameters[name], bool):
                    parser.add_argument(*option_names, help=help_text,
                                        action="store_true")
                else:
                    parser.add_argument(*option_names, help=help_text)

        args = parser.parse_args()

        for prompt in self.prompts:
            if prompt[1]:
                name = prompt[0]
                try:
                    value = getattr(args, name)
                except AttributeError:
                    value = None

                if value is not None:
                    value = self.coerce_type(value, parameters[name])
                    self.type_check(parameters, name, value)
                    parameters[name] = value

        return parameters


    def parse_parameter(self, parameters, n, name, value):
        """
        Parse a single parameter value, check it, and save it
        """
        if len(name) == 0:
            if n == 0:
                return
            else:
                raise ValueError("Error at line {} in option file".format(n+1))

        value = value.replace("\t", " ")
        if value.find(" ") >= 0:
            value = value.split()

        value = self.coerce_type(value, parameters[name])
        self.type_check(parameters, name, value)
        parameters[name] = value

        return


    def query_options(self, parameters):
        """
        Query the user for global option values
        """
        for prompt in self.prompts:
            if prompt[1]:
                default_choice = parameters[prompt[0]]
                choice = self.query_user(" ".join(prompt[:2]), default_choice)
                self.type_check(parameters, prompt[0], choice)
                parameters[prompt[0]] = choice
        return parameters


    def query_run(self, preamble, prompt):
        """
        Explain how to run the script and ask the user to continue
        """
        lines = preamble.split("\n")
        for line in lines:
            print(line.lstrip())

        choice = self.query_user(prompt, default_choice=False)
        return choice


    def query_user(self, prompt, default_choice=None):
        """
        Build the prompt and query the user
        """
        if default_choice is True:
            choices = " (Y|n)? "
        elif default_choice is False:
            choices = " (y|N)? "
        elif default_choice is None:
            choices = "? "
        elif len(default_choice) > 0:
            choices = " (%s)? " % default_choice
        else:
            choices = "? "

        choice = None
        prompt = prompt + choices
        while choice is None:
            choice = input(prompt)
            if len(choice) == 0:
                choice = default_choice
            else:
                choice = self.coerce_type(choice, default_choice)

        return choice


    def read(self, parameters):
        """
        Read the contents of an option file into a hash
        """
        n = 0
        name = ""
        value = ""
        with open(self.filename, "r") as fd:
            for line in fd:
                line = line.rstrip()
                if len(line) == 0:
                    name = ""
                    value = ""
                else:
                    if line.startswith(" ") or line.startswith("\t"):
                        value += line
                    else:
                        self.parse_parameter(parameters, n, name, value)

                        (name, eq, value) = line.partition("=")
                        name = name.strip()
                        value = value.strip()
                n += 1

            self.parse_parameter(parameters, n, name, value)

        return parameters


    def set_parameters(self, editor, parameters):
        """
        Set a field on the object after validating it
        """
        for name, value in parameters.items():
            if not hasattr(editor, name):
                invalid_msg = "{0} is not a valid parameter"
                raise ValueError(invalid_msg.format(name))

            value = self.coerce_type(value, parameters[name])
            self.type_check(parameters, name, value)
            setattr(editor, name, value)

    def type_check(self, parameters, name, value):
        """
        Check that the type of a new value agrees with the current type
        """
        badtype = "Invalid paramater value for %s (%s)"
        if isinstance(parameters[name], str):
            if not isinstance(value, str):
                raise ValueError(badtype % (name, value))

        elif isinstance(parameters[name], set):
            if isinstance(value, bool):
                raise ValueError(badtype % (name, value))

        else:
            ptype = type(parameters[name])
            if not isinstance(value, ptype):
                raise ValueError(badtype %  (name, value))


    def write(self, editor):
        """
        Write the contents of a hash as an option file

        Parameters
        ----------
        editor: The schema editor whose options are being saved to file
        """
        if self.filename is not None:
            with open(self.filename, "w") as fd:
                parameters = self.get_parameters(editor)

                for name in sorted(parameters.keys()):
                    value = parameters[name]
                    if isinstance(value, bool):
                        if value:
                            value = "yes"
                        else:
                            value = "no"
                    elif isinstance(value, set):
                        value = " " + "\n\t".join(list(value))

                    fd.write("{0} = {1}\n".format(name, value))


class Schema_editor:
    def __init__(self, **keywds):
        """
        Initialize the editor options

        Parameters
        ----------
        All arguments are keywords

        input  : directory name containing keyword database (string)

        output : directory name model schemas will be written to  (string)

        log    : file containing reported modifications to files

        add    : add fields in the keyword db but not in the model to the model (bool)

        delete : delete fields in the model but not in the keyword db from the model (bool)

        edit   : edit fields found in both to match the values in the keyword db (bool)

        rename : rename fields in the model to match the names in the keyword db (bool)

        list   : list changes without making them (bool)

        query  : query for approval of changes before they are made (bool)

        omit   : omit model schema keywords from editing (set)
        """
        # Default values of editor attributes
        self.query = False
        self.list = False
        self.add = False
        self.delete = False
        self.edit = False
        self.rename = False
        self.input = ""
        self.output = ""
        self.log = ""
        self.omit = set()
        self.omit_file = ''
        self.exclude = []
        self.exclude_file = ''
        self.model_db = None

        # Set attributes from keywds
        for name, value in keywds.items():
            if hasattr(self, name):
                setattr(self, name, value)
            else:
                invalid_msg = "{0} is not a valid keyword"
                raise ValueError(invalid_msg.format(name))

        # Initialize fields with info about current file
        self.current_file_name = None
        self.current_file_changed = False

        # Initialize the built in regex pattern
        _builtin_regexes = [
            '', 'NAXIS[0-9]{0,3}', 'BITPIX', 'XTENSION', 'PCOUNT', 'GCOUNT',
            'EXTEND', 'BSCALE', 'BUNIT', 'BZERO', 'BLANK', 'DATAMAX', 'DATAMIN',
            'EXTNAME', 'EXTVER', 'EXTLEVEL', 'GROUPS', 'PYTPE[0-9]',
            'PSCAL[0-9]', 'PZERO[0-9]', 'SIMPLE', 'TFIELDS',
            'TBCOL[0-9]{1,3}', 'TFORM[0-9]{1,3}', 'TTYPE[0-9]{1,3}',
            'TUNIT[0-9]{1,3}', 'TSCAL[0-9]{1,3}', 'TZERO[0-9]{1,3}',
            'TNULL[0-9]{1,3}', 'TDISP[0-9]{1,3}', 'HISTORY'
            ]

        self.builtin_regex = re.compile(
            '|'.join('(^{0}$)'.format(x) for x in _builtin_regexes))

    def add_subschema(self, model_schema, keyword_schema,
                      keyword_dict, fits_dict):
        """
        Add subschema from keyword db to model schema
        """
        for keyword_fits_name in keyword_dict:
            if not self.builtin_fits_keyword(keyword_fits_name):
                if keyword_fits_name not in fits_dict:
                    keyword_path = keyword_dict[keyword_fits_name]
                    keyword_path = keyword_path.split('.')

                    path = []
                    add = False
                    keyword_subschema = keyword_schema
                    model_subschema = model_schema
                    for name in keyword_path:
                        if "properties" in keyword_subschema:
                            keyword_subschema = keyword_subschema["properties"]
                        if "properties" in model_subschema:
                            model_subschema = model_subschema["properties"]

                        test = self.get_keyword_value(keyword_subschema, name)
                        if test is None:
                            break
                        else:
                            keyword_subschema = test

                        test = self.get_keyword_value(model_subschema, name)
                        path.append(name)

                        if test is None:
                            add = True
                            break
                        else:
                            model_subschema = test

                    if add:
                        self.schema_add_subschema(fits_dict,
                                                  model_subschema,
                                                  keyword_subschema,
                                                  path)


    def builtin_fits_keyword(self, fits_name):
        """
        Test if a fits keyword is built in
        """
        return self.builtin_regex.match(fits_name) is not None

    def change(self, options=None):
        """
        Change datamodels schema files to match info in the keyword db

        Parameters
        ----------
        options : An options object, used for querying for attribute values
                  interactively. Do not set if using in a program.
        """

        # Initialize the editor object from file and command line
        self.options = options
        if self.options is None:
            self.query = False
        else:
            if not self.options.get(self):
                return

        # If an omit file was specified, add the contents to the omit set.
        if self.omit_file:
            with open(self.omit_file) as fh:
                self.omit.update(yaml.safe_load(fh))

        # If an exclude file was specified, add the list of `DataModel` schema
        # to the exclusion list.
        if self.exclude_file:
            with open(self.exclude_file) as fh:
                self.exclude.extend(yaml.safe_load(fh))

        # Set output file for messages depending on list,
        # so output can be captured to a file
        if self.query or self.log == "":
            self.fd = sys.stdout
        else:
            self.fd = open(self.log, "w")

        # Is anything going to be done?
        if not any((
                self.add, self.delete, self.edit, self.rename,
        )):
            raise RuntimeError(
                'No operation has been requested. Set at least one of:\n'
                '\tadd, delete, edit, or rename'
            )

        # If not listing, get the output directory setup.
        if not self.list:
            output_dir = self.dated_directory(self.output, "schemas")

        # Parse the keyword database files
        keyword_db = Keyword_db(self.input)
        keyword_dict = keyword_db.create_dict()
        keyword_schema = keyword_db.schema

        # Loop over the model schema files, updating them from the keyword db

        fits_dict = {}
        model_db = Model_db(exclude=self.exclude)
        self.model_db = model_db
        for schema_file in model_db:
            model_path = []
            model_schema = model_db.read(schema_file)
            self.match_fits_keywords(schema_file, model_schema, keyword_schema,
                                     keyword_dict, fits_dict, model_path)

        for schema_file in model_db:
            self.current_file_name = schema_file
            self.current_file_changed = False

            model_path = []
            model_schema = model_db.read(schema_file)
            self.edit_schema(model_schema, keyword_schema,
                             keyword_dict, fits_dict,  model_path)

            if self.current_file_changed:
                # current_file_changed is set in report_and_query
                schema_file = os.path.join(output_dir, schema_file)
                model_db.save(schema_file, model_schema)

        schema_file = 'core.schema.yaml'
        self.current_file_changed = False
        model_schema = model_db.read(schema_file)
        self.add_subschema(model_schema, keyword_schema,
                           keyword_dict, fits_dict)

        if self.current_file_changed:
            schema_file = os.path.join(output_dir, schema_file)
            model_db.save(schema_file, model_schema)

        # Save the list of omitted parameters
        if not self.list and self.omit:
            omit_path = os.path.join(output_dir, 'aaa_omitted_keywords.yaml')
            with open(omit_path, 'w') as fh:
                yaml.safe_dump(list(self.omit), fh)

        # Write the opject attributes back to disk
        if self.options is not None:
            self.options.write(self)

    def check_type(self, model_schema, default_value):
        """
        Check if default value agrees with data type
        """
        valid = True
        model_type = model_schema.get("type")
        if model_type is not None:
            if model_type == "number":
                try:
                    float(default_value)
                except ValueError:
                    valid = False

            elif model_type == "integer":
                try:
                    int(default_value)
                except ValueError:
                    valid = False

        return valid


    def compare_schema_values(self, keyword_value, model_value):
        """
        Compare iwo values for type specific kind of equality
        """
        if (isinstance(keyword_value, str) and
            isinstance(model_value, str)):
            result = (keyword_value.lower().strip() ==
                      model_value.lower().strip())

        elif isinstance(keyword_value, float) or isinstance(model_value, float):
            try:
                result = float(keyword_value) == float(model_value)
            except ValueError:
                result = False

        elif isinstance(keyword_value, int) or isinstance(model_value, int):
            try:
                result = int(keyword_value) == int(model_value)
            except ValueError:
                result = False

        elif isinstance(keyword_value, list) and isinstance(model_value, list):
            result = set(keyword_value) == set(model_value)

        else:
            result = keyword_value == model_value

        return result


    def dated_directory(self, directory, prefix):
        """
        Create a dated directory in a specified directory
        """
        d = datetime.date.today()

        for i in range(100):
            subdirectory = ("%s_%04d_%02d_%02d_%02d" %
                            (prefix, d.year, d.month, d.day, i))

            subdirectory = os.path.join(directory, subdirectory, "")
            if not os.path.exists(subdirectory):
                os.makedirs(subdirectory)
                return subdirectory

        return None


    def edit_schema(self, model_schema, keyword_schema,
                    keyword_dict, fits_dict, model_path):
        """
        Edit model schema to match keywords in keyword db
        """
        properties = self.get_keyword_value(model_schema, "properties")
        if properties is not None:
            model_schema = properties

        model_names = list(model_schema.keys())
        for model_name in model_names:
            model_subschema = model_schema[model_name]
            new_model_path = new_path(model_path, model_name)
            if isinstance(model_subschema, dict):
                model_fits_name = self.get_keyword_value(model_subschema,
                                                         "fits_keyword")

                if model_fits_name is not None:
                    if not self.builtin_fits_keyword(model_fits_name):
                        # Special case for reference file keywords
                        p_name = model_fits_name[0:2] == "P_"

                        keyword_fits_name = fits_dict.get(model_fits_name)
                        if keyword_fits_name is None:
                            keyword_path = None
                        else:
                            keyword_path = keyword_dict.get(keyword_fits_name)

                        if keyword_path is None:
                            if not p_name:
                                self.schema_del_value(model_schema,
                                                            model_name,
                                                            new_model_path)

                        else:
                            keyword_subschema = self.find_subschema(keyword_schema,
                                                                    keyword_path)
                            if keyword_subschema is not None:
                                keyword_path = keyword_path.split('.')
                                self.update_schema_fields(keyword_subschema,
                                                            model_subschema,
                                                            new_model_path)

                                if not p_name and model_name != keyword_path[-1]:
                                    self.schema_rename_value(model_schema,
                                                            keyword_path[-1],
                                                            model_name,
                                                            new_model_path)

                                if model_fits_name != keyword_fits_name:
                                    fits_dict[keyword_fits_name] = \
                                        fits_dict[model_fits_name]
                                    del fits_dict[model_fits_name]

                else:
                    fits_hdu = self.get_keyword_value(model_subschema, "fits_hdu")

                    if fits_hdu is None:
                        self.edit_schema(model_subschema,
                                         keyword_schema,
                                         keyword_dict,
                                         fits_dict,
                                         new_model_path)

                        if not self.has_properties(model_subschema):
                            del model_schema[model_name]

            elif isinstance(model_subschema, list) and model_name == "allOf":
                for model_subsubschema in model_subschema:
                    if isinstance(model_subsubschema, dict):
                        self.edit_schema(model_subsubschema,
                                         keyword_schema,
                                         keyword_dict,
                                         fits_dict,
                                         model_path)



    def find_subschema(self, schema, path):
        """
        Trace a path thru the schema. If not found, return partial
        match if partial is True, None if partial is false
        """
        subschema = schema
        if type(path) != list:
            path = path.split('.')

        for name in path:
            if "properties" in subschema:
                subschema = subschema["properties"]

            subschema = self.get_keyword_value(subschema, name)
            if subschema is None:
                return None

        return subschema


    def get_keyword_value(self, submodel, keyword_name):
        """
        Lppk for keyword name in submodel, return value if found
        """
        keyword_value = None
        if isinstance(submodel, dict):
            if keyword_name in submodel:
                keyword_value = submodel[keyword_name]

            elif "allOf" in submodel:
                for subsubmodel in submodel["allOf"]:
                    if (isinstance(subsubmodel, dict) and
                        keyword_name in subsubmodel):
                        keyword_value = subsubmodel[keyword_name]
                        break

        return keyword_value


    def has_properties(self, model_schema):
        """
        Determine if schema dictionary has any real properties
        """
        properties = self.get_keyword_value(model_schema, "properties")
        if properties is not None:
            model_schema = properties

        return len(model_schema) > 0


    def match_fits_keywords(self, schema_file, model_schema, keyword_schema,
                            keyword_dict, fits_dict, model_path):
        """
        Match model names to keyword db names by fits keyword name
        """
        properties = self.get_keyword_value(model_schema, "properties")
        if properties is not None:
            model_schema = properties

        for model_name in model_schema:
            model_subschema = model_schema[model_name]
            new_model_path = new_path(model_path, model_name)
            if isinstance(model_subschema, dict):
                model_fits_name = self.get_keyword_value(model_subschema,
                                                         "fits_keyword")

                if model_fits_name is not None:
                    # Special case for reference file keywords
                    model_fits_name = model_fits_name.strip()
                    p_name = model_fits_name[0:2] == "P_"
                    if p_name:
                        model_name = model_name[2:]
                        keyword_path = keyword_dict.get(model_fits_name[2:])
                    else:
                        keyword_path = keyword_dict.get(model_fits_name)

                    if keyword_path is None:
                        keyword_subschema = self.find_subschema(keyword_schema,
                                                                new_model_path)
                    else:
                        keyword_subschema = self.find_subschema(keyword_schema,
                                                                keyword_path)

                    if keyword_subschema is not None:
                        keyword_fits_name = keyword_subschema.get("fits_keyword")
                        if keyword_fits_name is not None:
                            fits_dict[model_fits_name] = keyword_fits_name.strip()

                else:
                    fits_hdu = self.get_keyword_value(model_subschema,
                                                      "fits_hdu")
                    if fits_hdu is None:
                        self.match_fits_keywords(schema_file,
                                                 model_subschema,
                                                 keyword_schema,
                                                 keyword_dict,
                                                 fits_dict,
                                                 new_model_path)

            elif isinstance(model_subschema, list) and model_name == "allOf":
                for model_subsubschema in model_subschema:
                    if isinstance(model_subsubschema, dict):
                        self.match_fits_keywords(schema_file,
                                                 model_subsubschema,
                                                 keyword_schema,
                                                 keyword_dict,
                                                 fits_dict,
                                                 new_model_path)


    def pattern_from_enum(self, value):
        """
        Convert the value of an enum field to a regular expression pattern
        """
        return r'^((' + '|'.join(value) +  r')\\s*\\|\\s*)+$'

    def report_and_query(self, verb, title, keyword_value, model_value,
                         path, field):
        """
        Report a difference and optionally wait for use input
        """
        if ".".join(path) in self.omit:
            choice = False

        else:
            if self.query or self.list:
                self.report_schema_difference(title,
                                              keyword_value,
                                              model_value,
                                              path,
                                              field)
            if self.query and not self.list:
                choice = self.options.query_user(verb, default_choice=True)
                self.update_omit(choice, path)
            else:
                choice = not self.list

        if choice:
            self.current_file_changed = True

        return choice


    def report_schema_difference(self, title, keyword_value, model_value,
                                 path, field):
        """
        Report the differences between the values of a field in the
        datamodels schema and the keywords database
        """
        self.fd.write("=== {0} ===\n".format('.'.join(path)))
        if title is not None:
            self.fd.write(title + "\n")

        if isinstance(keyword_value, list) and isinstance(model_value, list):
            (keyword_value, model_value) = self.sort_for_report(keyword_value,
                                                                model_value)

        if field == "pattern":
            model_value  = self.pattern_from_enum(model_value)
            keyword_value = self.pattern_from_enum(keyword_value)

        self.report_schema_value("model schema", model_value, field)
        self.report_schema_value("keyword db", keyword_value, field)


    def report_schema_value(self, schema_name, schema_value, field):
        """
        Write a message for a single field value
        """

        if field is None:
            if schema_value is None:
                self.fd.write("%s not found\n" % schema_name)
            elif "properties" in schema_value:
                self.fd.write("%s is a group\n" % schema_name)
            else:
                schema_value = schema_value.get('fits_keyword')
                if schema_value is None:
                    self.fd.write("%s not found\n" % schema_name)
                else:
                    schema_value = schema_value.lower()
                    self.fd.write("%s found as %s\n" % (schema_name, schema_value))

        else:
            prefix = "%s %s = " % (schema_name, field)
            if type(schema_value) == list:
                self.fd.write(save_simple_list(prefix, schema_value))
            else:
                self.fd.write(save_scalar(prefix, schema_value))
            if field in ("enum", "pattern"):
                self.fd.write("\n")

    def schema_add_subschema(self, fits_dict, model_schema, keyword_schema,
                             keyword_path):
        """
        Add a subschema tree to the model schema
        """
        self.schema_add_value(model_schema, keyword_schema, keyword_path)

        if "properties" in keyword_schema:
            model_name = keyword_path[-1]
            model_subschema = self.find_subschema(model_schema, model_name)
            keyword_schema = keyword_schema["properties"]

            for keyword_name in keyword_schema:
                keyword_subschema = self.get_keyword_value(keyword_schema,
                                                           keyword_name)

                keyword_fits_name = \
                    self.get_keyword_value(keyword_subschema,
                                           "fits_keyword")

                if keyword_fits_name is None:
                    add = True
                elif keyword_fits_name in fits_dict:
                    add = False
                elif self.builtin_fits_keyword(keyword_fits_name):
                    add = False
                else:
                    add = True

                if add:
                    self.schema_add_subschema(fits_dict,
                                              model_subschema,
                                              keyword_subschema,
                                              new_path(keyword_path,
                                                       keyword_name))


    def schema_add_value(self, model_schema, keyword_schema, keyword_path):
        """
        Add a new keyword to the model schema
        """
        if self.add:
            title = keyword_schema.get('title')
            if self.report_and_query("Add", title, keyword_schema, None,
                                  keyword_path, None):
                keyword_name = keyword_path[-1]
                model_schema[keyword_name] = OrderedDict()

                for keyword_field in ('title', 'type', 'enum', 'default_value',
                                      'fits_keyword'):
                    keyword_value = \
                        self.strip_blanks(keyword_schema.get(keyword_field))
                    if keyword_value is not None and len(keyword_value) > 0:
                        if keyword_field == "default_value":
                            model_field = "default"
                        else:
                            model_field = keyword_field
                        model_schema[keyword_name][model_field] = keyword_value

                if "properties" in keyword_schema:
                    model_schema[keyword_name]["properties"] = OrderedDict()

    def schema_del_value(self, model_schema, model_name, path):
        """
        Delete a keyword from the model schema
        """
        if self.delete:
            title = model_schema.get('title')
            if self.report_and_query("Delete", title, None,
                                    model_schema[model_name], path, None):
                del model_schema[model_name]

    def schema_edit_value(self, model_schema, keyword_value, model_value,
                          path, field):
        """
        Edit a field within the model schema to match the value in the
        keyword database.
        """

        # Logic for performing updates
        perform = (self.edit and
                   keyword_value is not None and
                   model_value is not None)

        if perform:
            title = model_schema.get('title')
            if self.report_and_query("Change", title, keyword_value, model_value,
                                     path, field):
                if field == "pattern":
                    keyword_value = self.pattern_from_enum(keyword_value)
                model_schema[field] = keyword_value

    def schema_rename_value(self, model_schema, keyword_name, model_name, path):
        """
        Rename a field in the model schema to match the keyword database
        """
        if self.rename:
            title = model_schema.get('title')
            if self.report_and_query("Rename", title, keyword_name, model_name,
                                     path, "name"):
                model_schema[keyword_name] = model_schema[model_name]
                del model_schema[model_name]

    def sort_for_report(self, keyword_value, model_value):
        """
        Reorder the values in two lists of values so that the elements that
        differ are first in the list
        """
        common_values = list(set(keyword_value) & set(model_value))

        different_values = list(set(keyword_value) - set(common_values))
        sorted_keyword_value = different_values + common_values

        different_values = list(set(model_value) - set(common_values))
        sorted_model_value = different_values + common_values

        return (sorted_keyword_value, sorted_model_value)

    def strip_blanks(self, value):
        """
        Remove leading and trailing blanks from values
        """
        def strip_single(elem):
            if isinstance(elem, str):
                elem = elem.strip()
            return elem

        if isinstance(value, list):
            for i in range(len(value)):
                value[i] = strip_single(value[i])
        else:
            value = strip_single(value)

        return value

    def update_omit(self, choice, path):
        """
        Update the omit parameter based on user input
        """
        path = ".".join(path)
        if choice is True:
            self.omit.discard(path)
        elif choice is False:
            self.omit.add(path)


    def update_schema_fields(self, keyword_schema, model_schema, path):
        """
        Compare the fields of a sinlgle item bteween the datamodels schema
        and the keyword database and update when they differ
        """
        p_name = path[-1][0:2] == "p_"
        for model_field in ('fits_keyword', 'type', 'enum', 'pattern', 'default'):
            # Skip checking fits keyword for reference file fields
            if p_name and model_field == 'fits_keyword':
                continue

            # Bridge differences in terminology between
            # keyword db and model schema

            if model_field == "pattern":
                keyword_field = "enum"
            elif model_field == "default":
                keyword_field = "default_value"
            else:
                keyword_field = model_field

            keyword_value = self.strip_blanks(keyword_schema.get(keyword_field))
            model_value = self.strip_blanks(model_schema.get(model_field))

            if model_field == "pattern" and model_value is not None:
                # Pattern is inside innermost set of parentehses
                patstart = model_value.rfind('(') + 1
                patend = model_value.find(')')
                model_value = model_value[patstart:patend].split('|')

            if keyword_value == "float" and model_value == "number":
                keyword_value = "number"

            if keyword_value == "int" and model_value == "integer":
                keyword_value = "integer"

            if keyword_value == "" and model_value is None:
                keyword_value = None

            # Check if value in keyword db differs from value in model schema

            if not self.compare_schema_values(keyword_value, model_value):
                # Check if default value agrees with model schema type
                if model_field == "default":
                    valid = self.check_type(model_schema, keyword_value)
                else:
                    valid = True

                if valid:
                    self.schema_edit_value(model_schema, keyword_value,
                                           model_value, path, model_field)
