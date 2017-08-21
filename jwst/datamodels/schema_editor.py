from __future__ import absolute_import, division, print_function, unicode_literals

# The schema editor is desgned to be run as a command line script. It can be
# run either interactively or not. To run it non-inteactively, you must set
# the options when creating  a new editor object, as in the following
# example. The options are listed in the __doc__ comment for the Schema_editor
# class.
# 
# from jwst.datamodels import schema_editor
# editor = schema_editor.Schema_editor(add=True, edit=True)
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

import os
import sys
import os.path
import inspect
import datetime
import argparse
from collections import OrderedDict

from astropy.extern import six
from six.moves.urllib import parse as urlparse
from six.moves import input

from asdf import schema as aschema
from asdf import resolver as aresolver
from asdf import generic_io
from asdf import reference
from asdf import treeutil

from . import model_base
from . import schema as mschema


def dated_directory(directory, prefix):
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


def find_directory(directory, prefix):
    """
    Find a directory with a given prefix in the specified directory
    """
    if len(directory) == 0:
        directory = os.getcwd()

    directory_path = os.path.split(directory)
    if directory_path[1].startswith(prefix):
        chosen_dir = directory
    
    else:
        chosen_dir = None
        for subdirectory in os.listdir(directory):
            if os.path.isdir(subdirectory) and subdirectory.startswith(prefix):
                if chosen_dir is None:
                    chosen_dir = subdirectory
                elif subdirectory > chosen_dir:
                    # Naming convention means that more recent directories
                    # have names later in alphabetical order
                    chosen_dir = subdirectory

    if chosen_dir is not None:
        chosen_dir = os.path.join(os.path.abspath(chosen_dir), "")

    return chosen_dir


def merge_schemas(schema):
    """
    Merge the contents of two schemas, including enums in both
    """
    def merge_dictionaries(merged_subschema, dictionary):
        for name, subdictionary in six.iteritems(dictionary):
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

def most_similar(model_schema, keyword_name):
    """
    Find the most similar name in the datamodels schema to
    a name in the keyword database
    """
    from difflib import SequenceMatcher

    def similar(model_name):
        return (SequenceMatcher(None, keyword_name, model_name).ratio() /
                (len(keyword_name) + len(model_name)))

    model_names = list(model_schema.keys())
    if len(model_names) == 0:
        similar_name = None
        
    if len(model_names) == 1:
        similar_name = model_names[0]

    else:
        model_names.sort(key=similar, reverse=True)
        ratio = similar(model_names[0]) / similar(model_names[1])
        if ratio < 1.02:
            similar_name = None
        else:
            similar_name = model_names[0]

    return similar_name


def new_path(path, model_name):
    """
    Add  the model name to the path if it is not None
    """
    new_path = path[:]
    if model_name is not None:
        new_path.append(model_name)
    return new_path
    

class Keyword_db(object):
    def __init__(self, directory=""):
        """
        Combine all the top dbs in a directory into a single db
    
        Parameters
        ----------
        directory: The directory containing the downloaded schema files.
                   If blank, use the current directory
        """
        
        directory = find_directory(directory, "JWSTDP")
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
            
            for (name, other_subschema) in six.iteritems(other_schema):
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
                
            for keyword_name, keyword_subschema in \
                 six.iteritems(keyword_schema):
                if isinstance(keyword_subschema, dict):
                    if "fits_keyword" in keyword_subschema:
                        # Save the path to any dictionary
                        # with a fits_keyword field
                        fits_name = keyword_subschema["fits_keyword"]
                        keyword_dict[fits_name] = \
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
                    merged_subschema = merge_schemas(keyword_subschema)
                    recurse(merged_subschema, keyword_dict, path)
                
        path = []
        keyword_dict = {}
        keyword_schema = self.schema
        recurse(keyword_schema, keyword_dict, path)
        
        return keyword_dict
    
    
    def resolve_references(self, url):
        """
        Resolve urls in the schema
        """
        resolver = aresolver.default_url_mapping

        def resolve_refs(node, json_id):
            if json_id is None:
                json_id = url
            if isinstance(node, dict) and '$ref' in node:
                suburl = generic_io.resolve_uri(json_id, node['$ref'])
                parts = urlparse.urlparse(suburl)
                fragment = parts.fragment
                if len(fragment):
                    suburl_path = suburl[:-(len(fragment) + 1)]
                else:
                    suburl_path = suburl
                suburl_path = resolver(suburl_path)
                if suburl_path == url:
                    subschema = schema
                else:
                    try:
                        subschema = aschema.load_schema(suburl_path,
                                                        resolver,
                                                        True)
                    except IOError as err:
                        print("Could not read " + suburl_path)
                        subschema = OrderedDict()
                subschema_fragment = reference.resolve_fragment(
                    subschema, fragment)
                return subschema_fragment
            return node

        self.schema = treeutil.walk_and_modify(self.schema, resolve_refs)

        
class Model_db(object):    
    def __init__(self):
        """
        Load the list of datamodels schema files from the schema directory
        """
        source_file = os.path.abspath(inspect.getfile(model_base.DataModel))
        self.base_url = os.path.join(os.path.dirname(source_file),
                                     'schemas', '')
        self.schema_files = []

        for filename in os.listdir(self.base_url):
            if filename.endswith(".yaml"):
                self.schema_files.append(filename)


    def __iter__(self):
        """
        Return an iterator to the list of schema files
        """
        return iter(self.schema_files)

        
    def order(self, schema):
        """
        Preserve the order of the of the schema entries in a new field
        """
        index = 0
        for (name, subschema) in six.iteritems(schema):
            if name == "properties" or isinstance(subschema, dict):
                index += 100
                subschema["ordering"] = index
                self.order(subschema)

    def read(self, schema_file):
        """
        Read a schema file into memory
        """
        fname = os.path.join(self.base_url, schema_file)
        schema = aschema.load_schema(fname, resolve_references=False)
        self.order(schema)

        return schema

        
    def save(self, schema_file, schema):
        """
        Save the modified schema back to disk
        """
        def save_dictionary(fd, schema, leading=''):
            delayed = []
            for (name, value) in six.iteritems(schema):
                if name == "properties":
                    value = self.sort(value)

                if name == "$schema":
                    delayed.append(name)
                    
                elif name == "properties" or isinstance(value, dict):
                    fd.write(leading + name + ":\n")
                    save_dictionary(fd, value, leading=leading+'  ')
                    
                elif isinstance(value, list):
                    if len(value) > 0 and isinstance(value[0], dict):
                        save_complex_list(fd, name, value, leading)
                    else:
                        fd.write(save_simple_list(leading, name, value))
                        
                elif name != "ordering":
                    fd.write(save_scalar(leading, name, value))
                    
                leading = leading.replace('-', ' ')
                
            for name in delayed:
                value = schema[name]
                if isinstance(value, list):
                    fd.write(save_list(leading, name, value))        
                elif name != "ordering":
                    fd.write(save_scalar(leading, name, value))
                
        def save_complex_list(fd, name, values, leading):
            fd.write(leading + name + ":\n")
            for value in values:
                save_dictionary(fd, value, leading=leading+'- ')
            
        def is_long_line(leading, name, value, max_length=80):
            long_line = leading + name + ': ' + str(value)
            return len(long_line) > max_length
    
        def save_long_line(leading, name, value, sep1, sep2, max_length=80):
            long_line = leading + name + ': '
            values = value.split(sep1)
                
            start = True
            line_start = 0
            prefix = leading + '  '

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
                        long_line +=  sep2 + prefix + value
                    
            long_line += '\n'

            return long_line

        def save_short_line(leading, name, value):
            short_line = leading + name + ': ' + str(value) + '\n'
            return short_line
            
        def save_simple_list(leading, name, values, max_length=80):
            sep1 = ', '
            sep2 = sep1 + '\n'
            vstr = [str(v) for v in values]
            value = '[' + sep1.join(vstr) + ']'
            if is_long_line(leading, name, value):
                line = save_long_line(leading, name, value, sep1, sep2)
            else:
                line = save_short_line(leading, name, value)
            return line
    
        def save_scalar(leading, name, value):
            sep = '|'
            if is_long_line(leading, name, value):
                value = '"\n' + value + '"\\'
                line = save_long_line(leading, name, value, '|', '|\\\n')
            else:
                line = save_short_line(leading, name, value)
            return line
    
        with open(schema_file, mode='w') as fd:
            save_dictionary(fd, schema)

    def sort(self, schema):
        def by_ordering(value):
            if type(value[1]) == OrderedDict:
                ordering = value[1].get("ordering", 0)
            else:
                ordering = 0
            return ordering
        sorted_items = sorted(schema.items(), key=by_ordering)
        return OrderedDict(sorted_items)

        
class Options(object):
    """
    Get options for running the schema editor from file, command line arguments,
    and interactive user input.
    """

    preamble = \
        """
         Welcome to the schema editor. The editor changes datamodel schema
         files to match information in the keyword database. You have control
         over which changes are made when the editor is run.
         
         Before starting, download the keyword database files from
         
         https://iwjwdmsdauiwebv.stsci.edu/portal/Mashup/Clients/jwkeywords/
         
         and unpack them. Then run this script. First it will ask you for
         the name of the directory containing the keyword database and the
         output directory you want the changed model schemas written to. 
         Then it will ask you what kind of changes you wish to make to the
         schema files. (Additions to the schema, deletions, and so on.)
         Then it will determine the differences, display them one at a time,
         and ask if you want to make the change. If you say yes, the change
         will be made to the schema. Finally, it will create a new subdirectory
         in the output directory whose name starts with "schemas" and write
         all modified schemas to that directory. It will also write the options
         to a file in  your home directory so that the next time the script
         is run it will not ask you about the same changes twice. It is a 
         text file, you can edit it or delete it. You are seeing the message 
         because no options file was found in your home directory. If one 
         is found, this message  will not be displayed. 
         
         The first two inputs are the names of the input and output directories.
         If you leave them blank, this script will use the current directory.
         All the following input to the script should either be y (yes) or n (no).
         Each question has a default answer which will be displayed as a
         capital letter (Y or N). If you hit return, the script will use
         the default value.
        """

    run_script = "do you want to continue running the editor"

    prompts = (
               ("input", "directory name containing keyword database"),
               ("output", "directory name model schemas will be written to"),
               ("add", "fields in the keyword db but not in the model to the model") ,
               ("delete", "fields in the model but not in the keyword db from the model") ,
               ("edit", "fields found in both to match the values in the keyword db") ,
               ("rename", "fields in the model to match the names in the keyword db") ,
               ("list", "changes without making them"),
               ("query", "for approval of changes before they are made"),
               ("omit", "")
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
            if isinstance(value, six.string_types):
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
        for name, prompt in self.prompts:
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
                help_text = " ".join(prompt)
                full = "--" + name
                abbrev = full[1:3]
                if isinstance(parameters[name], bool):
                    parser.add_argument(abbrev, full, help=help_text,
                                        action="store_true")
                else:
                    parser.add_argument(abbrev, full, help=help_text)
            
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
                choice = self.query_user(" ".join(prompt), default_choice)
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
        if default_choice == True:
            choices = " (Y|n)? "
        elif default_choice == False:
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
        for name, value in six.iteritems(parameters):
            try:
                current_value = getattr(editor, name)
            except AttributeError:
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
        if isinstance(parameters[name], six.string_types):
            if not isinstance(value, six.string_types):
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
            

class Schema_editor(object):
    def __init__(self, **keywds):
        """
        Initialize the editor options

        Parameters
        ----------
        All arguments are keywords

        input  : directory name containing keyword database (string)

        output : directory name model schemas will be written to  (string)

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
        self.omit = set()

        # Set attributes from keywds
        for name, value in six.iteritems(keywds):
            if hasattr(self, name):
                setattr(self, name, value)
            else:
                invalid_msg = "{0} is not a valid keyword"
                raise ValueError(invalid_msg.format(name))

        # Initialize fields with info about current file
        self.current_file_name = None
        self.current_file_changed = False

        
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
        
        # Parse the keyword database files            
        keyword_db = Keyword_db(self.input)
        keyword_dict = keyword_db.create_dict()      
        keyword_schema = keyword_db.schema

        # Loop over the model schema files, updating them from the keyword db
        first = True
        model_db = Model_db()
        for schema_file in model_db:
            self.current_file_name = schema_file
            self.current_file_changed = False
            
            model_schema = model_db.read(schema_file)
            self.edit_schema(schema_file, model_schema,
                             keyword_schema, keyword_dict)
            
            if self.current_file_changed:
                # current_file_changed is set in report_and_query
                if first:
                    first = False
                    output_dir = dated_directory(self.output, "schemas")
                schema_file = os.path.join(output_dir, schema_file)
                model_db.save(schema_file, model_schema)
        
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
        if (isinstance(keyword_value, six.string_types) and
            isinstance(model_value, six.string_types)):
            result = (keyword_value.lower().strip() ==
                      model_value.lower().strip())
    
        elif isinstance(keyword_value, float) or isinstance(model_value, float):
            try:
                result = float(keyword_value) == float(model_value)
            except:
                result = False
    
        elif isinstance(keyword_value, int) or isinstance(model_value, int):
            try:
                result = int(keyword_value) == int(model_value)
            except:
                result = False
    
        elif isinstance(keyword_value, list) and isinstance(model_value, list):
            result = set(keyword_value) == set(model_value)
    
        else:
            result = keyword_value == model_value
    
        return result

       
    def create_model_dict(self, model_schema):
        """
        Create a mapping between fits keyword name and datamodel schema name
        """
        model_dict = {}
        for (submodel_name, submodel) in six.iteritems(model_schema):
            keyword_name = self.get_keyword_value(submodel, "fits_keyword")
            
            if keyword_name is not None:
                model_dict[keyword_name] = submodel_name
    
        return model_dict


    def edit_schema(self, schema_file, model_schema,
                    keyword_schema, keyword_dict):
        """
        Top level call to edit model schema to match keyword_db
        """
        path = []
        if schema_file == "core.schema.yaml":
            model_dict = {}
            self.update_model_schema(model_schema, keyword_schema, 
                                     model_dict, path)
        else:
            self.update_fits_keywords(model_schema, keyword_schema, 
                                      keyword_dict, path)            
        
        
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
        
        
    def report_and_query(self, keyword_value, model_value, path, field):
        """
        Report a difference and optionally wait for use input
        """
        if ".".join(path) in self.omit:
            choice = False
   
        else:
            if self.query or self.list:
                self.report_schema_difference(keyword_value,
                                              model_value,
                                              path,
                                              field)      
            if self.query and not self.list:
                choice = self.options.query_user("Change", default_choice=True)
                self.update_omit(choice, path)
            else:
                choice = not self.list

        if choice:
            self.current_file_changed = True

        return choice

    
    def report_schema_difference(self, keyword_value, model_value, path, field):
        """
        Report the differences between the values of a field in the 
        datamodels schema and the keywords database
        """
        print("=== {0} {1} ===".format(self.current_file_name, '.'.join(path)))
        
        if isinstance(keyword_value, list) and isinstance(model_value, list):
            (keyword_value, model_value) = self.sort_for_report(keyword_value,
                                                                model_value)
    
        self.report_schema_value("model schema", model_value, field)
        self.report_schema_value("keyword db", keyword_value, field)


    def report_schema_value(self, schema_name, schema_value, field):
        """
        Print a message for a single field value
        """
        if field is None:
            if schema_value is None:
                print("%s not found" % schema_name)
            else:
                print("%s found" % schema_name)
    
        else:
            print("%s %s = %s" % (schema_name, field, str(schema_value)))
            if field == "enum":
                print()
    
    def schema_add_value(self, model_schema, keyword_schema,
                         keyword_name, ordering, path):
        """
        Add a new keyword to the model schema
        """
        # Logic for performing updates
        
        done = False
        if self.add:
            if self.report_and_query(keyword_schema, None, path, None):
                done = True
                model_schema[keyword_name] = OrderedDict()
                for keyword_field in ("title", "type", "enum", "default_value",
                                      "fits_keyword"):
                    keyword_value = keyword_schema.get(keyword_field)
                    if keyword_value is not None and len(keyword_value) > 0:
                        if keyword_field == "default_value":
                            model_field = "default"
                        else:
                            model_field = keyword_field
                        model_schema[keyword_name][model_field] = keyword_value

                if "properties" in keyword_schema:
                    model_schema[keyword_name]["properties"] = OrderedDict()

                model_schema[keyword_name]["ordering"] = ordering
        return done

      
    def schema_del_value(self, model_schema, model_name, path):
        """
        Delete a keyword from the model schema
        """
        done = False
        if self.delete:
            if self.report_and_query(None, model_schema, path, None):
                done = True
                del model_schema[model_name]
        return done
      
    def schema_edit_value(self, model_schema, keyword_value, model_value,
                          path, field):
        """
        Edit a field within the model schema to match the value in the
        keyword database.
        """
    
        # Logic for performing updates
        done = False
        perform = (self.edit and
                   keyword_value is not None and
                   model_value is not None)
    
        if perform:
            if self.report_and_query(keyword_value, model_value, path, field):
                done = True
                if field == "pattern":
                    keyword_value = (r'^((' +
                                     '|'.join(keyword_value) +
                                     r')\\s*\\|\\s*)+$')
                model_schema[field] = keyword_value

        return done
    
    def schema_rename_value(self, model_schema, keyword_name, model_name, path):
        """
        Rename a field in the model schema to match the keyword database
        """
        done = False
        if self.rename:
            if self.report_and_query(keyword_name, model_name, path, "name"):
                done = True
                model_schema[keyword_name] = model_schema[model_name]
                del model_schema[model_name]
        return done


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
    
    def update_fits_keywords(self, model_schema, keyword_schema, 
                             keyword_dict, model_path):
        """
        Update model schema to match keywords in keyword db
        """
        properties = self.get_keyword_value(model_schema, "properties")
        if properties is not None:
            model_schema = properties

        for model_name, model_subschema in six.iteritems(model_schema):
            if isinstance(model_subschema, dict):
                fits_name = self.get_keyword_value(model_subschema,
                                                   "fits_keyword")
                if fits_name is not None:
                    # Special case for reference file keywords
                    p_name = fits_name[0:2] == "P_"
                    import pdb ## DBG
                    pdb.set_trace() ## DBG
                    if p_name:
                        fits_name = fits_name[2:]

                    keyword_path = keyword_dict.get(fits_name)

                    if keyword_path is None:
                        if not p_name:
                            done = self.schema_del_value(model_subschema,
                                                         model_name,
                                                         new_path(model_path,
                                                                  model_name))
                    else:
                        found = True
                        keyword_subschema = keyword_schema
                        keyword_path = keyword_path.split(".")
                        for keyword_name in keyword_path:
                            if "properties" in keyword_subschema:
                                keyword_subschema = \
                                    keyword_subschema["properties"]
                                
                            if keyword_name in keyword_subschema:
                                keyword_subschema = \
                                    keyword_subschema[keyword_name]

                            else:
                                found = False
                                break
                            
                        if found:
                            done = self.update_schema_fields(keyword_subschema,
                                                             model_subschema,
                                                             new_path(model_path,
                                                                      model_name))
                            
                            if model_name != keyword_path[-1] and not p_name:
                                done = self.schema_rename_value(model_schema,
                                                                keyword_path[-1],
                                                                model_name,
                                                                new_path(model_path,
                                                                         model_name))

                else:
                    fits_hdu = self.get_keyword_value(model_subschema,
                                                   "fits_hdu")
                    if fits_hdu is not None:
                        # No info on hdu in keyword database
                        pass
                
                    else:
                        self.update_fits_keywords(model_subschema,
                                                  keyword_schema,
                                                  keyword_dict,                                          
                                                  new_path(model_path,
                                                           model_name))
                
            elif isinstance(model_subschema, list) and model_name == "allOf":
                for model_subsubschema in model_subschema:
                    if isinstance(model_subsubschema, dict):
                        self.update_fits_keywords(model_subsubschema,
                                                  keyword_schema,
                                                  keyword_dict,
                                                  model_path)
                    
                
    def update_model_schema(self, model_schema, keyword_schema, 
                            model_dict, path):
        """
        Update the model file to agree with the keyword db
        """
        properties = self.get_keyword_value(keyword_schema, "properties")
        if properties is not None:
            keyword_schema = properties

        properties = self.get_keyword_value(model_schema, "properties")
        if properties is not None:
            model_schema = properties
            model_dict = self.create_model_dict(model_schema)
    
        ordering = 0
        for keyword_name, keyword_subschema in six.iteritems(keyword_schema):
            # Find the item in the model schema (model_name and model_subschema) 
            # that corresponds to an item in the keyword database

            fits_name = self.get_keyword_value(keyword_subschema, "fits_keyword")
            if fits_name is not None:
                model_name = model_dict.get(fits_name)
        
                if model_name is None:
                    model_subschema = None
                    # Value in keyword db but not in model schema
                    ordering += 1
                    done = self.schema_add_value(model_schema,
                                                 keyword_subschema,
                                                 keyword_name,
                                                 ordering,
                                                 new_path(path,
                                                          keyword_name))
                    
                else:
                    del model_dict[fits_name]
                    model_subschema = model_schema[model_name]
                    if keyword_name != model_name:
                        # Name in keyword db differs from name in model schema
                        done = self.schema_rename_value(model_schema,
                                                        keyword_name,
                                                        model_name,
                                                        new_path(path,
                                                                 model_name))
            else:
                model_subschema = self.get_keyword_value(model_schema,
                                                         keyword_name)
                if model_subschema is not None:
                    model_name = keyword_name
                    
                else:
                    if keyword_name in ("standard", "misc", "basic", "allOf"):
                        model_name = None
                        
                    else:
                        # Name in keyword db differs from name in model schema
                        model_name = most_similar(model_schema, keyword_name)
                        if model_name is None:
                            ordering += 1
                            done = self.schema_add_value(model_schema,
                                                         keyword_subschema,
                                                         keyword_name,
                                                         ordering,
                                                         new_path(path,
                                                                  keyword_name))
                            if done:
                                model_name = keyword_name
                                
                        else:
                            done = self.schema_rename_value(model_schema,
                                                            keyword_name,
                                                            model_name,
                                                            new_path(path,
                                                                     model_name))
                                
                    if model_name:
                        if model_name in model_schema:
                            model_subschema = model_schema[model_name]
                        else:
                            model_subschema = model_schema[keyword_name]
                    else:
                        model_subschema = model_schema

            if model_subschema:
                # If a subschema of the model schema was found that matches
                # the suschema of the keyword database, descend into both
                # to check for further updates

                ordering = model_subschema.get("ordering", ordering)

                if isinstance(keyword_subschema, dict):
                    if "properties" in keyword_subschema:
                        # Skip "standard" and "misc",
                        # they are not in the model schema
                        
                        if keyword_name not in ("standard", "misc"):
                            self.update_model_schema(model_subschema,
                                                     keyword_subschema,
                                                     model_dict,
                                                     new_path(path, model_name))
                            
                    else:
                        done = self.update_schema_fields(keyword_subschema,
                                                        model_subschema,
                                                        new_path(path,
                                                                 model_name))
    
                elif (isinstance(keyword_subschema,list) and
                      keyword_name == "allOf"):
                    # Need to merge keyword db schemas if they are under
                    # an "allOf"
                    merged_subschema = merge_schemas(keyword_subschema)
                    
                    self.update_model_schema(model_subschema,
                                             merged_subschema,
                                             model_dict,
                                             path)
    
        for model_name in model_dict.values():
            # Name in model schema but not in keyword database
            done = self.schema_del_value(model_schema,
                                         model_name,
                                         new_path(path, model_name))
    

    def update_omit(self, choice, path):
        """
        Update the omit parameter based on user input
        """
        path = ".".join(path)
        if choice == True:
            self.omit.discard(path)
        elif choice == False:
            self.omit.add(path)


    def update_schema_fields(self, keyword_schema, model_schema, path):
        """
        Compare the fields of a sinlgle item bteween the datamodels schema
        and the keyword database and update when they differ
        """
        done = False
        for model_field in ('type', 'enum', 'pattern', 'default'):

            # Bridge differences in terminology between
            # keyword db and model schema

            if model_field == "pattern":
                keyword_field = "enum"
            elif model_field == "default":
                keyword_field = "default_value"
            else:
                model_field = keyword_field
                
            keyword_value = keyword_schema.get(keyword_field)
            model_value = model_schema.get(model_field)
  
            if model_field == "pattern":
                # Pattern is inside innermost set of parentehses
                patstart = model_value.rfind('(') + 1
                patend = model_value.find(')')
                model_value = model_value.split[patstart:patend]('|')

            if keyword_value == "float" and model_value == "number":
                keyword_value = "number"
                
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
                    if self.schema_edit_value(model_schema, keyword_value, 
                                              model_value, path, model_field):
                        done = True
        return done
