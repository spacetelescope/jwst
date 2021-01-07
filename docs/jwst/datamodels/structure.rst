The Structure of DataModels
===========================

Datamodels allows for the creation of a new model though the usual
method of calling the __init__ method. Each type of model has its own
class and schema. The schema is specified through the class variable
schema_url. The schema gives the binding between the FITS header
keyword and/or extension and the datamodels attribute name. The chief
distinction between the two is that FITS has a flat model, datamodels
supports a hierarchical data model. The typical structure of a
datamodels class is that it first calls the __init__ method of the
base class and then initializes the required arrays of the models with
lines that look like they shouldn't do anything, for example

    self.dq = self.dq

The reason why a line like the above initializes the array is that the
access to the array on the right side of the assignment will
initialize the array to a default value if it is not already defined
and one is found in the schema. The sequence of calls is that the dot
notation invokes __getattr__, which calls _make_default if the
attribute is not defined, which in turn calls _make_default_array if
the schema says the attribute is an array. All these methods can be
found in properties.py.

The base class for Datamodels isDataModel,  in model_base.py. It takes
several arguments, the most important of which is init, which as the
name suggests, specifies how to initialize the primary data array of
the model. Init is most usually the name of a file, but can be an
already opened fits or asdf file, a numpy array, a shape tuple, or
None. If init is a shape tuple the primary data array is initialized
to its default value.

Optional arguments to __init__  can give a schema which overrides the
class schema, extensions to the schema, two flages pass_invalid_values
and strict_validation, which control the data validation, and numpy arrays
which are used to initialized the model arrays by using parameters of the
same name.

As an alternative to creating a model by initializing an object of the
specific class, you can call the open function, which is in
util.py. This function takes the same arguments as the __init_
method. If it is called with the name of a FITS file, it looks in the
primary header for a keyword named DATAMODL that conains the name of
the class to use to open the model. If that keyword is not found,
checks the dimensionality of the image and uses a generic model type
to open the image.

The base class for Datamodels loads the schema from the a file in the
schemas subdirectory. If the base class is passed a descriptor of an
already open model, it returns a shallow copy of the already open
image. This is done to speed the code, as re-opening already open
models is a common operation in the pipeline. If it is passed the
name of a file, it peeks at the first several bytes of the file to
determine the file type. This test is in filetype.py.

If the file type is a FITS file, it calls from_fits in fits_support.py
to open the file. From_fits first reads the serialized version of the
asdf tree stored in the asdf extension of the FITS file. It then walks
through the schema, which has a tree structure and uses the fields
fits_keyword and fits_hdu to locate and read items in the fits file
into the asdf tree. So items in the rest of the FITS file override
items in the asdf extension. It keeps track of the names of these
items and then makes a pass over the FITS file and writes any keywords
and hdus to another area of the asdf tree called
extra_fits. Extra_fits has subtrees for each hdu. Keywords in each hdu
not found in the schema are placed in the header subtree and data is
placed in the data subtree.  Finally, it reads the history keywords
and places them in a history structure.

To write a model back to a file, call the save method on the file. It
first calls validate_required to check the schema to see if all the
required fields are present in the model. Then it calls the function
to_fits in fits_support.py. It first creates an empty fits file and
then calls a custom validator to write the contents of the asdf tree
into this file. The functions called are defined by the dictionary
FITS_VALIDATORS found in fits_support.py. Since the validator uses the
schema as a guide, only items in the schema are added to the FITS file
at this time. After saving items mentioned in the schema, to_fits then
saves the contents of extra_fits, and then the history. Finally, it
serializes the asdf tree and writes it to the asdf extension.

Items within a model are accessed as attribute, that is, with dot
motation. The code which handles getting and setting attributes is
found in properties.py. Datamodels distinguishes between items at the
endpoints of the asdf tree and subtrees within the asdf tree. The
former are returned as scalars or numpy arrays, depending on whether
the endpoint represents a FITS keyword or data array. Subtrees are
returned as nodes. A node is an object containing the subtree as well
as the subschema which describes the subtree.  If one tries to get an
attribute that does not exist in the asdf tree, one of several things
may happen. If the attribute is not mentioned in the schema, the value
of the attribute is set to None. If it is in the schema and the schema
has a default value, the code creates the item with the default value
and then returns it. The functions that do this are _make_default and
_make_default_array, which it calls. If not only the item, but the
subtree containg the item is missing, the code throws an
AttributeError. When an attribute representing an array is accessed,
the type of the array is compared to the type in the schema and if
they are different, the array is cast to the type in the schema. The
same is true for numpy records, which represent rows in a FITS
table. The casting is done by the function gentle_asarray in util.py.

When setting or deleting an attribute, the code validates the
change. The code which does the validation can be found in
validate.py. The validator checks the values of pass_invalid_values,
which allows values inconsistent with the schema to be set, and
strict_validation, which throws an exception if the value does not
match the schema.
