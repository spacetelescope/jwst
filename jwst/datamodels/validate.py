"""
Functions that support validation of model changes
"""

import warnings
import jsonschema
from asdf import AsdfFile
from asdf import schema as asdf_schema
from asdf.util import HashableDict
from asdf import yamlutil


class ValidationWarning(Warning):
    pass


def value_change(path, value, schema, pass_invalid_values,
                 strict_validation):
    """
    Validate a change in value against a schema.
    Trap error and return a flag.
    """
    try:
        _check_value(value, schema)
        update = True

    except jsonschema.ValidationError as error:
        update = False
        errmsg = _error_message(path, error)
        if pass_invalid_values:
            update = True
        if strict_validation:
            raise jsonschema.ValidationError(errmsg)
        else:
            warnings.warn(errmsg, ValidationWarning)
    return update


def _check_type(validator, types, instance, schema):
    """
    Callback to check data type. Skips over null values.
    """
    if instance is None:
        errors = []
    else:
        errors = asdf_schema.validate_type(validator, types,
                                           instance, schema)
    return errors


validator_callbacks = HashableDict(asdf_schema.YAML_VALIDATORS)
validator_callbacks.update({'type': _check_type})


def _check_value(value, schema):
    """
    Perform the actual validation.
    """
    if value is None:
        if schema.get('fits_required'):
            name = schema.get("fits_keyword") or schema.get("fits_hdu")
            raise jsonschema.ValidationError("%s is a required value"
                                              % name)
    else:
        validator_context = AsdfFile()
        validator_resolver = validator_context.resolver

        temp_schema = {
            '$schema':
            'http://stsci.edu/schemas/asdf-schema/0.1.0/asdf-schema'}
        temp_schema.update(schema)
        validator = asdf_schema.get_validator(temp_schema,
                                              validator_context,
                                              validator_callbacks,
                                              validator_resolver)

        value = yamlutil.custom_tree_to_tagged_tree(value, validator_context)
        validator.validate(value, _schema=temp_schema)
        validator_context.close()


def _error_message(path, error):
    """
    Add the path to the attribute as context for a validation error
    """
    if isinstance(path, list):
        spath = [str(p) for p in path]
        name = '.'.join(spath)
    else:
        name = str(path)

    error = str(error)
    if len(error) > 2000:
        error = error[0:1996] + " ..."
    errfmt = "While validating {} the following error occurred:\n{}"
    errmsg = errfmt.format(name, error)
    return errmsg
