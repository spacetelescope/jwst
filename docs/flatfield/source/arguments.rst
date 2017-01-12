Step Arguments
==============

The flat_field step has one step-specific argument, and it is only
relevant for NIRSpec data.

*  ``--flat_suffix``

``flat_suffix`` is a string, the suffix to use when constructing the
name of on optional output file for on-the-fly flat fields.  If
``flat_suffix`` is specified (and if the input data are NIRSpec),
the extracted and interpolated flat fields will be saved to a file with
this suffix.  The default (if ``flat_suffix`` was not specified) is to
not write this optional output file.
