Step Arguments
==============

The persistence step has three step-specific arguments.

*  ``--input_trapsfilled``

``input_trapsfilled`` is the name of the most recent trapsfilled file
for the current detector.  If this is not specified, an array of zeros
will be used as an initial value.  If this is specified, it will be used
to predict persistence for the input science file.
The step writes an output trapsfilled file, and that could be used
as input to the persistence step for a subsequent exposure.

*  ``--flag_pers_cutoff``

If this floating-point value is specified, pixels that receive a
persistence correction greater than or equal to ``flag_pers_cutoff`` DN
(the default is 40) are flagged in the PIXELDQ array of the
output file with the DQ value "PERSISTENCE".

*  ``--save_persistence``

If this boolean parameter is specified and is True (the default is False),
the persistence that was subtracted (group by group, integration by
integration) will be written to an output file with suffix "_output_pers".
