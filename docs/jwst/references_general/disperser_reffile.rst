:orphan:

.. _disperser_reffile:

DISPERSER Reference File (NIRSpec only)
---------------------------------------

:REFTYPE: DISPERSER
:Data model: `~jwst.datamodels.DisperserModel`

Reference Selection Keywords for DISPERSER
++++++++++++++++++++++++++++++++++++++++++
CRDS selects appropriate DISPERSER references based on the following keywords.
DISPERSER is not applicable for instruments not in the table.
All keywords used for file selection are *required*.

========== ======================================
Instrument Keywords
========== ======================================
NIRSpec    INSTRUME, EXP_TYPE, DATE-OBS, TIME-OBS
========== ======================================

.. include:: ../includes/standard_keywords.inc

Reference File Format
+++++++++++++++++++++
The disperser reference file contains reference data about the NIRSpec dispersers (gratings or the prism).

Files applicable to gratings have a field:

:groovedensity: Number of grooves per meter in a grating

The following fields are common for all gratings and the prism:

:grating: Name of grating
:gwa_tiltx:
    :temperatures: Temperatures measured where the GWA sensor is
    :zeroreadings: Value of GWA sensor reading which corresponds to disperser model parameters
    :tilt_model: Model of the relation between THETA_Y vs GWA_X sensor reading
:gwa_tilty:
    :temperatures: Temperatures measured where the GWA sensor is
    :zeroreadings: Value of GWA sensor reading which corresponds to disperser model parameters
    :tilt_model: Model of the relation between THETA_X vs GWA_Y sensor reading
:tilt_x: Angle (in degrees) between the grating surface and the reference surface (the mirror)
:tilt_y: Angle (in degrees) between the grating surface and the reference surface (the mirror)
:theta_x: Element alignment angle in x-axis (in degrees)
:theta_y: Element alignment angle in y-axis (in degrees)
:theta_z: Element alignment angle in z-axis (in degrees)

The prism reference file has in addition the following fields:

:angle: Angle between the front and back surface of the prosm (in degrees)
:kcoef: K coefficients of Selmeir equation, describing the material
:lcoef: L coeffficients describing the material
:tcoef: Six constants, describing the thermal behavior of the glass
:tref: Temperature (in K), used to compute the change in temperature relative to the reference temperature of the glass
:pref: Reference pressure (in ATM)
:wbound: Min and Max wavelength (in meters) for which the model is valid

