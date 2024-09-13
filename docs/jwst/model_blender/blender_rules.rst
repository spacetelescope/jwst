.. _blender_rules:

Model Blender Rules
====================

Blending models relies on rules to define how to evaluate all the input values
for a model attribute in order to determine the final output value. These rules
are derived from the model schema for each attribute.

The rules are applied to a collection of all values to be blended (for example
if blending models that have different exposure times, the blended exposure
time will be the sum of all exposure times).

Supported rule names and corresponding functions are listed in
`jwst.model_blender.rules.RULE_FUNCTIONS`.

.. automodapi:: jwst.model_blender.rules
  :no-inheritance-diagram:
  :include-all-objects:
