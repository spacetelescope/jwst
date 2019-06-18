.. _blender_rules:

Model Blender Rules
====================

Blending models relies on rules to define how to evaluate all the input values
for a model attribute in order to determine the final output value. These rules
then get specified in the model schema for each attribute.

The rules get interpreted and applied as list or array operations that work on
the set of input values for each attribute.  The full set of pre-defined rules
includes

.. code-block:: python

    import numpy as np
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
                      'stddev': np.std,
                      'mintime': mintime,
                      'maxtime': maxtime,
                      'mindate': mindate,
                      'maxdate': maxdate,
		      'mindatetime': mindatetime,
		      'maxdatetime': maxdatetime}

The rules that should be referenced in the model schema definition are the
keys defined for `jwst.model_blender.blender_rules.blender_funcs` listed
above.  This definition illustrates how several rules are simply interfaces for
numpy array operations, while others are defined internally to `model_blender`.

.. automodapi:: jwst.model_blender.blendrules
   :skip: OrderedDict
   :skip: Time
   :skip: time
