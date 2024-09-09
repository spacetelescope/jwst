from datetime import time

import numpy as np


__all__ = ['RULE_FUNCTIONS', 'Blender', 'make_blender']


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


def _isotime(time_str):
    hms = [float(i) for i in time_str.split(':')]
    sec_ms = hms[2] - int(hms[2])
    isotime = time(int(hms[0]), int(hms[1]), int(hms[2]), int(sec_ms * 1000000))
    return isotime


RULE_FUNCTIONS = {
    'multi': multi,
    'mean': np.mean,
    'sum': np.sum,
    'max': np.max,
    'min': np.min,

     # retained date/time names for backwards compatibility
     # as these all assume ISO8601 format the lexical and
     # chronological sorting match
    'mintime': min,
    'maxtime': max,
    'mindate': min,
    'maxdate': max,
    'mindatetime': min,
    'maxdatetime': max
}


class Blender:
    def __init__(self, blend_function):
        self.blend_function = blend_function
        self.values = []

    def accumulate(self, value):
        self.values.append(value)

    def finalize(self):
        if not self.values:
            return None
        return self.blend_function(self.values)


def make_blender(rule):
    return Blender(RULE_FUNCTIONS[rule])
