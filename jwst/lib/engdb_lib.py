"""Engineering DB common library"""
from collections import namedtuple


# Define the returned value tuple.
EngDB_Value = namedtuple('EngDB_Value', ['obstime', 'value'])
