import sys
if sys.version_info[0] >= 3:
    __builtins__['unicode'] = str
    __builtins__['basestring'] = str

from .step import Step
from .pipeline import Pipeline
from .linear_pipeline import LinearPipeline

__version__ = '0.9.3'
