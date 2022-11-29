import re
import sys
from pkg_resources import get_distribution, DistributionNotFound

__version_commit__ = ''
_regex_git_hash = re.compile(r'.*\+g(\w+)')

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = 'dev'

if '+' in __version__:
    commit = _regex_git_hash.match(__version__).groups()
    if commit:
        __version_commit__ = commit[0]

if sys.version_info < (3, 8):
    raise ImportError("JWST requires Python 3.8 and above.")
