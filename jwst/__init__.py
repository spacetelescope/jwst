import re

from importlib.metadata import version

__version__ = version(__name__)

_regex_git_hash = re.compile(r".*\+g(\w+)")
__version__commit__ = ""
if "+" in __version__:
    commit = _regex_git_hash.match(__version__).groups()
    if commit:
        __version_commit__ = commit[0]
