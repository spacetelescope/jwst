from __future__ import absolute_import

try:
    import version_noop
except ImportError:
    __version__ = 'unknown'
    __svn_revision__ = 'unknown'
    __svn_full_info__ = 'unknown'
    __setup_datetime__ = 'unknown'
