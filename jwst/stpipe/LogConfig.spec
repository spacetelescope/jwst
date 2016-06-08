#
# Logger configuration file
#

# Each section name is a `fnmatch` expression that matches against a
# logger name.  For example, to match all loggers involving the
# MyStep, use ``[*.MyStep]``.

[__many__]
# The minimum level at which to emit log records.  May be one of the
# standard names (DEBUG, INFO, WARNING, CRITICAL, ERROR), or a numeric
# value.  The default of 0 will emit all log records.
level = level(default=0)

# The minimum level at which to treat log records as exceptions.  May
# be one of the standard names (DEBUG, INFO, WARNING, CRITICAL,
# ERROR), or a numeric value.  The default of 0 will treat no log
# messages as exceptions.
break_level = level(default=0)

# A list of handlers to add to the logger.  Each item may be one of:
#
#   - "file:filename" Save log messages to "filename"
#   - "append:filename" Append log messages to "filename"
#   - "stdout" Send log messages to stdout
#   - "stderr" Send log messages to stderr
#
handler = string(default='')

# The format of the log messages.  See
# http://docs.python.org/library/logging.html#logrecord-attributes for
# information on creating a log format string.
format = string(default=None)
