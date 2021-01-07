# Copyright (C) 2010 Association of Universities for Research in Astronomy(AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
"""
Logging setup etc.
"""
import fnmatch
import io
import logging
import os
import sys
import threading

from ..extern.configobj.configobj import ConfigObj
from ..extern.configobj import validate

from . import config_parser

STPIPE_ROOT_LOGGER = 'stpipe'
DEFAULT_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
DEFAULT_CONFIGURATION = b"""
[*]
handler = stderr
level = INFO
"""

MAX_CONFIGURATION = b"""
[*]
handler = stderr
level = DEBUG
"""


###########################################################################
# LOGS AS EXCEPTIONS

class LoggedException(Exception):
    """
    This is an exception used when a log record is converted into an
    exception.

    Use its `record` member to get the original `logging.LogRecord`.
    """
    def __init__(self, record):
        self.record = record
        Exception.__init__(self, record.getMessage())


class BreakHandler(logging.Handler):
    """
    A handler that turns logs of a certain severity or higher into
    exceptions.
    """
    _from_config = True

    def emit(self, record):
        raise LoggedException(record)


#########################################################################
# LOGGING CONFIGURATION

# A dictionary mapping patterns to
log_config = {}


class LogConfig():
    """
    Stores a single logging configuration.

    Parameters
    ----------
    name : str
        The `fnmatch` pattern used to match the logging class

    handler, level, break_level, format : str
        See LogConfig.spec for a description of these values.
    """
    def __init__(self, name, handler=None, level=logging.NOTSET,
                 break_level=logging.NOTSET, format=None):
        if name in ('', '.', 'root'):
            name = '*'
        self.name = name
        self.handler = handler
        if not isinstance(self.handler, list):
            if self.handler.strip() == '':
                self.handler = []
            else:
                self.handler = [x.strip() for x in self.handler.split(',')]
        self.level = level
        self.break_level = break_level
        if format is None:
            format = DEFAULT_FORMAT
        self.format = format

    def match(self, log_name):
        """
        Returns `True` if `log_name` matches the pattern of this
        configuration.
        """
        if log_name.startswith(STPIPE_ROOT_LOGGER):
            log_name = log_name[len(STPIPE_ROOT_LOGGER) + 1:]
            if fnmatch.fnmatchcase(log_name, self.name):
                return True
        return False

    def get_handler(self, handler_str):
        """
        Given a handler string, returns a `logging.Handler` object.
        """
        if handler_str.startswith("file:"):
            return logging.FileHandler(handler_str[5:], 'w', 'utf-8', True)
        elif handler_str.startswith("append:"):
            return logging.FileHandler(handler_str[7:], 'a', 'utf-8', True)
        elif handler_str == 'stdout':
            return logging.StreamHandler(sys.stdout)
        elif handler_str == 'stderr':
            return logging.StreamHandler(sys.stderr)
        else:
            raise ValueError("Can't parse handler {0!r}".format(handler_str))

    def apply(self, log):
        """
        Applies the configuration to the given `logging.Logger`
        object.
        """
        for handler in log.handlers[:]:
            if hasattr(handler, '_from_config'):
                log.handlers.remove(handler)

        # Set a handler
        for handler_str in self.handler:
            handler = self.get_handler(handler_str)
            handler._from_config = True
            handler.setLevel(self.level)
            log.addHandler(handler)

        # Set the log level
        log.setLevel(self.level)

        # Set the break level
        if self.break_level != logging.NOTSET:
            log.addHandler(BreakHandler(self.break_level))

        formatter = logging.Formatter(self.format)
        for handler in log.handlers:
            if (isinstance(handler, logging.Handler) and
                hasattr(handler, '_from_config')):
                handler.setFormatter(formatter)

    def match_and_apply(self, log):
        """
        If the given `logging.Logger` object matches the pattern of
        this configuration, it applies the configuration to it.
        """
        if self.match(log.name):
            self.apply(log)


def load_configuration(config_file):
    """
    Loads a logging configuration file.  The format of this file is
    defined in LogConfig.spec.

    Parameters
    ----------
    config_file : str or readable file-like object
    """
    def _level_check(value):
        try:
            value = int(value)
        except ValueError:
            pass
        try:
            value = logging._checkLevel(value)
        except ValueError:
            raise validate.VdtTypeError(value)
        return value

    spec = config_parser.load_spec_file(LogConfig)
    config = ConfigObj(config_file, raise_errors=True, interpolation=False)
    val = validate.Validator()
    val.functions['level'] = _level_check
    config_parser.validate(config, spec, validator=val)

    log_config.clear()

    for key, val in config.items():
        log_config[key] = LogConfig(key, **val)

    for log in logging.Logger.manager.loggerDict.values():
        if isinstance(log, logging.Logger):
            for cfg in log_config.values():
                cfg.match_and_apply(log)


def getLogger(name=None):
    log = logging.getLogger(name)

    return log


def _find_logging_config_file():
    files = [
        "stpipe-log.cfg",
        "~/.stpipe-log.cfg",
        "/etc/stpipe-log.cfg"
        ]

    for file in files:
        file = os.path.expanduser(file)
        if os.path.exists(file):
            return os.path.abspath(file)

    buffer = io.BytesIO(DEFAULT_CONFIGURATION)
    return buffer


###########################################################################
# LOGGING DELEGATION


class DelegationHandler(logging.Handler):
    """
    A handler that delegates messages along to the currently active
    `Step` logger.  It only delegates messages that come from outside
    of the `stpipe` heirarchy, in order to prevent infinite recursion.

    Since we could be multi-threaded and each thread may be running a
    different thread, we need to manage a dictionary mapping the
    current thread to the Step's logger on that thread.
    """
    def __init__(self, *args, **kwargs):
        self._logs = {}
        logging.Handler.__init__(self, *args, **kwargs)

    def emit(self, record):
        log = self.log
        if (log is not None and
            not record.name.startswith(STPIPE_ROOT_LOGGER)):
            record.name = log.name
            log.handle(record)

    @property
    def log(self):
        return self._logs.get(threading.current_thread(), None)

    @log.setter
    def log(self, log):
        assert (log is None or
                (isinstance(log, logging.Logger) and
                 log.name.startswith(STPIPE_ROOT_LOGGER)))
        self._logs[threading.current_thread()] = log

# Install the delegation handler on the root logger.  The Step class
# uses the `delegator` instance to change what the current Step logger
# is.
log = getLogger()
delegator = DelegationHandler()
delegator.log = getLogger(STPIPE_ROOT_LOGGER)
log.addHandler(delegator)

logging_config_file = _find_logging_config_file()
if logging_config_file is not None:
    load_configuration(logging_config_file)

logging.captureWarnings(True)
