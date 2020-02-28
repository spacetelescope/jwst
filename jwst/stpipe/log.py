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


class LogConfig:
    """
    Stores a single logging configuration.

    Parameters
    ----------
    name : str
        The `fnmatch` pattern used to match the logging class

    handler : str, list

    level, break_level, format : str
        See LogConfig.spec for a description of these values.
    """
    def __init__(self, name, handler=None, level=None, break_level=None, format=None):
        if name in ('', '.', 'root'):
            name = '*'

        self.name = name
        self.handler = handler

        if isinstance(self.handler, str) and self.handler is not None:
            self.handler = [x.strip() for x in self.handler.split(',')] if self.handler.strip() else []

        self.level = level if level is not None else logging.NOTSET
        self.break_level = break_level if break_level is not None else logging.NOTSET
        self.format = format if format is not None else DEFAULT_FORMAT

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

    @staticmethod
    def get_handler(handler_str):
        """
        Given a handler string, returns a `logging.Handler` object.
        """
        if handler_str.startswith("file:"):
            return logging.FileHandler(handler_str[5:], 'w', 'utf-8', True)

        if handler_str.startswith("append:"):
            return logging.FileHandler(handler_str[7:], 'a', 'utf-8', True)

        if handler_str == 'stdout':
            return logging.StreamHandler(sys.stdout)

        if handler_str == 'stderr':
            return logging.StreamHandler(sys.stderr)

        raise ValueError(f"Can't parse handler {handler_str!r}")

    def apply(self, logger):
        """
        Applies the configuration to the given `logging.Logger`
        object.
        """
        for handler in logger.handlers[:]:
            if hasattr(handler, '_from_config'):
                logger.handlers.remove(handler)

        # Set a handler
        if self.handler is not None:
            for handler_str in self.handler:
                handler = self.get_handler(handler_str)
                handler._from_config = True
                handler.setLevel(self.level)
                logger.addHandler(handler)

        # Set the log level
        logger.setLevel(self.level)

        # Set the break level
        if self.break_level != logging.NOTSET:
            logger.addHandler(BreakHandler(self.break_level))

        formatter = logging.Formatter(self.format)
        for handler in logger.handlers:
            if isinstance(handler, logging.Handler) and hasattr(handler, '_from_config'):
                handler.setFormatter(formatter)

    def match_and_apply(self, logger):
        """
        If the given `logging.Logger` object matches the pattern of
        this configuration, it applies the configuration to it.
        """
        if self.match(logger.name):
            self.apply(logger)


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
    logger = logging.getLogger(name)

    return logger


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
        logger = self.log

        if logger is not None and not record.name.startswith(STPIPE_ROOT_LOGGER):
            record.name = logger.name
            logger.handle(record)

    @property
    def log(self):
        return self._logs.get(threading.current_thread(), None)

    @log.setter
    def log(self, logger):
        assert (logger is None or
                (isinstance(logger, logging.Logger) and
                 logger.name.startswith(STPIPE_ROOT_LOGGER)))
        self._logs[threading.current_thread()] = logger


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
