"""Configure Association Logging."""

import sys
import logging
from logging.config import dictConfig
from collections import defaultdict


__all__ = ["log_config"]

DMS_DEFAULT_FORMAT = "%(asctime)s %(levelname)s pid=%(process)d src=%(name)s.%(funcName)s"


class ContextFilter:
    """Set Association Generator logging context."""

    def __init__(self):
        self.context = {}

    def __call__(self):
        return self

    def filter(self, record):
        record._context = self.context  # noqa: SLF001
        return True

    def set(self, key, value):
        self.context[key] = value


class LogLevelFilter:
    """Filter on a specific level."""

    def __init__(self, level):
        self.__level = level

    def filter(self, logrecord):
        return logrecord.levelno == self.__level


class DMSFormatter(logging.Formatter):
    """DMS-specific formatting."""

    def format(self, record):
        log_parts = []
        log_parts.append(super(DMSFormatter, self).format(record))
        for key in record._context:  # noqa: SLF001
            log_parts.append(f"{key}={record._context[key]}")  # noqa: SLF001
        log_parts.append(f'msg="{record.getMessage()}"')
        log_line = " ".join(log_parts)
        return log_line


# Define the common logging configuration
# Basic logger has the following config
base_logger = {
    "handlers": [
        "info",
        "debug",
        "default",
    ],
    "propagate": False,
}

# Basic, user-centric logging
context = ContextFilter()
base_config = {
    "version": 1,
    "disable_existing_loggers": False,
    "datefmt": "%Y%m%d%H%M",
    "formatters": {
        "info": {"format": "%(message)s", "datefmt": "cfg://datefmt"},
        "debug": {
            "format": "%(asctime)s:%(levelname)s:%(name)s.%(funcName)s:%(message)s",
            "datefmt": "cfg://datefmt",
        },
    },
    "handlers": {
        "info": {
            "class": "logging.StreamHandler",
            "formatter": "info",
            "level": logging.INFO,
            "stream": sys.stdout,
            "filters": ["info", "context"],
        },
        "debug": {
            "class": "logging.StreamHandler",
            "formatter": "debug",
            "level": logging.DEBUG,
            "stream": sys.stderr,
            "filters": ["debug", "context"],
        },
        "default": {
            "class": "logging.StreamHandler",
            "formatter": "debug",
            "level": logging.WARNING,
            "stream": sys.stderr,
            "filters": ["context"],
        },
    },
    "filters": {
        "context": {"()": context},
        "info": {"()": LogLevelFilter, "level": logging.INFO},
        "debug": {"()": LogLevelFilter, "level": logging.DEBUG},
    },
    "DMS": {
        "datefmt": "%Y%m%d%H%M",
        "logformat": ("%(asctime)s %(levelname)s pid=%(process)d src=%(name)s.%(funcName)s"),
    },
}

# DMS-specific configuration
DMS_config = {
    "formatters": {
        "info": {
            "()": DMSFormatter,
            "format": "cfg://DMS.logformat",
            "datefmt": "cfg://DMS.datefmt",
        },
        "debug": {
            "()": DMSFormatter,
            "format": "cfg://DMS.logformat",
            "datefmt": "cfg://DMS.datefmt",
        },
    }
}


def log_config(name=None, user_name=None, logger_config=None, config=None, merge=True):
    """
    Set up logging with defaults.

    logging.dictConfig is used with optional default
    configuration dict.

    Parameters
    ----------
    name : str
        The name of the logger to instantiate.
        If None, the root logger will be modified.

    user_name : str
        User-understandable name. If not specified, it will be the same
        as `name`.

    logger_config : dict
        The dict to use to setup the logger. This is used
        as the value to the key of `name`.

    config : dict
        User-specified logging configuration as specified by

    merge : bool
        Merge the user-specified config in with `base_config`.
        If `False`, just use the user-specified config.

    Returns
    -------
    Logger
        The `logging.Logger` instance.

    Notes
    -----
    Internally, the logging configuration is placed into a dict defined by
    docs.python.org/3.5/library/logging.config.html#configuration-dictionary-schema
    `base_config` defines a set of handlers, formatters, and filters that
    separate out the different levels as follows:

        `default`
            Handles `logging.WARNING` and above.

        `debug`
            Handles `logging.DEBUG` exclusively.

        `info`
            Handles `logging.INFO` exclusively.

    Given the defaults, a logger is created with the given `name`,
    with the levels as described above.

    To make modifications, the rest of this note
    describes how the parameters interact
    with the configuration dict.

    The pair of parameters `name` and `logger_config` together make up the
    dict that is the value to the top-level key `loggers`

    The parameter `config` defines the rest of the top-level keys for
    configuration dict.

    When the parameter `merge` is `True`, the following happens:
    For the `formatters` dict, the dict for the key `name` is created by
    first creating the default dict, and then doing an
    `dict.update(logger_config)` to create the logger configuration.

    For the rest of the config, and dict.update is done on each of
    the top-level keys from the `config`.
    """
    global context

    if user_name is None:
        user_name = "root" if name is None else name
    if logger_config is None:
        logger_config = {}
    if config is None:
        config = {}
    use_logger = base_logger.copy()
    use_config = defaultdict(dict, base_config)
    if merge:
        use_logger.update(logger_config)
        for key in config:
            use_config[key].update(config[key])
        use_config["loggers"] = {name: use_logger}

    dictConfig(use_config)
    logger = logging.getLogger(name)
    logger.context = context
    return logger
