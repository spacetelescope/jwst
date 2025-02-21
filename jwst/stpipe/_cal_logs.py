"""Cal log handling code."""

import getpass
import logging
import re
import socket
import time

_IPV4_REGEX = re.compile(r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}")
_IPV6_REGEX = re.compile(
    r"([0-9a-fA-F]{1,4}:){1,7}:|:(:[0-9a-fA-F]{1,4}){1,7}|([0-9a-fA-F]{1,4}:){1,7}[0-9a-fA-F]{0,4}(:[0-9a-fA-F]{1,4}){1,7}"
)
_HOSTNAME = socket.gethostname()
_USER = getpass.getuser()


def _scrub(msg):
    """
    Scrub sensitive information from a message.

    Parameters
    ----------
    msg : str
        The string to scrub

    Returns
    -------
    scrubbed : str
        The scrubbed string
    """
    if _USER in msg:
        return ""
    if _HOSTNAME in msg:
        return ""
    if re.search(_IPV4_REGEX, msg):
        return ""
    # if re.search(_IPV6_REGEX, msg):
    #    return ""
    return msg


class _ScrubbingFormatter(logging.Formatter):
    """Formatter that removes sensitive information."""

    def format(self, record):
        return _scrub(super().format(record))


_LOG_FORMATTER = _ScrubbingFormatter(
    "%(asctime)s.%(msecs)03d - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
_LOG_FORMATTER.converter = time.gmtime
