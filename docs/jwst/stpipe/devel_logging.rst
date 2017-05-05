=======
Logging
=======

The logging in stpipe is built on the Python standard library’s
`logging` module.  For detailed information about logging, refer to
the documentation there.  This document basically outlines some simple
conventions to follow so that the configuration mechanism described in
:ref:`user-logging` works.

Logging from a Step or Pipeline
===============================

Each Step instance (and thus also each Pipeline instance) has a `log`
member, which is a Python `logging.Logger` instance.  All messages
from the Step should use this object to log messages.  For example,
from a `process` method::

    self.log.info("This Step wants to say something")

Logging from library code
=========================

Often, you may want to log something from code that is oblivious to
the concept of stpipe Steps.  In that case, stpipe has special code
that allows library code to use any logger and have those messages
appear as if they were coming from the step that used the library.
All the library code has to do is use a Python `logging.Logger` as
normal::

    import logging

    # ...
    log = logging.getLogger()

    # If the log on its own won’t emit, neither will it in the
    # context of an stpipe step, so make sure the level is set to
    # allow everything through
    log.setLevel(logging.DEBUG)

    def my_library_call():
        # ...
        log.info("I want to make note of something")
        # ...
