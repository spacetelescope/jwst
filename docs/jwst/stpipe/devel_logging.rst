=======
Logging
=======

The logging in stpipe is built on the Python standard library’s
`logging` module.  For detailed information about logging, refer to
the documentation there.

By default, stpipe will pick up and configure any logger used in
the `jwst` library code, so developers do not need to do anything
special to log messages from their modules: just get a logger by
name, and log to it.  By convention, loggers should be named
for the module they are used in.

All the library code has to do is use a Python `logging.Logger` as
normal::

    import logging

    log = logging.getLogger(__name__)

    def my_library_call():
        # ...
        log.info("I want to make note of something")
        # ...
