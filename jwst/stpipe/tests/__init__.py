import io


def setup():
    """Configure the logging system to turn off default logging when running tests."""
    from jwst.stpipe import log

    # Turn off default logging when running tests
    buffer = io.BytesIO(b"[*]\n")

    log.load_configuration(buffer)
