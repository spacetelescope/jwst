import io


def setup():
    from jwst.stpipe import log

    # Turn off default logging when running tests
    buffer = io.BytesIO(b"[*]\n")

    log.load_configuration(buffer)
