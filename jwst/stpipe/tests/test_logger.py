import io
import contextlib
import tempfile
import os

import pytest

from .. import log as stpipe_log


@contextlib.contextmanager
def get_tempfile():
    with tempfile.NamedTemporaryFile(delete=False) as fd:
        filename = fd.name

    yield filename

    os.remove(filename)


def test_configuration():
    with get_tempfile() as logfilename:

        configuration = """
[.]
handler = file:{0}
break_level = ERROR
level = WARNING
format = '%(message)s'
""".format(logfilename)

        fd = io.BytesIO()
        fd.write(configuration.encode('latin-1'))
        fd.seek(0)
        stpipe_log.load_configuration(fd)

        log = stpipe_log.getLogger(stpipe_log.STPIPE_ROOT_LOGGER)

        log.info("Hidden")
        log.warn("Shown")

        with pytest.raises(stpipe_log.LoggedException):
            log.critical("Breaking")

        with open(logfilename, 'r') as fd:
            lines = [x.strip() for x in fd.readlines()]

        assert lines == ['Shown', 'Breaking']
