import io
import logging

import pytest

from .. import log as stpipe_log


@pytest.fixture(autouse=True)
def clean_up_logging():
    yield
    logging.shutdown()
    stpipe_log.load_configuration(io.BytesIO(stpipe_log.DEFAULT_CONFIGURATION))


def test_configuration(tmpdir):
    logfilename = tmpdir.join('output.log')

    configuration = """
[.]
handler = file:{0}
break_level = ERROR
level = WARNING
format = '%(message)s'
""".format(logfilename)

    fd = io.StringIO()
    fd.write(configuration)
    fd.seek(0)
    stpipe_log.load_configuration(fd)
    fd.close()

    log = stpipe_log.getLogger(stpipe_log.STPIPE_ROOT_LOGGER)

    log.info("Hidden")
    log.warning("Shown")

    with pytest.raises(stpipe_log.LoggedException):
        log.critical("Breaking")

    logging.shutdown()

    with open(logfilename, 'r') as fd:
        lines = [x.strip() for x in fd.readlines()]

    assert lines == ['Shown', 'Breaking']
