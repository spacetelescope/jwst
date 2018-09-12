import io

import pytest

from .. import log as stpipe_log


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

    log = stpipe_log.getLogger(stpipe_log.STPIPE_ROOT_LOGGER)

    log.info("Hidden")
    log.warn("Shown")

    with pytest.raises(stpipe_log.LoggedException):
        log.critical("Breaking")

    with open(logfilename, 'r') as fd:
        lines = [x.strip() for x in fd.readlines()]

    assert lines == ['Shown', 'Breaking']
