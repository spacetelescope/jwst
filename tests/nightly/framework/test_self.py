import pytest


pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


def test_bigdata(_bigdata):
    with open('bigdata.txt', 'w+') as fp:
        fp.write(_bigdata + '\n')


def test_jail(_bigdata):
    with open('output.txt', 'w+') as fp:
        fp.write('hello world!\n')
