import pytest


@pytest.mark.bigdata
def test_bigdata():
    with open('bigdata.txt', 'w+') as fp:
        fp.write("bigdata defined" + '\n')


def test_jail(_jail):
    with open('output.txt', 'w+') as fp:
        fp.write('hello world!\n')
