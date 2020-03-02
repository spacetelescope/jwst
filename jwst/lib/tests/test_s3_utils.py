import pytest

from jwst.lib import s3_utils
from . import helpers


@pytest.fixture
def s3_text_file(s3_root_dir):
    path = str(s3_root_dir.join("test.txt"))
    with open(path, "w") as text_file:
        print("foo", file=text_file)

    return path


def test_object_exists(s3_text_file):
    assert s3_utils.object_exists("s3://test-s3-data/test.txt") is True
    assert s3_utils.object_exists("s3://test-s3-data/missing.fits") is False
    assert s3_utils.object_exists("s3://missing-bucket/test.txt") is False


def test_get_object(s3_text_file):
    assert s3_utils.get_object("s3://test-s3-data/test.txt").read() == b"foo\n"


def test_get_client(s3_text_file):
    assert isinstance(s3_utils.get_client(), helpers.MockS3Client)


def test_is_s3_uri(s3_text_file):
    assert s3_utils.is_s3_uri("s3://test-s3-data/test.fits") is True
    assert s3_utils.is_s3_uri("some/filesystem/path") is False


def test_split_uri(s3_text_file):
    assert s3_utils.split_uri("s3://test-s3-data/key") == ("test-s3-data", "key")
    assert s3_utils.split_uri("s3://test-s3-data/some/longer/key") == ("test-s3-data", "some/longer/key")
