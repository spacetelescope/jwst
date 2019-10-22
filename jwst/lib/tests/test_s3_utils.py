from .. import s3_utils
from . import helpers


def test_object_exists():
    assert s3_utils.object_exists("s3://test-s3-data/test.txt") is True
    assert s3_utils.object_exists("s3://test-s3-data/missing.fits") is False
    assert s3_utils.object_exists("s3://missing-bucket/test.txt") is False


def test_get_object():
    assert s3_utils.get_object("s3://test-s3-data/test.txt").read() == b"foo"


def test_get_client():
    assert isinstance(s3_utils.get_client(), helpers.MockS3Client)


def test_is_s3_uri():
    assert s3_utils.is_s3_uri("s3://test-s3-data/test.fits") is True
    assert s3_utils.is_s3_uri("some/filesystem/path") is False


def test_split_uri():
    assert s3_utils.split_uri("s3://test-s3-data/key") == ("test-s3-data", "key")
    assert s3_utils.split_uri("s3://test-s3-data/some/longer/key") == ("test-s3-data", "some/longer/key")
