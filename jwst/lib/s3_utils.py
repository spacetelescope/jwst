"""
Experimental support for reading reference files from S3.  Use of these functions
requires installing the [aws] extras (but this module can be safely imported without
them).
"""
import atexit


__all__ = ["object_exists", "get_object", "get_client", "is_s3_uri", "split_uri"]


_CLIENT = None


def object_exists(uri):
    """
    Determine if an object exists on S3.

    Parameters
    ----------
    uri : str
        S3 URI (s3://bucket-name/some/key)

    Returns
    -------
    bool
        `True` if object exists, `False` if not.
    """
    bucket_name, key = split_uri(uri)
    return get_client().object_exists(bucket_name, key)


def get_object(uri):
    """
    Fetch the content of an object from S3.

    Parameters
    ----------
    uri : str
        S3 URI (s3://bucket-name/some/key)

    Returns
    -------
    io.BytesIO
        The content of the object.
    """
    bucket_name, key = split_uri(uri)
    return get_client().get_object(bucket_name, key)


def get_client():
    """
    Get the shared instance of ConcurrentS3Client.

    Returns
    -------
    stsci_aws_utils.s3.ConcurrentS3Client
    """
    global _CLIENT
    if _CLIENT is None:
        from stsci_aws_utils.s3 import ConcurrentS3Client
        _CLIENT = ConcurrentS3Client()
        atexit.register(_CLIENT.close)
    return _CLIENT


def is_s3_uri(value):
    """
    Determine if a value represents an S3 URI.

    Parameters
    ----------
    value : str
        Value to test.

    Returns
    -------
    bool
        `True` if value is an S3 URI, `False` if not.
    """
    return value.startswith("s3://")


def split_uri(uri):
    """
    Split an S3 URI into bucket name and key components.

    Parameters
    ----------
    uri : str
        S3 URI (s3://bucket-name/some/key)

    Returns
    -------
    str
        Bucket name URI component
    str
        Key URI component
    """
    if not uri.startswith("s3://"):
        raise ValueError("Expected S3 URI")

    bucket_name, key = uri.replace("s3://", "").split("/", 1)
    return bucket_name, key
