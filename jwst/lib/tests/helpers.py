import os
import glob
import io


S3_BUCKET_NAME = "test-s3-data"


class MockS3Client:
    def __init__(self, s3_test_data_path):
        self.s3_test_data_path = s3_test_data_path

    def get_object(self, bucket_name, key):
        assert self.object_exists(bucket_name, key)

        with open(self._get_path(key), "rb") as f:
            return io.BytesIO(f.read())

    def object_exists(self, bucket_name, key):
        if bucket_name != S3_BUCKET_NAME:
            return False

        return os.path.isfile(self._get_path(key))

    def prefix_exists(self, bucket_name, key_prefix):
        return any(self.iterate_keys(bucket_name, key_prefix))

    def iterate_keys(self, bucket_name, key_prefix):
        if bucket_name != S3_BUCKET_NAME:
            return

        for k in self._list_keys():
            if k.startswith(key_prefix):
                yield k

    def _get_path(self, key):
        return os.path.join(self.s3_test_data_path, key)

    def _list_keys(self):
        paths = glob.glob(self.s3_test_data_path + "/**", recursive=True)
        paths = [p for p in paths if os.path.isfile(p)]
        return [p.replace(self.s3_test_data_path, "") for p in paths]
