import os
import tempfile
from contextlib import contextmanager
from urllib.parse import urlparse

from google.cloud import storage


_GOOGLE_STORAGE_SCHEMA = "gs"

## Things that can go here ##
class GoogleStorageEnvironmentError(EnvironmentError): pass
class GoogleStorageValueError(ValueError): pass
class GoogleStorageIOError(IOError): pass
## End of things that can go wrong


class GoogleStorageIO:
    """
    This is a simple, light-weight wrapper aound Google Cloud's SDK
    that comes with pretty low-level REST-like operations. Given that
    we mostly depend on files sitting either on the NFS or the GS,
    I thought this could facilitate a few things and make it easier
    to switch between NFS and GCloud worlds.

    Note that Google's SDK tries its best to locate a configuration or
    authentication file around for it to access the data; but if
    you might need to use the following command once (if you haven't done so
    before) so that the auth part goes much smoother (basically hands-off):

    ```
        # Beta channel
        $ gcloud beta auth application-default login
        # Stable channel
        $ gcloud auth application-default login
    ```

    More information on the auth procedures:
        https://cloud.google.com/docs/authentication/
    """
    def __init__(self, *args, **kwargs):
        """
        You can pass authentication related configuration directly down to
        gcloud client, so if the default automated auth method doesn't work
        this will help you go with other options: CREDs, service accounts, etc.
        """
        try:
            self._client = storage.Client(*args, **kwargs)
        except Exception as e:
            err = ("Couldn't initialize Google Storage client, possibly due "
                   "to authentication issues. Please check your gcloud "
                   "auth setup: {}".format(e))
            raise GoogleStorageEnvironmentError(err)


    @property
    def client(self):
        """
        A standard Google Storage client. Exposing this carefully, because
        there might be helpful storage-based operations that are not
        implemented in this class; and this might come handy in such
        cases.
        """
        return self._client


    def list_files(self, gsuri, **kwargs):
        bucket_name, gs_rel_path = self.parse_uri(gsuri)
        bucket = self._client.get_bucket(bucket_name)
        return bucket.list_blobs(**kwargs)


    @client.setter
    def client(self, new_client):
        self._client = new_client


    def parse_uri(self, gsuri):
        gs_parts = urlparse(gsuri)
        assert gs_parts.scheme == _GOOGLE_STORAGE_SCHEMA
        bucket_name = gs_parts.netloc
        gs_rel_path = gs_parts.path
        # And apparently gcloud doesn't like relative URIs
        # that starts with a `/`, so let's get rid of that:
        gs_rel_path = gs_rel_path if gs_rel_path[0] != '/' else gs_rel_path[1::]
        return (bucket_name, gs_rel_path)


    def upload_file(self, localpath, gsuri):
        # And now request the handles for bucket and the file
        bucket_name, rel_path = self.parse_uri(gsuri)
        bucket = self.client.get_bucket(bucket_name)
        ublob = storage.Blob(rel_path, bucket)
        ublob.upload_from_filename(localpath)


    def download_to_path(self, gsuri, localpath, binary_mode=False, tmpdir=None):
        """
        This method is analogous to "gsutil cp gsuri localpath", but in a
        programatically accesible way. The only difference is that we
        have to make a guess about the encoding of the file to not upset
        downstream file operations. If you are downloading a VCF, then
        "False" is great. If this is a BAM file you are asking for, you
        should enable the "binary_mode" to make sure file doesn't get
        corrupted.

        gsuri: full GS-based URI, e.g. gs://cohorts/rocks.txt
        localpath: the path for the downloaded file, e.g. /mnt/cohorts/yep.txt
        binary_mode: (logical) if yes, the binary file operations will be
                     used; if not, standard ascii-based ones.
        """
        bucket_name, gs_rel_path = self.parse_uri(gsuri)
        # And now request the handles for bucket and the file
        bucket = self._client.get_bucket(bucket_name)
        # Just assignment, no downloading (yet)
        ablob = bucket.get_blob(gs_rel_path)
        if not ablob:
            raise GoogleStorageIOError(
                "No such file on Google Storage: '{}'".format(gs_rel_path))

        # A tmp file to serve intermediate phase
        # should be on same filesystem as localpath
        tmp_fid, tmp_file_path = tempfile.mkstemp(text=(not binary_mode),
                                                  dir=tmpdir)
        # set chunk_size to reasonable default
        # https://github.com/GoogleCloudPlatform/google-cloud-python/issues/2222
        ablob.chunk_size = 1<<30
        # Download starts in a sec....
        ablob.download_to_filename(client=self._client, filename=tmp_file_path)
        # ... end download ends. Let's move our finished file over.

        # You will see that below, instead of directly writing to a file
        # we are instead first using a different file and then move it to
        # its final location. We are doing this because we don't want
        # corrupted/incomplete data to be around as much as possible.
        return os.rename(tmp_file_path, localpath)

    def get_blob(self, gsuri):
        bucket_name, gs_rel_path = self.parse_uri(gsuri)
        bucket = self._client.get_bucket(bucket_name)
        ablob = bucket.get_blob(gs_rel_path)
        return ablob


class GoogleStorageFile(object):
    """
    This class is a context manager and a context decorator (future) that
    is meant to make working on GS files feel similar to local ones.
    Here is a one way via the `with` keyword:
        ```
        gsuri = 'gs://cohorts/biomarker.csv'
        with GoogleStorageFile(gsuri, 'r+') as fgs:
            contents = fgs.read()
            # since `fgs` is a local one, it is easy to
            # to work against it.
            # ...
            gfs.write(...)
        # And boom, all the tmp files are cleaned and our file is
        # uploaded back to the GS in its new form.
        ```
    """
    def __init__(self, gcio, gsuri, mode):
        self.gcio = gcio
        self.gsuri = gsuri
        self.is_read = 'r' in mode
        self.is_write = ('w' in mode) or ('r+' in mode)
        self.is_binary = 'b' in mode
        self.mode = mode
        self._localfile = None

    def __enter__(self):
        # Get a temp file and use it as a local file to immitate
        # direct GS access by downloading contents onto that file.
        tmp_fid, tmp_file_path = tempfile.mkstemp(text=(not self.is_binary))
        self.gcio.download_to_path(self.gsuri, tmp_file_path)
        self._localfile = open(tmp_file_path, self.mode)
        self._localfile_path = tmp_file_path
        return self._localfile

    def __exit__(self, *args):
        # We are done with the open file, so let's close it
        self._localfile.close()
        # If write mode is on, then upload the altered one to GS
        if self.is_write:
            bucket_name, rel_path = self.gcio.parse_uri(self.gsuri)
            bucket = self.gcio.client.get_bucket(bucket_name)
            ublob = storage.Blob(rel_path, bucket)
            ublob.upload_from_filename(self._localfile_path)
        # And because of the type of the temp file we created,
        # we are responsible for deleting it when we are finshed:
        os.remove(self._localfile_path)
