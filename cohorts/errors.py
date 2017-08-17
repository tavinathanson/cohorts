
class BamFileNotFound(Exception):

    def __init__(self, filepath=None, patient_id=None, filetype=None):
        self.filepath = filepath
        self.patient_id = patient_id
        self.filetype = 'Bamfile' if filetype is None else filetype

    def __str__(self):
        if self.patient_id is not None:
            patient_str = 'Patient {}'.format(self.patient_id)
        else:
            patient_str = 'Patient'

        if self.filepath is not None:
            print_str = 'The {} ({}) for {} was not found.'.format(self.filetype, repr(self.filepath), patient_str)
        else:
            print_str = '{} has no {}.'.format(patient_str, self.filetype)
        return print_str

class RNABamFileNotFound(BamFileNotFound):
    def __init__(self, filetype='tumor RNA Bamfile', *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class TumorBamFileNotFound(BamFileNotFound):
    def __init__(self, filetype='tumor sample bamfile', *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class NormalBamFileNotFound(BamFileNotFound):
    def __init__(self, filetype='normal sample bamfile', *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

