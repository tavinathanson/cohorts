
class MissingData(Exception):

    def __init__(self, filepath=None, patient_id=None, filetype='Data'):
        self.filepath = filepath
        self.patient_id = patient_id
        self.filetype = filetype

    def __str__(self):
        if self.patient_id is not None:
            patient_str = "Patient {}".format(self.patient_id)
        else:
            patient_str = "Patient"

        if self.filepath is not None:
            print_str = "The {} ({}) for {} was not found.".format(self.filetype, repr(self.filepath), patient_str)
        else:
            print_str = "{} has no {}.".format(patient_str, self.filetype)
        return print_str

class MissingBamFile(MissingData):
    def __init__(self, filetype='Bamfile', *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class MissingRNABamFile(MissingBamFile):
    def __init__(self, filetype="tumor RNA Bamfile", *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class MissingTumorBamFile(MissingBamFile):
    def __init__(self, filetype="tumor sample bamfile", *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class MissingNormalBamFile(MissingBamFile):
    def __init__(self, filetype="normal sample bamfile", *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class MissingHLAType(MissingData):
    def __init__(self, filetype="HLA alleles", *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)

class MissingVariantFile(MissingData):
    def __init__(self, filetype="vcf/maf files", *args, **kwargs):
        super().__init__(filetype=filetype, *args, **kwargs)
