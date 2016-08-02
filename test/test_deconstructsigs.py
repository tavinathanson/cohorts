from __future__ import print_function

from . import data_path, generated_data_path
from .data_generate import generate_vcfs
from .test_count import make_cohort

from os import path, makedirs

FILE_FORMAT_1 = "patient_format1_%s.vcf"
FILE_FORMAT_2 = "patient_format2_%s.vcf"


def test_deconstructsigs_output():
    """
    Run on a subset (3 patients) of TCGA
    -> TCGA
    ->
    """

def test_input_generation():
    """
    A test of make_deconstructsigs_patient_inputs()

    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])
        variants = cohort.load_variants()
        inputs_path = generated_data_path("deconstructsigs-inputs")
        print "inputs_path ", inputs_path
        inputs_dir = path.dirname(inputs_path)
        print "inputs_dir", inputs_dir
        if not path.exists(inputs_dir):
            makedirs(inputs_dir)

        make_deconstructsigs_patient_inputs(inputs_dir, variants)

        # count




