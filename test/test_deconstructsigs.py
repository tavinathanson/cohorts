from __future__ import print_function

from . import data_path, generated_data_path
from .data_generate import generate_vcfs
from .test_count import make_cohort
from cohorts.signatures.deconstructsigs_utils \
    import make_deconstructsigs_patient_inputs, run_deconstructsigs

from shutil import rmtree
from os import path, makedirs, listdir
from nose.tools import raises, eq_, ok_
import pandas as pd

FILE_FORMAT_1 = "patient_format1_%s.vcf"

def test_deconstructsigs_output():
    """
    Run on a subset (3 patients)
    Test that output csv is non-empty
    """
    return

def test_input_generation():
    """
    A test of make_deconstructsigs_patient_inputs
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        variants = cohort.load_variants() # dict of patient to variantcollection

        truth_file_lengths_dict= { "{}.input".format(patient_id) : len(variant_collection) \
                for patient_id, variant_collection in variants.iteritems()
            }

        inputs_path = generated_data_path("deconstructsigs-inputs")
        inputs_dir = path.dirname(inputs_path)
        if not path.exists(inputs_dir):
            makedirs(inputs_dir)

        make_deconstructsigs_patient_inputs(inputs_path, variants)

        # count the files and their lengths that were generated
        files_list = listdir(inputs_path)
        num_input_files = len(files_list)
        eq_(num_input_files, 3)

        file_to_length_dict = {}
        for input_file in files_list:
            df = pd.read_csv("{}/{}".format(inputs_path, input_file))
            file_to_length_dict[input_file] = len(df)

        eq_(file_to_length_dict, truth_file_lengths_dict)

    finally:
        if inputs_dir is not None and path.exists(inputs_dir):
            rmtree(inputs_dir)
        if cohort is not None:
            cohort.clear_caches()
