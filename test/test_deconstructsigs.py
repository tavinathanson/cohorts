from __future__ import print_function

from . import data_path, generated_data_path
from .data_generate import generate_vcfs
from .test_count import make_cohort
from cohorts.signatures.deconstructsigs_utils \
    import make_deconstructsigs_patient_inputs, run_deconstructsigs

from shutil import rmtree
from os import path, makedirs, listdir
from nose.tools import eq_, assert_not_equal
import pandas as pd
import numpy as np

FILE_FORMAT_1 = "patient_format1_%s.vcf"

def test_deconstructsigs_output():
    """
    Run on a subset (3 patients)
    Test that output csv is non-empty
    """

    try:
        input_dir = data_path("signatures")
        sample_id_list = ["1", "2", "3"]
        output_dir = generated_data_path("deconstructsigs-outputs")
        path_to_write_r_script = generated_data_path("tmp_r")
        if not path.exists(output_dir):
            makedirs(output_dir)

        run_deconstructsigs(input_dir,
                            sample_id_list,
                            output_dir=output_dir,
                            path_to_write_r_script=path_to_write_r_script
                            )
        signatures_path = "{}/cohort_signatures.csv".format(output_dir)
        variant_counts_path = "{}/signature_variant_counts.csv".format(output_dir)

        # check that vals in the csvs are non-zero
        sigs_matrix = _sigs_csv_to_matrix(signatures_path, sample_id_list)
        sum = np.sum(sigs_matrix)
        assert_not_equal(sum, 0)

        variant_counts = _sigs_csv_to_matrix(variant_counts_path, sample_id_list)
        counts_sum = np.sum(variant_counts)
        assert_not_equal(counts_sum, 0)

    finally:
        if output_dir is not None and path.exists(output_dir):
            rmtree(output_dir)

def _sigs_csv_to_matrix(path, samples_ordering):
    df = pd.read_csv(path)
    return df[df.columns[1:-1]].values

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
