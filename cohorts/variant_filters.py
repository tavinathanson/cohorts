# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from .variant_stats import variant_stats_from_variant

import numpy as np
from varcode import Substitution, Variant
import re
import pandas as pd
from os import path

DEFAULT_DEPTH = 30
DEFAULT_NORMAL_VAF = 0.02
DEFAULT_ALT_DEPTH = 5

def variant_qc_filter(variant, variant_metadata, depth=DEFAULT_DEPTH,
                      normal_vaf=DEFAULT_NORMAL_VAF, alt_depth=DEFAULT_ALT_DEPTH):
    somatic_stats = variant_stats_from_variant(variant, variant_metadata)

    # Filter variant with depth < 30
    if (somatic_stats.tumor_stats.depth < depth or
        somatic_stats.normal_stats.depth < depth):
        return False

    # Filter based on normal evidence
    if somatic_stats.normal_stats.variant_allele_frequency > normal_vaf:
        return False

    if somatic_stats.tumor_stats.alt_depth < alt_depth:
        return False

    return True

def effect_qc_filter(effect, variant_metadata, **kwargs):
    return variant_qc_filter(effect.variant, variant_metadata, **kwargs)

def neoantigen_qc_filter(row, variants, **kwargs):
    # We need to re-create the Variant object from the neoeantigen DataFrame.
    # TODO: Just pickle the Neoantigen collection rather than dealing with
    # DataFrames.
    genome = variants[0].ensembl
    variant = Variant(contig=row["chr"],
                      ref=row["ref"],
                      alt=row["alt"],
                      start=row["start"],
                      ensembl=genome)
    return variant_qc_filter(variant, variants.metadata[variant], **kwargs)

def polyphen_qc_filter(row, variants, **kwargs):
    genome = variants[0].ensembl
    variant = Variant(contig=row["chrom"],
                      ref=row["ref"],
                      alt=row["alt"],
                      start=row["pos"],
                      ensembl=genome)
    return variant_qc_filter(variant, variants.metadata[variant], **kwargs)

def load_ensembl_coverage(cohort, coverage_path, min_depth=30):
    """
    Load in Pageant CoverageDepth results with Ensembl loci.

    coverage_path is a path to Pageant CoverageDepth output directory, with
    one subdirectory per patient and a `cdf.csv` file inside each patient subdir.
    """
    columns = [
        "NormalDepth",
        "TumorDepth",
        "Normal BP",
        "Tumor BP",
        "Num Loci",
        "% Normal BP",
        "% Tumor BP",
        "% Loci",
        "Off-target Normal BP",
        "Off-target Tumor BP",
        "Off-target Num Loci",
        "Off-target % Normal BP",
        "Off-target % Tumor BP",
        "Off-target % Loci",
    ]
    ensembl_loci_dfs = []
    for patient in cohort:
        patient_ensembl_loci_df = pd.read_csv(
            path.join(coverage_path, patient.id, "cdf.csv"),
            names=columns)
        # pylint: disable=no-member
        # pylint gets confused by read_csv
        patient_ensembl_loci_df = patient_ensembl_loci_df[(
            (patient_ensembl_loci_df.NormalDepth == min_depth) &
            (patient_ensembl_loci_df.TumorDepth == min_depth))]
        assert len(patient_ensembl_loci_df) == 1, (
            "Incorrect number of %d depth loci results: %d" % (
                min_depth, len(patient_ensembl_loci_df)))
        patient_ensembl_loci_df["patient_id"] = patient.id
        ensembl_loci_dfs.append(patient_ensembl_loci_df)
    ensembl_loci_df = pd.concat(ensembl_loci_dfs)
    ensembl_loci_df["MB"] = ensembl_loci_df["Num Loci"] / 1000000.0
    return ensembl_loci_df[["patient_id", "Num Loci", "MB"]]
