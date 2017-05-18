# Copyright (c) 2017. Mount Sinai School of Medicine
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
from .utils import get_logger

from varcode import Variant
from varcode.common import memoize
import pandas as pd
from os import path

logger = get_logger(__name__)

def no_filter(filterable_variant):
    return True

def variant_qc_filter(filterable_variant,
                      min_tumor_depth,
                      min_normal_depth,
                      min_tumor_vaf,
                      max_normal_vaf,
                      min_tumor_alt_depth):
    logger.debug('Applying variant_qc_filter with params: min_tumor_depth={}, min_normal_depth={}, min_tumor_vaf={}, max_normal_vaf={}, min_tumor_alt_depth={}'.format(min_tumor_depth, min_normal_depth, min_tumor_vaf, max_normal_vaf, min_tumor_alt_depth))

    somatic_stats = variant_stats_from_variant(filterable_variant.variant,
                                               filterable_variant.variant_metadata)

    # Filter variant with depth < depth
    if (somatic_stats.tumor_stats.depth < min_tumor_depth or
        somatic_stats.normal_stats.depth < min_normal_depth):
        return False

    # Filter based on normal evidence
    if somatic_stats.normal_stats.variant_allele_frequency > max_normal_vaf:
        return False

    if somatic_stats.tumor_stats.variant_allele_frequency < min_tumor_vaf:
        return False

    if somatic_stats.tumor_stats.alt_depth < min_tumor_alt_depth:
        return False

    return True

@memoize
def expressed_variant_set(cohort, patient, variant_collection):
    # Warning: we previously had an issue where we used the same
    # memoized set across different cohorts.
    # TODO: we're currently using the same isovar cache that we use for expressed
    # neoantigen prediction; so we pass in the same epitope lengths.
    # This is hacky and should be addressed.
    df_isovar = patient.cohort.load_single_patient_isovar(
        patient=patient,
        variants=variant_collection,
        epitope_lengths=[8, 9, 10, 11])
    expressed_variant_set = set()
    for _, row in df_isovar.iterrows():
        expressed_variant = Variant(contig=row["chr"],
                          start=row["pos"],
                          ref=row["ref"],
                          alt=row["alt"],
                          ensembl=variant_collection[0].ensembl)
        expressed_variant_set.add(expressed_variant)
    return expressed_variant_set

def variant_expressed_filter(filterable_variant, **kwargs):
    expressed_variants = expressed_variant_set(
        cohort=filterable_variant.patient.cohort,
        patient=filterable_variant.patient,
        variant_collection=filterable_variant.variant_collection)
    return filterable_variant.variant in expressed_variants

def effect_expressed_filter(filterable_effect, **kwargs):
    return variant_expressed_filter(filterable_effect, **kwargs)

def load_ensembl_coverage(cohort, coverage_path, min_tumor_depth, min_normal_depth=0,
                          pageant_dir_fn=None):
    """
    Load in Pageant CoverageDepth results with Ensembl loci.

    coverage_path is a path to Pageant CoverageDepth output directory, with
    one subdirectory per patient and a `cdf.csv` file inside each patient subdir.

    If min_normal_depth is 0, calculate tumor coverage. Otherwise, calculate
    join tumor/normal coverage.

    pageant_dir_fn is a function that takes in a Patient and produces a Pageant
    dir name.

    Last tested with Pageant CoverageDepth version 1ca9ed2.
    """
    # Function to grab the pageant file name using the Patient
    if pageant_dir_fn is None:
        pageant_dir_fn = lambda patient: patient.id

    columns_both = [
        "depth1", # Normal
        "depth2", # Tumor
        "onBP1",
        "onBP2",
        "numOnLoci",
        "fracBPOn1",
        "fracBPOn2",
        "fracLociOn",
        "offBP1",
        "offBP2",
        "numOffLoci",
        "fracBPOff1",
        "fracBPOff2",
        "fracLociOff",
    ]
    columns_single = [
        "depth",
        "onBP",
        "numOnLoci",
        "fracBPOn",
        "fracLociOn",
        "offBP",
        "numOffLoci",
        "fracBPOff",
        "fracLociOff"
    ]
    if min_normal_depth < 0:
        raise ValueError("min_normal_depth must be >= 0")
    use_tumor_only = (min_normal_depth == 0)
    columns = columns_single if use_tumor_only else columns_both
    ensembl_loci_dfs = []
    for patient in cohort:
        patient_ensembl_loci_df = pd.read_csv(
            path.join(coverage_path, pageant_dir_fn(patient), "cdf.csv"),
            names=columns,
            header=1)
        # pylint: disable=no-member
        # pylint gets confused by read_csv
        if use_tumor_only:
            depth_mask = (patient_ensembl_loci_df.depth == min_tumor_depth)
        else:
            depth_mask = (
                (patient_ensembl_loci_df.depth1 == min_normal_depth) &
                (patient_ensembl_loci_df.depth2 == min_tumor_depth))
        patient_ensembl_loci_df = patient_ensembl_loci_df[depth_mask]
        assert len(patient_ensembl_loci_df) == 1, (
            "Incorrect number of tumor={}, normal={} depth loci results: {} for patient {}".format(
                min_tumor_depth, min_normal_depth, len(patient_ensembl_loci_df), patient))
        patient_ensembl_loci_df["patient_id"] = patient.id
        ensembl_loci_dfs.append(patient_ensembl_loci_df)
    ensembl_loci_df = pd.concat(ensembl_loci_dfs)
    ensembl_loci_df["MB"] = ensembl_loci_df.numOnLoci / 1000000.0
    return ensembl_loci_df[["patient_id", "numOnLoci", "MB"]]
