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

def variant_qc_filter(variant, variant_metadata):
    somatic_stats = variant_stats_from_variant(variant, variant_metadata)

    # Filter variant with depth < 30
    if somatic_stats.tumor_stats.depth < 30 or somatic_stats.normal_stats.depth < 30:
        return False

    # Filter based on normal evidence
    if somatic_stats.normal_stats.variant_allele_frequency > .02:
        return False

    if somatic_stats.tumor_stats.alt_depth < 5:
        return False

    return True

def effect_qc_filter(effect, variant_metadata):
    return variant_qc_filter(effect.variant, variant_metadata)

def neoantigen_qc_filter(row, metadata):
    variant = variant_string_to_variant(row["variant"])
    return variant_qc_filter(variant, metadata[variant])

# TODO: Remove this hack and properly store the objects themselves.
genomic_coord_regex = r"g\.([0-9]+)([A|C|T|G])>([A|C|T|G])"
def variant_string_to_variant(variant_str, reference="grch37"):
    chrom, coord = variant_str.split()
    m = re.match(genomic_coord_regex, coord)
    start = m.group(1)
    ref = m.group(2)
    alt = m.group(3)
    variant = Variant(contig=chrom, start=start, ref=ref, alt=alt, ensembl=reference)
    return variant
