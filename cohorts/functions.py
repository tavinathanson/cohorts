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

from __future__ import print_function

from .variant_filters import variant_qc_filter, effect_qc_filter, neoantigen_qc_filter

import numpy as np
from varcode.effects import Substitution
from varcode.common import memoize

def snv_count(row, cohort, filter_fn=variant_qc_filter,
              normalized_per_mb=True, **kwargs):
    patient_id = row["patient_id"]
    patient_variants = cohort.load_variants(
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=filter_fn,
        **kwargs)
    if patient_id in patient_variants:
        count = len(patient_variants[patient_id])
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan

def nonsynonymous_snv_count(row, cohort, filter_fn=effect_qc_filter,
                            normalized_per_mb=True, **kwargs):
    patient_id = row["patient_id"]
    patient_nonsynonymous_effects = cohort.load_effects(
        only_nonsynonymous=True,
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=filter_fn,
        **kwargs)
    if patient_id in patient_nonsynonymous_effects:
        count = len(patient_nonsynonymous_effects[patient_id])
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan

def missense_snv_count(row, cohort, filter_fn=effect_qc_filter,
                       normalized_per_mb=True, **kwargs):
    patient_id = row["patient_id"]
    def missense_filter_fn(effect, variant_metadata):
        if filter_fn is not None:
            return type(effect) == Substitution and filter_fn(effect, variant_metadata)
        return type(effect) == Substitution
    patient_missense_effects = cohort.load_effects(
        only_nonsynonymous=True,
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=missense_filter_fn,
        **kwargs)
    if patient_id in patient_missense_effects:
        count = len(patient_missense_effects[patient_id])
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan

def neoantigen_count(row, cohort, filter_fn=neoantigen_qc_filter,
                     normalized_per_mb=True, **kwargs):
    patient_id = row["patient_id"]
    patient = cohort.patient_from_id(row["patient_id"])
    patient_neoantigens = cohort.load_neoantigens(patients=[patient],
                                                  filter_fn=filter_fn,
                                                  **kwargs)
    if patient_id in patient_neoantigens:
        patient_neoantigens_df = patient_neoantigens[patient_id]
        count = len(patient_neoantigens_df)
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan

def expressed_missense_snv_count(row, cohort, **kwargs):
    return missense_snv_count(row, cohort, only_expressed=True, **kwargs)

def expressed_neoantigen_count(row, cohort, **kwargs):
    return neoantigen_count(row, cohort, only_expressed=True, **kwargs)

@memoize
def get_patient_to_mb(cohort):
    patient_to_mb = dict(cohort.as_dataframe(join_with="ensembl_coverage")[["patient_id", "MB"]].to_dict("split")["data"])
    return patient_to_mb
