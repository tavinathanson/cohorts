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

import numpy as np
import pandas as pd

from varcode import EffectCollection
from varcode.effects import Substitution

def snv_count(row, cohort, **kwargs):
    patient_id = row["patient_id"]
    patient_variants = cohort.load_variants(
        patients=[cohort.patient_from_id(patient_id)], **kwargs)
    if patient_id in patient_variants:
        return len(patient_variants[patient_id])
    return np.nan

def nonsynonymous_snv_count(row, cohort, **kwargs):
    patient_id = row["patient_id"]
    patient_nonsynonymous_effects = cohort.load_effects(
        only_nonsynonymous=True, patients=[cohort.patient_from_id(patient_id)], **kwargs)
    if patient_id in patient_nonsynonymous_effects:
        return len(patient_nonsynonymous_effects[patient_id])
    return np.nan

def missense_snv_count(row, cohort, **kwargs):
    patient_id = row["patient_id"]
    patient_missense_effects = cohort.load_effects(
        only_nonsynonymous=True, 
        patients=[cohort.patient_from_id(patient_id)], 
        filter_fn=lambda effect, variant_metadata: type(effect) == Substitution,
        **kwargs)
    if patient_id in patient_missense_effects:
        return len(patient_missense_effects[patient_id])
    return np.nan

def neoantigen_count(row, cohort, **kwargs):
    patient_id = row["patient_id"]
    patient_neoantigens = cohort.load_neoantigens(patients=[cohort.patient_from_id(patient_id)], **kwargs)
    if patient_id in patient_neoantigens["patient_id"].unique():
        return len(patient_neoantigens[patient_neoantigens["patient_id"] == patient_id])
    return np.nan

def expressed_neoantigen_count(row, cohort, **kwargs):
    return neoantigen_count(row, cohort, only_expressed=True)
