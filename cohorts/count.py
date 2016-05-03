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

def snv_count(cohort, **kwargs):
    patient_variants = cohort.load_variants(**kwargs)
    def count_func(patient_id):
        if patient_id in patient_variants:
            return len(patient_variants[patient_id])
        return np.nan
    return count(cohort, count_func, count_col="snv_count")

def nonsynonymous_snv_count(cohort, **kwargs):
    patient_nonsynonymous_effects = cohort.load_effects(only_nonsynonymous=True, **kwargs)
    def count_func(patient_id):
        if patient_id in patient_nonsynonymous_effects:
            return len(patient_nonsynonymous_effects[patient_id])
        return np.nan
    return count(cohort, count_func, count_col="nonsynonymous_snv_count")

def missense_snv_count(cohort, **kwargs):
    patient_nonsynonymous_effects = cohort.load_effects(only_nonsynonymous=True, **kwargs)
    patient_missense_effects = dict(
        [(patient_id,
          EffectCollection(
              [effect for effect in effects if type(effect) == Substitution]))
         for (patient_id, effects) in patient_nonsynonymous_effects.items()])
    def count_func(patient_id):
        if patient_id in patient_missense_effects:
            return len(patient_missense_effects[patient_id])
        return np.nan
    return count(cohort, count_func, count_col="missense_snv_count")

def neoantigen_count(cohort, **kwargs):
    patient_neoantigens = cohort.load_neoantigens(**kwargs)
    def count_func(patient_id):
        if patient_id in patient_neoantigens["patient_id"].unique():
            return len(patient_neoantigens[patient_neoantigens["patient_id"] == patient_id])
        return np.nan
    return count(cohort, count_func, count_col="neoantigen_count")

def count(cohort, count_func, count_col):
    df = cohort.as_dataframe()
    df[count_col] = df["patient_id"].map(count_func)
    original_len = len(df)
    df = df[~df[count_col].isnull()]
    updated_len = len(df)
    if updated_len < original_len:
        print("Missing count for %d patients: from %d to %d" % (original_len - updated_len, original_len, updated_len))
    return count_col, df
