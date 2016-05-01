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

from .load import verify_group_by

def snv_count(cohort, **kwargs):
    sample_variants = cohort.load_variants(**kwargs)
    def count_func(sample):
        if sample in sample_variants:
            return len(sample_variants[sample])
        return np.nan
    return count(cohort, count_func, count_col="snv_count")

def nonsynonymous_snv_count(cohort, **kwargs):
    sample_nonsynonymous_effects = cohort.load_effects(only_nonsynonymous=True, **kwargs)
    def count_func(sample):
        if sample in sample_nonsynonymous_effects:
            return len(sample_nonsynonymous_effects[sample])
        return np.nan
    return count(cohort, count_func, count_col="nonsynonymous_snv_count")

def missense_snv_count(cohort, **kwargs):
    sample_nonsynonymous_effects = cohort.load_effects(only_nonsynonymous=True, **kwargs)
    sample_missense_effects = dict(
        [(sample,
          EffectCollection(
              [effect for effect in effects if type(effect) == Substitution]))
         for (sample, effects) in sample_nonsynonymous_effects.items()])
    def count_func(sample):
        if sample in sample_missense_effects:
            return len(sample_missense_effects[sample])
        return np.nan
    return count(cohort, count_func, count_col="missense_snv_count")

def neoantigen_count(cohort, **kwargs):
    sample_neoantigens = cohort.load_neoantigens(**kwargs)
    def count_func(sample):
        if sample in sample_neoantigens["sample_id"].unique():
            return len(sample_neoantigens[sample_neoantigens["sample_id"] == sample])
        return np.nan
    return count(cohort, count_func, count_col="neoantigen_count")

def count(cohort, count_func, count_col, group_by="patient"):
    verify_group_by(group_by)

    id_col = "patient_id" if group_by == "patient" else "sample_id"

    df = cohort.as_dataframe(group_by=group_by)
    df[count_col] = df[id_col].map(count_func)
    original_len = len(df)
    df = df[~df[count_col].isnull()]
    updated_len = len(df)
    if updated_len < original_len:
        print("Missing count for %d samples: from %d to %d" % (original_len - updated_len, original_len, updated_len))
    return count_col, df
