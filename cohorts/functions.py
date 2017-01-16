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

from .variant_filters import no_filter, effect_expressed_filter
from .utils import first_not_none_param

from functools import wraps
import numpy as np
from varcode.effects import Substitution
from varcode.common import memoize

def use_defaults(func):
    """
    Decorator for functions that should automatically fall back to the Cohort-default filter_fn and
    normalized_per_mb if not specified.
    """
    @wraps(func)
    def wrapper(row, cohort, filter_fn=None, normalized_per_mb=None, **kwargs):
        filter_fn = first_not_none_param([filter_fn, cohort.filter_fn], no_filter)
        normalized_per_mb = first_not_none_param([normalized_per_mb, cohort.normalized_per_mb], False)
        return func(row=row,
                    cohort=cohort,
                    filter_fn=filter_fn,
                    normalized_per_mb=normalized_per_mb,
                    **kwargs)
    return wrapper

def count_function(func):
    """
    Decorator for functions that return a collection (technically a dict of collections) that should be
    counted up. Also automatically falls back to the Cohort-default filter_fn and normalized_per_mb if
    not specified.
    """
    # Fall back to Cohort-level defaults.
    @use_defaults
    @wraps(func)
    def wrapper(row, cohort, filter_fn=None, normalized_per_mb=None, **kwargs):
        per_patient_data = func(row=row,
                                cohort=cohort,
                                filter_fn=filter_fn,
                                normalized_per_mb=normalized_per_mb,
                                **kwargs)
        patient_id = row["patient_id"]
        if patient_id in per_patient_data:
            count = len(per_patient_data[patient_id])
            if normalized_per_mb:
                count /= float(get_patient_to_mb(cohort)[patient_id])
            return count
        return np.nan
    return wrapper

@memoize
def get_patient_to_mb(cohort):
    patient_to_mb = dict(cohort.as_dataframe(join_with="ensembl_coverage")[["patient_id", "MB"]].to_dict("split")["data"])
    return patient_to_mb

@count_function
def snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    patient_id = row["patient_id"]
    return cohort.load_variants(
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=filter_fn,
        **kwargs)

@count_function
def nonsynonymous_snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    # This only loads one effect per variant.
    patient_id = row["patient_id"]
    return cohort.load_effects(
        only_nonsynonymous=True,
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=filter_fn,
        **kwargs)

@count_function
def missense_snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    def missense_filter_fn(filterable_effect):
        assert filter_fn is not None, "filter_fn should never be None, but it is."
        return (type(filterable_effect.effect) == Substitution and
                filter_fn(filterable_effect))
    # This only loads one effect per variant.
    patient_id = row["patient_id"]
    return cohort.load_effects(
        only_nonsynonymous=True,
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=missense_filter_fn,
        **kwargs)

@count_function
def neoantigen_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    patient = cohort.patient_from_id(row["patient_id"])
    return cohort.load_neoantigens(patients=[patient],
                                   filter_fn=filter_fn,
                                   **kwargs)

@use_defaults
def expressed_missense_snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    def expressed_filter_fn(filterable_effect):
        assert filter_fn is not None, "filter_fn should never be None, but it is."
        return filter_fn(filterable_effect) and effect_expressed_filter(filterable_effect)
    return missense_snv_count(row=row,
                              cohort=cohort,
                              filter_fn=expressed_filter_fn,
                              normalized_per_mb=normalized_per_mb)

@use_defaults
def expressed_neoantigen_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    return neoantigen_count(row=row,
                            cohort=cohort,
                            filter_fn=filter_fn,
                            normalized_per_mb=normalized_per_mb,
                            only_expressed=True,
                            **kwargs)
