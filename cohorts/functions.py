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

from .variant_filters import no_filter, effect_expressed_filter
from .varcode_utils import FilterableVariant
from .utils import first_not_none_param
from .variant_stats import variant_stats_from_variant

from functools import wraps
import numpy as np
import pandas as pd
from varcode.effects import Substitution, FrameShift
from varcode.common import memoize
from varcode.effects.effect_classes import Exonic
import inspect

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

def count_variants_function_builder(function_name, filterable_variant_function=None):
    """
    Creates a function that counts variants that are filtered by the provided filterable_variant_function.
    The filterable_variant_function is a function that takes a filterable_variant and returns True or False.

    Users of this builder need not worry about applying e.g. the Cohort's default `filter_fn`. That will be applied as well.
    """
    @count_function
    def count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
        def count_filter_fn(filterable_variant, **kwargs):
            assert filter_fn is not None, "filter_fn should never be None, but it is."
            return ((filterable_variant_function(filterable_variant) if filterable_variant_function is not None else True) and
                    filter_fn(filterable_variant, **kwargs))
        patient_id = row["patient_id"]
        return cohort.load_variants(
            patients=[cohort.patient_from_id(patient_id)],
            filter_fn=count_filter_fn,
            **kwargs)
    count.__name__ = function_name
    count.__doc__ = str("".join(inspect.getsourcelines(filterable_variant_function)[0])) if filterable_variant_function is not None else ""
    return count

def count_effects_function_builder(function_name, only_nonsynonymous, filterable_effect_function=None):
    """
    Create a function that counts effects that are filtered by the provided filterable_effect_function.
    The filterable_effect_function is a function that takes a filterable_effect and returns True or False.

    Users of this builder need not worry about applying e.g. the Cohort's default `filter_fn`. That will be applied as well.
    """
    @count_function
    def count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
        def count_filter_fn(filterable_effect, **kwargs):
            assert filter_fn is not None, "filter_fn should never be None, but it is."
            return ((filterable_effect_function(filterable_effect) if filterable_effect_function is not None else True) and
                    filter_fn(filterable_effect, **kwargs))
        # This only loads one effect per variant.
        patient_id = row["patient_id"]
        return cohort.load_effects(
            only_nonsynonymous=only_nonsynonymous,
            patients=[cohort.patient_from_id(patient_id)],
            filter_fn=count_filter_fn,
            **kwargs)
    count.__name__ = function_name
    count.__doc__ = (("only_nonsynonymous=%s\n" % only_nonsynonymous) +
                     str("".join(inspect.getsourcelines(filterable_effect_function)[0])) if filterable_effect_function is not None else "")
    return count

variant_count = count_variants_function_builder("variant_count")

snv_count = count_variants_function_builder(
    "snv_count",
    filterable_variant_function=lambda filterable_variant: (
        filterable_variant.variant.is_snv))

indel_count = count_variants_function_builder(
    "indel_count",
    filterable_variant_function=lambda filterable_variant: (
        filterable_variant.variant.is_indel))

deletion_count = count_variants_function_builder(
    "deletion_count",
    filterable_variant_function=lambda filterable_variant: (
        filterable_variant.variant.is_deletion))

insertion_count = count_variants_function_builder(
    "insertion_count",
    filterable_variant_function=lambda filterable_variant: (
        filterable_variant.variant.is_insertion))

effect_count = count_effects_function_builder(
    "effect_count",
    only_nonsynonymous=False)

nonsynonymous_snv_count = count_effects_function_builder(
    "nonsynonymous_snv_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        filterable_effect.variant.is_snv))

missense_snv_count = count_effects_function_builder(
    "missense_snv_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        type(filterable_effect.effect) == Substitution and
        filterable_effect.variant.is_snv))

nonsynonymous_indel_count = count_effects_function_builder(
    "nonsynonymous_indel_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        filterable_effect.variant.is_indel))

nonsynonymous_deletion_count = count_effects_function_builder(
    "nonsynonymous_deletion_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        filterable_effect.variant.is_deletion))

nonsynonymous_insertion_count = count_effects_function_builder(
    "nonsynonymous_insertion_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        filterable_effect.variant.is_insertion))

exonic_variant_count = count_effects_function_builder(
    "exonic_variant_count",
    only_nonsynonymous=False,
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, Exonic)))

exonic_snv_count = count_effects_function_builder(
    "exonic_snv_count",
    only_nonsynonymous=False,
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, Exonic) and
        filterable_effect.variant.is_snv))

exonic_indel_count = count_effects_function_builder(
    "exonic_indel_count",
    only_nonsynonymous=False,
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, Exonic) and
        filterable_effect.variant.is_indel))

exonic_deletion_count = count_effects_function_builder(
    "exonic_deletion_count",
    only_nonsynonymous=False,
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, Exonic) and
        filterable_effect.variant.is_deletion))

exonic_insertion_count = count_effects_function_builder(
    "exonic_insertion_count",
    only_nonsynonymous=False,
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, Exonic) and
        filterable_effect.variant.is_insertion))

frameshift_count = count_effects_function_builder(
    "frameshift_count",
    only_nonsynonymous=False, # Should not matter, because FrameShift extends NonsilentCodingMutation
    filterable_effect_function=lambda filterable_effect: (
        isinstance(filterable_effect.effect, FrameShift)))

missense_snv_and_nonsynonymous_indel_count = count_effects_function_builder(
    "missense_snv_and_nonsynonymous_indel_count",
    only_nonsynonymous=True,
    filterable_effect_function=lambda filterable_effect: (
        (filterable_effect.variant.is_indel) or
         (type(filterable_effect.effect) == Substitution and
          filterable_effect.variant.is_snv)))

@count_function
def neoantigen_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    patient = cohort.patient_from_id(row["patient_id"])
    return cohort.load_neoantigens(patients=[patient],
                                   filter_fn=filter_fn,
                                   **kwargs)

@use_defaults
def expressed_missense_snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    def expressed_filter_fn(filterable_effect, **kwargs):
        assert filter_fn is not None, "filter_fn should never be None, but it is."
        return filter_fn(filterable_effect) and effect_expressed_filter(filterable_effect)
    return missense_snv_count(row=row,
                              cohort=cohort,
                              filter_fn=expressed_filter_fn,
                              normalized_per_mb=normalized_per_mb, **kwargs)

@use_defaults
def expressed_exonic_snv_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    def expressed_filter_fn(filterable_effect, **kwargs):
        assert filter_fn is not None, "filter_fn should never be None, but it is."
        return filter_fn(filterable_effect) and effect_expressed_filter(filterable_effect)
    return exonic_snv_count(row=row,
                              cohort=cohort,
                              filter_fn=expressed_filter_fn,
                              normalized_per_mb=normalized_per_mb, **kwargs)

@use_defaults
def expressed_neoantigen_count(row, cohort, filter_fn, normalized_per_mb, **kwargs):
    return neoantigen_count(row=row,
                            cohort=cohort,
                            filter_fn=filter_fn,
                            normalized_per_mb=normalized_per_mb,
                            only_expressed=True,
                            **kwargs)

def median_vaf_purity(row, cohort, **kwargs):
    """
    Estimate purity based on 2 * median VAF.

    Even if the Cohort has a default filter_fn, ignore it: we want to use all variants for
    this estimate.
    """
    patient_id = row["patient_id"]
    patient = cohort.patient_from_id(patient_id)
    variants = cohort.load_variants(patients=[patient], filter_fn=no_filter)[patient_id]
    def grab_vaf(variant):
        filterable_variant = FilterableVariant(variant, variants, patient)
        return variant_stats_from_variant(variant, filterable_variant.variant_metadata).tumor_stats.variant_allele_frequency
    vafs = [grab_vaf(variant) for variant in variants]
    return 2 * pd.Series(vafs).median()
