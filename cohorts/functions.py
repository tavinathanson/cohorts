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
from .utils import first_not_none_param, get_logger
from .variant_stats import variant_stats_from_variant

from functools import wraps
import numpy as np
import pandas as pd
from varcode.effects import Substitution, FrameShift
from varcode.common import memoize
from varcode.effects.effect_classes import Exonic, StopLoss
import inspect

logger = get_logger(__name__)

def use_defaults(func):
    """
    Decorator for functions that should automatically fall back to the Cohort-default filter_fn and
    normalized_per_mb if not specified.
    """
    @wraps(func)
    def wrapper(row, cohort, filter_fn=None, normalized_per_mb=None, **kwargs):
        logger.debug("applying defaults")
        filter_fn = first_not_none_param([filter_fn, cohort.filter_fn], no_filter)
        normalized_per_mb = first_not_none_param([normalized_per_mb, cohort.normalized_per_mb], False)
        logger.debug("filter_fn set to: {}".format(str(filter_fn.__name__)))
        logger.debug("filter_fn has content: {}".format(str(filter_fn)))
        logger.debug("normalized_per_mb set to: {}".format(str(normalized_per_mb)))
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
    # Keep track of these to be able to query the returned function for these attributes
    count.only_nonsynonymous = only_nonsynonymous
    count.filterable_effect_function = filterable_effect_function
    return count

def create_effect_filter(effect_name, effect_filter):
    """
    Return a composable filter function that applies the provided effect_filter 
    to variants/effects of the provided function.

    For example:
        only_snv = create_effect_filter("snv", is_snv)
        snv_variant_count = only_snv(variant_count)
        snv_read_count = only_snv(read_count)

    Parameters
    ---------
    effect_name: (str)
      description of the filter applied. Used to construct a default "__doc__" & "__name__" 
      of the resulting function.
    effect_filter: (function)
      function taking a single input (filterable_effect) & returning a boolean
      only effects for which the filter returns True will be utilized.

    Returns
    -------
    A function that takes two parameters: func & (optionally) name
    """
    def filter_by_type(func, name=None):
        @use_defaults
        def filtered(row, cohort, filter_fn, normalized_per_mb, **kwargs):
            def new_filter_fn(filterable_effect, **kwargs):
                assert filter_fn is not None, "filter_fn should never be None, but it is."
                return filter_fn(filterable_effect, **kwargs) and effect_filter(filterable_effect)
            return func(row=row,
                        cohort=cohort,
                        filter_fn=new_filter_fn,
                        normalized_per_mb=normalized_per_mb,
                        **kwargs)
        if name is None:
            name = "_".join([effect_name, func.__name__])
        filtered.__name__ = name
        filtered.__doc__ = func.__doc__ + "; Only {} effects.".format(effect_name)
        return filtered
    filter_by_type.__doc__ = "Return a new count function limited to {} effects".format(effect_name)
    return filter_by_type

def is_exonic(filterable_effect):
    return isinstance(filterable_effect.effect, Exonic)
def is_frameshift(filterable_effect):
    return isinstance(filterable_effect.effect, FrameShift)
def is_snv(filterable_effect):
    return filterable_effect.variant.is_snv
def is_indel(filterable_effect):
    return filterable_effect.variant.is_indel
def is_missense(filterable_effect):
    return type(filterable_effect.effect) == Substitution
def is_insertion(filterable_effect):
    return filterable_effect.variant.is_insertion
def is_deletion(filterable_effect):
    return filterable_effect.variant.is_deletion
def is_stoploss(filterable_effect):
    return isinstance(filterable_effect.effect, StopLoss)
def is_expressed(filterable_effect):
    return effect_expressed_filter(filterable_effect)
def is_nonsynonymous(filterable_effect):
    return filterable_effect.effect.modifies_protein_sequence

only_exonic = create_effect_filter("exonic", is_exonic)
only_frameshift = create_effect_filter("frameshift", is_frameshift)
only_snv = create_effect_filter("snv", is_snv)
only_indel = create_effect_filter("indel", is_indel)
only_missense = create_effect_filter("missense", is_missense)
only_insertion = create_effect_filter("insertion", is_insertion)
only_deletion = create_effect_filter("deletion", is_deletion)
only_stoploss = create_effect_filter("stoploss", is_stoploss)
only_expressed = create_effect_filter("expressed", is_expressed)
only_nonsynonymous = create_effect_filter("nonsynonymous", is_nonsynonymous)

## base count functions
effect_count = count_effects_function_builder("effect_count", only_nonsynonymous=False)
nonsynonymous_effect_count = count_effects_function_builder("nonsynonymous_effect_count", only_nonsynonymous=True)

# for effect counts, make equivalent functions named "variants" (for backwards compat)
nonsynonymous_variant_count = nonsynonymous_effect_count
nonsynonymous_variant_count.__name__ = "nonsynonymous_variant_count"
variant_count = effect_count
variant_count.__name__ = "variant_count"

# main count functions, with custom names
insertion_count = only_insertion(variant_count, name="insertion_count")
snv_count = only_snv(variant_count, name="snv_count")
indel_count = only_indel(variant_count, name="indel_count")
deletion_count = only_deletion(variant_count, name="deletion_count")

# other common count functions
missense_snv_count = only_missense(snv_count)
nonsynonymous_snv_count = only_nonsynonymous(snv_count)
nonsynonymous_indel_count = only_nonsynonymous(indel_count)
nonsynonymous_deletion_count = only_nonsynonymous(deletion_count)
nonsynonymous_insertion_count = only_nonsynonymous(insertion_count)

# exonic counts
exonic_snv_count = only_exonic(snv_count)
exonic_variant_count = only_exonic(variant_count)
exonic_missense_snv_count = only_exonic(missense_snv_count)
exonic_indel_count = only_exonic(indel_count)
exonic_deletion_count = only_exonic(deletion_count)
exonic_insertion_count = only_exonic(insertion_count)
exonic_frameshift_deletion_count = only_exonic(only_frameshift(deletion_count))
exonic_frameshift_insertion_count = only_exonic(only_frameshift(insertion_count))
exonic_frameshift_indel_count = only_exonic(only_frameshift(indel_count))

# expressed counts
expressed_missense_snv_count = only_expressed(missense_snv_count)
expressed_exonic_variant_count = only_expressed(exonic_variant_count)
expressed_exonic_snv_count = only_expressed(exonic_snv_count)
expressed_exonic_indel_count = only_expressed(exonic_indel_count)
expressed_exonic_insertion_count = only_expressed(exonic_insertion_count)
expressed_exonic_deletion_count = only_expressed(exonic_deletion_count)

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
    variants = cohort.load_variants(patients=[patient], filter_fn=no_filter)
    if patient_id in variants.keys():
        variants = variants[patient_id]
    else:
        return np.nan
    def grab_vaf(variant):
        filterable_variant = FilterableVariant(variant, variants, patient)
        return variant_stats_from_variant(variant, filterable_variant.variant_metadata).tumor_stats.variant_allele_frequency
    vafs = [grab_vaf(variant) for variant in variants]
    return 2 * pd.Series(vafs).median()
