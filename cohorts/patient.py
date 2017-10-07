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

from .variant_filters import no_filter
from .varcode_utils import FilterableVariant
from .variant_stats import variant_stats_from_variant
from .utils import require_id_str, set_attributes
from collections import defaultdict, namedtuple
import pandas as pd
from varcode.effects import EffectCollection, effect_sort_key
from varcode.effects.effect_classes import KnownAminoAcidChange, FrameShift, Substitution
import numpy as np

FilterFnCol = namedtuple("FilterFnCol",
                         ["filter_fn", "only_nonsynonymous", "col_name"])

def build_filter_fn_col_from_count_fn(count_fn):
    if "_count" not in count_fn.__name__:
        raise ValueError("Cannot build a column from a count function without _count in its name: {}".format(count_fn.__name__))

    name_prefix = count_fn.__name__.split("_count")[0]
    return FilterFnCol(
        filter_fn=count_fn.filterable_effect_function,
        only_nonsynonymous=count_fn.only_nonsynonymous,
        col_name=name_prefix)

def build_filter_fn_cols_from_count_fns(count_fns):
    cols = []
    for count_fn in count_fns:
        cols.append(build_filter_fn_col_from_count_fn(count_fn=count_fn))
    return cols

transitions = set(
    ["C>T", "T>C", "G>A", "A>G"]
)
def is_transition(variant):
    change = variant.ref + ">" + variant.alt
    return change in transitions

transversions = set(
    ["C>G", "T>G", "C>A", "T>A", "G>C", "G>T", "A>C", "A>T"]
)
def is_transversion(variant):
    change = variant.ref + ">" + variant.alt
    return change in transversions

def min_five_percent_tumor_vaf_filter(filterable_variant):
    somatic_stats = variant_stats_from_variant(filterable_variant.variant,
                                               filterable_variant.variant_metadata)
    return (somatic_stats.tumor_stats.variant_allele_frequency >= 0.05)

class Patient(object):
    """
    Represents a patient, which contains zero or more of: normal `Sample`,
    tumor `Sample`.

    Parameters
    __________
    id : str
        ID of the patient.
    os : int
        Overall survival in days.
    pfs : int
        Progression-free survival in days.
    deceased : bool
        Is the patient deceased?
    progressed : bool
        Has the patient progressed?
    progressed_or_deceased : bool
        Has the patient either progressed or passed away?
    benefit : bool
        Has the patient seen a durable clinical benefit?
    variants : str or VariantCollection or list
        VCF or MAF path, VariantCollection, or list of some combination of those.
    normal_sample : Sample
        This patient's normal `Sample`.
    tumor_sample: Sample
        This patient's tumor `Sample`.
    hla_alleles : list
        A list of this patient's HLA class I alleles.
    additional_data : dict
        A dictionary of additional data: name of datum mapping to value.
        Will create these attributes in the Patient object.
    """
    def __init__(self,
                 id,
                 os,
                 pfs,
                 deceased,
                 progressed=None,
                 progressed_or_deceased=None,
                 benefit=None,
                 variants=[],
                 normal_sample=None,
                 tumor_sample=None,
                 hla_alleles=None,
                 additional_data=None,
                 cohort=None):
        require_id_str(id)
        self.id = id
        self.os = os
        self.pfs = pfs
        self.deceased = deceased
        self.progressed = progressed
        self.progressed_or_deceased = progressed_or_deceased
        self.benefit = benefit
        self.variants = variants
        self.normal_sample = normal_sample
        self.tumor_sample = tumor_sample
        self.hla_alleles = hla_alleles
        self.additional_data = additional_data

        if self.additional_data is not None:
            set_attributes(self, self.additional_data)

        # TODO: This can be removed once all patient-specific functions are
        # removed from Cohort.
        self.cohort = cohort

        self.add_progressed_or_deceased()

    @property
    def variants_list(self):
        return self.variants if type(self.variants) == list else [self.variants]

    def variants_dataframe(self, filter_fn_cols=[]):
        effects_sets = {}

        # Generate a set of unfiltered variants.
        variants_no_filter = self.cohort.load_variants(
            patients=[self],
            filter_fn=no_filter)[self.id]
        df_variants_no_filter = variants_no_filter.to_dataframe()
        assert len(df_variants_no_filter) == len(variants_no_filter.to_dataframe()[["chr", "start", "ref", "alt"]].drop_duplicates())

        # Also keep track of all effects.
        all_effects = self.cohort.load_effects(
                only_nonsynonymous=True,
                all_effects=True,
                patients=[self],
                filter_fn=no_filter)[self.id]
        variant_to_effects = all_effects.groupby_variant()

        # First, add depth/VAF columns.
        col_values = defaultdict(list)
        for variant in variants_no_filter:
            filterable_variant = FilterableVariant(variant, variants_no_filter, self.cohort)
            col_values["vaf_filter"].append(min_five_percent_tumor_vaf_filter(filterable_variant))
            col_values["transition"].append(is_transition(variant))
            col_values["transversion"].append(is_transversion(variant))
            somatic_stats = variant_stats_from_variant(filterable_variant.variant, filterable_variant.variant_metadata)
            if somatic_stats.tumor_stats is not None:
                col_values["tumor_vaf"].append(somatic_stats.tumor_stats.variant_allele_frequency)
                col_values["tumor_depth"].append(somatic_stats.tumor_stats.depth)
                col_values["tumor_alt_depth"].append(somatic_stats.tumor_stats.alt_depth)
            if somatic_stats.normal_stats is not None:
                col_values["normal_vaf"].append(somatic_stats.normal_stats.variant_allele_frequency)
                col_values["normal_depth"].append(somatic_stats.normal_stats.depth)
                col_values["normal_alt_depth"].append(somatic_stats.normal_stats.alt_depth)

            assert len(filterable_variant.variant_metadata.keys()) == 1
            inner_metadata_dict = filterable_variant.variant_metadata[list(filterable_variant.variant_metadata.keys())[0]]
            dbnsfp_pred = inner_metadata_dict.get("i_dbNSFP_LR_pred", ".") # Sometimes the key isn't even there
            if pd.isnull(dbnsfp_pred):
                dbnsfp_pred = "." # Or the key could be there, but value NaN
            col_values["possibly_deleterious"].append("D" in dbnsfp_pred)
            # Now look at effects.
            variant_effects = variant_to_effects.get(variant, EffectCollection([]))

            # Make sure the effect collection is sorted in priority order now that we're printing pieces of it out.
            # Also, throw in all alternate_effect's as well.
            expanded_effects = []
            for effect in variant_effects:
                expanded_effects.append(effect)
                if hasattr(effect, "alternate_effect"):
                    expanded_effects.append(effect.alternate_effect)
            variant_effects = EffectCollection(expanded_effects, sort_key=effect_sort_key, distinct=True)

            is_nonsynonymous = len(variant_effects.drop_silent_and_noncoding()) > 0
            col_values["nonsynonymous"].append(is_nonsynonymous)
            col_values["nonsynonymous_snv"].append(is_nonsynonymous and variant.is_snv)
            col_values["missense_snv"].append(is_nonsynonymous and variant.is_snv and (len(variant_effects.filter(lambda effect: type(effect) == Substitution)) > 0))
            col_values["nonsynonymous_indel"].append(is_nonsynonymous and variant.is_indel)
            col_values["frameshift"].append(is_nonsynonymous and (len(variant_effects.filter(lambda effect: isinstance(effect, FrameShift))) > 0))

            top_priority_effect = variant_effects.top_priority_effect() if len(variant_effects) > 0 else None
            col_values["priority_nonsynonymous_effect"].append(top_priority_effect.short_description if top_priority_effect is not None else "")
            col_values["priority_nonsynonymous_gene_id"].append(top_priority_effect.gene_id if top_priority_effect is not None else "")
            col_values["priority_nonsynonymous_gene_name"].append(top_priority_effect.gene_name if top_priority_effect is not None else "")

        # If this column has no information, remove it (e.g. not using MAFs, etc.).
        if set(col_values["possibly_deleterious"]) == set([False]):
            del col_values["possibly_deleterious"]

        # Put new values back into the unfiltered variant DataFrame.
        for col_name, col_value in col_values.items():
            if len(col_value) != len(df_variants_no_filter):
                raise ValueError("Trying to add a column that isn't the same length as the variants DataFrame: {} is {}".format(
                    col_name, len(col_value)))
            df_variants_no_filter[col_name] = col_value

        return df_variants_no_filter

    def add_progressed_or_deceased(self):
        assert self.progressed is not None or self.progressed_or_deceased is not None, (
            "Need at least one of progressed and progressed_or_deceased")
        if self.progressed_or_deceased is None:
            self.progressed_or_deceased = self.progressed or self.deceased

        # If we have both of these, ensure that they're in sync
        if self.progressed_or_deceased is not None and self.progressed is not None:
            assert self.progressed_or_deceased == self.progressed or self.deceased, (
                    "progressed_or_deceased should equal progressed || deceased")

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.__class__ == other.__class__ and
            self.id == other.id)

    def __str__(self):
        def str_if_not_none(s):
            return s if s is not None else "None"

        return "Patient(id=\"{}\", os={:.2f}, pfs={:.2f}, deceased={}, progressed={}, benefit={})".format(
            self.id,
            str_if_not_none(self.os),
            str_if_not_none(self.pfs),
            str_if_not_none(self.deceased),
            str_if_not_none(self.progressed),
            str_if_not_none(self.benefit))

    def __repr__(self):
        return str(self)
