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

from varcode import EffectCollection, Variant
from .errors import MissingBamFile
import logging

logger = logging.getLogger(__name__)


def genome(variant_collection):
    return variant_collection[0].ensembl

class FilterableVariant(object):
    def __init__(self, variant, variant_collection, patient):
        self.variant = variant
        self.variant_collection = variant_collection
        self.patient = patient

    @property
    def variant_metadata(self):
        source_to_metadata = self.variant_collection.source_to_metadata_dict
        return dict([
                (source, source_to_metadata[source][self.variant])
                for source in self.variant_collection.sources
                if self.variant in source_to_metadata[source]
            ])

    @property
    def genome(self):
        return genome(self.variant_collection)

class FilterableEffect(FilterableVariant):
    def __init__(self, effect, variant_collection, patient):
        self.effect = effect
        FilterableVariant.__init__(self,
                                   variant=effect.variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

class FilterableNeoantigen(FilterableVariant):
    def __init__(self, neoantigen_row, variant_collection, patient):
        self.neoantigen_row = neoantigen_row
        def build_variant(row, genome):
            return Variant(
                contig=row["chr"],
                ref=row["ref"],
                alt=row["alt"],
                start=row["start"],
                ensembl=genome)
        variant = build_variant(neoantigen_row, genome(variant_collection))
        FilterableVariant.__init__(self,
                                   variant=variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

class FilterablePolyphen(FilterableVariant):
    def __init__(self, polyphen_row, variant_collection, patient):
        self.polyphen_row = polyphen_row
        def build_variant(row, genome):
            return Variant(
                contig=row["chrom"],
                ref=row["ref"],
                alt=row["alt"],
                start=row["pos"],
                ensembl=genome)
        variant = build_variant(polyphen_row, genome(variant_collection))
        FilterableVariant.__init__(self,
                                   variant=variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

def filter_variants(variant_collection, patient, filter_fn, **kwargs):
    """Filter variants from the Variant Collection

    Parameters
    ----------
    variant_collection : varcode.VariantCollection
    patient : cohorts.Patient
    filter_fn: function
        Takes a FilterableVariant and returns a boolean. Only variants returning True are preserved.

    Returns
    -------
    varcode.VariantCollection
        Filtered variant collection, with only the variants passing the filter
    """
    if filter_fn:
        try:
            return variant_collection.clone_with_new_elements([
                variant
                for variant in variant_collection
                if filter_fn(FilterableVariant(
                        variant=variant,
                        variant_collection=variant_collection,
                        patient=patient,
                        ), **kwargs)
            ])
        except MissingBamFile as e:
            if patient.cohort.fail_on_missing_bams:
                raise
            else:
                logger.info(str(e))
                return None
    else:
        return variant_collection

def filter_effects(effect_collection, variant_collection, patient, filter_fn, all_effects, **kwargs):
    """Filter variants from the Effect Collection

    Parameters
    ----------
    effect_collection : varcode.EffectCollection
    variant_collection : varcode.VariantCollection
    patient : cohorts.Patient
    filter_fn : function
        Takes a FilterableEffect and returns a boolean. Only effects returning True are preserved.
    all_effects : boolean
        Return the single, top-priority effect if False. If True, return all effects (don't filter to top-priority).

    Returns
    -------
    varcode.EffectCollection
        Filtered effect collection, with only the variants passing the filter
    """
    def top_priority_maybe(effect_collection):
        """
        Always (unless all_effects=True) take the top priority effect per variant
        so we end up with a single effect per variant.
        """
        if all_effects:
            return effect_collection
        return EffectCollection(list(effect_collection.top_priority_effect_per_variant().values()))

    def apply_filter_fn(filter_fn, effect):
        """
        Return True if filter_fn is true for the effect or its alternate_effect.
        If no alternate_effect, then just return True if filter_fn is True.
        """
        applied = filter_fn(FilterableEffect(
            effect=effect,
            variant_collection=variant_collection,
            patient=patient), **kwargs)
        if hasattr(effect, "alternate_effect"):
            applied_alternate = filter_fn(FilterableEffect(
                effect=effect.alternate_effect,
                variant_collection=variant_collection,
                patient=patient), **kwargs)
            return applied or applied_alternate
        return applied

    if filter_fn:
        try:
            return top_priority_maybe(EffectCollection([
                effect
                for effect in effect_collection
                if apply_filter_fn(filter_fn, effect)]))
        except MissingBamFile as e:
            logger.warning(str(e))
            return None
    else:
        return top_priority_maybe(effect_collection)

def filter_neoantigens(neoantigens_df, variant_collection, patient, filter_fn):
    if filter_fn:
        try:
            filter_mask = neoantigens_df.apply(
                lambda row: filter_fn(
                    FilterableNeoantigen(neoantigen_row=row,
                                         variant_collection=variant_collection,
                                         patient=patient)),
                axis=1,
                # reduce ensures that an empty result is a Series vs. a DataFrame
                reduce=True)
            return neoantigens_df[filter_mask]
        except MissingBamFile as e:
            logger.warning(str(e))
            return None
    else:
        return neoantigens_df

def filter_polyphen(polyphen_df, variant_collection, patient, filter_fn):
    if filter_fn:
        filter_mask = polyphen_df.apply(
            lambda row: filter_fn(
                FilterablePolyphen(polyphen_row=row,
                                   variant_collection=variant_collection,
                                   patient=patient)),
            axis=1,
            # reduce ensures that an empty result is a Series vs. a DataFrame
            reduce=True)
        return polyphen_df[filter_mask]
    else:
        return polyphen_df
