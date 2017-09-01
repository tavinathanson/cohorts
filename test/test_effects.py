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

from . import generated_data_path
from cohorts import Patient, Cohort
from varcode import Variant, VariantCollection
from varcode.effects import Substitution, FrameShift, IntronicSpliceSite
import pandas as pd
from nose.tools import eq_, ok_

def test_splice_filtering_substitution():
    """
    Make sure that ExonicSpliceSite mutations with Substitution alternates are kept even when we filter to
    only Substitution effects.
    """
    cohort = None
    try:
        variant = Variant(contig=10, start=124340409, ref="C", alt="A", ensembl=75)
        patient = Patient(id="patient", os=3, pfs=2, deceased=False, progressed=False, variants=VariantCollection([variant]))
        cohort_cache_path = generated_data_path("cache")
        cohort = Cohort(
            patients=[patient],
            cache_dir=cohort_cache_path)

        def missense_snv_filter(filterable_effect):
            return (type(filterable_effect.effect) == Substitution and
                    filterable_effect.variant.is_snv)

        effects = cohort.load_effects(filter_fn=missense_snv_filter)[patient.id]
        eq_(len(effects), 1)
        all_effects = cohort.load_effects(filter_fn=missense_snv_filter,
                                          all_effects=True)[patient.id]
        eq_(len(all_effects), 7)
    finally:
        if cohort is not None:
            cohort.clear_caches()

def test_splice_filtering_frameshift():
    """
    Make sure that ExonicSpliceSite mutations with FrameShift alternates are kept even when we filter to
    only FrameShift effects.
    """
    cohort = None
    try:
        variant = Variant(contig=8, start=145617535, ref="GGGGGTGCAAGGTGA", alt="", ensembl=75)
        patient = Patient(id="patient", os=3, pfs=2, deceased=False, progressed=False, variants=VariantCollection([variant]))
        cohort_cache_path = generated_data_path("cache")
        cohort = Cohort(
            patients=[patient],
            cache_dir=cohort_cache_path)

        def frameshift_filter(filterable_effect):
            return (type(filterable_effect.effect) == FrameShift)

        effects = cohort.load_effects(filter_fn=frameshift_filter)[patient.id]
        eq_(len(effects), 1)
        all_effects = cohort.load_effects(filter_fn=frameshift_filter,
                                          all_effects=True)[patient.id]
        eq_(len(all_effects), 1)
    finally:
        if cohort is not None:
            cohort.clear_caches()

def test_effects_priority_caching():
    """
    Make sure that effects are cached such that they are not filtered
    prematurely. See https://github.com/hammerlab/cohorts/issues/252.
    """
    cohort = None
    try:
        # This variant has IntronicSpliceSite, Subsitution effects, and more.
        variant = Variant(contig=3, start=20212211, ref="C", alt="T", ensembl=75)
        patient = Patient(id="patient", os=3, pfs=2, deceased=False, progressed=False, variants=VariantCollection([variant]))
        cohort_cache_path = generated_data_path("cache")
        cohort = Cohort(
            patients=[patient],
            cache_dir=cohort_cache_path)

        # All of the effects.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects(all_effects=True)[patient.id]
            eq_(len(effects), 15)

        # Top priority effect.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects()[patient.id]
            eq_(len(effects), 1)
            eq_(type(effects[0]), IntronicSpliceSite)

        def missense_snv_filter(filterable_effect):
            return (type(filterable_effect.effect) == Substitution and
                    filterable_effect.variant.is_snv)

        # All missense SNV effects, from the large cache.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects(all_effects=True, filter_fn=missense_snv_filter)[patient.id]
            eq_(len(effects), 6)

        # Top missense SNV effect, from the large cache.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects(filter_fn=missense_snv_filter)[patient.id]
            eq_(len(effects), 1)
            eq_(type(effects[0]), Substitution)

        # Top missense SNV effects, from the small nonsynonymous cache.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects(only_nonsynonymous=True, filter_fn=missense_snv_filter)[patient.id]
            eq_(len(effects), 1)
            eq_(type(effects[0]), Substitution)

        # All nonsynonymous effects, from the small nonsynonymous cache.
        cohort.clear_caches()
        for i in range(2):
            effects = cohort.load_effects(all_effects=True, only_nonsynonymous=True)[patient.id]
            eq_(len(effects), 6)
    finally:
        if cohort is not None:
            cohort.clear_caches()
