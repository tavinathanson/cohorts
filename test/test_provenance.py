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

from nose.tools import eq_, ok_
import pandas as pd
from os import path
import warnings

from .test_basic import make_simple_cohort

def test_provenance():
    cohort = None
    try:
        # First, let's just make a cohort and cache a file. Before we do this,
        # we have no provenance.
        cohort = make_simple_cohort()
        cohort.check_provenance = True
        df_empty = pd.DataFrame({"a": [1]})
        patient_id = "1"
        cache_name = cohort.cache_names["variant"]
        cache_dir = path.join(cohort.cache_dir, cache_name)
        patient_cache_dir = path.join(cache_dir, str(patient_id))
        provenance_path = path.join(patient_cache_dir, "PROVENANCE")
        ok_(not path.exists(provenance_path))

        # After we cache the file, we should have provenance.
        cohort.save_to_cache(df_empty, cache_name, patient_id, "cached_file.csv")
        ok_(path.exists(provenance_path))

        # Now let's mess with the cached provenance file.
        df_empty = cohort.load_from_cache(cache_name, patient_id, "cached_file.csv")
        provenance = cohort.load_provenance(patient_cache_dir)
        provenance["hello"] = "1.0.1"
        cohort.save_provenance(patient_cache_dir, provenance)

        # We should get a warning because our saved provenance is different
        # from our current environment.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cohort.load_from_cache(cache_name, patient_id, "cached_file.csv")
            eq_(len(w), 1)

        # But, when check_provenance is off, we shouldn't get a warning.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cohort.check_provenance = False
            cohort.load_from_cache(cache_name, patient_id, "cached_file.csv")
            eq_(len(w), 0)
    finally:
        if cohort is not None:
            cohort.clear_caches()
