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

import pandas as pd
from numpy.random import choice, randint, seed

from . import Cohort, Patient

def random_cohort(size, cache_dir, seed_val=1234):
    seed(seed_val)
    d = {}
    d["id"] = [str(id) for id in range(size)]
    d["age"] = choice([10, 15, 28, 32, 59, 62, 64, 66, 68], size)
    d["OS"] = [os + randint(10) for os in choice([10, 100, 500, 1000], size)]
    # Note: these values are not currently consistent with each other.
    d["PFS"] = [int(os * 0.6) for os in d["OS"]]
    d["benefit"] = choice([False, True], size)
    d["random"] = [randint(100) for i in range(size)]
    d["random_boolean"] = choice([False, True], size)
    d["benefit_correlate"] = [randint(50) if benefit else randint(20) for benefit in d["benefit"]]
    d["benefit_correlate_boolean"] = [True if corr > 10 else False for corr in d["benefit_correlate"]]
    d["deceased"] = choice([False, True], size)
    d["progressed_or_deceased"] = choice([False, True], size)
    df = pd.DataFrame(d)
    patients = []
    for i, row in df.iterrows():
         patient = Patient(
             id=row["id"],
             os=row["OS"],
             pfs=row["PFS"],
             benefit=row["benefit"],
             deceased=row["deceased"],
             progressed_or_deceased=row["progressed_or_deceased"],
             additional_data=row)
         patients.append(patient)
    return Cohort(
        patients=patients,
        cache_dir=cache_dir)
