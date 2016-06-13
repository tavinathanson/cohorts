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

from cohorts.coxph import *

import pandas as pd
import numpy as np
import re
import statsmodels as sm
from nose.tools import raises, eq_, ok_

def prep_tcga_coxph_data():
    d = pd.DataFrame(np.genfromtxt('nationwidechildrens.org_clinical_patient_blca.txt', delimiter='\t', skip_header = 1, names = True, dtype=None, comments="CDE_ID", autostrip=True))
    d['status'] = (d.vital_status == 'Dead').astype(int)
    def mk_int(s, conv = int):
        non_decimal = re.compile(r'[^\d.]+')
        s = non_decimal.sub('', s)
        return conv(s) if s else np.nan
    d['days_to_death'] = d.apply(lambda row: mk_int(row.days_to_death), axis = 1)
    d['days_to_last_followup'] = d.apply(lambda row: mk_int(row.days_to_last_followup), axis = 1)
    d['weight'] = d.apply(lambda row: mk_int(row.weight, conv = float), axis = 1)
    d['height'] = d.apply(lambda row: mk_int(row.height, conv = float), axis = 1)
    d['time'] = d.apply(lambda row: min(row.days_to_death, row.days_to_last_followup), axis = 1)
    return(d)

def prep_alt_coxph_data():
    d = sm.datasets.get_rdataset('cancer','survival').data
    d['status'] = (d.status == 2).astype(int)
    d['ecog'] = d['ph.ecog']
    d['karno'] = d['ph.karno']
    d['weightloss'] = d['wt.loss']
    return(d)

def test_coxph_model_tcga():
    d = prep_tcga_coxph_data()
    cv, data, concordance  = coxph_model(data = d,
                         formula = '~ gender + height + weight + age_at_initial_pathologic_diagnosis',
                         time_col = 'time',
                         event_col = 'status',
                         )
    return(cv)

def test_coxph_model_alt():
    d = prep_alt_coxph_data()
    cv, data, concordance = coxph_model(data = d,
                                        formula = '~ age + sex + ecog + karno',
                                        time_col = 'time',
                                        event_col = 'status',
                                        )
    return(cv)

## example with missing data for some covariates
def test_coxph_model_alt_missing():
    d = prep_alt_coxph_data()
    cv, data, concordance = coxph_model(data = d,
                                        formula = '~ age + sex + ecog + karno + weightloss',
                                        time_col = 'time',
                                        event_col = 'status',
                                        )
    return(cv)

## test bootstrap code
def test_coxph_bootstrap():
    d = prep_alt_coxph_data()
    models = {'first_model': '~ age + sex', 'second_model': '~ age + sex + ecog'}
    results = bootstrap_coxph(data = d, models = models, time_col = 'time', event_col = 'status')
    return(results) 
