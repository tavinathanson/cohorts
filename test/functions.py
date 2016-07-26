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

from cohorts.functions import (
    snv_count as cohorts_snv_count,
    missense_snv_count as cohorts_missense_snv_count,
    expressed_missense_snv_count as cohorts_expressed_missense_snv_count,
    neoantigen_count as cohorts_neoantigen_count,
    expressed_neoantigen_count as cohorts_expressed_neoantigen_count
)

# Create functions that don't use QC filtering or exon MB normalization
# because those aren't tested at the moment.
# TODO: Tests should simply work with these options set to True.

def snv_count(row, cohort, **kwargs):
    return cohorts_snv_count(row, cohort, normalized_per_mb=False, **kwargs)

def missense_snv_count(row, cohort, **kwargs):
    return cohorts_missense_snv_count(row, cohort, normalized_per_mb=False, **kwargs)

def neoantigen_count(row, cohort, **kwargs):
    return cohorts_neoantigen_count(row, cohort, normalized_per_mb=False, **kwargs)

def expressed_missense_snv_count(row, cohort, **kwargs):
    return cohorts_expressed_missense_snv_count(row, cohort, normalized_per_mb=False, **kwargs)

def expressed_neoantigen_count(row, cohort, **kwargs):
    return cohorts_expressed_neoantigen_count(row, cohort, normalized_per_mb=False, **kwargs)
