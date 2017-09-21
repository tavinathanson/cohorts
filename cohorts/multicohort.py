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

from .cohort import Cohort

class MultiCohort(Cohort):
    def __init__(self,
                 cohorts):
        self.cohorts = cohorts
        patient_id_to_cache_dir = {}
        patients = []
        for cohort in cohorts:
            for patient in cohort:
                patient_id_to_cache_dir[patient.id] = cohort.cache_dir
                patients.append(patient)
        if len(patient_id_to_cache_dir) != len(patients):
            raise ValueError("Some patients are repeated across Cohorts, which MultiCohort does not currently support.")
        self.patient_id_to_cache_dir = patient_id_to_cache_dir

        Cohort.__init__(self,
                        patients=patients,
                        cache_dir="dummy",
                        cache_root_dir=None,
                        cache_dir_kwargs={},
                        show_progress=cohorts[0].show_progress,
                        kallisto_ensembl_version=cohorts[0].kallisto_ensembl_version,
                        cache_results=cohorts[0].cache_results,
                        extra_df_loaders=[],
                        join_with=cohorts[0].join_with,
                        filter_fn=cohorts[0].filter_fn,
                        mhc_class=cohorts[0].mhc_class,
                        normalized_per_mb=cohorts[0].normalized_per_mb,
                        min_coverage_normal_depth=cohorts[0].min_coverage_normal_depth,
                        min_coverage_tumor_depth=cohorts[0].min_coverage_tumor_depth,
                        responder_pfs_equals_os=cohorts[0].responder_pfs_equals_os,
                        check_provenance=cohorts[0].check_provenance,
                        print_provenance=True,
                        print_filter=True,
                        polyphen_dump_path=cohorts[0].polyphen_dump_path,
                        pageant_coverage_path=cohorts[0].pageant_coverage_path,
                        pageant_dir_fn=cohorts[0].pageant_dir_fn,
                        additional_maf_cols=cohorts[0].additional_maf_cols,
                        benefit_plot_name=cohorts[0].benefit_plot_name,
                        merge_type=cohorts[0].merge_type)

    def select_cache_dir(self, patient_id):
        return self.patient_id_to_cache_dir[patient_id]

    def clear_cache(self, cache):
        for cache_dir in set(self.patient_id_to_cache_dir.values()):
            Cohort.clear_cache_helper(self, cache, cache_dir)
