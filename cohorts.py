from __future__ import print_function

from os import path
import pandas as pd

import varcode
from varcode import VariantCollection

class Cohort(object):
    """Represents a cohort of patients."""

    def __init__(self,
                 data_dir,
                 sample_ids,
                 normal_bam_ids,
                 tumor_bam_ids,
                 snv_file_format_funcs=None,
                 snv_merge_type=None,
                 indel_file_format_funcs=None,
                 indel_merge_type=None):
        self.data_dir = data_dir
        self.normal_bam_ids = normal_bam_ids
        self.tumor_bam_ids = tumor_bam_ids
        self.snv_file_format_funcs = snv_file_format_funcs
        self.snv_merge_type = snv_merge_type
        self.indel_file_format_funcs = indel_file_format_funcs
        self.indel_merge_type = indel_merge_type
        self.sample_ids = sample_ids

    def load_mutations(self, variant_type="snv"):
        assert variant_type in ["snv", "indel"]
        sample_variants = {}

        for i, sample_id in enumerate(self.sample_ids):
            try:
                if variant_type == "snv":
                    variants = self._load_single_sample_mutations(
                        sample_id, self.normal_bam_ids[i], self.tumor_bam_ids[i],
                        self.snv_file_format_funcs, self.snv_merge_type)
                elif variant_type == "indel":
                    variants = self._load_single_sample_mutations(
                        sample_id, self.normal_bam_ids[i], self.tumor_bam_ids[i],
                        self.indel_file_format_funcs, self.indel_merge_type) 
            except:
                print("Variants did not exist for %s" % sample_id)
                continue

            sample_variants[sample_id] = variants
        return sample_variants

    def _load_single_sample_mutations(self, sample_id, normal_bam_id, tumor_bam_id,
                                    file_format_funcs, merge_type):
        combined_variants = []
        for file_format_func in file_format_funcs:
            file_name = file_format_func(
                sample_id, normal_bam_id, tumor_bam_id)
            variants = varcode.load_vcf_fast(path.join(self.data_dir, file_name))
            combined_variants.append(set(variants.elements))

        if len(combined_variants) == 1:
            assert merge_type == None
            merged_variants =  VariantCollection(combined_variants)
        else:
            assert merge_type in ["union", "intersection"]
            if merge_type == "union":
                merged_variants = VariantCollection(set.union(*combined_variants))
            elif merge_type == "intersection":
                merged_variants = VariantCollection(set.intersection(*combined_variants))
        return merged_variants
