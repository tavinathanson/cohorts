from __future__ import print_function

from os import path, makedirs
import pandas as pd
import pickle

import varcode
from varcode import VariantCollection

class Cohort(object):
    """Represents a cohort of patients."""

    def __init__(self,
                 data_dir,
                 cache_dir,
                 sample_ids,
                 normal_bam_ids,
                 tumor_bam_ids,
                 cache_results=True,
                 snv_file_format_funcs=None,
                 indel_file_format_funcs=None):
        self.data_dir = data_dir
        self.cache_dir = cache_dir
        self.normal_bam_ids = normal_bam_ids
        self.tumor_bam_ids = tumor_bam_ids
        self.cache_results = cache_results
        self.snv_file_format_funcs = snv_file_format_funcs
        self.indel_file_format_funcs = indel_file_format_funcs
        self.sample_ids = sample_ids

    def load_mutations(self, variant_type="snv", merge_type="union"):
        assert variant_type in ["snv", "indel"], "Unknown variant type: %s" % variant_type
        sample_variants = {}

        # Create a file to save cached_results
        cache_dir = path.join(self.cache_dir, "cached-mutations")

        for i, sample_id in enumerate(self.sample_ids):
            sample_cache_dir = path.join(cache_dir, str(sample_id))
            variants_cache_file = path.join(sample_cache_dir,
                                            "%s-%s-variants.pkl" % (variant_type, merge_type))

            # If we are caching results and we have the final final files we load them
            if self.cache_results and path.exists(variants_cache_file):
                sample_variants[sample_id] = pickle.load(open(variants_cache_file, "rb"))
            else:
                if self.cache_results and not path.exists(sample_cache_dir):
                    # Make the directory to save the results
                    makedirs(sample_cache_dir)
                try:
                    if variant_type == "snv":
                        variants = self._load_single_sample_mutations(
                            sample_id, self.normal_bam_ids[i], self.tumor_bam_ids[i],
                            self.snv_file_format_funcs, merge_type)
                    elif variant_type == "indel":
                        variants = self._load_single_sample_mutations(
                            sample_id, self.normal_bam_ids[i], self.tumor_bam_ids[i],
                            self.indel_file_format_funcs, merge_type) 
                except IOError:
                    print("Variants did not exist for %s" % sample_id)
                    continue

                sample_variants[sample_id] = variants

                if self.cache_results:
                    with open(variants_cache_file, "wb") as cache_file:
                        pickle.dump(variants, cache_file)

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
            assert merge_type == None, "Cannot specify a merge type when there is nothing to merge"
            merged_variants =  VariantCollection(combined_variants)
        else:
            assert merge_type in ["union", "intersection"], "Unknown merge type: %s" % merge_type
            if merge_type == "union":
                merged_variants = VariantCollection(set.union(*combined_variants))
            elif merge_type == "intersection":
                merged_variants = VariantCollection(set.intersection(*combined_variants))
        return merged_variants
