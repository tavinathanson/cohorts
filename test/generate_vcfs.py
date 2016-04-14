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

from . import data_path, generated_data_path

import vcf
from os import path, makedirs
from numpy.random import choice, seed

BASES = ["A", "C", "T", "G"]

def generate_test_vcfs(id_to_mutation_count):
    """
    Generate fake, random VCFs per-sample

    Parameters
    ----------
    id_to_mutation_count : dict
        sample ID to number of mutations we want to generate for that sample

    Returns
    -------
    str
        Path to the generated VCF directory
    """
    seed(1234)
    for sample_id in id_to_mutation_count.keys():
        vcf_reader = vcf.Reader(filename=data_path("vcf_template.vcf"))
        file_path = generated_data_path(
            path.join("vcfs", "sample_%d.vcf" % sample_id))
        file_dir = path.dirname(file_path)
        if not path.exists(file_dir):
            makedirs(file_dir)
        with open(file_path, "w") as f:
            vcf_writer = vcf.Writer(f, vcf_reader)
            record = list(vcf_reader)[0]
            for i in range(id_to_mutation_count[sample_id]):
                record.CHROM = choice(range(1, 20))
                record.POS = choice(range(1000))
                record.REF = choice(BASES)
                bases_no_ref = [base for base in BASES if base != record.REF]
                record.ALT = choice(bases_no_ref)
                vcf_writer.write_record(record)
    return path.dirname(f.name)
