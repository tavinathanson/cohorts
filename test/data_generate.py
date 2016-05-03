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
from varcode import load_vcf_fast
from os import path, makedirs

def generate_vcfs(id_to_mutation_count, file_format, template_name):
    """
    Generate cropped VCFs from a template, for each sample.

    Parameters
    ----------
    id_to_mutation_count : dict
        sample ID to number of mutations we want to generate for that sample

    Returns
    -------
    str
        Path to the generated VCF directory
    """
    for sample_id in id_to_mutation_count.keys():
        template_path = data_path(template_name)
        vcf_reader = vcf.Reader(filename=template_path)
        file_path = generated_data_path(
            path.join("vcfs", file_format % sample_id))
        file_dir = path.dirname(file_path)
        if not path.exists(file_dir):
            makedirs(file_dir)
        with open(file_path, "w") as f:
            vcf_writer = vcf.Writer(f, vcf_reader)
            i = 0
            num_records_in_template = len(load_vcf_fast(template_path))
            num_records_to_generate = id_to_mutation_count[sample_id]
            assert num_records_in_template >= num_records_to_generate, (
                "Cannot generate more records than exist in the template: %d is less than %d" % (
                    num_records_in_template, num_records_to_generate))
            for record in vcf_reader:
                if i < id_to_mutation_count[sample_id]:
                    vcf_writer.write_record(record)
                    i += 1
                else:
                    break

    return path.dirname(f.name)
