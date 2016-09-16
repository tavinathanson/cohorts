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

class Sample(object):
    """
    Represents a single tumor or normal sample. It can point to DNA and/or
    RNA reads.

    Parameters
    __________
    is_tumor : bool
        Does this `Sample` represent a tumor sample?
    bam_path_dna : str
        Path to the DNA BAM file.
    bam_path_rna : str
        Path to the RNA BAM file.
    cufflinks_path : str
        Path to the Cufflinks output file.
    """
    def __init__(self,
                 is_tumor,
                 bam_path_dna=None,
                 bam_path_rna=None,
                 cufflinks_path=None,
                 kallisto_path=None):
        self.is_tumor = is_tumor
        self.bam_path_dna = bam_path_dna
        self.bam_path_rna = bam_path_rna
        self.cufflinks_path = cufflinks_path
        self.kallisto_path = kallisto_path
