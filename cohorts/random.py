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
from os import path
from numpy.random import choice, randint, seed
from varcode import Variant, VariantCollection
from varcode.effects import Substitution
from mhctools.random_predictor import RandomBindingPredictor

from . import Cohort, Patient

def random_cohort(size, cache_dir, data_dir=None,
                  min_random_variants=None,
                  max_random_variants=None,
                  seed_val=1234):
    """
    Parameters
    ----------
    min_random_variants: optional, int
        Minimum number of random variants to be generated per patient.
    max_random_variants: optional, int
        Maximum number of random variants to be generated per patient.
    """
    seed(seed_val)
    d = {}
    d["id"] = [str(id) for id in range(size)]
    d["age"] = choice([10, 15, 28, 32, 59, 62, 64, 66, 68], size)
    d["smoker"] = choice([False, True], size)
    d["OS"] = [randint(10, 1000) for i in range(size)]
    d["PFS"] = [int(os * 0.6) for os in d["OS"]]
    d["benefit"] = [pfs < 50 for pfs in d["PFS"]]
    d["random"] = [randint(100) for i in range(size)]
    d["random_boolean"] = choice([False, True], size)
    d["benefit_correlate"] = [randint(50) if benefit else randint(20) for benefit in d["benefit"]]
    d["benefit_correlate_boolean"] = [True if corr > 10 else False for corr in d["benefit_correlate"]]
    d["deceased"] = choice([False, True], size)
    d["progressed_or_deceased"] = [deceased or choice([False, True]) for deceased in d["deceased"]]
    df = pd.DataFrame(d)
    patients = []
    for i, row in df.iterrows():
        snv_vcf_paths = None
        if max_random_variants is not None and min_random_variants is not None:
            if data_dir is None:
                raise ValueError("data_dir must be provided if random variants are being generated.")
            vcf_path = path.join(data_dir, "patient_%s_mutect.vcf" % row["id"])
            generate_simple_vcf(
                vcf_path, generate_random_missense_variants(num_variants=randint(min_random_variants, max_random_variants)))
            snv_vcf_paths = [vcf_path]
        patient = Patient(
            id=row["id"],
            os=row["OS"],
            pfs=row["PFS"],
            benefit=row["benefit"],
            deceased=row["deceased"],
            progressed_or_deceased=row["progressed_or_deceased"],
            hla_alleles=["HLA-A02:01"],
            variants={"snv": snv_vcf_paths},
            additional_data=row)
        patients.append(patient)
    return Cohort(
        patients=patients,
        cache_dir=cache_dir,
        mhc_class=RandomBindingPredictor)

def generate_random_missense_variants(num_variants=10, max_search=100000, reference="GRCh37"):
    """
    Generate a random collection of missense variants by trying random variants repeatedly.
    """
    variants = []
    for i in range(max_search):
        bases = ["A", "C", "T", "G"]
        random_ref = choice(bases)
        bases.remove(random_ref)
        random_alt = choice(bases)
        random_contig = choice(["1", "2", "3", "4", "5"])
        random_variant = Variant(contig=random_contig, start=randint(1, 1000000),
                                 ref=random_ref, alt=random_alt, ensembl=reference)
        try:
            effects = random_variant.effects()
            for effect in effects:
                if isinstance(effect, Substitution):
                    variants.append(random_variant)
                    break
        except:
            continue
        if len(variants) == num_variants:
            break
    return VariantCollection(variants)

def generate_simple_vcf(filename, variant_collection):
    """
    Output a very simple metadata-free VCF for each variant in a variant_collection.
    """
    contigs = []
    positions = []
    refs = []
    alts = []
    for variant in variant_collection:
        contigs.append("chr" + variant.contig)
        positions.append(variant.start)
        refs.append(variant.ref)
        alts.append(variant.alt)
    df = pd.DataFrame()
    df["contig"] = contigs
    df["position"] = positions
    df["id"] = ["."] * len(variant_collection)
    df["ref"] = refs
    df["alt"] = alts
    df["qual"] = ["."] * len(variant_collection)
    df["filter"] = ["."] * len(variant_collection)
    df["info"] = ["."] * len(variant_collection)
    df["format"] = ["GT:AD:DP"] * len(variant_collection)
    normal_ref_depths = [randint(1, 10) for v in variant_collection]
    normal_alt_depths = [randint(1, 10) for v in variant_collection]
    df["n1"] = ["0:%d,%d:%d" % (normal_ref_depths[i], normal_alt_depths[i],
                                normal_ref_depths[i] + normal_alt_depths[i])
                for i in range(len(variant_collection))]
    tumor_ref_depths = [randint(1, 10) for v in variant_collection]
    tumor_alt_depths = [randint(1, 10) for v in variant_collection]
    df["t1"] = ["0/1:%d,%d:%d" % (tumor_ref_depths[i], tumor_alt_depths[i], tumor_ref_depths[i] + tumor_alt_depths[i])
                for i in range(len(variant_collection))]

    with open(filename, "w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write("##reference=file:///projects/ngs/resources/gatk/2.3/ucsc.hg19.fasta\n")

    with open(filename, "a") as f:
        df.to_csv(f, sep="\t", index=None, header=None)
