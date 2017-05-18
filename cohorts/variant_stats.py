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

from collections import namedtuple

VariantStats = namedtuple("VariantStats",
                          ["depth", "alt_depth", "variant_allele_frequency"])
SomaticVariantStats = namedtuple("SomaticVariantStats",
                          ["tumor_stats", "normal_stats"])

def strelka_somatic_variant_stats(variant, variant_metadata):
    """Parse out the variant calling statistics for a given variant from a Strelka VCF

    Parameters
    ----------
    variant : varcode.Variant
    sample_info : dict
        Dictionary of sample to variant calling statistics, corresponds to the sample columns
        in a Strelka VCF

    Returns
    -------
    SomaticVariantStats
    """

    sample_info = variant_metadata["sample_info"]
    # Ensure there are exactly two samples in the VCF, a tumor and normal
    assert len(sample_info) == 2, "More than two samples found in the somatic VCF"
    tumor_stats = _strelka_variant_stats(variant, sample_info["TUMOR"])
    normal_stats = _strelka_variant_stats(variant, sample_info["NORMAL"])
    return SomaticVariantStats(tumor_stats=tumor_stats, normal_stats=normal_stats)

def _strelka_variant_stats(variant, sample_info):
    """Parse a single sample"s variant calling statistics based on Strelka VCF output

    Parameters
    ----------
    variant : varcode.Variant
    sample_info : dict
        Dictionary of Strelka-specific variant calling fields

    Returns
    -------
    VariantStats
    """
    
    if variant.is_deletion or variant.is_insertion:
        # ref: https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output
        ref_depth = int(sample_info['TAR'][0]) # number of reads supporting ref allele (non-deletion)
        alt_depth = int(sample_info['TIR'][0]) # number of reads supporting alt allele (deletion)
        depth = ref_depth + alt_depth
    else:
        # Retrieve the Tier 1 counts from Strelka
        ref_depth = int(sample_info[variant.ref+"U"][0])
        alt_depth = int(sample_info[variant.alt+"U"][0])
        depth = alt_depth + ref_depth
    if depth > 0:
        vaf = float(alt_depth) / depth
    else:
        # unclear how to define vaf if no reads support variant
        # up to user to interpret this (hopefully filtered out in QC settings)
        vaf = None

    return VariantStats(depth=depth, alt_depth=alt_depth, variant_allele_frequency=vaf)

def mutect_somatic_variant_stats(variant, variant_metadata):
    """Parse out the variant calling statistics for a given variant from a Mutect VCF

    Parameters
    ----------
    variant : varcode.Variant
    sample_info : dict
        Dictionary of sample to variant calling statistics, corresponds to the sample columns
        in a Mutect VCF

    Returns
    -------
    SomaticVariantStats
    """

    sample_info = variant_metadata["sample_info"]
    # Ensure there are exactly two samples in the VCF, a tumor and normal
    assert len(sample_info) == 2, "More than two samples found in the somatic VCF"

    # Find the sample with the genotype field set to variant in the VCF
    tumor_sample_infos = [info for info in sample_info.values() if info["GT"] == "0/1"]

    # Ensure there is only one such sample
    assert len(tumor_sample_infos) == 1, "More than one tumor sample found in the VCF file"

    tumor_sample_info = tumor_sample_infos[0]
    normal_sample_info = [info for info in sample_info.values() if info["GT"] != "0/1"][0]

    tumor_stats = _mutect_variant_stats(variant, tumor_sample_info)
    normal_stats = _mutect_variant_stats(variant, normal_sample_info)
    return SomaticVariantStats(tumor_stats=tumor_stats, normal_stats=normal_stats)

def _mutect_variant_stats(variant, sample_info):
    """Parse a single sample"s variant calling statistics based on Mutect"s (v1) VCF output

    Parameters
    ----------
    variant : varcode.Variant
    sample_info : dict
        Dictionary of Mutect-specific variant calling fields

    Returns
    -------
    VariantStats
    """

    # Parse out the AD (or allele depth field), which is an array of [REF_DEPTH, ALT_DEPTH]
    ref_depth, alt_depth = sample_info["AD"]
    depth = int(ref_depth) + int(alt_depth)
    vaf = float(alt_depth) / depth

    return VariantStats(depth=depth, alt_depth=alt_depth, variant_allele_frequency=vaf)

def _maf_variant_stats(variant, variant_metadata, prefix="t"):
    ref_depth = variant_metadata["%s_ref_count" % prefix]
    alt_depth = variant_metadata["%s_alt_count" % prefix]
    depth = int(ref_depth) + int(alt_depth)
    vaf = float(alt_depth) / depth
    return VariantStats(depth=depth, alt_depth=alt_depth, variant_allele_frequency=vaf)

def maf_somatic_variant_stats(variant, variant_metadata):
    """
    Parse out the variant calling statistics for a given variant from a MAF file

    Assumes the MAF format described here: https://www.biostars.org/p/161298/#161777

    Parameters
    ----------
    variant : varcode.Variant
    variant_metadata : dict
        Dictionary of metadata for this variant

    Returns
    -------
    SomaticVariantStats
    """
    tumor_stats = None
    normal_stats = None
    if "t_ref_count" in variant_metadata:
        tumor_stats = _maf_variant_stats(variant, variant_metadata, prefix="t")
    if "n_ref_count" in variant_metadata:
        normal_stats = _maf_variant_stats(variant, variant_metadata, prefix="n")
    return SomaticVariantStats(tumor_stats=tumor_stats, normal_stats=normal_stats)

def variant_stats_from_variant(variant,
                               metadata,
                               merge_fn=(lambda all_stats: \
                                max(all_stats, key=(lambda stats: stats.tumor_stats.depth)))):
    """Parse the variant calling stats from a variant called from multiple variant files. The stats are merged
    based on `merge_fn`

    Parameters
    ----------
    variant : varcode.Variant
    metadata : dict
        Dictionary of variant file to variant calling metadata from that file
    merge_fn : function
        Function from list of SomaticVariantStats to single SomaticVariantStats.
        This is used if a variant is called by multiple callers or appears in multiple VCFs.
        By default, this uses the data from the caller that had a higher tumor depth.

    Returns
    -------
    SomaticVariantStats
    """
    all_stats = []
    for (variant_file, variant_metadata) in metadata.items():
        if ".maf" in variant_file.lower():
            stats = maf_somatic_variant_stats(variant, variant_metadata)
        elif "strelka" in variant_file.lower():
            stats = strelka_somatic_variant_stats(variant, variant_metadata)
        elif "mutect" in variant_file.lower():
            stats = mutect_somatic_variant_stats(variant, variant_metadata)
        else:
            raise ValueError("Cannot parse sample fields, variant file {} is from an unsupported caller.".format(variant_file))
        all_stats.append(stats)
    return merge_fn(all_stats)
