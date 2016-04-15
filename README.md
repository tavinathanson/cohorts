[![Build Status](https://travis-ci.org/tavinathanson/cohorts.svg?branch=master)](https://travis-ci.org/tavinathanson/cohorts) [![Coverage Status](https://coveralls.io/repos/tavinathanson/cohorts/badge.svg?branch=master&service=github)](https://coveralls.io/github/tavinathanson/cohorts?branch=master)

Cohorts
=======

Cohorts is a library for analyzing and plotting clinical data, mutations and neoepitopes in patient cohorts.

It calls out to external libraries like [topiary](https://github.com/hammerlab/topiary) and caches the results for easy manipulation.

Installation
------------

You can install Cohorts using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```bash
pip install cohorts
```

Usage Examples
--------------

```python
cohort = Cohort(
    data_dir="/my/input/data",
    cache_dir="/where/cohorts/results/get/saved",
    sample_ids=["sample_1", "sample_2"],
    clinical_dataframe=pandas_dataframe_with_clinical_data,
    clinical_dataframe_id_col="sample_id_in_dataframe",
    os_col="Overall Survival",
    pfs_col="Progression-Free Survival",
    deceased_col="deceased",
    progressed_or_deceased_col="progressed_or_deceased"
)

cohort.plot_survival(how="os")
```

```python
def mutect_snv_file_format_func(sample_id, normal_bam_id, tumor_bam_id):
    return "Mutect-%d-normal=%s.bam-tumor=%s.bam-merged.vcf" % (
        sample_id, normal_bam_id, tumor_bam_id)

def strelka_snv_file_format_func(...):
    ...

cohort = Cohort(
    ...
    benefit_col="patient_durable_benefit",
    snv_file_format_funcs=[
        mutect_snv_file_format_func,
        strelka_snv_file_format_func
    ]
)

# Comparison plot of missense mutation counts between benefit and no-benefit patients
cohort.plot_benefit(on=missense_snv_count)

# Raw missense mutations counts
missense_snv_col, updated_dataframe = missense_snv_count(cohort)
```
