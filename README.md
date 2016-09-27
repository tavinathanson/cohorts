[![PyPI](https://img.shields.io/pypi/v/cohorts.svg?maxAge=2592000)]() [![Build Status](https://travis-ci.org/hammerlab/cohorts.svg?branch=master)](https://travis-ci.org/hammerlab/cohorts) [![Coverage Status](https://coveralls.io/repos/hammerlab/cohorts/badge.svg?branch=master&service=github)](https://coveralls.io/github/hammerlab/cohorts?branch=master)

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
patient_1 = Patient(
    id="patient_1",
    os=70,
    pfs=24,
    deceased=True,
    progressed=True,
    benefit=False
)
    
patient_2 = Patient(
    id="patient_2",
    os=100,
    pfs=50,
    deceased=False,
    progressed=True,
    benefit=False
)

cohort = Cohort(
    patients=[patient_1, patient_2],
    cache_dir="/where/cohorts/results/get/saved"
)

cohort.plot_survival(on="os")
```

```python
sample_1_tumor = Sample(
    is_tumor=True,
    bam_path_dna="/path/to/dna/bam",
    bam_path_rna="/path/to/rna/bam"
)

patient_1 = Patient(
    id="patient_1",
    ...
    snv_vcf_paths=["/where/my/mutect/vcfs/live",
                   "/where/my/strelka/vcfs/live"]
    indel_vcfs_paths=[...],
    tumor_sample=sample_1_tumor,
    ...
)

cohort = Cohort(
    ...
    patients=[patient_1]
)

# Comparison plot of missense mutation counts between benefit and no-benefit patients
cohort.plot_benefit(on=missense_snv_count)

# Raw missense mutations counts
missense_snv_col, dataframe = cohort.as_dataframe(missense_snv_count)
(col_1, col_2), dataframe = cohort.as_dataframe([missense_snv_count, neoantigen_count])
```
