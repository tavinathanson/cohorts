[![PyPI](https://img.shields.io/pypi/v/cohorts.svg?maxAge=21600)]() [![Build Status](https://travis-ci.org/hammerlab/cohorts.svg?branch=master)](https://travis-ci.org/hammerlab/cohorts) [![Coverage Status](https://coveralls.io/repos/hammerlab/cohorts/badge.svg?branch=master&service=github)](https://coveralls.io/github/hammerlab/cohorts?branch=master)

Cohorts
=======

Cohorts is a library for analyzing and plotting clinical data, mutations and neoepitopes in patient cohorts.

It calls out to external libraries like [topiary](https://github.com/hammerlab/topiary) and caches the results for easy manipulation.

Cohorts requires Python 3 (3.3+). We are no longer maintaining compatability with Python 2. For context, see this [Python 3 statement](www.python3statement.org).

Installation
------------

You can install Cohorts using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```bash
pip install cohorts
```

Features
--------

* Data management: construct a `Cohort` consisting of `Patient`s with `Sample`s.
* Use `varcode` and `topiary` to generate and cache variant effects and predicted neoantigens.
* Provenance: track the state of the world (package and data versions) for a given analysis.
* Aggregation functions: built-in functions such as `missense_snv_count`, `neoantigen_count`, `expressed_neoantigen_count`; or create your own functions.
* Plotting: survival curves via `lifelines`, response/no response plots (with Mann-Whitney and Fisher's Exact results), ROC curves. Example: `cohort.plot_survival(on=missense_snv_count, how="pfs")`.
* Filtering: filter collections of variants/effects/neoantigens by, for example, variant statistics.
* Pre-define data sets to work with. Example: `cohort.as_dataframe(join_with=["tcr", "pdl1"])`.

In addition, several other libraries make use of `cohorts`:
* [pygdc](http://github.com/hammerlab/pygdc)
* [query_tcga](http://github.com/jburos/query_tcga)

Quick Start
---------------

One way to get started using Cohorts is to use it to analyze TCGA data.

As an example, we can create a cohort using [query_tcga](http://github.com/jburos/query_tcga):

```python
from query_tcga import cohort, config

# provide authentication token
config.load_config('config.ini')

# load patient data
blca_patients = cohort.prep_patients(project_name='TCGA-BLCA',
                                     project_data_dir='data')

# create cohort
blca_cohort = cohort.prep_cohort(patients=blca_patients,
                                 cache_dir='data-cache')
```

Then, use `plot_survival()` to summarize a potential biomarker (e.g. `snv_count`) by survival:.

```python
from cohorts.functions import snv_count
blca_cohort.plot_survival(snv_count, how='os', threshold='median')
```

Which should produce a summary of results including this plot:

![Survival plot example](/docs/survival_plot_example.png)

We could alternatively use `plot_benefit()` to summarize OS>12mo instead of survival:

```python
blca_cohort.plot_benefit(snv_count)
```

![Benefit plot example](/docs/benefit_plot_example.png)


See the full example in the [quick-start notebook](http://nbviewer.jupyter.org/github/hammerlab/tcga-blca/blob/master/Quick-start%20-%20using%20Cohorts%20with%20TCGA%20data.ipynb)

Building from Scratch
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

```
