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

Usage
-----

```python
cohort = Cohort(
    data_dir="/my/clinical/data",
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
