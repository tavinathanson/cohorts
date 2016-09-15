
# Using cohorts

[Cohorts](http://github.com/hammerlab/cohorts) is a library for analyzing and plotting clinical data, mutations and neoepitopes in patient cohorts.

It calls out to external libraries like topiary and caches the results for easy manipulation.


# Motivating use case

The motivating use case for Cohorts is a clinical study for checkpoint blockade (ie cancer immunotherapy).

These studies typically collect a wide range of data for each patient, such as:

- clinical data (e.g. outcome, treatment assignment, gender, risk factors, etc)
- whole-exome sequencing (WES) or whole-genome sequencing (WGS) data
   - typically for tumor & normal samples
- expression data (e.g. RNA-seq)
- immune infiltrate data 
- targeted sequencing data (e.g. IMPACT calls or foundation panel data)
- etc.

Data inputs are processed into clinically-relevant features, such as:

- somatic mutations
- predicted neoantigens
- mutational signatures
- tumor heterogeneity estimates
- etc.

Because of the computational overhead frequently involved in computing the above-mentioned 
features, Cohorts caches these computed features by default. The [caching strategy](#caching-strategy)
will be described in some detail below, but in general Cohorts uses file-based caching.

Cohorts then provides utilities to aid in filtering, summarizing & evaluating these features and 
their correlation with a variety of clinical outcomes in either a Survival-analysis context or 
a logistic regression. 

Clinical outcomes evaluated typically include survival, progression-free survival, durable clinical 
benefit, RECIST, or mRECIST criteria.  

We expect that the sets of input data, features & clinical outcome measures will not only 
vary from cohort to cohort, but will also likely grow in variety over time. [Cohorts](http://github.com/hammerlab/cohorts) is designed with this flexibility in mind.


# Design philosophy

A [Cohorts](http://github.com/hammerlab/cohorts) is, at its most basic level, a collection of Patients.

## Minimal example

At a minimum, each Patient has the following data elements: 

## Example with genetic data

# Concepts

## Filters

## Provenance

# Caching strategy

