import pandas as pd
import path

def mean_ensembl_coverage(cohort, coverage_path=None,
                          pageant_dir_fn=None):
    """
    calculate mean coverage using Pageant CoverageDepth results with Ensembl loci.

    coverage_path is a path to Pageant CoverageDepth output directory, with
    one subdirectory per patient and a `cdf.csv` file inside each patient subdir.

    pageant_dir_fn is a function that takes in a Patient and produces a Pageant
    dir name.

    Last tested with Pageant CoverageDepth version 1ca9ed2.
    """
    # Function to grab the pageant file name using the Patient
    if pageant_dir_fn is None:
        pageant_dir_fn = cohort.pageant_dir_fn
    if pageant_dir_fn is None:
        pageant_dir_fn = lambda patient: patient.id
    if coverage_path is None:
        coverage_path = cohort.coverage_path

    columns_both = [
        "depth1", # Normal
        "depth2", # Tumor
        "onBP1",
        "onBP2",
        "numOnLoci",
        "fracBPOn1",
        "fracBPOn2",
        "fracLociOn",
        "offBP1",
        "offBP2",
        "numOffLoci",
        "fracBPOff1",
        "fracBPOff2",
        "fracLociOff",
    ]
    columns_single = [
        "depth",
        "onBP",
        "numOnLoci",
        "fracBPOn",
        "fracLociOn",
        "offBP",
        "numOffLoci",
        "fracBPOff",
        "fracLociOff"
    ]
    columns = columns_single if use_tumor_only else columns_both
    ensembl_loci_dfs = []
    for patient in cohort:
        patient_ensembl_loci_df = pd.read_csv(
            path.join(coverage_path, pageant_dir_fn(patient), "cdf.csv"),
            names=columns,
            header=1)
        patient_ensembl_loci_df["patient_id"] = patient.id
        ensembl_loci_dfs.append(patient_ensembl_loci_df)
    ensembl_loci_df = pd.concat(ensembl_loci_dfs)
    ensembl_loci_df["MB"] = ensembl_loci_df.numOnLoci / 1000000.0
    return ensembl_loci_df.groupby("patient_id").apply(wavg, "depth", "numOnLoci").mean()

def wavg(group, avg_name, weight_name):
    """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
    In rare instance, we may not have weights, so just return the mean. Customize this if your business case
    should return otherwise.
    """
    d = group[avg_name]
    w = group[weight_name]
    try:
        return (d * w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()

