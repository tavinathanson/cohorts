import pandas as pd
import sys
import os
from cohorts.r_utils import call_R, write_r_template_to_file
from shutil import rmtree

R_TEMPLATE = """
library('deconstructSigs')

result <- list()
counts <- list()

ids <- c({sample_ids})

for(sample in ids) {{
  print(sample)
  data <- read.csv(sprintf("{input_dir}/%s.input", sample), colClasses=c("Sample"="character"));

  data$Sample <- as.character(data$Sample)
  sigs.input <- mut.to.sigs.input(mut.ref = data,
                                  sample.id = "Sample",
                                  chr = "Chr",
                                  pos = "Position",
                                  ref = "Ref",
                                  alt = "Alt")
  output <- whichSignatures(tumor.ref = sigs.input,
                            signatures.ref = signatures.cosmic,
                            sample.id = sample, contexts.needed = TRUE);

  # Add the sample name to a column named sample
  output$weights['Sample'] <- sample

  # Save the weights dataframe with a key/index of the `sample`
  result[[sample]] <- output$weights
  counts[[sample]] <- sigs.input
}}
# Merge all variant counts into a single dataframe
finalCounts <- reshape::merge_all(counts)
write.csv(file = '{output_dir}/signature_variant_counts.csv', finalCounts)
# Merge all the weight vectors into a single dataframe
signatureDeconvolutions <- reshape::merge_all(result)
write.csv(file = '{output_dir}/cohort_signatures.csv', signatureDeconvolutions)
"""

def run_deconstructsigs(input_dir,
                        sample_id_list,
                        output_dir=None,
                        path_to_write_r_script=None,
                        write_variant_counts=True):
    """
    Uses subprocess to call the R script.
    Patient IDs are loaded (how?)
    Output variant counts
    """
    # Create the default dir for signatures_path and variants_counts_path
    if output_dir is None:
        home_dir = os.path.expanduser('~')
        output_dir = "{}/cohort_deconstructsigs_output".format(home_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    signatures_path = "{}/cohort_signatures.csv".format(output_dir)
    variant_counts_path = "{}/cohort_signature_variant_counts.csv".format(output_dir)

    sample_ids = "'" + "','".join(sample_id_list) + "'"
    r_script = R_TEMPLATE.format(sample_ids=sample_ids,
                                 input_dir=input_dir,
                                 output_dir=output_dir)
    #import pdb
    #pdb.set_trace()
    r_path = write_r_template_to_file(r_script, path_to_write_r_script)
    call_R(r_path)
    rmtree(path_to_write_r_script)

def make_deconstructsigs_patient_inputs(destination_dir, variants):
    """
    Parameters
    ___________
    destination_dir : str
        directory to create / put per-patient files
    variants : dict {str : VariantCollection}
        output of cohort.load_variants()
    """
    if not os.path.exists(destination_dir):
        os.mkdir(destination_dir)
        print "created dir"
    for sample in variants:
        df = variant_collection_to_df(variants[sample], sample)
        df.to_csv('{}/{}.input'.format(destination_dir, sample), index=False)
    print "completed"

def variant_collection_to_df(collection, sample_id):
    """
    Creates a DataFrame for the sample_id from the VariantCollection with the below format:
    Alt,Chr,Position,Ref,Sample
    A,chr1,985349,G,CI6743
    T,chr1,4235925,G,CI6743

    Skips variants on chromosomes that have b37-specific names.

    Parameters
    ______
    collection : VariantCollection
    sample_id : str
    """
    rows = []
    for variant in collection:
        # Chromosome nomenclature from b37 (deconstructSigs assumes hg19)
        if variant.contig.startswith('GL') or variant.contig == 'MT' or variant.contig == 'HS37D5':
            continue
        row = {}
        row['Sample'] = sample_id
        row['Chr'] = 'chr' + variant.contig
        row['Position'] = variant.start
        row['Ref'] = variant.ref
        row['Alt'] = variant.alt
        rows.append(row)
    return pd.DataFrame.from_records(rows)

