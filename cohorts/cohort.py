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

from os import path, makedirs
from shutil import rmtree
import pandas as pd
import seaborn as sb
import json
import warnings
import pprint
from copy import copy
import dill
import hashlib
import inspect
import logging
import pickle
import numpy as np

# pylint doesn't like this line
# pylint: disable=no-name-in-module
from types import FunctionType

import vap  ## vcf-annotate-polyphen
from sqlalchemy import create_engine

from pyensembl import cached_release

import varcode
from varcode import EffectCollection, VariantCollection
from mhctools import NetMHCcons, EpitopeCollection
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe
from topiary.sequence_helpers import contains_mutant_residues
from isovar.allele_reads import reads_overlapping_variants
from isovar.protein_sequences import reads_generator_to_protein_sequences_generator, protein_sequences_generator_to_dataframe
from pysam import AlignmentFile
from scipy.stats import pearsonr
from collections import defaultdict
from tqdm import tqdm

from .dataframe_loader import DataFrameLoader
from .utils import DataFrameHolder, first_not_none_param, filter_not_null, InvalidDataError, strip_column_names as _strip_column_names, get_logger, get_cache_dir
from .utils import get_function_name, make_hash, hash_function
from .provenance import compare_provenance
from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot, roc_curve_plot, stripboxplot, CorrelationResults
from .model import cohort_coxph, cohort_bootstrap_auc, cohort_mean_bootstrap_auc
from .collection import Collection
from .varcode_utils import (filter_variants, filter_effects,
                            filter_neoantigens, filter_polyphen)
from .variant_filters import no_filter
from .styling import set_styling
from . import variant_filters

logger = get_logger(__name__)
cache_logger = get_logger(__name__ + '.caching')

class Cohort(Collection):
    """
    Represents a cohort of `Patient`s.

    Parameters
    __________
    patients : List
        A list of `Patient`s for this cohort.
    cache_dir : str
        Path to store cached results, e.g. cached variant effects.
    cache_root_dir : str
        (optional) directory in which cache_dir should be created
    cache_dir_kwargs : dict
        (optional) dictionary of name=value data to use when formatting cache_dir str
    show_progress : bool
        Whether or not to show DataFrame application progress as an increasing percentage.
    kallisto_ensembl_version : int
        Cached release version to use from pyensembl to annotate Kallisto data
    cache_results : bool
        Whether or not to cache results.
    extra_df_loaders : List
        List of `DataFrameLoader`s to include as join options `Cohort`.
    join_with : str or List
        The name of one or more `DataFrameLoader`s to join with by default.
    join_how : str
        What type of default join to use for joining `DataFrameLoader`s.
    filter_fn : Function
        Specify a default filter function for `load_variants`, `load_effects`,
        `functions.missense_snv_count`, etc.
    mhc_class : mhctools.BaseCommandlinePredictor, defaults to NetMHCcons
        What MHC binding predictor to use for neoantigen calling.
    normalized_per_mb : bool
        Whether or not to normalize by number of loci.
    min_coverage_depth_normal : int
        When counting number of exonic loci, only count loci with at least this much normal depth.
    min_coverage_depth_tumor : int
        When counting number of exonic loci, only count loci with at least this much tumor depth.
    responder_pfs_equals_os : bool
        Ensure that the PFS values for responders (not progressed) are equal to
        OS values.
    check_provenance : bool
        Verify that the cached provenance is equal to the current environment.
    print_provenance : bool
        Print a summary of cache file provenance.
    print_filter : bool
        Print the name of the default `filter_fn` on initialization.
    polyphen_dump_path : str
        Path to a Polyphen database dump.
    pageant_coverage_path : str
        Path to Pageant CoverageDepth output.
    pageant_dir_fn : function
        Function from patient to a specific Pageant CoverageDepth directory within the path. Defaults to the Patietn ID.
    additional_maf_cols : list
       If loading variants from MAFs, specify any additional columns to pull in from the MAFs.
    benefit_plot_name : str
        What word to use for "benefit" when plotting.
    merge_type : {"union", "intersection"}, optional
        Use this method to merge multiple variant sets for a single patient, default "union"
    """
    def __init__(self,
                 patients,
                 cache_dir,
                 cache_root_dir=None,
                 cache_dir_kwargs=dict(),
                 show_progress=True,
                 kallisto_ensembl_version=None,
                 cache_results=True,
                 extra_df_loaders=[],
                 join_with=None,
                 join_how="inner",
                 filter_fn=None,
                 mhc_class=NetMHCcons,
                 normalized_per_mb=False,
                 min_coverage_normal_depth=0,
                 min_coverage_tumor_depth=0,
                 responder_pfs_equals_os=False,
                 check_provenance=False,
                 print_provenance=True,
                 print_filter=True,
                 polyphen_dump_path=None,
                 pageant_coverage_path=None,
                 pageant_dir_fn=None,
                 additional_maf_cols=None,
                 benefit_plot_name="Benefit",
                 merge_type="union"):
        Collection.__init__(
            self,
            elements=patients)

        # TODO: Patients shouldn't actually need to reference their Cohort; remove
        # this when patient-specific functions all live in Patient.
        for patient in patients:
            patient.cohort = self
        self.cache_dir = get_cache_dir(cache_dir=cache_dir, cache_root_dir=cache_root_dir, **cache_dir_kwargs)
        self.cache_root_dir = cache_root_dir
        self.show_progress = show_progress
        self.cache_results = cache_results
        self.kallisto_ensembl_version = kallisto_ensembl_version

        df_loaders = [
            DataFrameLoader("kallisto", self.load_kallisto),
            DataFrameLoader("cufflinks", self.load_cufflinks),
            DataFrameLoader("ensembl_coverage", self.load_ensembl_coverage)]
        df_loaders.extend(extra_df_loaders)
        self.df_loaders = df_loaders
        self.join_with = join_with
        self.join_how = join_how
        self.filter_fn = filter_fn
        self.mhc_class = mhc_class
        self.normalized_per_mb = normalized_per_mb
        self.min_coverage_normal_depth = min_coverage_normal_depth
        self.min_coverage_tumor_depth = min_coverage_tumor_depth
        self.responder_pfs_equals_os = responder_pfs_equals_os
        self.check_provenance = check_provenance
        self.polyphen_dump_path = polyphen_dump_path
        self.pageant_coverage_path = pageant_coverage_path
        self.pageant_dir_fn = pageant_dir_fn
        self.additional_maf_cols = additional_maf_cols
        self.benefit_plot_name = benefit_plot_name
        self.merge_type = merge_type
        self._genome = None

        self.verify_id_uniqueness()
        self.verify_survival()
        self.dataframe_hash = None

        self.cache_names = {"variant": "cached-variants",
                            "effect": "cached-effects",
                            "nonsynonymous_effect": "cached-nonsynonymous-effects",
                            "neoantigen": "cached-neoantigens",
                            "expressed_neoantigen": "cached-expressed-neoantigens",
                            "polyphen": "cached-polyphen-annotations",
                            "isovar": "cached-isovar-output"}

        if print_filter:
            print("Applying %s filter by default" % self.filter_fn.__name__ if
                  self.filter_fn is not None else "None")
        if print_provenance:
            pprint.pprint(self.summarize_data_sources())

        set_styling()

    def filter(self, filter_fn):
        new_cohort = copy(self)
        new_cohort.elements = [patient for patient in self if filter_fn(patient)]
        return new_cohort

    @property
    def genome(self):
        if self._genome is None:
            self._genome = self.get_genome()
        return self._genome

    def get_genome(self):
        variants = self.load_variants()
        for patient_id, variant_collection in variants.items():
            if len(variant_collection) > 0:
                return variant_collection[0].ensembl
        raise ValueError("No variants to derive genome from")

    def verify_id_uniqueness(self):
        patient_ids = set([patient.id for patient in self])
        if len(patient_ids) != len(self):
            raise ValueError("Non-unique patient IDs")

    def verify_survival(self):
        cohort_dataframe = self.as_dataframe()
        if not (cohort_dataframe["pfs"] <=
                cohort_dataframe["os"]).all():
            raise InvalidDataError("PFS should be <= OS, but PFS is larger than OS for some patients.")

        def func(row):
            if row["pfs"] < row["os"]:
                if not row["progressed_or_deceased"]:
                    raise InvalidDataError(
                        "A patient did not progress despite PFS being less than OS. "
                        "Full row: %s" % row)
        if self.responder_pfs_equals_os:
            cohort_dataframe.apply(func, axis=1)

    def _as_dataframe_unmodified(self, join_with=None, join_how=None):
        # Use join_with if specified, otherwise fall back to what is defined in the class
        join_with = first_not_none_param([join_with, self.join_with], default=[])
        if type(join_with) == str:
            join_with = [join_with]

        # Convert strings to DataFrameLoader objects
        df_loaders = [df_loader for df_loader in self.df_loaders if df_loader.name in join_with]

        # Use join_how if specified, otherwise fall back to what is defined in the class
        join_how = first_not_none_param([join_how, self.join_how], default="inner")

        patient_rows = []
        for patient in self:
            row = {} if patient.additional_data is None else patient.additional_data.copy()
            row["patient_id"] = patient.id
            for clinical_col in ["benefit", "os", "pfs", "deceased",
                                "progressed", "progressed_or_deceased"]:
                row[clinical_col] = getattr(patient, clinical_col)

            patient_rows.append(row)
        df = pd.DataFrame.from_records(patient_rows)

        # Are any columns duplicated in the DataFrame(s) to be joined?
        # If so, rename those columns to be suffixed by the DataFrameLoader
        # name.
        df_loader_dfs = {}
        col_counts = defaultdict(int)
        for df_loader in df_loaders:
            df_loader_dfs[df_loader] = df_loader.load_dataframe()
            for col in df_loader_dfs[df_loader].columns:
                col_counts[col] += 1
        for col, count in col_counts.items():
            # Don't rename columns that are not duplicated.
            if count > 1:
                for df_loader, loaded_df in df_loader_dfs.items():
                    # Don't rename a column that will be joined on.
                    if col != "patient_id" and col != df_loader.join_on:
                        loaded_df.rename(columns={col: "%s_%s" % (col, df_loader.name)}, inplace=True)

        for df_loader, loaded_df in df_loader_dfs.items():
            old_len_df = len(df)
            df = df.merge(
                loaded_df,
                left_on=df_loader.join_on_left,
                right_on=df_loader.join_on_right,
                how=join_how)
            print("%s join with %s: %d to %d rows" % (
                join_how,
                df_loader.name,
                old_len_df,
                len(df)))

        self.dataframe_hash = make_hash(df.sort_values("patient_id"))
        return df

    def as_dataframe(self, on=None, join_with=None, join_how=None,
                     return_cols=False, rename_cols=False,
                     keep_paren_contents=True, **kwargs):
        """
        Return this Cohort as a DataFrame, and optionally include additional columns
        using `on`.

        on : str or function or list or dict, optional
            - A column name.
            - Or a function that creates a new column for comparison, e.g. count.snv_count.
            - Or a list of column-generating functions or column names.
            - Or a map of new column names to their column-generating functions or column names.

        If `on` is a function or functions, kwargs is passed to those functions.
        Otherwise kwargs is ignored.

        Other parameters
        ----------------
        `return_cols`: (bool)
            If True, return column names generated via `on` along with the `DataFrame`
            as a `DataFrameHolder` tuple.
        `rename_cols`: (bool)
            If True, then return columns using "stripped" column names
            ("stripped" means lower-case names without punctuation other than `_`)
            See `utils.strip_column_names` for more details
            defaults to False
        `keep_paren_contents`: (bool)
            If True, then contents of column names within parens are kept.
            If False, contents of column names within-parens are dropped.
            Defaults to True
        ----------

        Return : `DataFrame` (or `DataFrameHolder` if `return_cols` is True)
        """
        df = self._as_dataframe_unmodified(join_with=join_with, join_how=join_how)
        if on is None:
            return DataFrameHolder.return_obj(None, df, return_cols)

        if type(on) == str:
            return DataFrameHolder.return_obj(on, df, return_cols)

        def apply_func(on, col, df):
            """
            Sometimes we have functions that, by necessity, have more parameters
            than just `row`. We construct a function with just the `row` parameter
            so it can be sent to `DataFrame.apply`. We hackishly pass `cohort`
            (as `self`) along if the function accepts a `cohort` argument.
            """
            on_argnames = on.__code__.co_varnames
            if "cohort" not in on_argnames:
                func = lambda row: on(row=row, **kwargs)
            else:
                func = lambda row: on(row=row, cohort=self, **kwargs)

            if self.show_progress:
                tqdm.pandas(desc=col)
                df[col] = df.progress_apply(func, axis=1) ## depends on tqdm on prev line
            else:
                df[col] = df.apply(func, axis=1)
            return DataFrameHolder(col, df)

        def func_name(func, num=0):
            return func.__name__ if not is_lambda(func) else "column_%d" % num

        def is_lambda(func):
            return func.__name__ == (lambda: None).__name__

        if type(on) == FunctionType:
            return apply_func(on, func_name(on), df).return_self(return_cols)

        if len(kwargs) > 0:
            logger.warning("Note: kwargs used with multiple functions; passing them to all functions")

        if type(on) == dict:
            cols = []
            for key, value in on.items():
                if type(value) == str:
                    df[key] = df[value]
                    col = key
                elif type(value) == FunctionType:
                    col, df = apply_func(on=value, col=key, df=df)
                else:
                    raise ValueError("A value of `on`, %s, is not a str or function" % str(value))
                cols.append(col)
        if type(on) == list:
            cols = []
            for i, elem in enumerate(on):
                if type(elem) == str:
                    col = elem
                elif type(elem) == FunctionType:
                    col = func_name(elem, i)
                    col, df = apply_func(on=elem, col=col, df=df)
                cols.append(col)

        if rename_cols:
            rename_dict = _strip_column_names(df.columns, keep_paren_contents=keep_paren_contents)
            df.rename(columns=rename_dict, inplace=True)
            cols = [rename_dict[col] for col in cols]
        return DataFrameHolder(cols, df).return_self(return_cols)

    def load_dataframe(self, df_loader_name):
        """
        Instead of joining a DataFrameJoiner with the Cohort in `as_dataframe`, sometimes
        we may want to just directly load a particular DataFrame.
        """
        logger.debug("loading dataframe: {}".format(df_loader_name))
        # Get the DataFrameLoader object corresponding to this name.
        df_loaders = [df_loader for df_loader in self.df_loaders if df_loader.name == df_loader_name]

        if len(df_loaders) == 0:
            raise ValueError("No DataFrameLoader with name %s" % df_loader_name)
        if len(df_loaders) > 1:
            raise ValueError("Multiple DataFrameLoaders with name %s" % df_loader_name)

        return df_loaders[0].load_dataframe()

    def generate_provenance(self):
        module_names = ["cohorts", "pyensembl", "varcode", "mhctools", "topiary", "isovar", "scipy", "numpy", "pandas"]
        module_versions = [__import__(module_name).__version__ for module_name in module_names]
        return dict(zip(module_names, module_versions))

    def load_provenance(self, patient_cache_dir):
        with open(path.join(patient_cache_dir, "PROVENANCE"), "r") as f:
            return json.load(f)

    def save_provenance(self, patient_cache_dir, provenance):
        with open(path.join(patient_cache_dir, "PROVENANCE"), "w") as f:
            json.dump(provenance, f)

    def load_from_cache(self, cache_name, patient_id, file_name):
        if not self.cache_results:
            return None

        cache_logger.debug("loading patient {} data from {} cache: {}".format(patient_id, cache_name, file_name))

        cache_dir = path.join(self.cache_dir, cache_name)
        patient_cache_dir = path.join(cache_dir, str(patient_id))
        cache_file = path.join(patient_cache_dir, file_name)

        if not path.exists(cache_file):
            cache_logger.debug("cache file does not exist. Checking for older format.")
            # We removed variant_type from the cache name. Eventually remove this notification.
            if (path.exists(path.join(patient_cache_dir, "snv-" + file_name)) or
                path.exists(path.join(patient_cache_dir, "indel-" + file_name))):
                raise ValueError("Cache is in an older format (with variant_type). Please re-generate it.")
            return None

        if self.check_provenance:
            cache_logger.debug("Checking cache provenance")
            num_discrepant = compare_provenance(
                this_provenance = self.generate_provenance(),
                other_provenance = self.load_provenance(patient_cache_dir),
                left_outer_diff = "In current environment but not cached in %s for patient %s" % (cache_name, patient_id),
                right_outer_diff = "In cached %s for patient %s but not current" % (cache_name, patient_id)
                )
        try:
            if path.splitext(cache_file)[1] == ".csv":
                cache_logger.debug("Loading cache as csv file")
                return pd.read_csv(cache_file, dtype={"patient_id": object})
            else:
                cache_logger.debug("Loading cache as pickled file")
                with open(cache_file, "rb") as f:
                    return pickle.load(f)
        except IOError:
            return None

    def save_to_cache(self, obj, cache_name, patient_id, file_name):
        if not self.cache_results:
            return

        cache_logger.debug("saving patient {} data to {} cache: {}".format(patient_id, cache_name, file_name))

        cache_dir = path.join(self.cache_dir, cache_name)
        patient_cache_dir = path.join(cache_dir, str(patient_id))
        cache_file = path.join(patient_cache_dir, file_name)

        if not path.exists(patient_cache_dir):
            makedirs(patient_cache_dir)

        if type(obj) == pd.DataFrame:
            obj.to_csv(cache_file, index=False)
        else:
            with open(cache_file, "wb") as f:
                # Protocol=2 for compatability with Py 2 and 3
                pickle.dump(obj, f)

        provenance = self.generate_provenance()
        self.save_provenance(patient_cache_dir, provenance)

    def iter_patients(self, patients):
        if patients is None:
            return self
        return iter(patients)

    def patient_from_id(self, id):
        for patient in self:
            if patient.id == id:
                return patient
        raise ValueError("No patient with ID %s found" % id)

    def load_variants(self, patients=None, filter_fn=None, **kwargs):
        """Load a dictionary of patient_id to varcode.VariantCollection

        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        filter_fn : function
            Takes a FilterableVariant and returns a boolean. Only variants returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        merged_variants
            Dictionary of patient_id to VariantCollection
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        filter_fn_name = get_function_name(filter_fn)
        logger.debug("loading variants with filter_fn: {}".format(filter_fn_name))
        patient_variants = {}

        for patient in self.iter_patients(patients):
            variants = self._load_single_patient_variants(patient, filter_fn, **kwargs)
            if variants is not None:
                patient_variants[patient.id] = variants
        return patient_variants


    def _load_single_patient_variants(self, patient, filter_fn, use_cache=True, **kwargs):
        """ Load filtered, merged variants for a single patient, optionally using cache

            Note that filtered variants are first merged before filtering, and
                each step is cached independently. Turn on debug statements for more
                details about cached files.

            Use `_load_single_patient_merged_variants` to see merged variants without filtering.
        """
        if filter_fn is None:
            use_filtered_cache = False
        else:
            filter_fn_name = get_function_name(filter_fn)
            logger.debug("loading variants for patient {} with filter_fn {}".format(patient.id, filter_fn_name))
            use_filtered_cache = use_cache

        ## confirm that we can get cache-name (else don't use filtered cache)
        if use_filtered_cache:
            cache_logger.debug("identifying filtered-cache file name")
            try:
                ## try to load filtered variants from cache
                filtered_cache_file_name = "%s-variants.%s.pkl" % (self.merge_type,
                                                                   hash_function(filter_fn, **kwargs))
            except:
                cache_logger.warning("error identifying filtered-cache file name for patient {}: {}".format(
                        patient.id, filter_fn_name))
                use_filtered_cache = False
            else:
                cache_logger.debug("trying to load filtered variants from cache: {}".format(filtered_cache_file_name))
                try:
                    cached = self.load_from_cache(self.cache_names["variant"], patient.id, filtered_cache_file_name)
                    if cached is not None:
                        return cached
                except:
                    cache_logger.warning("Error loading variants from cache for patient: {}".format(patient.id))
                    pass

        ## get merged variants
        logger.debug("... getting merged variants for: {}".format(patient.id))
        merged_variants = self._load_single_patient_merged_variants(patient, use_cache=use_cache)

        # Note None here is different from 0. We want to preserve None
        if merged_variants is None:
            logger.info("Variants did not exist for patient %s" % patient.id)
            return None

        logger.debug("... applying filters to variants for: {}".format(patient.id))
        filtered_variants = filter_variants(variant_collection=merged_variants,
                                            patient=patient,
                                            filter_fn=filter_fn,
                                            **kwargs)
        if use_filtered_cache:
            cache_logger.debug("saving filtered variants to cache: {}".format(filtered_cache_file_name))
            self.save_to_cache(filtered_variants, self.cache_names["variant"], patient.id, filtered_cache_file_name)
        return filtered_variants

    def _load_single_patient_merged_variants(self, patient, use_cache=True):
        """ Load merged variants for a single patient, optionally using cache

            Note that merged variants are not filtered.
            Use `_load_single_patient_variants` to get filtered variants
        """
        logger.debug("loading merged variants for patient {}".format(patient.id))
        no_variants = False
        try:
            # get merged-variants from cache
            if use_cache:
                ## load unfiltered variants into list of collections
                variant_cache_file_name = "%s-variants.pkl" % (self.merge_type)
                merged_variants = self.load_from_cache(self.cache_names["variant"], patient.id, variant_cache_file_name)
                if merged_variants is not None:
                    return merged_variants
            # get variant collections from file
            variant_collections = []
            optional_maf_cols = ["t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count"]
            if self.additional_maf_cols is not None:
                optional_maf_cols.extend(self.additional_maf_cols)
            for patient_variants in patient.variants_list:
                if type(patient_variants) == str:
                    if ".vcf" in patient_variants:
                        try:
                            variant_collections.append(varcode.load_vcf_fast(patient_variants))
                        # StopIteration is thrown for empty VCFs. For an empty VCF, don't append any variants,
                        # and don't throw an error. But do record a warning, in case the StopIteration was
                        # thrown for another reason.
                        except StopIteration as e:
                            logger.warning("Empty VCF (or possibly a VCF error) for patient {}: {}".format(
                                patient.id, str(e)))
                    elif ".maf" in patient_variants:
                        # See variant_stats.maf_somatic_variant_stats
                        variant_collections.append(
                            varcode.load_maf(
                                patient_variants,
                                optional_cols=optional_maf_cols,
                                encoding="latin-1"))
                    else:
                        raise ValueError("Don't know how to read %s" % patient_variants)
                elif type(patient_variants) == VariantCollection:
                    variant_collections.append(patient_variants)
                else:
                    raise ValueError("Don't know how to read %s" % patient_variants)
            # merge variant-collections
            if len(variant_collections) == 0:
                no_variants = True
            elif len(variant_collections) == 1:
                # There is nothing to merge
                variants = variant_collections[0]
                merged_variants = variants
            else:
                merged_variants = self._merge_variant_collections(variant_collections, self.merge_type)
        except IOError:
            logger.error("Error reading from variant files for patient {}".format(patient.id))
            no_variants = True

        # Note that this is the number of variant collections and not the number of
        # variants. 0 variants will lead to 0 neoantigens, for example, but 0 variant
        # collections will lead to NaN variants and neoantigens.
        if no_variants:
            print("Variants did not exist for patient %s" % patient.id)
            merged_variants = None

        # save merged variants to file
        if use_cache:
            self.save_to_cache(merged_variants, self.cache_names["variant"], patient.id, variant_cache_file_name)
        return merged_variants

    def _merge_variant_collections(self, variant_collections, merge_type):
        logger.debug("Merging variants using merge type: {}".format(merge_type))
        assert merge_type in ["union", "intersection"], "Unknown merge type: %s" % merge_type
        head = variant_collections[0]
        if merge_type == "union":
            merged_variants = head.union(*variant_collections[1:])
        elif merge_type == "intersection":
            merged_variants = head.intersection(*variant_collections[1:])

        return merged_variants

    def load_polyphen_annotations(self, as_dataframe=False,
                                  filter_fn=None):
        """Load a dataframe containing polyphen2 annotations for all variants

        Parameters
        ----------
        database_file : string, sqlite
            Path to the WHESS/Polyphen2 SQLite database.
            Can be downloaded and bunzip2"ed from http://bit.ly/208mlIU
        filter_fn : function
            Takes a FilterablePolyphen and returns a boolean.
            Only annotations returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        annotations
            Dictionary of patient_id to a DataFrame that contains annotations
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        patient_annotations = {}
        for patient in self:
            annotations = self._load_single_patient_polyphen(
                patient,
                filter_fn=filter_fn)
            if annotations is not None:
                annotations["patient_id"] = patient.id
                patient_annotations[patient.id] = annotations
        if as_dataframe:
            return pd.concat(patient_annotations.values())
        return patient_annotations

    def _load_single_patient_polyphen(self, patient, filter_fn):
        cache_name = self.cache_names["polyphen"]
        cached_file_name = "polyphen-annotations.csv"

        # Don't filter here, as these variants are used to generate the
        # PolyPhen cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient,
                                                      filter_fn=None)
        if variants is None:
            return None

        cached = self.load_from_cache(cache_name, patient.id, cached_file_name)
        if cached is not None:
            return filter_polyphen(polyphen_df=cached,
                                   variant_collection=variants,
                                   patient=patient,
                                   filter_fn=filter_fn)

        engine = create_engine("sqlite:///{}".format(self.polyphen_dump_path))
        conn = engine.connect()

        df = pd.DataFrame(columns=["chrom", "pos", "ref", "alt",
                                   "annotation_found", "gene", "protein",
                                   "aa_change", "hvar_pred", "hvar_prob",
                                   "hdiv_pred", "hdiv_prob"])
        for variant in variants:
            chrom = "chr{}".format(getattr(variant, "contig", None))
            pos = getattr(variant, "start", None)
            ref = getattr(variant, "ref", None)
            alt = getattr(variant, "alt",  None)
            annotation = vap.annotate_variant(conn, chrom, pos, ref, alt)
            datum = {"chrom": chrom,
                     "pos": pos,
                     "ref": ref,
                     "alt": alt,
                     "annotation_found": annotation is not None}
            attributes = ["gene", "protein", "aa_change",
                          "hvar_pred", "hvar_prob",
                          "hdiv_pred", "hdiv_prob"]
            for attr in attributes:
                datum[attr] = getattr(annotation, attr, None)
            df = df.append(datum, ignore_index=True)
        df["pos"] = df["pos"].astype("int")
        df["annotation_found"] = df["annotation_found"].astype("bool")
        self.save_to_cache(df, cache_name, patient.id, cached_file_name)
        return filter_polyphen(polyphen_df=df,
                               variant_collection=variants,
                               patient=patient,
                               filter_fn=filter_fn)

    def load_effects(self, patients=None, only_nonsynonymous=False,
                     all_effects=False, filter_fn=None, **kwargs):
        """Load a dictionary of patient_id to varcode.EffectCollection

        Note that this only loads one effect per variant.

        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        only_nonsynonymous : bool, optional
            If true, load only nonsynonymous effects, default False
        all_effects : bool, optional
            If true, return all effects rather than only the top-priority effect per variant
        filter_fn : function
            Takes a FilterableEffect and returns a boolean. Only effects returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        effects
             Dictionary of patient_id to varcode.EffectCollection
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        filter_fn_name = get_function_name(filter_fn)
        logger.debug("loading effects with filter_fn {}".format(filter_fn_name))
        patient_effects = {}
        for patient in self.iter_patients(patients):
            effects = self._load_single_patient_effects(
                patient, only_nonsynonymous, all_effects, filter_fn, **kwargs)
            if effects is not None:
                patient_effects[patient.id] = effects
        return patient_effects

    def _load_single_patient_effects_unfiltered(self, patient, only_nonsynonymous):
        """
        Load single patient effects without applying a filter-function
        """
        cached_file_name = "%s-effects.pkl" % self.merge_type

        # get result from cache
        if only_nonsynonymous:
            cached = self.load_from_cache(self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["effect"], patient.id, cached_file_name)
        if cached is not None:
            return cached

        # prepare various cached effects
        # Don't filter here, as these variants are used to generate the
        # effects cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient, filter_fn=None)
        if variants is None:
            return None

        ### all effects per variant (not the top priority)
        effects = variants.effects()
        self.save_to_cache(effects, self.cache_names["effect"], patient.id, cached_file_name)

        # Save all nonsynonymous effects, rather than top priority only.
        nonsynonymous_effects = effects.drop_silent_and_noncoding()
        self.save_to_cache(nonsynonymous_effects, self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)
        
        if only_nonsynonymous:
            return nonsynonymous_effects
        else:
            return effects

    def _load_single_patient_effects(self, patient, only_nonsynonymous, all_effects, filter_fn, **kwargs):
        filter_fn_name = get_function_name(filter_fn)
        logger.debug("loading effects for patient {} with filter_fn {}".format(patient.id, filter_fn_name))

        cached_file_name = "{}-effects.{}.pkl".format(self.merge_type, hash_function(filter_fn, **kwargs))
        cache_logger.debug("filtered effects cache name set to: {}".format(cached_file_name))

        # try to get filtered result from cache
        cached = self.load_from_cache(self.cache_names["effect"], patient.id, cached_file_name)
        if cached is not None:
            return cached

        # get effects in absence of filter
        effects = self._load_single_patient_effects_unfiltered(
            patient=patient,
            only_nonsynonymous=only_nonsynonymous)

        # (needed for metadata - should ideally be included in the `effects` object)
        variants = self._load_single_patient_variants(
            patient, filter_fn=None)
        
        if effects is None:
            return None

        # filter effects
        filtered_effects = filter_effects(
            effect_collection=(effects),
            patient=patient,
            variant_collection=variants,
            filter_fn=filter_fn,
            all_effects=all_effects,
            **kwargs)

        # save to cache
        self.save_to_cache(filtered_effects, self.cache_names["effect"], patient.id, cached_file_name)
        return filtered_effects

    def load_kallisto(self):
        """
        Load Kallisto transcript quantification data for a cohort

        Parameters
        ----------

        Returns
        -------
        kallisto_data : Pandas dataframe
            Pandas dataframe with Kallisto data for all patients
            columns include patient_id, gene_name, est_counts
        """
        kallisto_data = pd.concat(
            [self._load_single_patient_kallisto(patient) for patient in self],
            copy=False
        )

        if self.kallisto_ensembl_version is None:
            raise ValueError("Required a kallisto_ensembl_version but none was specified")

        ensembl_release = cached_release(self.kallisto_ensembl_version)

        kallisto_data["gene_name"] = \
            kallisto_data["target_id"].map(lambda t: ensembl_release.gene_name_of_transcript_id(t))

        # sum counts across genes
        kallisto_data = \
            kallisto_data.groupby(["patient_id", "gene_name"])[["est_counts"]].sum().reset_index()

        return kallisto_data

    def _load_single_patient_kallisto(self, patient):
        """
        Load Kallisto gene quantification given a patient

        Parameters
        ----------
        patient : Patient

        Returns
        -------
        data: Pandas dataframe
            Pandas dataframe of sample's Kallisto data
            columns include patient_id, target_id, length, eff_length, est_counts, tpm
        """
        data = pd.read_csv(patient.tumor_sample.kallisto_path, sep="\t")
        data["patient_id"] = patient.id
        return data

    def load_cufflinks(self, filter_ok=True):
        """
        Load a Cufflinks gene expression data for a cohort

        Parameters
        ----------
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == "OK"

        Returns
        -------
        cufflinks_data : Pandas dataframe
            Pandas dataframe with Cufflinks data for all patients
            columns include patient_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        return \
            pd.concat(
                [self._load_single_patient_cufflinks(patient, filter_ok) for patient in self],
                copy=False
        )

    def _load_single_patient_cufflinks(self, patient, filter_ok):
        """
        Load Cufflinks gene quantification given a patient

        Parameters
        ----------
        patient : Patient
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == "OK"

        Returns
        -------
        data: Pandas dataframe
            Pandas dataframe of sample's Cufflinks data
            columns include patient_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        data = pd.read_csv(patient.tumor_sample.cufflinks_path, sep="\t")
        data["patient_id"] = patient.id

        if filter_ok:
            # Filter to OK FPKM counts
            data = data[data["FPKM_status"] == "OK"]
        return data

    def load_neoantigens(self, patients=None, only_expressed=False,
                         epitope_lengths=[8, 9, 10, 11], ic50_cutoff=500,
                         process_limit=10, max_file_records=None,
                         filter_fn=None):
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)

        dfs = {}
        for patient in self.iter_patients(patients):
            df_epitopes = self._load_single_patient_neoantigens(
                patient=patient,
                only_expressed=only_expressed,
                epitope_lengths=epitope_lengths,
                ic50_cutoff=ic50_cutoff,
                process_limit=process_limit,
                max_file_records=max_file_records,
                filter_fn=filter_fn)
            if df_epitopes is not None:
                dfs[patient.id] = df_epitopes
        return dfs

    def _load_single_patient_neoantigens(self, patient, only_expressed, epitope_lengths,
                                         ic50_cutoff, process_limit, max_file_records,
                                         filter_fn):
        cached_file_name = "%s-neoantigens.csv" % self.merge_type

        # Don't filter here, as these variants are used to generate the
        # neoantigen cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient, filter_fn=None)
        if variants is None:
            return None

        if patient.hla_alleles is None:
            print("HLA alleles did not exist for patient %s" % patient.id)
            return None

        if only_expressed:
            cached = self.load_from_cache(self.cache_names["expressed_neoantigen"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["neoantigen"], patient.id, cached_file_name)
        if cached is not None:
            return filter_neoantigens(neoantigens_df=cached,
                                      variant_collection=variants,
                                      patient=patient,
                                      filter_fn=filter_fn)

        try:
            mhc_model = self.mhc_class(
                alleles=patient.hla_alleles,
                epitope_lengths=epitope_lengths,
                max_file_records=max_file_records,
                process_limit=process_limit)
        except TypeError:
            # The class may not support max_file_records and process_limit.
            mhc_model = self.mhc_class(
                alleles=patient.hla_alleles,
                epitope_lengths=epitope_lengths)

        if only_expressed:
            df_isovar = self.load_single_patient_isovar(patient=patient,
                                                         variants=variants,
                                                         epitope_lengths=epitope_lengths)

            # Map from isovar rows to protein sequences
            isovar_rows_to_protein_sequences = dict([
                (frozenset(row.to_dict().items()), row["amino_acids"]) for (i, row) in df_isovar.iterrows()])

            # MHC binding prediction
            epitopes = mhc_model.predict(isovar_rows_to_protein_sequences)

            # Call `get_filtered_isovar_epitopes` in order to only include peptides that
            # overlap a variant; without this filter, when we use
            # protein_sequence_length above, some 8mers generated from a 21mer source will
            # not overlap a variant.
            df_epitopes = self.get_filtered_isovar_epitopes(
                epitopes, ic50_cutoff=ic50_cutoff).dataframe()
            # Store chr/pos/ref/alt in the cached DataFrame so we can filter based on
            # the variant later.
            for variant_column in ["chr", "pos", "ref", "alt"]:
                # Be consistent with Topiary's output of "start" rather than "pos".
                # Isovar, on the other hand, outputs "pos".
                # See https://github.com/hammerlab/topiary/blob/5c12bab3d47bd86d396b079294aff141265f8b41/topiary/converters.py#L50
                df_column = "start" if variant_column == "pos" else variant_column
                df_epitopes[df_column] = df_epitopes.source_sequence_key.apply(
                    lambda key: dict(key)[variant_column])
            df_epitopes["patient_id"] = patient.id

            self.save_to_cache(df_epitopes, self.cache_names["expressed_neoantigen"], patient.id, cached_file_name)
        else:
            epitopes = predict_epitopes_from_variants(
                variants=variants,
                mhc_model=mhc_model,
                ic50_cutoff=ic50_cutoff,
                # Only include peptides with a variant
                only_novel_epitopes=True)
            df_epitopes = epitopes_to_dataframe(epitopes)
            df_epitopes["patient_id"] = patient.id

            self.save_to_cache(df_epitopes, self.cache_names["neoantigen"], patient.id, cached_file_name)

        return filter_neoantigens(neoantigens_df=df_epitopes,
                                  variant_collection=variants,
                                  patient=patient,
                                  filter_fn=filter_fn)

    def get_filtered_isovar_epitopes(self, epitopes, ic50_cutoff):
        """
        Mostly replicates topiary.build_epitope_collection_from_binding_predictions

        Note: topiary needs to do fancy stuff like subsequence_protein_offset + binding_prediction.offset
        in order to figure out whether a variant is in the peptide because it only has the variant's
        offset into the full protein; but isovar gives us the variant's offset into the protein subsequence
        (dictated by protein_sequence_length); so all we need to do is map that onto the smaller 8-11mer
        peptides generated by mhctools.
        """
        mutant_binding_predictions = []
        for binding_prediction in epitopes:
            peptide = binding_prediction.peptide
            peptide_offset = binding_prediction.offset
            isovar_row = dict(binding_prediction.source_sequence_key)
            is_mutant = contains_mutant_residues(
                peptide_start_in_protein=peptide_offset,
                peptide_length=len(peptide),
                mutation_start_in_protein=isovar_row["variant_aa_interval_start"],
                mutation_end_in_protein=isovar_row["variant_aa_interval_end"])
            if is_mutant and binding_prediction.value <= ic50_cutoff:
                mutant_binding_predictions.append(binding_prediction)
        return EpitopeCollection(mutant_binding_predictions)

    def load_single_patient_isovar(self, patient, variants, epitope_lengths):
        # TODO: different epitope lengths, and other parameters, should result in
        # different caches
        isovar_cached_file_name = "%s-isovar.csv" % self.merge_type
        df_isovar = self.load_from_cache(self.cache_names["isovar"], patient.id, isovar_cached_file_name)
        if df_isovar is not None:
            return df_isovar

        import logging
        logging.disable(logging.INFO)
        if patient.tumor_sample is None:
            raise ValueError("Patient %s has no tumor sample" % patient.id)
        if patient.tumor_sample.bam_path_rna is None:
            raise ValueError("Patient %s has no tumor RNA BAM path" % patient.id)
        rna_bam_file = AlignmentFile(patient.tumor_sample.bam_path_rna)

        # To ensure that e.g. 8-11mers overlap substitutions, we need at least this
        # sequence length: (max peptide length * 2) - 1
        # Example:
        # 123456789AB
        #           123456789AB
        # AAAAAAAAAAVAAAAAAAAAA
        protein_sequence_length = (max(epitope_lengths) * 2) - 1
        allele_reads_generator = reads_overlapping_variants(
            variants=variants,
            samfile=rna_bam_file,
            min_mapping_quality=1)
        protein_sequences_generator = reads_generator_to_protein_sequences_generator(
            allele_reads_generator,
            protein_sequence_length=protein_sequence_length,
            # Per Alex R.'s suggestion; equivalent to min_reads_supporting_rna_sequence previously
            min_variant_sequence_coverage=3,
            max_protein_sequences_per_variant=1, # Otherwise we might have too much neoepitope diversity
            variant_sequence_assembly=False)
        df_isovar = protein_sequences_generator_to_dataframe(protein_sequences_generator)
        self.save_to_cache(df_isovar, self.cache_names["isovar"], patient.id, isovar_cached_file_name)
        return df_isovar

    def load_ensembl_coverage(self):
        if self.pageant_coverage_path is None:
            raise ValueError("Need a Pageant CoverageDepth path to load ensembl coverage values")
        return variant_filters.load_ensembl_coverage(
            cohort=self,
            coverage_path=self.pageant_coverage_path,
            min_normal_depth=self.min_coverage_normal_depth,
            min_tumor_depth=self.min_coverage_tumor_depth,
            pageant_dir_fn=self.pageant_dir_fn)

    def clear_caches(self):
        for cache in self.cache_names.keys():
            self.clear_cache(cache)

    def clear_cache(self, cache):
        cache_path = path.join(self.cache_dir, self.cache_names[cache])
        if path.exists(cache_path):
            rmtree(cache_path)

    def cohort_columns(self):
        cohort_dataframe = self.as_dataframe()
        column_types = [cohort_dataframe[col].dtype for col in cohort_dataframe.columns]
        return dict(zip(list(cohort_dataframe.columns), column_types))

    def plot_col_from_cols(self, cols, only_allow_one=False, plot_col=None):
        if type(cols) == str:
            if plot_col is not None:
                raise ValueError("plot_col is specified when it isn't needed because there is only one col.")
            plot_col = cols
        elif type(cols) == list:
            # If e.g. an `on` dictionary is provided, that'll result in a list of cols.
            # But if there is just one col, we can use it as the plot_col.
            if len(cols) == 0:
                raise ValueError("Empty list of `on` cols: %s" % str(cols))
            elif len(cols) == 1:
                plot_col = cols[0]
            else:
                if only_allow_one:
                    raise ValueError("`on` has multiple columns, which is not allowed here.")
                if plot_col is None:
                    raise ValueError("plot_col must be specified when multiple `on`s are present.")
        else:
            raise ValueError("cols need to be a str or a list, but cols are %s" % str(cols))
        return plot_col

    def plot_roc_curve(self, on, bootstrap_samples=100, ax=None, **kwargs):
        """Plot an ROC curve for benefit and a given variable

        Parameters
        ----------
        on : str or function or list or dict
            See `cohort.load.as_dataframe`
        bootstrap_samples : int, optional
            Number of boostrap samples to use to compute the AUC
        ax : Axes, default None
            Axes to plot on

        Returns
        -------
        (mean_auc_score, plot): (float, matplotlib plot)
            Returns the average AUC for the given predictor over `bootstrap_samples`
            and the associated ROC curve
        """
        plot_col, df = self.as_dataframe(on, return_cols=True, **kwargs)
        df = filter_not_null(df, "benefit")
        df = filter_not_null(df, plot_col)
        df.benefit = df.benefit.astype(bool)
        return roc_curve_plot(df, plot_col, "benefit", bootstrap_samples, ax=ax)

    def plot_benefit(self, on, benefit_col="benefit", label="Response", ax=None,
                     alternative="two-sided", boolean_value_map={},
                     order=None, **kwargs):
        """Plot a comparison of benefit/response in the cohort on a given variable
        """
        no_benefit_plot_name = "No %s" % self.benefit_plot_name
        boolean_value_map = boolean_value_map or {True: self.benefit_plot_name, False: no_benefit_plot_name}
        order = order or [no_benefit_plot_name, self.benefit_plot_name]

        return self.plot_boolean(on=on,
                                 boolean_col=benefit_col,
                                 alternative=alternative,
                                 boolean_label=label,
                                 boolean_value_map=boolean_value_map,
                                 order=order,
                                 ax=ax,
                                 **kwargs)

    def plot_boolean(self,
                     on,
                     boolean_col,
                     plot_col=None,
                     boolean_label=None,
                     boolean_value_map={},
                     order=None,
                     ax=None,
                     alternative="two-sided",
                     **kwargs):
        """Plot a comparison of `boolean_col` in the cohort on a given variable via
        `on` or `col`.

        If the variable (through `on` or `col`) is binary this will compare
        odds-ratios and perform a Fisher's exact test.

        If the variable is numeric, this will compare the distributions through
        a Mann-Whitney test and plot the distributions with box-strip plot

        Parameters
        ----------
        on : str or function or list or dict
            See `cohort.load.as_dataframe`
        plot_col : str, optional
            If on has many columns, this is the one whose values we are plotting.
            If on has a single column, this is unnecessary.
            We might want many columns if, e.g. we're generating boolean_col from a
            function as well.
        boolean_col : str
            Column name of boolean column to plot or compare against.
        boolean_label : None, optional
            Label to give boolean column in the plot
        boolean_value_map : dict, optional
            Map of conversions for values in the boolean column, i.e. {True: 'High', False: 'Low'}
        order : None, optional
            Order of the labels on the x-axis
        ax : None, optional
            Axes to plot on
        alternative : str, optional
            Choose the sidedness of the mannwhitneyu or Fisher's Exact test.

        Returns
        -------
        (Test statistic, p-value): (float, float)

        """
        cols, df = self.as_dataframe(on, return_cols=True, **kwargs)
        plot_col = self.plot_col_from_cols(cols=cols, plot_col=plot_col)
        df = filter_not_null(df, boolean_col)
        df = filter_not_null(df, plot_col)

        if boolean_label:
            df[boolean_label] = df[boolean_col]
            boolean_col = boolean_label

        condition_value = None
        if boolean_value_map:
            assert set(boolean_value_map.keys()) == set([True, False]), \
                "Improper mapping of boolean column provided"
            df[boolean_col] = df[boolean_col].map(lambda v: boolean_value_map[v])
            condition_value = boolean_value_map[True]

        if df[plot_col].dtype == "bool":
            results = fishers_exact_plot(
                data=df,
                condition1=boolean_col,
                condition2=plot_col,
                condition1_value=condition_value,
                alternative=alternative,
                order=order,
                ax=ax)
        else:
            results = mann_whitney_plot(
                data=df,
                condition=boolean_col,
                distribution=plot_col,
                condition_value=condition_value,
                alternative=alternative,
                order=order,
                ax=ax)
        return results

    def plot_survival(self,
                      on,
                      how="os",
                      survival_units="Days",
                      strata=None,
                      ax=None,
                      ci_show=False,
                      with_condition_color="#B38600",
                      no_condition_color="#A941AC",
                      with_condition_label=None,
                      no_condition_label=None,
                      color_map=None,
                      label_map=None,
                      color_palette="Set2",
                      threshold=None, **kwargs):
        """Plot a Kaplan Meier survival curve by splitting the cohort into two groups
        Parameters
        ----------
        on : str or function or list or dict
            See `cohort.load.as_dataframe`
        how : {"os", "pfs"}, optional
            Whether to plot OS (overall survival) or PFS (progression free survival)
        survival_units : str
            Unit of time for the survival measure, i.e. Days or Months
        strata : str
            (optional) column name of stratifying variable
        ci_show : bool
            Display the confidence interval around the survival curve
        threshold : int, "median", "median-per-strata" or None (optional)
            Threshold of `col` on which to split the cohort
        """
        assert how in ["os", "pfs"], "Invalid choice of survival plot type %s" % how
        cols, df = self.as_dataframe(on, return_cols=True, **kwargs)
        plot_col = self.plot_col_from_cols(cols=cols, only_allow_one=True)
        df = filter_not_null(df, plot_col)
        results = plot_kmf(
            df=df,
            condition_col=plot_col,
            xlabel=survival_units,
            ylabel="Overall Survival (%)" if how == "os" else "Progression-Free Survival (%)",
            censor_col="deceased" if how == "os" else "progressed_or_deceased",
            survival_col=how,
            strata_col=strata,
            threshold=threshold,
            ax=ax,
            ci_show=ci_show,
            with_condition_color=with_condition_color,
            no_condition_color=no_condition_color,
            with_condition_label=with_condition_label,
            no_condition_label=no_condition_label,
            color_palette=color_palette,
            label_map=label_map,
            color_map=color_map,
        )
        return results

    def plot_correlation(self, on, x_col=None, plot_type="jointplot", stat_func=pearsonr, show_stat_func=True, plot_kwargs={}, **kwargs):
        """Plot the correlation between two variables.

        Parameters
        ----------
        on : list or dict of functions or strings
            See `cohort.load.as_dataframe`
        x_col : str, optional
            If `on` is a dict, this guarantees we have the expected ordering.
        plot_type : str, optional
            Specify "jointplot", "regplot", "boxplot", or "barplot".
        stat_func : function, optional.
            Specify which function to use for the statistical test.
        show_stat_func : bool, optional
            Whether or not to show the stat_func result in the plot itself.
        plot_kwargs : dict, optional
            kwargs to pass through to plotting functions.
        """
        if plot_type not in ["boxplot", "barplot", "jointplot", "regplot"]:
            raise ValueError("Invalid plot_type %s" % plot_type)
        plot_cols, df = self.as_dataframe(on, return_cols=True, **kwargs)
        if len(plot_cols) != 2:
            raise ValueError("Must be comparing two columns, but there are %d columns" % len(plot_cols))
        for plot_col in plot_cols:
            df = filter_not_null(df, plot_col)
        if x_col is None:
            x_col = plot_cols[0]
            y_col = plot_cols[1]
        else:
            if x_col == plot_cols[0]:
                y_col = plot_cols[1]
            else:
                y_col = plot_cols[0]
        series_x = df[x_col]
        series_y = df[y_col]
        coeff, p_value = stat_func(series_x, series_y)
        if plot_type == "jointplot":
            plot = sb.jointplot(data=df, x=x_col, y=y_col,
                                stat_func=stat_func if show_stat_func else None,
                                **plot_kwargs)
        elif plot_type == "regplot":
            plot = sb.regplot(data=df, x=x_col, y=y_col,
                              **plot_kwargs)
        elif plot_type == "boxplot":
            plot = stripboxplot(data=df, x=x_col, y=y_col, **plot_kwargs)
        else:
            plot = sb.barplot(data=df, x=x_col, y=y_col, **plot_kwargs)
        return CorrelationResults(coeff=coeff, p_value=p_value, stat_func=stat_func,
                                  series_x=series_x, series_y=series_y, plot=plot)

    def coxph(self, on, formula=None, how="pfs"):
        return cohort_coxph(self, on, formula=formula, how=how)

    def bootstrap_auc(self, on, pred_col="is_benefit", n_bootstrap=1000, **kwargs):
        return cohort_bootstrap_auc(self, on, pred_col=pred_col, n_bootstrap=n_bootstrap)

    def mean_bootstrap_auc(self, on, pred_col="is_benefit", n_bootstrap=1000, **kwargs):
        return cohort_mean_bootstrap_auc(self, on, pred_col=pred_col, n_bootstrap=n_bootstrap)

    def _list_patient_ids(self):
        """ Utility function to return a list of patient ids in the Cohort
        """
        results = []
        for patient in self:
            results.append(patient.id)
        return(results)

    def summarize_provenance_per_cache(self):
        """Utility function to summarize provenance files for cached items used by a Cohort,
        for each cache_dir that exists. Only existing cache_dirs are summarized.

        This is a summary of provenance files because the function checks to see whether all
        patients data have the same provenance within the cache dir. The function assumes
        that it will be desireable to have all patients data generated using the same
        environment, for each cache type.

        At the moment, most PROVENANCE files contain details about packages used to generat
        e the cached data file. However, this function is generic & so it summarizes the
        contents of those files irrespective of their contents.

        Returns
        ----------
        Dict containing summarized provenance for each existing cache_dir, after checking
        to see that provenance files are identical among all patients in the data frame for
        that cache_dir.

        If conflicting PROVENANCE files are discovered within a cache-dir:
         - a warning is generated, describing the conflict
         - and, a value of `None` is returned in the dictionary for that cache-dir

        See also
        -----------
        * `?cohorts.Cohort.summarize_provenance` which summarizes provenance files among
        cache_dirs.
        * `?cohorts.Cohort.summarize_dataframe` which hashes/summarizes contents of the data
        frame for this cohort.
        """
        provenance_summary = {}
        df = self.as_dataframe()
        for cache in self.cache_names:
            cache_name = self.cache_names[cache]
            cache_provenance = None
            num_discrepant = 0
            this_cache_dir = path.join(self.cache_dir, cache_name)
            if path.exists(this_cache_dir):
                for patient_id in self._list_patient_ids():
                    patient_cache_dir = path.join(this_cache_dir, patient_id)
                    try:
                        this_provenance = self.load_provenance(patient_cache_dir = patient_cache_dir)
                    except:
                        this_provenance = None
                    if this_provenance:
                        if not(cache_provenance):
                            cache_provenance = this_provenance
                        else:
                            num_discrepant += compare_provenance(this_provenance, cache_provenance)
                if num_discrepant == 0:
                    provenance_summary[cache_name] = cache_provenance
                else:
                    provenance_summary[cache_name] = None
        return(provenance_summary)

    def summarize_dataframe(self):
        """Summarize default dataframe for this cohort using a hash function.
        Useful for confirming the version of data used in various reports, e.g. ipynbs
        """
        if self.dataframe_hash:
            return(self.dataframe_hash)
        else:
            df = self._as_dataframe_unmodified()
            return(self.dataframe_hash)

    def summarize_provenance(self):
        """Utility function to summarize provenance files for cached items used by a Cohort.

        At the moment, most PROVENANCE files contain details about packages used to
        generate files. However, this function is generic & so it summarizes the contents
        of those files irrespective of their contents.

        Returns
        ----------
        Dict containing summary of provenance items, among all cache dirs used by the Cohort.

        IE if all provenances are identical across all cache dirs, then a single set of
        provenances is returned. Otherwise, if all provenances are not identical, the provenance
        items per cache_dir are returned.

        See also
        ----------
        `?cohorts.Cohort.summarize_provenance_per_cache` which is used to summarize provenance
        for each existing cache_dir.
        """
        provenance_per_cache = self.summarize_provenance_per_cache()
        summary_provenance = None
        num_discrepant = 0
        for cache in provenance_per_cache:
            if not(summary_provenance):
                ## pick arbitrary provenance & call this the "summary" (for now)
                summary_provenance = provenance_per_cache[cache]
                summary_provenance_name = cache
            ## for each cache, check equivalence with summary_provenance
            num_discrepant += compare_provenance(
                provenance_per_cache[cache],
                summary_provenance,
                left_outer_diff = "In %s but not in %s" % (cache, summary_provenance_name),
                right_outer_diff = "In %s but not in %s" % (summary_provenance_name, cache)
            )
        ## compare provenance across cached items
        if num_discrepant == 0:
            prov = summary_provenance ## report summary provenance if exists
        else:
            prov = provenance_per_cache ## otherwise, return provenance per cache
        return(prov)

    def summarize_data_sources(self):
        """Utility function to summarize data source status for this Cohort, useful for confirming
        the state of data used for an analysis

        Returns
        ----------
        Dictionary with summary of data sources

        Currently contains
        - dataframe_hash: hash of the dataframe (see `?cohorts.Cohort.summarize_dataframe`)
        - provenance_file_summary: summary of provenance file contents (see `?cohorts.Cohort.summarize_provenance`)
        """
        provenance_file_summary = self.summarize_provenance()
        dataframe_hash = self.summarize_dataframe()
        results = {
            "provenance_file_summary": provenance_file_summary,
            "dataframe_hash": dataframe_hash
        }
        return(results)
