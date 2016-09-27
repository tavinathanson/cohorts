from varcode import EffectCollection, Variant

def genome(variant_collection):
    return variant_collection[0].ensembl

class FilterableVariant(object):
    def __init__(self, variant, variant_collection, patient):
        self.variant = variant
        self.variant_collection = variant_collection
        self.patient = patient

    @property
    def variant_metadata(self):
        source_to_metadata = self.variant_collection.source_to_metadata_dict
        return dict([
                (source, source_to_metadata[source][self.variant])
                for source in self.variant_collection.sources
                if self.variant in source_to_metadata[source]
            ])

    @property
    def genome(self):
        return genome(self.variant_collection)

class FilterableEffect(FilterableVariant):
    def __init__(self, effect, variant_collection, patient):
        self.effect = effect
        FilterableVariant.__init__(self,
                                   variant=effect.variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

class FilterableNeoantigen(FilterableVariant):
    def __init__(self, neoantigen_row, variant_collection, patient):
        self.neoantigen_row = neoantigen_row
        def build_variant(row, genome):
            return Variant(
                contig=row["chr"],
                ref=row["ref"],
                alt=row["alt"],
                start=row["start"],
                ensembl=genome)
        variant = build_variant(neoantigen_row, genome(variant_collection))
        FilterableVariant.__init__(self,
                                   variant=variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

class FilterablePolyphen(FilterableVariant):
    def __init__(self, polyphen_row, variant_collection, patient):
        self.polyphen_row = polyphen_row
        def build_variant(row, genome):
            return Variant(
                contig=row["chrom"],
                ref=row["ref"],
                alt=row["alt"],
                start=row["pos"],
                ensembl=genome)
        variant = build_variant(polyphen_row, genome(variant_collection))
        FilterableVariant.__init__(self,
                                   variant=variant,
                                   variant_collection=variant_collection,
                                   patient=patient)

def filter_variants(variant_collection, patient, filter_fn):
    """Filter variants from the Variant Collection

    Parameters
    ----------
    variant_collection : varcode.VariantCollection
    patient : cohorts.Patient
    filter_fn: function
        Takes a FilterableVariant and returns a boolean. Only variants returning True are preserved.

    Returns
    -------
    varcode.VariantCollection
        Filtered variant collection, with only the variants passing the filter
    """
    if filter_fn:
        return variant_collection.clone_with_new_elements([
            variant
            for variant in variant_collection
            if filter_fn(FilterableVariant(
                    variant=variant,
                    variant_collection=variant_collection,
                    patient=patient))
        ])
    else:
        return variant_collection

def filter_effects(effect_collection, variant_collection, patient, filter_fn):
    """Filter variants from the Effect Collection

    Parameters
    ----------
    effect_collection : varcode.EffectCollection
    variant_collection : varcode.VariantCollection
    patient : cohorts.Patient
    filter_fn: function
        Takes a FilterableEffect and returns a boolean. Only effects returning True are preserved.

    Returns
    -------
    varcode.EffectCollection
        Filtered effect collection, with only the variants passing the filter
    """
    if filter_fn:
        return EffectCollection([
            effect
            for effect in effect_collection
            if filter_fn(FilterableEffect(
                    effect=effect,
                    variant_collection=variant_collection,
                    patient=patient))])
    else:
        return effect_collection

def filter_neoantigens(neoantigens_df, variant_collection, patient, filter_fn):
    if filter_fn:
        filter_mask = neoantigens_df.apply(
            lambda row: filter_fn(
                FilterableNeoantigen(neoantigen_row=row,
                                     variant_collection=variant_collection,
                                     patient=patient)), axis=1)
        return neoantigens_df[filter_mask]
    else:
        return neoantigens_df

def filter_polyphen(polyphen_df, variant_collection, patient, filter_fn):
    if filter_fn:
        filter_mask = polyphen_df.apply(
            lambda row: filter_fn(
                FilterablePolyphen(polyphen_row=row,
                                   variant_collection=variant_collection,
                                   patient=patient)), axis=1)
        return polyphen_df[filter_mask]
    else:
        return polyphen_df
