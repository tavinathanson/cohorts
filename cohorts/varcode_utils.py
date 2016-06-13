from varcode import VariantCollection

def filter_variants_with_metadata(variant_collection, filter_fn):
    """Filter variants from the Variant Collections

    Parameters
    ----------
    variant_collection : varcode.VariantCollection
    filter_fn: function
        Takes a variant and it's metadata and returns a boolean. Only variants returning True are preserved.

    Returns
    -------
    varcode.VariantCollection
        Filtered variant collection, with only the variants passing the filter
    """
    if filter_fn:
        return VariantCollection([variant
                                  for variant in variant_collection
                                  if filter_fn(variant, variant_collection.metadata[variant])
                                  ])
    else:
        return variant_collection
