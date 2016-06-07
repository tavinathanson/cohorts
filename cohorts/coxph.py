
import patsy
import pandas
import lifelines
import numpy

def coxph_model(data, formula, time_col, event_col, penalizer = 0.1, normalize = False):
    sdata = patsy.dmatrix(
                formula,
                 data = data,
                 return_type = 'dataframe'
                 ).join(data[[time_col, event_col]])
    sdata = sdata.ix[:, sdata.columns != 'Intercept']
    sdata = sdata.dropna()
    
    cf = lifelines.CoxPHFitter(penalizer = penalizer, normalize = normalize)
    cf.fit(sdata,
           time_col,
           event_col
           )
    ##cf.print_summary()
    try:
        concordance = lifelines.utils.concordance_index(
            cf.durations,
            -cf.predict_partial_hazard(cf.data).values.ravel(),
            cf.event_observed
            )
    except:
        print('Error computing concordance')
        concordance = numpy.nan
    
    return(cf, sdata, concordance)

def bootstrap_coxph(models, data, time_col, event_col, bootstrap_samples = 1000, frac = 0.9, **kwargs):
    scores = {}
    for m in models:
        scores[m] = numpy.zeros(bootstrap_samples)
    for i in range(bootstrap_samples):
        sample_data = data.sample(frac = frac, replace = True)

        for m in models:
            a = coxph_model(data = sample_data,
                            formula = models[m],
                            time_col = time_col,
                            event_col = event_col,
                            **kwargs
                           )
            scores[m][i] = a[3]
    results = pandas.DataFrame(scores)
    results['metric'] = 'concordance'
    results = pandas.melt(results, id_vars=['metric'], var_name = 'model', value_name='concordance')
    return(results)


