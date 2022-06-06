from scipy import stats
import numpy as np
from matplotlib import pyplot as plt


def generate_quantile_quantile_plot(observed_tvalues, permutation_tvalues, n_subjects, n_OLS_params):
    """
    Quantile-Quantile modeling is technique commonly used in genome wide association studies (Ehret 2011) to simultaneously assess whether a collection of many statistical tests deviates from the expected null (schweder, spjotvoll 1982).
    Schematically, the procedure is as follows: 
        - compute an association test for each edge of the connectome by fitting a multiple regression model with the phenotype of interest as the outcome variable and edge connectivity weight as the predictor of interest
        - any confound variables that we wish to control for should be included in the edgewise multiple regression models above
        - the observed p-values from the steps above are then negative-log transformed, and then sorted from smallest to largest. 
        - the rank-ordered negative-log transformed p-values are plotted versus negative-log transformed linearly scaled points on the x-axis 
          - ie: the x-y coordinates for the plotted points are (-log(1), -log(largest p-value)), (-log(n-1/n), -log(2nd largest p-value)), … (-log(1/n), -log(smallest p-value))
        - the corresponding null distribution can either be generated theoretically or using permutation t_values using the NCA function network_contingency_analysis
        
    Parameters
    ----------
    observed_tvalues : ndarray
        1-dim numpy array of shape (number_of_edges, ) with observed t-values for phenotype beta, for each edge
    permutation_tvalues : ndarray or None
        2-dim numpy array of shape (n_perms, number_of_edges) with permutation t-values for phenotype beta, for each edge
        if None, then the theoretical null distribution is plotted
    n_subjects: integer
        Number of subjects used in edgewise association tests
    n_OLS_params: integer
        Total number of predictors used in edgewise association tests (should be: intercept (1) + phenotype of interest (1) + number of confound variables)
        
    Returns
    -------
    None:
        Outputs matplotlib plot to stdout. Best to call this function from a jupyter notebook or IDE
    """
    
    def t_to_p(t, n_subjects, n_params):
        return scipy.stats.t.sf(np.abs(t), n_subjects-n_params)*2 

    def qbeta(p, shape1, shape2):
        """
        Calculates the cumulative of the Beta-distribution
        """
        result = stats.beta.ppf(q=p, a=shape1, b=shape2, loc=0, scale=1)
        return result

    observed_pvals = t_to_p(observed_tvalues, n_subjects, n_OLS_params)
    pval = -np.log10(sorted(observed_pvals)) 
    n = len(observed_pvals)
    a = np.arange(1, n+1)
    x = -np.log10(a/n)
    
    if permutation_tvalues is None:
        upper = qbeta(0.025, a, sorted(a)[::-1])
        lower = qbeta(0.975, a, sorted(a)[::-1])
    else:
        print('Converting permutation t-values to p-values. This can take a while for large number of perms...')
        permutation_pvals = t_to_p(permutation_tvalues, n_subjects, n_OLS_params)
        print('Sorting permutation p-values...')
        permutation_pvals = np.sort(permutation_pvals, axis=1)
        upper = np.percentile(permutation_pvals, q=2.5, axis=0)
        lower = np.percentile(permutation_pvals, q=97.5, axis=0)
        
    plt.figure(figsize=(11, 9))
    plt.xlabel(r"$-log_{10}$(Expected P-value)", fontsize=30, fontweight='bold')
    plt.ylabel(r"$-log_{10}$(Observed P-value)", fontsize=30, fontweight='bold')
    plt.plot(np.arange(min(x), max(x), 0.01), np.arange(min(x), max(x), 0.01), c='r', linewidth=2)
    plt.fill_between(x, -np.log10(upper), -np.log10(lower), color='lightgray', zorder=1)
    plt.scatter(x, pval, alpha=1, s=25, zorder=2)
    plt.xticks(np.arange(min(x), max(x)+1, 1), fontsize=20)
    plt.yticks(np.arange(min(pval), max(pval)+1, 1), fontsize=20)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.show()
    
    return None
