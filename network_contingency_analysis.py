import numpy as np
import scipy
import itertools
import statsmodels


def suprathreshold_count_cell_summary_fn(t_vals_per_edge, t_threshold):
    """ Returns the number of suprathreshold edges in t_vals_per_edge. This is the cell summary statistic used in Sripada (2021).
    """
    return np.sum(np.abs(t_vals_per_edge) > abs(t_threshold))


def network_contingency_analysis(edges, phenotype, covariates, n_perms, t_threshold, net_order_file, 
                                 perm_order=None, cell_summary_fn=suprathreshold_count_cell_summary_fn, random_seed=42):
    
    """
    Network contingency analysis (NCA) assesses whether the observed cell summary statistic (eg: count of suprathreshold edges) in a cell is higher than what is expected by chance. 
    Schematically, the procedure is as follows: 
        - per edge linear regression is conducted with the following formula: edge_i ~ intercept + phenotype + covariates
          the t-value for the beta associated with the phenotype in the formula above is stored, for each edge
        - the 2 steps above are repeated with Freedman and Lane permutations (1983) to determine the null/chance t-value distribution for the phenotype-beta each edge
          a schematic for Freedman and Lane permuatations is provided in the code comments below
        - summary statistic is computed for each cell (pair of networks) by aggregating the t-values for the edges in the cell, for both observed and permutation t-values
        - permutation based p-values are computed for each cell using the true (observed) summary statistic for each cell
        

    Parameters
    ----------
    edges : ndarray
        2-dimensional numpy array containing connectivity data for each subject. Shape must be (number_of_subjects, number_of_edges)
        Each row is the flattened upper-triangular matrix of pearson connectomes per subject
    phenotype : ndarray
        2-dimensional numpy array containing the phenotypic data for each subject. Shape must be (number_of_subjects, 1)
    covariates : ndarray
        2-dimensional numpy array containing the covariate/nuisance data for each subject. Shape must be (number_of_subjects, number_of_covariates)
    n_perms : int
        Number of permutations that define the null/chance distribution
    t_threshold : float
        T-value threshold used to determine suprathreshold edges. Generally t_threshold=1.96 approximates p=0.05 for large samples
    net_order_file : str
        location of file with network assignments for each edge. Must be single columned csv, with no header
    perm_order : ndarray, optional
        A precomputed order of permutations can be provided. If perm_order is not None, it must be a numpy array of shape [number of subject, n_perms],
        where each row, i, is the permuted order of subjects for the ith observation of null distribution
        If perm_order is not provided, random shuffle orders are generated using np.random with generator-seed=random_seed
    cell_summary_fn : function, optional
        An alternative function can be supplied to generate the per-cell summary statistic used in NCA. 
        inputs of cell_summary_fn must be:
            t_vals_per_edge: numpy array of T-values per edge
            t_threshold: T-value threshold used to determine suprathreshold edges
        note - cell_summary_fn doesn't necessary need to only consider superthreshold edges (eg: NCA can also be run on cell means of t-values per edge)
    random_seed : int, optional
        Random number generator seed for numpy, only used if perm_order is not provided
    
        
    Returns
    -------
    ndarray : 
        1-dim numpy array of shape (number_of_edges, ) with observed t-values for phenotype beta, for each edge
        
    ndarray : 
        2-dim numpy array of shape (n_perms, number_of_edges) with permutation t-values for phenotype beta, for each edge
        
    ndarray :
        1-dim numpy array of shape (number_of_cells, ) with observed cell-wise summary statistics for observed t-values
        
    ndarray :
        2-dim numpy array of shape (n_perms, number_of_cells) with cell-wise summary statistics for permutation t-values

    ndarray :
        1-dim numpy array of shape (number_of_cells, ) with raw p-values for observed cell-wise summary statistics versus null/chance
        
    ndarray :
        1-dim numpy array of shape (number_of_cells, ) with FDR corrected p-values for observed cell-wise summary statistics versus null/chance
    """

    
    assert len(edges.shape) == 2, "edges must be 2 dimensional array with shape (n_subjects, n_edges)"
    assert len(phenotype.shape) == 2, "phenotype must be 2 dimensional array with shape (n_subjects, 1)"
    assert len(covariates.shape) == 2, "covariates must be 2 dimensional array with shape (n_subjects, n_covariates)"
    assert edges.shape[0] == phenotype.shape[0], "number of rows in edges and phenotype arrays must match"
    assert covariates.shape[0] == phenotype.shape[0], "number of rows in edges and covariates arrays must match"
    if perm_order is not None:
        assert perm_order.shape[0] == edges.shape[0], "number of rows in edges and perm_order arrays must match"
        assert perm_order.shape[1] == n_perms, "number of columns in perm_order must match n_perms"
    else:
        np.random.seed(random_seed)
        n_subjects = edges.shape[0]
        perm_order = np.zeros((n_subjects, n_perms))
        for i in range(n_perms):
            perm_order[:, i] = np.random.choice(np.arange(n_subjects), replace=False, size=n_subjects)
	perm_order = perm_order.astype(np.int32)  # needs to be int array to be used for row-indexing later
    
    
    # first step of schematic above: 
    #     per edge linear regression is conducted with the following formula: edge_i ~ phenotype + covariates
    print('running step 1 - per edge linear regression is conducted with the following formula: edge_i ~ phenotype + covariates...')
    X = np.hstack((phenotype,                         # phenotype
                   np.ones((phenotype.shape[0], 1)),  # intercept term 
                   covariates))                       # covariate block
    n_subjects, n_predictors = X.shape
    X_pinv = scipy.linalg.pinv2(X)
    inv_vcv = scipy.linalg.inv(np.dot(X.T, X))               # shape: (n_predictors, n_predictors)
    inv_vcv = np.sqrt(np.diag(inv_vcv))
    betas = np.dot(X_pinv, edges)                            # shape: (n_predictors, n_edges)
    fit = np.dot(X, betas)                                   # shape: (n_subjects, n_edges)
    residuals = edges - fit                                  # shape: (n_subjects, n_edges)
    eps_std = np.std(residuals, axis=0, ddof=n_predictors)   # shape: (n_edges, )
    pheno_beta_std = inv_vcv[0] * eps_std                    # shape: (n_edges, )
    observed_tvalues =  betas[0, :]/pheno_beta_std           # shape: (n_edges, ) array of observed t-values for the phenotype-betas
    
    # second step of schematic above: 
    #     generate null/chance distribution for phenotype betas and t-values using Freedman and Lane permutations
    # note - perm_betas are standardized if edges, phenotypes and covariates are column-wise zscored (where appropriate)
    print('running step 2 - generate null/chance distribution for phenotype betas and t-values using Freedman and Lane permutations...')
    perm_tvalues, perm_betas = freedman_lane_permutations(edges, phenotype, covariates, n_perms, perm_order)
    
    # third step of schematic above:
    #     summary statistics and pvalues computed for each cell, for both observed and permuted tvalues
    print('running steps 3/4 - summary statistic and pvalues are computed for each cell, for both observed and permuted tvalues...')
    cell_edges_dict = generate_cell_edges_dict(net_order_file)
    
    observed_cell_summaries = np.zeros(len(cell_edges_dict))
    perm_cell_summaries = np.zeros((len(cell_edges_dict), n_perms))
    cell_raw_pvals = np.zeros(len(cell_edges_dict))
    cell_order = []
    for cell_i, cell in enumerate(cell_edges_dict):
        cell_edges_idxs = cell_edges_dict[cell]
        observed_cell_summary = cell_summary_fn(observed_tvalues[cell_edges_idxs], t_threshold)
        perms_cell_summary = np.apply_along_axis(lambda x: cell_summary_fn(x, t_threshold),
                                                axis=1,
                                                arr=perm_tvalues[:, cell_edges_idxs])
        observed_cell_summaries[cell_i] = observed_cell_summary
        perm_cell_summaries[cell_i, :] = perms_cell_summary
        cell_p_val = np.sum(perms_cell_summary >= observed_cell_summary) / n_perms  # one-sided hypothesis test
        cell_raw_pvals[cell_i] = cell_p_val
	cell_order.append(cell)
    _, fdr_corrected_cell_pvals = statsmodels.stats.multitest.fdrcorrection(cell_raw_pvals, alpha=0.05)
    
    return observed_tvalues, perm_tvalues, observed_cell_summaries, perm_cell_summaries, cell_raw_pvals, fdr_corrected_cell_pvals, cell_order
   
    
def freedman_lane_permutations(edges, phenotype, covariates, n_perms, perm_order):
    """
    third step of schematic described in network_contingency_analysis: 
         generate Freedman and Lane permutations (1983) to determine the null/chance t-value distribution for the phenotype-beta each edge
         the steps for Freedman and Lane permutation are as follows:
              1. fit edge ~ intercept + covariates, for each edge
              2. retain residuals R and predictions edge-hat
              3. permute R (call it R-tilde) and then make edge-tilde = R-tilde + edge-hat
              4. fit edge-tilde ~ phenotype + covariates, and store the tvalue for phenotype-beta
    
    Parameters
    ----------
    edges : ndarray
        2-dimensional numpy array containing connectivity data for each subject. Shape must be (number_of_subjects, number_of_edges)
        Each row is generally the flattened upper/lower-triangular matrix of pearson connectomes per subject
    phenotype : ndarray
        2-dimensional numpy array containing the phenotypic data for each subject. Shape must be (number_of_subjects, 1)
    covariates : ndarray
        2-dimensional numpy array containing the covariate/nuisance data for each subject. Shape must be (number_of_subjects, number_of_covariates)
    n_perms : int
        Number of permutations that define the null/chance distribution
    perm_order : ndarray, optional
        A precomputed order of permutations can be provided. If perm_order is not None, it must be a numpy array of shape [number of subject, n_perms],
        where each row, i, is the permuted order of subjects for the ith observation of null distribution
        If perm_order is not provided, random shuffle orders are generated using np.random with generator-seed=random_seed

    Returns
    -------
    ndarray 
        Numpy array of shape (number_of_permutations, number_of_edges), containing t-values per edge for phenotype beta, for each permutation iteration
        
    ndarray 
        Numpy array of shape (number_of_permutations, number_of_edges), containing phenotype-betas per edge, for each permutation iteration

    """
    # steps 1, 2 in function header
    X = np.hstack((np.ones((edges.shape[0], 1)).astype(np.float),
                   covariates))
    pinv_X = scipy.linalg.pinv2(X)
    betas = np.dot(pinv_X, edges)
    edges_hat = np.dot(X, betas)
    Rs = edges - edges_hat

    # precompute the pseudo-inv and estimator variance for permutation models
    X = np.hstack((phenotype,                         # phenotype
                   np.ones((phenotype.shape[0], 1)),  # intercept term 
                   covariates))                       # covariate block
    pinv = scipy.linalg.pinv2(X)
    inv_vcv = scipy.linalg.inv(np.dot(X.T, X))
    inv_vcv = np.sqrt(np.diag(inv_vcv))
    n_subjects, n_predictors = X.shape
    
    # generate permutation distribution
    perm_tvalues = np.zeros((n_perms, edges.shape[1]))
    perm_betas = np.zeros((n_perms, edges.shape[1]))
    for perm_idx in range(n_perms):
        Rs_tilde = Rs[perm_order[:, perm_idx], :]   # edges are shuffled row-wise, to maintain across-edge covariance structure (within subject)
        edges_tilde = Rs_tilde + edges_hat
        # compute betas and tvalues for edge-tilde ~ phenotype + covariates, for each edge
        betas = np.dot(pinv, edges_tilde)    
        fit = np.dot(X, betas)
        residuals = edges_tilde - fit
        eps_std = np.std(residuals, axis=0, ddof=n_predictors)
        beta_std = inv_vcv[0] * eps_std
        perm_tvalues[perm_idx, :] = betas[0, :]/beta_std
        perm_betas[perm_idx, :] = betas[0, :]
        
        if (perm_idx+1) % (n_perms // min(n_perms, 10)) == 0:
            print(f'     freedman-lane perms {int(100*((perm_idx+1) / n_perms))}% complete...')
            
    return perm_tvalues, perm_betas
	
	
def generate_cell_edges_dict(net_order_file):
    """ 
    Generates dictionary that maps cell_id -> (indices of edges) for cell
        cell_id is sorted(network_i_id, network_j_id), where network ids are provided in net_order_file
        (indices of edges) are locations in flattened upper-triangular connectome for each cell
        
    Parameters
    ----------
    net_order_file : str
        location of file with network assignments for each edge. Must be single columned csv, with no header
        
    Returns
    -------
    dict : (int/str, int/str) -> ndarray
        key is cell_id, where cell_id is sorted(network_i_id, network_j_id)
        value is indices of edges for cell for network_i and network_j in flattened upper-triangular connectome

    """
    net_order = pd.read_csv(net_order_file, header=None).values.flatten()
    utri_idxs = np.triu_indices(net_order.shape[0], k=1)
    unique_nets = np.unique(net_order)
    cell_edges_dict = dict()
    
    for net_i in unique_nets:
        net_i_idxs = np.argwhere(net_order == net_i).flatten()
        for net_j in unique_nets:
            net_j_idxs = np.argwhere(net_order == net_j).flatten()

            netmask = np.zeros((net_order.shape[0], net_order.shape[0]))
            for x, y in itertools.product(net_i_idxs, net_j_idxs):
                netmask[x, y] = 1
                netmask[y, x] = 1
            for x, y in itertools.product(net_j_idxs, net_i_idxs):
                netmask[y, x] = 1
                netmask[x, y] = 1

            cell_edges_dict[tuple(sorted((net_i, net_j)))] = np.argwhere(netmask[utri_idxs] == 1).flatten()
    
    return cell_edges_dict
