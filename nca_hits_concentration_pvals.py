import numpy as np
import pandas as pd

# roi to network assigned for our Gordon (333) + subcortical (65) + cerebellar (66) atlas
# total number of ROIs = 418
roi_to_network_mapping = [ 4,  0,  1,  4,  6,  4,  7,  6,  7,  3, 13, 14, 15, 15,  6,  6,  6,
                          13, 13,  6,  2,  2, 10,  7,  4,  4,  2,  2,  8,  0,  0,  0,  0,  2,
                          0,  0,  0,  0,  1,  2, 11, 11, 11,  4,  0,  0,  0,  0, 11,  0, 11,
                          11,  1,  0, 11,  0,  0,  0,  1, 10, 10, 10,  2,  3,  3,  3,  3,  3,
                          3,  3,  2,  2, 13, 11, 10,  2,  3,  7, 10, 10,  2,  2,  8,  2, 10,
                          10, 11, 11, 14,  6, 11, 11, 14,  4, 11,  7,  6,  6,  6, 11,  2,  3,
                          2,  3,  2, 11, 11,  7,  7, 11,  2,  2, 11,  4, 13,  4,  4, 13, 13,
                          13, 13, 13, 13, 13, 13,  4,  4, 13, 13, 15,  6,  6, 13, 13, 13,  6,
                          6,  6,  6,  6,  6, 13, 15, 13,  4,  4,  2,  7,  7,  4,  4,  4,  2,
                          4, 11,  4,  4, 10, 13,  3, 10,  4,  0,  1,  4,  6,  7,  7,  6,  7,
                          3, 13, 14, 15,  6,  6,  6, 13, 13,  2,  2,  7,  8,  4,  2,  4,  2,
                          2, 11,  0,  0,  2,  0,  0,  0,  2,  1,  2, 11,  4,  0,  0, 11,  0,
                          0,  0,  0, 11,  0,  0, 11,  1,  0,  0,  0,  0,  0,  1,  2,  4, 10,
                          10,  2,  3,  4, 10,  3, 10, 10,  3, 10,  3,  3,  2,  2, 11, 10,  2,
                          3,  7, 10, 10, 10,  3,  2,  2,  8,  2,  2, 11,  6, 11, 11, 14,  6,
                          6,  4,  6,  4,  7,  7, 11,  6,  6,  6, 11,  6,  3,  3,  0, 11,  7,
                          7,  2, 11,  7,  7,  4,  4, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                          4, 13, 13,  6, 15, 15, 13, 13,  6,  6, 13, 13, 13, 13, 13, 13, 13,
                          6,  6,  6,  6,  6, 13, 15, 13,  4,  4,  2,  2,  7,  7,  4,  4,  4,
                          4,  4,  4,  7,  7,  3,  3,  4, 10, 10,  9,  9,  9,  9,  9,  9,  9,
                          9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                          9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                          9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 12, 12, 12, 12,
                          12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                          12, 12, 12, 12, 12, 12, 12, 12, 12, 12]

# network names corresponding to the indices above
net_names = ['SMhand', 'SMmouth', 'CinguloOperc', 'Auditory', 'Default', 'N/A', 'Visual', 
             'FrontoParietal', 'Salience', 'Subcortical', 'VentralAttn', 'DorsalAttn', 'Cerebellum', 
             'None', 'CinguloParietal', 'RetrosplenialTemporal']
net_names_no_NA = [x for x in net_names if x != 'N/A']

# functional clusters used for concentration analysis
# somatosensory/subcortical 
group_a = ['SMhand', 'SMmouth', 'Cerebellum', 'Auditory', 'Visual', 'Subcortical']
# default/control
group_b = ['Default', 'FrontoParietal', 'DorsalAttn', 'VentralAttn', 'CinguloOperc']
# other
group_c = ['Salience', 'None', 'CinguloParietal', 'RetrosplenialTemporal']


def nca_hits_concentration_pvals(fdr_corrected_cell_pvals, cell_order, n_perms=10000):
    """
    Concentration analysis to assess whether cells significantly associated with outcome variables of interest (based on the previously described NCA analysis) were concentrated within three functional clusters    Schematically, the procedure is as follows: 
        - use function network_contingency_analysis to generate cell-wise NCA results (fdr_corrected_cell_pvals, cell_order are part of the outputs of network_contingency_analysis)
        - we calculate the concentration of significant cells within each cluster (our statistic of interest), defined as the average number of NCA significant cells for a functional cluster normalized by the total number of NCA significant cells across all three clusters
        - the non-parametric permutation null distribution of this statistic of interest is constructed by repeatedly shuffling cell labels across the entire connectome and recomputing the statistic of interest
        - The p-value of the observed concentration of NCA results per cluster is calculated as the proportion of the permutation distribution greater than observed
    
    Parameters
    ----------
    fdr_corrected_cell_pvals : ndarray
        1-dim numpy array of shape (number_of_cells, ) with FDR corrected p-values for observed cell-wise summary statistics versus null/chance
    cell_order : ndarray
    	1-dim numpy array of shape (number_of_cells, ) where each element is a tuple that uniquely identifies each cell (cell id tuple corresponds to network indices that make up the cell)
    n_perms: integer
        Number of permutations that define the null/chance distribution
        
    Returns
    -------
    float:
        permutation based p-value of concentration statistic for group A (somatosensory/subcortical)
        
    float:
        permutation based p-value of concentration statistic for group B (default/control)
        
    float:
        permutation based p-value of concentration statistic for group C (other)
    """

    
    cell_p_val_corrected_dict = dict()
    for cell, pval in zip(cell_order, fdr_corrected_cell_pvals):
        cell_p_val_corrected_dict[cell] = pval
        
    unique_nets = np.unique(roi_to_network_mapping)
    res_mat = np.zeros((len(unique_nets), len(unique_nets)))
    res_utri = utri_idxs = np.triu_indices(len(unique_nets), k=0)
    for cell in cell_p_val_corrected_dict:
        if 5 in cell: 
            raise Exception('should not happen')
        x = cell[0] if cell[0] < 5 else cell[0] - 1
        y = cell[1] if cell[1] < 5 else cell[1] - 1
        if cell_p_val_corrected_dict[cell] < 0.05:
            res_mat[x, y] = 1
            res_mat[y, x] = 1

    assert np.all(res_mat == res_mat.T)
    hits_vec = res_mat[res_utri].flatten()
    true_cell_sums = res_mat.sum(0)
    
    print(f'running {n_perms} permutations...')
    perm_cell_sums = np.zeros((n_perms, len(unique_nets)))
    for perm_i in range(n_perms):
        random_res_mat = np.zeros((len(unique_nets), len(unique_nets)))
        random_hits_vec = np.random.choice(hits_vec, size=len(hits_vec), replace=False)
        random_res_mat[res_utri] = random_hits_vec
        random_res_mat += random_res_mat.T
        random_cell_sums = random_res_mat.sum(0)
        perm_cell_sums[perm_i, :] = random_cell_sums

    df = pd.DataFrame()
    df['network'] = net_names_no_NA
    df['hits_count'] = true_cell_sums
    true_group_counts = []
    for group in [group_a, group_b, group_c]:
        group_df = df[df.network.isin(group)]
        true_group_counts.append(group_df.hits_count.values.mean())
    true_group_counts = np.array(true_group_counts) / np.sum(true_group_counts)

    perm_group_means = np.zeros((n_perms, 3))
    df = pd.DataFrame()
    df['network'] = net_names_no_NA
    for perm_i in range(n_perms):
        df['hits_count'] = perm_cell_sums[perm_i, :]
        perm_group_counts = []
        for group in [group_a, group_b, group_c]:
            group_df = df[df.network.isin(group)]
            perm_group_counts.append(group_df.hits_count.values.mean())

        perm_group_counts = np.array(perm_group_counts) / np.sum(perm_group_counts)
        perm_group_means[perm_i, :] = perm_group_counts
    
    group_1_pval = np.sum(perm_group_means[:, 0] > true_group_counts[0]) / n_perms
    group_2_pval = np.sum(perm_group_means[:, 1] > true_group_counts[1]) / n_perms
    group_3_pval = np.sum(perm_group_means[:, 2] > true_group_counts[2]) / n_perms
    return group_1_pval, group_2_pval, group_3_pval