import numpy as np
import pandas as pd
import scanpy as sc

def manifold_alignment(Wx, Wy, Wxy = None, mu = 0.9, d = 30, normlized = True, eps = 1e-5):
    Wx = np.array(Wx)
    Wy = np.array(Wy)
    n = Wx.shape[0]

    if normlized == True:
        Wx = Wx/np.max(np.abs(Wx))
        Wy = Wy/np.max(np.abs(Wy))
        # Wx = Wx + 1
        # Wy = Wy + 1

    if Wxy == None:
        Wxy = np.diag(np.repeat(1, n))

    Wxy = mu * (np.sum(Wx) + np.sum(Wy)) / (2 * np.sum(Wxy)) * Wxy
    W = np.asarray(np.bmat(((Wx, Wxy), (Wxy.T, Wy))))
    W = (W + W.T)/2
    D = np.diag(np.sum(W, axis = 0))
    L = D - W

    vals, vecs = np.linalg.eig(L)
    idx = np.argsort(vals)
    for i in range(len(idx)):
        if vals[idx[i]] >= eps:
            break

    alignedNet = vecs.real[:,idx[i:(i + d)]]
    for i in range(alignedNet.shape[1]):
        alignedNet[:,i] /= np.linalg.norm(alignedNet[:,i])

    return alignedNet

def data_process(dta, min_genes = 100, min_cells = 10, mt_pct = 10, npcs = 50, oversd = None):
    adata = sc.AnnData(dta)
    adata.var_names_make_unique()
    # quality control
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_genes(adata, min_cells = min_cells)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)

    if oversd is not None:
        mu = np.mean(adata.obs.n_genes_by_counts)
        sd = np.std(adata.obs.n_genes_by_counts)
        thres = mu + oversd * sd
        adata = adata[adata.obs.n_genes_by_counts < thres, :]
        
    adata = adata[adata.obs.pct_counts_mt < mt_pct, :]
    # normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # find highly variable gene
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    # pca
    sc.tl.pca(adata, svd_solver='arpack', n_comps = npcs)
    
    return adata