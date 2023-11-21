import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import os

adata = ad.read_h5ad('/mnt/plummergrp/QuPath_0.4.4/BrCa/session_file.h5ad')
adata.obs['batch'] = pd.Categorical(adata.obs['ImageID'])
# adata.obs['batch'] = adata.obs['batch'].cat.rename_categories({'0': 'sample1', '1': ' sample2'})
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, 'batch')

adata.write_h5ad('/mnt/plummergrp/QuPath_0.4.4/BrCa/session_file_harmony.h5ad')
