import numpy as np
import anndata as ad
import scanpy as sc
import scimap as sm
import pandas as pd
import sys

path_adata = str(sys.argv[1])
path_workflow = str(sys.argv[2])
path_output = str(sys.argv[3])

# for each marker clip value to mean of top 20 values
def mediantop20(subAD):
    outAD = subAD.copy()
    for ix,x in enumerate(outAD.var_names):
        aX = subAD[:,x].X.flatten()
        top20 = np.sort(aX)[-20:]
        outAD.X[:,ix] = np.clip(subAD[:,x].X.flatten(),0,np.median(top20))
    return outAD

# remove expression outliers from the data
def removeOutliers(ad):

    # separate each sample
    s = {}
    for sID in ad.obs.ImageID.sort_values().unique():
        print(sID)
        s[sID] = ad[ad.obs['ImageID']==sID]
        print(s[sID].X.max(axis=0))
        s[sID] = mediantop20(s[sID])
        print(s[sID].X.max(axis=0))
    outAD = sc.concat(s.values())
    return outAD

# rescale and run cell phenotyping
adata = ad.read_h5ad(path_adata)

newAD = removeOutliers(adata)
newAD = sm.pp.rescale(newAD,imageid='ImageID', method='by_image')

# phenotype the cells
phenoDF = pd.read_csv(path_workflow)
# phenoDF.drop(inplace=True, columns=malGenes)
sm.tl.phenotype_cells(newAD,phenoDF,label='phenotype')
newAD.obs.groupby('ImageID').phenotype.value_counts()
newAD.write_h5ad(path_output)