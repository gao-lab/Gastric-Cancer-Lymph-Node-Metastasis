#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import networkx as nx
import scipy.io
import anndata
import scanpy as sc
import os
from networkx.algorithms.bipartite import biadjacency_matrix

import scglue


# In[2]:


# os.chdir("ArchR/")
os.getcwd()


# # collect data

# ## scRNA-seq

# In[81]:


rna_matrix = scipy.io.mmread("seRNA/expr_counts.mtx").T.tocsr()
rna_obs = pd.read_csv("seRNA/metadata.csv", header=0, index_col=0)
rna_var = pd.read_csv("seRNA/genes.csv", index_col=0)
rna_obs.index.name, rna_var.index.name = "cells", "genes"
rna = anndata.AnnData(X=rna_matrix, obs=rna_obs, var=rna_var)
rna


# In[82]:


rna.obs["domain"] = "scRNA-seq"
rna.obs["dataset"] = "GC_RNA"


# In[83]:


scglue.data.get_gene_annotation(
    rna, gtf="/home/weil/GENOME/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
    gtf_by="gene_name"
)
rna.var["genome"] = "hg38"


# In[84]:


rna


# In[85]:


rna.obs.head()


# ## scATAC-seq

# In[3]:


atac_matrix = scipy.io.mmread("archr_results/peak_matrix/peak_matrix.mtx").T.tocsr()
atac_obs = pd.read_csv("archr_results/peak_matrix/cellname.txt", header=0, index_col=0)
atac_var = pd.read_csv("archr_results/peak_matrix/peakinfo.txt", header=None, index_col=0)

atac_obs.index.name, atac_var.index.name = "cells", "peaks"
atac = anndata.AnnData(X=atac_matrix, obs=atac_obs, var=atac_var)
atac


# In[4]:


atac.var.head()


# In[5]:


atac.obs.head()


# In[6]:


atac.obs["domain"] = "scATAC-seq"
# atac.obs["protocol"] = "SNARE-seq"
atac.obs["dataset"] = "GC_ATAC"


# In[7]:


atac.var["chrom"] = np.vectorize(lambda x: x.split("_")[0])(atac.var_names)
atac.var["chromStart"] = np.vectorize(lambda x: x.split("_")[1])(atac.var_names).astype(int)
atac.var["chromEnd"] = np.vectorize(lambda x: x.split("_")[2])(atac.var_names).astype(int)
atac.var["genome"] = "hg38"


# ## clean data

# In[86]:


retained_genes = rna.var.dropna(subset=["chrom", "chromStart", "chromEnd"]).index
rna = rna[:, retained_genes]
rna.var = rna.var.astype({"chromStart": int, "chromEnd": int})
rna


# In[87]:


sc.pp.filter_genes(rna, min_counts=1)
rna


# In[8]:


# ArchR has removed peaks from blacklist
atac.var = atac.var.astype({"chromStart": int, "chromEnd": int})
atac


# In[9]:


sc.pp.filter_genes(atac, min_counts=1)
atac


# In[ ]:


# rna.write("GLUE/clean_RNA.h5ad", compression="gzip")
# atac.write("GLUE/clean_ATAC.h5ad", compression="gzip")


# # preprocessing

# In[10]:


import anndata
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams


# In[11]:


scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# In[88]:


# rna = anndata.read_h5ad("GLUE/GC_B2_RNA.h5ad")
rna


# In[12]:


# atac = anndata.read_h5ad("GLUE/GC_B2_ATAC.h5ad")
atac


# In[16]:


# columns_to_drop = ['cell_type1_glue', 'subtype1_glue', 'LN_condition', 'subtype_tex', 'subtype_plasma']

# atac.obs=atac.obs.drop(columns=columns_to_drop)


# In[17]:


# atac = anndata.read_h5ad("GLUE/GC_B2_ATAC.h5ad")
atac


# ## preprocess RNA data

# In[89]:


rna.X, rna.X.data


# In[90]:


rna.layers["counts"] = rna.X.copy()


# In[91]:


from skmisc.loess import loess


# In[92]:


sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")


# In[93]:


genes_excluded = pd.read_csv("../../GENOME/scrna/genes_excluded.csv",index_col=0)["x"]


# In[94]:


rna


# In[95]:


rna.var["gene"]=rna.var.index
rna.var["gene"].isin(genes_excluded).value_counts()


# In[96]:


(rna.var["gene"].isin(genes_excluded) & rna.var["highly_variable"]).value_counts()


# In[97]:


rna.var.loc[rna.var["gene"].isin(genes_excluded), "highly_variable"] = False


# In[98]:


(rna.var["gene"].isin(genes_excluded) & rna.var["highly_variable"]).value_counts()


# In[99]:


rna


# In[ ]:





# In[100]:


sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")


# In[101]:


sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.umap(rna, color="cell_type")


# In[102]:


sc.pl.umap(rna, color="sample")


# In[103]:


sc.pl.umap(rna, color="subtype")


# In[106]:


sc.pl.umap(rna, color="patient")


# In[107]:


sc.pl.umap(rna, color="Fresh")


# In[108]:


sc.pl.umap(rna, color="cell_type1")


# In[ ]:


rna


# In[ ]:


# rna = anndata.read_h5ad("GLUE/rna_preprocessed.h5ad")


# ## re analysis scATAC-seq data with harmony

# In[35]:


# atac = anndata.read_h5ad("GLUE/atac_glue.h5ad")
# atac


# In[20]:


atac


# In[19]:


atac.X, atac.X.data


# In[21]:


scglue.data.lsi(atac, n_components=100, n_iter=15)


# In[22]:


atac


# In[23]:


import scanpy.external as sce


# In[24]:


sce.pp.harmony_integrate(atac, "patient", basis='X_lsi')


# In[32]:


sc.pp.neighbors(atac, use_rep="X_pca_harmony", metric="cosine")
sc.tl.umap(atac)


# In[34]:


atac


# In[36]:


sc.pl.umap(atac, color=["cell_type", 
                        "predictedGroup_Un"], wspace=0.65,ncols=2)


# In[37]:


sc.pl.umap(atac, color=["Sample", 
                        "patient"], wspace=0.65,ncols=2)


# In[38]:


sc.pl.umap(atac, color="region")


# In[ ]:





# ## Construct prior regulatory graph

# In[109]:


# atac=sc.read("GLUE/atac_preprocessed.h5ad")
# atac


# In[110]:


# atac.var = atac.var.drop(["highly_variable"], axis=1)


# In[111]:


atac


# ### Obtain genomic coordinates

# In[112]:


# done in collect data section


# ### Graph construction

# In[113]:


graph = scglue.genomics.rna_anchored_prior_graph(rna, atac)


# In[114]:


graph


# In[115]:


graph.number_of_nodes(), graph.number_of_edges()


# In[116]:


# Graph node covers all omic features
all(graph.has_node(gene) for gene in rna.var_names), all(graph.has_node(peak) for peak in atac.var_names)


# In[117]:


# Edge attributes contain weights and signs
for _, e in zip(range(5), graph.edges):
    print(f"{e}: {graph.edges[e]}")


# In[118]:


# Each node has a self-loop
all(graph.has_edge(gene, gene) for gene in rna.var_names), all(graph.has_edge(peak, peak) for peak in atac.var_names)


# In[119]:


# Graph is symmetric
all(graph.has_edge(j, i) for i, j, _ in graph.edges)


# In[120]:


atac.var.head()


# In[ ]:





# ## Save preprocessed data files

# In[121]:


rna.write("GLUE/rna_preprocessed.h5ad", compression="gzip")
atac.write("GLUE/atac_preprocessed.h5ad", compression="gzip")
nx.write_graphml(graph, "GLUE/prior.graphml.gz")


# In[122]:


# rna = anndata.read_h5ad("GLUE/rna_preprocessed.h5ad")
# atac = anndata.read_h5ad("GLUE/atac_preprocessed.h5ad")


# In[123]:


atac


# In[124]:


rna


# In[ ]:




