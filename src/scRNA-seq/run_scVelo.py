#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 整合loom文件和meta-data数据
import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import os
import scanpy as sc

import matplotlib.pyplot as plt


# In[2]:


scv.set_figure_params()


# In[15]:


# os.chdir("/home/weil/GC")
# os.getcwd()


# In[10]:


loom_data = scv.read('02_cr_results_raw/LIB_for_scVelo/B2-19_combined.loom', cache=False)
loom_data.obs


# In[11]:


# barcode名字去重后缀x，使其与seurat导出的barcode名称一致
loom_data.obs = loom_data.obs.rename(index = lambda x: 
                                     x.replace(':', '_').replace('x', '').replace('-5LIB', ''))
loom_data.obs = loom_data.obs.rename(index = lambda x: 
                                     x.replace('G14', 'B14').replace('16', 'B16').replace('B13L3-2', 'B13L3'))
loom_data.obs.index += "-1"
loom_data.obs.head()


# ## run scVelo for CD8+ T cell

# ### subset CD8+ T cells

# In[22]:


# 读取seurat中的meta信息
meta=pd.read_csv("processed_data/data_B2-19/metadata_all_ct_subtype_condition_noTLB.csv", 
                 header=0, index_col=0)
meta = meta[meta['cell_type1'] == 'T-CD8']
meta['subtype1'].value_counts()

cell_umap= pd.read_csv("processed_data/data_B2-19/cell_typing/tcell/cd8_tcell_cell_embeddings.csv", header=0, 
                       names=["Cell ID", "UMAP_1", "UMAP_2"])


# In[24]:


# 对细胞文件和RNA剪切速率文件取交集，保留关注的细胞类型
sample_cd8 = loom_data[np.isin(loom_data.obs.index, meta.index)]
sample_cd8.obs.head()


# In[25]:


sample_cd8.obs.shape


# In[26]:


cell_umap.head()


# In[27]:


meta_ordered=meta.loc[sample_cd8.obs.index.values,:]
cell_umap.set_index('Cell ID', inplace=True)
cell_umap_ordered=cell_umap.loc[sample_cd8.obs.index.values,:]


# In[28]:


cell_umap_ordered.shape


# In[29]:


np.all(sample_cd8.obs.index==meta_ordered.index)


# In[30]:


np.all(sample_cd8.obs.index==cell_umap_ordered.index)


# In[31]:


sample_cd8.obs=pd.concat([sample_cd8.obs, meta_ordered], axis=1)
sample_cd8.obsm['X_umap'] = cell_umap_ordered.values


# In[32]:


len(sample_cd8.obs["sample"].unique())


# In[33]:


sample_cd8.obs["subtype1"].unique()


# In[34]:


adata_cd8 = sample_cd8
# some gene labels are duplicated (Ensembl IDs are still unique!!)
adata_cd8.var_names_make_unique()

# save model to file
adata_cd8.write('processed_data/data_B2-19/RNA_velocity/cd8_tcell_dynamicModel.h5ad', compression = 'gzip')


# ### run scvelo for HLA-Tex

# In[35]:


# read h5ad file
adata_Tex= scv.read('processed_data/data_B2-19/RNA_velocity/cd8_tcell_dynamicModel.h5ad')


# In[36]:


adata_Tex = adata_Tex[np.isin(adata_Tex.obs["subtype1"], ["CD8_LAYN+ Tex","CD8_CCL5+ Tex","CD8_HLA-DRA+ Tem"]),:]


# In[37]:


adata_Tex


# In[38]:


# Running RNA Velocity
scv.pp.filter_and_normalize(adata_Tex,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata_Tex, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_Tex, mode = "stochastic")
scv.tl.velocity_graph(adata_Tex)


# In[46]:


plt.rcParams['figure.figsize'] = [6, 6]
scv.pl.velocity_embedding_stream(adata_Tex, basis='umap',color = "subtype1",title=None,dpi=300,
                                 save="plots/figures/figS5E_velocity_stream.png")
# figure cannot be saved as pdf, using png instead.
# saving figure to file plots/figures/figS5E_velocity_stream.png


# In[47]:


plt.rcParams['figure.figsize'] = [8,7]
scv.pl.velocity_embedding_grid(adata_Tex, basis='umap',
                               color = "subtype1",legend_loc='right margin',arrow_length=2,
                               title=None,dpi=300,
                                 save="plots/figures/figS5E_velocity_grid.png")


# ### run scvelo for all CD8 cells

# In[12]:


# read h5ad file
adata_cd8= scv.read('processed_data/data_B2-19/RNA_velocity/cd8_tcell_dynamicModel.h5ad')


# In[60]:


np.unique(adata_cd8.obs["subtype1"])


# In[61]:


adata_cd8


# In[62]:


# Running RNA Velocity
scv.pp.filter_and_normalize(adata_cd8,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata_cd8, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_cd8, mode = "stochastic")
scv.tl.velocity_graph(adata_cd8)


# In[69]:


plt.rcParams['figure.figsize'] = [12, 9]
scv.pl.velocity_embedding_stream(adata_cd8, basis='umap',
                               color = "subtype1")#,legend_loc='right margin')


# In[71]:


plt.rcParams['figure.figsize'] = [12, 9]
scv.pl.velocity_embedding_grid(adata_cd8, basis='umap',
                               color = "subtype1",legend_loc='right margin',arrow_length=3)# TODO: check Tgd, is it NK?


# In[ ]:





# ## integrating data (run only once)

# In[3]:


# files=os.listdir("02_cr_results_raw/LIB_for_scVelo")
# files.remove("B2-19_combined.loom")
# files.remove("B2-17_combined.loom")
# files.remove('4T-5LIB.loom')
# files.remove('B3T-5LIB.loom')
# files.remove('B8T-5LIB.loom')
# len(files)


# In[5]:


# os.chdir("02_cr_results_raw/LIB_for_scVelo")
# output_filename='B2-19_combined.loom'
# loompy.combine(files, output_filename, key="Accession")

