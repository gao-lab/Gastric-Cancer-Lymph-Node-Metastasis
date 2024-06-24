#!/usr/bin/env python
# coding: utf-8

# # blast model

# In[26]:


import time
import warnings
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import Cell_BLAST as cb
print(cb.__version__)

warnings.filterwarnings("ignore")
cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0


# ## read data and split train and test data

# In[2]:


cd8_tcell = anndata.read_h5ad("cd8_tcell_raw.h5ad")


# In[3]:


cd8_tcell.obs.head()


# In[4]:


sc.pl.umap(cd8_tcell, color="subtype1", palette="tab20")


# In[5]:


sc.pl.umap(cd8_tcell, color="Fresh", palette="tab10")


# In[6]:


np.unique(cd8_tcell.obs["subtype1"])


# In[10]:


cd8_tcell.obs["subtype0"]=cd8_tcell.obs["subtype1"].astype("object")
cd8_tcell.obs["subtype0"][np.isin(cd8_tcell.obs["subtype1"], ['CD8_CCL5+ Tex', 'CD8_LAYN+ Tex','CD8_MKI67+ Tex'])]="CD8_Tex"
cd8_tcell.obs["subtype0"][np.isin(cd8_tcell.obs["subtype1"], ['CD8_HLA-DRA+ Tem', 'CD8_HSPA1B+ Tem',
        'CD8_Tcm', 'CD8_Tem', 'CD8_Trm'])]="CD8_Tm"
cd8_tcell.obs["subtype0"][np.isin(cd8_tcell.obs["subtype1"], ['CD8_TIM3+ Trm', 'CD8_TIM3+MKI67+ Trm'])]="CD8_TIM3+ Trm"


# In[11]:


np.unique(cd8_tcell.obs["subtype0"])


# ## find_variable_genes

# In[12]:


get_ipython().run_cell_magic('capture', '', 'axes = cb.data.find_variable_genes(cd8_tcell, grouping="Fresh")')


# In[13]:


cd8_tcell.var["variable_genes"].sum()


# In[14]:


genes_excluded_all = pd.read_csv("../genes_excluded_all.csv",index_col=0)["x"].values
genes_excluded_all.shape


# In[15]:


# number of HVG in excluded 
np.isin(cd8_tcell.var["variable_genes"].index[cd8_tcell.var["variable_genes"]].values,
        genes_excluded_all).sum()


# In[16]:


# modify HVG
cd8_tcell.var["variable_genes"][np.isin(cd8_tcell.var["variable_genes"].index,genes_excluded_all)]=False


# In[17]:


cd8_tcell.var["variable_genes"].sum()


# In[18]:


cd8_tcell.var["variable_genes"].index[cd8_tcell.var["variable_genes"]].values


# ## one model for visualization (batch correction(fresh/frozen))

# In[19]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodel_rmbatch = cb.directi.fit_DIRECTi(\n    cd8_tcell, genes=cd8_tcell.var.query("variable_genes").index.to_numpy(),\n    batch_effect="Fresh", latent_dim=10, cat_dim=20\n)\ncd8_tcell.obsm["X_latent"] = model_rmbatch.inference(cd8_tcell)')


# In[20]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[21]:


sc.pp.neighbors(cd8_tcell, use_rep="X_latent")
sc.tl.umap(cd8_tcell)


# In[22]:


sc.pl.umap(cd8_tcell, color="subtype0", palette="tab10")


# In[23]:


sc.pl.umap(cd8_tcell, color="subtype1", palette="tab20b")


# In[24]:


sc.pl.umap(cd8_tcell, color="Fresh", palette="tab10")


# In[25]:


model_rmbatch.save("./cd8_tcell_rmbatch_model_new")


# In[ ]:





# ## multi models and blast model

# In[27]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodels = []\nfor i in range(4):\n    models.append(cb.directi.fit_DIRECTi(\n        cd8_tcell, genes=cd8_tcell.var.query("variable_genes").index.to_numpy(),\n        batch_effect="Fresh", latent_dim=10, cat_dim=20, random_seed=i\n    ))\n# cd8_tcell.obsm["X_latent"] = model_rmbatch.inference(cd8_tcell)')


# In[28]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[29]:


blast = cb.blast.BLAST(models, cd8_tcell)


# In[30]:


blast.save("./cd8_tcell_blast_subtype0")


# In[31]:


blast.save("./cd8_tcell_blast_new")

