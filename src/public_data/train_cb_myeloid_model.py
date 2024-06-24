#!/usr/bin/env python
# coding: utf-8

# In[1]:


import time
import warnings
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import Cell_BLAST as cb
print(cb.__version__)

warnings.simplefilter("ignore")
cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0


# ## read data

# In[2]:


myeloid = anndata.read_h5ad("./myeloid_raw.h5ad")


# In[3]:


sc.pl.umap(myeloid, color="subtype1", palette="tab20")


# In[4]:


sc.pl.umap(myeloid, color="batch", palette="tab10")


# In[5]:


myeloid.obs.head()


# ## find_variable_genes

# In[6]:


get_ipython().run_cell_magic('capture', '', 'axes = cb.data.find_variable_genes(myeloid, grouping="batch")')


# In[7]:


myeloid.var["variable_genes"].sum()


# In[8]:


genes_excluded_all = pd.read_csv("../genes_excluded_all.csv",index_col=0)["x"].values
genes_excluded_all.shape


# In[9]:


# number of HVG in excluded 
np.isin(myeloid.var["variable_genes"].index[myeloid.var["variable_genes"]].values,
        genes_excluded_all).sum()


# In[10]:


# modify HVG
myeloid.var["variable_genes"][np.isin(myeloid.var["variable_genes"].index,genes_excluded_all)]=False


# In[11]:


myeloid.var["variable_genes"].sum()


# ## one model for visualization (batch correction(fresh/frozen))

# In[12]:


time.time()


# In[13]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodel_rmbatch = cb.directi.fit_DIRECTi(\n    myeloid, genes=myeloid.var.query("variable_genes").index.to_numpy(),\n    batch_effect="batch", latent_dim=10, cat_dim=20\n)\nmyeloid.obsm["X_latent"] = model_rmbatch.inference(myeloid)')


# In[14]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[15]:


sc.pp.neighbors(myeloid, use_rep="X_latent")
sc.tl.umap(myeloid)


# In[16]:


sc.pl.umap(myeloid, color="subtype1", palette="tab20")


# In[17]:


sc.pl.umap(myeloid, color="subtype1", palette="tab20b")


# In[18]:


sc.pl.umap(myeloid, color="cell_type", palette="tab20c")


# In[19]:


sc.pl.umap(myeloid, color="batch", palette="tab10")


# In[20]:


model_rmbatch.save("./myeloid_rmbatch_model")


# ## multi models and blast model

# In[21]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodels = []\nfor i in range(4):\n    models.append(cb.directi.fit_DIRECTi(\n        myeloid, genes=myeloid.var.query("variable_genes").index.to_numpy(),\n        batch_effect="batch", latent_dim=10, cat_dim=20, random_seed=i\n    ))')


# In[22]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[23]:


blast = cb.blast.BLAST(models, myeloid)


# In[24]:


blast.save("./myeloid_blast")


# In[ ]:




