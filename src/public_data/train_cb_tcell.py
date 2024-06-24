#!/usr/bin/env python
# coding: utf-8

# # blast model

# In[1]:


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


# In[2]:


import scipy
print(scipy.__version__)


# In[3]:


import numba
print(numba.__version__)


# ## read data and split train and test data

# In[4]:


tcell = anndata.read_h5ad("tcell_raw.h5ad")


# In[3]:


tcell.obs.head()


# In[4]:


sc.pl.umap(tcell, color="subtype1", palette="tab20")


# In[5]:


sc.pl.umap(tcell, color="Fresh", palette="tab10")


# In[6]:


np.unique(tcell.obs["subtype1"])


# In[ ]:





# In[7]:


# tcell.obs["subtype0"]=tcell.obs["subtype1"].astype("object")
# tcell.obs["subtype0"][np.isin(tcell.obs["subtype1"], ['CD4_CCL5+ Tex', 'CD4_LAYN+ Tex','CD4_MKI67+ Tex'])]="CD4_Tex"
# tcell.obs["subtype0"][np.isin(tcell.obs["subtype1"], ['CD4_HLA-DRA+ Tem', 'CD4_HSPA1B+ Tem',
#         'CD4_Tcm', 'CD4_Tem', 'CD4_Trm'])]="CD4_Tm"
# tcell.obs["subtype0"][np.isin(tcell.obs["subtype1"], ['CD4_TIM3+ Trm', 'CD4_TIM3+MKI67+ Trm'])]="CD4_TIM3+ Trm"


# In[8]:


# np.unique(tcell.obs["subtype0"])


# ## find_variable_genes

# In[9]:


get_ipython().run_cell_magic('capture', '', 'axes = cb.data.find_variable_genes(tcell, grouping="Fresh")')


# In[10]:


tcell.var["variable_genes"].sum()


# In[11]:


genes_excluded_all = pd.read_csv("../genes_excluded_all.csv",index_col=0)["x"].values
genes_excluded_all.shape


# In[12]:


# number of HVG in excluded 
np.isin(tcell.var["variable_genes"].index[tcell.var["variable_genes"]].values,
        genes_excluded_all).sum()


# In[13]:


# modify HVG
tcell.var["variable_genes"][np.isin(tcell.var["variable_genes"].index,genes_excluded_all)]=False


# In[14]:


tcell.var["variable_genes"].sum()


# In[15]:


tcell.var["variable_genes"].index[tcell.var["variable_genes"]].values


# ## one model for visualization (batch correction(fresh/frozen))

# In[16]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodel_rmbatch = cb.directi.fit_DIRECTi(\n    tcell, genes=tcell.var.query("variable_genes").index.to_numpy(),\n    batch_effect="Fresh", latent_dim=10, cat_dim=20\n)\ntcell.obsm["X_latent"] = model_rmbatch.inference(tcell)')


# In[17]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[18]:


sc.pp.neighbors(tcell, use_rep="X_latent")
sc.tl.umap(tcell)


# In[19]:


sc.pl.umap(tcell, color="subtype1", palette="tab20b")


# In[20]:


sc.pl.umap(tcell, color="Fresh", palette="tab10")


# In[21]:


model_rmbatch.save("tcell_one_test/tcell_rmbatch_model_new")


# ### visualization

# In[5]:


model_rmbatch = cb.directi.DIRECTi.load("./tcell_one_test/tcell_rmbatch_model_new")


# In[6]:


tcell.obsm["X_latent"] = model_rmbatch.inference(tcell)


# In[7]:


sc.pp.neighbors(tcell, use_rep="X_latent")
sc.tl.umap(tcell)


# In[8]:


sc.pl.umap(tcell, color="subtype1", palette="tab20b")


# In[ ]:





# ## multi models and blast model

# In[24]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodels = []\nfor i in range(4):\n    models.append(cb.directi.fit_DIRECTi(\n        tcell, genes=tcell.var.query("variable_genes").index.to_numpy(),\n        batch_effect="Fresh", latent_dim=10, cat_dim=20, random_seed=i\n    ))\n# tcell.obsm["X_latent"] = model_rmbatch.inference(tcell)')


# In[25]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[26]:


blast = cb.blast.BLAST(models, tcell)


# In[ ]:


blast.save("./tcell_blast")


# In[ ]:





# In[ ]:





# In[ ]:




