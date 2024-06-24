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


# In[4]:


import anndata
print(anndata.__version__)


# In[5]:



print(sc.__version__)


# ## read data and split train and test data

# In[6]:


cd4_tcell = anndata.read_h5ad("cd4_tcell_raw.h5ad")


# In[3]:


cd4_tcell.obs.head()


# In[4]:


sc.pl.umap(cd4_tcell, color="subtype1", palette="tab20")


# In[5]:


sc.pl.umap(cd4_tcell, color="Fresh", palette="tab10")


# In[6]:


np.unique(cd4_tcell.obs["subtype1"])


# In[ ]:





# In[7]:


# cd4_tcell.obs["subtype0"]=cd4_tcell.obs["subtype1"].astype("object")
# cd4_tcell.obs["subtype0"][np.isin(cd4_tcell.obs["subtype1"], ['CD4_CCL5+ Tex', 'CD4_LAYN+ Tex','CD4_MKI67+ Tex'])]="CD4_Tex"
# cd4_tcell.obs["subtype0"][np.isin(cd4_tcell.obs["subtype1"], ['CD4_HLA-DRA+ Tem', 'CD4_HSPA1B+ Tem',
#         'CD4_Tcm', 'CD4_Tem', 'CD4_Trm'])]="CD4_Tm"
# cd4_tcell.obs["subtype0"][np.isin(cd4_tcell.obs["subtype1"], ['CD4_TIM3+ Trm', 'CD4_TIM3+MKI67+ Trm'])]="CD4_TIM3+ Trm"


# In[8]:


# np.unique(cd4_tcell.obs["subtype0"])


# ## find_variable_genes

# In[9]:


get_ipython().run_cell_magic('capture', '', 'axes = cb.data.find_variable_genes(cd4_tcell, grouping="Fresh")')


# In[10]:


cd4_tcell.var["variable_genes"].sum()


# In[11]:


genes_excluded_all = pd.read_csv("../genes_excluded_all.csv",index_col=0)["x"].values
genes_excluded_all.shape


# In[12]:


# number of HVG in excluded 
np.isin(cd4_tcell.var["variable_genes"].index[cd4_tcell.var["variable_genes"]].values,
        genes_excluded_all).sum()


# In[13]:


# modify HVG
cd4_tcell.var["variable_genes"][np.isin(cd4_tcell.var["variable_genes"].index,genes_excluded_all)]=False


# In[14]:


cd4_tcell.var["variable_genes"].sum()


# In[15]:


cd4_tcell.var["variable_genes"].index[cd4_tcell.var["variable_genes"]].values


# ## one model for visualization (batch correction(fresh/frozen))

# In[16]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodel_rmbatch = cb.directi.fit_DIRECTi(\n    cd4_tcell, genes=cd4_tcell.var.query("variable_genes").index.to_numpy(),\n    batch_effect="Fresh", latent_dim=10, cat_dim=20\n)\ncd4_tcell.obsm["X_latent"] = model_rmbatch.inference(cd4_tcell)')


# In[17]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[18]:


sc.pp.neighbors(cd4_tcell, use_rep="X_latent")
sc.tl.umap(cd4_tcell)


# In[19]:


sc.pl.umap(cd4_tcell, color="subtype1", palette="tab20b")


# In[20]:


sc.pl.umap(cd4_tcell, color="Fresh", palette="tab10")


# In[21]:


model_rmbatch.save("cd4_tcell_one_test/cd4_tcell_rmbatch_model_new")


# ### visualization

# In[7]:


model_rmbatch = cb.directi.DIRECTi.load("./cd4_tcell_one_test/cd4_tcell_rmbatch_model_new")


# In[8]:


cd4_tcell.obsm["X_latent"] = model_rmbatch.inference(cd4_tcell)


# In[9]:


sc.pp.neighbors(cd4_tcell, use_rep="X_latent")
sc.tl.umap(cd4_tcell)


# In[10]:


sc.pl.umap(cd4_tcell, color="subtype1", palette="tab20b")


# In[11]:


sc.pl.umap(cd4_tcell, color="Fresh", palette="tab10")


# In[ ]:





# ## multi models and blast model

# In[22]:


get_ipython().run_cell_magic('capture', '', 'start_time=time.time()\nmodels = []\nfor i in range(4):\n    models.append(cb.directi.fit_DIRECTi(\n        cd4_tcell, genes=cd4_tcell.var.query("variable_genes").index.to_numpy(),\n        batch_effect="Fresh", latent_dim=10, cat_dim=20, random_seed=i\n    ))\n# cd4_tcell.obsm["X_latent"] = model_rmbatch.inference(cd4_tcell)')


# In[23]:


print("Time elapsed: %.1fs" % (time.time() - start_time))


# In[24]:


blast = cb.blast.BLAST(models, cd4_tcell)


# In[25]:


blast.save("./cd4_tcell_blast")


# In[ ]:




