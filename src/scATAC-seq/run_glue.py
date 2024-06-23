#!/usr/bin/env python
# coding: utf-8

# In[4]:


import anndata
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import numpy as np


# In[5]:


scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# In[6]:


import os
# os.chdir("ATAC_with_B4/")
os.getcwd()


# ## Read preprocessed data

# In[7]:


rna = anndata.read_h5ad("GLUE/rna_preprocessed.h5ad")
atac = anndata.read_h5ad("GLUE/atac_preprocessed.h5ad")
graph = nx.read_graphml("GLUE/prior.graphml.gz")


# ## Configure data

# In[8]:


scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca", use_batch="Fresh"
)


# In[9]:


scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi", use_batch="patient"
)


# In[10]:


graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))


# ## Build and train GLUE model

# In[11]:


# import importlib; importlib.reload(scglue)


# In[12]:


glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "GLUE/glue"}
)


# In[13]:


glue.save("GLUE/glue.dill")


# In[14]:


# glue = scglue.models.load_model("GLUE/glue.dill")


# ## Check integration diagnostics

# In[15]:


dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, graph,
    count_layers={"rna": "counts"}
)
dx


# In[16]:


_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")


# ## Apply model for cell and feature embedding

# In[17]:


rna.obsm["X_glue"] = glue.encode_data("rna", rna)# domain name means the key name when training model
atac.obsm["X_glue"] = glue.encode_data("atac", atac)


# In[18]:


combined = anndata.concat([rna, atac], join="outer")


# In[19]:


sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)


# In[20]:


combined


# In[21]:


# combined.obs["cell_type"]=combined.obs["cell_type"].astype("str")


# In[22]:


# combined.obs.loc[combined.obs["domain"]=="scATAC-seq", "cell_type"]="NA"


# In[23]:


combined.shape


# In[24]:


atac


# In[25]:


# sc.pl.umap(combined[combined.obs["domain"]=='scRNA-seq', :], color=["cell_type", "sample"], wspace=0.65)


# In[26]:


sc.pl.umap(combined, color=["cell_type1", "domain"], wspace=0.65)


# In[27]:


sc.pl.umap(combined, color=["subtype1"], wspace=0.65)


# In[28]:


# sc.pl.umap(combined[combined.obs["domain"]=='scRNA-seq', :], color=["subtype"], wspace=0.65)


# ## transfer labels (cell_type1)

# In[29]:


atac.obs=atac.obs.rename(columns={"cell_type": "cell_type_archr"})
atac.obs["sample"]=atac.obs["Sample"]


# In[30]:


atac


# In[31]:


# sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
# sc.tl.umap(atac)


# In[32]:


scglue.data.transfer_labels(rna, atac, "cell_type1", n_neighbors=15, use_rep="X_glue")


# In[33]:


atac


# In[34]:


sc.pl.umap(atac, color=["cell_type1","cell_type_archr", "predictedGroup_Un", "sample"], wspace=0.65,ncols=2)


# In[ ]:





# In[35]:


# 用自己帮忙修正注释结果
scglue.data.transfer_labels(atac, atac, "cell_type1", n_neighbors=15, use_rep="X_lsi")


# In[36]:


sc.pl.umap(atac, color=["cell_type1","cell_type_archr", "predictedGroup_Un", "sample"], wspace=0.65,ncols=2)


# In[ ]:





# ## transfer labels (cell_type)

# In[49]:


scglue.data.transfer_labels(rna, atac, "cell_type", n_neighbors=15, use_rep="X_glue")


# In[50]:


atac


# In[51]:


sc.pl.umap(atac, color=["cell_type","cell_type_archr", "predictedGroup_Un", "sample"], wspace=0.65,ncols=2)


# In[ ]:





# In[52]:


# 用自己帮忙修正注释结果
scglue.data.transfer_labels(atac, atac, "cell_type", n_neighbors=15, use_rep="X_lsi")


# In[53]:


sc.pl.umap(atac, color=["cell_type","cell_type_archr", "predictedGroup_Un", "sample"], wspace=0.65,ncols=2)


# ## transfer labels-subtype

# In[37]:


scglue.data.transfer_labels(rna, atac, "subtype1",n_neighbors=15, use_rep="X_glue")


# In[38]:


atac


# In[39]:


sc.pl.umap(atac, color=["subtype1","sample"], wspace=0.65)


# In[ ]:





# In[40]:


scglue.data.transfer_labels(atac, atac, "subtype1",use_rep="X_lsi",n_neighbors=15)


# In[41]:


np.unique(atac.obs["subtype1"], return_counts=True)


# In[42]:


sc.pl.umap(atac, color=["subtype1","sample"], wspace=0.65)


# In[43]:


rna.write("GLUE/rna_glue.h5ad", compression="gzip")
atac.write("GLUE/atac_glue.h5ad", compression="gzip")
atac.obs.to_csv("GLUE/atac_ct_glue.csv")


# In[54]:


atac.obs.to_csv("GLUE/atac_ct_glue.csv")


# ## graph

# In[44]:


# To obtain feature embeddings
feature_embeddings = glue.encode_graph(graph)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]


# In[45]:


feature_embeddings.to_csv(f"GLUE/graph/feature_embeddings.csv", index=True, header=False)


# In[46]:


genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)


# In[47]:


dist_graph = scglue.genomics.window_graph(
    promoters, peaks, 150000,
    attr_fn=lambda l, r, d: {
        "dist": abs(d),
        "weight": scglue.genomics.dist_power_decay(abs(d)),
        "type": "dist"
    }
)
dist_graph = nx.DiGraph(dist_graph)
dist_graph.number_of_edges()


# In[48]:


nx.write_graphml(dist_graph, f"GLUE/graph/dist.graphml.gz")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




