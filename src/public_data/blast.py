###!/usr/bin/python3

import time
import os
import sys
import warnings
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import Cell_BLAST as cb

from numba import njit

njit(cache=True)


# 从命令行读取变量并赋值给gsm_name
gsm_name = sys.argv[1]
path=sys.argv[2]
blast_model = sys.argv[3] #cd8_tcell_blast
column=sys.argv[4] # subtype1
majority_threshold=float(sys.argv[5])
n_jobs=int(sys.argv[6]) # config.N_JOBS

warnings.filterwarnings("ignore")
cb.config.N_JOBS = n_jobs
cb.config.RANDOM_SEED = 0

# +
print("Load model...")
start_time = time.time()
blast = cb.blast.BLAST.load("blast_model/"+blast_model)
print("Time elapsed for model loading: %.1fs" % (time.time() - start_time))

# 读取gsm_name文件内容为一个matrix
print("Read query data...")
start_time = time.time()
try:
    query_data_0 = anndata.read_h5ad(path)
    #query_data_0 = anndata.read_h5ad("data/" + gsm_name + ".h5ad")
except FileNotFoundError:
    print("File not found.")
    sys.exit(1)
print("Time elapsed for query data loading: %.1fs" % (time.time() - start_time))
    
print(gsm_name, "contains", query_data_0.shape[0], "cells")

query_data = query_data_0.copy()
query_data.X = query_data_0.raw.X
query_data.var_names = query_data.var["gene_symbol"]
query_data.var_names = query_data.var_names.astype('object')
query_data.var_names_make_unique()
query_data

tmp = query_data.X[0:50,100:150].todense()
# 判断matrix是否全为整数，如果包含小数，则停止运行
if not np.all(np.equal(tmp, tmp.astype(int))):
    print("Data contains non-integer values. Stopping the program.")
    sys.exit(1)

print("Querying...")
start_time = time.time()
query_data_hits = blast.query(query_data)
print("Time per query: %.1fms" % (
    (time.time() - start_time) * 1000 / query_data.shape[0]
))
print("Time elapsed: %.1fs" % (time.time() - start_time))

query_data_hits = query_data_hits.reconcile_models().filter(by="pval", cutoff=0.05)

query_data_predictions = query_data_hits.annotate(column,majority_threshold=majority_threshold)

out_path = "predictions/"+blast_model + "_" +column + "_" + str(majority_threshold)
if not os.path.exists(out_path):
    os.mkdir(out_path)

query_data_predictions.to_csv(out_path +"/" + gsm_name + "_pred.csv")


print("Time elapsed: %.1fs" % (time.time() - start_time))

