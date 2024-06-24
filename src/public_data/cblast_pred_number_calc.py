#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os


# ### calc number for all subtypes in a blast model

# In[2]:


# pred_dir = "myeloid_blast_subtype1_0.4"
# pred_dir = "myeloid_blast_subtype1_0.5"
# pred_dir = "cd8_tcell_blast_subtype1_0.4"
# pred_dir = "cd8_tcell_blast_subtype1_0.5"
# pred_dir = "cd4_tcell_blast_subtype1_0.4"
pred_dir = "tcell_blast_subtype1_0.4"


# In[27]:


back_cell = pred_dir.split("_blast")[0]

result_df = pd.DataFrame(columns=["GSM", "Subtype","Count"])

csv_directory = "../cblast/eACA_predictions/"+pred_dir

for csv_file in os.listdir(csv_directory):
    if csv_file.endswith("_pred.csv"):
        file_path = os.path.join(csv_directory, csv_file)
        
        gsm_name = csv_file.split("_")[0]
        print(gsm_name)
            
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            df = pd.read_csv(file_path)
            value_counts = df['subtype1'].value_counts()

            for value, count in value_counts.items():
                result_df = pd.concat([result_df, pd.DataFrame([{"GSM": gsm_name, 
                                      'Subtype': value, 
                                      'Count': count}])], ignore_index=True)

out_csv = "subtype_cblast_num/"+pred_dir+"_num.csv"
result_df.to_csv(out_csv, index=False)

