#%%
import numpy as np
import os
import fnmatch
import pandas as pd

# %%
file_names = ['PcrA_Myc-PcrA_Sub_Myc_control_median.npy',
              'PcrA_Myc-PcrA_Sub_Myc_control_maxci.npy',
              'PcrA_Myc-PcrA_Sub_Myc_control_mean.npy',
              'PcrA_Myc-PcrA_Sub_Myc_control_mad.npy',
              'PcrA_Myc-PcrA_Sub_Myc_control_minci.npy',
              'PcrA_Myc-PcrA_Sub_Myc_control_var.npy']
# %%
arr_list = []
for fname in file_names:
    arr_list.append(np.load(fname))
dat = np.array(arr_list)
# %%
dat = np.transpose(dat)

# %%
dat.shape
# %%
df = pd.DataFrame(dat)

# %%
df.columns = ['median','maxci','mean','mad','minci','var']

# %%
df['position'] = np.arange(0, dat.shape[1]) * 10 + 5

# %%
df.to_csv('PcrA_over_Myc_control_enrichment_stats.csv')