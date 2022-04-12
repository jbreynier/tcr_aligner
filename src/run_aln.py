#%%
from alignment_functions import pairwise2_list, HashedAln
import pandas as pd
from Bio import Seq
import numpy as np 
import time 


df = pd.read_csv("../data/GSE179994/GSE179994_PBMC.bulkTCR.tsv", sep = '\t')
# %%
df = df.assign(clone_nt = df['TRB'].str.extract(r"(_\w+_)").iloc[:,0].str.replace("_", "").tolist())

df = df.assign(clone_aa =  [Seq.translate(i) for i in df.clone_nt.tolist()])
#%%
## check cannonicla seqs
df = df[df.clone_aa.str.contains(r"^C")].pipe(lambda x: x[x.clone_aa.str.contains(r"F$")])

##

df_pre = df[df['sample'].str.contains("pre")]
df_post = df[df['sample'].str.contains("post")]
#%%

#5817 exact matches 


#%%

clones_pre_aa = df_pre[~df_pre.clone_aa.isin(df_post['clone_aa'])]['clone_aa'].tolist()
clones_post_aa = df_post[~df_post.clone_aa.isin(df_pre['clone_aa'])]['clone_aa'].tolist()


# %%
# import time 

# tdiff = []
# for i in range(100):
#     a= time.time()
#     _ = pairwise2_list(clones_pre_aa, clones_post_aa[i])
#     b=time.time()
#     tdiff.append(b-a)

# # %%
# import numpy as np 
# np.mean(tdiff)
# # %%
# np.mean(tdiff) * len(clones_post_aa) / 3600 # ~ 19 hours for basic pairwise  
# %%
clones_pre_ids = df_pre[~df_pre.clone_aa.isin(df_post['clone_aa'])]['cloneId'].tolist()

# %%
alnr = HashedAln(clones_pre_aa, clones_pre_ids, 6)
# %%
alnr.build_index()
# %%
hits, hit_size = alnr.align_clones(clones_pre_aa[:100])# %%

# %%
n_missing = []
hit_sizes = []
for i in range(5, 11):
    alnr = HashedAln(clones_pre_aa, clones_pre_ids, i)
    istart = time.time()
    alnr.build_index()
    hits, hit_size = alnr.align_clones(clones_pre_aa[:1000])
    hit_sizes.append(np.mean(hit_size))
    n_missing.append(np.sum([1 for x in hits if x is None ]))

## at a kmer size of 9, miss about 10% of samples, and 

# %%
