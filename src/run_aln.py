from tcraligner.alignment_functions import HashedAln
import pandas as pd
from Bio import Seq
import numpy as np 
import time 


df = pd.read_csv("../data/GSE179994/GSE179994_PBMC.bulkTCR.tsv", sep = '\t')
print(df.head())

df = df.assign(clone_nt = df['TRB'].str.extract(r"(_\w+_)").iloc[:,0].str.replace("_", "").tolist())

df = df.assign(clone_aa =  [Seq.translate(i) for i in df.clone_nt.tolist()])

## check cannonical seqs
df = df[df.clone_aa.str.contains(r"^C")].pipe(lambda x: x[x.clone_aa.str.contains(r"F$")])

df_pre = df[df['sample'].str.contains("pre")]
df_post = df[df['sample'].str.contains("post")]

clones_pre_aa = df_pre[~df_pre.clone_aa.isin(df_post['clone_aa'])]['clone_aa'].tolist()
clones_post_aa = df_post[~df_post.clone_aa.isin(df_pre['clone_aa'])]['clone_aa'].tolist()
clones_pre_ids = df_pre[~df_pre.clone_aa.isin(df_post['clone_aa'])]['cloneId'].tolist()
alnr = HashedAln(clones_pre_aa, clones_pre_ids, 6)
alnr.build_index()
res = alnr.align_clones(clones_pre_aa[:100])
