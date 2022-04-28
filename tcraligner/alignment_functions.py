# Code is very inefficient, just initial framework
#%%
import pandas as pd

from Bio import pairwise2
import math

#%%
# def verify_canonical(row):
#     match = re.search(r"C[^_]*[FW]", row['cdr3aa'])
#     if not match:
#         row['cdr3aa'] = None
#         row['cdr3dna'] = None
#     else:
#         row['cdr3aa'] = match.group()
#         row['cdr3dna'] = row['cdr3dna'][match.span()[0]*3:match.span()[1]*3]
#     return row

def _pairwise2_list(list_prev, seq):
    align_score = -math.inf
    best_sequence = ""
    for prev_seq in list_prev:
        current_score = pairwise2.align.localxs(prev_seq, seq, -1, -1, score_only=True)
        if current_score > align_score:
            align_score = current_score
            best_sequence = prev_seq
    return (align_score, best_sequence)





def _multiple_pairwise2(row_prev_seq, list_curr, list_columns):
    align_score = -math.inf
    best_prev = ""
    best_curr = ""
    list_curr = set(list_curr)
    for column in list_columns:
        # quick fix: check if instance of list
        if isinstance(row_prev_seq[column], list):
            list_prev = set(row_prev_seq[column])
            for prev_seq in list_prev:
                for curr_seq in list_curr:
                    current_score = pairwise2.align.localxs(prev_seq, curr_seq, -1, -1, score_only=True)
                    if current_score > align_score:
                        align_score = current_score
                        best_prev = prev_seq
                        best_curr = curr_seq
    if align_score > -math.inf:
        return (align_score, best_prev, best_curr)
    else:
        return float("NaN")


class HashedAln:
    def __init__(self, ref_clones, clone_ids, kmer_size):
        ## decided which seq to index( one with more seqs)
        self.reference_clones = {k:v for k,v in zip(clone_ids, ref_clones)}
        self._reference_clones_seq = ref_clones
        self.kmer_size = kmer_size
    def build_index(self):
        self.index = {}
        for clone_id, seq in self.reference_clones.items():
            kmers = [seq[i:(i + self.kmer_size)] for i in range(0,len(seq), self.kmer_size)]
            for km in kmers:
                if km in self.index:
                    self.index[km].append(clone_id)
                else:
                    self.index[km]=[clone_id]
    def _align_clones_hash(self, query_clones):
        aligned_clones = []
        unaligned_clones = []
        for clone in query_clones:
            hits = set()
            kmers = [clone[i:(i + self.kmer_size)] for i in range(0,len(clone), self.kmer_size)]
            for km in kmers:
                if km in self.index:
                    hits.update(self.index[km])
            if len(hits) >0:
                hits_seq = [self.reference_clones[k] for k in hits ]
                aligned_clones.append(_pairwise2_list(hits_seq, clone))
            else:
                unaligned_clones.append(clone)
        return aligned_clones, unaligned_clones
    def _align_clones_pairwise(self, query_clones):
        aligned_clones = [ _pairwise2_list(self._reference_clones_seq, clone ) 
            for clone in query_clones ]
        return aligned_clones
    def align_clones(self, query_clones):
        print("Aligning with hash")
        aligned, unaligned = self._align_clones_hash(query_clones)
        print(f"{len(unaligned)} failed to align with hash, running standard pairwise alignment")
        aligned += self._align_clones_pairwise(unaligned)
        print("Done")
        return aligned
# 

# %%