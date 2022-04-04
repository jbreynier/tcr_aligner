# Code is very inefficient, just initial framework

import pandas as pd
import re
from Bio import pairwise2
import math

def verify_canonical(row):
    match = re.search(r"C[^_]*[FW]", row['cdr3aa'])
    if not match:
        row['cdr3aa'] = None
        row['cdr3dna'] = None
    else:
        row['cdr3aa'] = match.group()
        row['cdr3dna'] = row['cdr3dna'][match.span()[0]*3:match.span()[1]*3]
    return row

def pairwise2_list(list_prev, seq):
    align_score = -math.inf
    best_sequence = ""
    for prev_seq in list_prev:
        current_score = pairwise2.align.localxs(prev_seq, seq, -1, -1, score_only=True)
        if current_score > align_score:
            align_score = current_score
            best_sequence = prev_seq
    return (align_score, best_sequence)

def multiple_pairwise2(row_prev_seq, list_curr, list_columns):
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