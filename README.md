#  Rapid pairwise alignment of TCR clones

## Installation
```
git clone https://github.com/jbreynier/tcr_aligner.git
cd tcr_aligner
pip install .
```

## Use cases

Need reference set of clones (with separate list for ID and sequences), similar for query set of clones

```
alnr = HashedAln(ref_clone_seq, ref_clone_ids, n_mismatches_allowed)
alnr.build_index()
res = alnr.align_clones(query_clone_seq)
```

A complete example script is available: `src/run_aln.py`

## Description of Files

- `src` : auxiliary scripts to run example alignment 
- `data` : example data obtained from GEO
- `tcraligner` : main code base for hashTCR
- `env.yml` : Required python libraries
- `test_seqs.ipynb` : notebook for in depth benchmarking of hash-based approach
- `vdjdb_analysis.ipynb` : notebook for substitution analysis (amino acid similarity + confidence score)
