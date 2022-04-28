#  Rapid pairwise alignment of TCR clones

To install
```
git clone https://github.com/jbreynier/tcr_aligner.git
cd tcr_aligner
pip install .
```

to run:

Need Reference set of clones with separate list for ID, and sequnces, and similar for query set of clones 

```
alnr = HashedAln(ref_clone_seq ref_clone_ids, n_misatches_allowed)
alnr.build_index()
res = alnr.align_clones(query_clone_seq)
```