# Other manuscript plots

Figure-generation scripts retained for revised-manuscript plots that use compact derived inputs.

## UpSet and Venn plots

`upset_plot.py` and `venn_symmetric.py` use a local tab-separated file:

`upset_intersections.tsv`

Required header:

```text
fail_16S	fail_FastANI	fail_Mash	n	degree	combination
```

The three `fail_*` columns are Boolean membership indicators for literature-threshold discordance; `n` is the mutually exclusive intersection count and `degree` is the number of active sets.

## Complete-genome 16S copy-number plot

`plot_complete_genome_16S_copy_number.py` uses:

`CompleteGenome_per_GCA_full_length16S.tsv`

Required columns:

```text
Domain	Full-length 16S count
```

`Domain` should contain `Bacteria` or `Archaea`. Large derived input tables are not committed to GitHub; these schemas document the inputs needed to reproduce the figures.
