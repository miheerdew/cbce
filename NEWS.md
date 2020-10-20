# cbce 1.0.0
Moving to a stable verison!

# cbce 0.8.2.9000
The filtering process would produce an error for a single bimodule. Fixed that. 

# cbce 0.8.1

Filter bimodules with scores lower than the BH multiple testing threshold. But this doesn't help with scrambled much. 

# cbce 0.8.0

Th new function `half_permutation_fdr` estimates CBCE's false discovery rate across different values of alpha using half-permutation.

The function `filter_bimodules` filters bimodules for overlaps.


# cbce 0.7.0

Allow singleton bimodules (the KISS principle).

# cbce 0.6.0

Change the output format. `result$comms` will be all unqiue communities, `result$comms.fil` will be the filtered communities. `result$summary` gives a summary.

# cbce 0.5.0

Brought back the heursitic search!

# cbce 0.4.0

Allow covariates directly to be passed to our method.

# cbce 0.3.0

To correct for covariate removal, added a new parameter n.eff. Removal of each covariate, reduces the degrees of freedom and that must be accounted for in the p-value calculations. This doesn't matter when the sample size is large, but it does make a difference when it is small.

# cbce 0.2.5

Use conservative multiple testing correction (BY) by default.

# cbce 0.2.4

Use deterministic starts by default to avoid surpirses

# cbce 0.2.3

Randomly choose initial sets to starts with.

# cbce 0.2.2

Parameter to stop extractions from growing very large.

# cbce 0.2.1

Use local numbering for the extraction index.

# cbce 0.2.0.9000

* Added a `NEWS.md` file to track changes to the package.
