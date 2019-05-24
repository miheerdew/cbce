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
