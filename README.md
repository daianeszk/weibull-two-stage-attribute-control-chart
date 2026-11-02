# A Two-Stage Attribute Control Chart for Monitoring the Mean of Weibull-Distributed Lifetime Processes

This repository provides R scripts associated with the research article **"A Two-Stage Attribute Control Chart for Monitoring the Mean of Weibull-Distributed Lifetime Processes"** by Daiane S. Guimarães, Linda Lee Ho, and Roberto C. Quinino.

The scripts investigate threshold-based attribute control charts for detecting increases and decreases in the mean of Weibull-distributed lifetime processes. They evaluate two-stage conditional sampling and a one-sampling special operational case, and compare out-of-control average run lengths (ARL1) with published empirical Weibull X-bar chart benchmarks.

## Scripts

| File | Purpose |
| --- | --- |
| `Xbar_UCL_Comparison.R` | Two-stage attribute-chart design and comparison for increases in the Weibull process mean. |
| `Xbar_LCL_Comparison.R` | Two-stage attribute-chart design and comparison for decreases in the Weibull process mean. |
| `SpecialCaseUCL.R` | One-sampling special case (first-stage sample unavailable) for mean increases. |
| `SpecialCaseLCL.R` | One-sampling special case (first-stage sample unavailable) for mean decreases. |

## Requirements

- R (version used in the original analysis has not yet been documented).
- The scripts use functions from base R and its standard `stats` package, including `qbeta`, `qweibull`, `pweibull`, and `pbinom`.

## Running the scripts

Run each script separately from a writable working directory, for example:

```bash
Rscript Xbar_UCL_Comparison.R
Rscript Xbar_LCL_Comparison.R
Rscript SpecialCaseUCL.R
Rscript SpecialCaseLCL.R
```

Each script writes CSV results to the current working directory. These generated CSV files are not included in this repository.

## Numerical design and reproducibility

The numerical study was conducted in multiple runs, changing the global parameters at the beginning of the two-stage comparison scripts to evaluate manageable subsets of candidate configurations and avoid excessively large output files. The uploaded scripts retain one set of parameter values per file; they are **not** a single-command reproduction of every table in the article.

Relevant adjustable parameters include `ASS_min`, `ASS_max`, `ASS0_allowed`, `n_b_max`, `alpha_a_vec`, `k0_vec`, and `d_vec` (as applicable). The two-stage scripts calculate the in-control expected sample size as `ASS0 = n_a + n_b * alpha_a` and evaluate out-of-control ARL for the specified Weibull shape parameters and proportional mean shifts. They also calculate out-of-control expected sample size internally, although their selected CSV output columns emphasize ASS0 and ARL1.

For the special-case scripts, the first-stage sample is not collected, the remaining sample is mandatory, and its false-alarm probability is calibrated directly to 0.0027. The special-case sample sizes are 3, 5, 7, 10, 15, and 20.

All four scripts use an overall in-control false-alarm probability of 0.0027, corresponding to a target ARL0 of approximately 370.37. The original comparison benchmarks are embedded in the scripts. Please consult the associated article for the full study design, selection criteria, and interpretation of the published tables.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). A version-specific Zenodo DOI will be added after the first archived release.

## License

A software license has not yet been selected by the authors. Confirm the intended license before public release.
