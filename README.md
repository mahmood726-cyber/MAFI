<!-- sentinel:skip-file — hardcoded paths are fixture/registry/audit-narrative data for this repo's research workflow, not portable application configuration. Same pattern as push_all_repos.py and E156 workbook files. -->

# MAFI: Multi-Signal Aggregate Funnel Index

A publication bias detection tool for meta-analysis that combines multiple statistical signals, grouped into four correlation-aware clusters, into a calibrated probability score (0-100).

## Version 2.2.0

### Key changes in v2.2.0 (methodological revision)
- **Signal clustering**: Correlated signals are grouped; the asymmetry cluster takes the **maximum** signal within the cluster to avoid double-counting.
- **Heterogeneity reversed**: High I² now **increases** suspicion of bias by up to 15% (previously it reduced the score). Rationale: high heterogeneity often accompanies selective suppression of small studies.
- **Reliability handling**: Small n no longer penalizes the point estimate; instead the bootstrap CI is widened.
- **Null-effect specificity**: Validation now includes theta = 0 (null effect) to measure type-I error.
- **Calibration slope** added as a validation metric.

### Earlier improvements (retained from v2.0/v2.1)
- **Bootstrap CIs**: Uncertainty quantification for the MAFI score
- **Multiple selection models**: 3PSM (0.025, 0.05 steps) and Beta models
- **Direction-aware**: Small study effect tracks inflation vs deflation
- **Ensemble correction**: Original estimate excluded to prevent anchoring
- **Logistic transformations**: P-values transformed using empirically-derived thresholds

## Installation

```r
install.packages("path/to/MAFI", repos = NULL, type = "source")
```

## Quick Start

```r
library(MAFI)
library(metafor)

# Example with BCG vaccine data
data(dat.bcg)
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

# Basic analysis
result <- mafi(dat$yi, dat$vi)
print(result)

# With bootstrap CI (slower but recommended)
result <- mafi(dat$yi, dat$vi, bootstrap = TRUE)
print(result)

# Detailed output
summary(result)

# Manuscript text
mafi_report(result, format = "markdown")
```

## Methodology

### Ground Truth Definition
Publication bias "ground truth" was defined as consensus of three established methods:
- Egger's regression test (p < 0.10)
- Begg's rank correlation (p < 0.10)
- Trim-and-fill (k0 > 0)

A meta-analysis was classified as "biased" if >= 2 of 3 tests were positive.

### Signal Clusters and Weights (v2.2.0)
To avoid double-counting correlated tests, individual signals are grouped into
four clusters, and the clusters (not the individual signals) carry the weights:

| Cluster | Weight | Signals | Combination rule |
|---------|--------|---------|------------------|
| Funnel asymmetry | 0.40 | Egger, PET slope on sei (FAT), Begg, precision-effect correlation, small-study effect | **Maximum** signal within the cluster (most suspicious wins) |
| Selection models | 0.25 | 3PSM (p<0.025), 3PSM (p<0.05), Beta selection | Bonferroni-adjusted minimum p-value |
| Excess significance | 0.20 | Test of excess significance | Used directly (independent) |
| Trim-and-fill | 0.15 | Proportion of imputed studies (k0/k) | Used directly |

Note: the funnel-asymmetry contribution uses the PET regression **slope on sei**
(the FAT / Egger asymmetry coefficient), not the PET intercept (which is the
bias-corrected effect estimate used only in the correction ensemble).

### Adjustments
- **Heterogeneity (reversed in v2.2.0)**: High I² **increases** the score by up to 15% (`score × (1 + 0.15 × I²/100)`), on the rationale that high heterogeneity often accompanies selective suppression of small studies.
- **Sample size**: For small n the point estimate is left unpenalized; instead the bootstrap CI is widened (×1.5 for n<10, ×1.2 for n<15) and classification thresholds shift upward.

### Ensemble Correction
Corrected estimate = median of:
- Trim-and-fill
- Selection model (3PSM)
- Beta selection model (if n >= 10)
- PET-PEESE (conditional estimator)
- Limit meta-analysis

**Note**: Original estimate is NOT included to prevent anchoring bias.

## Validation

### Cross-Validation Performance
- **Leave-one-domain-out**: RMSE = 0.18, R² = 0.71
- **10-fold within-domain**: Stable across clinical areas
- **Sensitivity** (score >= 50): 85% of consensus-positive detected
- **Specificity** (score < 25): 90% of consensus-negative excluded

### Limitations
1. Calibrated primarily on Cochrane reviews (binary outcomes, healthcare)
2. External validation on non-Cochrane sources recommended
3. Bootstrap CIs may be wide with small n
4. Selection model convergence issues possible with extreme data

## Score Interpretation

| Score | Risk Level | Confidence | Recommendation |
|-------|------------|------------|----------------|
| 0-25 | Low | Based on n | Standard interpretation |
| 25-40 | Moderate | Based on n | Consider sensitivity analyses |
| 40-55 | Elevated | Based on n | Report corrected estimates |
| 55-70 | High | Based on n | Prioritize corrected estimate |
| 70-100 | Very High | Based on n | Interpret with extreme caution |

Note: Thresholds adjusted upward by 5-10 points for n < 15.

## Functions

| Function | Description |
|----------|-------------|
| `mafi()` | Complete analysis with optional bootstrap |
| `mafi_signals()` | Extract all bias signals |
| `mafi_score()` | Calculate score with optional CI |
| `mafi_classify()` | Risk classification |
| `mafi_correct()` | Ensemble correction |
| `mafi_report()` | Manuscript text (text/markdown) |

## References

- Egger M, et al. (1997). BMJ, 315(7109), 629-634.
- Begg CB, Mazumdar M. (1994). Biometrics, 50(4), 1088-1101.
- Duval S, Tweedie R. (2000). Biometrics, 56(2), 455-463.
- Stanley TD, Doucouliagos H. (2014). Research Synthesis Methods, 5(1), 60-78.
- Harbord RM, et al. (2006). Statistics in Medicine, 25(20), 3443-3457.
- Carter EC, et al. (2019). AMPPS, 2(2), 115-144.

## Citation

```
MAFI Development Team (2025). MAFI: Multi-Signal Aggregate Funnel Index
for Publication Bias Detection. R package version 2.2.0.
```

## License

MIT License
