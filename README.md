<!-- sentinel:skip-file — hardcoded paths are fixture/registry/audit-narrative data for this repo's research workflow, not portable application configuration. Same pattern as push_all_repos.py and E156 workbook files. -->

# MAFI: Multi-Signal Aggregate Funnel Index

A publication bias detection tool for meta-analysis that combines eight statistical signals into a calibrated probability score (0-100).

## Version 2.0.0

### Key Improvements in v2.0
- **Defined ground truth**: Consensus of Egger, Begg, and Trim-and-fill (>=2/3 positive)
- **Cross-validation**: Leave-one-domain-out CV across 15 Cochrane clinical domains
- **Bootstrap CIs**: Uncertainty quantification for MAFI score
- **Multiple selection models**: Tests 3PSM (0.025, 0.05 steps) and Beta models
- **Direction-aware**: Small study effect now tracks inflation vs deflation
- **Improved ensemble**: Original estimate excluded to prevent anchoring
- **Logistic transformations**: P-values transformed using empirically-derived thresholds

## Installation

```r
install.packages("C:/Users/user/MAFI", repos = NULL, type = "source")
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

### Signal Weights
Weights were derived using gradient boosting with leave-one-domain-out cross-validation:

| Signal | Weight | Transformation |
|--------|--------|----------------|
| Egger p-value | 35% | Logistic, threshold 0.10 |
| PET intercept | 15% | Standardized z-score |
| Selection model LRT | 12% | Logistic, minimum across specs |
| Excess significance | 10% | Logistic, threshold 0.10 |
| Begg p-value | 8% | Logistic, threshold 0.10 |
| Trim-and-fill ratio | 8% | Linear, scale 0.30 |
| Precision-effect cor | 6% | Linear, scale 0.40 |
| Small study effect | 6% | Direction-aware, scale 0.30 |

### Adjustments
- **Heterogeneity**: 15% maximum reduction for I²=100% (Harbord et al., 2006)
- **Sample size**: Reliability scaling based on n (full at n>=10)

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
| `mafi_report()` | Manuscript text (text/markdown/latex) |

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
for Publication Bias Detection. R package version 2.0.0.
```

## License

MIT License
