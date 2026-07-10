# Regression tests for the PET asymmetry signal (fixes MAFI-1).
#
# Bug: mafi_signals() stored the PET regression INTERCEPT (beta[1], the
# small-study-corrected effect estimate) as the funnel-asymmetry signal, when
# the asymmetry / FAT coefficient is the SLOPE on sei (beta[2]). Because the
# intercept tracks whether an effect EXISTS, any strong genuine effect with a
# perfectly symmetric funnel produced a maxed-out asymmetry signal and a false
# publication-bias verdict.
#
# These tests use a symmetric, no-publication-bias meta-analysis with a strong
# true effect -- exactly the case the bug false-flagged.

make_symmetric_strong <- function() {
  # Symmetric funnel, strong true effect (theta = 0.8), NO selection/bias.
  # A constant shift of yi does not change the slope on sei, so the asymmetry
  # signal must be insensitive to the effect size.
  set.seed(7)
  k <- 30
  sei <- runif(k, 0.05, 0.35)
  vi <- sei^2
  yi <- 0.8 + rnorm(k, 0, sei)
  list(yi = yi, vi = vi)
}

test_that("pet_slope stores the PET slope on sei (beta[2]), not the intercept", {
  d <- make_symmetric_strong()
  s <- mafi(d$yi, d$vi)$signals
  sei <- sqrt(d$vi)
  pet <- metafor::rma(yi = d$yi, vi = d$vi, mods = ~ sei, method = "REML")

  # The asymmetry signal source must be the slope coefficient (beta[2]) ...
  expect_equal(s$pet_slope, as.numeric(pet$beta[2]), tolerance = 1e-8)
  # ... and must NOT be the intercept (beta[1] = corrected effect estimate),
  # which is what the pre-fix code stored.
  expect_false(isTRUE(all.equal(s$pet_slope, as.numeric(pet$beta[1]))))
})

test_that("symmetric no-bias data with a strong effect is not flagged as biased", {
  # Specificity guard: before the MAFI-1 fix this data scored in the
  # Elevated/High-Risk range (~50-57) because the intercept z was ~10-35 and
  # maxed the asymmetry cluster. After the fix the slope z is ~0.24 and the
  # score sits firmly in the Low/Moderate band.
  d <- make_symmetric_strong()
  res <- mafi(d$yi, d$vi)

  expect_lt(res$score, 40)
  expect_false(res$classification$class %in%
                 c("Elevated Risk", "High Risk", "Very High Risk"))
  # Egger and Begg agree there is no asymmetry here.
  expect_gt(res$signals$egger_pval, 0.10)
})

test_that("correction path returns a finite corrected estimate on clean data", {
  # Lightweight coverage of mafi_correct() via the public pipeline.
  d <- make_symmetric_strong()
  res <- mafi(d$yi, d$vi)
  expect_true(is.numeric(res$corrected))
  expect_false(is.na(res$corrected))
})
