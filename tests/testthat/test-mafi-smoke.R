# Minimal smoke tests for the MAFI package.
# These exercise the public pipeline end-to-end on the bundled BCG example data
# and assert basic structural / range invariants. They do NOT pin exact numeric
# values (those depend on metafor internals); they guard against crashes,
# regressions in the public API, and out-of-range scores.

test_that("package version is v2.2.0 (v2.1 legacy is not sourced)", {
  # mafi.R defines version = "2.2.0"; the superseded mafi_v21.R now lives under
  # inst/legacy/ and must NOT shadow it.
  res <- mafi(c(0.2, 0.3, -0.1, 0.4, 0.1, 0.25), c(0.04, 0.03, 0.05, 0.02, 0.06, 0.03))
  expect_equal(res$version, "2.2.0")
})

test_that("mafi() runs on the BCG example data and returns a valid score", {
  data(mafi_example_bcg)
  res <- mafi(mafi_example_bcg$yi, mafi_example_bcg$vi)

  expect_s3_class(res, "mafi")
  expect_true(is.numeric(res$score))
  expect_false(is.na(res$score))
  expect_gte(res$score, 0)
  expect_lte(res$score, 100)
  expect_equal(res$n_studies, length(mafi_example_bcg$yi))
})

test_that("classification returns a non-NA confidence in the non-bootstrap path", {
  # Regression guard: previously the non-bootstrap caller passed ci_width = NA,
  # which made mafi_classify() return confidence = NA instead of an n-based level.
  data(mafi_example_bcg)
  res <- mafi(mafi_example_bcg$yi, mafi_example_bcg$vi)
  expect_true(res$classification$confidence %in% c("Low", "Moderate", "High"))
  expect_false(is.na(res$classification$confidence))
})

test_that("fewer than 5 studies yields NA score, not an error", {
  res_score <- mafi_score(c(0.1, 0.2, 0.3), c(0.04, 0.03, 0.05))
  expect_true(is.na(res_score))
})

test_that("mafi_classify handles NA score (insufficient data) gracefully", {
  cls <- mafi_classify(NA_real_, n_studies = 3)
  expect_equal(cls$class, "Unknown")
  expect_equal(cls$confidence, "None")
})

test_that("mafi_report produces text and markdown without error", {
  data(mafi_example_bcg)
  res <- mafi(mafi_example_bcg$yi, mafi_example_bcg$vi)
  expect_type(mafi_report(res, "text"), "character")
  expect_type(mafi_report(res, "markdown"), "character")
  expect_error(mafi_report(res, "latex"))
})
