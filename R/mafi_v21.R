#' MAFI: Multi-Signal Aggregate Funnel Index v2.1.0
#'
#' A publication bias detection tool that combines multiple statistical signals
#' into a calibrated probability score (0-100). Version 2.1 uses simulation-based
#' validation with known ground truth, avoiding circular validation.
#'
#' @name MAFI-package
#' @docType package
#' @author MAFI Development Team
#' @keywords package
NULL

# ============================================================
# MAFI v2.1.0 - EDITORIAL REVISION
# ============================================================
# KEY CHANGES from v2.0:
# 1. Ground truth: Simulation-based (known bias), NOT consensus-based
# 2. Logistic k: Changed from 5 to 4 (justified via ROC optimization)
# 3. Bonferroni correction: Applied for 3 selection model tests
# 4. Direction treatment: Non-binary (0.25*magnitude for deflation)
# 5. PET-PEESE: Always use PEESE when k >= 20, otherwise weighted
# 6. "Limit meta-analysis": Renamed to WLS extrapolation
# 7. Bootstrap: Increased minimum, added failure rate tracking
# 8. Confidence: Now considers CI width, not just n
# ============================================================

.mafi_weights <- list(
  egger_pval = 0.25,
  pet_intercept = 0.15,
  sel_lrt_pval = 0.15,
  excess_sig_pval = 0.12,
  begg_pval = 0.10,
  taf_k0_ratio = 0.10,
  precision_effect_cor = 0.07,
  small_study_effect = 0.06
)

.mafi_params <- list(
  egger_threshold = 0.10,
  begg_threshold = 0.10,
  sel_threshold = 0.10,
  excess_threshold = 0.10,
  logistic_k = 4,
  sel_bonferroni_n = 3,
  taf_ratio_scale = 0.30,
  excess_ratio_scale = 1.50,
  precision_cor_scale = 0.40,
  small_effect_scale = 0.30,
  deflation_weight = 0.25,
  i2_penalty_max = 0.15,
  n_reliability_threshold = 10,
  min_successful_boots = 100,
  boot_failure_tolerance = 0.30
)

#' @export
mafi_signals <- function(yi, vi, sei = NULL) {
  if (is.null(sei)) sei <- sqrt(vi)
  n <- length(yi)

  signals <- list(
    n_studies = n,
    re_estimate = NA, re_se = NA, re_pval = NA,
    tau2 = NA, i_squared = NA,
    egger_z = NA, egger_pval = NA,
    pet_intercept = NA, pet_intercept_se = NA,
    begg_tau = NA, begg_pval = NA,
    taf_k0 = NA, taf_k0_ratio = NA, taf_estimate = NA, taf_side = NA,
    sel_estimate = NA, sel_lrt_pval = NA,
    sel_estimate_05 = NA, sel_lrt_pval_05 = NA,
    sel_beta_estimate = NA, sel_beta_pval = NA,
    sel_min_pval = NA, sel_min_pval_adjusted = NA,
    obs_sig = NA, exp_sig = NA, excess_sig_ratio = NA, excess_sig_pval = NA,
    small_study_effect = NA, small_study_direction = NA, small_study_magnitude = NA,
    precision_effect_cor = NA
  )

  if (n < 5) {
    signals$error <- "Fewer than 5 studies"
    class(signals) <- c("mafi_signals", "list")
    return(signals)
  }

  tryCatch({
    re <- metafor::rma(yi = yi, vi = vi, method = "REML")
    theta <- as.numeric(re$beta)
    signals$re_estimate <- theta
    signals$re_se <- re$se
    signals$re_pval <- re$pval
    signals$tau2 <- re$tau2
    signals$i_squared <- re$I2

    tryCatch({
      egger <- metafor::regtest(re, model = "lm")
      signals$egger_z <- as.numeric(egger$zval)
      signals$egger_pval <- egger$pval
    }, error = function(e) {})

    tryCatch({
      pet <- metafor::rma(yi = yi, vi = vi, mods = ~ sei, method = "REML")
      signals$pet_intercept <- as.numeric(pet$beta[1])
      signals$pet_intercept_se <- pet$se[1]
    }, error = function(e) {})

    tryCatch({
      begg <- metafor::ranktest(re)
      signals$begg_tau <- as.numeric(begg$tau)
      signals$begg_pval <- begg$pval
    }, error = function(e) {})

    tryCatch({
      taf <- metafor::trimfill(re)
      signals$taf_k0 <- taf$k0
      signals$taf_k0_ratio <- taf$k0 / n
      signals$taf_estimate <- as.numeric(taf$beta)
      signals$taf_side <- taf$side
    }, error = function(e) {})

    tryCatch({
      re_ml <- metafor::rma(yi = yi, vi = vi, method = "ML")
      sel_pvals <- c()

      tryCatch({
        sel1 <- metafor::selmodel(re_ml, type = "stepfun", steps = c(0.025))
        signals$sel_estimate <- as.numeric(sel1$beta)
        signals$sel_lrt_pval <- sel1$LRTp
        sel_pvals <- c(sel_pvals, sel1$LRTp)
      }, error = function(e) {})

      tryCatch({
        sel2 <- metafor::selmodel(re_ml, type = "stepfun", steps = c(0.05))
        signals$sel_estimate_05 <- as.numeric(sel2$beta)
        signals$sel_lrt_pval_05 <- sel2$LRTp
        sel_pvals <- c(sel_pvals, sel2$LRTp)
      }, error = function(e) {})

      if (n >= 10) {
        tryCatch({
          sel3 <- metafor::selmodel(re_ml, type = "beta")
          signals$sel_beta_estimate <- as.numeric(sel3$beta)
          signals$sel_beta_pval <- sel3$LRTp
          sel_pvals <- c(sel_pvals, sel3$LRTp)
        }, error = function(e) {})
      }

      if (length(sel_pvals) > 0) {
        signals$sel_min_pval <- min(sel_pvals)
        signals$sel_min_pval_adjusted <- min(1, min(sel_pvals) * length(sel_pvals))
      }
    }, error = function(e) {})

    tryCatch({
      z_crit <- qnorm(0.975)
      power <- pnorm(abs(theta) / sei - z_crit)
      power[power < 0.05] <- 0.05
      power[power > 0.95] <- 0.95
      obs_sig <- sum(abs(yi / sei) > z_crit)
      exp_sig <- sum(power)
      signals$obs_sig <- obs_sig
      signals$exp_sig <- exp_sig
      signals$excess_sig_ratio <- obs_sig / max(exp_sig, 1)
      if (exp_sig >= 1 && exp_sig <= n - 1) {
        bt <- binom.test(obs_sig, n, p = exp_sig / n, alternative = "greater")
        signals$excess_sig_pval <- bt$p.value
      }
    }, error = function(e) {})

    tryCatch({
      tercile <- quantile(sei, c(1/3, 2/3))
      small_idx <- sei > tercile[2]
      large_idx <- sei < tercile[1]
      if (sum(small_idx) >= 2 && sum(large_idx) >= 2) {
        small_eff <- weighted.mean(yi[small_idx], 1/vi[small_idx])
        large_eff <- weighted.mean(yi[large_idx], 1/vi[large_idx])
        signals$small_study_effect <- small_eff - large_eff
        signals$small_study_magnitude <- abs(small_eff - large_eff)
        if (theta > 0) {
          signals$small_study_direction <- ifelse(small_eff > large_eff, "inflated", "deflated")
        } else if (theta < 0) {
          signals$small_study_direction <- ifelse(small_eff < large_eff, "inflated", "deflated")
        } else {
          signals$small_study_direction <- "neutral"
        }
      }
    }, error = function(e) {})

    tryCatch({
      if (theta >= 0) {
        signals$precision_effect_cor <- cor(1/sei, yi, method = "spearman")
      } else {
        signals$precision_effect_cor <- cor(1/sei, -yi, method = "spearman")
      }
    }, error = function(e) {})

  }, error = function(e) {
    signals$error <- as.character(e)
  })

  class(signals) <- c("mafi_signals", "list")
  return(signals)
}

#' @export
mafi_score <- function(yi, vi, signals = NULL, bootstrap = FALSE, n_boot = 200) {
  if (is.null(signals)) signals <- mafi_signals(yi, vi)
  n <- signals$n_studies
  if (n < 5) {
    if (bootstrap) return(list(score = NA, ci_lower = NA, ci_upper = NA, boot_failure_rate = NA))
    return(NA_real_)
  }

  params <- .mafi_params
  k <- params$logistic_k

  compute_score <- function(sig) {
    probs <- list()

    if (!is.na(sig$egger_pval)) {
      probs$egger <- 1 / (1 + exp(k * (sig$egger_pval - params$egger_threshold)))
    } else probs$egger <- 0.5

    if (!is.na(sig$pet_intercept) && !is.na(sig$pet_intercept_se) && sig$pet_intercept_se > 0) {
      pet_z <- abs(sig$pet_intercept) / sig$pet_intercept_se
      probs$pet <- pmin(1, pet_z / 2)
    } else probs$pet <- 0.5

    if (!is.na(sig$sel_min_pval_adjusted)) {
      probs$sel <- 1 / (1 + exp(k * (sig$sel_min_pval_adjusted - params$sel_threshold)))
    } else if (!is.na(sig$sel_lrt_pval)) {
      probs$sel <- 1 / (1 + exp(k * (sig$sel_lrt_pval - params$sel_threshold)))
    } else probs$sel <- 0.5

    if (!is.na(sig$excess_sig_pval)) {
      probs$excess_sig <- 1 / (1 + exp(k * (sig$excess_sig_pval - params$excess_threshold)))
    } else probs$excess_sig <- 0.5

    if (!is.na(sig$begg_pval)) {
      probs$begg <- 1 / (1 + exp(k * (sig$begg_pval - params$begg_threshold)))
    } else probs$begg <- 0.5

    if (!is.na(sig$taf_k0_ratio)) {
      probs$taf <- pmin(1, sig$taf_k0_ratio / params$taf_ratio_scale)
    } else probs$taf <- 0.5

    if (!is.na(sig$precision_effect_cor)) {
      probs$prec_cor <- pmin(1, pmax(0, -sig$precision_effect_cor / params$precision_cor_scale))
    } else probs$prec_cor <- 0.5

    if (!is.na(sig$small_study_effect) && !is.na(sig$small_study_direction)) {
      effect_magnitude <- pmin(1, abs(sig$small_study_effect) / params$small_effect_scale)
      if (sig$small_study_direction == "inflated") {
        probs$small <- effect_magnitude
      } else if (sig$small_study_direction == "deflated") {
        probs$small <- params$deflation_weight * effect_magnitude
      } else probs$small <- 0.5
    } else probs$small <- 0.5

    weights <- c(egger = .mafi_weights$egger_pval, pet = .mafi_weights$pet_intercept,
                 sel = .mafi_weights$sel_lrt_pval, excess_sig = .mafi_weights$excess_sig_pval,
                 begg = .mafi_weights$begg_pval, taf = .mafi_weights$taf_k0_ratio,
                 prec_cor = .mafi_weights$precision_effect_cor, small = .mafi_weights$small_study_effect)

    prob_vec <- unlist(probs)
    weight_vec <- weights[names(probs)]
    raw_score <- sum(prob_vec * weight_vec) / sum(weight_vec)
    i2_adj <- pmax(1 - params$i2_penalty_max, 1 - params$i2_penalty_max * (sig$i_squared / 100))
    n_adj <- pmin(1, sig$n_studies / params$n_reliability_threshold)
    reliability <- sqrt(n_adj)
    score <- 100 * raw_score * i2_adj * reliability
    return(pmax(0, pmin(100, score)))
  }

  main_score <- round(compute_score(signals), 1)
  if (!bootstrap) return(main_score)

  boot_scores <- numeric(n_boot)
  n_failed <- 0
  for (i in 1:n_boot) {
    idx <- sample(n, n, replace = TRUE)
    boot_signals <- tryCatch(mafi_signals(yi[idx], vi[idx]), error = function(e) NULL)
    if (!is.null(boot_signals) && is.null(boot_signals$error)) {
      boot_scores[i] <- compute_score(boot_signals)
    } else {
      boot_scores[i] <- NA
      n_failed <- n_failed + 1
    }
  }

  boot_failure_rate <- n_failed / n_boot
  boot_scores <- boot_scores[!is.na(boot_scores)]

  if (length(boot_scores) < params$min_successful_boots) {
    return(list(score = main_score, ci_lower = NA, ci_upper = NA,
                boot_failure_rate = boot_failure_rate, warning = "Insufficient bootstrap samples"))
  }

  ci <- quantile(boot_scores, c(0.025, 0.975))
  result <- list(score = main_score, ci_lower = round(ci[1], 1), ci_upper = round(ci[2], 1),
                 boot_failure_rate = round(boot_failure_rate, 3))
  if (boot_failure_rate > params$boot_failure_tolerance) {
    result$warning <- sprintf("High bootstrap failure rate (%.1f%%)", boot_failure_rate * 100)
  }
  return(result)
}

#' @export
mafi_classify <- function(score, n_studies = 15, ci_width = NULL) {
  base_thresholds <- c(25, 40, 55, 70)
  thresholds <- if (n_studies < 10) base_thresholds + 10 else if (n_studies < 15) base_thresholds + 5 else base_thresholds

  if (is.null(ci_width)) {
    confidence <- ifelse(n_studies >= 20, "High", ifelse(n_studies >= 10, "Moderate", "Low"))
  } else {
    ci_conf <- ifelse(ci_width < 15, "High", ifelse(ci_width < 30, "Moderate", "Low"))
    n_conf <- ifelse(n_studies >= 20, "High", ifelse(n_studies >= 10, "Moderate", "Low"))
    conf_order <- c("Low" = 1, "Moderate" = 2, "High" = 3)
    confidence <- names(conf_order)[min(conf_order[n_conf], conf_order[ci_conf])]
  }

  if (is.na(score)) return(list(class = "Unknown", color = "#808080",
    description = "Insufficient data", recommendation = "Cannot evaluate with < 5 studies", confidence = "None"))

  if (score < thresholds[1]) return(list(class = "Low Risk", color = "#28a745",
    description = "Little evidence of publication bias", recommendation = "Effect estimate appears robust", confidence = confidence))
  if (score < thresholds[2]) return(list(class = "Moderate Risk", color = "#ffc107",
    description = "Some signals of potential bias", recommendation = "Consider sensitivity analyses", confidence = confidence))
  if (score < thresholds[3]) return(list(class = "Elevated Risk", color = "#fd7e14",
    description = "Multiple signals suggest publication bias", recommendation = "Report corrected estimates alongside original", confidence = confidence))
  if (score < thresholds[4]) return(list(class = "High Risk", color = "#dc3545",
    description = "Strong evidence of publication bias", recommendation = "Corrected estimate should be prioritized", confidence = confidence))
  return(list(class = "Very High Risk", color = "#721c24",
    description = "Substantial publication bias detected", recommendation = "Effect likely inflated; interpret with extreme caution", confidence = confidence))
}

#' @export
mafi_correct <- function(yi, vi, signals = NULL) {
  sei <- sqrt(vi)
  n <- length(yi)
  if (n < 5) return(list(corrected = NA, prediction_lower = NA, prediction_upper = NA,
                         methods = list(), n_methods = 0, agreement = NA, original = NA, attenuation_pct = NA))

  if (is.null(signals)) signals <- mafi_signals(yi, vi)
  original <- signals$re_estimate
  estimates <- list()

  if (!is.na(signals$taf_estimate)) estimates$trim_fill <- signals$taf_estimate
  if (!is.na(signals$sel_estimate)) estimates$selection_3psm <- signals$sel_estimate
  if (!is.na(signals$sel_beta_estimate)) estimates$selection_beta <- signals$sel_beta_estimate

  tryCatch({
    pet <- metafor::rma(yi = yi, vi = vi, mods = ~ sei, method = "REML")
    peese <- metafor::rma(yi = yi, vi = vi, mods = ~ vi, method = "REML")
    pet_est <- as.numeric(pet$beta[1])
    peese_est <- as.numeric(peese$beta[1])
    if (n >= 20) {
      estimates$peese <- peese_est
    } else {
      pet_weight <- pmax(0, pmin(1, (20 - n) / 10))
      estimates$pet_peese_weighted <- pet_weight * pet_est + (1 - pet_weight) * peese_est
    }
  }, error = function(e) {})

  tryCatch({
    fit <- lm(yi ~ sei, weights = 1/vi)
    estimates$wls_extrapolation <- as.numeric(coef(fit)[1])
  }, error = function(e) {})

  est_vec <- unlist(estimates)
  est_vec <- est_vec[!is.na(est_vec)]

  if (length(est_vec) == 0) return(list(corrected = original, prediction_lower = NA, prediction_upper = NA,
                                        methods = estimates, n_methods = 0, agreement = NA, original = original, attenuation_pct = NA))

  corrected <- median(est_vec)
  if (length(est_vec) >= 2) {
    cv <- sd(est_vec) / (abs(mean(est_vec)) + 0.01)
    agreement <- ifelse(cv < 0.25, "High", ifelse(cv < 0.50, "Moderate", "Low"))
  } else agreement <- "Insufficient"

  attenuation <- NA
  if (!is.na(original) && abs(original) > 0.001) {
    attenuation <- round(100 * (original - corrected) / original, 1)
    attenuation <- pmax(-200, pmin(200, attenuation))
  }

  return(list(corrected = round(corrected, 4), prediction_lower = round(min(est_vec), 4),
              prediction_upper = round(max(est_vec), 4), methods = estimates, n_methods = length(est_vec),
              agreement = agreement, original = round(original, 4), attenuation_pct = attenuation))
}

#' @export
mafi <- function(yi, vi, sei = NULL, measure = "Effect Size", bootstrap = FALSE, n_boot = 200) {
  if (is.null(sei)) sei <- sqrt(vi)
  signals <- mafi_signals(yi, vi, sei)

  if (bootstrap) {
    score_result <- mafi_score(yi, vi, signals, bootstrap = TRUE, n_boot = n_boot)
    score <- score_result$score
    score_ci <- c(score_result$ci_lower, score_result$ci_upper)
    boot_failure_rate <- score_result$boot_failure_rate
    ci_width <- if (!is.na(score_ci[2]) && !is.na(score_ci[1])) score_ci[2] - score_ci[1] else NA
  } else {
    score <- mafi_score(yi, vi, signals)
    score_ci <- c(NA, NA)
    boot_failure_rate <- NA
    ci_width <- NA
  }

  classification <- mafi_classify(score, signals$n_studies, ci_width)
  correction <- mafi_correct(yi, vi, signals)

  result <- list(score = score, score_ci = score_ci, boot_failure_rate = boot_failure_rate,
                 classification = classification, original = signals$re_estimate, original_se = signals$re_se,
                 original_pval = signals$re_pval, corrected = correction$corrected,
                 correction_interval = c(correction$prediction_lower, correction$prediction_upper),
                 attenuation = correction$attenuation_pct, correction_agreement = correction$agreement,
                 n_studies = signals$n_studies, i_squared = signals$i_squared, tau2 = signals$tau2,
                 signals = signals, correction_methods = correction$methods, measure = measure, version = "2.1.0")
  class(result) <- "mafi"
  return(result)
}

#' @export
print.mafi <- function(x, ...) {
  cat("\n======================================================\n")
  cat("  MAFI: Multi-Signal Aggregate Funnel Index v2.1\n")
  cat("======================================================\n\n")
  cat(sprintf("Studies: %d | I^2: %.1f%% | Tau^2: %.4f\n\n", x$n_studies, x$i_squared, x$tau2))

  if (!is.na(x$score_ci[1])) {
    cat(sprintf("MAFI Score: %.1f [95%% CI: %.1f to %.1f]\n", x$score, x$score_ci[1], x$score_ci[2]))
    if (!is.na(x$boot_failure_rate) && x$boot_failure_rate > 0.1) {
      cat(sprintf("  (Bootstrap failure rate: %.1f%%)\n", x$boot_failure_rate * 100))
    }
  } else cat(sprintf("MAFI Score: %.1f / 100\n", x$score))

  cat(sprintf("Risk Level: %s (Confidence: %s)\n\n", x$classification$class, x$classification$confidence))
  cat("Effect Estimates:\n")
  cat(sprintf("  Original:  %.4f (SE: %.4f, p = %.4f)\n", x$original, x$original_se, x$original_pval))
  cat(sprintf("  Corrected: %.4f [%.4f to %.4f]\n", x$corrected, x$correction_interval[1], x$correction_interval[2]))

  if (!is.na(x$attenuation)) {
    direction <- ifelse(x$attenuation > 0, "toward null", "away from null")
    cat(sprintf("  Change: %.1f%% %s (Agreement: %s)\n", abs(x$attenuation), direction, x$correction_agreement))
  }
  cat("\nRecommendation:\n")
  cat(strwrap(x$classification$recommendation, width = 60, prefix = "  "), sep = "\n")
  cat("\n")
  invisible(x)
}

#' @export
summary.mafi <- function(object, ...) {
  cat("\n======================================================\n")
  cat("  MAFI v2.1 DETAILED SUMMARY\n")
  cat("======================================================\n\n")
  cat(sprintf("META-ANALYSIS: %d studies, effect = %.4f (SE: %.4f, p = %.4f)\n",
              object$n_studies, object$original, object$original_se, object$original_pval))
  cat(sprintf("Heterogeneity: I^2 = %.1f%%, Tau^2 = %.4f\n\n", object$i_squared, object$tau2))
  cat(sprintf("MAFI SCORE: %.1f / 100 [%s]\n\n", object$score, object$classification$class))

  s <- object$signals
  cat("SIGNALS:\n")
  cat(sprintf("  Egger: z = %.2f, p = %.4f\n", ifelse(is.na(s$egger_z), NA, s$egger_z), ifelse(is.na(s$egger_pval), NA, s$egger_pval)))
  cat(sprintf("  Begg: tau = %.3f, p = %.4f\n", ifelse(is.na(s$begg_tau), NA, s$begg_tau), ifelse(is.na(s$begg_pval), NA, s$begg_pval)))
  cat(sprintf("  T&F: %d imputed\n", ifelse(is.na(s$taf_k0), 0, s$taf_k0)))
  if (!is.na(s$sel_min_pval_adjusted)) {
    cat(sprintf("  Selection: min p = %.4f (Bonferroni: %.4f)\n", s$sel_min_pval, s$sel_min_pval_adjusted))
  }
  cat(sprintf("  Small study: %.3f (%s)\n\n", ifelse(is.na(s$small_study_effect), NA, s$small_study_effect),
              ifelse(is.na(s$small_study_direction), "NA", s$small_study_direction)))

  cat(sprintf("CORRECTION: Original = %.4f, Corrected = %.4f (%.1f%% change)\n",
              object$original, object$corrected, abs(object$attenuation)))
  cat("Methods:", paste(names(object$correction_methods), collapse = ", "), "\n\n")
  invisible(object)
}

#' @export
mafi_report <- function(x, format = "text") {
  score_text <- sprintf("%.1f/100", x$score)
  risk_text <- tolower(x$classification$class)
  if (format == "text") {
    sprintf("MAFI v2.1 score was %s (%s; k=%d, I^2=%.1f%%). Original estimate: %.3f, corrected: %.3f (%.1f%% change).",
            score_text, risk_text, x$n_studies, x$i_squared, x$original, x$corrected, abs(x$attenuation))
  } else if (format == "markdown") {
    sprintf("**MAFI v2.1:** %s (%s) | Original: %.3f | Corrected: %.3f", score_text, x$classification$class, x$original, x$corrected)
  } else stop("format must be 'text' or 'markdown'")
}

#' @export
simulate_biased_ma <- function(k = 20, theta = 0.3, tau = 0.1, bias_strength = 0.5,
                                n_per_study = c(20, 200), bias_type = "significance") {
  n_gen <- k * 3
  ni <- if (length(n_per_study) == 2) round(runif(n_gen, n_per_study[1], n_per_study[2])) else rep(n_per_study, n_gen)
  theta_i <- rnorm(n_gen, theta, tau)
  sei <- 1 / sqrt(ni)
  yi <- rnorm(n_gen, theta_i, sei)
  vi <- sei^2
  z <- yi / sei
  p <- 2 * pnorm(-abs(z))

  prob_select <- if (bias_type == "significance") {
    ifelse(p < 0.05, 1, 1 - bias_strength * (1 - exp(-5 * p)))
  } else if (bias_type == "direction") {
    ifelse(yi > 0, 1, 1 - bias_strength * 0.8)
  } else if (bias_type == "both") {
    ifelse(p < 0.05, 1, 1 - bias_strength * 0.6) * ifelse(yi > 0, 1, 1 - bias_strength * 0.4)
  } else rep(1, n_gen)

  selected <- runif(n_gen) < prob_select
  idx <- which(selected)[1:min(k, sum(selected))]
  if (length(idx) < 5) idx <- sample(n_gen, min(k, n_gen))

  data.frame(yi = yi[idx], vi = vi[idx], sei = sei[idx], true_theta = theta,
             true_tau = tau, bias_strength = bias_strength, bias_type = bias_type,
             naive_estimate = mean(yi[idx]), true_bias = mean(yi[idx]) - theta)
}

#' @export
validate_mafi <- function(n_sims = 100, k_range = c(10, 20, 50), bias_levels = c(0, 0.25, 0.5, 0.75, 1.0),
                          seed = 42, verbose = TRUE) {
  set.seed(seed)
  results <- data.frame()
  total <- length(k_range) * length(bias_levels) * n_sims
  count <- 0

  for (k in k_range) {
    for (bias in bias_levels) {
      for (i in 1:n_sims) {
        count <- count + 1
        if (verbose && count %% 50 == 0) cat(sprintf("Progress: %d/%d\n", count, total))
        dat <- tryCatch(simulate_biased_ma(k = k, theta = 0.3, tau = 0.1, bias_strength = bias), error = function(e) NULL)
        if (is.null(dat)) next
        mafi_result <- tryCatch(mafi(dat$yi, dat$vi), error = function(e) NULL)
        if (is.null(mafi_result)) next
        results <- rbind(results, data.frame(k = k, true_bias = bias, true_bias_binary = as.numeric(bias > 0.25),
                                              mafi_score = mafi_result$score, original_estimate = mafi_result$original,
                                              corrected_estimate = mafi_result$corrected, true_theta = dat$true_theta[1],
                                              naive_bias = dat$naive_estimate[1] - dat$true_theta[1]))
      }
    }
  }

  if (nrow(results) == 0) return(list(results = results, error = "No successful simulations"))

  results$score_bin <- cut(results$mafi_score, breaks = c(0, 20, 40, 60, 80, 100), include.lowest = TRUE)
  calibration <- aggregate(true_bias ~ score_bin, data = results, FUN = mean)

  pos <- results$mafi_score[results$true_bias_binary == 1]
  neg <- results$mafi_score[results$true_bias_binary == 0]
  auc <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
  rmse <- sqrt(mean((results$mafi_score/100 - results$true_bias)^2))

  results$bias_reduction <- abs(results$naive_bias) - abs(results$corrected_estimate - results$true_theta)
  correction_effectiveness <- mean(results$bias_reduction > 0, na.rm = TRUE)

  list(results = results, calibration = calibration, auc = round(auc, 3), rmse = round(rmse, 3),
       correction_effectiveness = round(correction_effectiveness, 3), n_sims_successful = nrow(results),
       summary = sprintf("MAFI Validation: AUC=%.3f, RMSE=%.3f, Correction=%.1f%%", auc, rmse, correction_effectiveness * 100))
}
