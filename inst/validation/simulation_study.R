# ============================================================
# MAFI Simulation Validation Study
# ============================================================
# Compares MAFI to traditional methods under known publication bias

library(MAFI)
library(metafor)

set.seed(42)

# Simulation parameters
n_sims <- 200
k_values <- c(10, 20, 50)
bias_levels <- c(0, 0.3, 0.6, 0.9)
theta_true <- 0.3
tau_true <- 0.15

cat("\n")
cat("============================================================\n")
cat("  MAFI SIMULATION VALIDATION STUDY\n")
cat("============================================================\n\n")

# Storage
results <- data.frame()

for (k in k_values) {
  for (bias in bias_levels) {
    cat(sprintf("Simulating: k=%d, bias=%.1f\n", k, bias))

    for (sim in 1:n_sims) {
      # Generate meta-analysis with publication bias
      dat <- .simulate_biased_ma(
        k = k,
        theta = theta_true,
        tau = tau_true,
        bias_strength = bias
      )

      if (nrow(dat) < 5) next

      yi <- dat$yi
      vi <- dat$vi

      # Run MAFI
      mafi_result <- tryCatch({
        suppressWarnings(mafi(yi, vi))
      }, error = function(e) NULL)

      if (is.null(mafi_result)) next

      # Traditional methods
      egger_p <- mafi_result$signals$egger_pval
      begg_p <- mafi_result$signals$begg_pval
      taf_k0 <- mafi_result$signals$taf_k0

      # Ground truth consensus
      consensus <- sum(c(
        !is.na(egger_p) && egger_p < 0.10,
        !is.na(begg_p) && begg_p < 0.10,
        !is.na(taf_k0) && taf_k0 > 0
      ), na.rm = TRUE)

      # Correction accuracy
      original <- mafi_result$original
      corrected <- mafi_result$corrected
      orig_error <- abs(original - theta_true)
      corr_error <- abs(corrected - theta_true)

      results <- rbind(results, data.frame(
        k = k,
        bias = bias,
        sim = sim,
        mafi_score = mafi_result$score,
        egger_p = egger_p,
        begg_p = begg_p,
        taf_k0 = taf_k0,
        consensus = consensus,
        original = original,
        corrected = corrected,
        orig_error = orig_error,
        corr_error = corr_error,
        correction_improved = corr_error < orig_error
      ))
    }
  }
}

cat("\n")
cat("============================================================\n")
cat("  RESULTS: DETECTION PERFORMANCE\n")
cat("============================================================\n\n")

# Summarize by condition
summary_table <- aggregate(
  cbind(mafi_score, egger_sig = egger_p < 0.10, begg_sig = begg_p < 0.10,
        taf_detected = taf_k0 > 0) ~ k + bias,
  data = results,
  FUN = function(x) mean(x, na.rm = TRUE)
)

cat("Detection rates by condition:\n")
cat(sprintf("%-5s %-6s %-12s %-12s %-12s %-12s\n",
            "k", "Bias", "MAFI(mean)", "Egger(%)", "Begg(%)", "T&F(%)"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (i in 1:nrow(summary_table)) {
  cat(sprintf("%-5d %-6.1f %-12.1f %-12.1f %-12.1f %-12.1f\n",
              summary_table$k[i],
              summary_table$bias[i],
              summary_table$mafi_score[i],
              100 * summary_table$egger_sig[i],
              100 * summary_table$begg_sig[i],
              100 * summary_table$taf_detected[i]))
}

cat("\n")
cat("============================================================\n")
cat("  RESULTS: CORRECTION PERFORMANCE\n")
cat("============================================================\n\n")

correction_summary <- aggregate(
  cbind(orig_error, corr_error, correction_improved) ~ k + bias,
  data = results,
  FUN = function(x) mean(x, na.rm = TRUE)
)

cat("Correction accuracy (true theta = 0.30):\n")
cat(sprintf("%-5s %-6s %-12s %-12s %-15s\n",
            "k", "Bias", "Orig.Error", "Corr.Error", "Improved(%)"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (i in 1:nrow(correction_summary)) {
  cat(sprintf("%-5d %-6.1f %-12.3f %-12.3f %-15.1f\n",
              correction_summary$k[i],
              correction_summary$bias[i],
              correction_summary$orig_error[i],
              correction_summary$corr_error[i],
              100 * correction_summary$correction_improved[i]))
}

cat("\n")
cat("============================================================\n")
cat("  RESULTS: MAFI SCORE CALIBRATION\n")
cat("============================================================\n\n")

# Check if MAFI score correlates with actual bias
cat("MAFI score by true bias level:\n")
for (bias in bias_levels) {
  scores <- results$mafi_score[results$bias == bias]
  cat(sprintf("  Bias = %.1f: Mean MAFI = %.1f (SD = %.1f)\n",
              bias, mean(scores, na.rm = TRUE), sd(scores, na.rm = TRUE)))
}

# Correlation
cat(sprintf("\nCorrelation (MAFI score vs true bias): r = %.3f\n",
            cor(results$mafi_score, results$bias, use = "complete.obs")))

cat("\n")
cat("============================================================\n")
cat("  KEY FINDINGS\n")
cat("============================================================\n\n")

# No bias condition (false positive rate)
no_bias <- results[results$bias == 0, ]
fp_mafi <- mean(no_bias$mafi_score > 50, na.rm = TRUE)
fp_egger <- mean(no_bias$egger_p < 0.10, na.rm = TRUE)
fp_begg <- mean(no_bias$begg_p < 0.10, na.rm = TRUE)

cat("1. False Positive Rate (no true bias):\n")
cat(sprintf("   MAFI (score > 50): %.1f%%\n", 100 * fp_mafi))
cat(sprintf("   Egger (p < 0.10):  %.1f%%\n", 100 * fp_egger))
cat(sprintf("   Begg (p < 0.10):   %.1f%%\n", 100 * fp_begg))

# Strong bias condition (true positive rate)
strong_bias <- results[results$bias == 0.9, ]
tp_mafi <- mean(strong_bias$mafi_score > 50, na.rm = TRUE)
tp_egger <- mean(strong_bias$egger_p < 0.10, na.rm = TRUE)
tp_begg <- mean(strong_bias$begg_p < 0.10, na.rm = TRUE)

cat("\n2. True Positive Rate (strong bias = 0.9):\n")
cat(sprintf("   MAFI (score > 50): %.1f%%\n", 100 * tp_mafi))
cat(sprintf("   Egger (p < 0.10):  %.1f%%\n", 100 * tp_egger))
cat(sprintf("   Begg (p < 0.10):   %.1f%%\n", 100 * tp_begg))

# Correction improvement
strong_bias_corr <- mean(strong_bias$correction_improved, na.rm = TRUE)
cat(sprintf("\n3. Correction improved estimate in %.1f%% of strongly biased MAs\n",
            100 * strong_bias_corr))

cat("\n")
cat("============================================================\n")
cat("  SIMULATION STUDY COMPLETE\n")
cat("============================================================\n\n")
