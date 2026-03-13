# ==============================================================================
# Comprehensive Test Suite for cptac_correlation()
# Tests all 7 scenarios with real CPTAC data
# ==============================================================================

library(SLCPTAC)

cat("\n")
cat("================================================================================\n")
cat("COMPREHENSIVE TEST: cptac_correlation()\n")
cat("================================================================================\n")
cat("\n")

# Store all results
all_results <- list()

# ==============================================================================
# Scenario 1: 1 continuous vs 1 continuous (mRNA-Protein correlation)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 1: Scenario 1 - mRNA vs Protein correlation\n")
cat("Research Question: Does TP53 mRNA level correlate with protein level?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result1 <- cptac_correlation(
      var1 = "TP53",
      var1_modal = "RNAseq",
      var1_cancers = "BRCA",
      var2 = "TP53",
      var2_modal = "Protein",
      var2_cancers = "BRCA",
      method = "pearson"
    )

    runtime1 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 1 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime1))
    cat(sprintf("  Correlation: %.3f\n", result1$stats$r[1]))
    cat(sprintf("  P-value: %.2e\n", result1$stats$p[1]))
    cat(sprintf("  Sample size: %d\n", nrow(result1$raw_data)))
    cat(sprintf("  Plot type: %s\n", class(result1$plot)[1]))

    all_results$scenario1 <- list(
      success = TRUE,
      runtime = runtime1,
      stats = result1$stats,
      n_samples = nrow(result1$raw_data),
      plot_class = class(result1$plot)[1]
    )
  },
  error = function(e) {
    cat("\n✗ TEST 1 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario1 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 2: 1 continuous vs multiple continuous (Protein vs Phospho sites)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 2: Scenario 2 - Protein vs multiple Phospho sites\n")
cat("Research Question: Which AKT1 phospho sites correlate with protein level?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result2 <- cptac_correlation(
      var1 = "AKT1",
      var1_modal = "Protein",
      var1_cancers = "BRCA",
      var2 = "AKT1",
      var2_modal = "Phospho",
      var2_cancers = "BRCA",
      method = "pearson"
    )

    runtime2 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 2 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime2))
    cat(sprintf("  Number of phospho sites: %d\n", nrow(result2$stats)))
    cat(sprintf(
      "  Top correlation: %.3f (p=%.2e)\n",
      max(abs(result2$stats$r)),
      min(result2$stats$p)
    ))
    cat(sprintf("  Sample size: %d\n", nrow(result2$raw_data)))

    all_results$scenario2 <- list(
      success = TRUE,
      runtime = runtime2,
      n_phospho_sites = nrow(result2$stats),
      top_r = max(abs(result2$stats$r)),
      n_samples = nrow(result2$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 2 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario2 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 3: Multiple continuous vs multiple continuous (Phospho correlation matrix)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 3: Scenario 3 - Phospho correlation matrix\n")
cat("Research Question: How do AKT1 phospho sites correlate with each other?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result3 <- cptac_correlation(
      var1 = "AKT1",
      var1_modal = "Phospho",
      var1_cancers = "BRCA",
      var2 = "AKT1",
      var2_modal = "Phospho",
      var2_cancers = "BRCA",
      method = "pearson"
    )

    runtime3 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 3 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime3))
    cat(sprintf("  Number of correlations: %d\n", nrow(result3$stats)))
    cat(sprintf("  Mean correlation: %.3f\n", mean(abs(result3$stats$r))))
    cat(sprintf("  Sample size: %d\n", nrow(result3$raw_data)))

    all_results$scenario3 <- list(
      success = TRUE,
      runtime = runtime3,
      n_correlations = nrow(result3$stats),
      mean_r = mean(abs(result3$stats$r)),
      n_samples = nrow(result3$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 3 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario3 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 4: 1 categorical vs 1 continuous (Mutation impact on expression)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 4: Scenario 4 - Mutation impact on protein\n")
cat("Research Question: Does PIK3CA mutation affect AKT1 protein level?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result4 <- cptac_correlation(
      var1 = "PIK3CA",
      var1_modal = "Mutation",
      var1_cancers = "BRCA",
      var2 = "AKT1",
      var2_modal = "Protein",
      var2_cancers = "BRCA"
    )

    runtime4 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 4 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime4))
    cat(sprintf("  P-value: %.2e\n", result4$stats$p_value[1]))
    cat(sprintf("  Test method: %s\n", result4$stats$test_method[1]))
    cat(sprintf("  Effect size: %.3f\n", result4$stats$effect_size[1]))
    cat(sprintf("  Sample size: %d\n", nrow(result4$raw_data)))

    all_results$scenario4 <- list(
      success = TRUE,
      runtime = runtime4,
      pvalue = result4$stats$p_value[1],
      test_method = result4$stats$test_method[1],
      n_samples = nrow(result4$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 4 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario4 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 5: 1 continuous vs multiple categorical (Protein vs multiple mutations)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 5: Scenario 5 - Protein vs multiple mutations\n")
cat("Research Question: Which mutations affect EGFR protein level?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result5 <- cptac_correlation(
      var1 = "EGFR",
      var1_modal = "Protein",
      var1_cancers = "LUAD",
      var2 = c("KRAS", "EGFR", "TP53"),
      var2_modal = "Mutation",
      var2_cancers = "LUAD"
    )

    runtime5 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 5 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime5))
    cat(sprintf("  Number of tests: %d\n", nrow(result5$stats)))
    cat(sprintf("  Significant (p<0.05): %d\n", sum(result5$stats$p_value < 0.05)))
    cat(sprintf("  Sample size: %d\n", nrow(result5$raw_data)))

    all_results$scenario5 <- list(
      success = TRUE,
      runtime = runtime5,
      n_tests = nrow(result5$stats),
      n_significant = sum(result5$stats$p_value < 0.05),
      n_samples = nrow(result5$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 5 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario5 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 6: Multiple continuous vs 1 categorical (Multiple proteins vs mutation)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 6: Scenario 6 - Multiple proteins vs mutation\n")
cat("Research Question: How does TP53 mutation affect multiple pathway proteins?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result6 <- cptac_correlation(
      var1 = c("AKT1", "MTOR", "RPS6"),
      var1_modal = "Protein",
      var1_cancers = "BRCA",
      var2 = "TP53",
      var2_modal = "Mutation",
      var2_cancers = "BRCA"
    )

    runtime6 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 6 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime6))
    cat(sprintf("  Number of tests: %d\n", nrow(result6$stats)))
    cat(sprintf("  Significant (p<0.05): %d\n", sum(result6$stats$p_value < 0.05)))
    cat(sprintf("  Sample size: %d\n", nrow(result6$raw_data)))

    all_results$scenario6 <- list(
      success = TRUE,
      runtime = runtime6,
      n_tests = nrow(result6$stats),
      n_significant = sum(result6$stats$p_value < 0.05),
      n_samples = nrow(result6$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 6 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario6 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Scenario 7: Categorical vs categorical (Co-mutation analysis)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 7: Scenario 7 - Co-mutation analysis\n")
cat("Research Question: Which mutations co-occur or are mutually exclusive?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result7 <- cptac_correlation(
      var1 = c("KRAS", "EGFR"),
      var1_modal = "Mutation",
      var1_cancers = "LUAD",
      var2 = c("TP53", "STK11"),
      var2_modal = "Mutation",
      var2_cancers = "LUAD"
    )

    runtime7 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 7 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime7))
    cat(sprintf("  Number of pairs: %d\n", nrow(result7$stats)))
    cat(sprintf("  Significant (p<0.05): %d\n", sum(result7$stats$p_value < 0.05)))
    cat(sprintf("  Sample size: %d\n", nrow(result7$raw_data)))
    if ("log2_or" %in% colnames(result7$stats)) {
      cat(sprintf("  Mean log2(OR): %.3f\n", mean(result7$stats$log2_or, na.rm = TRUE)))
    }

    all_results$scenario7 <- list(
      success = TRUE,
      runtime = runtime7,
      n_pairs = nrow(result7$stats),
      n_significant = sum(result7$stats$p_value < 0.05),
      n_samples = nrow(result7$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 7 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$scenario7 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Multi-Cancer Test (Scenario 1 extended)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("BONUS TEST: Multi-cancer comparison\n")
cat("Research Question: Is mRNA-protein correlation consistent across cancers?\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result_multi <- cptac_correlation(
      var1 = "TP53",
      var1_modal = "RNAseq",
      var1_cancers = c("BRCA", "LUAD", "COAD"),
      var2 = "TP53",
      var2_modal = "Protein",
      var2_cancers = c("BRCA", "LUAD", "COAD"),
      method = "pearson"
    )

    runtime_multi <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ MULTI-CANCER TEST PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime_multi))
    cat(sprintf("  Number of cancer types: 3\n"))
    cat(sprintf("  Total comparisons: %d\n", nrow(result_multi$stats)))
    cat(sprintf("  Mean correlation: %.3f\n", mean(result_multi$stats$r)))
    cat(sprintf("  Sample size: %d\n", nrow(result_multi$raw_data)))

    all_results$multi_cancer <- list(
      success = TRUE,
      runtime = runtime_multi,
      n_comparisons = nrow(result_multi$stats),
      mean_r = mean(result_multi$stats$r),
      n_samples = nrow(result_multi$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ MULTI-CANCER TEST FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$multi_cancer <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Summary Report
# ==============================================================================
cat("\n")
cat("================================================================================\n")
cat("TEST SUMMARY\n")
cat("================================================================================\n")
cat("\n")

success_count <- sum(sapply(all_results, function(x) x$success))
total_count <- length(all_results)

cat(sprintf("Tests Passed: %d / %d\n\n", success_count, total_count))

if (success_count > 0) {
  cat("Performance Summary:\n")
  cat("-------------------\n")

  for (scenario_name in names(all_results)) {
    result <- all_results[[scenario_name]]
    if (result$success) {
      cat(sprintf(
        "%-20s: %.2f sec",
        toupper(gsub("_", " ", scenario_name)),
        result$runtime
      ))

      if (!is.null(result$n_samples)) {
        cat(sprintf(", %d samples", result$n_samples))
      }

      if (!is.null(result$n_phospho_sites)) {
        cat(sprintf(", %d phospho sites", result$n_phospho_sites))
      } else if (!is.null(result$n_correlations)) {
        cat(sprintf(", %d correlations", result$n_correlations))
      } else if (!is.null(result$n_tests)) {
        cat(sprintf(", %d tests", result$n_tests))
      } else if (!is.null(result$n_pairs)) {
        cat(sprintf(", %d pairs", result$n_pairs))
      }

      cat("\n")
    }
  }

  cat("\n")
  cat("All test results saved in 'all_results' object\n")
  cat("\n")
}

if (success_count < total_count) {
  cat("\nFailed Tests:\n")
  cat("-------------\n")
  for (scenario_name in names(all_results)) {
    result <- all_results[[scenario_name]]
    if (!result$success) {
      cat(sprintf(
        "%-20s: %s\n",
        toupper(gsub("_", " ", scenario_name)),
        result$error
      ))
    }
  }
}

cat("\n")
cat("================================================================================\n")
cat("TEST COMPLETED\n")
cat("================================================================================\n")
cat("\n")

# Return results for further analysis
invisible(all_results)
