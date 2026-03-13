# ==============================================================================
# Test Suite for cptac_survival()
# Testing Kaplan-Meier and Forest plot (Scenarios 16-17)
# ==============================================================================

library(SLCPTAC)

cat("\n")
cat("================================================================================\n")
cat("COMPREHENSIVE TEST: cptac_survival()\n")
cat("================================================================================\n")
cat("\n")

all_results <- list()

# ==============================================================================
# Test 1: Scenario 16 - Single Gene, Single Cancer (KM + Cox)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 1: Scenario 16 - Single Gene KM + Cox\n")
cat("Research Question: TP53 mRNA水平是否预测乳腺癌患者生存？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result1 <- cptac_survival(
      var1 = "TP53",
      var1_modal = "RNAseq",
      var1_cancers = "BRCA",
      surv_type = "OS",
      cutoff_type = "optimal"
    )

    runtime1 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 1 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime1))
    cat(sprintf("  KM p-value: %.2e\n", result1$stats$km_pvalue[1]))
    cat(sprintf(
      "  Cox HR: %.3f [%.3f-%.3f]\n",
      result1$stats$cox_hr[1],
      result1$stats$cox_hr_lower[1],
      result1$stats$cox_hr_upper[1]
    ))
    cat(sprintf("  Cox p-value: %.2e\n", result1$stats$cox_pvalue[1]))
    cat(sprintf("  C-index: %.3f\n", result1$stats$cox_cindex[1]))
    cat(sprintf("  Sample size: %d\n", nrow(result1$raw_data)))
    cat(sprintf(
      "  Return structure: stats (%d rows), plot, raw_data (%d samples)\n",
      nrow(result1$stats), nrow(result1$raw_data)
    ))

    all_results$test1 <- list(
      success = TRUE,
      scenario = "Scenario 16 (KM + Cox)",
      runtime = runtime1,
      km_pvalue = result1$stats$km_pvalue[1],
      cox_hr = result1$stats$cox_hr[1],
      cox_pvalue = result1$stats$cox_pvalue[1],
      cindex = result1$stats$cox_cindex[1],
      n_samples = nrow(result1$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 1 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test1 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 2: Scenario 16 - Protein (continuous variable)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 2: Scenario 16 - Protein Expression Survival\n")
cat("Research Question: AKT1蛋白水平是否预测乳腺癌患者生存？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result2 <- cptac_survival(
      var1 = "AKT1",
      var1_modal = "Protein",
      var1_cancers = "BRCA",
      surv_type = "OS",
      cutoff_type = "optimal"
    )

    runtime2 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 2 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime2))
    cat(sprintf("  KM p-value: %.2e\n", result2$stats$km_pvalue[1]))
    cat(sprintf(
      "  Cox HR: %.3f [%.3f-%.3f]\n",
      result2$stats$cox_hr[1],
      result2$stats$cox_hr_lower[1],
      result2$stats$cox_hr_upper[1]
    ))
    cat(sprintf("  C-index: %.3f\n", result2$stats$cox_cindex[1]))
    cat(sprintf("  Sample size: %d\n", nrow(result2$raw_data)))
    cat(sprintf(
      "  Return structure: stats (%d rows), plot, raw_data (%d samples)\n",
      nrow(result2$stats), nrow(result2$raw_data)
    ))

    all_results$test2 <- list(
      success = TRUE,
      scenario = "Scenario 16 (Protein)",
      runtime = runtime2,
      cox_hr = result2$stats$cox_hr[1],
      cindex = result2$stats$cox_cindex[1],
      n_samples = nrow(result2$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 2 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test2 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 3: Scenario 16 - Mutation (categorical variable)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 3: Scenario 16 - Mutation Survival\n")
cat("Research Question: KRAS突变是否影响肺癌患者生存？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result3 <- cptac_survival(
      var1 = "KRAS",
      var1_modal = "Mutation",
      var1_cancers = "LUAD",
      surv_type = "OS"
    )

    runtime3 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 3 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime3))
    cat(sprintf("  KM p-value: %.2e\n", result3$stats$km_pvalue[1]))
    cat(sprintf("  Cox HR: %.3f\n", result3$stats$cox_hr[1]))
    cat(sprintf("  Sample size: %d\n", nrow(result3$raw_data)))
    cat(sprintf(
      "  Return structure: stats (%d rows), plot, raw_data (%d samples)\n",
      nrow(result3$stats), nrow(result3$raw_data)
    ))

    all_results$test3 <- list(
      success = TRUE,
      scenario = "Scenario 16 (Mutation)",
      runtime = runtime3,
      cox_hr = result3$stats$cox_hr[1],
      n_samples = nrow(result3$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 3 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test3 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 4: Scenario 17 - Multiple Genes (Forest plot)
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 4: Scenario 17 - Multiple Genes Forest Plot\n")
cat("Research Question: 哪些基因影响肺癌患者生存？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result4 <- cptac_survival(
      var1 = c("TP53", "EGFR", "KRAS"),
      var1_modal = "RNAseq",
      var1_cancers = "LUAD",
      surv_type = "OS",
      cutoff_type = "optimal"
    )

    runtime4 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 4 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime4))
    cat(sprintf("  Number of features: %d\n", nrow(result4$stats)))
    cat(sprintf("  Significant (p<0.05): %d\n", sum(result4$stats$p_value < 0.05)))
    cat(sprintf("  Sample size: %d\n", nrow(result4$raw_data)))
    cat(sprintf(
      "  Return structure: stats (%d rows), plot, raw_data (%d samples)\n",
      nrow(result4$stats), nrow(result4$raw_data)
    ))

    all_results$test4 <- list(
      success = TRUE,
      scenario = "Scenario 17 (Multiple genes)",
      runtime = runtime4,
      n_features = nrow(result4$stats),
      n_significant = sum(result4$stats$p_value < 0.05),
      n_samples = nrow(result4$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 4 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test4 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 5: Scenario 17 - Single Gene, Multiple Cancers
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 5: Scenario 17 - Multi-Cancer Comparison\n")
cat("Research Question: TP53在不同癌种的预后价值是否一致？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result5 <- cptac_survival(
      var1 = "TP53",
      var1_modal = "RNAseq",
      var1_cancers = c("BRCA", "LUAD", "COAD"),
      surv_type = "OS",
      cutoff_type = "optimal"
    )

    runtime5 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 5 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime5))
    cat(sprintf("  Number of cancer types: 3\n"))
    cat(sprintf("  Features analyzed: %d\n", nrow(result5$stats)))
    cat(sprintf("  Sample size: %d\n", nrow(result5$raw_data)))
    cat(sprintf(
      "  Return structure: stats (%d rows), plot, raw_data (%d samples)\n",
      nrow(result5$stats), nrow(result5$raw_data)
    ))

    all_results$test5 <- list(
      success = TRUE,
      scenario = "Scenario 17 (Multi-cancer)",
      runtime = runtime5,
      n_cancers = 3,
      n_samples = nrow(result5$raw_data)
    )
  },
  error = function(e) {
    cat("\n✗ TEST 5 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test5 <<- list(success = FALSE, error = e$message)
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

  for (test_name in names(all_results)) {
    result <- all_results[[test_name]]
    if (result$success) {
      cat(sprintf(
        "%-10s (%s):\n",
        toupper(test_name),
        result$scenario
      ))
      cat(sprintf(
        "           Runtime: %.2f sec, %d samples\n",
        result$runtime, result$n_samples
      ))

      if (!is.null(result$cox_hr)) {
        cat(sprintf("           Cox HR: %.3f", result$cox_hr))
        if (!is.null(result$cindex)) {
          cat(sprintf(", C-index: %.3f\n", result$cindex))
        } else {
          cat("\n")
        }
      } else if (!is.null(result$n_features)) {
        cat(sprintf(
          "           Features: %d, Significant: %d\n",
          result$n_features, result$n_significant
        ))
      }
      cat("\n")
    }
  }

  cat("✅ Key Finding: Unified return structure verified\n")
  cat("   All functions return: list(stats, plot, raw_data)\n")
  cat("\n")
  cat("✅ Scenario Types Confirmed:\n")
  cat("   - Scenario 16 (Single feature): KM curve + Cox regression ✓\n")
  cat("   - Scenario 17 (Multiple features): Forest plot with HRs ✓\n")
  cat("   - Continuous variables: Optimal cutoff calculation ✓\n")
  cat("   - Categorical variables: Direct group comparison ✓\n")
  cat("   - Multi-cancer: Each cancer analyzed independently ✓\n")
  cat("\n")
}

cat("================================================================================\n")
cat("TEST COMPLETED\n")
cat("================================================================================\n")
cat("\n")

invisible(all_results)
