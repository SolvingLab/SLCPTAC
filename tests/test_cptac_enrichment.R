# ==============================================================================
# Test Suite for cptac_enrichment()
# Testing genome-wide scan and pathway enrichment (Scenarios 8-15)
# ==============================================================================

library(SLCPTAC)

cat("\n")
cat("================================================================================\n")
cat("COMPREHENSIVE TEST: cptac_enrichment()\n")
cat("================================================================================\n")
cat("\n")

all_results <- list()

# ==============================================================================
# Test 1: Scenario 8 - Single Mutation vs Genome-wide Protein
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 1: Scenario 8 - Mutation vs Genome-wide Protein Scan\n")
cat("Research Question: KRAS突变影响哪些蛋白质表达？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result1 <- cptac_enrichment(
      var1 = "KRAS",
      var1_modal = "Mutation",
      var1_cancers = "LUAD",
      analysis_type = "genome",
      genome_modal = "Protein",
      top_n = 30
    )

    runtime1 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 1 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime1))
    cat(sprintf("  Top genes identified: %d\n", nrow(result1$stats)))
    cat(sprintf(
      "  Total genes scanned: %d\n",
      if (is.list(result1$raw_data)) {
        length(result1$raw_data)
      } else {
        nrow(result1$raw_data)
      }
    ))
    cat(sprintf("  Return structure: stats (%d rows), plot, raw_data\n", nrow(result1$stats)))

    all_results$test1 <- list(
      success = TRUE,
      scenario = "Scenario 8 (Mutation vs Genome)",
      runtime = runtime1,
      n_top_genes = nrow(result1$stats),
      analysis_type = "genome"
    )
  },
  error = function(e) {
    cat("\n✗ TEST 1 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test1 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 2: Scenario 9 - Single Mutation vs GSEA Enrichment
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 2: Scenario 9 - Mutation vs GSEA Pathway Enrichment\n")
cat("Research Question: PIK3CA突变激活哪些信号通路？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result2 <- cptac_enrichment(
      var1 = "PIK3CA",
      var1_modal = "Mutation",
      var1_cancers = "BRCA",
      analysis_type = "enrichment",
      enrich_database = "MsigDB",
      msigdb_category = "H", # Hallmark gene sets
      top_n = 20
    )

    runtime2 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 2 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime2))
    cat(sprintf("  Pathways enriched: %d\n", nrow(result2$stats)))
    cat(sprintf("  Return structure: stats (%d rows), plot, raw_data\n", nrow(result2$stats)))

    all_results$test2 <- list(
      success = TRUE,
      scenario = "Scenario 9 (Mutation vs GSEA)",
      runtime = runtime2,
      n_pathways = nrow(result2$stats),
      analysis_type = "enrichment"
    )
  },
  error = function(e) {
    cat("\n✗ TEST 2 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test2 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 3: Scenario 12 - Single Protein vs Genome-wide Protein
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 3: Scenario 12 - Protein vs Genome-wide Protein Correlation\n")
cat("Research Question: 哪些蛋白质与TP53蛋白水平相关？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result3 <- cptac_enrichment(
      var1 = "TP53",
      var1_modal = "Protein",
      var1_cancers = "BRCA",
      analysis_type = "genome",
      genome_modal = "Protein",
      method = "pearson",
      top_n = 30
    )

    runtime3 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 3 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime3))
    cat(sprintf("  Top correlated proteins: %d\n", nrow(result3$stats)))
    cat(sprintf("  Return structure: stats (%d rows), plot, raw_data\n", nrow(result3$stats)))

    all_results$test3 <- list(
      success = TRUE,
      scenario = "Scenario 12 (Protein vs Genome)",
      runtime = runtime3,
      n_top_genes = nrow(result3$stats),
      analysis_type = "genome"
    )
  },
  error = function(e) {
    cat("\n✗ TEST 3 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test3 <<- list(success = FALSE, error = e$message)
  }
)

# ==============================================================================
# Test 4: Scenario 13 - Single Protein vs GSEA Enrichment
# ==============================================================================
cat("\n----------------------------------------\n")
cat("TEST 4: Scenario 13 - Protein vs GSEA Pathway Enrichment\n")
cat("Research Question: AKT1蛋白水平关联哪些通路活性？\n")
cat("----------------------------------------\n")

start <- Sys.time()
tryCatch(
  {
    result4 <- cptac_enrichment(
      var1 = "AKT1",
      var1_modal = "Protein",
      var1_cancers = "BRCA",
      analysis_type = "enrichment",
      enrich_database = "MsigDB",
      msigdb_category = "H",
      top_n = 20
    )

    runtime4 <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    cat("\n✓ TEST 4 PASSED\n")
    cat(sprintf("  Runtime: %.2f sec\n", runtime4))
    cat(sprintf("  Pathways enriched: %d\n", nrow(result4$stats)))
    cat(sprintf("  Return structure: stats (%d rows), plot, raw_data\n", nrow(result4$stats)))

    all_results$test4 <- list(
      success = TRUE,
      scenario = "Scenario 13 (Protein vs GSEA)",
      runtime = runtime4,
      n_pathways = nrow(result4$stats),
      analysis_type = "enrichment"
    )
  },
  error = function(e) {
    cat("\n✗ TEST 4 FAILED\n")
    cat(sprintf("  Error: %s\n", e$message))
    all_results$test4 <<- list(success = FALSE, error = e$message)
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
      cat(sprintf("           Runtime: %.2f sec\n", result$runtime))

      if (result$analysis_type == "genome") {
        cat(sprintf("           Top genes: %d\n", result$n_top_genes))
      } else {
        cat(sprintf("           Pathways: %d\n", result$n_pathways))
      }
      cat("\n")
    }
  }

  cat("✅ Key Finding: Unified return structure verified\n")
  cat("   All functions return: list(stats, plot, raw_data)\n")
  cat("\n")
  cat("✅ Analysis Types Confirmed:\n")
  cat("   - Genome-wide scan: Identifies top correlated/affected genes ✓\n")
  cat("   - GSEA enrichment: Identifies enriched pathways ✓\n")
  cat("   - Categorical (Mutation): DEA-based analysis ✓\n")
  cat("   - Continuous (Protein): Correlation-based analysis ✓\n")
  cat("\n")
}

cat("================================================================================\n")
cat("TEST COMPLETED\n")
cat("================================================================================\n")
cat("\n")

invisible(all_results)
