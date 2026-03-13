cat("=== SLCPTAC Modified Package Verification ===\n\n")

library(SLCPTAC)
passed <- 0
failed <- 0

run_test <- function(name, expr_fn) {
  cat(sprintf("Test: %s\n", name))
  tryCatch({
    expr_fn()
    passed <<- passed + 1
  }, error = function(e) {
    cat(sprintf("  FAIL: %s\n", e$message))
    failed <<- failed + 1
  })
}

# Test 1: No BioEnricher in namespace
run_test("No BioEnricher in namespace", function() {
  ns_imports <- getNamespaceImports("SLCPTAC")
  if ("BioEnricher" %in% names(ns_imports)) {
    stop("BioEnricher still in namespace imports!")
  }
  cat("  PASS: No BioEnricher dependency\n")
})

# Test 2: limma DEA with synthetic data
run_test("limma DEA with synthetic data", function() {
  set.seed(42)
  n_s <- 50; n_g <- 100
  expr <- matrix(rnorm(n_g * n_s), nrow = n_g, ncol = n_s)
  rownames(expr) <- paste0("Gene", 1:n_g)
  colnames(expr) <- paste0("S", 1:n_s)
  groups <- factor(rep(c("WildType", "Mutation"), each = 25),
                   levels = c("WildType", "Mutation"))
  expr[1:10, groups == "Mutation"] <- expr[1:10, groups == "Mutation"] + 2

  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  cm <- limma::makeContrasts(contrasts = "Mutation-WildType", levels = design)
  fit <- limma::lmFit(expr, design)
  fit <- limma::contrasts.fit(fit, cm)
  fit <- limma::eBayes(fit, trend = TRUE)
  res <- limma::topTable(fit, number = Inf, sort.by = "none")

  n_sig <- sum(res$adj.P.Val < 0.05)
  cat(sprintf("  PASS: %d genes, %d significant\n", nrow(res), n_sig))
  if (n_sig < 5) stop("Too few DE genes detected")
})

# Test 3: geneset package KEGG loading
run_test("geneset package KEGG loading", function() {
  gs <- SLCPTAC:::.fetch_raw_geneset(
    enrich_type = "KEGG", GO_ont = "BP", kegg_category = "pathway",
    msigdb_category = "H", hgdisease_source = "do", mesh_method = "gendoo",
    mesh_category = "A", enrichrdb_library = "Cancer_Cell_Line_Encyclopedia"
  )
  stopifnot(is.list(gs))
  stopifnot("geneset" %in% names(gs))
  stopifnot(nrow(gs$geneset) > 1000)
  cat(sprintf("  PASS: KEGG loaded, %d gene-pathway pairs\n", nrow(gs$geneset)))
})

# Test 4: Full geneset_df pipeline (MsigDB Hallmark)
run_test("Full geneset_df (MsigDB Hallmark)", function() {
  df <- SLCPTAC:::.load_geneset_df(enrich_type = "MsigDB", msigdb_category = "H")
  stopifnot(is.data.frame(df))
  stopifnot(all(c("id", "term", "gene") %in% colnames(df)))
  n_t <- length(unique(df$id))
  cat(sprintf("  PASS: %d terms, %d rows, sample genes: %s\n",
              n_t, nrow(df), paste(head(unique(df$gene), 3), collapse = ", ")))
  is_symbol <- !all(grepl("^[0-9]+$", head(unique(df$gene), 20)))
  if (!is_symbol) stop("Got Entrez IDs instead of symbols")
})

# Test 5: GO BP loading
run_test("GO BP gene set loading", function() {
  df <- SLCPTAC:::.load_geneset_df(enrich_type = "GO", GO_ont = "BP")
  n_t <- length(unique(df$id))
  cat(sprintf("  PASS: GO BP, %d terms, %d rows\n", n_t, nrow(df)))
})

# Test 6: Reactome loading
run_test("Reactome gene set loading", function() {
  df <- SLCPTAC:::.load_geneset_df(enrich_type = "Reactome")
  n_t <- length(unique(df$id))
  cat(sprintf("  PASS: Reactome, %d terms\n", n_t))
})

# Test 7: Full GSEA pipeline
run_test("Full GSEA pipeline (MsigDB H)", function() {
  df <- SLCPTAC:::.load_geneset_df(enrich_type = "MsigDB", msigdb_category = "H")
  all_genes <- unique(df$gene)
  set.seed(42)
  ranked <- sort(setNames(rnorm(length(all_genes)), all_genes), decreasing = TRUE)

  result <- SLCPTAC:::.perform_gsea(
    ranked_genes = ranked,
    enrich_type = "MsigDB",
    msigdb_category = "H",
    n_workers = 2,
    minSize = 10,
    maxSize = 500
  )
  stopifnot(is.data.frame(result))
  stopifnot(nrow(result) > 0)
  stopifnot(all(c("ID", "Description", "NES", "pvalue", "qvalue") %in% colnames(result)))
  cat(sprintf("  PASS: %d pathways, top: %s (NES=%.2f)\n",
              nrow(result), result$Description[1], result$NES[1]))
})

cat(sprintf("\n=== Results: %d passed, %d failed ===\n", passed, failed))
if (failed > 0) quit(status = 1)
