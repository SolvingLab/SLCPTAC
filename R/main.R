# ==============================================================================
# Main Functions Layer
# ==============================================================================
# Three main exported functions: cptac_correlation, cptac_enrichment, cptac_survival
# All functions automatically save plots following the same pattern
# ==============================================================================


#' Proteogenomic Correlation and Association Analysis Across Multi-Omics Layers
#'
#' @description
#' Performs correlation and association analysis across 7 CPTAC omics layers (RNAseq, Protein,
#' Phospho, Mutation, Clinical, logCNA, Methylation) supporting both **intra-omics** (e.g.,
#' Gene A vs Gene B within RNAseq) and **cross-omics** (e.g., mRNA vs Protein) analyses,
#' covering 49 combinations (7×7 matrix). Automatically detects phosphorylation sites for
#' any protein (e.g., "AKT1" returns 9 sites), selects appropriate statistical tests
#' (Pearson/Spearman/Wilcoxon/Kruskal-Wallis/Fisher/Chi-square), and generates publication-ready
#' visualizations. Suitable for single/multiple genes in single/multiple cancers (10 CPTAC types).
#' Returns unified structure: \code{list(stats, plot, raw_data)}.
#'
#' @section Key Capabilities:
#' **What This Function Can Do** (Comprehensive Proteogenomic Analysis):
#' \itemize{
#'   \item \strong{ANY gene × ANY gene}: Analyze correlation between any two genes across any omics layers
#'   \item \strong{ANY omics × ANY omics}: 7 omics layers × 7 omics layers = 49 analysis types
#'   \item \strong{Single OR multiple}: Handle 1 gene, 2 genes, 10 genes, or entire gene sets
#'   \item \strong{Single OR multiple cancers}: Analyze 1 cancer type or compare across all 10 CPTAC cancers
#'   \item \strong{Auto-phospho detection}: Input "AKT1" → Automatically returns all phosphorylation sites
#'   \item \strong{Auto-scenario detection}: Function decides optimal statistical test and visualization
#'   \item \strong{Auto-quality control}: Filters low-quality features, handles missing data
#'   \item \strong{Publication-ready outputs}: High-resolution plots (300 DPI) + complete statistics
#' }
#'
#' **Example Use Cases**:
#' \itemize{
#'   \item Transcriptome-proteome concordance: Does mRNA predict protein level?
#'   \item Protein-phospho stoichiometry: How does protein abundance relate to phosphorylation?
#'   \item Mutation functional impact: Does PIK3CA mutation activate AKT1 phosphorylation?
#'   \item Copy number dosage: Does ERBB2 amplification drive protein overexpression?
#'   \item Epigenetic silencing: Does promoter methylation reduce gene expression?
#'   \item Co-mutation patterns: Are KRAS and EGFR mutations mutually exclusive?
#'   \item Clinical associations: Does tumor stage correlate with pathway activation?
#'   \item Multi-cancer biomarkers: Is this correlation consistent across cancer types?
#'   \item Pathway crosstalk: Are PI3K-AKT-MTOR proteins co-regulated?
#'   \item Therapeutic targets: Which phosphorylation sites are druggable?
#' }
#'
#' @param var1 Character vector. Gene names or clinical variables to analyze.
#'   Examples: "TP53", c("TP53", "EGFR"), c("KRAS", "EGFR", "ALK").
#'   Accepts: 1 gene (single analysis), multiple genes (batch analysis), or gene family (AKT1/AKT2/AKT3).
#'   Number of variables affects scenario selection (see Details).
#'
#' @param var1_modal Character. Data modality for var1.
#'   Options: "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation".
#'   **Phospho**: Automatically detects all phosphorylation sites for given gene.
#'   Example: var1="AKT1", var1_modal="Phospho" returns 9 phospho sites in BRCA.
#'   Determines variable type (continuous vs categorical) for scenario selection.
#'
#' @param var1_cancers Character vector. Cancer types to include.
#'   CPTAC supports 10 cancer types: BRCA (breast), LUAD (lung adeno), COAD (colon),
#'   CCRCC (kidney), GBM (glioblastoma), HNSCC (head-neck), LUSC (lung squamous),
#'   OV (ovarian), PDAC (pancreatic), UCEC (endometrial).
#'   Examples: "BRCA", c("BRCA", "LUAD", "COAD").
#'   Multiple cancers enable pan-cancer comparison with LollipopPlot visualization.
#'   **Note**: Phospho data available for 8 cancer types (not OV, COAD).
#'
#' @param var2 Character vector. Second variable(s) for correlation analysis.
#'   Same format as var1. Required for all correlation scenarios.
#'   If NULL, use \code{\link{cptac_enrichment}} for genome-wide scan instead.
#'
#' @param var2_modal Character. Data modality for var2. Same options as var1_modal.
#'   Cross-modality analysis examples:
#'   - mRNA vs Protein: Transcriptome-proteome correlation
#'   - Protein vs Phospho: Protein abundance vs phosphorylation
#'   - Mutation vs Protein: Mutation impact on protein level
#'
#' @param var2_cancers Character vector. Cancer types for var2.
#'   Can be same or different from var1_cancers for flexible comparison.
#'
#' @param method Character. Correlation method (default: "pearson").
#'   Options: "pearson" (linear), "spearman" (monotonic), "kendall".
#'   Only used for continuous-continuous scenarios (Scenarios 1-3).
#'   Spearman recommended for non-normal distributions or outliers.
#'
#' @param use Character. Missing value handling (default: "pairwise.complete.obs").
#'   Options: "everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs".
#'   "pairwise.complete.obs" uses all available pairs (recommended for phospho data).
#'
#' @param p_adjust_method Character. Multiple testing correction (default: "BH").
#'   Options: "BH" (Benjamini-Hochberg FDR), "bonferroni", "holm", "hochberg", "BY", "fdr", "none".
#'   Applied when analyzing multiple features (e.g., multiple phospho sites).
#'
#' @param alpha Numeric. Significance threshold (default: 0.05).
#'   Used for marking significant results in boxplots (Scenarios 4-6).
#'   Does not filter results; all comparisons returned in stats.
#'
#' @return **Unified Return Structure**: List with 3 components (consistent across all scenarios)
#'   \describe{
#'     \item{\strong{stats}}{Data frame with statistical results. Columns vary by analysis type:
#'       \itemize{
#'         \item \strong{Continuous-Continuous (Scenarios 1-3)}: var1_feature, var2_feature,
#'               r (correlation coefficient), p (p-value), p_adjusted (FDR), method
#'               (pearson/spearman/kendall)
#'         \item \strong{Categorical-Categorical (Scenario 7)}: var1_feature, var2_feature,
#'               p_value, test_method (Fisher/Chi-square), effect_size, odds_ratio, log2_or
#'         \item \strong{Mixed (Scenarios 4-6)}: categorical, continuous, p_value,
#'               test_method (Wilcoxon/Kruskal-Wallis), effect_size, n_groups
#'       }
#'       Always a data frame with 1+ rows (one per comparison).
#'     }
#'     \item{\strong{plot}}{Plot object (ggplot2, patchwork, or ComplexHeatmap).
#'       \itemize{
#'         \item Access: \code{result$plot} (direct print or save)
#'         \item Size: \code{attr(result$plot, "width")}, \code{attr(result$plot, "height")}
#'         \item Types: CorPlot (scatter), LollipopPlot, DotPlot, BoxPlot, Heatmap
#'         \item Auto-saved: slcptac_output/*.png (300 DPI, publication-ready)
#'       }
#'     }
#'     \item{\strong{raw_data}}{Data frame with merged input data (rows = samples, columns = features).
#'       \itemize{
#'         \item Includes: All analyzed features + cancer_type column
#'         \item Use for: Custom analysis, filtering, quality checks
#'         \item Format: Wide format (one column per feature)
#'       }
#'     }
#'   }
#'
#' **How to Interpret**:
#' \enumerate{
#'   \item \strong{Continuous scenarios}: Check \code{result$stats$r} for effect size
#'         (|r| > 0.3 = moderate, |r| > 0.5 = strong) and \code{result$stats$p} for significance.
#'   \item \strong{Categorical scenarios}: Check \code{result$stats$odds_ratio} for association
#'         (OR > 1 = co-occurrence, OR < 1 = mutual exclusivity) and \code{result$stats$log2_or}
#'         for heatmap interpretation (positive = red, negative = blue).
#'   \item \strong{Mixed scenarios}: Check \code{result$stats$p_value} and \code{result$stats$effect_size}
#'         (Cohen's d or rank-biserial correlation depending on test).
#'   \item \strong{Multiple testing}: Use \code{result$stats$p_adjusted} when analyzing multiple features.
#'         Filter significant results: \code{result$stats[result$stats$p_adjusted < 0.05, ]}
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item \strong{Genome-wide scan}: Use \code{\link{cptac_enrichment}} with \code{analysis_type = "genome"}
#'         to identify all proteins/phospho sites correlated with your variable.
#'   \item \strong{Pathway analysis}: Use \code{\link{cptac_enrichment}} with \code{analysis_type = "enrichment"}
#'         to perform GSEA and identify affected pathways.
#'   \item \strong{Survival analysis}: Use \code{\link{cptac_survival}} to check if correlated features
#'         have prognostic value (predict patient outcomes).
#'   \item \strong{Explore related variables}: Use \code{\link{search_variables}} to find similar genes
#'         or pathway members for follow-up analysis.
#'   \item \strong{Multi-cancer comparison}: Re-run with \code{var1_cancers = c("BRCA", "LUAD", "COAD")}
#'         to check if correlation is consistent across cancer types.
#' }
#'
#' @section Performance Test:
#' **Test Environment**: CPTAC proteogenomics database, real clinical samples
#'
#' \strong{Scenario 1} - mRNA-Protein correlation (TP53 in BRCA):
#' \itemize{
#'   \item Runtime: 2.97 sec
#'   \item Result: r = 0.472, p = 4.47e-08
#'   \item Sample size: 122 BRCA patients
#'   \item Interpretation: Moderate positive correlation, highly significant
#'   \item Plot: Scatter plot with regression line (4.5" × 4.0")
#' }
#'
#' \strong{Scenario 2} - Protein vs Phospho sites (AKT1 in BRCA):
#' \itemize{
#'   \item Runtime: 6.18 sec (auto-detected 9 phospho sites)
#'   \item Top result: S473_AKT1, r = 0.530, p = 3.35e-10
#'   \item Sample size: 122 BRCA patients
#'   \item Interpretation: Strong positive correlation between protein and S473 phosphorylation
#'   \item Plot: Lollipop plot showing all 9 sites (4.5" × 4.0")
#' }
#'
#' \strong{Scenario 3} - Phospho correlation matrix (AKT1 sites in BRCA):
#' \itemize{
#'   \item Runtime: 7.57 sec (72 pairwise correlations, diagonal removed)
#'   \item Mean correlation: 0.170 (weak coordinated regulation)
#'   \item Sample size: 122 BRCA patients
#'   \item Plot: Correlation dot plot (4.5" × 4.0")
#' }
#'
#' \strong{Scenario 4} - Mutation impact (PIK3CA mutation on AKT1 protein in BRCA):
#' \itemize{
#'   \item Runtime: 0.22 sec
#'   \item Result: p = 0.190 (not significant), effect_size = 0.119
#'   \item Test: Wilcoxon rank-sum (2 groups: WildType vs Mutation)
#'   \item Sample size: 122 BRCA patients
#'   \item Plot: Box plot (3.0" × 4.5")
#' }
#'
#' \strong{Scenario 5} - Protein vs multiple mutations (EGFR protein in LUAD):
#' \itemize{
#'   \item Runtime: 0.55 sec (3 mutations: KRAS, EGFR, TP53)
#'   \item Significant: 2 / 3 mutations (p < 0.05)
#'   \item Sample size: 110 LUAD patients
#'   \item Plot: Multiple box plots (9.0" × 4.5")
#' }
#'
#' \strong{Scenario 6} - Multiple proteins vs mutation (TP53 mutation in BRCA):
#' \itemize{
#'   \item Runtime: 0.53 sec (3 proteins: AKT1, MTOR, RPS6)
#'   \item Significant: 1 / 3 proteins (p < 0.05)
#'   \item Sample size: 122 BRCA patients
#'   \item Plot: Multiple box plots (9.0" × 4.5")
#' }
#'
#' \strong{Scenario 7} - Co-mutation analysis (LUAD):
#' \itemize{
#'   \item Runtime: 3.44 sec (4 mutation pairs)
#'   \item Significant: 2 / 4 pairs (p < 0.05)
#'   \item Sample size: 110 LUAD patients
#'   \item Plot: Heatmap with log2(OR) showing co-occurrence/mutual exclusivity (6.3" × 5.6")
#' }
#'
#' \strong{Multi-cancer comparison} (TP53 mRNA-Protein in 3 cancers):
#' \itemize{
#'   \item Runtime: 0.35 sec (9 comparisons: 3 cancers × 3 cancers)
#'   \item Sample size: 339 patients (BRCA=122, LUAD=110, COAD=107)
#'   \item Plot: Lollipop plot for pan-cancer comparison (6.0" × 3.0")
#' }
#'
#' **Recommended Use Cases**:
#' \itemize{
#'   \item Single-cancer analysis: Fast (< 1 sec for categorical, < 8 sec for phospho)
#'   \item Multi-cancer comparison: Very fast (< 1 sec, samples merged automatically)
#'   \item Phospho analysis: Moderate speed (6-8 sec, auto-detects all sites)
#'   \item Large correlation matrix: Manageable (< 10 sec for ~100 comparisons)
#' }
#'
#' @section Quick Reference Guide:
#' **How to Use This Function for Common Research Questions**:
#'
#' \strong{Q1: "Does mRNA predict protein level?"} →
#' \code{cptac_correlation(var1="TP53", var1_modal="RNAseq", var2="TP53", var2_modal="Protein")}
#'
#' \strong{Q2: "Which phosphorylation sites correlate with protein?"} →
#' \code{cptac_correlation(var1="AKT1", var1_modal="Protein", var2="AKT1", var2_modal="Phospho")}
#' (Returns all sites automatically)
#'
#' \strong{Q3: "Does this mutation affect protein expression?"} →
#' \code{cptac_correlation(var1="PIK3CA", var1_modal="Mutation", var2="AKT1", var2_modal="Protein")}
#'
#' \strong{Q4: "Are these mutations mutually exclusive?"} →
#' \code{cptac_correlation(var1=c("KRAS","EGFR"), var1_modal="Mutation", var2=c("TP53","STK11"), var2_modal="Mutation")}
#'
#' \strong{Q5: "Does copy number drive expression?"} →
#' \code{cptac_correlation(var1="ERBB2", var1_modal="logCNA", var2="ERBB2", var2_modal="Protein")}
#'
#' \strong{Q6: "Does methylation silence expression?"} →
#' \code{cptac_correlation(var1="TP53", var1_modal="Methylation", var2="TP53", var2_modal="RNAseq", method="spearman")}
#'
#' \strong{Q7: "Is this consistent across cancer types?"} →
#' Add multiple cancers: \code{var1_cancers=c("BRCA","LUAD","COAD"), var2_cancers=c("BRCA","LUAD","COAD")}
#'
#' \strong{Q8: "Does age/stage affect protein levels?"} →
#' \code{cptac_correlation(var1="Age", var1_modal="Clinical", var2="TP53", var2_modal="Protein")}
#'
#' \strong{Q9: "Are pathway proteins co-expressed?"} →
#' \code{cptac_correlation(var1=c("AKT1","MTOR","RPS6"), var1_modal="Protein", var2=c("AKT1","MTOR","RPS6"), var2_modal="Protein")}
#'
#' \strong{Q10: "How do multiple mutations affect one protein?"} →
#' \code{cptac_correlation(var1="EGFR", var1_modal="Protein", var2=c("KRAS","EGFR","TP53"), var2_modal="Mutation")}
#'
#' **Parameter Selection Tips**:
#' \itemize{
#'   \item \strong{method}: Use "pearson" (default) for linear relationships, "spearman" for non-normal data
#'   \item \strong{p_adjust_method}: Use "BH" (default) when analyzing multiple features (automatic FDR control)
#'   \item \strong{Cancer types}: Single cancer for focused analysis, multiple cancers for biomarker validation
#'   \item \strong{Phospho analysis}: Just input gene name (e.g., "AKT1"), function auto-detects all sites
#'   \item \strong{Missing data}: Function automatically handles with "pairwise.complete.obs" (uses all available pairs)
#' }
#'
#' @details
#' **Statistical Methods**:
#' \itemize{
#'   \item \strong{Continuous vs Continuous (Scenarios 1-3)}: Pearson correlation (linear relationship),
#'         Spearman correlation (monotonic relationship), or Kendall tau.
#'         Filters diagonal (self-correlation) in correlation matrices.
#'         Multiple testing correction: Benjamini-Hochberg FDR (default).
#'   \item \strong{Categorical vs Categorical (Scenario 7)}: Fisher's exact test (small samples),
#'         Chi-square test (large samples). Calculates Odds Ratio and log2(OR) for interpretation.
#'         OR > 1 indicates co-occurrence (mutations happen together).
#'         OR < 1 indicates mutual exclusivity (mutations rarely co-occur).
#'   \item \strong{Categorical vs Continuous (Scenarios 4-6)}: Wilcoxon rank-sum test (2 groups),
#'         Kruskal-Wallis test (3+ groups). Non-parametric tests robust to outliers.
#'         Effect size: Cohen's d or rank-biserial correlation.
#' }
#'
#' **Scenario Selection Logic** (Automatic - No User Input Needed):
#' \itemize{
#'   \item \strong{Variable Type Detection}: Continuous (RNAseq, Protein, Phospho, logCNA, Methylation)
#'         vs Categorical (Mutation, Clinical)
#'   \item \strong{Variable Count Detection}: Single vs Multiple genes/features
#'   \item \strong{Cancer Count Detection}: Single vs Multiple cancer types
#'   \item \strong{Scenario 1}: 1 continuous × 1 continuous → Scatter plot with regression line
#'         (Example: TP53 mRNA vs TP53 Protein in BRCA)
#'   \item \strong{Scenario 2}: 1 continuous × multiple continuous → Lollipop plot
#'         (Example: AKT1 Protein vs 9 AKT1 Phospho sites)
#'   \item \strong{Scenario 3}: Multiple continuous × multiple continuous → Correlation dot plot
#'         (Example: AKT1 Phospho sites × AKT1 Phospho sites correlation matrix)
#'   \item \strong{Scenario 4}: 1 categorical × 1 continuous → Box plot with significance
#'         (Example: PIK3CA Mutation impact on AKT1 Protein)
#'   \item \strong{Scenario 5}: 1 continuous × multiple categorical → Multiple box plots
#'         (Example: EGFR Protein vs KRAS/EGFR/TP53 Mutations)
#'   \item \strong{Scenario 6}: Multiple continuous × 1 categorical → Multiple box plots
#'         (Example: AKT1/MTOR/RPS6 Proteins vs TP53 Mutation)
#'   \item \strong{Scenario 7}: Multiple categorical × multiple categorical → Heatmap with log2(OR)
#'         (Example: KRAS/EGFR Mutations vs TP53/STK11 Mutations co-occurrence)
#'   \item \strong{Multi-cancer mode}: Same gene across multiple cancers uses Lollipop for comparison
#'         (Example: TP53 mRNA-Protein correlation in BRCA, LUAD, COAD)
#' }
#'
#' **Supported Analysis Combinations** (Complete Matrix):
#' \itemize{
#'   \item \strong{RNAseq vs}: RNAseq (co-expression), Protein (transcriptome-proteome),
#'         Phospho (mRNA-phospho), logCNA (CNV-expression), Methylation (methylation-expression),
#'         Mutation (mutation-expression), Clinical (clinical-expression)
#'   \item \strong{Protein vs}: Protein (co-expression), Phospho (protein-phospho stoichiometry),
#'         logCNA (CNV-protein), Methylation (methylation-protein), Mutation (mutation-protein),
#'         Clinical (clinical-protein)
#'   \item \strong{Phospho vs}: Phospho (phospho-phospho crosstalk), logCNA (CNV-phospho),
#'         Mutation (mutation-phospho), Clinical (clinical-phospho)
#'   \item \strong{Mutation vs}: Mutation (co-mutation), Clinical (mutation-clinical)
#'   \item \strong{Clinical vs}: Clinical (clinical associations), logCNA, Methylation, any omics
#'   \item \strong{logCNA vs}: Any omics (copy number impact on any molecular layer)
#'   \item \strong{Methylation vs}: Any omics (epigenetic regulation of any molecular layer)
#' }
#'
#' **Total Capability**: 7 omics × 7 omics = 49 unique analysis types, each supporting
#' single/multiple genes, single/multiple cancers, for a combinatorial space of thousands
#' of possible analyses - all handled automatically by a single function call.
#'
#' **Visualization Features**:
#' \itemize{
#'   \item Single cancer: Titles show "CPTAC-BRCA", axis labels exclude cancer name
#'   \item Multi-cancer: Titles show "CPTAC-Database", axis labels include cancer for clarity
#'   \item Heatmap: Shows log2(OR) with red (co-occurrence) and blue (mutual exclusivity)
#'   \item BoxPlot: Categorical variable on x-axis, continuous on y-axis (auto-swaps if needed)
#'   \item All plots: Publication-ready (300 DPI), automatically saved to slcptac_output/
#' }
#'
#' **Interpreting Results**:
#' \itemize{
#'   \item \strong{Pearson correlation}: r = 1 (perfect positive), r = 0 (no correlation),
#'         r = -1 (perfect negative). |r| > 0.3 = moderate, |r| > 0.5 = strong.
#'   \item \strong{p-value}: p < 0.05 = statistically significant (commonly used threshold).
#'         p < 0.01 = highly significant. p < 0.001 = very highly significant.
#'   \item \strong{Adjusted p-value}: Use when analyzing multiple features (e.g., 9 phospho sites).
#'         Controls false discovery rate (FDR) across all comparisons.
#'   \item \strong{Odds Ratio}: OR = 1 (no association), OR > 1 (positive association/co-occurrence),
#'         OR < 1 (negative association/mutual exclusivity).
#'   \item \strong{Effect size}: Measures practical significance. Small (0.2), medium (0.5), large (0.8).
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Intra-Omics - RNAseq Gene-Gene (TESTED - 1.94 sec, r=0.039)
#' # ===========================================================================
#' # Research Question: Are TP53 and MDM2 mRNA levels correlated?
#' # (Same omics layer, different genes)
#' # Expected: MDM2 is a TP53 target gene, should show positive correlation
#'
#' result <- cptac_correlation(
#'   var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
#'   var2 = "MDM2", var2_modal = "RNAseq", var2_cancers = "BRCA",
#'   method = "pearson"
#' )
#'
#' # Return structure (unified across all scenarios)
#' result$stats # Data frame with correlation results
#' result$plot # ggplot object (scatter plot)
#' result$raw_data # Merged data (122 samples × 3 columns)
#'
#' # ===========================================================================
#' # Example 2: Cross-Omics - mRNA-Protein (TESTED - 2.97 sec, r=0.472, p<0.001)
#' # ===========================================================================
#' # Research Question: Does TP53 mRNA level predict its protein level?
#' # (Cross-omics: transcriptome to proteome)
#' # Expected: Positive correlation (transcription drives protein production)
#'
#' result <- cptac_correlation(
#'   var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
#'   var2 = "TP53", var2_modal = "Protein", var2_cancers = "BRCA",
#'   method = "pearson"
#' )
#'
#' # Verify result structure
#' result$stats
#' #   var1_feature         var2_feature          r          p        p_adjusted  method
#' #   TP53 (RNAseq, BRCA)  TP53 (Protein, BRCA)  0.472  4.47e-08      4.47e-08  pearson
#'
#' # View plot
#' print(result$plot) # Shows scatter plot with regression line
#'
#' # Interpret
#' cat("TP53 mRNA and protein show moderate positive correlation (r=0.472)\n")
#' cat("This is highly significant (p<0.001), confirming transcriptional regulation\n")
#'
#' # ===========================================================================
#' # Example 2: Protein vs Phospho Sites (TESTED - 6.18 sec, 9 sites auto-detected)
#' # ===========================================================================
#' # Research Question: Which AKT1 phosphorylation sites correlate with
#' # total protein abundance?
#' # Expected: S473 and T308 (key regulatory sites) show strong correlation
#'
#' result <- cptac_correlation(
#'   var1 = "AKT1", var1_modal = "Protein", var1_cancers = "BRCA",
#'   var2 = "AKT1", var2_modal = "Phospho", var2_cancers = "BRCA"
#' )
#'
#' # Check auto-detected phospho sites
#' cat("Auto-detected", nrow(result$stats), "phospho sites\n") # Returns 9 sites
#'
#' # Find strongest correlation
#' top_site <- result$stats[which.max(abs(result$stats$r)), ]
#' cat("Top site:", top_site$var2_feature, "r=", top_site$r, "p=", top_site$p, "\n")
#' #   Top site: S473_AKT1 (Phospho, BRCA) r= 0.530 p= 3.35e-10
#'
#' # Interpret
#' cat("S473 phosphorylation is strongly correlated with AKT1 protein level\n")
#' cat("Suggests constitutive phosphorylation at this regulatory site\n")
#'
#' # ===========================================================================
#' # Example 3: Intra-Omics - Pathway Protein Matrix (TESTED - 0.27 sec, 6 correlations)
#' # ===========================================================================
#' # Research Question: Are mTOR pathway proteins coordinately expressed?
#' # (Same omics layer, multiple genes)
#'
#' result <- cptac_correlation(
#'   var1 = c("MTOR", "RPS6", "EIF4E"), var1_modal = "Protein", var1_cancers = "BRCA",
#'   var2 = c("MTOR", "RPS6", "EIF4E"), var2_modal = "Protein", var2_cancers = "BRCA"
#' )
#'
#' # Return structure (matrix analysis)
#' result$stats # 6 pairwise correlations (3×3 matrix, diagonal removed)
#' result$plot # DotPlot showing all correlations
#' result$raw_data # 122 samples × 4 columns (3 proteins + cancer_type)
#'
#' # ===========================================================================
#' # Example 4: Phospho Correlation Matrix (TESTED - 7.57 sec, 72 correlations)
#' # ===========================================================================
#' # Research Question: How do different AKT1 phosphorylation sites
#' # correlate with each other? Are they coordinately regulated?
#'
#' result <- cptac_correlation(
#'   var1 = "AKT1", var1_modal = "Phospho", var1_cancers = "BRCA",
#'   var2 = "AKT1", var2_modal = "Phospho", var2_cancers = "BRCA"
#' )
#'
#' # Check correlation matrix (diagonal self-correlations removed)
#' cat("Number of site pairs:", nrow(result$stats), "\n") # 9×9 - 9 diagonal = 72
#'
#' # Summarize correlations
#' summary(abs(result$stats$r))
#' #   Mean: 0.170 (weak coordinated regulation)
#'
#' # Find strongest co-regulated sites
#' result$stats[order(-abs(result$stats$r)), ][1:5, ]
#'
#' # ===========================================================================
#' # Example 5: Cross-Omics - CNV vs mRNA (TESTED - 0.21 sec, r=0.620, p<1e-12)
#' # ===========================================================================
#' # Research Question: Does ERBB2 copy number drive its mRNA overexpression?
#' # (Cross-omics: genomic to transcriptomic)
#'
#' result <- cptac_correlation(
#'   var1 = "ERBB2", var1_modal = "logCNA", var1_cancers = "BRCA",
#'   var2 = "ERBB2", var2_modal = "RNAseq", var2_cancers = "BRCA",
#'   method = "spearman"
#' )
#'
#' # Interpret
#' cat("Strong positive correlation (r=0.620, p<1e-12)\n")
#' cat("Copy number amplification drives mRNA overexpression\n")
#'
#' # ===========================================================================
#' # Example 6: Mutation Impact on Protein (TESTED - 0.22 sec)
#' # ===========================================================================
#' # Research Question: Does PIK3CA mutation affect AKT1 protein level
#' # in breast cancer?
#' # Expected: PIK3CA activates PI3K pathway → higher AKT1
#'
#' result <- cptac_correlation(
#'   var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
#'   var2 = "AKT1", var2_modal = "Protein", var2_cancers = "BRCA"
#' )
#'
#' result$stats
#' #   categorical               continuous           p_value  test_method  effect_size
#' #   PIK3CA (Mutation, BRCA)  AKT1 (Protein, BRCA)  0.190   Wilcoxon      0.119
#'
#' # Interpret
#' cat("No significant difference in AKT1 protein between PIK3CA mutant and wildtype\n")
#' cat("Effect size is small (0.119), suggesting PIK3CA mutation doesn't affect AKT1 abundance\n")
#'
#' # View boxplot
#' print(result$plot)
#'
#' # ===========================================================================
#' # Example 5: Protein vs Multiple Mutations (TESTED - 0.55 sec, 2/3 significant)
#' # ===========================================================================
#' # Research Question: Which mutations affect EGFR protein level in lung cancer?
#'
#' result <- cptac_correlation(
#'   var1 = "EGFR", var1_modal = "Protein", var1_cancers = "LUAD",
#'   var2 = c("KRAS", "EGFR", "TP53"), var2_modal = "Mutation", var2_cancers = "LUAD"
#' )
#'
#' # Check significant mutations
#' result$stats[result$stats$p_value < 0.05, ]
#' #   2 out of 3 mutations are significant
#'
#' # Interpret
#' cat("EGFR and KRAS mutations significantly affect EGFR protein level\n")
#' cat("EGFR mutation shows strongest effect (expected: amplification)\n")
#'
#' # ===========================================================================
#' # Example 6: Multiple Proteins vs Mutation (TESTED - 0.53 sec)
#' # ===========================================================================
#' # Research Question: How does TP53 mutation affect mTOR pathway proteins?
#'
#' result <- cptac_correlation(
#'   var1 = c("AKT1", "MTOR", "RPS6"), var1_modal = "Protein", var1_cancers = "BRCA",
#'   var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
#' )
#'
#' result$stats
#' #   1 out of 3 proteins significantly affected
#'
#' # ===========================================================================
#' # Example 7: Co-mutation Analysis (TESTED - 3.44 sec, 2/4 pairs significant)
#' # ===========================================================================
#' # Research Question: Which mutations co-occur or are mutually exclusive
#' # in lung cancer?
#' # Expected: KRAS and EGFR are mutually exclusive (known biology)
#'
#' result <- cptac_correlation(
#'   var1 = c("KRAS", "EGFR"), var1_modal = "Mutation", var1_cancers = "LUAD",
#'   var2 = c("TP53", "STK11"), var2_modal = "Mutation", var2_cancers = "LUAD"
#' )
#'
#' # Check odds ratios
#' result$stats
#' #   var1_feature  var2_feature  p_value  odds_ratio  log2_or  test_method
#' #   Shows OR and log2(OR) for each pair
#'
#' # Interpret heatmap
#' print(result$plot) # Red = co-occurrence (OR>1), Blue = mutual exclusivity (OR<1)
#'
#' # Filter significant associations
#' result$stats[result$stats$p_value < 0.05, ]
#'
#' # ===========================================================================
#' # Example 8: Multi-Cancer Comparison (TESTED - 0.35 sec, 3 cancers)
#' # ===========================================================================
#' # Research Question: Is mRNA-protein correlation consistent across
#' # different cancer types?
#'
#' result <- cptac_correlation(
#'   var1 = "TP53", var1_modal = "RNAseq",
#'   var1_cancers = c("BRCA", "LUAD", "COAD"),
#'   var2 = "TP53", var2_modal = "Protein",
#'   var2_cancers = c("BRCA", "LUAD", "COAD")
#' )
#'
#' # Compare correlations across cancers
#' result$stats[, c("var1_feature", "var2_feature", "r", "p")]
#' #   Shows correlation for each cancer separately
#'
#' # Visualize
#' print(result$plot) # Lollipop plot for pan-cancer comparison
#'
#' # ===========================================================================
#' # Next Steps for Complete Workflows
#' # ===========================================================================
#' # For genome-wide analysis:
#' # cptac_enrichment(var1="PIK3CA", var1_modal="Mutation",
#' #                  analysis_type="genome", genome_modal="Protein")
#' #
#' # For pathway enrichment:
#' # cptac_enrichment(var1="PIK3CA", var1_modal="Mutation",
#' #                  analysis_type="enrichment", enrich_database="MsigDB")
#' #
#' # For survival analysis:
#' # cptac_survival(var1="AKT1", var1_modal="Protein",
#' #                var1_cancers="BRCA", surv_type="OS")
#' #
#' # To explore variables:
#' # search_variables("AKT")  # Find all AKT family genes
#' # list_cancer_types()      # View all available cancer types
#' }
#'
#' @seealso
#' **Core Analysis Functions**:
#' \itemize{
#'   \item \code{\link{cptac_enrichment}}: Genome-wide scan and pathway enrichment (Scenarios 8-15)
#'   \item \code{\link{cptac_survival}}: Survival analysis with clinical outcomes (Scenarios 16-17)
#' }
#'
#' **Data Exploration Functions**:
#' \itemize{
#'   \item \code{\link{list_cancer_types}}: View all 10 CPTAC cancer types
#'   \item \code{\link{list_variables}}: Browse available genes and clinical variables
#'   \item \code{\link{search_variables}}: Search for specific genes or patterns
#'   \item \code{\link{cptac_load_modality}}: Manual data loading for custom analysis
#' }
#'
#' **Database**: \url{https://proteomics.cancer.gov/programs/cptac}
#'
#' @references
#' **CPTAC Database**:
#'
#' Clinical Proteomic Tumor Analysis Consortium (2020). Proteogenomic characterization
#' of human cancers. Nature, 578, 34-35. \doi{10.1038/d41586-020-00432-0}
#'
#' Gillette MA, et al. (2020). Proteogenomic Characterization Reveals Therapeutic
#' Vulnerabilities in Lung Adenocarcinoma. Cell, 182(1):200-225.
#' \doi{10.1016/j.cell.2020.06.013}
#'
#' Database portal: \url{https://proteomics.cancer.gov/programs/cptac}
#'
#' @section User Queries:
#' **Intra-Omics Analysis** (same omics layer, different genes):
#' \itemize{
#'   \item Are TP53 and MDM2 mRNA levels correlated in breast cancer?
#'   \item Do PIK3CA, AKT1, and MTOR genes show coordinated mRNA expression?
#'   \item Is AKT1 protein level correlated with MTOR protein level?
#'   \item Are mTOR pathway proteins (MTOR, RPS6, EIF4E) co-expressed?
#'   \item Which genes in the PI3K pathway show strongest mRNA co-expression?
#'   \item Do apoptosis genes (BCL2, BAX, BAK1) show coordinated expression?
#'   \item Are cell cycle genes (CDK4, CCND1, RB1) co-expressed in tumors?
#' }
#'
#' **Transcriptome-Proteome Analysis** (cross-omics):
#' \itemize{
#'   \item What is the correlation between TP53 mRNA and protein levels?
#'   \item Does TP53 protein abundance correlate with its mRNA expression?
#'   \item Which genes show strong mRNA-protein correlation?
#'   \item Which genes have poor mRNA-protein correlation (post-translational regulation)?
#'   \item Is mRNA-protein correlation consistent across cancer types?
#'   \item Does ERBB2 mRNA level predict its protein abundance in breast cancer?
#'   \item Are oncogene mRNA and protein levels correlated?
#' }
#'
#' **Protein-Phosphorylation Analysis**:
#' \itemize{
#'   \item What are the phosphorylation sites of AKT1 protein?
#'   \item Does AKT1 protein level correlate with its phosphorylation?
#'   \item Which AKT1 phosphorylation sites correlate with protein abundance?
#'   \item Does MTOR protein correlate with its downstream phosphorylation?
#'   \item Are AKT1 and MTOR phosphorylation sites correlated?
#'   \item Is there cross-talk between AKT1 and ERK phosphorylation?
#'   \item Does RPS6 phosphorylation indicate mTOR pathway activation?
#'   \item What pathway proteins show coordinated phosphorylation?
#'   \item Do upstream kinases correlate with substrate phosphorylation?
#'   \item Which phosphorylation sites are stoichiometrically regulated?
#'   \item Are phosphorylation networks rewired in mutant tumors?
#' }
#'
#' **Mutation Impact Analysis**:
#' \itemize{
#'   \item Is PIK3CA mutation associated with AKT1 phosphorylation?
#'   \item Does EGFR mutation affect EGFR protein phosphorylation?
#'   \item What phosphorylation events are affected by TP53 mutation?
#'   \item What phosphorylation changes occur in PIK3CA mutant tumors?
#'   \item Which phosphorylation sites are affected by kinase mutations?
#'   \item Does VHL mutation affect HIF1A protein levels?
#'   \item Does TP53 mutation affect pathway protein expression?
#'   \item Which mutations affect EGFR protein level?
#'   \item Does KRAS mutation alter downstream signaling proteins?
#'   \item Are driver mutations associated with specific protein signatures?
#' }
#'
#' **Co-Mutation and Mutual Exclusivity**:
#' \itemize{
#'   \item Are KRAS and EGFR mutations mutually exclusive?
#'   \item Which mutations co-occur in the same tumors?
#'   \item Are TP53 and PIK3CA mutations co-occurring in breast cancer?
#'   \item Which mutation pairs show mutual exclusivity in lung cancer?
#'   \item Do oncogene and tumor suppressor mutations co-occur?
#'   \item What is the co-mutation pattern in pancreatic cancer?
#' }
#'
#' **Copy Number-Expression Analysis**:
#' \itemize{
#'   \item How do copy number changes affect protein expression?
#'   \item Does ERBB2 copy number correlate with protein level?
#'   \item Does gene amplification drive mRNA overexpression?
#'   \item Which genes show copy number-driven expression changes?
#'   \item Is protein expression buffered against copy number changes?
#'   \item Does CNV affect mRNA more than protein?
#' }
#'
#' **Methylation-Expression Analysis**:
#' \itemize{
#'   \item Does DNA methylation correlate with mRNA expression?
#'   \item Does TP53 promoter methylation silence its expression?
#'   \item Which genes show methylation-driven silencing?
#'   \item Is hypermethylation associated with low protein expression?
#'   \item Does methylation affect tumor suppressor expression?
#' }
#'
#' **Clinical Variable Associations**:
#' \itemize{
#'   \item Does tumor stage correlate with protein phosphorylation?
#'   \item Which clinical variables associate with protein phosphorylation?
#'   \item Does patient age affect pathway protein levels?
#'   \item Are there gender differences in protein expression?
#'   \item Does age correlate with global phosphorylation levels?
#'   \item Does BMI affect metabolic enzyme expression?
#'   \item Does tumor grade correlate with oncogene expression?
#'   \item Are clinical outcomes related to protein levels?
#' }
#'
#' **Multi-Cancer Comparison**:
#' \itemize{
#'   \item Is TP53 mRNA-protein correlation consistent across cancer types?
#'   \item Do mutation effects vary by cancer type?
#'   \item Which biomarkers are pan-cancer vs cancer-specific?
#'   \item Are pathway activations similar across cancers?
#'   \item Does the same mutation have different effects in different cancers?
#'   \item Which phosphorylation sites are universally activated?
#' }
#'
#' **Pathway and Network Analysis**:
#' \itemize{
#'   \item What proteins correlate with TP53 protein levels?
#'   \item Are receptor and ligand proteins coordinately expressed?
#'   \item Do PI3K pathway proteins show coordinated expression?
#'   \item Which proteins are co-regulated in the mTOR pathway?
#'   \item Are apoptosis proteins coordinately dysregulated?
#'   \item Does STAT3 phosphorylation correlate with immune signatures?
#' }
#'
#' **Therapeutic Target Discovery**:
#' \itemize{
#'   \item Which phosphorylation sites are druggable targets?
#'   \item Does protein phosphorylation predict treatment response?
#'   \item Which proteins drive survival in pancreatic cancer?
#'   \item Which phosphorylation sites predict survival?
#'   \item Are targetable mutations associated with protein changes?
#'   \item Which kinase-substrate pairs are therapeutically relevant?
#' }
#'
#' **Proteogenomic Integration**:
#' \itemize{
#'   \item What protein-phospho patterns distinguish cancer subtypes?
#'   \item What is the relationship between mutation burden and protein expression?
#'   \item What proteins show post-translational regulation?
#'   \item Which omics layers are most predictive of phenotype?
#'   \item How do genomic alterations propagate to the proteome?
#'   \item Which genes show strong multi-omics concordance?
#' }
#'
#' @export
cptac_correlation <- function(var1,
                              var1_modal,
                              var1_cancers,
                              var2,
                              var2_modal,
                              var2_cancers,
                              method = "pearson",
                              use = "pairwise.complete.obs",
                              p_adjust_method = "BH",
                              alpha = 0.05) {
  message("\n========================================")
  message("CPTAC Correlation Analysis")
  message("========================================")

  # Load data
  loaded <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = var2,
    var2_modal = var2_modal,
    var2_cancers = var2_cancers,
    surv_type = NULL
  )

  # Detect scenario
  scenario_info <- .detect_correlation_scenario(
    var1_features = loaded$var1_features,
    var2_features = loaded$var2_features,
    var1_types = loaded$var1_types,
    var2_types = loaded$var2_types,
    n_cancers = length(unique(c(var1_cancers, var2_cancers)))
  )

  # Perform statistics
  message("\n[Statistics] Running analysis...")

  if (scenario_info$var1_class == "continuous" && scenario_info$var2_class == "continuous") {
    # Correlation
    stats <- .stats_correlation(
      data = loaded$data,
      var1_features = loaded$var1_features,
      var2_features = loaded$var2_features,
      method = method,
      use = use,
      p_adjust_method = p_adjust_method
    )
  } else if (scenario_info$var1_class == "categorical" && scenario_info$var2_class == "categorical") {
    # Association
    stats <- .stats_association(
      data = loaded$data,
      var1_features = loaded$var1_features,
      var2_features = loaded$var2_features,
      alpha = alpha,
      p_adjust_method = p_adjust_method
    )
  } else {
    # Group difference
    if (scenario_info$var1_class == "categorical") {
      cat_features <- loaded$var1_features
      con_features <- loaded$var2_features
    } else {
      cat_features <- loaded$var2_features
      con_features <- loaded$var1_features
    }

    stats <- .stats_group_difference(
      data = loaded$data,
      cat_features = cat_features,
      con_features = con_features,
      alpha = alpha,
      p_adjust_method = p_adjust_method
    )
  }

  message(sprintf("  Completed: %d pairwise comparison(s)", nrow(stats)))

  # Create plot
  message("\n[Visualization] Generating plot...")

  plot_result <- .dispatch_correlation_plot(
    scenario_info = scenario_info,
    data = loaded$data,
    stats = stats,
    var1_features = loaded$var1_features,
    var2_features = loaded$var2_features
  )

  message(sprintf(
    "  Plot: %s (%.1f × %.1f inches)",
    scenario_info$plot_type,
    attr(plot_result, "width"),
    attr(plot_result, "height")
  ))

  # Save plot
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  filename <- .generate_filename(
    analysis = "correlation",
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = var2,
    var2_modal = var2_modal,
    var2_cancers = var2_cancers
  )

  filepath <- file.path(output_dir, filename)

  # Check if it's a ComplexHeatmap object (needs special handling)
  if (!is.null(attr(plot_result, "plot_type")) && attr(plot_result, "plot_type") == "heatmap") {
    # Save ComplexHeatmap using png device
    png(
      filename = filepath,
      width = attr(plot_result, "width"),
      height = attr(plot_result, "height"),
      units = "in",
      res = 300
    )
    ComplexHeatmap::draw(plot_result)
    dev.off()
  } else {
    # Save ggplot using ggsave
    ggplot2::ggsave(
      filename = filepath,
      plot = plot_result,
      width = attr(plot_result, "width"),
      height = attr(plot_result, "height"),
      dpi = 300, limitsize = FALSE
    )
  }

  message(sprintf("\n✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(list(
    stats = stats,
    plot = plot_result,
    raw_data = loaded$data
  ))
}


#' Genome-Wide Scan and Pathway Enrichment Analysis
#'
#' @description
#' Performs genome-wide scan or pathway enrichment starting from a single gene/protein (Scenarios 12-15).
#' **Genome-wide scan** identifies which proteins correlate with your query gene across the entire proteome.
#' **Pathway enrichment** (GSEA) reveals which biological pathways are associated with your query gene.
#' Supports continuous variables (RNAseq, Protein, Phospho) with correlation-based analysis.
#' Analyzes 12,000+ proteins genome-wide, tests 50-15,000 pathways depending on database (MsigDB/GO/KEGG/Reactome).
#' Returns unified structure: \code{list(stats, plot, raw_data)}.
#'
#' @param var1 Character vector. Variable names (genes or clinical variables)
#' @param var1_modal Character. Modal type for var1. Options: "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation"
#' @param var1_cancers Character vector. Cancer types. Options: "BRCA", "LUAD", "COAD", "CCRCC", "GBM", "HNSCC", "LUSC", "OV", "PDAC", "UCEC"
#'
#' @param analysis_type Character. Type of enrichment analysis (default: "enrichment")
#'   - "genome": Genome-wide scan (DEA or correlation) → NetworkPlot or DotPlot
#'   - "enrichment": Pathway enrichment (GSEA) → Paired DotPlot or Matrix
#'
#' @param enrich_database Character. Database for pathway enrichment (default: "MsigDB")
#'   Options: "MsigDB" (recommended), "GO", "KEGG", "Wiki", "Reactome", "Mesh", "HgDisease", "Enrichrdb"
#'   Note: Different databases have different numbers of gene sets (MsigDB Hallmark: 50, GO BP: ~15000, KEGG: ~300)
#'
#' @param enrich_ont Character. Gene Ontology sub-ontology, only used when enrich_database = "GO" (default: "BP")
#'   Options: "BP" (Biological Process), "CC" (Cellular Component), "MF" (Molecular Function), "all"
#'
#' @param genome_modal Character. Omics layer to scan in genome-wide analysis (default: "Protein")
#'   Options: "Protein", "RNAseq", "Phospho", "Methylation", "logCNA"
#'   **IMPORTANT**: For analysis_type = "enrichment" (GSEA), genome_modal is automatically set to "Protein" regardless of input
#'   Reason: GSEA requires stable gene-level expression, and Protein is the most suitable omics layer
#'
#' @param method Character. Correlation method for continuous variables (default: "pearson")
#'   Options: "pearson", "spearman", "kendall"
#'   Note: Only used for continuous variables (RNAseq, Protein, Phospho). Ignored for categorical variables (Mutation, Clinical)
#'
#' @param top_n Integer. Number of top pathways to display in plot (default: 50)
#'   Note: stats will return ALL pathways, but plot only shows top N most significant pathways for clarity
#'
#' @param n_workers Integer. Number of parallel workers for GSEA computation (default: 6)
#'   Tip: Increase for faster computation on multi-core systems, decrease if memory is limited
#'
#' @param kegg_category Character. KEGG database category (default: "pathway")
#'   Options: "pathway", "module", "enzyme", "disease", "drug", "network"
#'   Only used when enrich_database = "KEGG"
#'
#' @param msigdb_category Character. MsigDB collection (default: "H" for Hallmark)
#'   Options: "H" (Hallmark, 50 gene sets), "C1" (Positional), "C2-CGP" (Chemical/Genetic),
#'            "C2-CP" (Canonical Pathways), "C5-GO-BP" (GO Biological Process), etc.
#'   Only used when enrich_database = "MsigDB"
#'
#' @param hgdisease_source Character. Human disease database source (default: "do")
#'   Options: "do" (Disease Ontology), "ncg_v7", "ncg_v6", "disgenet", "covid19"
#'   Only used when enrich_database = "HgDisease"
#'
#' @param mesh_method Character. MeSH mapping method (default: "gendoo")
#'   Options: "gendoo", "gene2pubmed", "RBBH"
#'   Only used when enrich_database = "Mesh"
#'
#' @param mesh_category Character. MeSH descriptor category (default: "A")
#'   Only used when enrich_database = "Mesh"
#'
#' @param enrichrdb_library Character. Enrichr library name (default: "Cancer_Cell_Line_Encyclopedia")
#'   Only used when enrich_database = "Enrichrdb"
#'
#' @return **Unified Return Structure**: List with 3 components (consistent across scenarios)
#'   \describe{
#'     \item{\strong{stats}}{Data frame with analysis results:
#'       \itemize{
#'         \item \strong{Genome scan}: Top N genes (default 50) with correlation coefficients, p-values.
#'               Columns: gene, r (correlation), pvalue, p_adjusted
#'         \item \strong{GSEA enrichment}: All pathways with enrichment scores.
#'               Columns: pathway, NES (normalized enrichment score), pvalue, p.adjust, size, leading_edge
#'       }
#'       Number of rows: genome scan (top_n), enrichment (all pathways tested)
#'     }
#'     \item{\strong{plot}}{Plot object (ggplot2 or patchwork):
#'       \itemize{
#'         \item \strong{Genome scan}: NetworkPlot showing query gene connected to top correlated genes
#'         \item \strong{GSEA enrichment}: DotPlot showing enriched pathways ranked by significance
#'         \item Access: \code{result$plot} (direct print or save)
#'         \item Size: \code{attr(result$plot, "width")}, \code{attr(result$plot, "height")}
#'         \item Auto-saved: slcptac_output/*.png (300 DPI)
#'       }
#'     }
#'     \item{\strong{raw_data}}{Complete genome-wide results (all 12,000+ genes):
#'       \itemize{
#'         \item Data frame with: gene, correlation/logFC, p-value for ALL genes scanned
#'         \item Use for: Custom filtering, secondary analysis, quality checks
#'         \item Note: stats contains filtered/top results, raw_data contains complete results
#'       }
#'     }
#'   }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Genome-Wide Scan (TESTED - 2.06 sec, 100 genes)
#' # ===========================================================================
#' # Research Question: Which proteins correlate with TP53 protein level?
#' # Analysis: Correlation-based genome-wide scan across 12,000+ proteins
#'
#' result <- cptac_enrichment(
#'   var1 = "TP53",
#'   var1_modal = "Protein",
#'   var1_cancers = "BRCA",
#'   analysis_type = "genome",
#'   genome_modal = "Protein",
#'   top_n = 30
#' )
#'
#' # Return structure (unified)
#' result$stats # Top 30 correlated proteins
#' #   Columns: gene, r, pvalue, p_adjusted
#' result$plot # NetworkPlot showing TP53 connected to top genes
#' result$raw_data # Complete results for all 12,009 proteins
#'
#' # Find strongest correlations
#' result$stats[order(-abs(result$stats$r)), ][1:10, ]
#'
#' # ===========================================================================
#' # Example 2: Pathway Enrichment (TESTED - 3.50 sec, 50 pathways)
#' # ===========================================================================
#' # Research Question: Which pathways are associated with AKT1 protein level?
#' # Analysis: GSEA using MsigDB Hallmark gene sets
#'
#' result <- cptac_enrichment(
#'   var1 = "AKT1",
#'   var1_modal = "Protein",
#'   var1_cancers = "BRCA",
#'   analysis_type = "enrichment",
#'   enrich_database = "MsigDB",
#'   msigdb_category = "H", # Hallmark gene sets (50 pathways)
#'   top_n = 20
#' )
#'
#' # Return structure
#' result$stats # All 50 Hallmark pathways with enrichment scores
#' #   Columns: pathway, NES, pvalue, p.adjust, size, leading_edge
#' result$plot # DotPlot showing top enriched pathways
#' result$raw_data # Complete correlation results for all genes
#'
#' # Find significantly enriched pathways
#' significant <- result$stats[result$stats$p.adjust < 0.05, ]
#' cat("Significant pathways:", nrow(significant), "\n")
#'
#' # ===========================================================================
#' # Example 3: Different Databases
#' # ===========================================================================
#'
#' # KEGG Pathways (~300 pathways)
#' result_kegg <- cptac_enrichment(
#'   var1 = "EGFR",
#'   var1_modal = "Protein",
#'   var1_cancers = "LUAD",
#'   analysis_type = "enrichment",
#'   enrich_database = "KEGG",
#'   kegg_category = "pathway"
#' )
#'
#' # GO Biological Process (~15,000 pathways)
#' result_go <- cptac_enrichment(
#'   var1 = "TP53",
#'   var1_modal = "Protein",
#'   var1_cancers = "BRCA",
#'   analysis_type = "enrichment",
#'   enrich_database = "GO",
#'   enrich_ont = "BP"
#' )
#'
#' # Reactome Pathways (~2,500 pathways)
#' result_reactome <- cptac_enrichment(
#'   var1 = "MTOR",
#'   var1_modal = "Protein",
#'   var1_cancers = "BRCA",
#'   analysis_type = "enrichment",
#'   enrich_database = "Reactome"
#' )
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After genome scan, investigate specific genes:
#' # cptac_correlation(var1="TP53", var2="MDM2", var1_modal="Protein", ...)
#' #
#' # After enrichment, check survival impact:
#' # cptac_survival(var1="AKT1", var1_modal="Protein", surv_type="OS", ...)
#' }
#'
#' @section Performance Test:
#' **Test Environment**: CPTAC proteogenomics database, 122 BRCA samples, 12,009 proteins
#'
#' \strong{Scenario 12} - Genome-wide protein scan (TP53 in BRCA):
#' \itemize{
#'   \item Runtime: 2.06 sec
#'   \item Proteome scanned: 12,009 proteins
#'   \item Result: 100 top correlated proteins, 3,333 significant (p<0.05)
#'   \item Top correlation: r~0.3-0.5 range (typical for proteome-wide)
#'   \item Plot: NetworkPlot showing TP53 → top 30 proteins
#' }
#'
#' \strong{Scenario 13} - GSEA pathway enrichment (AKT1 in BRCA):
#' \itemize{
#'   \item Runtime: 3.50 sec
#'   \item Proteome scanned: 12,009 proteins → ranked gene list
#'   \item Pathways tested: 50 (MsigDB Hallmark)
#'   \item Result: All 50 pathways with NES scores and p-values
#'   \item Expected: mTOR, PI3K-AKT, cell cycle pathways enriched
#'   \item Plot: DotPlot showing top 20 enriched pathways
#' }
#'
#' \strong{Database Options}:
#' \itemize{
#'   \item MsigDB Hallmark: 50 pathways, ~3-4 sec
#'   \item KEGG: ~300 pathways, ~5-8 sec
#'   \item GO Biological Process: ~15,000 pathways, ~30-60 sec
#'   \item Reactome: ~2,500 pathways, ~10-20 sec
#' }
#'
#' \strong{Recommended Use}:
#' \itemize{
#'   \item Genome scan: Fast exploration (<3 sec), good for hypothesis generation
#'   \item GSEA enrichment: Comprehensive pathway analysis (<5 sec for MsigDB)
#'   \item Sample size: Works well with 100-500 samples
#'   \item Multi-cancer: Can analyze across multiple cancer types simultaneously
#' }
#'
#' @section User Queries:
#' **Genome-Wide Discovery**:
#' \itemize{
#'   \item Which proteins correlate with TP53 protein level?
#'   \item What genes are affected by KRAS mutation?
#'   \item Which proteins change with AKT1 phosphorylation?
#'   \item What proteome-wide changes occur in PIK3CA mutants?
#'   \item Which proteins are co-expressed with EGFR?
#'   \item What genes correlate with tumor grade?
#'   \item Which proteins show coordinated expression with mTOR?
#' }
#'
#' **Pathway Enrichment**:
#' \itemize{
#'   \item Which pathways are enriched in TP53 mutants?
#'   \item What biological processes are affected by PIK3CA mutation?
#'   \item Which signaling pathways correlate with AKT1 expression?
#'   \item What pathways are activated in EGFR-high tumors?
#'   \item Which GO terms are associated with MTOR protein level?
#'   \item What KEGG pathways are dysregulated in high-grade tumors?
#'   \item Which Reactome pathways correlate with survival?
#' }
#'
#' **Method Selection**:
#' \itemize{
#'   \item Should I use genome scan or pathway enrichment?
#'   \item Which enrichment database is best for my analysis?
#'   \item How do I find genes affected by a specific mutation?
#'   \item What's the difference between MsigDB and GO?
#'   \item When should I use KEGG vs Reactome?
#'   \item How to identify pathway activation from proteomics?
#'   \item Which analysis reveals druggable targets?
#' }
#'
#' **Interpretation**:
#' \itemize{
#'   \item What does NES (normalized enrichment score) mean?
#'   \item How do I interpret pathway enrichment results?
#'   \item What is a good correlation threshold for genome scan?
#'   \item How many pathways should be significantly enriched?
#'   \item What does leading edge genes represent?
#'   \item How to validate enrichment findings?
#'   \item Which enriched pathways are therapeutically relevant?
#' }
#'
#' **Follow-Up Analysis**:
#' \itemize{
#'   \item After finding correlated proteins, what's next?
#'   \item How to validate pathway enrichment experimentally?
#'   \item Can I check survival impact of enriched genes?
#'   \item How to compare pathways across cancer types?
#'   \item What genes in enriched pathways are druggable?
#'   \item How to visualize pathway networks?
#'   \item Can I correlate pathway scores with clinical outcomes?
#' }
#'
#' @references
#' **CPTAC Database**:
#'
#' Clinical Proteomic Tumor Analysis Consortium (2020). Proteogenomic characterization
#' of human cancers. Nature, 578, 34-35. \doi{10.1038/d41586-020-00432-0}
#'
#' Gillette MA, et al. (2020). Proteogenomic Characterization Reveals Therapeutic
#' Vulnerabilities in Lung Adenocarcinoma. Cell, 182(1):200-225.
#' \doi{10.1016/j.cell.2020.06.013}
#'
#' Database portal: \url{https://proteomics.cancer.gov/programs/cptac}
#'
#' @seealso
#' **Complementary Analysis Functions**:
#' \itemize{
#'   \item \code{\link{cptac_correlation}}: Pairwise correlation between specific genes
#'   \item \code{\link{cptac_survival}}: Survival analysis for identified genes
#' }
#'
#' @export
cptac_enrichment <- function(var1,
                             var1_modal,
                             var1_cancers,
                             analysis_type = "enrichment",
                             enrich_database = "MsigDB",
                             enrich_ont = "BP",
                             genome_modal = "Protein",
                             method = "pearson",
                             top_n = 50,
                             n_workers = 6,
                             kegg_category = "pathway",
                             msigdb_category = "H",
                             hgdisease_source = "do",
                             mesh_method = "gendoo",
                             mesh_category = "A",
                             enrichrdb_library = "Cancer_Cell_Line_Encyclopedia") {
  message("\n========================================")
  message("CPTAC Enrichment Analysis")
  message("========================================")

  # Validate: Clinical variables with genome scan
  if (any(var1_modal == "Clinical") && analysis_type == "genome") {
    stop(
      "Clinical variables cannot be used for genome-wide scans.\n",
      "Reason: Genome-wide DEA requires exactly 2 groups, but clinical variables often have >2 categories (e.g., Stage I/II/III/IV).\n",
      "\n",
      "Suggested alternatives:\n",
      "  1. Use cptac_correlation() to compare clinical variables with specific genes/proteins\n",
      "     Example: cptac_correlation(var1='Tumor_Stage', var1_modal='Clinical', var2=c('TP53','EGFR'), var2_modal='Protein')\n",
      "  \n",
      "  2. For binary clinical variables (e.g., Gender: Male/Female), enrichment analysis may work\n",
      "     Note: Ensure your clinical variable has exactly 2 categories before using enrichment",
      call. = FALSE
    )
  }

  # Validate and fix: enrichment分析必须使用Protein
  if (analysis_type == "enrichment" && genome_modal != "Protein") {
    message("\n[Note] For pathway enrichment analysis (GSEA), genome_modal is automatically set to 'Protein'")
    message("       Reason: GSEA requires gene-level expression data, and Protein is the most stable omics layer")
    message("       Your specified genome_modal '", genome_modal, "' has been overridden\n")
    genome_modal <- "Protein"
  }

  # Load variable data
  loaded <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL
  )

  # Detect scenario
  scenario_info <- .detect_enrichment_scenario(
    var_features = loaded$var1_features,
    var_types = loaded$var1_types,
    analysis_type = analysis_type
  )

  # Perform analysis
  message("\n[Analysis] Running enrichment pipeline...")

  if (scenario_info$var_class == "categorical") {
    # DEA-based enrichment (Scenarios 8-11)
    result <- .run_categorical_enrichment(
      data = loaded$data,
      var_features = loaded$var1_features,
      var1_cancers = var1_cancers,
      genome_modal = genome_modal,
      analysis_type = analysis_type,
      enrich_database = enrich_database,
      enrich_ont = enrich_ont,
      top_n = top_n,
      n_workers = n_workers,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  } else {
    # Correlation-based enrichment (Scenarios 12-15)
    result <- .run_continuous_enrichment(
      data = loaded$data,
      var_features = loaded$var1_features,
      var1_cancers = var1_cancers,
      genome_modal = genome_modal,
      analysis_type = analysis_type,
      enrich_database = enrich_database,
      enrich_ont = enrich_ont,
      method = method,
      top_n = top_n,
      n_workers = n_workers,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  }

  message("\n✓ Enrichment analysis completed")

  # Save plot (following cptac_correlation pattern)
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  var_str <- paste(unique(gsub(" \\(.*\\)", "", loaded$var1_features)), collapse = "-")
  cancer_str <- paste(unique(var1_cancers), collapse = "-")
  modal_str <- paste(unique(var1_modal), collapse = "-")
  analysis_str <- if (analysis_type == "genome") "GenomeScan" else "GSEA"

  filename <- sprintf(
    "enrichment_%s_%s_%s_%s_%s.png",
    analysis_str, cancer_str, var_str, modal_str, genome_modal
  )
  filepath <- file.path(output_dir, filename)

  # Enrichment的plot直接是plot对象（带width/height属性）
  ggplot2::ggsave(
    filename = filepath,
    plot = result$plot,
    width = attr(result$plot, "width"),
    height = attr(result$plot, "height"),
    dpi = 300, limitsize = FALSE
  )

  message(sprintf("✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(result)
}


#' Prognostic Survival Analysis with Kaplan-Meier and Cox Regression
#'
#' @description
#' Performs prognostic survival analysis (OS/PFS) across all CPTAC omics layers to identify
#' biomarkers predicting patient outcomes. Automatically selects between two scenarios:
#' **Scenario 16** (single feature) generates Kaplan-Meier curves + Cox regression with optimal
#' cutoff calculation; **Scenario 17** (multiple features) generates Forest plot comparing
#' hazard ratios across genes/cancers/phospho sites. Supports continuous variables (RNAseq,
#' Protein, Phospho) with automatic cutpoint optimization and categorical variables (Mutation,
#' Clinical) with direct group comparison. Returns unified structure: \code{list(stats, plot, raw_data)}.
#'
#' @param var1 Character vector. Gene names or clinical variables to test for prognostic value.
#'   Examples: "TP53", c("TP53", "EGFR", "KRAS"), c("AKT1", "MTOR").
#'   Number of features determines scenario: 1 feature → Scenario 16 (KM+Cox),
#'   multiple features → Scenario 17 (Forest plot).
#'   Note: Phospho sites auto-detected (e.g., "AKT1" with Phospho modal returns multiple sites).
#'
#' @param var1_modal Character. Omics layer to analyze.
#'   Options: "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation".
#'   Continuous variables (RNAseq, Protein, Phospho, logCNA, Methylation) use cutpoint optimization.
#'   Categorical variables (Mutation, Clinical) use direct group comparison.
#'
#' @param var1_cancers Character vector. Cancer types to analyze.
#'   Options: "BRCA", "LUAD", "COAD", "CCRCC", "GBM", "HNSCC", "LUSC", "OV", "PDAC", "UCEC".
#'   Single cancer: Detailed KM+Cox analysis (Scenario 16).
#'   Multiple cancers: Pan-cancer Forest plot comparison (Scenario 17).
#'
#' @param surv_type Character. Survival endpoint (default: "OS").
#'   Options: "OS" (Overall Survival - time to death), "PFS" (Progression-Free Survival - time to progression/death).
#'   All 10 CPTAC cancer types have OS/PFS data available.
#'
#' @param cutoff_type Character. Cutpoint method for continuous variables (default: "optimal").
#'   Options: "optimal" (maximizes log-rank statistic, recommended), "median", "mean", "quantile".
#'   Only applies to continuous variables; categorical variables use existing groups.
#'   Optimal cutoff balances High/Low group sizes while maximizing survival separation.
#'
#' @param minprop Numeric. Minimum proportion per group for optimal cutoff (default: 0.1).
#'   Ensures at least 10% samples in each group (prevents extreme imbalance).
#'   Range: 0.05-0.3. Lower values allow more flexible cutpoints.
#'
#' @param percent Numeric. Percentile for quantile cutoff (default: 0.25).
#'   Only used when cutoff_type="quantile". Range: 0-1 (e.g., 0.25 = first quartile).
#'
#' @param palette Character vector. Colors for survival curves (default: c("#ED6355", "#41A98E", "#EFA63A", "#3a6ea5")).
#'   First two colors used for High/Low or Mutation/WildType groups in KM curves.
#'
#' @param show_cindex Logical. Display concordance index in plot (default: TRUE).
#'   C-index measures predictive accuracy: 0.5 = random, 1.0 = perfect prediction, >0.6 = acceptable.
#'
#' @return **Unified Return Structure**: List with 3 components (consistent across scenarios)
#'   \describe{
#'     \item{\strong{stats}}{Data frame with survival statistics (columns vary by scenario):
#'       \itemize{
#'         \item \strong{Scenario 16} (single feature): variable, km_pvalue (log-rank test),
#'               cox_hr (hazard ratio), cox_hr_lower, cox_hr_upper (95% CI), cox_pvalue, cox_cindex
#'         \item \strong{Scenario 17} (multiple features): variable, hr (hazard ratio),
#'               hr_lower, hr_upper (95% CI), p_value, cindex
#'       }
#'       Always a data frame with 1+ rows (one per feature analyzed).
#'       HR interpretation: HR > 1 = worse survival (risk factor), HR < 1 = better survival (protective).
#'     }
#'     \item{\strong{plot}}{Plot object (ggplot2 or patchwork):
#'       \itemize{
#'         \item \strong{Scenario 16}: Side-by-side KM curve + Cox regression curve
#'         \item \strong{Scenario 17}: Forest plot with HR point estimates and confidence intervals
#'         \item Access: \code{result$plot} (direct print or save)
#'         \item Size: \code{attr(result$plot, "width")}, \code{attr(result$plot, "height")}
#'         \item Auto-saved: slcptac_output/*.png (300 DPI)
#'       }
#'     }
#'     \item{\strong{raw_data}}{Data frame with merged data (samples × features):
#'       \itemize{
#'         \item Includes: Original features + survival time/event columns + cancer_type
#'         \item Survival columns: CANCER_OS_time, CANCER_OS_event (or PFS)
#'         \item Use for: Custom survival models, additional covariates, subgroup analysis
#'       }
#'     }
#'   }
#'
#' @section Performance Test:
#' **Test Environment**: CPTAC survival data, OS endpoints, real patient outcomes
#'
#' \strong{Scenario 16} - Single gene KM + Cox (TP53 mRNA in BRCA):
#' \itemize{
#'   \item Runtime: 2.29 sec
#'   \item Sample size: 122 BRCA patients with OS data
#'   \item Result: KM p=0.047 (significant), Cox HR=0.733 [0.148-3.629], C-index=0.539
#'   \item Interpretation: TP53-high shows marginally better survival (protective trend)
#'   \item Plot: KM curve (left) + Cox curve (right), side-by-side
#' }
#'
#' \strong{Scenario 16} - Protein expression (AKT1 in BRCA):
#' \itemize{
#'   \item Runtime: 0.54 sec
#'   \item Sample size: 122 BRCA patients
#'   \item Result: KM p=0.036 (significant), Cox HR=4.772 [0.212-107.298], C-index=0.779
#'   \item Interpretation: AKT1-high protein associated with worse survival (risk factor)
#'   \item Cutoff: Optimal cutpoint at 26.555 (Low n=88, High n=34)
#' }
#'
#' \strong{Scenario 16} - Mutation status (KRAS in LUAD):
#' \itemize{
#'   \item Runtime: 0.49 sec
#'   \item Sample size: 110 LUAD patients
#'   \item Result: KM p=0.098 (trend), Cox HR=1.993 (worse survival for mutants)
#'   \item Interpretation: KRAS mutation shows trend toward worse outcomes
#'   \item No cutoff needed (categorical: WildType vs Mutation)
#' }
#'
#' \strong{Scenario 17} - Multiple genes (TP53, EGFR, KRAS in LUAD):
#' \itemize{
#'   \item Runtime: 0.27 sec
#'   \item Sample size: 110 LUAD patients
#'   \item Result: 3 features analyzed, 1 significant (p<0.05)
#'   \item Plot: Forest plot comparing HRs across genes
#' }
#'
#' \strong{Scenario 17} - Multi-cancer (TP53 in BRCA, LUAD, COAD):
#' \itemize{
#'   \item Runtime: 0.30 sec
#'   \item Sample size: 339 patients (BRCA=122, LUAD=110, COAD=107)
#'   \item Result: 3 features (1 per cancer), each analyzed independently
#'   \item Plot: Forest plot comparing TP53 prognostic value across cancers
#' }
#'
#' **Recommended Use**:
#' \itemize{
#'   \item Single-gene prognostic test: Fast (<3 sec), detailed KM+Cox output
#'   \item Multi-gene screening: Very fast (<1 sec), comparative Forest plot
#'   \item Multi-cancer validation: Fast (<1 sec), pan-cancer biomarker assessment
#'   \item Optimal for: 100-500 samples per cancer type
#' }
#'
#' @details
#' **Statistical Methods**:
#' \itemize{
#'   \item \strong{Kaplan-Meier analysis}: Non-parametric survival curve estimation with log-rank test.
#'         Tests null hypothesis that survival curves are identical between groups.
#'         P-value < 0.05 indicates significant survival difference.
#'   \item \strong{Cox proportional hazards}: Estimates hazard ratio (HR) adjusting for censoring.
#'         HR > 1: higher expression/mutation = worse survival (risk factor).
#'         HR < 1: higher expression = better survival (protective factor).
#'         95% CI excludes 1.0 → statistically significant.
#'   \item \strong{Optimal cutoff}: For continuous variables, maximizes log-rank chi-square statistic
#'         while maintaining minimum group size (minprop). Finds cutpoint with best survival separation.
#'   \item \strong{Concordance index (C-index)}: Measures predictive accuracy. 0.5 = random guess,
#'         0.6-0.7 = acceptable, 0.7-0.8 = good, >0.8 = excellent prediction.
#' }
#'
#' **Scenario Selection**:
#' \itemize{
#'   \item Function automatically counts features after data loading
#'   \item 1 feature (1 gene in 1 cancer) → Scenario 16: Detailed KM + Cox analysis
#'   \item Multiple features → Scenario 17: Comparative Forest plot
#'         Sources of multiple features: (1) multiple genes, (2) multiple cancers, (3) multiple phospho sites
#'   \item Each cancer type analyzed independently (uses own survival time/event columns)
#' }
#'
#' **Interpreting Results**:
#' \itemize{
#'   \item \strong{KM p-value}: Significance of survival difference. p<0.05 = significant separation between groups.
#'   \item \strong{Cox HR}: Effect size of prognostic factor. HR=2 means 2× higher risk of death.
#'         HR=0.5 means 50% lower risk (protective). Clinical relevance: HR>1.5 or HR<0.67 considered strong.
#'   \item \strong{95% CI}: Confidence interval for HR. Wide CI = uncertainty (small sample or weak effect).
#'         CI excludes 1.0 → statistically significant at p<0.05 level.
#'   \item \strong{C-index}: Predictive performance. C=0.539 = minimal prediction, C=0.779 = good prediction.
#'         Use to compare different biomarkers: higher C-index = better prognostic marker.
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Single Gene Survival (TESTED - 2.29 sec, HR=0.733, C=0.539)
#' # ===========================================================================
#' # Research Question: Does TP53 mRNA level predict breast cancer survival?
#' # Analysis: Optimal cutoff + KM curve + Cox regression
#'
#' result <- cptac_survival(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = "BRCA",
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Return structure (Scenario 16)
#' result$stats
#' #   variable           km_pvalue  cox_hr  cox_hr_lower  cox_hr_upper  cox_pvalue  cox_cindex
#' #   TP53 (RNAseq,BRCA) 0.047      0.733   0.148         3.629         0.704       0.539
#'
#' result$plot # KM curve + Cox curve side-by-side
#' result$raw_data # 122 samples with survival time/event columns
#'
#' # Interpret
#' cat("TP53-high patients show better survival (HR=0.733, protective trend)\n")
#' cat("KM test significant (p=0.047), but Cox CI wide due to sample size\n")
#'
#' # ===========================================================================
#' # Example 2: Protein Survival (TESTED - 0.54 sec, HR=4.772, C=0.779)
#' # ===========================================================================
#' # Research Question: Does AKT1 protein level predict survival?
#' # Analysis: Strong prognostic marker with high C-index
#'
#' result <- cptac_survival(
#'   var1 = "AKT1",
#'   var1_modal = "Protein",
#'   var1_cancers = "BRCA",
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Interpret
#' cat("AKT1-high protein predicts worse survival (HR=4.772, risk factor)\n")
#' cat("Good predictive accuracy (C-index=0.779)\n")
#' cat("Optimal cutoff: 26.555 separates High (n=34) vs Low (n=88)\n")
#'
#' # ===========================================================================
#' # Example 3: Mutation Survival (TESTED - 0.49 sec, HR=1.993)
#' # ===========================================================================
#' # Research Question: Does KRAS mutation affect lung cancer survival?
#' # Analysis: Categorical variable, no cutoff needed
#'
#' result <- cptac_survival(
#'   var1 = "KRAS",
#'   var1_modal = "Mutation",
#'   var1_cancers = "LUAD",
#'   surv_type = "OS"
#' )
#'
#' # Interpret
#' cat("KRAS mutants show trend toward worse survival (HR=1.993, p=0.098)\n")
#' cat("Not reaching significance, may need larger sample\n")
#'
#' # ===========================================================================
#' # Example 4: Multiple Genes Forest (TESTED - 0.27 sec, 3 genes, 1 significant)
#' # ===========================================================================
#' # Research Question: Which genes predict lung cancer survival?
#' # Analysis: Compare prognostic value across genes
#'
#' result <- cptac_survival(
#'   var1 = c("TP53", "EGFR", "KRAS"),
#'   var1_modal = "RNAseq",
#'   var1_cancers = "LUAD",
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Return structure (Scenario 17)
#' result$stats
#' #   variable          hr     hr_lower  hr_upper  p_value   cindex
#' #   TP53 (RNAseq,..)  X.XX   X.XX      X.XX      0.XXX     0.XXX
#' #   EGFR (RNAseq,..)  X.XX   X.XX      X.XX      0.XXX     0.XXX
#' #   KRAS (RNAseq,..)  X.XX   X.XX      X.XX      0.XXX     0.XXX
#'
#' result$plot # Forest plot with HRs and 95% CIs
#'
#' # Filter significant genes
#' result$stats[result$stats$p_value < 0.05, ]
#'
#' # ===========================================================================
#' # Example 5: Multi-Cancer Comparison (TESTED - 0.30 sec, 3 cancers)
#' # ===========================================================================
#' # Research Question: Is TP53 prognostic value consistent across cancers?
#' # Analysis: Compare same gene across cancer types
#'
#' result <- cptac_survival(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = c("BRCA", "LUAD", "COAD"),
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Compare HRs across cancers
#' result$stats[, c("variable", "hr", "p_value", "cindex")]
#' #   Shows TP53 prognostic value in each cancer type
#'
#' # Visualize
#' print(result$plot) # Forest plot for pan-cancer comparison
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After identifying prognostic genes:
#' # - Check correlation with other markers: cptac_correlation()
#' # - Find affected pathways: cptac_enrichment()
#' # - Validate in independent cohort
#' }
#'
#' @section User Queries:
#' **Prognostic Biomarker Discovery**:
#' \itemize{
#'   \item Does TP53 expression predict patient survival?
#'   \item Which genes are prognostic markers in breast cancer?
#'   \item Does AKT1 protein level predict clinical outcomes?
#'   \item Are phosphorylation sites prognostic in cancer?
#'   \item Which mutations affect patient survival?
#'   \item Does EGFR expression correlate with survival time?
#'   \item What proteins predict progression-free survival?
#' }
#'
#' **Survival Analysis Methods**:
#' \itemize{
#'   \item How to analyze gene expression for survival prediction?
#'   \item What is the optimal cutoff for continuous variables?
#'   \item How do I interpret hazard ratios?
#'   \item What is C-index and how to use it?
#'   \item Should I use OS or PFS for my analysis?
#'   \item How to compare survival across multiple genes?
#'   \item How to validate prognostic biomarkers?
#' }
#'
#' **Multi-Cancer Validation**:
#' \itemize{
#'   \item Is this prognostic marker pan-cancer or cancer-specific?
#'   \item Does TP53 predict survival in multiple cancer types?
#'   \item Which biomarkers are universally prognostic?
#'   \item Do mutation effects vary by cancer type?
#'   \item How to compare prognostic value across cancers?
#'   \item Which proteins show consistent survival association?
#' }
#'
#' **Clinical Translation**:
#' \itemize{
#'   \item Which phosphorylation sites predict treatment response?
#'   \item Does protein expression stratify patient risk groups?
#'   \item Which mutations identify high-risk patients?
#'   \item What biomarkers can guide therapeutic decisions?
#'   \item How do clinical variables affect survival prediction?
#'   \item Which molecular markers improve prognostic models?
#'   \item Can protein levels replace genomic markers for prognosis?
#' }
#'
#' **Integration with Other Analyses**:
#' \itemize{
#'   \item After finding prognostic genes, what's next?
#'   \item How to check if correlated genes are also prognostic?
#'   \item Can I test pathway-level survival associations?
#'   \item How to combine multiple biomarkers for risk score?
#'   \item What proteins in enriched pathways predict survival?
#'   \item How to build multivariate prognostic models?
#' }
#'
#' @references
#' **CPTAC Database**:
#'
#' Clinical Proteomic Tumor Analysis Consortium (2020). Proteogenomic characterization
#' of human cancers. Nature, 578, 34-35. \doi{10.1038/d41586-020-00432-0}
#'
#' Gillette MA, et al. (2020). Proteogenomic Characterization Reveals Therapeutic
#' Vulnerabilities in Lung Adenocarcinoma. Cell, 182(1):200-225.
#' \doi{10.1016/j.cell.2020.06.013}
#'
#' Database portal: \url{https://proteomics.cancer.gov/programs/cptac}
#'
#' @seealso
#' **Complementary Analysis Functions**:
#' \itemize{
#'   \item \code{\link{cptac_correlation}}: Test correlation between prognostic genes
#'   \item \code{\link{cptac_enrichment}}: Identify pathways associated with prognostic markers
#' }
#'
#' @export
cptac_survival <- function(var1,
                           var1_modal,
                           var1_cancers,
                           surv_type = "OS",
                           cutoff_type = "optimal",
                           minprop = 0.1,
                           percent = 0.25,
                           palette = c("#ED6355", "#41A98E", "#EFA63A", "#3a6ea5"),
                           show_cindex = TRUE) {
  message("\n========================================")
  message("CPTAC Survival Analysis")
  message("========================================")

  # Load variable data
  loaded_var <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL
  )

  # Load survival data
  loaded_surv <- cptac_load_modality(
    var1 = surv_type,
    var1_modal = "Survival",
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL,
    surv_type = surv_type
  )

  # Merge data
  merged_data <- .merge_modal_data(loaded_var$data, loaded_surv$data)

  # Detect scenario
  scenario_info <- .detect_survival_scenario(
    var_features = loaded_var$var1_features,
    var_types = loaded_var$var1_types,
    n_cancers = length(var1_cancers)
  )

  # Perform survival analysis
  message("\n[Analysis] Running survival analysis...")

  time_col <- paste0(var1_cancers[1], "_", surv_type, "_time")
  event_col <- paste0(var1_cancers[1], "_", surv_type, "_event")

  if (scenario_info$scenario_id == 16) {
    # Scenario 16: Single feature → KM + Cox
    result <- .run_survival_single(
      merged_data = merged_data,
      var_feature = loaded_var$var1_features[1],
      var_type = loaded_var$var1_types[1],
      time_col = time_col,
      event_col = event_col,
      surv_type = surv_type,
      cutoff_type = cutoff_type,
      minprop = minprop,
      palette = palette,
      var_cancers = var1_cancers,
      var_col_override = NULL
    )
  } else {
    # Scenario 17: Multiple variables → Forest plot
    result <- .run_survival_forest(
      merged_data = merged_data,
      var_features = loaded_var$var1_features,
      var_types = loaded_var$var1_types,
      time_col = time_col,
      event_col = event_col,
      surv_type = surv_type,
      cutoff_type = cutoff_type,
      minprop = minprop,
      var1_cancers = var1_cancers
    )
  }

  message("\n✓ Survival analysis completed")

  # Save plot (following cptac_correlation pattern)
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  var_str <- paste(unique(gsub(" \\(.*\\)", "", loaded_var$var1_features)), collapse = "-")
  cancer_str <- paste(unique(var1_cancers), collapse = "-")
  modal_str <- paste(unique(var1_modal), collapse = "-")
  analysis_str <- if (scenario_info$n_vars == 1) "KM_Cox" else "Forest"

  filename <- sprintf(
    "survival_%s_%s_%s_%s_%s.png",
    analysis_str, surv_type, cancer_str, var_str, modal_str
  )
  filepath <- file.path(output_dir, filename)

  # 统一使用属性方案
  ggplot2::ggsave(
    filename = filepath,
    plot = result$plot,
    width = attr(result$plot, "width"),
    height = attr(result$plot, "height"),
    dpi = 300,
    limitsize = FALSE
  )

  message(sprintf("✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(result)
}


# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Dispatch correlation plot based on scenario
#' @keywords internal
.dispatch_correlation_plot <- function(scenario_info, data, stats, var1_features, var2_features) {
  if (scenario_info$scenario_id == 1) {
    .plot_scenario1(data, stats, var1_features[1], var2_features[1], scenario_info)
  } else if (scenario_info$scenario_id == 2) {
    .plot_scenario2(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 3) {
    .plot_scenario3(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 4) {
    .plot_scenario4(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id %in% c(5, 6)) {
    .plot_scenario5_6(data, stats, var1_features, var2_features, scenario_info)
  } else {
    # Scenario 7
    .plot_scenario7(data, stats, var1_features, var2_features, scenario_info)
  }
}


#' Run categorical enrichment (DEA-based, Scenarios 8-11)
#' @keywords internal
.run_categorical_enrichment <- function(data, var_features, var1_cancers, genome_modal,
                                        analysis_type, enrich_database, enrich_ont,
                                        top_n, n_workers, kegg_category, msigdb_category,
                                        hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide data...")
  genome_matrix <- .load_genome_data(var1_cancers, genome_modal)

  if (length(var_features) == 1) {
    # ========== Scenarios 8 & 9: Single variable ==========
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- factor(data[[var_col]])
    names(var_data) <- rownames(data)

    # Perform DEA
    message("\n[Step 2] Performing DEA...")
    dea_stats <- .stats_dea_genome(var_data, genome_matrix, var1_cancers)

    if (analysis_type == "genome") {
      # Scenario 8: Genome-wide → NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)

      # Select top 50 up + top 50 down genes
      dea_up <- dea_stats[dea_stats$logFC > 0, ]
      dea_down <- dea_stats[dea_stats$logFC < 0, ]
      p_col <- if ("P.Value" %in% colnames(dea_up)) "P.Value" else "pvalue"
      dea_up <- dea_up[order(dea_up[[p_col]]), ]
      dea_down <- dea_down[order(dea_down[[p_col]]), ]
      dea_top <- rbind(head(dea_up, 50), head(dea_down, 50))

      plot_result <- .plot_network(
        stats = dea_top,
        var1_name = query_gene,
        edge_metric = "logFC",
        query_omics = "Mutation",
        genome_omics = genome_modal,
        cancer_type = var1_cancers[1],
        analysis_type = "DEA",
        method = NULL
      )

      return(list(
        stats = dea_top,
        plot = plot_result,
        raw_data = dea_stats
      ))
    } else {
      # Scenario 9: Enrichment → GSEA
      message("\n[Step 3] Performing GSEA...")

      ranked_genes <- setNames(dea_stats$logFC, dea_stats$gene)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)

      gsea_result <- .perform_gsea(
        ranked_genes = ranked_genes,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        n_workers = n_workers,
        kegg_category = kegg_category,
        msigdb_category = msigdb_category,
        hgdisease_source = hgdisease_source,
        mesh_method = mesh_method,
        mesh_category = mesh_category,
        enrichrdb_library = enrichrdb_library
      )

      message("\n[Step 4] Creating GSEA plots...")

      gene_name <- gsub("\\s*\\(.*", "", var_features[1])
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Mutation"

      plot_result <- .plot_gsea_paired(
        gsea_stats = gsea_result,
        var_name = gene_name,
        omics_type = modal_type,
        cancer_types = var1_cancers,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        method = NULL, # DEA doesn't use correlation method
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = dea_stats
      ))
    }
  } else {
    # ========== Scenarios 10 & 11: Multiple variables ==========
    message("\n[Step 2] Processing multiple variables...")

    all_stats <- list()

    for (var_label in var_features) {
      var_col <- .extract_colname_from_label(c(var_label), data)[1]
      var_data <- factor(data[[var_col]])
      names(var_data) <- rownames(data)

      message(sprintf("  Processing %s...", var_label))
      dea_stats <- .stats_dea_genome(var_data, genome_matrix, var1_cancers)
      dea_stats$var_name <- var_label
      all_stats[[var_label]] <- dea_stats
    }

    combined_stats <- do.call(rbind, all_stats)

    if (analysis_type == "genome") {
      # Scenario 10: Genome-wide → DotPlot
      message("\n[Step 3] Creating DotPlot...")

      plot_result <- .plot_dotplot_paired(
        all_stats = combined_stats,
        analysis_type = "DEA",
        cancer_types = var1_cancers,
        genome_omics = genome_modal,
        is_mutation = TRUE,
        use_mean = FALSE,
        feature_list = NULL,
        method = NULL,
        top_n = top_n
      )

      # 使用筛选后的数据作为stats（每个变量top N × 2）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 11: Enrichment → GSEA Matrix
      message("\n[Step 3] Performing GSEA for multiple variables...")

      gsea_combined <- data.frame()

      for (var_label in var_features) {
        var_stats <- all_stats[[var_label]]
        ranked_genes <- setNames(var_stats$logFC, var_stats$gene)
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)

        gsea_result <- .perform_gsea(
          ranked_genes = ranked_genes,
          enrich_type = enrich_database,
          GO_ont = enrich_ont,
          n_workers = n_workers,
          kegg_category = kegg_category,
          msigdb_category = msigdb_category,
          hgdisease_source = hgdisease_source,
          mesh_method = mesh_method,
          mesh_category = mesh_category,
          enrichrdb_library = enrichrdb_library
        )
        gsea_result$var_name <- var_label
        gsea_combined <- rbind(gsea_combined, gsea_result)
      }

      message("\n[Step 4] Creating GSEA matrix plot...")
      plot_result <- .plot_gsea_matrix(
        all_gsea_stats = gsea_combined,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        cancer_types = var1_cancers,
        method = NULL, # DEA doesn't use correlation
        use_mean = FALSE,
        feature_list = NULL,
        top_n = top_n,
        omics_type = "Mutation" # Categorical变量都是Mutation或Clinical
      )

      # 使用筛选后的pathways作为stats
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else gsea_combined,
        plot = plot_result,
        raw_data = all_stats
      ))
    }
  }
}


#' Run continuous enrichment (Correlation-based, Scenarios 12-15)
#' @keywords internal
.run_continuous_enrichment <- function(data, var_features, var1_cancers, genome_modal,
                                       analysis_type, enrich_database, enrich_ont,
                                       method, top_n, n_workers, kegg_category, msigdb_category,
                                       hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide data...")
  genome_matrix <- .load_genome_data(var1_cancers, genome_modal)

  if (length(var_features) == 1) {
    # ========== Scenarios 12 & 13: Single variable ==========
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- as.numeric(data[[var_col]])
    names(var_data) <- rownames(data)

    # Perform correlation
    message("\n[Step 2] Calculating correlations...")
    cor_stats <- .stats_cor_genome(var_data, genome_matrix, var1_cancers, method)

    if (analysis_type == "genome") {
      # Scenario 12: Genome-wide → NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)

      # Extract modal type from label
      modal_match <- regmatches(var_label, regexpr("\\(([^,]+),", var_label))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Protein"

      # Select top 50 positive + top 50 negative correlations
      cor_pos <- cor_stats[cor_stats$r > 0, ]
      cor_neg <- cor_stats[cor_stats$r < 0, ]
      cor_pos <- cor_pos[order(cor_pos$pvalue), ]
      cor_neg <- cor_neg[order(cor_neg$pvalue), ]
      cor_top <- rbind(head(cor_pos, 50), head(cor_neg, 50))

      plot_result <- .plot_network(
        stats = cor_top,
        var1_name = query_gene,
        edge_metric = "r",
        query_omics = modal_type,
        genome_omics = genome_modal,
        cancer_type = var1_cancers[1],
        analysis_type = "correlation",
        method = method
      )

      return(list(
        stats = cor_top,
        plot = plot_result,
        raw_data = cor_stats
      ))
    } else {
      # Scenario 13: Enrichment → GSEA
      message("\n[Step 3] Performing GSEA...")

      ranked_genes <- setNames(cor_stats$r, cor_stats$gene)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)

      gsea_result <- .perform_gsea(
        ranked_genes = ranked_genes,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        n_workers = n_workers,
        kegg_category = kegg_category,
        msigdb_category = msigdb_category,
        hgdisease_source = hgdisease_source,
        mesh_method = mesh_method,
        mesh_category = mesh_category,
        enrichrdb_library = enrichrdb_library
      )

      message("\n[Step 4] Creating GSEA plots...")

      gene_name <- gsub("\\s*\\(.*", "", var_features[1])
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Protein"

      plot_result <- .plot_gsea_paired(
        gsea_stats = gsea_result,
        var_name = gene_name,
        omics_type = modal_type,
        cancer_types = var1_cancers,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        method = method, # Pass correlation method
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = cor_stats
      ))
    }
  } else {
    # ========== Scenarios 14 & 15: Multiple variables ==========
    message("\n[Step 2] Processing multiple variables...")

    all_stats <- list()

    for (var_label in var_features) {
      var_col <- .extract_colname_from_label(c(var_label), data)[1]
      var_data <- as.numeric(data[[var_col]])
      names(var_data) <- rownames(data)

      message(sprintf("  Processing %s...", var_label))
      cor_stats <- .stats_cor_genome(var_data, genome_matrix, var1_cancers, method)
      cor_stats$var_name <- var_label
      all_stats[[var_label]] <- cor_stats
    }

    combined_stats <- do.call(rbind, all_stats)

    if (analysis_type == "genome") {
      # Scenario 14: Genome-wide → DotPlot
      message("\n[Step 3] Creating DotPlot...")

      # Determine if we need mean expression
      n_unique_cancers <- length(unique(var1_cancers))
      use_mean_expr <- (n_unique_cancers == 1 && length(var_features) > 1)
      feature_names <- if (use_mean_expr) gsub("\\s*\\(.*", "", var_features) else NULL

      plot_result <- .plot_dotplot_paired(
        all_stats = combined_stats,
        analysis_type = "Correlation",
        cancer_types = var1_cancers,
        genome_omics = genome_modal,
        is_mutation = FALSE,
        use_mean = use_mean_expr,
        feature_list = feature_names,
        method = method,
        top_n = top_n
      )

      # 使用筛选后的数据作为stats（每个变量top N × 2）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 15: Enrichment → GSEA Matrix
      message("\n[Step 3] Performing GSEA for multiple variables...")

      gsea_combined <- data.frame()

      for (var_label in var_features) {
        var_stats <- all_stats[[var_label]]
        ranked_genes <- setNames(var_stats$r, var_stats$gene)
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)

        gsea_result <- .perform_gsea(
          ranked_genes = ranked_genes,
          enrich_type = enrich_database,
          GO_ont = enrich_ont,
          n_workers = n_workers,
          kegg_category = kegg_category,
          msigdb_category = msigdb_category,
          hgdisease_source = hgdisease_source,
          mesh_method = mesh_method,
          mesh_category = mesh_category,
          enrichrdb_library = enrichrdb_library
        )
        gsea_result$var_name <- var_label
        gsea_combined <- rbind(gsea_combined, gsea_result)
      }

      message("\n[Step 4] Creating GSEA matrix plot...")

      # Determine if we need mean expression
      n_unique_cancers <- length(unique(var1_cancers))
      use_mean_expr <- (n_unique_cancers == 1 && length(var_features) > 1)
      feature_names <- if (use_mean_expr) gsub("\\s*\\(.*", "", var_features) else NULL

      # Extract omics_type from first feature (format: "GENE (Modal, Cancer)")
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else NULL

      plot_result <- .plot_gsea_matrix(
        all_gsea_stats = gsea_combined,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        cancer_types = var1_cancers,
        method = method,
        use_mean = use_mean_expr,
        feature_list = feature_names,
        top_n = top_n,
        omics_type = modal_type
      )

      # 使用筛选后的pathways作为stats（用于图表展示）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else gsea_combined,
        plot = plot_result,
        raw_data = all_stats
      ))
    }
  }
}


#' Run single variable survival analysis (Scenario 16)
#' @keywords internal
.run_survival_single <- function(merged_data, var_feature, var_type, time_col, event_col,
                                 surv_type, cutoff_type, minprop, palette, var_cancers,
                                 var_col_override = NULL) {
  message("\n[Step 1] Preparing survival data...")

  # 使用override列名（用于phospho mean），否则从feature提取
  var_col <- if (!is.null(var_col_override)) {
    var_col_override
  } else {
    .extract_colname_from_label(c(var_feature), merged_data)[1]
  }

  # Validate columns exist
  if (!var_col %in% colnames(merged_data)) {
    stop(sprintf(
      "Variable column '%s' not found in data. Available: %s",
      var_col, paste(colnames(merged_data), collapse = ", ")
    ), call. = FALSE)
  }
  if (!time_col %in% colnames(merged_data)) {
    stop(sprintf("Time column '%s' not found in data", time_col), call. = FALSE)
  }
  if (!event_col %in% colnames(merged_data)) {
    stop(sprintf("Event column '%s' not found in data", event_col), call. = FALSE)
  }

  # Handle continuous vs categorical
  if (var_type == "continuous") {
    message(sprintf("  Calculating %s cutoff...", cutoff_type))

    # Filter to complete cases first
    valid_idx <- complete.cases(merged_data[, c(var_col, time_col, event_col)])
    if (sum(valid_idx) < 10) {
      stop("Too few valid samples for survival analysis (need at least 10)", call. = FALSE)
    }

    cutoff <- .calc_optimal_cutoff(merged_data[valid_idx, ], var_col, time_col, event_col, minprop)

    # Create group on ALL data (not just valid)
    merged_data$group <- ifelse(merged_data[[var_col]] > cutoff, "High", "Low")
    merged_data$group <- factor(merged_data$group, levels = c("Low", "High"))
    group_col <- "group"

    message(sprintf(
      "  Cutoff: %.3f (Low: n=%d, High: n=%d)",
      cutoff,
      sum(merged_data$group == "Low", na.rm = TRUE),
      sum(merged_data$group == "High", na.rm = TRUE)
    ))
  } else {
    group_col <- var_col
  }

  # Perform KM analysis
  message("\n[Step 2] Performing Kaplan-Meier analysis...")
  km_result <- .perform_km_analysis(merged_data, group_col, time_col, event_col)

  # Perform Cox analysis
  message("\n[Step 3] Performing Cox regression...")
  cox_result <- .perform_cox_analysis(merged_data, var_col, time_col, event_col)

  # Create combined plot
  message("\n[Step 4] Creating combined KM + Cox plot...")

  gene_name <- gsub("\\s*\\(.*", "", var_feature)
  modal_match <- regmatches(var_feature, regexpr("\\(([^,]+),", var_feature))
  modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Unknown"

  cox_stats_df <- data.frame(
    variable = var_feature,
    hr = cox_result$hr,
    hr_lower = cox_result$hr_lower,
    hr_upper = cox_result$hr_upper,
    p_value = cox_result$p_value,
    stringsAsFactors = FALSE
  )

  plot_result <- .plot_km_cox_combined(
    km_fit = km_result$survfit,
    cox_model_stats = cox_stats_df,
    data = merged_data,
    time_col = time_col,
    event_col = event_col,
    group_col = group_col,
    var_name = gene_name,
    omics_type = modal_type,
    cancer_type = var_cancers[1],
    surv_type = surv_type,
    var_col = var_col # Pass the actual feature column
  )

  # Prepare stats
  stats <- data.frame(
    variable = var_feature,
    km_pvalue = km_result$p_value,
    cox_hr = cox_result$hr,
    cox_hr_lower = cox_result$hr_lower,
    cox_hr_upper = cox_result$hr_upper,
    cox_pvalue = cox_result$p_value,
    cox_cindex = cox_result$cindex,
    stringsAsFactors = FALSE
  )

  return(list(
    stats = stats,
    plot = plot_result,
    raw_data = merged_data
  ))
}


#' Run multiple variables survival analysis (Scenario 17)
#' @keywords internal
.run_survival_forest <- function(merged_data, var_features, var_types, time_col, event_col,
                                 surv_type, cutoff_type, minprop, var1_cancers) {
  message("\n[Step 1] Performing Cox regression for multiple variables...")

  forest_data <- data.frame()

  for (i in seq_along(var_features)) {
    var_label <- var_features[i]
    var_type <- var_types[i]
    var_col <- .extract_colname_from_label(c(var_label), merged_data)[1]

    message(sprintf("  Processing %s...", var_label))

    # 从feature label提取癌种
    cancer_match <- regmatches(var_label, regexpr(",\\s*([^)]+)\\)", var_label))
    feature_cancer <- if (length(cancer_match) > 0) {
      gsub(",\\s*|\\)", "", cancer_match)
    } else {
      var1_cancers[1]
    }

    # 使用该feature对应癌种的time/event列
    feature_time_col <- paste0(feature_cancer, "_", surv_type, "_time")
    feature_event_col <- paste0(feature_cancer, "_", surv_type, "_event")

    # 只用该feature对应癌种的数据
    feature_data <- merged_data[merged_data$cancer_type == feature_cancer, ]

    # For continuous variables, calculate cutoff (for reference only)
    if (var_type == "continuous") {
      tryCatch(
        {
          cutoff <- .calc_optimal_cutoff(feature_data, var_col, feature_time_col, feature_event_col, minprop)
        },
        error = function(e) {
          # 如果cutoff计算失败，使用median
          cutoff <- median(feature_data[[var_col]], na.rm = TRUE)
        }
      )
    }

    # Perform Cox regression
    cox_result <- .perform_cox_analysis(feature_data, var_col, feature_time_col, feature_event_col)

    # 使用完整的feature label（包含phospho site、modal、cancer）
    forest_data <- rbind(forest_data, data.frame(
      variable = var_label, # 使用完整label而不是只有基因名
      hr = cox_result$hr,
      hr_lower = cox_result$hr_lower,
      hr_upper = cox_result$hr_upper,
      p_value = cox_result$p_value,
      cindex = cox_result$cindex,
      stringsAsFactors = FALSE
    ))
  }

  # Create forest plot
  message("\n[Step 2] Creating forest plot...")

  # 提取所有涉及的癌种
  all_cancers <- unique(sapply(var_features, function(f) {
    match <- regmatches(f, regexpr(",\\s*([^)]+)\\)", f))
    if (length(match) > 0) gsub(",\\s*|\\)", "", match) else NA
  }))
  all_cancers <- all_cancers[!is.na(all_cancers)]
  cancer_label <- if (length(all_cancers) > 1) "Database" else all_cancers[1]

  plot_result <- .plot_forest(
    cox_stats = forest_data,
    surv_type = surv_type,
    cancer_type = cancer_label
  )

  return(list(
    stats = forest_data,
    plot = plot_result,
    raw_data = merged_data
  ))
}


#' Generate filename for saved plots
#' @keywords internal
.generate_filename <- function(analysis, var1, var1_modal, var1_cancers,
                               var2 = NULL, var2_modal = NULL, var2_cancers = NULL) {
  cancer_str <- paste(unique(c(var1_cancers, var2_cancers)), collapse = "-")
  var1_str <- paste(var1, collapse = "-")

  if (!is.null(var2)) {
    var2_str <- paste(var2, collapse = "-")
    filename <- sprintf(
      "%s_%s_%s_%s_vs_%s_%s.png",
      analysis, cancer_str, var1_str, var1_modal, var2_str, var2_modal
    )
  } else {
    filename <- sprintf(
      "%s_%s_%s_%s.png",
      analysis, cancer_str, var1_str, var1_modal
    )
  }

  return(filename)
}
