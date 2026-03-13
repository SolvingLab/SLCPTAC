# =============================================================================
# astat internal functions (inlined to remove external dependency)
# Original package: astat by Zaoqu Liu
# License: GPL (>= 3)
# =============================================================================

#' @title Smart Parameter Matching
#' @description Intelligently matches user input to valid options with case-insensitive and partial matching support.
#' @param input User input string
#' @param valid_options Character vector of valid options
#' @param param_name Parameter name for error messages
#' @return Matched option from valid_options
#' @noRd
smart_match <- function(input, valid_options, param_name) {
  if (is.null(input) || is.na(input)) {
    return(NULL)
  }

  input_lower <- tolower(as.character(input))
  valid_lower <- tolower(valid_options)

  # Exact match (case-insensitive)
  exact_match <- valid_options[valid_lower == input_lower]
  if (length(exact_match) == 1) {
    return(exact_match)
  }
  if (length(exact_match) > 1) {
    return(exact_match[1])
  }

  # Partial match (prefix matching)
  partial_matches <- valid_options[startsWith(valid_lower, input_lower)]

  if (length(partial_matches) == 0) {
    stop(sprintf(
      "Invalid %s '%s'. Valid options: %s",
      param_name, input, paste(valid_options, collapse = ", ")
    ))
  } else if (length(partial_matches) == 1) {
    return(partial_matches)
  } else {
    # Multiple matches - apply priority rules

    # Rule 1: Prefer longer matches
    match_lengths <- nchar(valid_lower[valid_options %in% partial_matches])
    if (length(unique(match_lengths)) > 1) {
      longest <- max(match_lengths)
      longest_matches <- partial_matches[match_lengths == longest]
      if (length(longest_matches) == 1) {
        return(longest_matches)
      }
      partial_matches <- longest_matches
    }

    # Rule 2: Prefer matches without special characters
    no_special <- partial_matches[!grepl("[._]", partial_matches)]
    if (length(no_special) == 1) {
      return(no_special)
    } else if (length(no_special) > 1) {
      partial_matches <- no_special
    }

    # Rule 3: Prefer all-lowercase matches
    all_lower <- partial_matches[partial_matches == tolower(partial_matches)]
    if (length(all_lower) == 1) {
      return(all_lower)
    } else if (length(all_lower) > 1) {
      partial_matches <- all_lower
    }

    # Special case for p.adjust.method
    if (param_name == "p.adjust.method" && "bonferroni" %in% partial_matches) {
      return("bonferroni")
    }

    stop(sprintf(
      "Ambiguous %s '%s'. Could match: %s. Please be more specific.",
      param_name, input, paste(partial_matches, collapse = ", ")
    ))
  }
}


#' @title Correlation and P-value Calculation
#' @description Efficiently calculates correlation matrix and corresponding p-values for large matrices with varying numbers of missing data.
#' @param x Numeric vector, matrix or data frame
#' @param y Numeric vector, matrix or data frame, or NULL
#' @param use Method for handling missing values
#' @param alternative Alternative hypothesis: "two.sided", "less", or "greater"
#' @param ... Additional arguments passed to cor()
#' @return List with cor (correlation matrix), p (p-value matrix), and nObs (number of observations)
#' @noRd
CorPvalue <- function(x, y = NULL,
                      use = "pairwise.complete.obs",
                      alternative = "two.sided",
                      ...) {
  # Handle vector inputs - convert to single-column matrix with names
  if (is.null(dim(x))) {
    # x is a vector
    x <- matrix(x, ncol = 1)
    colnames(x) <- "x"
  } else {
    x <- as.matrix(x)
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("Var", seq_len(ncol(x)))
    }
  }

  if (!is.null(y)) {
    if (is.null(dim(y))) {
      # y is a vector
      y <- matrix(y, ncol = 1)
      colnames(y) <- "y"
    } else {
      y <- as.matrix(y)
      if (is.null(colnames(y))) {
        colnames(y) <- paste0("Var", seq_len(ncol(y)))
      }
    }
  }

  # Calculate correlation matrix
  cor_result <- cor(x, y, use = use, ...)

  # Ensure cor_result is a matrix with proper dimnames
  if (!is.matrix(cor_result)) {
    cor_result <- as.matrix(cor_result)
  }
  if (is.null(rownames(cor_result))) {
    rownames(cor_result) <- colnames(x)
  }
  if (is.null(colnames(cor_result))) {
    if (is.null(y)) {
      colnames(cor_result) <- colnames(x)
    } else {
      colnames(cor_result) <- colnames(y)
    }
  }

  # Calculate finite value matrix
  finMat <- !is.na(x)

  # Calculate number of observations for each pair
  if (is.null(y)) {
    np <- t(finMat) %*% finMat
  } else {
    np <- t(finMat) %*% (!is.na(y))
  }

  # Calculate t-statistic
  T <- sqrt(np - 2) * cor_result / sqrt(1 - cor_result^2)

  # Calculate p-values based on alternative hypothesis
  if (alternative == "two.sided") {
    p <- 2 * pt(abs(T), np - 2, lower.tail = FALSE)
  } else if (alternative == "less") {
    p <- pt(T, np - 2, lower.tail = TRUE)
  } else if (alternative == "greater") {
    p <- pt(T, np - 2, lower.tail = FALSE)
  }

  # Ensure p has the same dimnames as cor_result
  dimnames(p) <- dimnames(cor_result)

  list(cor = cor_result, p = p, nObs = np)
}


#' @title Matrix of Correlations and P-values
#' @description Calculates the correlation matrix along with p-values. Supports multiple correlation methods and handles missing values. Can adjust p-values for multiple comparisons.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param x A numeric matrix or data frame.
#' @param y A numeric matrix or data frame, or NULL. If NULL, correlations are calculated within x.
#' @param cor.method Correlation method: "pearson" (or "p"), "kendall" (or "k"), or "spearman" (or "s"). Case-insensitive. Default: "pearson"
#' @param p.adjust Logical. If TRUE, adjust p-values for multiple comparisons. Default: FALSE
#' @param p.adjust.method Method for p-value adjustment: "holm", "hochberg", "hommel", "bonferroni" (or "b"), "BH" (or "bh"), "BY" (or "by"), "fdr", "none". Case-insensitive. Default: "holm"
#' @param use Method for handling missing values: "everything", "all.obs", "complete.obs" (or "c"), "na.or.complete", "pairwise.complete.obs" (or "p"). Default: "everything"
#' @param alternative Alternative hypothesis for correlation test: "two.sided", "less", or "greater". Default: "two.sided"
#' @param ... Additional arguments passed to cor().
#' @return A data frame with columns: x (first variable), y (second variable), r (correlation), p (p-value).
# @keywords internal
#' @examples
#' # Basic usage
#' result <- stat_cor(mtcars[, 1:5])
#'
#' # With abbreviations
#' result <- stat_cor(mtcars[, 1:5], cor.method = "s", p.adjust = TRUE, p.adjust.method = "b")
#'
#' # Between two matrices
#' result <- stat_cor(mtcars[, 1:3], mtcars[, 4:6])
stat_cor <- function(x,
                     y = NULL,
                     cor.method = "pearson",
                     p.adjust = FALSE,
                     p.adjust.method = "holm",
                     use = "everything",
                     alternative = "two.sided",
                     ...) {
  # Check and install required packages
  required_packages <- c("tidyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing package '%s'...", pkg))
      utils::install.packages(pkg, quiet = TRUE)
    }
  }

  # Match parameters using smart matching
  cor.method <- smart_match(cor.method, c("pearson", "kendall", "spearman"), "cor.method")
  use <- smart_match(use, c(
    "everything", "all.obs", "complete.obs",
    "na.or.complete", "pairwise.complete.obs"
  ), "use")
  alternative <- smart_match(alternative, c("two.sided", "less", "greater"), "alternative")

  if (p.adjust) {
    p.adjust.method <- smart_match(
      p.adjust.method,
      c(
        "holm", "hochberg", "hommel", "bonferroni",
        "BH", "BY", "fdr", "none"
      ),
      "p.adjust.method"
    )
  }

  # Calculate correlations and p-values
  corr <- CorPvalue(x, y, use = use, alternative = alternative, method = cor.method)

  # Adjust p-values if requested
  if (p.adjust) {
    corr$p <- matrix(
      stats::p.adjust(as.vector(corr$p), method = p.adjust.method),
      nrow = nrow(corr$p),
      ncol = ncol(corr$p)
    )
    dimnames(corr$p) <- dimnames(corr$cor)
  }

  # Convert to tidy format
  d1 <- corr$cor %>%
    as.data.frame() %>%
    { d <- .; data.frame(x = rownames(d), d, check.names = FALSE, stringsAsFactors = FALSE) } %>%
    tidyr::pivot_longer(2:ncol(.), names_to = "y", values_to = "r") %>%
    as.data.frame()

  if (p.adjust) {
    d1$p <- as.vector(corr$p)
    d3 <- d1
  } else {
    d2 <- corr$p %>%
      as.data.frame() %>%
      { d <- .; data.frame(x = rownames(d), d, check.names = FALSE, stringsAsFactors = FALSE) } %>%
      tidyr::pivot_longer(2:ncol(.), names_to = "y", values_to = "p") %>%
      as.data.frame()
    d3 <- merge(d1, d2, by = c("x", "y"))
  }

  return(d3)
}
