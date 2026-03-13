# ==============================================================================
# Package Initialization
# ==============================================================================
# This file contains functions that run when the package is loaded
# ==============================================================================

#' Package Load Hook
#'
#' @description
#' Called when the package is loaded. Sets up default options.
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  .ggforge_init()
  # Set ragg max dimension to 100000px (default is 50000px)
  # This allows saving large plots with many variables/genes
  # Example: 36 phospho sites × 30 genes each = very large plot
  options(ragg.max_dim = 100000)
}
