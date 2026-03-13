# =============================================================================
# ggforge internal functions (inlined to remove external dependency)
# Original package: ggforge by Zaoqu Liu
# License: GPL (>= 3)
# =============================================================================

# Palette list cache (loaded from ggforge data on first use)
.ggforge_palette_cache <- new.env(parent = emptyenv())

.ggforge_palette_list_fn <- function() {
  if (is.null(.ggforge_palette_cache$data)) {
    tryCatch({
      e <- new.env(parent = emptyenv())
      utils::data("palette_list", package = "SLCPTAC", envir = e)
      .ggforge_palette_cache$data <- e$palette_list
    }, error = function(err) {
      .ggforge_palette_cache$data <- list(
        Paired = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
        Spectral = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
        RdBu = c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061"),
        Set1 = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999"),
        Dark2 = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
        stripe = c("#FFFFFF","#F0F0F0")
      )
    })
  }
  .ggforge_palette_cache$data
}



# --- Source: ggforge/R/00-config.R ---

#' ggforge Global Configuration
#'
#' @description
#' Configuration system for ggforge package. This provides a centralized
#' way to manage default settings and options.
#'
#' @name config
#' @keywords internal
NULL

#' Get or set ggforge options
#'
#' @param ... Named arguments to set options, or character strings to get options
#' @return If setting options, returns invisible NULL. If getting options, returns the option value.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Set options
#' ggforge_option(theme.base_size = 14, theme.font_family = "Arial")
#'
#' # Get options
#' ggforge_option("theme.base_size")
#' }
ggforge_option <- function(...) {
  args <- list(...)

  if (length(args) == 0) {
    # Return all ggforge options
    opts <- options()
    opts[grepl("^ggforge\\.", names(opts))]
  } else if (length(args) == 1 && is.null(names(args))) {
    # Get single option
    opt_name <- paste0("ggforge.", args[[1]])
    getOption(opt_name, default = .ggforge_defaults[[args[[1]]]])
  } else {
    # Set options
    opt_names <- paste0("ggforge.", names(args))
    names(args) <- opt_names
    do.call(options, args)
    invisible(NULL)
  }
}

#' Default configuration values
#' @keywords internal
.ggforge_defaults <- list(
  theme.base_size = 12,
  theme.font_family = NULL,
  theme.default = "theme_ggforge",
  palette.default = "Paired",
  seed = 8525,
  gglogger.enabled = FALSE
)

#' Initialize ggforge options on package load
#' @keywords internal
.ggforge_init <- function() {
  # Set default options if not already set
  for (opt_name in names(.ggforge_defaults)) {
    full_opt_name <- paste0("ggforge.", opt_name)
    if (is.null(getOption(full_opt_name))) {
      opt_list <- list(.ggforge_defaults[[opt_name]])
      names(opt_list) <- full_opt_name
      do.call(options, opt_list)
    }
  }
}

# --- Source: ggforge/R/00-specs.R ---

#' ggforge Style Specifications System
#'
#' @description
#' Centralized style specifications for all visual elements.
#' This is the authoritative source for all styling rules in ggforge.
#'
#' @name style-specs
#' @keywords internal
NULL

# =============================================================================
# STYLE SPECIFICATIONS - Single Source of Truth
# =============================================================================

#' Core style specifications
#' @keywords internal
.ggforge_style_specs <- list(
  # Font specifications
  font = list(
    # Plot title
    title = list(
      size = 13,
      face = "bold",
      hjust = 0.5,
      colour = "black"
    ),

    # Axis titles
    axis_title = list(
      size = 12,
      face = "bold",
      colour = "black"
    ),

    # Axis text - varies by variable type
    axis_text = list(
      continuous = list(
        size = 10,
        colour = "black"
      ),
      discrete = list(
        size = 12,
        colour = "black"
      )
    ),

    # Facet strip text
    facet_text = list(
      size = 12,
      face = "bold",
      colour = "black"
    ),

    # Legend title
    legend_title = list(
      size = 12,
      face = "bold",
      colour = "black"
    ),

    # Legend text
    legend_text = list(
      size = 10,
      colour = "black"
    )
  ),

  # Line specifications
  lines = list(
    base_size = 0.6,
    grid_major = list(
      colour = "grey80",
      linetype = 2,
      linewidth = 0.5
    ),
    grid_minor = list(
      colour = "grey90",
      linetype = 2,
      linewidth = 0.3
    )
  ),

  # Spacing specifications
  spacing = list(
    plot_margin = c(10, 10, 10, 10),
    panel_spacing = 5
  )
)

#' Get style specification by path
#'
#' @description
#' Retrieves a style specification using dot-notation path.
#' Automatically scales font sizes based on base_size.
#'
#' @param path Character. Dot-separated path (e.g., "font.axis_text.continuous")
#' @param base_size Numeric. Base font size for scaling (default: 12)
#' @return The specification value (list or atomic)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' get_style_spec("font.title")
#' get_style_spec("font.axis_text.continuous", base_size = 14)
#' }
get_style_spec <- function(path, base_size = 12) {
  # Split path by dots (fixed = TRUE means literal dot, no need for \\)
  parts <- strsplit(path, ".", fixed = TRUE)[[1]]

  # Get specs - direct access works in both installed and dev mode
  spec <- .ggforge_style_specs

  # Navigate to the specification
  for (part in parts) {
    if (is.null(spec[[part]])) {
      stop("Style specification not found: ", path,
        "\nAvailable at this level: ", paste(names(spec), collapse = ", "),
        call. = FALSE
      )
    }
    spec <- spec[[part]]
  }

  # Scale font sizes if present
  if (is.list(spec) && "size" %in% names(spec)) {
    spec$size <- spec$size * base_size / 12
  }

  return(spec)
}

#' Build element_text from style spec
#'
#' @description
#' Converts a style specification to ggplot2::element_text()
#'
#' @param spec_path Character. Path to style specification
#' @param base_size Numeric. Base font size
#' @param ... Additional arguments to override spec
#' @return A ggplot2::element_text object
#' @keywords internal
build_element_text <- function(spec_path, base_size = 12, ...) {
  spec <- get_style_spec(spec_path, base_size)

  # Merge with user overrides
  overrides <- list(...)
  if (length(overrides) > 0) {
    spec <- utils::modifyList(spec, overrides)
  }

  do.call(ggplot2::element_text, spec)
}

#' Build element_line from style spec
#'
#' @description
#' Converts a style specification to ggplot2::element_line()
#'
#' @param spec_path Character. Path to style specification
#' @param ... Additional arguments to override spec
#' @return A ggplot2::element_line object
#' @keywords internal
build_element_line <- function(spec_path, ...) {
  spec <- get_style_spec(spec_path)

  # Merge with user overrides
  overrides <- list(...)
  if (length(overrides) > 0) {
    spec <- utils::modifyList(spec, overrides)
  }

  do.call(ggplot2::element_line, spec)
}

# --- Source: ggforge/R/00-types.R ---

#' Variable Type Detection System
#'
#' @description
#' Unified system for detecting and classifying variable types.
#' This eliminates scattered type-checking logic across plot functions.
#'
#' @name variable-types
#' @keywords internal
NULL

# =============================================================================
# TYPE DETECTION
# =============================================================================

#' Detect variable type
#'
#' @description
#' Determines the type of a variable for styling purposes.
#'
#' Type hierarchy:
#' - continuous: numeric (not factor)
#' - discrete: factor, character, logical
#' - temporal: Date, POSIXct, POSIXlt
#' - ordered: ordered factor
#'
#' @param data Data frame
#' @param column Character. Column name to check
#' @return Character. One of "continuous", "discrete", "temporal", "ordered"
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   x = 1:10,
#'   y = factor(letters[1:10]),
#'   z = Sys.Date() + 1:10
#' )
#' detect_var_type(data, "x") # "continuous"
#' detect_var_type(data, "y") # "discrete"
#' detect_var_type(data, "z") # "temporal"
#' }
detect_var_type <- function(data, column) {
  if (is.null(column)) {
    return(NULL)
  }

  # Handle internal columns (like ".y", ".density")
  if (startsWith(column, ".")) {
    return("continuous")
  }

  col_data <- data[[column]]

  # Check for custom type attribute
  custom_type <- attr(col_data, "ggforge_type", exact = TRUE)
  if (!is.null(custom_type)) {
    return(custom_type)
  }

  # Type detection hierarchy
  if (is.ordered(col_data)) {
    return("ordered")
  }

  if (is.factor(col_data) || is.character(col_data) || is.logical(col_data)) {
    return("discrete")
  }

  if (inherits(col_data, c("Date", "POSIXct", "POSIXlt"))) {
    return("temporal")
  }

  if (is.numeric(col_data)) {
    return("continuous")
  }

  # Fallback
  return("discrete")
}

#' Check if variable is continuous
#'
#' @description
#' Helper to check if a variable should use continuous styling.
#'
#' @param data Data frame
#' @param column Character. Column name
#' @return Logical
#' @keywords internal
is_continuous_var <- function(data, column) {
  if (is.null(column)) {
    return(FALSE)
  }
  type <- detect_var_type(data, column)
  type %in% c("continuous", "temporal")
}

#' Check if variable is discrete
#'
#' @description
#' Helper to check if a variable should use discrete styling.
#'
#' @param data Data frame
#' @param column Character. Column name
#' @return Logical
#' @keywords internal
is_discrete_var <- function(data, column) {
  if (is.null(column)) {
    return(FALSE)
  }
  type <- detect_var_type(data, column)
  type %in% c("discrete", "ordered")
}

#' Annotate variable with explicit type
#'
#' @description
#' Allows users to override automatic type detection.
#'
#' @param x Vector. The variable to annotate
#' @param type Character. One of "continuous", "discrete", "temporal", "ordered"
#' @return The vector with type attribute attached
# @keywords internal
#'
#' @examples
#' # Force a numeric year column to be treated as discrete
#' data <- data.frame(year = 2020:2025, value = rnorm(6))
#' data$year <- var_type(data$year, "discrete")
#'
#' # Force a character column to be treated as continuous (rare)
#' data$score <- var_type(c(1.5, 2.3, 3.1), "continuous")
var_type <- function(x, type = c("continuous", "discrete", "temporal", "ordered")) {
  type <- match.arg(type)
  attr(x, "ggforge_type") <- type
  x
}

# --- Source: ggforge/R/01-validator.R ---

#' Parameter Validation System
#'
#' @description
#' A robust validation system for plot parameters with clear error messages
#' and automatic type coercion where appropriate.
#'
#' @name validator
#' @keywords internal
NULL

#' Column Validator
#'
#' @description
#' Validates that specified columns exist in a data frame and optionally
#' converts them to factors with proper level ordering.
#'
#' @param data A data frame
#' @param columns Column name(s) to validate
#' @param force_factor Whether to convert to factor
#' @param allow_multi Whether to allow multiple columns
#' @param concat_multi Whether to concatenate multiple columns
#' @param concat_sep Separator for concatenation
#' @return Validated column name(s)
#' @keywords internal
#' @importFrom tidyr unite
#' @importFrom rlang sym syms
#' @importFrom tidyr expand_grid
validate_columns <- function(
    data,
    columns,
    force_factor = FALSE,
    allow_multi = FALSE,
    concat_multi = FALSE,
    concat_sep = "_") {
  if (is.null(data)) {
    stop("Data is NULL", call. = FALSE)
  }

  if (is.null(columns)) {
    return(NULL)
  }

  # Get calling context for better error messages
  param_name <- deparse(substitute(columns))
  df_name <- deparse(substitute(data))

  # Check if multiple columns are allowed
  if (!allow_multi && length(columns) > 1) {
    stop(
      sprintf("Only one column allowed in '%s', got %d", param_name, length(columns)),
      call. = FALSE
    )
  }

  # Check all columns exist
  missing_cols <- setdiff(columns, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Column(s) not found in data: %s",
        paste0("'", missing_cols, "'", collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Handle multiple column concatenation
  data_modified <- FALSE
  if (allow_multi && concat_multi && length(columns) > 1) {
    message(
      sprintf(
        "Concatenating %d columns in '%s' with separator '%s'",
        length(columns), param_name, concat_sep
      )
    )

    new_col <- paste(columns, collapse = concat_sep)
    data <- unite(data, !!sym(new_col), !!!syms(columns), sep = concat_sep, remove = FALSE)

    # Preserve factor level ordering if needed
    if (force_factor) {
      all_levels <- lapply(columns, function(col) {
        if (is.factor(data[[col]])) {
          levels(data[[col]])
        } else {
          unique(data[[col]])
        }
      })
      all_levels <- do.call(expand_grid, all_levels)
      all_levels <- apply(all_levels, 1, paste, collapse = concat_sep)
      data[[new_col]] <- droplevels(factor(data[[new_col]], levels = unique(all_levels)))
    }

    columns <- new_col
    data_modified <- TRUE
  }

  # Force to factor if requested
  if (force_factor) {
    for (col in columns) {
      if (!is.factor(data[[col]])) {
        data[[col]] <- factor(data[[col]], levels = unique(data[[col]]))
        data_modified <- TRUE
      }
    }
  }

  # Propagate data modifications back to the caller via parent.frame().
  #
  # NOTE: We deliberately do NOT attach .updated_data as an attribute on the
  # returned column names. Attributed strings break tibble's [[<- indexing
  # (e.g., data[[x]] <- value fails when x has attributes). Instead we rely
  # on parent.frame() for data propagation, which works for all callers that
  # use simple variable names for the data argument.
  if (data_modified) {
    if (grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", df_name)) {
      parent_env <- parent.frame()
      if (exists(df_name, envir = parent_env, inherits = FALSE)) {
        parent_env[[df_name]] <- data
      }
    }
  }

  return(columns)
}

#' Validate and normalize common plot arguments
#'
#' @description
#' Validates common arguments used across all plot types
#'
#' @param seed Random seed
#' @param facet_by Faceting columns
#' @param split_by Split columns
#' @param group_by Grouping columns
#' @param facet_scales Facet scales type
#' @param theme Theme name or function
#' @param palette Palette name
#' @param alpha Transparency value
#' @param aspect.ratio Aspect ratio
#' @param legend.position Legend position
#' @param legend.direction Legend direction
#' @param ... Additional arguments (ignored)
#' @return Invisible NULL (sets seed as side effect)
#' @keywords internal
validate_common_args <- function(
    seed = 8525,
    facet_by = NULL,
    split_by = NULL,
    group_by = NULL,
    facet_scales = "fixed",
    theme = "theme_ggforge",
    palette = "Paired",
    alpha = 1,
    aspect.ratio = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    ...) {
  # Validate seed
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("'seed' must be a single numeric value", call. = FALSE)
  }
  set.seed(seed)

  # Validate facet_by
  if (length(facet_by) > 2) {
    stop("Maximum 2 columns allowed in 'facet_by'", call. = FALSE)
  }

  # Validate facet_scales
  valid_scales <- c("fixed", "free", "free_x", "free_y")
  if (!facet_scales %in% valid_scales) {
    stop(
      sprintf(
        "'facet_scales' must be one of: %s",
        paste(valid_scales, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Validate alpha
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a numeric value between 0 and 1", call. = FALSE)
  }

  # Validate aspect.ratio
  if (!is.null(aspect.ratio) && (!is.numeric(aspect.ratio) || aspect.ratio <= 0)) {
    stop("'aspect.ratio' must be a positive numeric value", call. = FALSE)
  }

  # Validate legend.position (handle vectors for split_by)
  valid_positions <- c("none", "left", "right", "bottom", "top")
  if (is.character(legend.position) && length(legend.position) == 1) {
    if (!legend.position %in% valid_positions) {
      stop(
        sprintf(
          "'legend.position' must be one of: %s",
          paste(valid_positions, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else if (is.character(legend.position) && length(legend.position) > 1) {
    # Validate each element for vectors (used with split_by)
    invalid <- legend.position[!legend.position %in% valid_positions]
    if (length(invalid) > 0) {
      stop(
        sprintf(
          "Invalid 'legend.position' value(s): %s. Must be one of: %s",
          paste(invalid, collapse = ", "),
          paste(valid_positions, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # Validate legend.direction (handle vectors for split_by)
  valid_directions <- c("horizontal", "vertical")
  if (length(legend.direction) == 1) {
    if (!legend.direction %in% valid_directions) {
      stop(
        sprintf(
          "'legend.direction' must be one of: %s",
          paste(valid_directions, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else if (length(legend.direction) > 1) {
    # Validate each element for vectors (used with split_by)
    invalid <- legend.direction[!legend.direction %in% valid_directions]
    if (length(invalid) > 0) {
      stop(
        sprintf(
          "Invalid 'legend.direction' value(s): %s. Must be one of: %s",
          paste(invalid, collapse = ", "),
          paste(valid_directions, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

#' Validate numeric range
#'
#' @param value Value to validate
#' @param param_name Parameter name for error message
#' @param min_val Minimum allowed value
#' @param max_val Maximum allowed value
#' @param allow_null Whether NULL is allowed
#' @return The validated value
#' @keywords internal
validate_numeric_range <- function(
    value,
    param_name,
    min_val = -Inf,
    max_val = Inf,
    allow_null = TRUE) {
  if (is.null(value)) {
    if (allow_null) {
      return(NULL)
    } else {
      stop(sprintf("'%s' cannot be NULL", param_name), call. = FALSE)
    }
  }

  if (!is.numeric(value)) {
    stop(sprintf("'%s' must be numeric", param_name), call. = FALSE)
  }

  if (any(value < min_val | value > max_val, na.rm = TRUE)) {
    stop(
      sprintf("'%s' must be between %s and %s", param_name, min_val, max_val),
      call. = FALSE
    )
  }

  return(value)
}

#' Validate choice from options
#'
#' @param value Value to validate
#' @param param_name Parameter name
#' @param choices Valid choices
#' @param allow_null Whether NULL is allowed
#' @return The validated value
#' @keywords internal
validate_choice <- function(
    value,
    param_name,
    choices,
    allow_null = TRUE) {
  if (is.null(value)) {
    if (allow_null) {
      return(NULL)
    } else {
      stop(sprintf("'%s' cannot be NULL", param_name), call. = FALSE)
    }
  }

  if (!value %in% choices) {
    stop(
      sprintf(
        "'%s' must be one of: %s",
        param_name,
        paste(choices, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  return(value)
}

# --- Source: ggforge/R/02-parameters.R ---

#' Common Plot Parameters Documentation
#'
#' @description
#' This file documents common parameters used across multiple plot functions.
#' By centralizing documentation, we maintain consistency and reduce duplication.
#'
#' @name parameters
#' @param data A data frame containing the data to plot
#' @param x Column name for x-axis variable
#' @param y Column name for y-axis variable
#' @param group_by Column name(s) for grouping data
#' @param group_by_sep Separator when concatenating multiple group_by columns
#' @param split_by Column name(s) to split data into multiple plots
#' @param split_by_sep Separator when concatenating multiple split_by columns
#' @param facet_by Column name(s) for faceting the plot
#' @param facet_scales Scales for facets: "fixed", "free", "free_x", "free_y"
#' @param facet_nrow Number of rows in facet layout
#' @param facet_ncol Number of columns in facet layout
#' @param facet_byrow Fill facets by row (TRUE) or column (FALSE)
#' @param keep_empty Keep empty factor levels
#' @param theme Theme name (string) or theme function
#' @param theme_args List of arguments passed to theme function
#' @param palette Color palette name
#' @param palcolor Custom colors for palette
#' @param alpha Transparency level (0-1)
#' @param x_text_angle Angle for x-axis text labels
#' @param aspect.ratio Aspect ratio of plot panel
#' @param legend.position Legend position: "none", "left", "right", "bottom", "top"
#' @param legend.direction Legend direction: "horizontal" or "vertical"
#' @param expand Expansion values for plot axes (CSS-like padding)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param combine Whether to combine split plots into one
#' @param nrow Number of rows when combining plots
#' @param ncol Number of columns when combining plots
#' @param byrow Fill combined plots by row
#' @param axes How to handle axes in combined plots ("keep", "collect", "collect_x", "collect_y")
#' @param axis_titles How to handle axis titles in combined plots
#' @param guides How to handle guides in combined plots ("collect", "keep", "auto")
#' @param design Custom layout design for combined plots
#' @param seed Random seed for reproducibility
#' @keywords internal
NULL

#' Plot Parameter Class
#'
#' @description
#' S3 class for storing and validating plot parameters
#'
#' @param ... Named parameters
#' @return An object of class "plot_params"
#' @keywords internal
new_plot_params <- function(...) {
  params <- list(...)
  class(params) <- c("plot_params", "list")
  params
}

# Note: merge_params() and get_default_params() were part of the original
# architecture plan but are not used - all plot functions use
# as.list(environment()) to capture parameters directly.

# --- Source: ggforge/R/03-themes.R ---

#' ggforge Theming System
#'
#' @description
#' A flexible and elegant theming system for ggforge plots
#'
#' @name themes
NULL

#' Build common ggforge theme arguments
#'
#' @description
#' Internal helper that generates the shared theme argument list used by
#' both \code{theme_ggforge} and \code{theme_ggforge_grid}.
#'
#' @param aspect.ratio Aspect ratio
#' @param text_size_scale Font size multiplier
#' @param font_family Font family
#' @param legend_text_size Legend text size (before scaling)
#' @return A named list of theme arguments
#' @keywords internal
.build_common_theme_args <- function(aspect.ratio, text_size_scale, font_family,
                                     legend_text_size = 10) {
  list(
    aspect.ratio = aspect.ratio,
    text = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      color = "black"
    ),
    plot.title = ggplot2::element_text(
      size = 13 * text_size_scale,
      family = font_family,
      colour = "black",
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      hjust = 0.5,
      margin = ggplot2::margin(b = 3),
      colour = "black"
    ),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(10, 10, 10, 10),
    axis.title = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      colour = "black",
      face = "bold"
    ),
    axis.text = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      colour = "black",
      face = "plain"
    ),
    strip.text = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      colour = "black",
      hjust = 0.5,
      face = "bold",
      margin = ggplot2::margin(3, 3, 3, 3)
    ),
    strip.background = ggplot2::element_rect(fill = "grey90", colour = "black", linewidth = 0.5),
    legend.title = ggplot2::element_text(
      size = 12 * text_size_scale,
      family = font_family,
      colour = "black",
      hjust = 0,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = legend_text_size * text_size_scale,
      family = font_family,
      colour = "black"
    ),
    legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
    legend.key.size = grid::unit(12, "pt"),
    legend.key.spacing.y = grid::unit(2, "pt"),
    legend.spacing.y = grid::unit(1, "pt"),
    legend.background = ggplot2::element_blank()
  )
}

#' Build a ggforge theme from base + overrides
#'
#' @description
#' Internal helper that combines a base theme with common args and user overrides.
#'
#' @param base_theme A ggplot2 base theme
#' @param theme_args Common theme arguments (from \code{.build_common_theme_args})
#' @param ... User overrides
#' @return A ggplot2 theme object
#' @keywords internal
.apply_theme_overrides <- function(base_theme, theme_args, ...) {
  user_args <- list(...)
  theme_args <- modifyList(theme_args, user_args)
  valid_args <- theme_args[names(theme_args) %in% methods::formalArgs(ggplot2::theme)]
  base_theme + do.call(ggplot2::theme, valid_args)
}

#' Main ggforge Theme
#'
#' @description
#' The default theme for ggforge, providing a clean and modern appearance
#'
#' @param aspect.ratio Aspect ratio of the plot panel
#' @param base_size Base font size (scales all text elements)
#' @param font_family Font family for all text
#' @param ... Additional arguments passed to \code{\link[ggplot2]{theme}}
#' @return A ggplot2 theme object
# @keywords internal
#' @importFrom ggplot2 theme element_text element_rect element_blank element_line margin unit
#' @importFrom methods formalArgs
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_ggforge()
theme_ggforge <- function(
    aspect.ratio = NULL,
    base_size = NULL,
    font_family = NULL,
    ...) {
  base_size <- base_size %||% ggforge_option("theme.base_size")
  font_family <- font_family %||% ggforge_option("theme.font_family")
  text_size_scale <- base_size / 12

  base_theme <- ggplot2::theme_classic(base_line_size = 0.6)
  theme_args <- .build_common_theme_args(aspect.ratio, text_size_scale, font_family,
    legend_text_size = 10
  )
  .apply_theme_overrides(base_theme, theme_args, ...)
}

#' ggforge Theme for Grid-based Plots (Heatmap-like)
#'
#' @description
#' A theme based on theme_bw, designed for grid-based plots like heatmaps,
#' dot plots, and correlation matrices
#'
#' @param aspect.ratio Aspect ratio of the plot panel
#' @param base_size Base font size (scales all text elements)
#' @param font_family Font family for all text
#' @param ... Additional arguments passed to \code{\link[ggplot2]{theme}}
#' @return A ggplot2 theme object
# @keywords internal
#' @importFrom ggplot2 theme_bw theme element_text element_rect element_blank element_line margin unit
#' @importFrom methods formalArgs
theme_ggforge_grid <- function(
    aspect.ratio = NULL,
    base_size = NULL,
    font_family = NULL,
    ...) {
  base_size <- base_size %||% ggforge_option("theme.base_size")
  font_family <- font_family %||% ggforge_option("theme.font_family")
  text_size_scale <- base_size / 12

  base_theme <- ggplot2::theme_bw(base_line_size = 0.6, base_rect_size = 1.2)
  theme_args <- .build_common_theme_args(aspect.ratio, text_size_scale, font_family,
    legend_text_size = 11
  )
  .apply_theme_overrides(base_theme, theme_args, ...)
}

#' Minimal Theme
#'
#' @description
#' A minimal theme with coordinate axes
#'
#' @param add_coord Whether to add coordinate arrows
#' @param xlen_npc Length of x-axis arrow (in npc units)
#' @param ylen_npc Length of y-axis arrow (in npc units)
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param lab_size Label size
#' @param ... Additional arguments passed to theme
#' @return List of ggplot2 theme components
# @keywords internal
#' @importFrom ggplot2 theme element_blank margin annotation_custom coord_cartesian
#' @importFrom grid grobTree gList linesGrob textGrob arrow gpar unit
theme_minimal_axes <- function(
    add_coord = TRUE,
    xlen_npc = 0.15,
    ylen_npc = 0.15,
    xlab = "",
    ylab = "",
    lab_size = 12,
    ...) {
  theme_args <- list(
    panel.border = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key.size = grid::unit(15, "pt"),
    legend.key.spacing.y = grid::unit(4, "pt"),
    legend.spacing.y = grid::unit(2, "pt"),
    plot.margin = ggplot2::margin(
      lab_size + 2,
      lab_size + 2,
      lab_size + 2,
      lab_size + 2,
      unit = "points"
    )
  )

  # Merge with user args
  user_args <- list(...)
  theme_args <- modifyList(theme_args, user_args)
  valid_args <- theme_args[names(theme_args) %in% methods::formalArgs(ggplot2::theme)]

  out <- do.call(ggplot2::theme, valid_args)

  if (add_coord) {
    g <- grid::grobTree(grid::gList(
      grid::linesGrob(
        x = grid::unit(c(0, xlen_npc), "npc"),
        y = grid::unit(c(0, 0), "npc"),
        arrow = grid::arrow(length = grid::unit(0.02, "npc")),
        gp = grid::gpar(lwd = 2)
      ),
      grid::textGrob(
        label = xlab,
        x = grid::unit(0, "npc"),
        y = grid::unit(0, "npc"),
        vjust = 4 / 3,
        hjust = 0,
        gp = grid::gpar(fontsize = lab_size)
      ),
      grid::linesGrob(
        x = grid::unit(c(0, 0), "npc"),
        y = grid::unit(c(0, ylen_npc), "npc"),
        arrow = grid::arrow(length = grid::unit(0.02, "npc")),
        gp = grid::gpar(lwd = 2)
      ),
      grid::textGrob(
        label = ylab,
        x = grid::unit(0, "npc"),
        y = grid::unit(0, "npc"),
        vjust = -2 / 3,
        hjust = 0,
        rot = 90,
        gp = grid::gpar(fontsize = lab_size)
      )
    ))

    return(list(
      ggplot2::annotation_custom(g),
      theme_ggforge() + out,
      ggplot2::coord_cartesian(clip = "off")
    ))
  } else {
    return(list(theme_ggforge() + out))
  }
}

#' Process theme argument
#'
#' @description
#' Converts theme name string to theme function
#'
#' @param theme Theme name or function
#' @return Theme function
#' @keywords internal
#' @importFrom utils getFromNamespace
process_theme <- function(theme) {
  if (is.function(theme)) {
    return(theme)
  }

  if (!is.character(theme)) {
    stop("Theme must be a character string or function", call. = FALSE)
  }

  # Handle namespace notation (e.g., "ggplot2::theme_minimal")
  if (grepl("::", theme)) {
    parts <- strsplit(theme, "::")[[1]]
    if (length(parts) != 2) {
      stop("Invalid theme specification: ", theme, call. = FALSE)
    }
    return(getFromNamespace(parts[2], parts[1]))
  }

  # Try to get from ggforge namespace first
  if (exists(theme, mode = "function", envir = environment())) {
    return(get(theme, envir = environment()))
  }

  # Try ggplot2 namespace
  if (exists(theme, mode = "function", envir = asNamespace("ggplot2"))) {
    return(get(theme, envir = asNamespace("ggplot2")))
  }

  stop("Cannot find theme: ", theme, call. = FALSE)
}

# =============================================================================
# STYLE APPLICATION SYSTEM
# =============================================================================

#' Apply data-driven styling to plot
#'
#' @description
#' Automatically applies axis and legend styling based on variable types.
#' This is the core function that eliminates the need for manual
#' get_axis_text_size() calls in every plot function.
#'
#' @param plot ggplot object
#' @param data Data frame used in the plot
#' @param x_var Character. X variable name (NULL if not applicable)
#' @param y_var Character. Y variable name (NULL if not applicable)
#' @param flip Logical. Whether axes are flipped
#' @param base_size Numeric. Base font size
#' @param ... Additional theme arguments to merge
#' @return ggplot object with styling applied
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = factor(cyl))) +
#'   ggplot2::geom_point()
#'
#' # Automatic styling based on variable types
#' p <- apply_style_theme(p, mtcars, x_var = "mpg", y_var = "cyl")
#' }
apply_style_theme <- function(
    plot,
    data,
    x_var = NULL,
    y_var = NULL,
    flip = FALSE,
    base_size = 12,
    ...) {
  # Detect variable types
  x_type <- detect_var_type(data, x_var)
  y_type <- detect_var_type(data, y_var)

  # Map types to style specs
  x_spec_path <- map_type_to_style(x_type)
  y_spec_path <- map_type_to_style(y_type)

  # Handle flip
  if (flip) {
    temp <- x_spec_path
    x_spec_path <- y_spec_path
    y_spec_path <- temp
  }

  # Build theme elements
  theme_updates <- list()

  if (!is.null(x_spec_path)) {
    theme_updates$axis.text.x <- build_element_text(x_spec_path, base_size)
  }

  if (!is.null(y_spec_path)) {
    theme_updates$axis.text.y <- build_element_text(y_spec_path, base_size)
  }

  # Apply legend styling
  theme_updates$legend.title <- build_element_text("font.legend_title", base_size)
  theme_updates$legend.text <- build_element_text("font.legend_text", base_size)

  # Apply standard grid lines (dashed grey) - centralized here to avoid
  # duplicating 'panel.grid.major = element_line(...)' in every atomic function
  theme_updates$panel.grid.major <- build_element_line("lines.grid_major")

  # Merge with additional arguments (user overrides take precedence)
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    theme_updates <- utils::modifyList(theme_updates, extra_args)
  }

  # Apply to plot
  plot + do.call(ggplot2::theme, theme_updates)
}

#' Map variable type to style specification path
#'
#' @param type Character. Variable type from detect_var_type()
#' @return Character. Style spec path
#' @keywords internal
map_type_to_style <- function(type) {
  if (is.null(type)) {
    return(NULL)
  }

  switch(type,
    continuous = "font.axis_text.continuous",
    discrete = "font.axis_text.discrete",
    ordered = "font.axis_text.discrete",
    temporal = "font.axis_text.continuous",
    NULL
  )
}

# --- Source: ggforge/R/04-palettes.R ---

#' Color Palette System for ggforge
#'
#' @description
#' A comprehensive color palette system with support for both discrete
#' and continuous palettes, custom colors, and intelligent color mapping.
#'
#' @name palettes
NULL

#' Get colors from palette
#'
#' @description
#' Main function for retrieving colors from palettes with intelligent
#' handling of discrete and continuous data.
#'
#' @param x Vector of values to map to colors (character, factor, or numeric)
#' @param n Number of colors for continuous palettes
#' @param palette Name of the palette to use
#' @param palcolor Custom colors (overrides palette)
#' @param type Type of palette: "auto", "discrete", or "continuous"
#' @param keep_names Whether to keep names on color vector
#' @param alpha Transparency level (0-1)
#' @param matched Return colors matched to x values
#' @param reverse Reverse the color order
#' @param NA_keep Include color for NA values
#' @param NA_color Color for NA values
#' @param transparent Use true transparency vs color blending
#' @return Named vector of colors
# @keywords internal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom stats setNames
#' @importFrom rlang "%||%"
#' @importFrom scales alpha
#' @examples
#' # Discrete palette
#' get_palette(c("A", "B", "C"), palette = "Paired")
#'
#' # Continuous palette
#' get_palette(1:100, palette = "Spectral", type = "continuous")
#'
#' # Custom colors
#' get_palette(c("A", "B", "C"), palcolor = c("red", "blue", "green"))
get_palette <- function(
    x,
    n = 100,
    palette = "Paired",
    palcolor = NULL,
    type = "auto",
    keep_names = TRUE,
    alpha = 1,
    matched = FALSE,
    reverse = FALSE,
    NA_keep = FALSE,
    NA_color = "grey80",
    transparent = TRUE) {
  # Access palette_list from package data (works with lazy loading)
  palette_list <- .ggforge_palette_list_fn()

  # Handle missing x
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }

  # Handle case where palette is a named color vector (should be palcolor)
  if (length(palette) > 1 || (!is.null(names(palette)) && length(palette) == 1)) {
    if (is.null(palcolor)) {
      palcolor <- palette
    }
    palette <- "Paired" # Use default palette name
  }

  # Validate palette
  if (length(palette) == 1 && !palette %in% names(palette_list)) {
    stop(
      sprintf(
        "Palette '%s' not found. Use show_palettes() to see available palettes",
        palette
      ),
      call. = FALSE
    )
  }

  # Handle custom colors
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }

  if (is.null(palcolor) || length(palcolor) == 0 || all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }

  # Handle named palcolors that match x values
  if (!is.null(names(palcolor))) {
    matched_colors <- palcolor[intersect(names(palcolor), x)]

    if (length(matched_colors) < length(x) && length(matched_colors) > 0) {
      # Partial match: fill in missing with palette
      palcolor <- get_palette(
        x = x, n = n, palette = palette, palcolor = NULL,
        type = type, keep_names = TRUE, alpha = 1,
        matched = matched, reverse = reverse,
        NA_keep = NA_keep, NA_color = NA_color,
        transparent = transparent
      )
      palcolor[names(matched_colors)] <- matched_colors
      reverse <- FALSE # Already reversed if needed
    } else if (length(matched_colors) == length(x)) {
      palcolor <- matched_colors
    }
  }

  pal_n <- length(palcolor)

  # Determine type
  if (type == "auto") {
    type <- if (is.numeric(x)) "continuous" else "discrete"
  }

  # Validate type
  if (!type %in% c("discrete", "continuous")) {
    stop("'type' must be 'auto', 'discrete', or 'continuous'", call. = FALSE)
  }

  # Generate colors based on type
  if (type == "discrete") {
    color <- .get_discrete_colors(
      x, palcolor, pal_n, matched, NA_color
    )
  } else {
    color <- .get_continuous_colors(
      x, palcolor, pal_n, n, matched, NA_color
    )
  }

  # Reverse if requested
  if (reverse) {
    if (!is.null(names(color))) {
      color <- setNames(rev(color), names(color))
    } else {
      color <- rev(color)
    }
  }

  # Handle NA
  if (!NA_keep) {
    color <- color[names(color) != "NA"]
  }

  # Remove names if requested
  if (!keep_names) {
    names(color) <- NULL
  }

  # Apply alpha
  if (alpha < 1) {
    color <- .apply_alpha(color, alpha, transparent)
  }

  return(color)
}

#' Get discrete colors
#' @keywords internal
.get_discrete_colors <- function(x, palcolor, pal_n, matched, NA_color) {
  if (!is.factor(x)) {
    x <- factor(x, levels = unique(x))
  }

  n_x <- nlevels(x)

  # Check if palette is continuous type
  if (isTRUE(attr(palcolor, "type") == "continuous")) {
    color <- colorRampPalette(palcolor)(n_x)
    names(color) <- levels(x)
  } else if (!is.null(names(palcolor))) {
    # Use named colors
    color <- palcolor[intersect(names(palcolor), levels(x))]
  } else {
    # Generate colors
    if (n_x <= pal_n) {
      color <- palcolor[1:n_x]
    } else {
      color <- colorRampPalette(palcolor)(n_x)
    }
    names(color) <- levels(x)
  }

  # Add NA color if needed
  if (any(is.na(x))) {
    color <- c(color, setNames(NA_color, "NA"))
  }

  # Match to x values
  if (matched) {
    color <- color[as.character(x)]
    color[is.na(color)] <- NA_color
  }

  return(color)
}

#' Get continuous colors
#' @keywords internal
.get_continuous_colors <- function(x, palcolor, pal_n, n, matched, NA_color) {
  if (!is.numeric(x) && all(!is.na(x))) {
    stop("x must be numeric for continuous palettes", call. = FALSE)
  }

  # Handle edge cases
  if (all(is.na(x))) {
    values <- as.factor(rep(0, n))
  } else if (length(unique(na.omit(as.numeric(x)))) == 1) {
    values <- as.factor(rep(unique(na.omit(as.numeric(x))), n))
  } else {
    if (matched) {
      values <- cut(
        x,
        breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1),
        include.lowest = TRUE
      )
    } else {
      values <- cut(
        1:100,
        breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1),
        include.lowest = TRUE
      )
    }
  }

  n_x <- nlevels(values)

  # Generate colors
  if (n_x <= pal_n) {
    color <- palcolor[1:n_x]
  } else {
    color <- colorRampPalette(palcolor)(n_x)
  }
  names(color) <- levels(values)

  # Add NA color
  if (any(is.na(x))) {
    color <- c(color, setNames(NA_color, "NA"))
  }

  # Match to values
  if (matched) {
    if (all(is.na(x))) {
      color <- NA_color
    } else if (length(unique(na.omit(x))) == 1) {
      color <- color[as.character(unique(na.omit(x)))]
      color[is.na(color)] <- NA_color
    } else {
      color <- color[as.character(values)]
      color[is.na(color)] <- NA_color
    }
  }

  return(color)
}

#' Apply alpha to colors
#' @keywords internal
#' @importFrom grDevices col2rgb rgb
.apply_alpha <- function(colors, alpha_val, transparent) {
  if (transparent) {
    # Use true transparency
    return(scales::alpha(colors, alpha_val))
  } else {
    # Blend with white
    has_names <- !is.null(names(colors))
    color_df <- as.data.frame(col2rgb(colors) / 255)

    colors_out <- sapply(color_df, function(color) {
      # Blend with white background
      color_rgb <- color * alpha_val + c(1, 1, 1) * (1 - alpha_val)
      rgb(color_rgb[1], color_rgb[2], color_rgb[3])
    })

    if (has_names) {
      names(colors_out) <- names(colors)
    }

    return(colors_out)
  }
}

#' Show available palettes
#'
#' @description
#' Display available color palettes visually
#'
#' @param palettes Custom palette list (NULL to use built-in)
#' @param type Type of palettes to show: "discrete", "continuous", or both
#' @param index Indices of palettes to show
#' @param palette_names Specific palette names to show
#' @param return_names Return palette names instead of plotting
#' @param return_palettes Return palette colors instead of plotting
#' @return Plot, palette names, or palette colors
# @keywords internal
#' @importFrom ggplot2 ggplot geom_col scale_fill_manual scale_x_continuous element_blank aes
show_palettes <- function(
    palettes = NULL,
    type = c("discrete", "continuous"),
    index = NULL,
    palette_names = NULL,
    return_names = TRUE,
    return_palettes = FALSE) {
  # Get palette list
  if (!is.null(palettes)) {
    palette_list <- palettes
  } else {
    # Access palette_list from package data (works with lazy loading)
    palette_list <- .ggforge_palette_list_fn()
    palette_list <- palette_list[
      unlist(lapply(palette_list, function(x) {
        isTRUE(attr(x, "type") %in% type)
      }))
    ]
  }

  # Filter by index
  if (!is.null(index)) {
    index <- index[index %in% seq_along(palette_list)]
    palette_list <- palette_list[index]
  }

  # Set names if missing
  if (is.null(names(palette_list))) {
    names(palette_list) <- seq_along(palette_list)
  }

  # Filter by names
  if (!is.null(palette_names)) {
    missing_names <- setdiff(palette_names, names(palette_list))
    if (length(missing_names) > 0) {
      stop("Cannot find palettes: ", paste(missing_names, collapse = ", "), call. = FALSE)
    }
    palette_list <- palette_list[palette_names]
  }

  # Return options
  if (return_palettes) {
    return(palette_list)
  }
  if (return_names) {
    return(names(palette_list))
  }

  # Create visualization
  df <- data.frame(
    palette = rep(names(palette_list), sapply(palette_list, length)),
    color = unlist(palette_list),
    stringsAsFactors = FALSE
  )
  df$palette <- factor(df$palette, levels = rev(unique(df$palette)))
  df$color_order <- factor(seq_len(nrow(df)), levels = seq_len(nrow(df)))
  df$proportion <- as.numeric(1 / table(df$palette)[df$palette])

  p <- ggplot2::ggplot(df, aes(y = .data$palette, x = .data$proportion, fill = .data$color_order)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = df$color) +
    scale_x_continuous(expand = c(0, 0), trans = "reverse") +
    theme_ggforge(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank()
    )

  print(p)
  invisible(names(palette_list))
}

#' Check and normalize palette argument
#' @keywords internal
check_palette <- function(palette, split_names) {
  palette <- as.list(palette)

  if (length(palette) == 0) {
    stop("'palette' must be specified", call. = FALSE)
  }

  # Special case: no split_by (split_names = "...")
  # If palette has names and doesn't match "...", it's meant for group_by
  # So just wrap it in a list and return
  if (length(split_names) == 1 && split_names[1] == "...") {
    if (!is.null(names(palette)) && !"..." %in% names(palette)) {
      # Named palette for group_by, wrap it
      result <- list(palette)
      names(result) <- "..."
      return(result)
    }
  }

  # Handle named palettes (partial naming allowed)
  if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
    missing <- setdiff(split_names, names(palette))
    if (length(missing) > 0) {
      # Fill missing values with default palette
      default_palette <- if (is.null(palette[[1]]) || palette[[1]] == "") "Paired" else palette[[1]]
      for (m in missing) {
        palette[[m]] <- default_palette
      }
      message(sprintf(
        "Using palette '%s' for split_by values without explicit palette: %s",
        default_palette, paste(missing, collapse = ", ")
      ))
    }
    # Reorder to match split_names
    palette <- palette[split_names]
    return(palette)
  }

  # Replicate if needed (unnamed palettes)
  if (length(palette) == 1 && length(split_names) > 1) {
    palette <- rep(palette, length(split_names))
  }

  if (length(palette) < length(split_names)) {
    stop(
      sprintf(
        "Palette length (%d) less than split_by values (%d)",
        length(palette), length(split_names)
      ),
      call. = FALSE
    )
  }

  # Set names
  if (is.null(names(palette))) {
    names(palette)[seq_along(split_names)] <- split_names
  }

  return(palette)
}

#' Check and normalize palcolor argument
#' @keywords internal
check_palcolor <- function(palcolor, split_names) {
  if (is.null(palcolor)) {
    return(NULL)
  }

  # Convert to list if needed
  if (!is.list(palcolor)) {
    palcolor <- list(palcolor)
  }

  # Special case: no split_by (split_names = "...")
  # If palcolor has names and doesn't match "...", it's meant for group_by
  # So just wrap it in a list and return
  if (length(split_names) == 1 && split_names[1] == "...") {
    if (!is.null(names(palcolor)) && length(palcolor) > 1 && !"..." %in% names(palcolor)) {
      # Named palcolor for group_by, wrap it
      result <- list(palcolor)
      names(result) <- "..."
      return(result)
    }
  }

  # Handle special case
  if (identical(split_names, "...") && !identical(names(palcolor), "...")) {
    names(palcolor) <- split_names
  }

  # Replicate if needed
  if (length(palcolor) == 1 && length(split_names) > 1) {
    palcolor <- rep(palcolor, length(split_names))
  }

  # Set names
  if (is.null(names(palcolor))) {
    names(palcolor)[seq_along(split_names)] <- split_names
  }

  return(palcolor)
}

# --- Source: ggforge/R/05-utils.R ---

#' Utility Functions for ggforge
#'
#' @description
#' Helper functions for data manipulation, layout, and plotting
#'
#' @name utils
#' @keywords internal
NULL

# NOTE: get_axis_text_size() has been replaced by apply_style_theme()
# See R/03-themes.R for the new implementation based on the style specs system

#' Apply column validation updates to data
#'
#' @description
#' Kept for backward compatibility. Data propagation from validate_columns()
#' is now handled entirely via parent.frame(). This function simply returns
#' the data unchanged.
#'
#' @param data The data frame
#' @param validated_cols The return value from \code{validate_columns()} (unused)
#' @return The data frame (unchanged)
#' @keywords internal
apply_column_update <- function(data, validated_cols) {
  data
}

#' Check linewidth parameter (compatibility shim for qqplotr)
#'
#' @description
#' This is a compatibility shim for the qqplotr package which expects
#' the check_linewidth() function from older versions of ggplot2.
#' This function handles the conversion of the old 'size' aesthetic to
#' the new 'linewidth' aesthetic.
#'
#' @param data Data frame containing aesthetic data
#' @param snake_class Class name in snake_case (unused, for compatibility)
#' @return Data frame with size converted to linewidth if needed
#' @keywords internal
check_linewidth <- function(data, snake_class = NULL) {
  # This is a compatibility shim for the old ggplot2 internal function
  # that was used to convert 'size' to 'linewidth'

  # If data has a 'size' column but no 'linewidth', convert it
  if (!is.null(data$size) && is.null(data$linewidth)) {
    data$linewidth <- data$size
  }

  # If data has a 'linewidth' column, validate it
  if (!is.null(data$linewidth)) {
    if (!is.numeric(data$linewidth)) {
      # Try to coerce to numeric
      data$linewidth <- as.numeric(data$linewidth)
      if (all(is.na(data$linewidth))) {
        stop("linewidth must be numeric", call. = FALSE)
      }
    }

    if (any(data$linewidth < 0, na.rm = TRUE)) {
      warning("linewidth values should be non-negative", call. = FALSE)
    }
  }

  return(data)
}

#' Get ggplot function (with gglogger support)
#'
#' @description
#' Centralized function to get the ggplot function with optional gglogger support.
#' This eliminates duplicate code across all plot functions.
#'
#' @return ggplot function (either from ggplot2 or gglogger)
#' @keywords internal
#' @importFrom utils getFromNamespace
get_ggplot <- function() {
  if (getOption("ggforge.gglogger.enabled", FALSE) &&
    requireNamespace("gglogger", quietly = TRUE)) {
    getFromNamespace("ggplot", "gglogger")
  } else {
    ggplot2::ggplot
  }
}

#' Normalize legend position (handle waiver)
#'
#' @description
#' Handles waiver objects for legend.position, converting them to actual positions
#' based on whether a grouping variable exists.
#'
#' @param position Legend position (can be waiver or string)
#' @param has_group Logical, whether a grouping variable exists
#' @param default_with_group Default position when group exists
#' @param default_no_group Default position when no group
#' @return Normalized position string
#' @keywords internal
#' @importFrom ggplot2 waiver
normalize_legend_position <- function(
    position,
    has_group,
    default_with_group = "right",
    default_no_group = "none") {
  if (inherits(position, "waiver")) {
    if (has_group) default_with_group else default_no_group
  } else {
    position
  }
}

#' Calculate text justification based on angle
#'
#' @param angle Text angle in degrees
#' @return List with h (hjust) and v (vjust) values
#' @keywords internal
calc_justification <- function(angle) {
  angle <- angle %% 360
  if (angle < 0) angle <- angle + 360

  if (angle < 10) {
    list(h = 0.5, v = 1)
  } else if (angle < 90) {
    list(h = 1, v = 1)
  } else if (angle < 180) {
    list(h = 1, v = 0.5)
  } else if (angle < 270) {
    list(h = 0, v = 0)
  } else if (angle < 315) {
    list(h = 0, v = 0.5)
  } else {
    list(h = 0, v = 1)
  }
}

#' Normalize expansion values (CSS-like padding)
#'
#' @param expand Numeric vector of expansion values
#' @param x_type Type of x-axis ("continuous" or "discrete")
#' @param y_type Type of y-axis ("continuous" or "discrete")
#' @param continuous_default Default for continuous axes
#' @param discrete_default Default for discrete axes
#' @return List with x and y expansion values
#' @keywords internal
#' @importFrom ggplot2 expansion
normalize_expansion <- function(
    expand,
    x_type,
    y_type,
    continuous_default = c(0.05, 0),
    discrete_default = c(0, 0.6)) {
  # Helper to get expansion by type
  expand_by_type <- function(ex, type, both = FALSE) {
    if (type == "continuous") {
      ret <- if (!is.null(ex)) c(ex, 0) else continuous_default
    } else {
      ret <- if (!is.null(ex)) c(0, ex) else discrete_default
    }
    if (both) c(ret, ret) else ret
  }

  if (is.null(expand)) {
    return(list(
      x = expand_by_type(NULL, x_type, both = TRUE),
      y = expand_by_type(NULL, y_type, both = TRUE)
    ))
  }

  # Parse CSS-like notation
  if (is.null(names(expand))) {
    expand <- switch(as.character(length(expand)),
      "1" = c(top = expand, right = expand, bottom = expand, left = expand),
      "2" = c(top = expand[1], right = expand[2], bottom = expand[1], left = expand[2]),
      "3" = c(top = expand[1], right = expand[2], bottom = expand[3], left = expand[2]),
      "4" = c(top = expand[1], right = expand[2], bottom = expand[3], left = expand[4]),
      stop("Invalid expand length", call. = FALSE)
    )
  }

  expand <- as.list(expand)

  # Validate conflicts
  if ("x" %in% names(expand) && any(c("left", "right") %in% names(expand))) {
    stop("Cannot have both 'x' and 'left'/'right' in expand", call. = FALSE)
  }
  if ("y" %in% names(expand) && any(c("top", "bottom") %in% names(expand))) {
    stop("Cannot have both 'y' and 'top'/'bottom' in expand", call. = FALSE)
  }

  # Expand x/y notation
  if ("x" %in% names(expand)) {
    expand$left <- expand$right <- expand$x
    expand$x <- NULL
  }
  if ("y" %in% names(expand)) {
    expand$bottom <- expand$top <- expand$y
    expand$y <- NULL
  }

  list(
    x = c(expand_by_type(expand$left, x_type), expand_by_type(expand$right, x_type)),
    y = c(expand_by_type(expand$bottom, y_type), expand_by_type(expand$top, y_type))
  )
}

#' Add faceting to a plot
#'
#' @param plot ggplot object
#' @param facet_by Column(s) to facet by
#' @param facet_scales Scale type for facets
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param byrow Fill by row
#' @param ... Additional arguments for facet functions
#' @return Faceted plot
#' @keywords internal
#' @importFrom ggplot2 facet_wrap facet_grid vars
#' @importFrom rlang sym
add_facets <- function(
    plot,
    facet_by,
    facet_scales,
    nrow,
    ncol,
    byrow,
    ...) {
  if (is.null(facet_by)) {
    return(plot)
  }

  if (length(facet_by) == 1) {
    plot + facet_wrap(
      facets = facet_by,
      scales = facet_scales,
      nrow = nrow,
      ncol = ncol,
      dir = if (byrow) "h" else "v",
      ...
    )
  } else {
    args <- list(...)
    args$strip.position <- NULL
    args$rows <- vars(!!sym(facet_by[1]))
    args$cols <- vars(!!sym(facet_by[2]))
    args$scales <- facet_scales

    plot + do.call(facet_grid, args)
  }
}

#' Combine multiple plots
#'
#' @param plots List of plot objects
#' @param combine Whether to combine
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param byrow Fill by row
#' @param axes Axis handling
#' @param axis_titles Axis title handling
#' @param guides Guide handling
#' @param design Custom design
#' @return Combined plot or list
#' @keywords internal
#' @importFrom patchwork wrap_plots
combine_plots <- function(
    plots,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = NULL,
    axes = NULL,
    axis_titles = NULL,
    guides = NULL,
    design = NULL) {
  if (!combine) {
    return(plots)
  }

  # Single plot - just return it
  if (length(plots) == 1 && !inherits(plots[[1]], "gTree")) {
    return(plots[[1]])
  }

  # Combine using patchwork
  wrap_plots(
    plots,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    axes = axes,
    axis_titles = axis_titles,
    guides = guides,
    design = design
  )
}

#' Split data by column
#'
#' @param data Data frame
#' @param split_by Column to split by
#' @return Named list of data frames
#' @keywords internal
split_data <- function(data, split_by) {
  if (is.null(split_by)) {
    return(list("..." = data))
  }

  split(data, data[[split_by]])
}

#' Create background layer for plots
#'
#' @param data Data frame
#' @param x X column
#' @param palette Palette name
#' @param palcolor Custom colors
#' @param alpha Alpha value
#' @param keep_empty Keep empty levels
#' @param facet_by Faceting columns
#' @param direction Background direction
#' @return ggplot layer
#' @keywords internal
#' @importFrom ggplot2 geom_rect aes
#' @importFrom dplyr distinct
#' @importFrom tidyr expand_grid
#' @importFrom rlang sym syms
create_bg_layer <- function(
    data,
    x,
    palette,
    palcolor,
    alpha,
    keep_empty,
    facet_by,
    direction = "vertical") {
  fct <- data[[x]]
  if (!keep_empty) {
    fct <- droplevels(fct)
  }

  bg_color <- get_palette(levels(fct), palette = palette, palcolor = palcolor)

  # Create background data
  bg_data <- data.frame(x = factor(levels(fct), levels = levels(fct)))
  bg_data$x_num <- as.numeric(bg_data$x)
  bg_data$xmin <- ifelse(bg_data$x_num == min(bg_data$x_num), -Inf, bg_data$x_num - 0.5)
  bg_data$xmax <- ifelse(bg_data$x_num == max(bg_data$x_num), Inf, bg_data$x_num + 0.5)
  bg_data$ymin <- -Inf
  bg_data$ymax <- Inf
  bg_data$fill <- bg_color[levels(fct)]

  # Add facet columns if needed
  if (!is.null(facet_by)) {
    unique_facets <- distinct(data, !!!syms(facet_by))
    bg_data <- expand_grid(bg_data, unique_facets)
    for (fb in facet_by) {
      bg_data[[fb]] <- factor(bg_data[[fb]], levels = levels(data[[fb]]))
    }
  }

  # Create appropriate geom based on direction
  if (direction == "vertical") {
    geom_rect(
      data = bg_data,
      aes(
        xmin = !!sym("xmin"), xmax = !!sym("xmax"),
        ymin = !!sym("ymin"), ymax = !!sym("ymax")
      ),
      fill = bg_data$fill,
      alpha = alpha,
      inherit.aes = FALSE
    )
  } else {
    geom_rect(
      data = bg_data,
      aes(
        xmin = !!sym("ymin"), xmax = !!sym("ymax"),
        ymin = !!sym("xmin"), ymax = !!sym("xmax")
      ),
      fill = bg_data$fill,
      alpha = alpha,
      inherit.aes = FALSE
    )
  }
}

#' Blend multiple colors
#'
#' @param colors Vector of colors
#' @param mode Blend mode
#' @return Blended color
#' @keywords internal
#' @importFrom grDevices col2rgb rgb
blend_colors <- function(colors, mode = c("blend", "average", "screen", "multiply")) {
  mode <- match.arg(mode)
  # "blend" is a legacy alias for "average"
  if (mode == "blend") mode <- "average"
  colors <- colors[!is.na(colors)]

  if (length(colors) == 0) {
    return(NA)
  }
  if (length(colors) == 1) {
    return(colors)
  }

  rgb_vals <- col2rgb(colors) / 255

  result <- switch(mode,
    average = apply(rgb_vals, 1, mean),
    screen = 1 - apply(1 - rgb_vals, 1, prod),
    multiply = apply(rgb_vals, 1, prod)
  )

  rgb(result[1], result[2], result[3])
}

# =============================================================================
# Text Angle Auto-Detection
# =============================================================================

#' Auto-detect optimal x-axis text angle based on label length
#'
#' @description
#' Automatically determines the best angle for x-axis labels based on
#' the maximum label length, number of labels, and plot orientation.
#' Useful for any plot with a discrete x-axis.
#'
#' @param data A data frame containing the plot data
#' @param x Character string specifying the x-axis column name
#' @param flip Logical; whether the plot coordinates are flipped
#' @param stack Logical; whether facets are stacked
#' @param n_groups Integer; number of groups (for dodged plots)
#'
#' @return Numeric value (0, 45, or 90) representing the optimal text angle
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   short_labels = factor(c("A", "B", "C")),
#'   long_labels = factor(c("Category One", "Category Two", "Category Three"))
#' )
#' auto_detect_text_angle(df, "short_labels", flip = FALSE, stack = FALSE)
#' auto_detect_text_angle(df, "long_labels", flip = FALSE, stack = FALSE)
#' }
#'
#' @keywords internal
auto_detect_text_angle <- function(data, x, flip = FALSE, stack = FALSE, n_groups = 1) {
  x_labels <- levels(data[[x]])
  if (is.null(x_labels)) {
    x_labels <- unique(as.character(data[[x]]))
  }

  max_length <- max(nchar(x_labels), na.rm = TRUE)
  n_labels <- length(x_labels)

  if (flip && stack) {
    if (max_length > 15) return(90)
    else if (max_length > 8) return(45)
    else return(0)
  } else if (flip) {
    if (max_length > 20) return(45)
    else return(0)
  } else {
    if (max_length > 12 || (max_length > 6 && n_labels > 8)) return(45)
    else if (max_length > 8 && n_labels > 5) return(45)
    else return(0)
  }
}

# =============================================================================
# Statistical Comparison Helpers
# =============================================================================

#' Check if a vector has zero variance
#'
#' @description
#' Determines whether a numeric vector has zero variance (all values identical
#' or all NA). Used to prevent statistical test failures.
#'
#' @param vec Numeric vector to check
#' @return Logical; TRUE if vector has zero variance, FALSE otherwise
#' @keywords internal
has_zero_variance <- function(vec) {
  all(is.na(vec)) || length(unique(vec[!is.na(vec)])) <= 1
}

#' Add minimal variance to vector
#'
#' @description
#' Adds tiny perturbations to a constant vector to enable statistical tests.
#'
#' @param yval Numeric vector to modify
#' @param base_value Optional base value for calculating epsilon
#' @return Modified numeric vector with minimal variance added
#' @keywords internal
#' @noRd
add_minimal_variance <- function(yval, base_value = NULL) {
  non_na_idx <- which(!is.na(yval))
  if (length(non_na_idx) < 2) return(yval)

  unique_y <- unique(yval[!is.na(yval)])
  if (length(unique_y) != 1) return(yval)

  if (is.null(base_value)) base_value <- unique_y[1]
  epsilon <- max(abs(base_value) * 1e-10, 1e-10)

  yval[non_na_idx[1]] <- unique_y[1] - epsilon
  yval[non_na_idx[2]] <- unique_y[1] + epsilon
  return(yval)
}

#' Check if data groups have variance issues for group_by comparisons
#' @param data_groups List of data frames, each representing a group
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @param group_by Character string; grouping column name
#' @return Logical; TRUE if variance issues detected
#' @keywords internal
check_variance_issues_grouped <- function(data_groups, x, y, group_by) {
  for (group_data in data_groups) {
    gs <- unique(as.character(group_data[[group_by]]))
    if (length(gs) >= 2) {
      yval1 <- group_data[[y]][group_data[[group_by]] == gs[1]]
      yval2 <- group_data[[y]][group_data[[group_by]] == gs[2]]
      if (has_zero_variance(yval1) && has_zero_variance(yval2)) return(TRUE)
    }
  }
  return(FALSE)
}

#' Check if data groups have variance issues for x comparisons
#' @param data_groups List of data frames, each representing a group
#' @param y Character string; y-axis column name
#' @return Logical; TRUE if variance issues detected
#' @keywords internal
check_variance_issues_simple <- function(data_groups, y) {
  for (group_data in data_groups) {
    if (has_zero_variance(group_data[[y]])) return(TRUE)
  }
  return(FALSE)
}

#' Fix variance issues in facet data for grouped comparisons
#' @param facet_data Data frame containing facet data
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @param group_by Character string; grouping column name
#' @return Modified data frame with variance issues fixed
#' @keywords internal
fix_facet_variance_grouped <- function(facet_data, x, y, group_by) {
  xdata <- split(facet_data, facet_data[[x]])
  all_gs <- unique(as.character(facet_data[[group_by]]))[1:2]

  for (xval in names(xdata)) {
    df <- xdata[[xval]]
    gs <- unique(as.character(df[[group_by]]))

    if (length(gs) < 2) {
      new_df <- data.frame(x_val = xval, y_val = c(0, 1), group_val = all_gs, stringsAsFactors = FALSE)
      colnames(new_df) <- c(x, y, group_by)
      if (is.factor(facet_data[[x]])) new_df[[x]] <- factor(new_df[[x]], levels = levels(facet_data[[x]]))
      if (is.factor(facet_data[[group_by]])) new_df[[group_by]] <- factor(new_df[[group_by]], levels = levels(facet_data[[group_by]]))
      xdata[[xval]] <- new_df
      next
    }

    yval1 <- df[[y]][df[[group_by]] == gs[1]]
    yval2 <- df[[y]][df[[group_by]] == gs[2]]
    if (all(is.na(yval1))) yval1 <- c(0, rep(NA, length(yval1) - 1))
    if (all(is.na(yval2))) yval2 <- c(1, rep(NA, length(yval2) - 1))
    if (has_zero_variance(yval1) && has_zero_variance(yval2)) {
      yval1 <- add_minimal_variance(yval1)
      yval2 <- add_minimal_variance(yval2)
    }
    df[[y]][df[[group_by]] == gs[1]] <- yval1
    df[[y]][df[[group_by]] == gs[2]] <- yval2
    xdata[[xval]] <- df
  }
  do.call(rbind, xdata)
}

#' Fix variance issues in facet data for simple comparisons
#' @param facet_data Data frame containing facet data
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @return Modified data frame with variance issues fixed
#' @keywords internal
fix_facet_variance_simple <- function(facet_data, x, y) {
  xdata <- split(facet_data, facet_data[[x]])
  for (xval in names(xdata)) {
    df <- xdata[[xval]]
    yval <- df[[y]]
    if (all(is.na(yval))) yval <- c(0, 1, rep(NA, max(0, length(yval) - 2)))
    if (has_zero_variance(yval)) yval <- add_minimal_variance(yval)
    df[[y]] <- yval
    xdata[[xval]] <- df
  }
  result <- do.call(rbind, xdata)
  if (is.factor(facet_data[[x]])) result[[x]] <- factor(result[[x]], levels = levels(facet_data[[x]]))
  return(result)
}

#' Preprocess data for pairwise comparisons with group_by
#'
#' @description
#' Prepares data for grouped pairwise statistical comparisons by detecting
#' and fixing variance issues that would cause test failures.
#'
#' @param data Data frame containing the plot data
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @param group_by Character string; grouping column name
#' @param facet_by Character vector; faceting column name(s), or NULL
#' @return Preprocessed data frame ready for statistical comparisons
#' @keywords internal
#' @importFrom dplyr group_by summarise add_count filter mutate
#' @importFrom rlang syms sym
preprocess_comparison_data_grouped <- function(data, x, y, group_by, facet_by = NULL) {
  split_cols <- c(x, y, group_by)
  grouping_vars <- x
  if (!is.null(facet_by)) {
    split_cols <- c(split_cols, facet_by)
    grouping_vars <- c(grouping_vars, facet_by)
  }

  if (length(grouping_vars) > 1) {
    split_key <- interaction(data[grouping_vars], drop = TRUE, sep = " // ")
  } else {
    split_key <- data[[grouping_vars]]
  }

  data_groups <- split(data[, split_cols, drop = FALSE], split_key)
  needs_fix <- check_variance_issues_grouped(data_groups, x, y, group_by)
  if (!needs_fix) return(data)

  warning("Some pairwise comparisons may fail due to insufficient variability. Adjusting data to ensure valid comparisons.", call. = FALSE)

  if (!is.null(facet_by)) {
    facet_key <- interaction(data[facet_by], drop = TRUE, sep = " // ")
    facet_splits <- split(data[, split_cols, drop = FALSE], facet_key)
  } else {
    facet_splits <- list(data[, split_cols, drop = FALSE])
  }

  fixed_data_list <- lapply(seq_along(facet_splits), function(i) {
    facet_data <- facet_splits[[i]]
    result <- fix_facet_variance_grouped(facet_data, x, y, group_by)
    if (!is.null(facet_by)) {
      facet_values <- unique(facet_data[, facet_by, drop = FALSE])
      if (nrow(facet_values) == 1) for (fb in facet_by) result[[fb]] <- facet_values[[fb]]
    }
    return(result)
  })
  do.call(rbind, fixed_data_list)
}

#' Preprocess data for pairwise comparisons without group_by
#'
#' @description
#' Prepares data for simple (non-grouped) pairwise statistical comparisons
#' by detecting and fixing variance issues.
#'
#' @param data Data frame containing the plot data
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @param facet_by Character vector; faceting column name(s), or NULL
#' @return Preprocessed data frame ready for statistical comparisons
#' @keywords internal
preprocess_comparison_data_simple <- function(data, x, y, facet_by = NULL) {
  split_cols <- c(x, y)
  grouping_vars <- x
  if (!is.null(facet_by)) {
    split_cols <- c(split_cols, facet_by)
    grouping_vars <- c(grouping_vars, facet_by)
  }

  if (length(grouping_vars) > 1) {
    split_key <- interaction(data[grouping_vars], drop = TRUE, sep = " // ")
  } else {
    split_key <- data[[grouping_vars]]
  }

  data_groups <- split(data[, split_cols, drop = FALSE], split_key)
  needs_fix <- check_variance_issues_simple(data_groups, y)
  if (!needs_fix) return(data)

  warning("Some pairwise comparisons may fail due to insufficient variability. Adjusting data to ensure valid comparisons.", call. = FALSE)

  if (!is.null(facet_by)) {
    facet_key <- interaction(data[facet_by], drop = TRUE, sep = " // ")
    facet_splits <- split(data[, split_cols, drop = FALSE], facet_key)
  } else {
    facet_splits <- list(data[, split_cols, drop = FALSE])
  }

  fixed_data_list <- lapply(seq_along(facet_splits), function(i) {
    facet_data <- facet_splits[[i]]
    result <- fix_facet_variance_simple(facet_data, x, y)
    if (!is.null(facet_by)) {
      facet_values <- unique(facet_data[, facet_by, drop = FALSE])
      if (nrow(facet_values) == 1) for (fb in facet_by) result[[fb]] <- facet_values[[fb]]
    }
    return(result)
  })
  do.call(rbind, fixed_data_list)
}

# =============================================================================
# Paired Data Validation Helpers
# =============================================================================

#' Format paired validation error message
#' @param problem_groups Data frame containing problematic group combinations
#' @param x Character string; x-axis column name
#' @param paired_by Character string; pairing column name
#' @param group_by Character string; grouping column name (can be NULL)
#' @param n_total_col Character string; name of the count column
#' @param with_group Logical; whether group_by is used
#' @return Character string containing the formatted error message
#' @keywords internal
format_paired_error <- function(problem_groups, x, paired_by, group_by = NULL,
                                n_total_col, with_group = TRUE) {
  cols_to_show <- if (with_group) {
    c(x, paired_by, group_by, ".n", n_total_col)
  } else {
    c(x, paired_by, ".n", n_total_col)
  }

  error_details <- apply(problem_groups[, cols_to_show], 1, function(row) {
    paste(paste(names(row), row, sep = "="), collapse = ", ")
  })

  base_msg <- if (with_group) {
    "When 'paired_by' and 'group_by' are both provided, each combination of 'x' and 'paired_by' must have exactly two observations, one for each group in 'group_by'."
  } else {
    "When 'paired_by' is provided without 'group_by', each combination of 'x' and 'paired_by' must have exactly two observations, one for each value of 'x'."
  }

  paste0(base_msg, " The following combinations do not satisfy this requirement:\n",
    paste0(error_details, collapse = "\n"))
}

#' Validate paired data with groups
#' @param data Data frame to validate
#' @param x Character string; x-axis column name
#' @param paired_by Character string; pairing column name
#' @param group_by Character string; grouping column name
#' @return NULL (invisibly); throws error if validation fails
#' @keywords internal
#' @importFrom dplyr group_by summarise add_count filter mutate
#' @importFrom rlang syms sym
validate_paired_groups <- function(data, x, paired_by, group_by) {
  n_total_col <- paste0(".n_total_", paired_by)

  problem_groups <- data %>%
    dplyr::group_by(!!!syms(c(x, paired_by, group_by))) %>%
    dplyr::summarise(.n = dplyr::n(), .groups = "drop") %>%
    dplyr::add_count(!!!syms(c(x, paired_by)), name = n_total_col) %>%
    dplyr::filter(!!sym(".n") != 1 | !!sym(n_total_col) != 2) %>%
    dplyr::mutate(
      .n = ifelse(!!sym(".n") == 1, !!sym(".n"), paste0(!!sym(".n"), " (expecting 1)")),
      !!sym(n_total_col) := ifelse(!!sym(n_total_col) == 2, !!sym(n_total_col),
        paste0(!!sym(n_total_col), " (expecting 2)"))
    )

  if (nrow(problem_groups) > 0) {
    stop(format_paired_error(problem_groups, x, paired_by, group_by, n_total_col, TRUE),
      call. = FALSE)
  }
}

#' Validate paired data without groups
#' @param data Data frame to validate
#' @param x Character string; x-axis column name
#' @param paired_by Character string; pairing column name
#' @return NULL (invisibly); throws error if validation fails
#' @keywords internal
#' @importFrom dplyr group_by summarise add_count filter mutate n_distinct
#' @importFrom rlang syms sym
validate_paired_simple <- function(data, x, paired_by) {
  if (dplyr::n_distinct(data[[x]], na.rm = TRUE) != 2) {
    stop("Exactly two unique values of 'x' are required when 'paired_by' is provided without 'group_by'.",
      call. = FALSE)
  }

  n_total_col <- paste0(".n_total_", paired_by)

  problem_groups <- data %>%
    dplyr::group_by(!!!syms(c(x, paired_by))) %>%
    dplyr::summarise(.n = dplyr::n(), .groups = "drop") %>%
    dplyr::add_count(!!!syms(paired_by), name = n_total_col) %>%
    dplyr::filter(!!sym(".n") != 1 | !!sym(n_total_col) != 2) %>%
    dplyr::mutate(
      .n = ifelse(!!sym(".n") == 1, !!sym(".n"), paste0(!!sym(".n"), " (expecting 1)")),
      !!sym(n_total_col) := ifelse(!!sym(n_total_col) == 2, !!sym(n_total_col),
        paste0(!!sym(n_total_col), " (expecting 2)"))
    )

  if (nrow(problem_groups) > 0) {
    stop(format_paired_error(problem_groups, x, paired_by, NULL, n_total_col, FALSE),
      call. = FALSE)
  }
}

#' Validate paired data structure
#'
#' @description
#' Main validation function for paired data. Handles NA values in paired_by,
#' validates structure based on whether group_by is provided, and sorts
#' data for proper pairing.
#'
#' @param data Data frame to validate
#' @param x Character string; x-axis column name
#' @param y Character string; y-axis column name
#' @param paired_by Character string; pairing column name
#' @param group_by Character string; grouping column name (can be NULL)
#' @return Validated and sorted data frame
#' @keywords internal
#' @importFrom dplyr arrange
#' @importFrom rlang syms
validate_paired_data <- function(data, x, y, paired_by, group_by = NULL) {
  if (any(is.na(data[[paired_by]]))) {
    warning("'paired_by' contains missing values, removing corresponding rows.", call. = FALSE)
    data <- data[!is.na(data[[paired_by]]), , drop = FALSE]
  }

  if (!is.null(group_by)) {
    validate_paired_groups(data, x, paired_by, group_by)
  } else {
    validate_paired_simple(data, x, paired_by)
  }

  data %>% dplyr::arrange(!!!syms(unique(c(paired_by, x, group_by))))
}

# --- Source: ggforge/R/06-plot-builder.R ---

#' Plot Builder Base Class
#'
#' @description
#' Abstract base system for building consistent plot functions with
#' split/facet/combine support
#'
#' @name plot-builder
#' @keywords internal
NULL

#' Build a plot with standard workflow
#'
#' @description
#' This is the core function that handles the standard plot workflow:
#' 1. Validate parameters
#' 2. Split data if requested
#' 3. Build atomic plots
#' 4. Combine or facet plots
#'
#' @param data Data frame
#' @param atomic_fn Function to create atomic plot
#' @param params List of parameters
#' @param split_by Column to split by
#' @param facet_by Column to facet by
#' @param combine Whether to combine plots
#' @param ... Additional arguments
#' @return Plot or list of plots
#' @keywords internal
build_plot <- function(
    data,
    atomic_fn,
    params,
    split_by = NULL,
    facet_by = NULL,
    combine = TRUE,
    ...) {
  # Validate and prepare split_by column
  # Most wrappers validate split_by before calling build_plot, but some
  # (e.g., LinePlot, DensityPlot) pass it through directly. This ensures
  # the split column exists and is properly factored in all cases.
  if (!is.null(split_by)) {
    split_by <- validate_columns(
      data, split_by,
      force_factor = TRUE,
      allow_multi = TRUE,
      concat_multi = TRUE,
      concat_sep = params$split_by_sep %||% "_"
    )
    data <- apply_column_update(data, split_by)
  }
  # Remove split_by_sep from params to prevent it leaking through ...
  # to ggplot2 geom layers (which would warn "Ignoring unknown parameters")
  params$split_by_sep <- NULL

  # Split data if requested
  if (!is.null(split_by)) {
    data_list <- split_data(data, split_by)
    split_names <- names(data_list)

    # Normalize palette and palcolor for splits
    params$palette <- check_palette(params$palette, split_names)
    params$palcolor <- check_palcolor(params$palcolor, split_names)
    params$legend.position <- check_legend_param(
      params$legend.position, split_names, "legend.position"
    )
    params$legend.direction <- check_legend_param(
      params$legend.direction, split_names, "legend.direction"
    )

    # Build plots for each split
    plots <- lapply(split_names, function(name) {
      split_params <- params
      split_params$palette <- params$palette[[name]]
      split_params$palcolor <- params$palcolor[[name]]
      split_params$legend.position <- params$legend.position[[name]]
      split_params$legend.direction <- params$legend.direction[[name]]
      split_params$title <- format_title(
        params$title, name, paste0(split_by, ": ", name)
      )

      do.call(atomic_fn, c(list(data = data_list[[name]]), split_params))
    })
    names(plots) <- split_names

    # Combine plots
    return(combine_plots(
      plots,
      combine = combine,
      nrow = params$nrow,
      ncol = params$ncol,
      byrow = params$byrow,
      axes = params$axes,
      axis_titles = params$axis_titles,
      guides = params$guides,
      design = params$design
    ))
  } else {
    # Single plot
    params$title <- format_title(params$title, NULL, NULL)
    return(do.call(atomic_fn, c(list(data = data), params)))
  }
}

#' Format title with dynamic content
#'
#' @param title Title (string or function)
#' @param split_name Name of current split
#' @param default_title Default title
#' @return Formatted title
#' @keywords internal
format_title <- function(title, split_name, default_title) {
  if (is.null(title)) {
    return(default_title)
  }

  if (is.function(title)) {
    return(title(default_title))
  }

  return(title)
}

#' Check and normalize legend parameter
#'
#' @param legend_param Legend position or direction
#' @param split_names Names of splits
#' @param param_name Parameter name for error messages
#' @return Normalized legend parameter
#' @keywords internal
#' @importFrom ggplot2 waiver
check_legend_param <- function(legend_param, split_names, param_name) {
  # Handle waiver first before converting to list
  if (inherits(legend_param, "waiver")) {
    legend_param <- list(legend_param)
  } else {
    legend_param <- as.list(legend_param)
  }

  # Check if first element is waiver (in case it was already a list)
  if (length(legend_param) > 0 && inherits(legend_param[[1]], "waiver")) {
    legend_param <- list(legend_param[[1]])
  }

  # Replicate if needed
  if (length(legend_param) == 1 && length(split_names) > 1) {
    legend_param <- rep(legend_param, length(split_names))
  }

  if (length(legend_param) < length(split_names)) {
    stop(
      sprintf(
        "%s length (%d) less than split_by values (%d)",
        param_name, length(legend_param), length(split_names)
      ),
      call. = FALSE
    )
  }

  # Set names
  if (is.null(names(legend_param))) {
    names(legend_param)[seq_along(split_names)] <- split_names
  } else {
    missing <- setdiff(split_names, names(legend_param))
    if (length(missing) > 0) {
      stop(
        sprintf("Missing %s for split_by values: %s", param_name, paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  return(legend_param)
}

#' Template for creating standard plot functions
#'
#' @description
#' This is a template showing the standard structure for plot functions.
#' All plot functions should follow this pattern:
#'
#' 1. A main exported function that:
#'    - Validates common arguments
#'    - Processes theme
#'    - Handles splits via build_plot()
#'
#' 2. An atomic function that:
#'    - Creates the actual ggplot
#'    - Takes simple, validated inputs
#'    - Returns a ggplot object
#'
#' @name plot-template
#' @keywords internal
NULL

# Example structure (not actual code):
#
# ExamplePlot <- function(
#     data, x, y,
#     split_by = NULL, split_by_sep = "_",
#     theme = "theme_ggforge", theme_args = list(),
#     palette = "Paired", palcolor = NULL,
#     ...) {
#
#   # Validate common arguments
#   validate_common_args(...)
#
#   # Process theme
#   theme <- process_theme(theme)
#
#   # Get default parameters
#   params <- get_default_params("example")
#   params <- merge_params(as.list(environment()), params)
#
#   # Build plot using standard workflow
#   build_plot(
#     data = data,
#     atomic_fn = ExamplePlotAtomic,
#     params = params,
#     split_by = split_by,
#     ...
#   )
# }
#
# ExamplePlotAtomic <- function(
#     data, x, y,
#     theme = "theme_ggforge", theme_args = list(),
#     palette = "Paired", palcolor = NULL,
#     ...) {
#
#   # Validate columns
#   x <- validate_columns(data, x)
#   y <- validate_columns(data, y)
#
#   # Create plot
#   p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
#     geom_point() +
#     do.call(theme, theme_args)
#
#   return(p)
# }

# --- Source: ggforge/R/corplot.R ---

#' Calculate Annotation Item
#'
#' @description
#' Helper function to calculate and format individual annotation items
#' for correlation plots. Supports equation, R², p-value, and correlation coefficients.
#'
#' @param item Annotation item type: "eq", "r2", "p", "spearman", "pearson", "kendall", "n"
#' @param dat Data frame with x and y variables
#' @param x X variable name
#' @param y Y variable name
#' @param m Linear model object (lm result)
#' @return Character string of formatted annotation expression (plotmath format)
#' @keywords internal
#' @importFrom stats coef cor
calculate_annotation <- function(item, dat, x, y, m) {
  if (item == "eq") {
    coefs <- stats::coef(m)
    if (is.na(coefs[2])) {
      return(as.character(as.expression(substitute(italic(y) == "NaN"))))
    }
    a <- format(as.numeric(coefs[1]), digits = 2)
    b <- format(as.numeric(abs(coefs[2])), digits = 2)
    if (coefs[2] >= 0) {
      anno_eq <- substitute(italic(y) == a + b %.% italic(x), list(a = a, b = b))
    } else {
      anno_eq <- substitute(italic(y) == a - b %.% italic(x), list(a = a, b = b))
    }
    return(as.character(as.expression(anno_eq)))
  } else if (item == "r2") {
    r2 <- format(summary(m)$r.squared, digits = 2)
    anno_r2 <- substitute(italic(R)^2 ~ "=" ~ r2, list(r2 = r2))
    return(as.character(as.expression(anno_r2)))
  } else if (item == "p") {
    coefs <- summary(m)$coefficients
    if (nrow(coefs) < 2 || all(is.na(coefs[2, ]))) {
      anno_p <- substitute(italic(P) ~ "=" ~ "NA")
    } else {
      pval <- format(coefs[2, 4], digits = 2)
      anno_p <- substitute(italic(P) ~ "=" ~ pvalue, list(pvalue = pval))
    }
    return(as.character(as.expression(anno_p)))
  } else if (item == "spearman") {
    rho <- stats::cor(dat[[x]], dat[[y]], method = "spearman", use = "complete.obs")
    anno_rho <- substitute(
      italic("Spearman's") ~ italic(rho) ~ "=" ~ value,
      list(value = format(rho, digits = 2))
    )
    return(as.character(as.expression(anno_rho)))
  } else if (item == "pearson") {
    r <- stats::cor(dat[[x]], dat[[y]], method = "pearson", use = "complete.obs")
    anno_r <- substitute(
      italic("Pearson's") ~ italic(r) ~ "=" ~ value,
      list(value = format(r, digits = 2))
    )
    return(as.character(as.expression(anno_r)))
  } else if (item == "kendall") {
    tau <- stats::cor(dat[[x]], dat[[y]], method = "kendall", use = "complete.obs")
    anno_tau <- substitute(
      italic("Kendall's") ~ italic(tau) ~ "=" ~ value,
      list(value = format(tau, digits = 2))
    )
    return(as.character(as.expression(anno_tau)))
  } else if (item == "n") {
    n <- sum(stats::complete.cases(dat[[x]], dat[[y]]))
    anno_n <- substitute(italic(N) ~ "=" ~ value, list(value = n))
    return(as.character(as.expression(anno_n)))
  } else {
    stop(
      "Unknown annotation item: ", item,
      ". Expected: eq, r2, p, spearman, pearson, kendall, n",
      call. = FALSE
    )
  }
}

#' Add Point Layers (Normal + Highlight)
#'
#' @description
#' Unified function to add point layers with optional highlighting.
#' Handles both regular geom_point and raster mode (scattermore) for large datasets.
#' Highlighted points are rendered on top with custom styling.
#'
#' @param p ggplot object
#' @param data Data frame with .highlight column
#' @param x X variable name
#' @param y Y variable name
#' @param group_by Grouping variable for coloring
#' @param pt_size Point size for normal points
#' @param pt_shape Point shape (0-25)
#' @param alpha Alpha transparency for normal points
#' @param raster Use raster graphics (scattermore)
#' @param raster_dpi DPI for raster mode
#' @param highlight_color Color for highlighted points
#' @param highlight_size Size for highlighted points
#' @param highlight_alpha Alpha for highlighted points
#' @param highlight_stroke Stroke width for highlighted points
#' @return ggplot object with point layers added
#' @keywords internal
#' @importFrom rlang sym
#' @importFrom ggplot2 geom_point aes
add_point_layers <- function(p, data, x, y, group_by,
                             pt_size, pt_shape, alpha,
                             raster, raster_dpi,
                             highlight_color, highlight_size,
                             highlight_alpha, highlight_stroke) {
  normal_data <- data[!data$.highlight, , drop = FALSE]
  highlight_data <- data[data$.highlight, , drop = FALSE]

  if (isTRUE(raster)) {
    # Raster mode
    if (nrow(normal_data) > 0) {
      p <- p + scattermore::geom_scattermore(
        data = normal_data,
        aes(color = !!sym(group_by)),
        pointsize = ceiling(pt_size),
        alpha = alpha,
        pixels = raster_dpi
      )
    }

    if (nrow(highlight_data) > 0) {
      p <- p +
        scattermore::geom_scattermore(
          data = highlight_data,
          aes(x = !!sym(x), y = !!sym(y)),
          color = highlight_color,
          pointsize = floor(highlight_size) + highlight_stroke,
          alpha = highlight_alpha,
          pixels = raster_dpi,
          inherit.aes = FALSE
        ) +
        scattermore::geom_scattermore(
          data = highlight_data,
          aes(color = !!sym(group_by)),
          pointsize = floor(highlight_size),
          alpha = highlight_alpha,
          pixels = raster_dpi
        )
    }
  } else {
    # Normal mode
    if (nrow(normal_data) > 0) {
      p <- p + geom_point(
        data = normal_data,
        aes(color = !!sym(group_by)),
        size = pt_size,
        shape = pt_shape,
        alpha = alpha
      )
    }

    if (nrow(highlight_data) > 0) {
      p <- p +
        geom_point(
          data = highlight_data,
          aes(x = !!sym(x), y = !!sym(y)),
          color = highlight_color,
          size = highlight_size + highlight_stroke,
          shape = pt_shape,
          alpha = highlight_alpha,
          inherit.aes = FALSE
        ) +
        geom_point(
          data = highlight_data,
          aes(color = !!sym(group_by)),
          size = highlight_size,
          shape = pt_shape,
          alpha = highlight_alpha
        )
    }
  }

  return(p)
}

#' Correlation Plot Atomic
#'
#' @description
#' Creates a single scatter correlation plot for two variables without splitting.
#'
#' @inheritParams parameters
#' @param x Column name for x-axis (numeric)
#' @param y Column name for y-axis (numeric)
#' @param group_by Column name for grouping. Different groups will be colored differently.
#' @param group_by_sep Separator for concatenating multiple columns in group_by
#' @param group_name Name for the group legend
#' @param pt_size Size of the points
#' @param pt_shape Shape of the points (0-25)
#' @param raster Whether to use raster graphics (faster for large datasets)
#' @param raster_dpi DPI for raster graphics as c(width, height)
#' @param highlight Items to highlight. Can be:
#'   - A vector of row indices
#'   - A vector of rownames
#'   - An expression to filter (e.g., "Species == 'setosa'")
#' @param highlight_color Color for highlighted points
#' @param highlight_size Size for highlighted points
#' @param highlight_alpha Alpha for highlighted points
#' @param highlight_stroke Stroke width for highlighted points
#' @param anno_items Annotation items to display. Options: "eq", "r2", "p", "spearman", "pearson", "kendall", "n"
#' @param anno_size Size of annotation text
#' @param anno_fg Foreground color of annotation text
#' @param anno_bg Background color of annotation text
#' @param anno_bg_r Radius of annotation background
#' @param anno_position Position of annotations. Options: "auto", "topleft", "topright", "bottomleft", "bottomright" (or shortcuts: "tl", "tr", "bl", "br")
#' @param add_smooth Whether to add linear regression line
#' @param smooth_color Color of regression line
#' @param smooth_width Width of regression line
#' @param smooth_se Whether to show standard error band
#'
#' @return A ggplot object
#' @keywords internal
#' @importFrom stats cor lm coef complete.cases
#' @importFrom rlang syms sym "%||%"
#' @importFrom dplyr group_by group_modify mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual labs theme waiver
#' @importFrom ggrepel geom_text_repel
CorPlotAtomic <- function(
    data, x, y,
    group_by = NULL,
    group_by_sep = "_",
    group_name = NULL,
    pt_size = 2,
    pt_shape = 16,
    alpha = 1,
    raster = FALSE,
    raster_dpi = c(512, 512),
    highlight = NULL,
    highlight_color = "black",
    highlight_size = 1,
    highlight_alpha = 1,
    highlight_stroke = 0.8,
    anno_items = c("n", "p", "pearson"),
    anno_size = 3.5,
    anno_fg = "black",
    anno_bg = "white",
    anno_bg_r = 0.1,
    anno_position = "auto",
    add_smooth = TRUE,
    smooth_color = "red2",
    smooth_width = 1.5,
    smooth_se = FALSE,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    aspect.ratio = NULL,
    legend.position = waiver(),
    legend.direction = "vertical",
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    ...) {
  # Get ggplot function with gglogger support
  ggplot <- get_ggplot()

  # Validate anno_position
  valid_positions <- c("auto", "topleft", "topright", "bottomleft", "bottomright", "tl", "tr", "bl", "br")
  if (!anno_position %in% valid_positions) {
    stop(
      sprintf("'anno_position' must be one of: %s", paste(valid_positions, collapse = ", ")),
      call. = FALSE
    )
  }

  # Convert shortcuts
  anno_position <- switch(anno_position,
    tl = "topleft",
    tr = "topright",
    bl = "bottomleft",
    br = "bottomright",
    anno_position
  )

  # Normalize raster_dpi
  if (length(raster_dpi) == 1) {
    raster_dpi <- rep(raster_dpi, 2)
  }

  # Validate columns
  x <- validate_columns(data, x)
  y <- validate_columns(data, y)
  group_by <- validate_columns(
    data, group_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = group_by_sep
  )

  # Store whether we have a real group_by for legend positioning
  has_real_group <- !is.null(group_by)

  # Handle NULL group_by by creating a dummy group for consistent behavior
  if (is.null(group_by)) {
    group_by <- ".group"
    data[[group_by]] <- factor("")
  }

  # Normalize legend position based on whether a real group was provided
  legend.position <- normalize_legend_position(
    legend.position,
    has_group = has_real_group
  )

  # Calculate text size scale
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")
  text_size_scale <- base_size / 12

  # Calculate annotations per facet group (needed for faceted plots)
  annodata <- data %>%
    dplyr::group_by(!!!syms(facet_by)) %>%
    dplyr::group_modify(function(dat, g) {
      m <- stats::lm(dat[[y]] ~ dat[[x]])

      # Keep each annotation as a separate row for proper spacing
      anno_list <- lapply(anno_items, function(item) {
        calculate_annotation(item, dat, x, y, m)
      })

      # Calculate correlation for auto positioning (positive -> topleft, negative -> topright)
      corr_val <- stats::cor(dat[[x]], dat[[y]], method = "pearson", use = "complete.obs")

      data.frame(
        anno = unlist(anno_list),
        corr = corr_val,
        stringsAsFactors = FALSE
      )
    })

  # Auto position: positive correlation -> topleft, negative -> topright
  if (anno_position == "auto") {
    first_corr <- annodata$corr[1]
    if (!is.na(first_corr) && first_corr >= 0) {
      anno_position <- "topleft"
    } else {
      anno_position <- "topright"
    }
  }

  # Process highlight: supports TRUE, expression strings, or row indices/names
  data$.highlight <- FALSE
  if (!is.null(highlight)) {
    if (isTRUE(highlight)) {
      # Highlight all points
      data$.highlight <- TRUE
    } else if (length(highlight) == 1 && is.character(highlight)) {
      # Expression filter (e.g., 'Species == "setosa"')
      data <- data %>% dplyr::mutate(.highlight = !!rlang::parse_expr(highlight))
    } else {
      # Row indices or names
      all_inst <- rownames(data) %||% seq_len(nrow(data))
      if (!any(highlight %in% all_inst)) {
        stop("No highlight items found in the data (rownames).", call. = FALSE)
      }
      if (!all(highlight %in% all_inst)) {
        warning("Some highlight items not found in the data (rownames).", call. = FALSE)
      }
      data$.highlight <- all_inst %in% highlight
      rm(all_inst)
    }
  }

  # Initialize plot
  p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y)))

  # Add points
  p <- add_point_layers(
    p, data, x, y, group_by,
    pt_size, pt_shape, alpha,
    raster, raster_dpi,
    highlight_color, highlight_size,
    highlight_alpha, highlight_stroke
  )

  # Add smooth line if requested
  if (add_smooth) {
    p <- p + geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = smooth_se,
      color = smooth_color,
      linewidth = smooth_width,
      alpha = 0.5
    )
  }

  # Add annotations (each as separate row, repel will space them vertically)
  anno_x <- if (grepl("left", anno_position)) -Inf else Inf
  anno_y <- if (grepl("top", anno_position)) Inf else -Inf
  anno_hjust <- if (grepl("left", anno_position)) 0 else 1

  p <- p +
    ggrepel::geom_text_repel(
      data = annodata,
      parse = TRUE,
      hjust = anno_hjust,
      direction = "y",
      aes(label = !!sym("anno")),
      x = anno_x,
      y = anno_y,
      size = text_size_scale * anno_size,
      bg.color = anno_bg,
      bg.r = anno_bg_r,
      color = anno_fg,
      min.segment.length = 0,
      max.overlaps = 100,
      force = 0.5,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "transparent"
    )

  # Add color scale
  p <- p + scale_color_manual(
    name = group_name %||% group_by,
    values = get_palette(
      levels(data[[group_by]]),
      palette = palette,
      palcolor = palcolor
    )
  )

  # Add labels
  p <- p + labs(
    title = title,
    subtitle = subtitle,
    x = xlab %||% x,
    y = ylab %||% y
  )

  # Apply theme
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")
  p <- p +
    do.call(theme, theme_args) +
    ggplot2::theme(aspect.ratio = aspect.ratio)

  # Apply data-driven styling (both x and y are continuous numeric variables)
  p <- apply_style_theme(
    plot = p,
    data = data,
    x_var = x,
    y_var = y,
    flip = FALSE,
    base_size = base_size,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Add faceting if requested
  p <- add_facets(
    p,
    facet_by = facet_by,
    facet_scales = facet_scales,
    nrow = facet_nrow,
    ncol = facet_ncol,
    byrow = facet_byrow
  )

  # Set plot dimensions as attributes
  height <- width <- 4.5
  if (!identical(legend.position, "none")) {
    if (legend.position %in% c("right", "left")) {
      width <- width + 1
    } else if (legend.direction == "horizontal") {
      height <- height + 1
    } else {
      width <- width + 2
    }
  }

  attr(p, "height") <- height
  attr(p, "width") <- width

  return(p)
}

#' Correlation Plot
#'
#' @description
#' Generate scatter correlation plot for two variables with optional linear
#' regression line, annotations, and highlighting.
#'
#' @inheritParams parameters
#' @inheritParams CorPlotAtomic
#' @param ... Additional arguments passed to atomic plotting functions.
#'
#' @return A ggplot object or list/combined plots if split_by is used
# @keywords internal
#' @examples
#' # Basic correlation plot with grouping
#' data(iris)
#' CorPlot(iris, x = "Sepal.Length", y = "Sepal.Width", group_by = "Species")
#'
#' # With custom annotations and positioning
#' CorPlot(iris,
#'   x = "Sepal.Length", y = "Sepal.Width",
#'   group_by = "Species",
#'   anno_items = c("n", "eq", "r2", "pearson"),
#'   anno_position = "bottomright"
#' )
#'
#' # With highlighting specific points
#' CorPlot(iris,
#'   x = "Sepal.Length", y = "Sepal.Width",
#'   group_by = "Species",
#'   highlight = 'Species == "setosa"',
#'   highlight_color = "red",
#'   highlight_size = 3,
#'   highlight_stroke = 1.5
#' )
#'
#' # With faceting by groups
#' CorPlot(iris,
#'   x = "Sepal.Length", y = "Sepal.Width",
#'   facet_by = "Species",
#'   facet_scales = "free",
#'   add_smooth = TRUE,
#'   smooth_color = "blue"
#' )
#'
#' # With splitting and custom palettes
#' CorPlot(iris,
#'   x = "Sepal.Length", y = "Sepal.Width",
#'   split_by = "Species",
#'   palette = c(setosa = "Set1", versicolor = "Dark2", virginica = "Paired"),
#'   combine = TRUE
#' )
CorPlot <- function(
    data, x, y,
    group_by = NULL,
    group_by_sep = "_",
    group_name = NULL,
    split_by = NULL,
    split_by_sep = "_",
    pt_size = 2,
    pt_shape = 16,
    raster = FALSE,
    alpha = 1,
    raster_dpi = c(512, 512),
    highlight = NULL,
    highlight_color = "black",
    highlight_size = 1,
    highlight_alpha = 1,
    highlight_stroke = 0.8,
    anno_items = c("n", "p", "pearson"),
    anno_size = 3.5,
    anno_fg = "black",
    anno_bg = "white",
    anno_bg_r = 0.1,
    anno_position = "auto",
    add_smooth = TRUE,
    smooth_color = "red2",
    smooth_width = 1.5,
    smooth_se = FALSE,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    aspect.ratio = NULL,
    legend.position = waiver(),
    legend.direction = "vertical",
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    seed = 8525,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    axes = NULL,
    axis_titles = NULL,
    guides = NULL,
    design = NULL,
    ...) {
  # Validate common arguments
  validate_common_args(
    seed = seed,
    facet_by = facet_by,
    split_by = split_by,
    theme = theme,
    palette = palette,
    alpha = alpha,
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Process theme
  theme <- process_theme(theme)

  # Validate split_by column
  split_by <- validate_columns(
    data, split_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = split_by_sep
  )

  # Collect all parameters for passing to atomic function
  params <- as.list(environment())
  params$data <- NULL # Remove data from params

  # Build plot using standard workflow
  build_plot(
    data = data,
    atomic_fn = CorPlotAtomic,
    params = params,
    split_by = split_by,
    facet_by = facet_by,
    combine = combine
  )
}

# --- Source: ggforge/R/dotplot.R ---

#' Dot Plot / Scatter Plot
#'
#' @description
#' Creates a dot plot where both X and Y axes can be numeric or categorical.
#' When both are numeric, the plot functions as a scatter plot.
#' Supports sizing dots by a numeric column and filling by another numeric column.
#'
#' @inheritParams parameters
#' @param x A character string specifying the column to use for the x-axis.
#'   Can be numeric or factor/character. When multiple columns are provided,
#'   they will be concatenated with `x_sep`.
#' @param y A character string specifying the column to use for the y-axis.
#'   Can be numeric or factor/character. When multiple columns are provided,
#'   they will be concatenated with `y_sep`.
#' @param x_sep A character string to concatenate multiple columns in x.
#' @param y_sep A character string to concatenate multiple columns in y.
#' @param size_by Which column to use as the size of the dots (numeric column).
#'   If not provided, the size will be the count of instances for each (x, y) pair.
#'   Can also be a single numeric value to specify a fixed size.
#' @param fill_by Which column to use to fill the dots (numeric column).
#'   If not provided, all dots will be filled with the middle color of the palette.
#' @param fill_cutoff A numeric value specifying the cutoff for the fill column.
#'   Values below (or above if `fill_reverse = TRUE`) will be shown in grey.
#' @param fill_reverse Whether to reverse the fill direction. If FALSE (default),
#'   values < cutoff are grey. If TRUE, values > cutoff are grey.
#' @param size_name A character string to name the size legend.
#' @param fill_name A character string to name the fill legend.
#' @param fill_cutoff_name A character string to name the fill cutoff legend.
#' @param flip Whether to flip the x and y axes.
#' @param add_bg Whether to add a background color to the plot.
#' @param bg_palette Palette for the background color.
#' @param bg_palcolor Custom colors for the background.
#' @param bg_alpha Alpha value for the background color.
#' @param bg_direction Direction for background stripes ("vertical" or "horizontal").
#' @param lollipop Whether to create a lollipop plot (requires numeric x and factor y).
#' @param x_text_angle Angle for x-axis text.
#' @param keep_empty Whether to keep empty factor levels.
#' @param ... Additional arguments passed to atomic plotting functions.
#'
#' @return A ggplot object or wrap_plots object or a list of ggplot objects
# @keywords internal
#' @importFrom ggplot2 waiver
#' @examples
#' \donttest{
#' mtcars <- datasets::mtcars
#' mtcars$carb <- factor(mtcars$carb)
#' mtcars$gear <- factor(mtcars$gear)
#' DotPlot(mtcars,
#'   x = "carb", y = "gear", size_by = "wt",
#'   fill_by = "mpg", fill_cutoff = 18
#' )
#' DotPlot(mtcars,
#'   x = "carb", y = "gear", size_by = "wt",
#'   fill_by = "mpg", fill_cutoff = 18, add_bg = TRUE
#' )
#' # Scatter plot (both axes numeric)
#' DotPlot(mtcars,
#'   x = "qsec", y = "drat", size_by = "wt",
#'   fill_by = "mpg", fill_cutoff = 18
#' )
#' }
DotPlot <- function(
    data, x, y,
    x_sep = "_",
    y_sep = "_",
    flip = FALSE,
    split_by = NULL,
    split_by_sep = "_",
    size_by = NULL,
    fill_by = NULL,
    fill_cutoff = NULL,
    fill_reverse = FALSE,
    size_name = NULL,
    fill_name = NULL,
    fill_cutoff_name = NULL,
    add_bg = FALSE,
    bg_palette = "stripe",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    bg_direction = c("vertical", "horizontal", "v", "h"),
    lollipop = FALSE,
    theme = "theme_ggforge_grid",
    theme_args = list(),
    palette = "Spectral",
    palcolor = NULL,
    alpha = 1,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    x_text_angle = 0,
    seed = 8525,
    aspect.ratio = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    keep_empty = FALSE,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    axes = NULL,
    axis_titles = NULL,
    guides = NULL,
    design = NULL,
    ...) {
  # Validate common arguments
  validate_common_args(
    seed = seed,
    facet_by = facet_by,
    split_by = split_by,
    theme = theme,
    palette = palette,
    alpha = alpha,
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Validate bg_direction
  bg_direction <- match.arg(bg_direction)
  if (bg_direction %in% c("h", "horizontal")) {
    bg_direction <- "horizontal"
  } else {
    bg_direction <- "vertical"
  }

  # Process theme
  theme <- process_theme(theme)

  # Validate split_by column
  split_by <- validate_columns(
    data, split_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = split_by_sep
  )

  # Collect all parameters for passing to atomic function
  params <- as.list(environment())
  params$data <- NULL # Remove data from params

  # Build plot using standard workflow
  build_plot(
    data = data,
    atomic_fn = DotPlotAtomic,
    params = params,
    split_by = split_by,
    facet_by = facet_by,
    combine = combine
  )
}

#' Dot Plot Atomic
#'
#' @description
#' Creates a single dot plot without splitting
#'
#' @inheritParams DotPlot
#' @keywords internal
#' @importFrom ggplot2 aes geom_point scale_size_area scale_fill_gradientn
#' @importFrom ggplot2 scale_color_manual scale_x_discrete scale_y_discrete
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous coord_flip
#' @importFrom ggplot2 guide_colorbar guide_legend guide_none guides labs theme
#' @importFrom ggplot2 element_line element_text geom_segment scale_color_gradientn
#' @importFrom ggnewscale new_scale_color
#' @importFrom rlang sym syms "%||%"
DotPlotAtomic <- function(
    data, x, y,
    x_sep = "_",
    y_sep = "_",
    flip = FALSE,
    lollipop = FALSE,
    size_by = NULL,
    fill_by = NULL,
    fill_cutoff = NULL,
    fill_reverse = FALSE,
    size_name = NULL,
    fill_name = NULL,
    fill_cutoff_name = NULL,
    theme = "theme_ggforge_grid",
    theme_args = list(),
    palette = "Spectral",
    palcolor = NULL,
    alpha = 1,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    x_text_angle = 0,
    aspect.ratio = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    add_bg = FALSE,
    bg_palette = "stripe",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    bg_direction = c("vertical", "horizontal", "v", "h"),
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    keep_empty = FALSE,
    ...) {
  # Determine if x and y are numeric BEFORE validation
  x_is_numeric <- length(x) == 1 && !is.character(data[[x]]) && !is.factor(data[[x]])
  y_is_numeric <- length(y) == 1 && !is.character(data[[y]]) && !is.factor(data[[y]])

  # Validate and process columns
  if (!x_is_numeric) {
    x <- validate_columns(
      data, x,
      force_factor = TRUE,
      allow_multi = TRUE,
      concat_multi = TRUE,
      concat_sep = x_sep
    )
  }

  if (!y_is_numeric) {
    y <- validate_columns(
      data, y,
      force_factor = TRUE,
      allow_multi = TRUE,
      concat_multi = TRUE,
      concat_sep = y_sep
    )
  }

  # Validate fill_cutoff usage
  if (!is.null(fill_cutoff) && is.null(fill_by)) {
    stop("'fill_by' must be provided when 'fill_cutoff' is specified.", call. = FALSE)
  }

  # Validate facet_by
  facet_by <- validate_columns(
    data, facet_by,
    force_factor = TRUE,
    allow_multi = TRUE
  )

  # Handle size_by
  if (!is.numeric(size_by)) {
    size_by <- validate_columns(data, size_by)
  }

  if (is.null(size_by)) {
    if (is.null(fill_by)) {
      data <- data |>
        dplyr::group_by(!!!syms(unique(c(x, y, facet_by)))) |>
        dplyr::summarise(.size = dplyr::n(), .groups = "drop")
    } else {
      # Aggregate data and check if there are duplicates
      data <- data |>
        dplyr::group_by(!!!syms(unique(c(x, y, facet_by)))) |>
        dplyr::summarise(!!sym(fill_by) := dplyr::first(!!sym(fill_by)), .size = dplyr::n(), .groups = "drop")

      # Only warn if there are actually duplicates (count > 1)
      if (any(data$.size > 1)) {
        warning("Using the first value of fill_by as size_by is calculated from count of instances.", immediate. = TRUE)
      }
    }
    size_by <- ".size"
    # Set default size_name if not provided
    if (is.null(size_name)) {
      size_name <- "Count"
    }
  }

  # Handle fill_by
  fill_by <- validate_columns(data, fill_by)

  if (!is.null(fill_by) && !is.null(fill_cutoff)) {
    # Add a column to indicate the fill cutoff
    if (isFALSE(fill_reverse)) {
      fill_cutoff_label <- paste0(fill_by, " < ", fill_cutoff)
      data[[fill_by]][data[[fill_by]] < fill_cutoff] <- NA
    } else {
      fill_cutoff_label <- paste0(fill_by, " > ", fill_cutoff)
      data[[fill_by]][data[[fill_by]] > fill_cutoff] <- NA
    }
  }

  if (is.null(fill_by)) {
    data$.fill_by <- 1
    fill_by <- ".fill_by"
    fill_legend <- FALSE
  } else {
    fill_legend <- TRUE
  }

  # Get ggplot (support gglogger)
  ggplot <- get_ggplot()

  # Calculate text justification
  just <- calc_justification(x_text_angle)

  # Build plot
  p <- ggplot(data, ggplot2::aes(x = !!sym(x), y = !!sym(y)))

  # Add background layer if requested
  if (isTRUE(add_bg)) {
    if (bg_direction %in% c("vertical", "v")) {
      if (x_is_numeric) {
        stop("Vertical 'bg_direction' is not supported when 'x' is numeric.", call. = FALSE)
      }
      p <- p + create_bg_layer(data, x, bg_palette, bg_palcolor, bg_alpha, keep_empty, facet_by, "vertical")
    } else {
      if (y_is_numeric) {
        stop("Horizontal 'bg_direction' is not supported when 'y' is numeric.", call. = FALSE)
      }
      p <- p + create_bg_layer(data, y, bg_palette, bg_palcolor, bg_alpha, keep_empty, facet_by, "horizontal")
    }
  }

  # Add scale for x and y
  if (!x_is_numeric) {
    p <- p + ggplot2::scale_x_discrete(drop = !keep_empty)
  }
  if (!y_is_numeric) {
    p <- p + ggplot2::scale_y_discrete(drop = !keep_empty)
  }

  # Add lollipop if requested
  if (isTRUE(lollipop)) {
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = !!sym(x), yend = !!sym(y)), color = "black", linewidth = 2) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = !!sym(x), yend = !!sym(y), color = !!sym(fill_by)), linewidth = 1) +
      ggplot2::scale_x_continuous(expand = c(0, 0, 0.05, 0)) +
      ggplot2::scale_color_gradientn(
        n.breaks = 5,
        colors = get_palette(
          1:100,
          palette = palette,
          palcolor = palcolor,
          type = "continuous",
          reverse = fill_reverse
        ),
        na.value = "grey80",
        guide = "none"
      ) +
      ggnewscale::new_scale_color()
  }

  # Add points with size and fill
  if (is.numeric(size_by)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = !!sym(fill_by), color = ""),
      size = size_by,
      shape = 21,
      alpha = alpha
    )
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(size = !!sym(size_by), fill = !!sym(fill_by), color = ""),
        shape = 21,
        alpha = alpha
      ) +
      ggplot2::scale_size_area(max_size = 6, n.breaks = 4) +
      ggplot2::guides(size = ggplot2::guide_legend(
        title = size_name %||% size_by,
        override.aes = list(fill = "transparent", shape = 21, colour = "black"),
        order = 1
      ))
  }

  # Add fill scale
  p <- p +
    ggplot2::scale_fill_gradientn(
      n.breaks = 5,
      colors = get_palette(
        1:100,
        palette = palette,
        palcolor = palcolor,
        type = "continuous",
        reverse = fill_reverse
      ),
      na.value = "grey80",
      guide = if (isTRUE(fill_legend)) {
        ggplot2::guide_colorbar(
          title = fill_name %||% fill_by,
          frame.colour = "black",
          ticks.colour = "black",
          frame.linewidth = 0.3,
          ticks.linewidth = 0.3,
          title.hjust = 0,
          order = 2
        )
      } else {
        ggplot2::guide_none()
      }
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab %||% x,
      y = ylab %||% y
    )

  # Apply theme
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")
  p <- p +
    do.call(theme, theme_args) +
    ggplot2::theme(aspect.ratio = aspect.ratio)

  # Apply data-driven styling (with x-axis angle override)
  p <- apply_style_theme(
    plot = p,
    data = data,
    x_var = if (flip) y else x,
    y_var = if (flip) x else y,
    flip = FALSE, # Already handled above
    base_size = base_size,
    legend.position = legend.position,
    legend.direction = legend.direction,
    axis.text.x = ggplot2::element_text(
      angle = x_text_angle,
      hjust = just$h,
      vjust = just$v
    )
  )

  # Add color scale and fill cutoff guide
  p <- p + ggplot2::scale_color_manual(values = "black", na.value = "black", guide = "none")

  if (!is.null(fill_by) && !is.null(fill_cutoff) && anyNA(data[[fill_by]])) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(
      title = fill_cutoff_name %||% fill_cutoff_label,
      override.aes = list(colour = "black", fill = "grey80", size = 3),
      order = 3
    ))
  }

  # Flip coordinates if requested
  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }

  # Add faceting if requested
  if (!is.null(facet_by)) {
    p <- add_facets(
      p, facet_by, facet_scales,
      facet_nrow, facet_ncol, facet_byrow
    )
  }

  return(p)
}

#' Lollipop Plot
#'
#' @description
#' Creates a lollipop plot with a numeric x-axis and categorical y-axis.
#' This is a convenience wrapper around `DotPlot` with `lollipop = TRUE`.
#'
#' @inheritParams DotPlot
#' @param x A character string specifying the column to use for the x-axis
#'   (numeric column expected).
#' @param y A character string specifying the column to use for the y-axis
#'   (factor/character column expected).
#' @param ... Additional arguments passed to atomic plotting functions.
#'
#' @return A ggplot object or wrap_plots object or a list of ggplot objects
# @keywords internal
#' @examples
#' \donttest{
#' mtcars <- datasets::mtcars
#' LollipopPlot(mtcars,
#'   x = "qsec", y = "gear", size_by = "wt",
#'   fill_by = "mpg"
#' )
#' LollipopPlot(mtcars,
#'   x = "qsec", y = "gear", size_by = "wt",
#'   fill_by = "mpg", fill_cutoff = 18
#' )
#' }
LollipopPlot <- function(
    data, x, y,
    y_sep = "_",
    flip = FALSE,
    split_by = NULL,
    split_by_sep = "_",
    size_by = NULL,
    fill_by = NULL,
    fill_cutoff = NULL,
    fill_reverse = FALSE,
    size_name = NULL,
    fill_name = NULL,
    fill_cutoff_name = NULL,
    theme = "theme_ggforge_grid",
    theme_args = list(),
    palette = "Spectral",
    palcolor = NULL,
    alpha = 1,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    x_text_angle = 0,
    seed = 8525,
    aspect.ratio = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    keep_empty = FALSE,
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    axes = NULL,
    axis_titles = NULL,
    guides = NULL,
    design = NULL,
    ...) {
  DotPlot(
    data = data,
    x = x,
    y = y,
    x_sep = NULL,
    y_sep = y_sep,
    flip = flip,
    lollipop = TRUE,
    split_by = split_by,
    split_by_sep = split_by_sep,
    size_by = size_by,
    fill_by = fill_by,
    fill_cutoff = fill_cutoff,
    fill_reverse = fill_reverse,
    size_name = size_name,
    fill_name = fill_name,
    fill_cutoff_name = fill_cutoff_name,
    theme = theme,
    theme_args = theme_args,
    palette = palette,
    palcolor = palcolor,
    alpha = alpha,
    facet_by = facet_by,
    facet_scales = facet_scales,
    facet_ncol = facet_ncol,
    facet_nrow = facet_nrow,
    facet_byrow = facet_byrow,
    x_text_angle = x_text_angle,
    seed = seed,
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    keep_empty = keep_empty,
    combine = combine,
    nrow = nrow,
    ncol = ncol,
    byrow = byrow,
    axes = axes,
    axis_titles = axis_titles,
    guides = guides,
    design = design,
    ...
  )
}

# --- Source: ggforge/R/heatmap.R ---

#' Join the meta data to the main data frame for heatmap
#'
#' @description
#' Helper function to join metadata (rows_data or columns_data) to the main heatmap data.
#' Preserves factor levels after joining.
#'
#' @param data A data frame containing the main data for the heatmap.
#' @param meta_data A data frame containing the meta data to be joined.
#' @param by A character string specifying the column name in `meta_data` to join on.
#' Either `rows_by` or `columns_by` should be specified in `data`.
#' @param cr_split_by A character string specifying the column name in `data` to join on.
#' Either `rows_split_by` or `columns_split_by` should be specified in `data`.
#' @param split_by A character string specifying the column name in `data` to join on.
#' Used to split the data into multiple heatmaps.
#' @param which A character string specifying whether to join on rows or columns.
#' Can be either `"row"` or `"column"`.
#' @importFrom dplyr left_join
#' @return A data frame with the meta data joined to the main data.
#' @keywords internal
join_heatmap_meta <- function(data, meta_data, by, cr_split_by, split_by, which) {
  if (!by %in% colnames(meta_data)) {
    stop(
      sprintf("[Heatmap] '%ss_by' (%s) must be a column in '%ss_data'", which, by, which),
      call. = FALSE
    )
  }

  join_by <- by
  if (!is.null(cr_split_by) && cr_split_by %in% colnames(data) &&
    cr_split_by %in% colnames(meta_data)) {
    # if split_by is in both data and meta_data, we join by both
    join_by <- c(by, cr_split_by)
  }
  if (!is.null(split_by) && split_by %in% colnames(data) &&
    split_by %in% colnames(meta_data)) {
    # if split_by is in both data and meta_data, we join by both
    join_by <- c(join_by, split_by)
  }

  # Join metadata while preserving factor levels
  out <- dplyr::left_join(data, meta_data, by = join_by, suffix = c("", paste0(".", which)))

  # Restore factor levels lost during join
  out[[by]] <- factor(out[[by]], levels = levels(data[[by]]))
  if (!is.null(cr_split_by) && cr_split_by %in% join_by) {
    if (is.factor(data[[cr_split_by]])) {
      out[[cr_split_by]] <- factor(out[[cr_split_by]], levels = levels(data[[cr_split_by]]))
    } else if (is.factor(meta_data[[cr_split_by]])) {
      # if the split_by is a factor in meta_data, we need to restore the levels
      out[[cr_split_by]] <- factor(out[[cr_split_by]], levels = levels(meta_data[[cr_split_by]]))
    } else {
      # if the split_by is not a factor, we convert it to a factor
      out[[cr_split_by]] <- factor(out[[cr_split_by]], levels = unique(out[[cr_split_by]]))
    }
  }
  if (!is.null(split_by) && split_by %in% join_by) {
    if (is.factor(data[[split_by]])) {
      out[[split_by]] <- factor(out[[split_by]], levels = levels(data[[split_by]]))
    } else if (is.factor(meta_data[[split_by]])) {
      # if the split_by is a factor in meta_data, we need to restore the levels
      out[[split_by]] <- factor(out[[split_by]], levels = levels(meta_data[[split_by]]))
    } else {
      # if the split_by is not a factor, we convert it to a factor
      out[[split_by]] <- factor(out[[split_by]], levels = unique(out[[split_by]]))
    }
  }

  out
}

#' Process/normalize data passed to [Heatmap()]
#'
#' This function is used to process the data passed to [Heatmap()].
#' @param data A data frame or matrix containing the data to be plotted.
#' Based on the `in_form`, the data can have the following formats:
#' * `matrix`: A matrix with rows and columns directly representing the heatmap.
#' * `long`: A data frame in long format with columns for values, rows, and columns.
#' * `wide-rows`: A data frame in wide format with columns for heatmap rows and values,
#'    and a single column for heatmap columns.
#' * `wide-columns`: A data frame in wide format with columns for heatmap columns and values,
#'    and a single column for heatmap rows.
#' * `auto`: Automatically inferred from the data format.
#'    When `data` is a matrix, `in_form` is set to `"matrix"`. When `columns_by` has more than one column,
#'    `in_form` is set to `"wide-columns"`. When `rows_by` has more than one column,
#'    `in_form` is set to `"wide-rows"`. Otherwise, it is set to `"long"`.
#' @param in_form The format of the data. Can be one of `"matrix"`, `"long"`, `"wide-rows"`, `"wide-columns"`, or `"auto"`.
#' Defaults to `"auto"`.
#' @param values_by A character of column name in `data` that contains the values to be plotted.
#' This is required when `in_form` is `"long"`. For other formats, the values are pivoted into a column named by `values_by`.
#' @param name A character string to name the heatmap (will be used to rename `values_by`).
#' @param split_by A character of column name in `data` that contains the split information to split into multiple heatmaps.
#' This is used to create a list of heatmaps, one for each level of the split.
#' Defaults to `NULL`, meaning no split.
#' @param split_by_sep A character string to concat multiple columns in `split_by`.
#' @param rows_by A vector of column names in `data` that contains the row information.
#' This is used to create the rows of the heatmap.
#' When `in_form` is `"long"` or `"wide-columns"`, this is requied, and multiple columns can be specified,
#' which will be concatenated by `rows_by_sep` into a single column.
#' @param rows_by_sep A character string to concat multiple columns in `rows_by`.
#' @param rows_name A character string to rename the column created by `rows_by`, which will be reflected in the name of the annotation or legend.
#' @param rows_split_by A character of column name in `data` that contains the split information for rows.
#' @param rows_split_by_sep A character string to concat multiple columns in `rows_split_by`.
#' @param rows_split_name A character string to rename the column created by `rows_split_by`, which will be reflected in the name of the annotation or legend.
#' @param columns_by A vector of column names in `data` that contains the column information.
#' This is used to create the columns of the heatmap.
#' When `in_form` is `"long"` or `"wide-rows"`, this is required, and multiple columns can be specified,
#' which will be concatenated by `columns_by_sep` into a single column.
#' @param columns_by_sep A character string to concat multiple columns in `columns_by`.
#' @param columns_name A character string to rename the column created by `columns_by`, which will be reflected in the name of the annotation or legend.
#' @param columns_split_by A character of column name in `data` that contains the split information for columns.
#' @param columns_split_by_sep A character string to concat multiple columns in `columns_split_by`.
#' @param columns_split_name A character string to rename the column created by `columns_split_by`, which will be reflected in the name of the annotation or legend.
#' @param pie_group_by A character of column name in `data` that contains the group information for pie charts.
#' This is used to create pie charts in the heatmap when `cell_type` is `"pie"`.
#' @param pie_group_by_sep A character string to concat multiple columns in `pie_group_by`.
#' @param pie_name A character string to rename the column created by `pie_group_by`, which will be reflected in the name of the annotation or legend.
#' @param rows_data A data frame containing additional data for rows, which can be used to add annotations to the heatmap.
#' It will be joined to the main data by `rows_by` and `split_by` if `split_by` exists in `rows_data`.
#' This is useful for adding additional information to the rows of the heatmap.
#' @param columns_data A data frame containing additional data for columns, which can be used to add annotations to the heatmap.
#' It will be joined to the main data by `columns_by` and `split_by` if `split_by` exists in `columns_data`.
#' This is useful for adding additional information to the columns of the heatmap.
#' @return A list containing the processed data and metadata:
#' * `data`: A list of data frames, one for each level of `split_by`. If no `split_by` is provided, the name will be `"..."`.
#'    Each data frame is in the long format.
#' * `values_by`: The name of the column containing the values to be plotted.
#' * `rows_by`: The name of the column containing the row information.
#' * `rows_split_by`: The name of the column containing the row split information.
#' * `columns_by`: The name of the column containing the column information.
#' * `columns_split_by`: The name of the column containing the column split information.
#' * `pie_group_by`: The name of the column containing the pie group information.
#' @importFrom rlang sym %||%
#' @keywords internal
process_heatmap_data <- function(
    data, in_form, values_by, name,
    split_by, split_by_sep,
    rows_by, rows_by_sep, rows_name,
    rows_split_by, rows_split_by_sep, rows_split_name,
    columns_by, columns_by_sep, columns_name,
    columns_split_by, columns_split_by_sep, columns_split_name,
    pie_group_by, pie_group_by_sep, pie_name,
    rows_data, columns_data) {
  if (identical(rows_by, columns_by) && !is.null(rows_by)) {
    stop("[Heatmap] 'rows_by' and 'columns_by' cannot be the same", call. = FALSE)
  }
  # Auto-detect data format if not specified
  # Priority: matrix > wide-rows > wide-columns > long (default)
  if (in_form == "auto") {
    if (is.matrix(data)) {
      in_form <- "matrix"
    } else if (length(rows_by) > 1) {
      in_form <- "wide-rows"
    } else if (length(columns_by) > 1) {
      in_form <- "wide-columns"
    } else {
      in_form <- "long"
    }
  }


  if (in_form == "matrix") {
    if (!is.null(split_by)) {
      stop("[Heatmap] 'split_by' is not supported when 'in_form = \"matrix\"'", call. = FALSE)
    }
    if (!is.null(rows_by)) {
      stop("[Heatmap] 'rows_by' is not supported when 'in_form = \"matrix\"'", call. = FALSE)
    }
    if (!is.null(columns_by)) {
      stop("[Heatmap] 'columns_by' is not supported when 'in_form = \"matrix\"'", call. = FALSE)
    }
    if (!is.null(pie_group_by)) {
      stop("[Heatmap] 'pie_group_by' is not supported when 'in_form = \"matrix\"'", call. = FALSE)
    }

    rows_name <- rows_name %||% "rows"
    data <- as.data.frame(data)
    columns_by <- colnames(data)
    data[rows_name] <- rownames(data)
    rows_by <- rows_name

    in_form <- "wide-columns"
  }

  if (identical(rows_name %||% rows_by, columns_name %||% columns_by)) {
    if (!is.null(columns_name)) {
      # consider flip and names_side?
      columns_name <- paste0(columns_name, " ")
    } else {
      rows_name <- paste0(" ", rows_name)
    }
  }

  # pie_group_by should be always in the main data
  pie_group_by <- validate_columns(
    data, pie_group_by,
    force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = pie_group_by_sep
  )

  split_by <- validate_columns(
    data, split_by,
    force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = split_by_sep
  )

  if (in_form == "long") {
    # values_by
    values_by <- validate_columns(data, values_by)
    if (is.null(values_by)) {
      stop("[Heatmap] 'values_by' must be specified when 'in_form = \"long\"'", call. = FALSE)
    }

    # rows_by/rows_split_by
    rows_by <- validate_columns(
      data, rows_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = rows_by_sep
    )
    if (is.null(rows_by)) {
      stop("[Heatmap] 'rows_by' must be specified when 'in_form = \"long\"'", call. = FALSE)
    }
    data[[rows_by]] <- droplevels(data[[rows_by]])

    # columns_by/columns_split_by
    columns_by <- validate_columns(
      data, columns_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = columns_by_sep
    )
    if (is.null(columns_by)) {
      stop("[Heatmap] 'columns_by' must be specified when 'in_form = \"long\"'", call. = FALSE)
    }
    data[[columns_by]] <- droplevels(data[[columns_by]])

    # join rows_data/columns_data
    if (!is.null(rows_data)) {
      data <- join_heatmap_meta(
        data, rows_data,
        by = rows_by, cr_split_by = rows_split_by,
        split_by = split_by, which = "row"
      )
    }
    rows_split_by <- validate_columns(
      data, rows_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = rows_split_by_sep
    )
    if (!is.null(rows_split_by)) {
      data[[rows_split_by]] <- droplevels(data[[rows_split_by]])
    }

    if (!is.null(columns_data)) {
      data <- join_heatmap_meta(
        data, columns_data,
        by = columns_by, cr_split_by = columns_split_by,
        split_by = split_by, which = "column"
      )
    }
    columns_split_by <- validate_columns(
      data, columns_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = columns_split_by_sep
    )
    if (!is.null(columns_split_by)) {
      data[[columns_split_by]] <- droplevels(data[[columns_split_by]])
    }

    # rename
    if (!is.null(rows_name)) {
      data <- dplyr::rename(data, !!sym(rows_name) := rows_by)
      rows_by <- rows_name
    }
    if (!is.null(columns_name)) {
      data <- dplyr::rename(data, !!sym(columns_name) := columns_by)
      columns_by <- columns_name
    }
  } else if (in_form == "wide-rows") {
    # columns_split_by columns_by pie_group_by rows1 rows2 ...
    # csb1             cb1        pgb1         0.1   0.2   ...
    # csb2             cb2        pgb2         0.3   0.4   ...
    #                    ... ...
    # rows_by
    rows_by <- rows_by %||% setdiff(colnames(data), c(columns_by, columns_split_by, pie_group_by))
    rows_name <- rows_name %||% ifelse("rows" %in% colnames(data), "rows.1", "rows")
    values_by <- values_by %||% ifelse("value" %in% colnames(data), "value.1", "value")
    data <- tidyr::pivot_longer(data, cols = rows_by, names_to = rows_name, values_to = values_by)
    data[[rows_name]] <- factor(data[[rows_name]], levels = unique(rows_by))
    data <- data[order(data[[rows_name]]), , drop = FALSE]
    rows_by <- rows_name

    # columns_by/columns_split_by
    columns_by <- validate_columns(
      data, columns_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = columns_by_sep
    )
    if (is.null(columns_by)) {
      stop("[Heatmap] 'columns_by' must be specified when 'in_form = \"wide-rows\"'", call. = FALSE)
    }
    data[[columns_by]] <- droplevels(data[[columns_by]])

    # rows_data/columns_data
    if (!is.null(rows_data)) {
      data <- join_heatmap_meta(
        data, rows_data,
        by = rows_by, cr_split_by = rows_split_by,
        split_by = split_by, which = "row"
      )
    }
    row_split_by <- validate_columns(
      data, rows_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = rows_split_by_sep
    )
    if (!is.null(row_split_by)) {
      data[[row_split_by]] <- droplevels(data[[row_split_by]])
    }
    if (!is.null(columns_data)) {
      data <- join_heatmap_meta(
        data, columns_data,
        by = columns_by, cr_split_by = columns_split_by,
        split_by = split_by, which = "column"
      )
    }
    columns_split_by <- validate_columns(
      data, columns_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = columns_split_by_sep
    )
    if (!is.null(columns_split_by)) {
      data[[columns_split_by]] <- droplevels(data[[columns_split_by]])
    }

    if (!is.null(columns_name)) {
      data <- dplyr::rename(data, !!sym(columns_name) := columns_by)
      columns_by <- columns_name
    }
  } else { # wide-columns
    # rows_split_by rows_by pie_group_by columns1 columns2 ...
    # rsb1          cb1        pgb1         0.1   0.2   ...
    # rsb2          cb2        pgb2         0.3   0.4   ...
    #                    ... ...
    # rows_by/rows_split_by
    rows_by <- validate_columns(
      data, rows_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = rows_by_sep
    )
    if (is.null(rows_by)) {
      stop("[Heatmap] 'rows_by' must be specified when 'in_form = \"wide-columns\"'", call. = FALSE)
    }

    # columns_by
    columns_by <- columns_by %||% setdiff(colnames(data), c(rows_by, rows_split_by, pie_group_by))
    columns_name <- columns_name %||% ifelse("columns" %in% colnames(data), "columns.1", "columns")
    values_by <- values_by %||% ifelse("value" %in% colnames(data), "value.1", "value")
    data <- tidyr::pivot_longer(data, cols = columns_by, names_to = columns_name, values_to = values_by)
    data[[columns_name]] <- factor(data[[columns_name]], levels = columns_by)
    data <- data[order(data[[columns_name]]), , drop = FALSE]
    columns_by <- columns_name

    # rows_data/columns_data
    if (!is.null(rows_data)) {
      data <- join_heatmap_meta(
        data, rows_data,
        by = rows_by, cr_split_by = rows_split_by,
        split_by = split_by, which = "row"
      )
    }
    rows_split_by <- validate_columns(
      data, rows_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = rows_split_by_sep
    )
    if (!is.null(rows_split_by)) {
      data[[rows_split_by]] <- droplevels(data[[rows_split_by]])
    }
    if (!is.null(columns_data)) {
      data <- join_heatmap_meta(
        data, columns_data,
        by = columns_by, cr_split_by = columns_split_by,
        split_by = split_by, which = "column"
      )
    }
    columns_split_by <- validate_columns(
      data, columns_split_by,
      force_factor = TRUE, allow_multi = TRUE, concat_multi = TRUE, concat_sep = columns_split_by_sep
    )
    if (!is.null(columns_split_by)) {
      data[[columns_split_by]] <- droplevels(data[[columns_split_by]])
    }

    if (!is.null(rows_name)) {
      data <- dplyr::rename(data, !!sym(rows_name) := rows_by)
      rows_by <- rows_name
    }
  }

  if (!is.null(rows_split_name) && !is.null(rows_split_by)) {
    data <- dplyr::rename(data, !!sym(rows_split_name) := rows_split_by)
    rows_split_by <- rows_split_name
  }
  if (!is.null(columns_split_name) && !is.null(columns_split_by)) {
    data <- dplyr::rename(data, !!sym(columns_split_name) := columns_split_by)
    columns_split_by <- columns_split_name
  }
  if (!is.null(pie_name) && !is.null(pie_group_by)) {
    data <- dplyr::rename(data, !!sym(pie_name) := pie_group_by)
    pie_group_by <- pie_name
  }
  if (!is.null(name)) {
    data <- dplyr::rename(data, !!sym(name) := !!sym(values_by))
    values_by <- name
  }

  list(
    data = if (is.null(split_by)) {
      stats::setNames(list(data), "...")
    } else {
      split(select(data, -!!sym(split_by)), data[[split_by]])
    },
    values_by = values_by,
    rows_by = rows_by,
    rows_split_by = rows_split_by,
    columns_by = columns_by,
    columns_split_by = columns_split_by,
    pie_group_by = pie_group_by
  )
}

#' Get the grid.draw-able ggplot grob
#' The output from ggplotGrob can not be directly used in grid.draw, the position can not be set.
#' This function extracts the gTree from the ggplot grob.
#' @param p A ggplot object
#' @param void If TRUE, the theme_void will be added to the ggplot object
#' @param nolegend If TRUE, the legend will be removed from the ggplot object
#' @return A gTree object
#' @importFrom ggplot2 ggplotGrob theme_void theme
#' @keywords internal
gggrob <- function(p, void = TRUE, nolegend = TRUE) {
  if (isTRUE(void)) {
    p <- p + theme_void()
  }
  if (isTRUE(nolegend)) {
    p <- p + theme(legend.position = "none")
  }
  for (g in ggplotGrob(p)$grobs) {
    if (inherits(g, "gTree") && !inherits(g, "zeroGrob") && !inherits(g, "absoluteGrob")) {
      return(g)
    }
  }
  stop("No gTree found in the ggplot grob.", call. = FALSE)
}

#' Heatmap annotation function for categorical data
#' @param x A data frame
#' @param split_by A character string of the column name to split the data
#' @param group_by A character string of the column name to group the data
#' @param column A character string of the column name to plot
#' @param title A character string to name the legend
#' @param which A character string specifying the direction of the annotation. Default is "row".
#'  Other options are "column".
#' @param palette A character string specifying the palette of the annotation
#' @param palcolor A character vector of colors to override the palette
#' @param border A logical value indicating whether to draw the border of the annotation
#' @param legend.direction A character string specifying the direction of the legend. Default is "vertical".
#'  Other options are "horizontal".
#' @param show_legend A logical value indicating whether to show the legend
#' @param .plotting A function to create the plot for each split and each group
#' @param ... Other arguments passed to `ComplexHeatmap::AnnotationFunction`
#' @keywords internal
#' @importFrom grid grid.draw grid.lines viewport gpar
#' @importFrom tidyr unite
.anno_ggcat <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                        palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, .plotting, ...) {
  column <- validate_columns(x, column, force_factor = TRUE)
  clevels <- levels(x[[column]])

  if (!is.null(split_by)) {
    x <- x %>% unite("..split", split_by, group_by, remove = FALSE)
    plots <- .plotting(data = x, column = column, group_by = "..split", palette = palette, palcolor = palcolor)
    plots <- lapply(plots, gggrob)
  } else {
    plots <- .plotting(data = x, column = column, group_by = group_by, palette = palette, palcolor = palcolor)
    plots <- lapply(plots, gggrob)
  }
  # add legend
  if (isTRUE(show_legend)) {
    lgd <- ComplexHeatmap::Legend(
      title = title,
      labels = clevels,
      legend_gp = grid::gpar(fill = get_palette(clevels, palette = palette, palcolor = palcolor)),
      border = TRUE, nrow = if (legend.direction == "horizontal") 1 else NULL
    )
  } else {
    lgd <- NULL
  }
  anno <- ComplexHeatmap::AnnotationFunction(
    fun = function(index, k, n) {
      if (which == "row") {
        index <- rev(index)
      }
      grobs <- grobs[index]
      total <- length(index)
      # draw border
      grid.lines(x = c(0, 0), y = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(1, 1), y = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(0, 1), y = c(0, 0), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(0, 1), y = c(1, 1), gp = grid::gpar(col = "black", lwd = 1))
      for (i in seq_along(grobs)) {
        if (which == "row") {
          grobs[[i]]$vp <- viewport(x = 0.5, y = (i - 1) * 1 / total + 1 / (2 * total), width = 0.95, height = 1 / total)
        } else {
          grobs[[i]]$vp <- viewport(x = (i - 1) * 1 / total + 1 / (2 * total), y = 0.5, width = 1 / total, height = 1)
        }
        grid.draw(grobs[[i]])
      }
    },
    var_import = list(grobs = plots, which = which),
    n = length(plots),
    which = which,
    subsettable = TRUE,
    ...
  )
  list(anno = anno, legend = lgd)
}

#' Heatmap annotation functions
#'
#' @rdname heatmap-anno
#' @param x A data frame
#' @param split_by A character string of the column name to split the data
#' @param group_by A character string of the column name to group the data
#' @param column A character string of the column name to plot
#' @param title A character string to name the legend
#' @param which A character string specifying the direction of the annotation. Default is "row".
#'  Other options are "column".
#' @param palette A character string specifying the palette of the annotation
#' @param palcolor A character vector of colors to override the palette
#' @param border A logical value indicating whether to draw the border of the annotation
#' @param legend.direction A character string specifying the direction of the legend. Default is "vertical".
#'  Other options are "horizontal".
#' @param show_legend A logical value indicating whether to show the legend
#' @param .plotting A function to create the plot for each split and each group
#' @param ... Other arguments passed to `ComplexHeatmap::AnnotationFunction`
#' @keywords internal
#' @importFrom grid grid.draw grid.lines viewport gpar
#' @importFrom tidyr unite
.anno_ggseries <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                           palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, .plotting, ...) {
  x$.x <- x[[group_by]]
  glevels <- levels(x[[group_by]])
  colors <- get_palette(glevels, palette = palette, palcolor = palcolor)
  if (!is.null(split_by)) {
    x <- x %>% unite("..split", split_by, group_by, remove = FALSE)
    # We need to use show the original group_by in the legend
    # but here we have united the split_by and group_by
    x <- x %>% mutate(..palcolors = colors[as.numeric(x[[group_by]])])
    palcolors <- x$..palcolors
    names(palcolors) <- x$..split
    plots <- .plotting(data = x, column = column, group_by = "..split", palette = palette, palcolor = as.list(palcolors))
    plots <- lapply(plots, gggrob)
  } else {
    plots <- .plotting(data = x, column = column, group_by = group_by, palette = palette, palcolor = as.list(colors))
    plots <- lapply(plots, gggrob)
  }

  # add legend
  if (isTRUE(show_legend)) {
    lgd <- ComplexHeatmap::Legend(
      title = title,
      labels = glevels,
      legend_gp = grid::gpar(fill = get_palette(glevels, palette = palette, palcolor = colors)),
      border = TRUE, nrow = if (legend.direction == "horizontal") 1 else NULL
    )
  } else {
    lgd <- NULL
  }

  anno <- ComplexHeatmap::AnnotationFunction(
    fun = function(index, k, n) {
      if (identical(which, "row")) {
        index <- rev(index)
      }
      grobs <- grobs[index]
      total <- length(index)
      # draw border
      grid.lines(x = c(0, 0), y = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(1, 1), y = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(0, 1), y = c(0, 0), gp = grid::gpar(col = "black", lwd = 1))
      grid.lines(x = c(0, 1), y = c(1, 1), gp = grid::gpar(col = "black", lwd = 1))
      for (i in seq_along(grobs)) {
        if (which == "row") {
          grobs[[i]]$vp <- viewport(x = 0.5, y = (i - 1) * 1 / total + 1 / (2 * total), width = 0.95, height = 1 / total)
        } else {
          grobs[[i]]$vp <- viewport(x = (i - 1) * 1 / total + 1 / (2 * total), y = 0.5, width = 1 / total, height = 0.95)
        }
        grid.draw(grobs[[i]])
      }
    },
    var_import = list(grobs = plots, which = which),
    n = length(plots),
    which = which,
    subsettable = TRUE,
    ...
  )

  list(anno = anno, legend = lgd)
}

#' @rdname heatmap-anno
#' @param x A data frame
#' @param split_by A character string of the column name to split the data (heatmap)
#' @param group_by A character string of the column name to group the data (rows or columns of the heatmap)
#' @param column A character string of the column name of the data `x` to plot
#' @param title A character string to name the legend
#' @param which A character string specifying the direction of the annotation. Default is "row".
#'  Other options are "column".
#' @param palette A character string specifying the palette of the annotation
#' @param palcolor A character vector of colors to override the palette
#' @param border A logical value indicating whether to draw the border of the annotation
#' @param legend.direction A character string specifying the direction of the legend. Default is "vertical".
#'  Other options are "horizontal".
#' @param show_legend A logical value indicating whether to show the legend
#' @param ... Other arguments passed to `ComplexHeatmap::AnnotationFunction`
#'  The parameters passed to `row_annotation_params` and `column_annotation_params` will be passed here.
anno_pie <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                     palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggcat(
    x = x, split_by = split_by, group_by = group_by, column = column, title = title, which = which,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      PieChart(data, x = column, split_by = group_by, palette = palette, palcolor = palcolor, combine = FALSE)
    }, ...
  )
}

#' @rdname heatmap-anno
anno_ring <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                      palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggcat(
    x = x, split_by = split_by, group_by = group_by, column = column, which = which, title = title,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      RingPlot(data, group_by = column, split_by = group_by, palette = palette, palcolor = palcolor, combine = FALSE)
    }, ...
  )
}

#' @rdname heatmap-anno
anno_bar <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                     palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggcat(
    x = x, split_by = split_by, group_by = group_by, column = column, which = which, title = title,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      BarPlot(data,
        x = column, split_by = group_by, expand = c(0.05, 1),
        palette = palette, palcolor = palcolor, combine = FALSE
      )
    }, ...
  )
}

#' @rdname heatmap-anno
anno_violin <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                        palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggseries(
    x = x, split_by = split_by, group_by = group_by, column = column, which = which, title = title,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      ViolinPlot(data,
        x = ".x", y = column, split_by = group_by, combine = FALSE,
        palette = palette, palcolor = palcolor, flip = which == "row"
      )
    }, ...
  )
}

#' @rdname heatmap-anno
anno_boxplot <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                         palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggseries(
    x = x, split_by = split_by, group_by = group_by, column = column, which = which, title = title,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      BoxPlot(data,
        x = ".x", y = column, split_by = group_by, combine = FALSE,
        palette = palette, palcolor = palcolor, flip = which == "row"
      )
    }, ...
  )
}

#' @rdname heatmap-anno
anno_density <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                         palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, ...) {
  .anno_ggseries(
    x = x, split_by = split_by, group_by = group_by, column = column, title = title, which = which,
    palette = palette, palcolor = palcolor, border = border, legend.direction = legend.direction, show_legend = show_legend,
    .plotting = function(data, column, group_by, palette, palcolor) {
      DensityPlot(data,
        x = column, split_by = group_by, combine = FALSE,
        palette = palette, palcolor = palcolor, flip = which == "row"
      )
    }, ...
  )
}

#' @rdname heatmap-anno
#' @param alpha A numeric value between 0 and 1 specifying the transparency of the annotation
anno_simple <- function(x, split_by = NULL, group_by = NULL, column = NULL, title, which = "row", palette,
                        palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, alpha = 1, ...) {
  if (!is.null(split_by)) {
    x <- do.call(rbind, split(x, x[[split_by]]))
  }
  if (!is.null(column)) {
    x <- x[[column]]
  }
  is_cont <- is.numeric(x)
  if (isFALSE(is_cont) && !is.factor(x)) {
    x <- factor(x, levels = unique(x))
  }
  # add legend
  if (is_cont) {
    col_fun <- circlize::colorRamp2(
      seq(min(x), max(x), length = 100),
      get_palette(palette = palette, palcolor = palcolor, alpha = alpha)
    )
    lgd <- if (isTRUE(show_legend)) {
      ComplexHeatmap::Legend(
        title = title,
        col_fun = col_fun,
        border = TRUE, direction = legend.direction
      )
    }
    anno <- ComplexHeatmap::anno_simple(x, col = col_fun, which = which, border = border, ...)
  } else {
    colors <- get_palette(levels(x), palette = palette, palcolor = palcolor, alpha = alpha)
    lgd <- if (isTRUE(show_legend)) {
      ComplexHeatmap::Legend(
        title = title,
        labels = levels(x),
        legend_gp = grid::gpar(fill = colors),
        border = TRUE, nrow = if (legend.direction == "horizontal") 1 else NULL
      )
    }
    anno <- ComplexHeatmap::anno_simple(as.character(x), col = colors, which = which, border = border, ...)
  }

  list(anno = anno, legend = lgd)
}

#' @rdname heatmap-anno
anno_points <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                        palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, alpha = 1, ...) {
  if (!is.null(split_by)) {
    x <- do.call(rbind, split(x, x[[split_by]]))
  }

  anno <- ComplexHeatmap::anno_points(
    x[[column]],
    which = which, border = border, ...
  )
  list(anno = anno, legend = NULL)
}

#' @rdname heatmap-anno
#' @param add_points A logical value indicating whether to add points to the annotation
anno_lines <- function(x, split_by = NULL, group_by, column, title, which = "row", palette,
                       palcolor = NULL, border = TRUE, legend.direction, show_legend = TRUE, alpha = 1, add_points = TRUE, ...) {
  anno <- ComplexHeatmap::anno_lines(
    x[[column]],
    which = which, border = border, add_points = add_points, ...
  )
  list(anno = anno, legend = NULL)
}

#' Heatmap layer functions used to draw on the heatmap cells
#'
#' @rdname heatmap-layer
#' @param j An integer specifying the column index
#' @param i An integer specifying the row index
#' @param x A numeric vector specifying the x position
#' @param y A numeric vector specifying the y position
#' @param w A numeric vector specifying the width
#' @param h A numeric vector specifying the height
#' @param fill A character vector specifying the fill color
#' @keywords internal
layer_white_bg <- function(j, i, x, y, w, h, fill) {
  grid.rect(x = x, y = y, width = w, height = h, gp = grid::gpar(fill = "white", col = "white", lwd = 1))
}

#' @rdname heatmap-layer
#' @param alpha A numeric value between 0 and 1 specifying the transparency of the fill color
#' @keywords internal
layer_bg <- function(j, i, x, y, w, h, fill, alpha) {
  grid.rect(x = x, y = y, width = w, height = h, gp = grid::gpar(col = fill, lwd = 1, fill = adjcolors(fill, alpha)))
}

#' @rdname heatmap-layer
#' @param color A character vector specifying the color of the reticle
#' @keywords internal
layer_reticle <- function(j, i, x, y, w, h, fill, color) {
  ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
  ind_top <- ind_mat[1, ]
  ind_left <- ind_mat[, 1]
  for (col in seq_len(ncol(ind_mat))) {
    grid.lines(
      x = unit(rep(x[ind_top[col]], each = 2), "npc"), y = unit(c(0, 1), "npc"),
      gp = grid::gpar(col = color, lwd = 1.5)
    )
  }
  for (row in seq_len(nrow(ind_mat))) {
    grid.lines(
      x = unit(c(0, 1), "npc"), y = unit(rep(y[ind_left[row]], each = 2), "npc"),
      gp = grid::gpar(col = color, lwd = 1.5)
    )
  }
}

#' @rdname heatmap-layer
#' @param data A dataframe used to create the annotation. Different from the data used to
#'  create the heatmap itself, which is aggregated data. This dataframe is the original data,
#'  where each cell could have multiple values.
#' @param dot_size A numeric value specifying the size of the dot or a function to calculate the size
#'  from the values in the cell. The function can take 1, 3, or 5 arguments: the first argument is
#'  the values in the cell before aggregation; the 2nd and 3rd arguments are the row and column
#'  indices; the 4th and 5th arguments are the row and column names.
#' @param row_names Row names from the heatmap matrix.
#' @param col_names Column names from the heatmap matrix.
#' @keywords internal
layer_dot <- function(j, i, x, y, w, h, fill, data, dot_size, alpha, row_names = NULL, col_names = NULL) {
  if (is.numeric(dot_size) && length(dot_size) == 1) {
    # Simple numeric size
    grid.points(x, y,
      pch = 21, size = unit(dot_size, "mm"),
      gp = grid::gpar(col = "black", lwd = 1, fill = adjcolors(fill, alpha))
    )
  } else {
    # dot_size is a pre-computed vector/list or will be computed from function
    grid.points(x, y,
      pch = 21, size = unit(scales::rescale(unlist(dot_size), to = c(.5, 12)), "mm"),
      gp = grid::gpar(col = "black", lwd = 1, fill = adjcolors(fill, alpha))
    )
  }
}

#' @rdname heatmap-layer
#' @param col_fun A function to calculate the color of the bars
#' @keywords internal
layer_bars <- function(j, i, x, y, w, h, fill, flip, col_fun, data, alpha) {
  # colors be like [1,1.9]: '#A6CEE3' (1.9,2.8]: '#1F78B4' (2.8,3.7]: '#B2DF8A'
  indices <- paste(i, j, sep = "-")
  data <- data[indices]
  ns <- lengths(data)
  if (flip) {
    # rep(w / ns, ns) can't keep the unit
    bw <- rep(w, lengths(data))
    bh <- rep(sapply(seq_along(h), function(m) h[m] / ns[m]), ns)
    by <- unlist(lapply(seq_along(y), function(m) {
      y[m] - h[m] / 2 + seq_along(data[[m]]) * h[m] / length(data[[m]])
    })) - bh / 2
    bx <- rep(x, lengths(data))
  } else {
    # rep(w / ns, ns) can't keep the unit
    bw <- rep(sapply(seq_along(w), function(m) w[m] / ns[m]), ns)
    bh <- rep(h, lengths(data))
    bx <- unlist(lapply(seq_along(x), function(m) {
      x[m] - w[m] / 2 + seq_along(data[[m]]) * w[m] / length(data[[m]])
    })) - bw / 2
    by <- rep(y, lengths(data))
  }
  bf <- unlist(lapply(data, col_fun))
  grid.rect(x = bx, y = by, width = bw, height = bh, gp = grid::gpar(fill = bf, col = "transparent"))
}

#' @rdname heatmap-layer
#' @keywords internal
layer_pie <- function(j, i, x, y, w, h, fill, palette, palcolor, data, pie_size) {
  indices <- paste(i, j, sep = "-")
  data <- data[indices]
  pies <- lapply(indices, function(idx) {
    p <- PieChart(data[[idx]], x = "Var", y = "Freq", label = NULL, palette = palette, palcolor = palcolor)
    ggplotGrob(p + theme_void() + theme(legend.position = "none"))
  })

  if (!is.function(pie_size)) {
    pie_sizes <- rep(pie_size %||% 1, length(pies))
  } else {
    pie_sizes <- sapply(data, function(d) pie_size(sum(d$Freq, na.rm = TRUE)))
    pie_sizes <- scales::rescale(pie_sizes, to = c(0.2, 1))
  }
  idx <- which(sapply(pies[[1]]$grobs, function(g) inherits(g, "gTree") && !inherits(g, "zeroGrob") && !inherits(g, "absoluteGrob")))[1]
  for (m in seq_along(pies)) {
    pies[[m]]$grobs[[idx]]$vp <- viewport(x = x[m], y = y[m], width = pie_sizes[m] * w[m], height = pie_sizes[m] * h[m])
    grid.draw(pies[[m]]$grobs[[idx]])
  }
}

#' @rdname heatmap-layer
#' @param colors A character vector specifying the fill color of the violin plot.
#'  If not provided, the fill color of row/column annotation will be used
#' @keywords internal
layer_boxviolin <- function(j, i, x, y, w, h, fill, flip, data, colors, fn) {
  vlndata <- data[paste(i, j, sep = "-")]
  vlnplots <- lapply(seq_along(vlndata), function(m) {
    p <- fn(data.frame(x = 1, y = vlndata[[m]]), x = "x", y = "y", palcolor = colors %||% fill[m], flip = flip)
    ggplotGrob(p + theme_void() + theme(legend.position = "none"))
  })
  idx <- which(sapply(vlnplots[[1]]$grobs, function(g) inherits(g, "gTree") && !inherits(g, "zeroGrob") && !inherits(g, "absoluteGrob")))[1]
  for (m in seq_along(vlnplots)) {
    wm <- if (flip) w[m] * 0.95 else w[m]
    hm <- if (flip) h[m] * 0.95 else h[m]
    vlnplots[[m]]$grobs[[idx]]$vp <- viewport(x = x[m], y = y[m], width = wm, height = hm)
    grid.draw(vlnplots[[m]]$grobs[[idx]])
  }
}

#' Atomic heatmap without split
#'
#' @inheritParams process_heatmap_data
#' @inheritParams parameters
#' @param data A data frame used to create the heatmap.
#'  The data should be in a long form where each row represents a instance in the heatmap.
#' @param values_fill A value to fill in the missing values in the heatmap.
#' When there is missing value in the data, the cluster_rows and cluster_columns will fail.
#' @param border A logical value indicating whether to draw the border of the heatmap.
#'  If TRUE, the borders of the slices will be also drawn.
#' @param title The global (column) title of the heatmap
#' @param lower_quantile,upper_quantile,lower_cutoff,upper_cutoff Vector of minimum and maximum cutoff values or quantile values for each feature.
#'  It's applied to aggregated values when aggregated values are used (e.g. plot_type tile, label, etc).
#'  It's applied to raw values when raw values are used (e.g. plot_type bars, etc).
#' @param rows_palette A character string specifying the palette of the row group annotation.
#'  The default is "Paired".
#' @param rows_palcolor A character vector of colors to override the palette of the row group annotation.
#' @param columns_palette A character string specifying the palette of the column group annotation.
#'  The default is "Paired".
#' @param columns_palcolor A character vector of colors to override the palette of the column group annotation.
#' @param columns_split_palette A character string specifying the palette of the column split annotation.
#'  The default is "simspec".
#' @param columns_split_palcolor A character vector of colors to override the palette of the column split annotation.
#' @param rows_split_palette A character string specifying the palette of the row split annotation.
#'  The default is "simspec".
#' @param rows_split_palcolor A character vector of colors to override the palette of the row split annotation.
#' @param cluster_columns A logical value indicating whether to cluster the columns.
#'  If TRUE and columns_split_by is provided, the clustering will only be applied to the columns within the same split.
#' @param cluster_rows A logical value indicating whether to cluster the rows.
#'  If TRUE and rows_split_by is provided, the clustering will only be applied to the rows within the same split.
#' @param legend_items A numeric vector with names to specify the items in the main legend.
#'  The names will be working as the labels of the legend items.
#' @param legend_discrete A logical value indicating whether the main legend is discrete.
#' @param show_row_names A logical value indicating whether to show the row names.
#'  If TRUE, the legend of the row group annotation will be hidden.
#' @param show_column_names A logical value indicating whether to show the column names.
#'  If TRUE, the legend of the column group annotation will be hidden.
#' @param column_title A character string/vector of the column name(s) to use as the title of the column group annotation.
#' @param row_title A character string/vector of the column name(s) to use as the title of the row group annotation.
#' @param na_col A character string specifying the color for missing values.
#'  The default is "grey85".
#' @param row_names_side A character string specifying the side of the row names.
#'  The default is "right".
#' @param column_names_side A character string specifying the side of the column names.
#'  The default is "bottom".
#' @param bars_sample An integer specifying the number of samples to draw the bars.
#' @param flip A logical value indicating whether to flip the heatmap.
#' The idea is that, you can simply set `flip = TRUE` to flip the heatmap.
#' You don't need to swap the arguments related to rows and columns, except those you specify via `...`
#' that are passed to `ComplexHeatmap::Heatmap()` directly.
#' @param pie_palette A character string specifying the palette of the pie chart.
#' @param pie_palcolor A character vector of colors to override the palette of the pie chart.
#' @param pie_values A function or character that can be converted to a function by [match.arg()]
#' to calculate the values for the pie chart. Default is "length".
#' The function should take a vector of values as the argument and return a single value, for each
#' group in `pie_group_by`.
#' @param pie_size A numeric value or a function specifying the size of the pie chart.
#'  If it is a function, the function should take `count` as the argument and return the size.
#' @param pie_size_name A character string specifying the name of the legend for the pie size.
#' @param label_size A numeric value specifying the size of the labels when `cell_type = "label"`.
#' @param label A function to calculate the labels for the heatmap cells.
#' It can take either 1, 3, or 5 arguments. The first argument is the aggregated values.
#' If it takes 3 arguments, the second and third arguments are the row and column indices.
#' If it takes 5 arguments, the second and third arguments are the row and column indices,
#' the fourth and fifth arguments are the row and column names.
#' The function should return a character vector of the same length as the aggregated values.
#' If the function returns NA, no label will be shown for that cell.
#' For the indices, if you have the same dimension of data (same order of rows and columns) as the heatmap, you need to use `ComplexHeatmap::pindex()` to get the correct values.
#' @param layer_fun_callback A function to add additional layers to the heatmap.
#'  The function should have the following arguments: `j`, `i`, `x`, `y`, `w`, `h`, `fill`, `sr` and `sc`.
#'  Please also refer to the `layer_fun` argument in `ComplexHeatmap::Heatmap`.
#' @param cell_type A character string specifying the type of the heatmap cells.
#'  The default is values. Other options are "bars", "label", "dot", "violin", "boxplot".
#'  Note that for pie chart, the values under columns specified by `rows` will not be used directly. Instead, the values
#'  will just be counted in different `pie_group_by` groups. `NA` values will not be counted.
#' @param cell_agg A function to aggregate the values in the cell, for the cell type "tile" and "label".
#'  The default is `mean`.
#' @param add_bg A logical value indicating whether to add a background to the heatmap.
#'  Does not work with `cell_type = "bars"` or `cell_type = "tile"`.
#' @param bg_alpha A numeric value between 0 and 1 specifying the transparency of the background.
#' @param violin_fill A character vector of colors to override the fill color of the violin plot.
#'  If NULL, the fill color will be the same as the annotion.
#' @param boxplot_fill A character vector of colors to override the fill color of the boxplot.
#'  If NULL, the fill color will be the same as the annotion.
#' @param dot_size A numeric value specifying the size of the dot or a function to calculate the size
#'  from the values in the cell or a function to calculate the size from the values in the cell.
#' @param dot_size_name A character string specifying the name of the legend for the dot size.
#' If NULL, the dot size legend will not be shown.
#' @param column_name_annotation A logical value indicating whether to add the column annotation for the column names.
#'  which is a simple annotaion indicating the column names.
#' @param column_name_legend A logical value indicating whether to show the legend of the column name annotation.
#' @param column_annotation A character string/vector of the column name(s) to use as the column annotation.
#'  Or a list with the keys as the names of the annotation and the values as the column names.
#' @param column_annotation_side A character string specifying the side of the column annotation.
#'  Could be a list with the keys as the names of the annotation and the values as the sides.
#' @param column_annotation_palette A character string specifying the palette of the column annotation.
#'  The default is "Paired".
#'  Could be a list with the keys as the names of the annotation and the values as the palettes.
#' @param column_annotation_palcolor A character vector of colors to override the palette of the column annotation.
#'  Could be a list with the keys as the names of the annotation and the values as the palcolors.
#' @param column_annotation_type A character string specifying the type of the column annotation.
#'  The default is "auto". Other options are "simple", "pie", "ring", "bar", "violin", "boxplot", "density".
#'  Could be a list with the keys as the names of the annotation and the values as the types.
#'  If the type is "auto", the type will be determined by the type and number of the column data.
#' @param column_annotation_params A list of parameters passed to the annotation function.
#'  Could be a list with the keys as the names of the annotation and the values as the parameters passed to the annotation function. For the parameters for names (columns_by, rows_by, columns_split_by, rows_split_by), the key should be "name.(name)", where `(name)` is the name of the annotation.
#'  See [anno_pie()], [anno_ring()], [anno_bar()], [anno_violin()], [anno_boxplot()], [anno_density()], [anno_simple()], [anno_points()] and [anno_lines()] for the parameters of each annotation function.
#' @param column_annotation_agg A function to aggregate the values in the column annotation.
#' @param row_name_annotation A logical value indicating whether to add the row annotation for the row names.
#'  which is a simple annotaion indicating the row names.
#' @param row_name_legend A logical value indicating whether to show the legend of the row name annotation.
#' @param row_annotation A character string/vector of the column name(s) to use as the row annotation.
#' Or a list with the keys as the names of the annotation and the values as the column names.
#' @param row_annotation_side A character string specifying the side of the row annotation.
#' Could be a list with the keys as the names of the annotation and the values as the sides.
#' @param row_annotation_palette A character string specifying the palette of the row annotation.
#' The default is "Paired".
#' Could be a list with the keys as the names of the annotation and the values as the palettes.
#' @param row_annotation_palcolor A character vector of colors to override the palette of the row annotation.
#' Could be a list with the keys as the names of the annotation and the values as the palcolors.
#' @param row_annotation_type A character string specifying the type of the row annotation.
#' The default is "auto". Other options are "simple", "pie", "ring", "bar", "violin", "boxplot", "density".
#' Could be a list with the keys as the names of the annotation and the values as the types.
#' If the type is "auto", the type will be determined by the type and number of the row data.
#' @param row_annotation_params A list of parameters passed to the annotation function.
#' Could be a list with the keys as the names of the annotation and the values as the parameters.
#' Same as `column_annotation_params`.
#' @param row_annotation_agg A function to aggregate the values in the row annotation.
#' @param add_reticle A logical value indicating whether to add a reticle to the heatmap.
#' @param reticle_color A character string specifying the color of the reticle.
#' @param cell_border A logical value indicating whether to draw borders around heatmap cells (only for `cell_type = "tile"`). Default is `FALSE`.
#' @param cell_border_color A character string specifying the color of the cell borders. Default is `"grey80"`.
#' @param cell_border_width A numeric value specifying the width of the cell borders. Default is `0.1`.
#' @param palette A character string specifying the palette of the heatmap cells.
#' @param palcolor A character vector of colors to override the palette of the heatmap cells.
#' @param scale A character string specifying the scaling method for the heatmap matrix.
#'  Options are "none" (no scaling, default), "row" (z-score by row), or "column" (z-score by column).
#'  Row scaling is useful for comparing values across columns (e.g., cell types),
#'  while column scaling is useful for comparing values across rows (e.g., genes/pathways).
#' @param alpha A numeric value between 0 and 1 specifying the transparency of the heatmap cells.
#' @param return_grob A logical value indicating whether to return the grob object of the heatmap.
#'  This is useful when merging multiple heatmaps using patchwork.
#' @param ... Other arguments passed to [ComplexHeatmap::Heatmap()]
#' When `row_names_max_width` is passed, a unit is expected. But you can also pass a numeric values,
#' with a default unit "inches", or a string like "5inches" to specify the number and unit directly.
#' @importFrom dplyr group_by across ungroup %>% all_of summarise first slice_sample everything group_map
#' @importFrom tidyr pivot_longer pivot_wider unite expand_grid
#' @importFrom ggplot2 ggplotGrob theme_void
#' @importFrom grid grid.rect grid.text grid.lines grid.points viewport gpar unit grid.draw is.unit
#' @importFrom grid convertUnit grid.grabExpr
#' @return A drawn HeatmapList object if `return_grob = FALSE`. Otherwise, a grob/gTree object.
#' @keywords internal
HeatmapAtomic <- function(
    data, values_by, values_fill = NA,
    # data definition
    rows_by = NULL, rows_split_by = NULL,
    columns_by = NULL, columns_split_by = NULL,
    # palettes
    palette = "RdBu", palcolor = NULL,
    rows_palette = "Paired", rows_palcolor = NULL, rows_split_palette = "simspec", rows_split_palcolor = NULL,
    columns_palette = "Paired", columns_palcolor = NULL, columns_split_palette = "simspec", columns_split_palcolor = NULL,
    # cell_type: pies
    pie_size_name = "size", pie_size = NULL, pie_values = "length",
    pie_group_by = NULL, pie_palette = "Spectral", pie_palcolor = NULL,
    # cell_type: bars
    bars_sample = 100,
    # cell_type: label
    label = identity, label_size = 10,
    # cell_type: violin
    violin_fill = NULL,
    # cell_type: boxplot
    boxplot_fill = NULL,
    # cell_type: dot
    dot_size = 8, dot_size_name = "size",
    # legend
    legend_items = NULL, legend_discrete = FALSE,
    legend.position = "right", legend.direction = "vertical",
    # values
    scale = c("none", "row", "column"),
    lower_quantile = 0, upper_quantile = 0.99, lower_cutoff = NULL, upper_cutoff = NULL,
    # bg
    add_bg = FALSE, bg_alpha = 0.5,
    # reticle
    add_reticle = FALSE, reticle_color = "grey",
    # cell border (for tile type)
    cell_border = FALSE, cell_border_color = "grey80", cell_border_width = 0.1,
    # passed to ComplexHeatmap::Heatmap
    column_name_annotation = TRUE, column_name_legend = NULL,
    row_name_annotation = TRUE, row_name_legend = NULL,
    cluster_columns = TRUE, cluster_rows = TRUE, show_row_names = !row_name_annotation, show_column_names = !column_name_annotation,
    border = TRUE, title = NULL, column_title = character(0), row_title = character(0), na_col = "grey85",
    row_names_side = "right", column_names_side = "bottom",
    column_annotation = NULL, column_annotation_side = "top", column_annotation_palette = "Paired", column_annotation_palcolor = NULL,
    column_annotation_type = "auto", column_annotation_params = list(), column_annotation_agg = NULL,
    row_annotation = NULL, row_annotation_side = "left", row_annotation_palette = "Paired", row_annotation_palcolor = NULL,
    row_annotation_type = "auto", row_annotation_params = list(), row_annotation_agg = NULL,
    # misc
    flip = FALSE, alpha = 1, seed = 8525, return_grob = FALSE,
    # cell customization
    layer_fun_callback = NULL, cell_type = "tile", cell_agg = NULL,
    ...) {
  # NOTE: Data validation is performed in parent Heatmap() function

  # Normalize title parameters for ComplexHeatmap compatibility
  # FALSE -> NULL, TRUE -> character(0)
  if (isFALSE(column_title)) column_title <- NULL
  if (isFALSE(row_title)) row_title <- NULL
  if (isTRUE(column_title)) column_title <- character(0)
  if (isTRUE(row_title)) row_title <- character(0)

  get_col_fun <- function(lower, upper, a = alpha) {
    # If the lower and upper cutoff are the same, we need to adjust the upper cutoff
    if (upper == lower) {
      if (upper == 0) {
        upper <- 1e-3
      } else {
        upper <- upper + upper * 1e-3
      }
    }
    circlize::colorRamp2(
      seq(lower, upper, length = 100),
      get_palette(palette = palette, palcolor = palcolor, alpha = a, transparent = FALSE)
    )
  }

  flip_side <- function(side) {
    match.arg(side, c("left", "right", "top", "bottom"))
    if (side == "left") {
      return("top")
    }
    if (side == "right") {
      return("bottom")
    }
    if (side == "top") {
      return("left")
    }
    if (side == "bottom") {
      return("right")
    }
  }

  # Initialize the heatmap arguments
  hmargs <- list(
    # name = name,   # error when name has irregular characters (e.g. "-")
    heatmap_legend_param = list(title = values_by),
    border = border, na_col = na_col,
    cluster_columns = if (flip) cluster_rows else cluster_columns,
    cluster_rows = if (flip) cluster_columns else cluster_rows,
    cluster_column_slices = FALSE, cluster_row_slices = FALSE, show_heatmap_legend = FALSE,
    show_row_names = if (flip) show_column_names else show_row_names,
    show_column_names = if (flip) show_row_names else show_column_names,
    row_names_side = if (flip) flip_side(column_names_side) else row_names_side,
    column_names_side = if (flip) flip_side(row_names_side) else column_names_side,
    column_title = column_title, row_title = row_title,
    ...
  )

  # Set the row_names_max_width based on the length of the row names
  if (!is.null(hmargs$row_names_max_width)) {
    if (is.character(hmargs$row_names_max_width)) {
      if (grepl("^[0-9]+(\\.[0-9]+)?$", hmargs$row_names_max_width)) {
        hmargs$row_names_max_width <- unit(as.numeric(hmargs$row_names_max_width), "inches")
      }
      if (!grepl("^[0-9]+(\\.[0-9]+\\s*)?(npc|cm|centimetre|centimeter|in|inch|inches|mm|points|picas|bigpts|cicero|scalepts|lines|char|native|snpc|strwidth|strheight|grobwidth|grobheight|null)$", hmargs$row_names_max_width)) {
        stop(
          "[Heatmap] 'row_names_max_width' should be in a format of '2inches' or '2cm' if given as a string",
          call. = FALSE
        )
      }

      # convert to grid::unit
      hmargs$row_names_max_width <- unit(
        as.numeric(gsub("[^0-9.]", "", hmargs$row_names_max_width)),
        gsub("[0-9.]", "", hmargs$row_names_max_width)
      )
    } else if (is.numeric(hmargs$row_names_max_width)) {
      hmargs$row_names_max_width <- unit(hmargs$row_names_max_width, "inches")
    }
  } else if (flip) {
    hmargs$row_names_max_width <- ComplexHeatmap::max_text_width(levels(data[[columns_by]]))
  } else {
    hmargs$row_names_max_width <- ComplexHeatmap::max_text_width(levels(data[[rows_by]]))
  }

  # Collect the legends
  legends <- list()

  # Set the default cell aggregation function for pie chart (will be plotted as the background)
  cell_agg <- cell_agg %||% ifelse(cell_type == "pie", "nansum", "nanmean")
  if (is.character(cell_agg)) {
    if (startsWith(cell_agg, "nan")) {
      fn <- match.fun(substring(cell_agg, 4))
      cell_agg <- function(x) fn(x[is.finite(x)])
    } else {
      cell_agg <- match.fun(cell_agg)
    }
  }

  # Extract the matrix for the heatmap (aggregated values, for e.g. tile, label, pie background, etc)
  # We also need it for bars, because ComplexHeatmap::Heatmap need the matrix to plot anyway
  # rows_split_by  rows_by  columns_split_by1::columns_by1 columns_split_by2::columns_by2 ...
  # rows_split_by1 rows_by1 0.1                            0.2                            ...
  # rows_split_by2 rows_by2 0.3                            0.4                            ...
  # ...
  hmargs$matrix <- data %>%
    group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by)))) %>%
    summarise(.value = cell_agg(!!sym(values_by)), .groups = "drop") %>%
    unite(".columns", !!!syms(unique(c(columns_split_by, columns_by))), sep = " // ") %>%
    unite(".rows", !!!syms(unique(c(rows_split_by, rows_by))), sep = " // ") %>%
    pivot_wider(
      names_from = ".columns",
      values_from = ".value",
      values_fill = values_fill
    ) %>%
    as.data.frame()

  rownames(hmargs$matrix) <- hmargs$matrix$.rows
  hmargs$matrix$.rows <- NULL
  hmargs$matrix <- as.matrix(hmargs$matrix)
  hmargs$matrix[is.na(hmargs$matrix)] <- values_fill

  columns_order <- data %>%
    tidyr::expand(!!!syms(unique(c(columns_split_by, columns_by)))) %>%
    unite(".columns", !!!syms(unique(c(columns_split_by, columns_by))), sep = " // ") %>%
    dplyr::pull(".columns") %>%
    unique() %>%
    intersect(colnames(hmargs$matrix))
  rows_order <- data %>%
    tidyr::expand(!!!syms(unique(c(rows_split_by, rows_by)))) %>%
    unite(".rows", !!!syms(unique(c(rows_split_by, rows_by))), sep = " // ") %>%
    dplyr::pull(".rows") %>%
    unique() %>%
    intersect(rownames(hmargs$matrix))

  hmargs$matrix <- hmargs$matrix[rows_order, columns_order, drop = FALSE]

  # Apply scaling if requested
  scale <- match.arg(scale)
  if (scale == "row") {
    hmargs$matrix <- t(scale(t(hmargs$matrix)))
    hmargs$matrix[is.nan(hmargs$matrix)] <- 0
  } else if (scale == "column") {
    hmargs$matrix <- scale(hmargs$matrix)
    hmargs$matrix[is.nan(hmargs$matrix)] <- 0
  }

  if (flip) {
    hmargs$matrix <- t(hmargs$matrix)
  }

  r_split_by <- if (flip) columns_split_by else rows_split_by
  c_split_by <- if (flip) rows_split_by else columns_split_by
  r_by <- if (flip) columns_by else rows_by
  c_by <- if (flip) rows_by else columns_by
  r_split_levels <- if (!is.null(r_split_by)) levels(data[[r_split_by]])
  c_split_levels <- if (!is.null(c_split_by)) levels(data[[c_split_by]])
  r_levels <- if (!is.null(r_by)) levels(data[[r_by]])
  c_levels <- if (!is.null(c_by)) levels(data[[c_by]])
  if (!is.null(r_split_by)) {
    row_split_labels <- strsplit(rownames(hmargs$matrix), " // ", fixed = TRUE)
    hmargs$row_split <- factor(sapply(row_split_labels, `[`, 1), levels = r_split_levels)
    hmargs$row_labels <- factor(sapply(row_split_labels, `[`, 2), levels = r_levels)
  } else {
    hmargs$row_labels <- factor(rownames(hmargs$matrix), levels = r_levels)
  }

  if (!is.null(c_split_by)) {
    column_split_labels <- strsplit(colnames(hmargs$matrix), " // ", fixed = TRUE)
    hmargs$column_split <- factor(sapply(column_split_labels, `[`, 1), levels = c_split_levels)
    hmargs$column_labels <- factor(sapply(column_split_labels, `[`, 2), levels = c_levels)
  } else {
    hmargs$column_labels <- factor(colnames(hmargs$matrix), levels = c_levels)
  }

  if (cell_type == "bars") { # where multiple values are used, operating on data
    lower_cutoff <- lower_cutoff %||% quantile(data[[values_by]][is.finite(data[[values_by]])], lower_quantile, na.rm = TRUE)
    upper_cutoff <- upper_cutoff %||% quantile(data[[values_by]][is.finite(data[[values_by]])], upper_quantile, na.rm = TRUE)
    data[[values_by]][data[[values_by]] < lower_cutoff] <- lower_cutoff
    data[[values_by]][data[[values_by]] > upper_cutoff] <- upper_cutoff
  } else { # where aggregated values are used
    lower_cutoff <- lower_cutoff %||% quantile(hmargs$matrix[is.finite(hmargs$matrix)], lower_quantile, na.rm = TRUE)
    upper_cutoff <- upper_cutoff %||% quantile(hmargs$matrix[is.finite(hmargs$matrix)], upper_quantile, na.rm = TRUE)
  }

  # Set the color function for the heatmap cells
  hmargs$col <- get_col_fun(lower_cutoff, upper_cutoff)

  # Indices for data in layer_fun
  indices <- if (flip) {
    # 1-1, 2-1, 1-2, 2-2, ...
    expand.grid(1:nrow(hmargs$matrix), 1:ncol(hmargs$matrix))
  } else {
    # 1-1, 1-2, 2-1, 2-2, ...
    expand_grid(1:nrow(hmargs$matrix), 1:ncol(hmargs$matrix))
  }
  indices <- paste(indices[[1]], indices[[2]], sep = "-")

  # Compose the main legend
  get_main_legend <- function(allow_discreate = TRUE) {
    if (identical(legend.position, "none")) {
      return(NULL)
    }
    if (!allow_discreate && isTRUE(legend_discrete)) {
      stop("[Heatmap] 'legend_discrete = TRUE' is not allowed", call. = FALSE)
    }

    if (isTRUE(legend_discrete)) {
      if (is.null(legend_items)) {
        lgd_items <- sort(unique(as.vector(hmargs$matrix)), decreasing = TRUE)
        names(lgd_items) <- as.character(lgd_items)
      } else {
        lgd_items <- unlist(legend_items)
      }
      ComplexHeatmap::Legend(
        title = values_by, at = lgd_items, labels = names(lgd_items),
        legend_gp = grid::gpar(fill = hmargs$col(lgd_items)), border = TRUE, direction = legend.direction
      )
    } else {
      ComplexHeatmap::Legend(title = values_by, col_fun = hmargs$col, border = TRUE, direction = legend.direction)
    }
  }

  nrow_multiplier <- ncol_multiplier <- 1
  if (cell_type == "pie") {
    if (is.null(pie_group_by)) {
      stop("[Heatmap] Parameter 'pie_group_by' is required when 'cell_type = \"pie\"'", call. = FALSE)
    }
    pie_values <- pie_values %||% "length"
    if (is.character(pie_values)) {
      if (startsWith(pie_values, "nan")) {
        fn <- match.fun(substring(pie_values, 4))
        pie_values <- function(x) fn(x[is.finite(x)])
      } else {
        pie_values <- match.fun(pie_values)
      }
    }

    pie_group_levels <- levels(data[[pie_group_by]])
    pie_data <- data %>%
      group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by, pie_group_by))), .drop = FALSE) %>%
      summarise(.value = pie_values(!!sym(values_by)), .groups = "drop") %>%
      group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by)))) %>%
      group_map(
        ~ data.frame(Var = .x[[pie_group_by]], Freq = .x$.value)
      )
    names(pie_data) <- indices

    pie_colors <- get_palette(pie_group_levels, palette = pie_palette, palcolor = pie_palcolor)
    if (is.character((pie_size))) {
      if (startsWith(pie_size, "nan")) {
        fn <- match.fun(substring(pie_size, 4))
        pie_size <- function(x) fn(x[is.finite(x)])
      } else {
        pie_size <- match.fun(pie_size)
      }
    }
    hmargs$layer_fun <- function(j, i, x, y, w, h, fill, sr, sc) {
      layer_white_bg(j, i, x, y, w, h, fill)
      if (isTRUE(add_bg)) {
        layer_bg(j, x, x, y, w, h, fill, alpha = bg_alpha)
      }
      if (isTRUE(add_reticle)) {
        layer_reticle(j, i, x, y, w, h, fill, color = reticle_color)
      }
      layer_pie(
        j, i, x, y, w, h, fill,
        palette = pie_palette, palcolor = pie_palcolor, data = pie_data, pie_size = pie_size
      )
      if (is.function(layer_fun_callback)) {
        layer_fun_callback(j, i, x, y, w, h, fill, sr, sc)
      }
    }

    if (!identical(legend.position, "none") && is.function(pie_size)) {
      pie_sizes <- sapply(pie_data, function(d) pie_size(sum(d$Freq, na.rm = TRUE)))
      pie_size_min <- min(pie_sizes, na.rm = TRUE)
      pie_size_max <- max(pie_sizes, na.rm = TRUE)
      legends$.pie_size <- ComplexHeatmap::Legend(
        title = pie_size_name,
        labels = scales::number((seq(pie_size_min, pie_size_max, length.out = ifelse(pie_size_max > pie_size_min, 4, 1)))),
        type = "points",
        pch = 21,
        size = unit(8, "mm") * seq(0.2, 1, length.out = 4),
        grid_height = unit(8, "mm") * seq(0.2, 1, length.out = 4) * 0.8,
        grid_width = unit(8, "mm"),
        legend_gp = grid::gpar(fill = "grey30"),
        border = FALSE,
        background = "transparent",
        direction = legend.direction
      )
    }
    if (isTRUE(add_bg)) {
      legends$.heatmap <- get_main_legend()
    }
    if (!identical(legend.position, "none")) {
      legends$.pie <- ComplexHeatmap::Legend(
        title = pie_group_by, direction = legend.direction,
        border = TRUE, labels = pie_group_levels, legend_gp = grid::gpar(fill = pie_colors)
      )
    }
    nrow_multiplier <- 4
    ncol_multiplier <- 6
  } else if (cell_type == "bars") {
    if (isTRUE(add_bg)) {
      stop("[Heatmap] Cannot use 'add_bg' with 'cell_type = \"bars\"'", call. = FALSE)
    }
    if (isTRUE(add_reticle)) {
      stop("[Heatmap] Cannot use 'add_reticle' with 'cell_type = \"bars\"'", call. = FALSE)
    }

    bars_data <- data %>%
      group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by)))) %>%
      group_map(~ .x[[values_by]])

    names(bars_data) <- indices

    # plot bars in each cell
    hmargs$layer_fun <- function(j, i, x, y, w, h, fill, sr, sc) {
      layer_white_bg(j, i, x, y, w, h, fill)
      layer_bars(
        j, i, x, y, w, h, fill,
        flip = flip,
        col_fun = hmargs$col, data = bars_data, alpha = alpha
      )
      if (is.function(layer_fun_callback)) {
        layer_fun_callback(j, i, x, y, w, h, fill, sr, sc)
      }
    }
    # Override the main legend
    legends$.heatmap <- get_main_legend(FALSE)
    nrow_multiplier <- 0.5
  } else if (cell_type == "dot") {
    if (is.character(dot_size)) {
      if (startsWith(dot_size, "nan")) {
        fn <- match.fun(substring(dot_size, 4))
        dot_size <- function(x) fn(x[is.finite(x)])
      } else {
        dot_size <- match.fun(dot_size)
      }
    }
    # Store raw values for each cell to pass to dot_size function later
    dot_data <- data %>%
      group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by)))) %>%
      summarise(.value = list(!!sym(values_by)), .groups = "drop") %>%
      unite(".columns", !!!syms(unique(c(columns_split_by, columns_by))), sep = " // ") %>%
      unite(".rows", !!!syms(unique(c(rows_split_by, rows_by))), sep = " // ") %>%
      pivot_wider(names_from = ".columns", values_from = ".value") %>%
      select(-!!sym(".rows")) %>%
      as.data.frame()

    if (flip) {
      dot_data <- t(dot_data)
    }

    if (!identical(legend.position, "none") && is.function(dot_size) && !is.null(dot_size_name)) {
      # Optimized: only compute min/max for legend, not all sizes
      nargs <- length(formalArgs(dot_size))
      dot_size_min <- Inf
      dot_size_max <- -Inf

      for (idx in seq_along(indices)) {
        cell_key <- indices[idx]
        # Parse indices from the key "i-j"
        ij <- as.integer(strsplit(cell_key, "-")[[1]])
        cell_values <- dot_data[ij[1], ij[2]][[1]]

        size_val <- if (nargs == 1 || is.primitive(dot_size)) {
          dot_size(cell_values)
        } else if (nargs == 3) {
          dot_size(cell_values, ij[1], ij[2])
        } else if (nargs == 5) {
          dot_size(
            cell_values, ij[1], ij[2],
            rownames(hmargs$matrix)[ij[1]],
            colnames(hmargs$matrix)[ij[2]]
          )
        } else {
          stop("[Heatmap] 'dot_size' function should take 1, 3, or 5 arguments", call. = FALSE)
        }

        if (is.finite(size_val)) {
          if (size_val < dot_size_min) dot_size_min <- size_val
          if (size_val > dot_size_max) dot_size_max <- size_val
        }
      }

      legends$.dot_size <- ComplexHeatmap::Legend(
        title = dot_size_name,
        labels = scales::number((seq(dot_size_min, dot_size_max, length.out = ifelse(dot_size_max > dot_size_min, 4, 1)))),
        type = "points",
        pch = 21,
        size = unit(8, "mm") * seq(0.2, 1, length.out = 4),
        grid_height = unit(8, "mm") * seq(0.2, 1, length.out = 4) * 0.8,
        grid_width = unit(8, "mm"),
        legend_gp = grid::gpar(fill = "grey30"),
        border = FALSE,
        background = "transparent",
        direction = legend.direction
      )
    }

    hmargs$layer_fun <- function(j, i, x, y, w, h, fill, sr, sc) {
      layer_white_bg(j, i, x, y, w, h, fill)
      if (isTRUE(add_bg)) {
        layer_bg(j, i, x, y, w, h, fill, alpha = bg_alpha)
      }
      if (isTRUE(add_reticle)) {
        layer_reticle(j, i, x, y, w, h, fill, color = reticle_color)
      }
      # Compute dot sizes based on function arguments
      if (is.function(dot_size)) {
        nargs <- length(formalArgs(dot_size))
        sizes <- numeric(length(i))
        for (idx in seq_along(i)) {
          cell_values <- dot_data[[i[idx], j[idx]]]
          if (nargs == 1 || is.primitive(dot_size)) {
            sizes[idx] <- dot_size(cell_values)
          } else if (nargs == 3) {
            sizes[idx] <- dot_size(cell_values, i[idx], j[idx])
          } else if (nargs == 5) {
            sizes[idx] <- dot_size(
              cell_values, i[idx], j[idx],
              rownames(hmargs$matrix)[i[idx]],
              colnames(hmargs$matrix)[j[idx]]
            )
          } else {
            stop("[Heatmap] 'dot_size' function should take 1, 3, or 5 arguments", call. = FALSE)
          }
        }
        layer_dot(
          j, i, x, y, w, h, fill,
          data = dot_data, dot_size = sizes, alpha = alpha
        )
      } else {
        layer_dot(
          j, i, x, y, w, h, fill,
          data = dot_data, dot_size = dot_size, alpha = alpha
        )
      }
      if (is.function(layer_fun_callback)) {
        layer_fun_callback(j, i, x, y, w, h, fill, sr, sc)
      }
    }
    # Override the main legend
    legends$.heatmap <- get_main_legend()
  } else if (cell_type %in% c("violin", "boxplot")) {
    # df with multiple values in each cell
    vdata <- data %>%
      group_by(!!!syms(unique(c(rows_split_by, rows_by, columns_split_by, columns_by)))) %>%
      group_map(~ .x[[values_by]])

    names(vdata) <- indices
    vcolors <- if (cell_type == "violin") violin_fill else boxplot_fill

    hmargs$layer_fun <- function(j, i, x, y, w, h, fill, sr, sc) {
      layer_white_bg(j, i, x, y, w, h, fill)
      if (isTRUE(add_bg)) {
        layer_bg(j, i, x, y, w, h, fill, alpha = bg_alpha)
      }
      if (isTRUE(add_reticle)) {
        layer_reticle(j, i, x, y, w, h, fill, color = reticle_color)
      }
      layer_fn <- if (cell_type == "violin") ViolinPlot else BoxPlot
      layer_boxviolin(
        j, i, x, y, w, h, fill,
        flip = flip, data = vdata, colors = vcolors, fn = layer_fn
      )
      if (is.function(layer_fun_callback)) {
        layer_fun_callback(j, i, x, y, w, h, fill, sr, sc)
      }
    }
    if (!identical(legend.position, "none")) {
      if (!is.null(legend_items)) {
        stop("[Heatmap] Cannot use 'legend_items' with 'cell_type = \"violin\"' or '\"boxplot\"'", call. = FALSE)
      }
      if (isTRUE(legend_discrete)) {
        stop("[Heatmap] Cannot use 'legend_discrete = TRUE' with 'cell_type = \"violin\"' or '\"boxplot\"'", call. = FALSE)
      }
      if (is.null(vcolors)) {
        legends$.heatmap <- ComplexHeatmap::Legend(
          title = values_by, col_fun = hmargs$col, border = TRUE,
          direction = legend.direction
        )
      } else if (isTRUE(add_bg)) {
        legends$.heatmap <- ComplexHeatmap::Legend(
          title = values_by, col_fun = get_col_fun(lower_cutoff, upper_cutoff, bg_alpha),
          border = TRUE, direction = legend.direction
        )
      }
    }
  } else if (cell_type == "tile") {
    if (isTRUE(add_bg)) {
      stop("[Heatmap] Cannot use 'add_bg' with 'cell_type = \"tile\"'", call. = FALSE)
    }
    if (isTRUE(add_reticle)) {
      stop("[Heatmap] Cannot use 'add_reticle' with 'cell_type = \"tile\"'", call. = FALSE)
    }
    if (isTRUE(cell_border)) {
      hmargs$rect_gp <- grid::gpar(col = cell_border_color, lwd = cell_border_width)
    } else {
      hmargs$rect_gp <- grid::gpar(col = NA)
    }
    hmargs$layer_fun <- layer_fun_callback
    legends$.heatmap <- get_main_legend()
  } else if (cell_type == "label") {
    if (isTRUE(add_bg)) {
      stop("[Heatmap] Cannot use 'add_bg' with 'cell_type = \"label\"'", call. = FALSE)
    }
    if (isTRUE(add_reticle)) {
      stop("[Heatmap] Cannot use 'add_reticle' with 'cell_type = \"label\"'", call. = FALSE)
    }
    hmargs$layer_fun <- function(j, i, x, y, w, h, fill, sr, sc) {
      labels <- ComplexHeatmap::pindex(hmargs$matrix, i, j)
      if (is.function(label)) {
        nargs <- length(formalArgs(label))
        if (nargs == 1) {
          labels <- label(labels)
        } else if (nargs == 3) {
          labels <- label(labels, i, j)
        } else if (nargs == 5) {
          labels <- label(labels, i, j, rownames(hmargs$matrix)[i], colnames(hmargs$matrix)[j])
        } else {
          stop("[Heatmap] 'label' function should take 1, 3, or 5 arguments", call. = FALSE)
        }
      }
      inds <- !is.na(labels)
      if (any(inds)) {
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        lapply(theta, function(a) {
          x_out <- x[inds] + unit(cos(a) * label_size / 30, "mm")
          y_out <- y[inds] + unit(sin(a) * label_size / 30, "mm")
          grid.text(labels[inds], x = x_out, y = y_out, gp = grid::gpar(fontsize = label_size, col = "white"))
        })
        grid.text(labels[inds], x[inds], y[inds], gp = grid::gpar(fontsize = label_size, col = "black"))
      }
      if (is.function(layer_fun_callback)) {
        layer_fun_callback(j, i, x, y, w, h, fill, sr, sc)
      }
    }
    legends$.heatmap <- get_main_legend()
  }

  .ncols <- nlevels(data[[columns_by]])
  if (!is.null(columns_split_by)) {
    .ncols <- .ncols * nlevels(data[[columns_split_by]])
  }
  .ncols <- .ncols * ncol_multiplier
  .nrows <- nlevels(data[[rows_by]])
  if (!is.null(rows_split_by)) {
    .nrows <- .nrows * nlevels(data[[rows_split_by]])
  }
  .nrows <- .nrows * nrow_multiplier
  nrows <- ifelse(flip, .ncols, .nrows)
  ncols <- ifelse(flip, .nrows, .ncols)

  ## Set up the top annotations
  setup_name_annos <- function(names_side, anno_title, show_names, params, which,
                               split_by, splits, split_palette, split_palcolor,
                               by, by_name_annotation, by_labels, by_palette, by_palcolor, by_name_legend) {
    annos <- list(
      annotation_name_side = names_side,
      show_annotation_name = list()
    )
    if (!is.null(split_by)) {
      param <- params[[paste0("name.", split_by)]] %||% list()
      param$x <- splits
      param$title <- split_by
      param$palette <- split_palette
      param$palcolor <- split_palcolor
      param$border <- param$border %||% TRUE
      param$legend.direction <- legend.direction
      param$which <- ifelse(flip, setdiff(c("column", "row"), which), which)
      param$show_legend <- is.null(anno_title) && !identical(legend.position, "none")
      worh <- ifelse(param$which == "row", "width", "height")
      param[[worh]] <- param[[worh]] %||% unit(2.5, "mm")
      if (is.numeric(param[[worh]]) && !grid::is.unit(param[[worh]])) {
        param[[worh]] <- unit(param[[worh]], "mm")
      }

      annos$show_annotation_name[[split_by]] <- TRUE
      anno_legend <- do.call(anno_simple, param)
      annos[[split_by]] <- anno_legend$anno
      if (isTRUE(param$show_legend)) {
        legends[[paste0("name.", split_by)]] <<- anno_legend$legend
      }
    }

    if (!is.null(by) && by_name_annotation) {
      param <- params[[paste0("name.", by)]] %||% list()
      param$x <- by_labels
      param$title <- by
      param$palette <- by_palette
      param$palcolor <- by_palcolor
      param$border <- param$border %||% TRUE
      param$legend.direction <- legend.direction
      param$which <- ifelse(flip, setdiff(c("column", "row"), which), which)
      param$show_legend <- !identical(legend.position, "none") &&
        (isTRUE(by_name_legend) || (is.null(by_name_legend) && !show_names))
      worh <- ifelse(param$which == "row", "width", "height")
      param[[worh]] <- param[[worh]] %||% unit(2.5, "mm")
      if (is.numeric(param[[worh]]) && !grid::is.unit(param[[worh]])) {
        param[[worh]] <- unit(param[[worh]], "mm")
      }

      annos$show_annotation_name[[by]] <- TRUE
      anno_legend <- do.call(anno_simple, param)
      annos[[by]] <- anno_legend$anno
      if (isTRUE(param$show_legend)) {
        legends[[paste0("name.", by)]] <<- anno_legend$legend
      }
    }

    return(annos)
  }

  setup_annos <- function(which, names_side,
                          annotation, annotation_type,
                          annotation_palette, annotation_palcolor,
                          annotation_agg, annotation_params,
                          split_by, by) {
    annos <- list()
    if (!is.list(annotation)) {
      annotation <- as.list(annotation)
      names(annotation) <- unlist(annotation)
    }
    if (is.character(annotation_type) && identical(annotation_type, "auto")) {
      annotation_type <- as.list(rep("auto", length(annotation)))
      names(annotation_type) <- names(annotation)
    }
    if (is.character(annotation_palette) && length(annotation_palette) == 1) {
      annotation_palette <- as.list(rep(annotation_palette, length(annotation)))
      names(annotation_palette) <- names(annotation)
    }
    if (!is.list(annotation_palcolor)) {
      annotation_palcolor <- list(annotation_palcolor)
      annotation_palcolor <- rep(annotation_palcolor, length(annotation))
      names(annotation_palcolor) <- names(annotation)
    }
    annotation_agg <- annotation_agg %||% list()
    annotation_params <- annotation_params %||% list()
    for (aname in names(annotation)) {
      if (aname %in% formalArgs(ComplexHeatmap::HeatmapAnnotation)) {
        annos[[aname]] <- annotation[[aname]]
        next
      }
      annocol <- annotation[[aname]]
      annoagg <- annotation_agg[[aname]]
      annotype <- annotation_type[[aname]] %||% "auto"
      param <- annotation_params[[aname]] %||% list()
      annodata <- param$x %||% data
      annocol <- validate_columns(annodata, annocol)
      if (annotype == "auto") {
        all_ones <- annodata %>%
          group_by(!!!syms(unique(c(split_by, by)))) %>%
          summarise(n = n(), .groups = "drop") %>%
          pull("n")
        all_ones <- all(all_ones == 1)
        if (is.character(annodata[[annocol]]) || is.factor(annodata[[annocol]]) || is.logical(annodata[[annocol]])) {
          annotype <- ifelse(all_ones, "pie", "simple")
        } else if (is.numeric(annodata[[annocol]])) {
          annotype <- ifelse(all_ones, "points", "violin")
        } else {
          stop(
            sprintf("[Heatmap] Don't know how to handle %s annotation type for column: %s", which, annocol),
            call. = FALSE
          )
        }
      }
      if (annotype %in% c("simple", "points", "lines") && is.null(annoagg)) {
        warning("[Heatmap] Assuming '", which, "_annotation_agg[\"", aname, "\"] = dplyr::first' for the simple annotation", call. = FALSE)
        annoagg <- dplyr::first
      }
      if (is.null(annoagg)) {
        annodata <- annodata %>% select(!!!syms(unique(c(split_by, by, annocol))))
      } else {
        annodata <- annodata %>%
          group_by(!!!syms(unique(c(split_by, by)))) %>%
          summarise(!!sym(annocol) := annoagg(!!sym(annocol)), .groups = "drop")
      }
      param$x <- annodata
      param$split_by <- split_by
      param$group_by <- by
      param$column <- annocol
      param$title <- aname
      param$which <- ifelse(flip, setdiff(c("column", "row"), which), which)
      param$palette <- annotation_palette[[aname]] %||% "Paired"
      param$palcolor <- annotation_palcolor[[aname]]
      param$legend.direction <- legend.direction
      # swap width and height if flip is TRUE
      if (flip) {
        pheight <- param$height
        param$height <- param$width
        param$width <- pheight
      }
      if (legend.position == "none") {
        param$show_legend <- FALSE
      }
      if (!exists(paste0("anno_", annotype))) {
        stop(
          sprintf("[Heatmap] Unsupported annotation type: %s", annotype),
          call. = FALSE
        )
      }
      anno <- do.call(paste0("anno_", annotype), param)
      annos[[aname]] <- anno$anno
      legends[[paste0(which, ".", aname)]] <<- anno$legend
    }

    if (length(annos) > 0) {
      annos$annotation_name_side <- names_side
    }

    return(annos)
  }

  ncol_annos <- sum(cluster_columns, show_column_names) * 4
  ncol_annos <- ncol_annos +
    ifelse(is.null(columns_split_by), 0, 1) +
    ifelse(is.null(columns_by) || !column_name_annotation, 0, 1)
  top_annos <- setup_name_annos(
    names_side = ifelse(flip, column_names_side, row_names_side), anno_title = column_title,
    show_names = show_column_names, params = column_annotation_params, which = "column",
    split_by = columns_split_by, splits = if (flip) hmargs$row_split else hmargs$column_split,
    split_palette = columns_split_palette, split_palcolor = columns_split_palcolor,
    by = columns_by, by_name_annotation = column_name_annotation,
    by_labels = if (flip) hmargs$row_labels else hmargs$column_labels, by_palette = columns_palette,
    by_palcolor = columns_palcolor, by_name_legend = column_name_legend
  )
  column_annos <- setup_annos(
    which = "column", names_side = ifelse(flip, column_names_side, row_names_side),
    annotation = column_annotation, annotation_type = column_annotation_type,
    annotation_palette = column_annotation_palette, annotation_palcolor = column_annotation_palcolor,
    annotation_agg = column_annotation_agg, annotation_params = column_annotation_params,
    split_by = columns_split_by, by = columns_by
  )

  if (column_annotation_side == "top") {
    column_annos$annotation_name_side <- NULL
    top_annos <- c(top_annos, column_annos)
  } else if (length(column_annos) > 0) {
    if (isTRUE(flip)) {
      hmargs$right_annotation <- do.call(ComplexHeatmap::rowAnnotation, column_annos)
    } else {
      hmargs$bottom_annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, column_annos)
    }
    rm(column_annos)
  }
  # Check if top_annos has any actual annotation (excluding metadata fields)
  has_top_anno <- length(setdiff(names(top_annos), c("annotation_name_side", "show_annotation_name"))) > 0
  if (length(top_annos$show_annotation_name) > 0 || has_top_anno) {
    if (isTRUE(flip)) {
      hmargs$left_annotation <- do.call(ComplexHeatmap::rowAnnotation, top_annos)
    } else {
      hmargs$top_annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, top_annos)
    }
  }
  rm(top_annos)

  nrow_annos <- sum(cluster_rows, show_row_names) * 4
  nrow_annos <- nrow_annos +
    ifelse(is.null(rows_split_by), 0, 1) +
    ifelse(is.null(rows_by) || !row_name_annotation, 0, 1)
  left_annos <- setup_name_annos(
    names_side = ifelse(flip, row_names_side, column_names_side), anno_title = row_title,
    show_names = show_row_names, params = row_annotation_params, which = "row",
    split_by = rows_split_by, splits = if (flip) hmargs$column_split else hmargs$row_split,
    split_palette = rows_split_palette, split_palcolor = rows_split_palcolor,
    by = rows_by, by_name_annotation = row_name_annotation,
    by_labels = if (flip) hmargs$column_labels else hmargs$row_labels,
    by_palette = rows_palette, by_palcolor = rows_palcolor, by_name_legend = row_name_legend
  )
  row_annos <- setup_annos(
    which = "row", names_side = ifelse(flip, row_names_side, column_names_side),
    annotation = row_annotation, annotation_type = row_annotation_type,
    annotation_palette = row_annotation_palette, annotation_palcolor = row_annotation_palcolor,
    annotation_agg = row_annotation_agg, annotation_params = row_annotation_params,
    split_by = rows_split_by, by = rows_by
  )

  if (row_annotation_side == "left") {
    row_annos$annotation_name_side <- NULL
    left_annos <- c(left_annos, row_annos)
  } else if (length(row_annos) > 0) {
    if (isTRUE(flip)) {
      hmargs$bottom_annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, row_annos)
    } else {
      hmargs$right_annotation <- do.call(ComplexHeatmap::rowAnnotation, row_annos)
    }
    rm(row_annos)
  }
  # Check if left_annos has any actual annotation (excluding metadata fields)
  has_left_anno <- length(setdiff(names(left_annos), c("annotation_name_side", "show_annotation_name"))) > 0
  if (length(left_annos$show_annotation_name) > 0 || has_left_anno) {
    if (isTRUE(flip)) {
      hmargs$top_annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, left_annos)
    } else {
      hmargs$left_annotation <- do.call(ComplexHeatmap::rowAnnotation, left_annos)
    }
  }
  rm(left_annos)

  ## Set up the heatmap
  rownames_width <- grid::convertUnit(hmargs$row_names_max_width, "inches", valueOnly = TRUE) * 0.5 - 0.2
  rownames_width <- max(rownames_width, 0)
  if (isTRUE(flip)) {
    width <- nrows * 0.25 + ncol_annos * 0.5 + rownames_width
    # How about column name length (nchars)?
    height <- ncols * 0.25 + nrow_annos * 0.5
  } else {
    width <- ncols * 0.25 + nrow_annos * 0.5 + rownames_width
    height <- nrows * 0.25 + ncol_annos * 0.5
  }
  if (cell_type == "pie") {
    width <- max(width, height)
    height <- max(width, height)
  }
  if (!identical(legend.position, "none")) {
    if (legend.position %in% c("right", "left")) {
      if (legend.direction == "horizontal") {
        width <- width + 3
      } else {
        width <- width + 1.5
      }
    } else if (legend.direction == "horizontal") {
      height <- height + 3
    } else {
      height <- height + 1.5
    }
  }

  p <- do.call(ComplexHeatmap::Heatmap, hmargs)
  mat <- p@matrix
  # Default title style: centered and bold

  title_gp <- grid::gpar(fontsize = 14, fontface = "bold", just = "center")

  if (isTRUE(return_grob)) {
    if (identical(legend.position, "none")) {
      p <- grid.grabExpr(ComplexHeatmap::draw(p,
        annotation_legend_list = legends,
        show_annotation_legend = FALSE,
        column_title = title,
        column_title_gp = title_gp
      ))
    } else {
      p <- grid.grabExpr(ComplexHeatmap::draw(p,
        annotation_legend_list = legends,
        annotation_legend_side = legend.position,
        column_title = title,
        column_title_gp = title_gp
      ))
    }
  } else {
    # When return_grob = FALSE (ggplot2 v4), ComplexHeatmap::draw() will render
    # to the graphics device. To prevent unwanted output during assignment while
    # still allowing proper display when explicitly printed, we capture to a
    # null device first, then return the HeatmapList which can be drawn later.
    current_dev <- grDevices::dev.cur()
    null_dev <- grDevices::pdf(NULL)
    on.exit(
      {
        grDevices::dev.off() # Close the null device
        if (current_dev > 1) {
          grDevices::dev.set(current_dev) # Restore the previous device if it wasn't null
        }
      },
      add = TRUE
    )

    if (identical(legend.position, "none")) {
      p <- ComplexHeatmap::draw(p,
        annotation_legend_list = legends,
        show_annotation_legend = FALSE,
        column_title = title,
        column_title_gp = title_gp
      )
    } else {
      p <- ComplexHeatmap::draw(p,
        annotation_legend_list = legends,
        annotation_legend_side = legend.position,
        column_title = title,
        column_title_gp = title_gp
      )
    }
  }

  attr(p, "height") <- max(height, 4)
  attr(p, "width") <- max(width, 4)
  attr(p, "data") <- mat
  p
}

#' Heatmap
#'
#' @description Heatmap is a popular way to visualize data in matrix format. It is widely used in biology to visualize gene expression data in microarray and RNA-seq data. The heatmap is a matrix where rows represent the samples and columns represent the features. The color of each cell represents the value of the feature in the sample. The color can be continuous or discrete. The heatmap can be split by the columns or rows to show the subgroups in the data. The heatmap can also be annotated by the columns or rows to show the additional information of the samples or features.
#' @inheritParams process_heatmap_data
#' @inheritParams HeatmapAtomic
#' @inheritParams parameters
# @keywords internal
#' @importFrom patchwork wrap_plots
#' @seealso \code{\link{anno_simple}}, \code{\link{anno_points}}, \code{\link{anno_lines}}, \code{\link{anno_pie}}, \code{\link{anno_violin}}, \code{\link{anno_boxplot}}, \code{\link{anno_density}}
#' @examples
#' \donttest{
#' set.seed(8525)
#'
#' matrix_data <- matrix(rnorm(60), nrow = 6, ncol = 10)
#' rownames(matrix_data) <- paste0("R", 1:6)
#' colnames(matrix_data) <- paste0("C", 1:10)
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(matrix_data)
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # use a different color palette
#'   # change the main legend title
#'   # show row names (legend will be hidden)
#'   # show column names
#'   # change the row name annotation name and side
#'   # change the column name annotation name
#'   Heatmap(matrix_data,
#'     palette = "viridis", values_by = "z-score",
#'     show_row_names = TRUE, show_column_names = TRUE,
#'     rows_name = "Features", row_names_side = "left",
#'     columns_name = "Samples"
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # flip the heatmap
#'   Heatmap(matrix_data,
#'     palette = "viridis", values_by = "z-score",
#'     show_row_names = TRUE, show_column_names = TRUE,
#'     rows_name = "Features", row_names_side = "left",
#'     columns_name = "Samples", flip = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # add annotations to the heatmap
#'   rows_data <- data.frame(
#'     rows = paste0("R", 1:6),
#'     group = sample(c("X", "Y", "Z"), 6, replace = TRUE)
#'   )
#'   Heatmap(matrix_data,
#'     rows_data = rows_data,
#'     row_annotation = list(Group = "group"),
#'     row_annotation_type = list(Group = "simple"),
#'     row_annotation_palette = list(Group = "Spectral")
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(matrix_data,
#'     rows_data = rows_data,
#'     rows_split_by = "group"
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # add labels to the heatmap
#'   Heatmap(matrix_data,
#'     rows_data = rows_data,
#'     rows_split_by = "group", cell_type = "label",
#'     label = function(x) {
#'       ifelse(
#'         x > 0, scales::number(x, accuracy = 0.01), NA
#'       )
#'     }
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # add labels based on an external data
#'   pvalues <- matrix(runif(60, 0, 0.5), nrow = 6, ncol = 10)
#'   Heatmap(matrix_data,
#'     rows_data = rows_data,
#'     rows_split_by = "group", cell_type = "label",
#'     label = function(x, i, j) {
#'       pv <- ComplexHeatmap::pindex(pvalues, i, j)
#'       ifelse(pv < 0.01, "***",
#'         ifelse(pv < 0.05, "**",
#'           ifelse(pv < 0.1, "*", NA)
#'         )
#'       )
#'     }
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # quickly simulate a GO board
#'   go <- matrix(sample(c(0, 1, NA), 81, replace = TRUE), ncol = 9)
#'
#'   Heatmap(
#'     go,
#'     # Do not cluster rows and columns and hide the annotations
#'     cluster_rows = FALSE, cluster_columns = FALSE,
#'     row_name_annotation = FALSE, column_name_annotation = FALSE,
#'     show_row_names = FALSE, show_column_names = FALSE,
#'     # Set the legend items
#'     values_by = "Players", legend_discrete = TRUE,
#'     legend_items = c("Player 1" = 0, "Player 2" = 1),
#'     # Set the pawns
#'     cell_type = "dot", dot_size = function(x) ifelse(is.na(x), 0, 10),
#'     dot_size_name = NULL, # hide the dot size legend
#'     palcolor = c("white", "black"),
#'     # Set the board
#'     add_reticle = TRUE,
#'     # Set the size of the board
#'     width = ggplot2::unit(105, "mm"), height = ggplot2::unit(105, "mm")
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # Make the row/column name annotation thicker
#'   Heatmap(matrix_data,
#'     # Use the "name." prefix
#'     column_annotation_params = list(name.columns = list(height = 5)),
#'     row_annotation_params = list(name.rows = list(width = 5))
#'   )
#' }
#'
#' # Use long form data
#' N <- 500
#' data <- data.frame(
#'   value = rnorm(N),
#'   c = sample(letters[1:8], N, replace = TRUE),
#'   r = sample(LETTERS[1:5], N, replace = TRUE),
#'   p = sample(c("x", "y"), N, replace = TRUE),
#'   q = sample(c("X", "Y", "Z"), N, replace = TRUE),
#'   a = as.character(sample(1:5, N, replace = TRUE)),
#'   p1 = runif(N),
#'   p2 = runif(N)
#' )
#'
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     rows_by = "r", columns_by = "c", values_by = "value",
#'     rows_split_by = "p", columns_split_by = "q", show_column_names = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # split into multiple heatmaps
#'   Heatmap(data,
#'     values_by = "value", columns_by = "c", rows_by = "r", split_by = "p",
#'     upper_cutoff = 2, lower_cutoff = -2, legend.position = c("none", "right"),
#'     design = "AAAAAA#BBBBBBB"
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # cell_type = "bars" (default is "tile")
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "bars"
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   p <- Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "dot", dot_size = length, dot_size_name = "data points",
#'     add_bg = TRUE, add_reticle = TRUE
#'   )
#'   p
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   dot_size_data <- p@data
#'   # Make it big so we can see if we get the right indexing
#'   # for dot_size function
#'   dot_size_data["A", "a"] <- max(dot_size_data) * 2
#'
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "dot", dot_size_name = "data points",
#'     dot_size = function(x, i, j) ComplexHeatmap::pindex(dot_size_data, i, j),
#'     show_row_names = TRUE, show_column_names = TRUE,
#'     add_bg = TRUE, add_reticle = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "pie", pie_group_by = "q", pie_size = sqrt,
#'     add_bg = TRUE, add_reticle = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "violin", add_bg = TRUE, add_reticle = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     cell_type = "boxplot", add_bg = TRUE, add_reticle = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     column_annotation = list(r1 = "p", r2 = "q", r3 = "p1"),
#'     column_annotation_type = list(r1 = "ring", r2 = "bar", r3 = "violin"),
#'     column_annotation_params = list(
#'       r1 = list(height = grid::unit(10, "mm"), show_legend = FALSE),
#'       r3 = list(height = grid::unit(18, "mm"))
#'     ),
#'     row_annotation = c("q", "p2", "a"),
#'     row_annotation_side = "right",
#'     row_annotation_type = list(q = "pie", p2 = "density", a = "simple"),
#'     row_annotation_params = list(q = list(width = grid::unit(12, "mm"))),
#'     show_row_names = TRUE, show_column_names = TRUE
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   Heatmap(data,
#'     values_by = "value", rows_by = "r", columns_by = "c",
#'     split_by = "p", palette = list(x = "Reds", y = "Blues")
#'   )
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # implies in_form = "wide-rows"
#'   Heatmap(data, rows_by = c("p1", "p2"), columns_by = "c")
#' }
#' if (requireNamespace("cluster", quietly = TRUE)) {
#'   # implies wide-columns
#'   Heatmap(data, rows_by = "r", columns_by = c("p1", "p2"))
#' }
#' }
Heatmap <- function(
    data, values_by = NULL, values_fill = NA, name = NULL,
    # data definition
    in_form = c("auto", "matrix", "wide-columns", "wide-rows", "long"),
    split_by = NULL, split_by_sep = "_",
    rows_by = NULL, rows_by_sep = "_", rows_split_by = NULL, rows_split_by_sep = "_",
    columns_by = NULL, columns_by_sep = "_", columns_split_by = NULL, columns_split_by_sep = "_",
    rows_data = NULL, columns_data = NULL,
    # names
    columns_name = NULL, columns_split_name = NULL,
    rows_name = NULL, rows_split_name = NULL,
    # palettes
    palette = "RdBu", palcolor = NULL,
    rows_palette = "Paired", rows_palcolor = NULL, rows_split_palette = "simspec", rows_split_palcolor = NULL,
    columns_palette = "Paired", columns_palcolor = NULL, columns_split_palette = "simspec", columns_split_palcolor = NULL,
    # cell_type: pies
    pie_size_name = "size", pie_size = NULL, pie_values = "length", pie_name = NULL,
    pie_group_by = NULL, pie_group_by_sep = "_", pie_palette = "Spectral", pie_palcolor = NULL,
    # cell_type: bars
    bars_sample = 100,
    # cell_type: label
    label = identity, label_size = 10,
    # cell_type: violin
    violin_fill = NULL,
    # cell_type: boxplot
    boxplot_fill = NULL,
    # cell_type: dot
    dot_size = 8, dot_size_name = "size",
    # legend
    legend_items = NULL, legend_discrete = FALSE,
    legend.position = "right", legend.direction = "vertical",
    # values
    scale = c("none", "row", "column"),
    lower_quantile = 0, upper_quantile = 0.99, lower_cutoff = NULL, upper_cutoff = NULL,
    # bg
    add_bg = FALSE, bg_alpha = 0.5,
    # reticle
    add_reticle = FALSE, reticle_color = "grey",
    # cell border (for tile type)
    cell_border = FALSE, cell_border_color = "grey80", cell_border_width = 0.1,
    # passed to ComplexHeatmap::Heatmap
    column_name_annotation = TRUE, column_name_legend = NULL,
    row_name_annotation = TRUE, row_name_legend = NULL,
    cluster_columns = TRUE, cluster_rows = TRUE, show_row_names = !row_name_annotation, show_column_names = !column_name_annotation,
    border = TRUE, title = NULL, column_title = character(0), row_title = character(0), na_col = "grey85",
    row_names_side = "right", column_names_side = "bottom",
    column_annotation = NULL, column_annotation_side = "top", column_annotation_palette = "Paired", column_annotation_palcolor = NULL,
    column_annotation_type = "auto", column_annotation_params = list(), column_annotation_agg = NULL,
    row_annotation = NULL, row_annotation_side = "left", row_annotation_palette = "Paired", row_annotation_palcolor = NULL,
    row_annotation_type = "auto", row_annotation_params = list(), row_annotation_agg = NULL,
    # misc
    flip = FALSE, alpha = 1, seed = 8525,
    # cell customization
    layer_fun_callback = NULL, cell_type = c("tile", "bars", "label", "dot", "violin", "boxplot", "pie"), cell_agg = NULL,
    # subplots
    combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, axes = NULL, axis_titles = axes, guides = NULL, design = NULL,
    ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for Heatmap(). Install from Bioconductor:\n",
      "  BiocManager::install('ComplexHeatmap')", call. = FALSE)
  }

  validate_common_args(
    seed = seed,
    alpha = alpha,
    legend.position = legend.position,
    legend.direction = legend.direction
  )
  in_form <- match.arg(in_form)
  cell_type <- match.arg(cell_type)

  hmdata <- process_heatmap_data(
    data,
    in_form = in_form, values_by = values_by, name = name,
    split_by = split_by, split_by_sep = split_by_sep,
    rows_by = rows_by, rows_by_sep = rows_by_sep, rows_name = rows_name,
    rows_split_by = rows_split_by, rows_split_by_sep = rows_split_by_sep, rows_split_name = rows_split_name,
    columns_by = columns_by, columns_by_sep = columns_by_sep, columns_name = columns_name,
    columns_split_by = columns_split_by, columns_split_by_sep = columns_split_by_sep, columns_split_name = columns_split_name,
    pie_group_by = pie_group_by, pie_group_by_sep = pie_group_by_sep, pie_name = pie_name,
    rows_data = rows_data, columns_data = columns_data
  )

  palette <- check_palette(palette, names(hmdata$data))
  palcolor <- check_palcolor(palcolor, names(hmdata$data))
  rows_palette <- check_palette(rows_palette, names(hmdata$data))
  rows_palcolor <- check_palcolor(rows_palcolor, names(hmdata$data))
  rows_split_palette <- check_palette(rows_split_palette, names(hmdata$data))
  rows_split_palcolor <- check_palcolor(rows_split_palcolor, names(hmdata$data))
  columns_palette <- check_palette(columns_palette, names(hmdata$data))
  columns_palcolor <- check_palcolor(columns_palcolor, names(hmdata$data))
  columns_split_palette <- check_palette(columns_split_palette, names(hmdata$data))
  columns_split_palcolor <- check_palcolor(columns_split_palcolor, names(hmdata$data))
  pie_palette <- check_palette(pie_palette, names(hmdata$data))
  pie_palcolor <- check_palcolor(pie_palcolor, names(hmdata$data))
  legend.direction <- check_legend_param(legend.direction, names(hmdata$data), "legend.direction")
  legend.position <- check_legend_param(legend.position, names(hmdata$data), "legend.position")

  ggplot2_v4 <- utils::compareVersion(as.character(utils::packageVersion("ggplot2")), "4") >= 0
  return_grob <- !ggplot2_v4 || length(hmdata$data) > 1

  plots <- lapply(
    names(hmdata$data), function(nm) {
      default_title <- if (length(hmdata$data) == 1 && identical(nm, "...")) NULL else nm
      if (is.function(title)) {
        title <- title(default_title)
      } else {
        title <- title %||% default_title
      }

      HeatmapAtomic(
        data = hmdata$data[[nm]], values_by = hmdata$values_by, values_fill = values_fill,
        rows_by = hmdata$rows_by, rows_split_by = hmdata$rows_split_by,
        columns_by = hmdata$columns_by, columns_split_by = hmdata$columns_split_by,
        palette = palette[[nm]], palcolor = palcolor[[nm]],
        rows_palette = rows_palette[[nm]], rows_palcolor = rows_palcolor[[nm]],
        rows_split_palette = rows_split_palette[[nm]], rows_split_palcolor = rows_split_palcolor[[nm]],
        columns_palette = columns_palette[[nm]], columns_palcolor = columns_palcolor[[nm]],
        columns_split_palette = columns_split_palette[[nm]], columns_split_palcolor = columns_split_palcolor[[nm]],
        pie_size_name = pie_size_name, pie_size = pie_size, pie_values = pie_values,
        pie_group_by = hmdata$pie_group_by, pie_palette = pie_palette[[nm]], pie_palcolor = pie_palcolor[[nm]],
        bars_sample = bars_sample,
        label = label, label_size = label_size,
        violin_fill = violin_fill,
        boxplot_fill = boxplot_fill,
        dot_size = dot_size, dot_size_name = dot_size_name,
        legend_items = legend_items, legend_discrete = legend_discrete,
        legend.position = legend.position[[nm]], legend.direction = legend.direction[[nm]],
        scale = scale,
        lower_quantile = lower_quantile, upper_quantile = upper_quantile,
        lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff,
        add_bg = add_bg, bg_alpha = bg_alpha,
        add_reticle = add_reticle, reticle_color = reticle_color,
        cell_border = cell_border, cell_border_color = cell_border_color, cell_border_width = cell_border_width,
        column_name_annotation = column_name_annotation, column_name_legend = column_name_legend,
        row_name_annotation = row_name_annotation, row_name_legend = row_name_legend,
        cluster_columns = cluster_columns, cluster_rows = cluster_rows, show_row_names = show_row_names, show_column_names = show_column_names,
        border = border, title = title, column_title = column_title, row_title = row_title, na_col = na_col,
        row_names_side = row_names_side, column_names_side = column_names_side,
        column_annotation = column_annotation, column_annotation_side = column_annotation_side,
        column_annotation_palette = column_annotation_palette, column_annotation_palcolor = column_annotation_palcolor,
        column_annotation_type = column_annotation_type, column_annotation_params = column_annotation_params,
        column_annotation_agg = column_annotation_agg,
        row_annotation = row_annotation, row_annotation_side = row_annotation_side,
        row_annotation_palette = row_annotation_palette, row_annotation_palcolor = row_annotation_palcolor,
        row_annotation_type = row_annotation_type, row_annotation_params = row_annotation_params,
        row_annotation_agg = row_annotation_agg,
        flip = flip, alpha = alpha, seed = seed, return_grob = return_grob,
        layer_fun_callback = layer_fun_callback, cell_type = cell_type, cell_agg = cell_agg,
        ...
      )
    }
  )

  p <- combine_plots(plots,
    combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,
    axes = axes, axis_titles = axis_titles, guides = guides, design = design
  )
  if (length(plots) == 1) {
    attr(p, "data") <- attr(plots[[1]], "data")
  }

  # Return the plot object
  # When return_grob = FALSE, p is a HeatmapList object with auto-printing behavior
  # The initial draw() in HeatmapAtomic was captured to a null device to prevent
  # printing during assignment, so this returned object can print normally when
  # called directly (e.g., in a Jupyter cell)
  p
}

# --- Source: ggforge/R/kmplot.R ---

#' Kaplan-Meier Survival Plot Atomic
#'
#' @description
#' Creates a single Kaplan-Meier survival plot without splitting.
#' This is the core plotting function that handles the actual ggplot construction.
#'
#' @inheritParams parameters
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_step geom_ribbon scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs annotate coord_cartesian geom_hline geom_vline scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous expansion geom_text geom_segment
#' @importFrom dplyr group_by summarise ungroup arrange filter mutate .data
#' @importFrom rlang sym syms := "%||%"
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom survival survfit Surv survdiff
#' @importFrom stats pchisq
KMPlotAtomic <- function(
    data,
    time,
    status,
    group_by = NULL,
    group_by_sep = "_",
    group_name = NULL,
    show_pval = TRUE,
    pval_method = "logrank",
    pval_digits = 4,
    pval_size = 4,
    pval_coord = c(0.05, 0.1),
    show_conf_int = FALSE,
    conf_alpha = 0.2,
    show_median_line = "none",
    median_linetype = 2,
    median_linewidth = 0.6,
    line_width = 1.3,
    show_risk_table = FALSE,
    risk_table_height = 0.25,
    risk_table_fontsize = 3.5,
    show_censors = TRUE,
    censor_shape = 3,
    censor_size = 4,
    censor_stroke = 0.5,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    aspect.ratio = NULL,
    x_breaks = NULL,
    y_breaks = waiver(),
    x_min = NULL,
    x_max = NULL,
    y_min = 0,
    y_max = 1,
    legend.position = "top",
    legend.direction = "horizontal",
    title = NULL,
    subtitle = NULL,
    xlab = "Time",
    ylab = "Survival Probability",
    facet_by = NULL,
    facet_scales = "fixed",
    facet_nrow = NULL,
    facet_ncol = NULL,
    facet_byrow = TRUE,
    ...) {
  # Get ggplot function (with gglogger support)
  ggplot <- get_ggplot()

  # Validate columns
  time <- validate_columns(data, time)
  status <- validate_columns(data, status)
  group_by <- validate_columns(
    data, group_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = group_by_sep
  )

  # Handle NULL group_by (single curve)
  has_groups <- !is.null(group_by)
  if (!has_groups) {
    group_by <- ".group"
    data[[group_by]] <- factor("All")
  }

  # Build survival fit
  if (has_groups && length(unique(data[[group_by]])) > 1) {
    formula_str <- sprintf("Surv(%s, %s) ~ %s", time, status, group_by)
  } else {
    formula_str <- sprintf("Surv(%s, %s) ~ 1", time, status)
  }
  sfit <- survival::survfit(as.formula(formula_str), data = data)

  # Extract survival data
  if (has_groups && length(unique(data[[group_by]])) > 1) {
    surv_data <- data.frame(
      time = sfit$time,
      n.risk = sfit$n.risk,
      n.event = sfit$n.event,
      n.censor = sfit$n.censor,
      surv = sfit$surv,
      std.err = sfit$std.err,
      lower = sfit$lower,
      upper = sfit$upper,
      strata = rep(names(sfit$strata), sfit$strata),
      stringsAsFactors = FALSE
    )
    surv_data[[group_by]] <- sub(".*=", "", surv_data$strata)
    surv_data[[group_by]] <- factor(surv_data[[group_by]], levels = sub(".*=", "", names(sfit$strata)))
  } else {
    surv_data <- data.frame(
      time = sfit$time,
      n.risk = sfit$n.risk,
      n.event = sfit$n.event,
      n.censor = sfit$n.censor,
      surv = sfit$surv,
      std.err = sfit$std.err,
      lower = sfit$lower,
      upper = sfit$upper,
      stringsAsFactors = FALSE
    )
    surv_data[[group_by]] <- factor("All")
  }

  # Get group levels
  group_levels <- levels(surv_data[[group_by]])
  n_groups <- length(group_levels)

  # Add time=0 point for step plot
  time0_list <- lapply(group_levels, function(grp) {
    data.frame(
      time = 0,
      n.risk = max(surv_data$n.risk[surv_data[[group_by]] == grp], na.rm = TRUE),
      n.event = 0,
      n.censor = 0,
      surv = 1,
      std.err = 0,
      lower = 1,
      upper = 1,
      stringsAsFactors = FALSE
    ) -> df
    df[[group_by]] <- factor(grp, levels = group_levels)
    if (has_groups && length(unique(data[[group_by]])) > 1) {
      df$strata <- names(sfit$strata)[which(group_levels == grp)]
    }
    df
  })
  time0 <- do.call(rbind, time0_list)

  surv_data <- rbind(time0, surv_data)
  surv_data <- surv_data %>% dplyr::arrange(!!sym(group_by), .data$time)

  # Setup colors
  colors <- get_palette(group_levels, palette = palette, palcolor = palcolor)

  # Calculate axis limits
  if (is.null(x_min)) x_min <- 0
  if (is.null(x_max)) x_max <- max(surv_data$time, na.rm = TRUE) * 1.05

  # Calculate x-axis breaks
  if (is.null(x_breaks)) {
    x_breaks <- pretty(c(x_min, x_max), n = 6)
    x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  }

  # Calculate p-value if requested
  pval_text <- NULL
  if (show_pval && has_groups && n_groups > 1) {
    tryCatch(
      {
        sdiff <- survival::survdiff(as.formula(formula_str), data = data)
        pval <- 1 - stats::pchisq(sdiff$chisq, length(sdiff$n) - 1)

        # Smart p-value formatting
        if (pval >= 1e-4) {
          # Regular format with 4 significant figures for p >= 0.0001
          pval_formatted <- signif(pval, 4)
          pval_text <- sprintf("P = %g", pval_formatted)
        } else {
          # Scientific notation with 2 decimal places for p < 0.0001
          pval_text <- sprintf("P = %.2e", pval)
        }
      },
      error = function(e) {
        warning("Failed to calculate p-value: ", e$message, call. = FALSE)
      }
    )
  }

  # Build main survival plot
  p <- ggplot(surv_data, aes(x = .data$time, y = .data$surv, color = !!sym(group_by)))

  # Add confidence interval ribbons
  if (show_conf_int) {
    p <- p + ggplot2::geom_ribbon(
      aes(ymin = .data$lower, ymax = .data$upper, fill = !!sym(group_by)),
      alpha = conf_alpha,
      linetype = 0,
      show.legend = FALSE
    )
  }

  # Add survival curves
  p <- p + ggplot2::geom_step(linewidth = line_width, alpha = alpha)

  # Add censoring marks
  if (show_censors) {
    censor_data <- surv_data[surv_data$n.censor > 0, ]
    if (nrow(censor_data) > 0) {
      p <- p + ggplot2::geom_point(
        data = censor_data,
        aes(x = .data$time, y = .data$surv, color = !!sym(group_by)),
        shape = censor_shape,
        size = censor_size,
        stroke = censor_stroke,
        show.legend = FALSE
      )
    }
  }

  # Add median survival lines
  if (show_median_line %in% c("h", "v", "hv")) {
    # Calculate median survival times
    median_surv <- summary(sfit)$table
    if (is.matrix(median_surv)) {
      median_times <- median_surv[, "median"]
      names(median_times) <- sub(".*=", "", rownames(median_surv))
    } else {
      median_times <- median_surv["median"]
      names(median_times) <- "All"
    }

    if (show_median_line %in% c("h", "hv")) {
      p <- p + ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = median_linetype,
        linewidth = median_linewidth,
        color = "grey50"
      )
    }

    if (show_median_line %in% c("v", "hv")) {
      for (i in seq_along(median_times)) {
        if (!is.na(median_times[i])) {
          p <- p + ggplot2::geom_vline(
            xintercept = median_times[i],
            linetype = median_linetype,
            linewidth = median_linewidth,
            color = colors[i]
          )
        }
      }
    }
  }

  # Add p-value annotation with background (mimicking corplot.R style)
  if (!is.null(pval_text)) {
    # Combine method and p-value with newline
    label_text <- paste0("Log-rank\n", pval_text)

    # Create annotation data frame
    pval_df <- data.frame(
      x = x_max * pval_coord[1],
      y = pval_coord[2],
      label = label_text,
      stringsAsFactors = FALSE
    )

    # Add combined label with background using ggrepel (exactly like corplot.R)
    p <- p + ggrepel::geom_text_repel(
      data = pval_df,
      aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0,
      vjust = 0,
      size = pval_size,
      color = "black",
      bg.color = "white",
      bg.r = 0.1,
      min.segment.length = Inf,
      max.overlaps = 100,
      force = 0,
      box.padding = 0,
      point.padding = 0,
      segment.color = "transparent",
      inherit.aes = FALSE
    )
  }

  # Setup legend labels
  if (has_groups) {
    legend_labels <- paste0(group_levels, " (n=", summary(sfit)$n, ")")
    names(legend_labels) <- group_levels
  } else {
    legend_labels <- paste0("All (n=", summary(sfit)$n, ")")
    names(legend_labels) <- "All"
  }

  # Add color and fill scales
  p <- p +
    ggplot2::scale_color_manual(
      name = group_name %||% group_by,
      values = colors,
      labels = legend_labels
    ) +
    ggplot2::scale_fill_manual(
      name = group_name %||% group_by,
      values = colors,
      labels = legend_labels,
      guide = "none"
    )

  # Add scales
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = c(x_min, x_max),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      limits = c(y_min, y_max),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    )

  # Add labels
  p <- p + ggplot2::labs(
    title = title,
    subtitle = subtitle,
    x = xlab,
    y = ylab
  )

  # Apply theme
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")
  text_size_scale <- base_size / 12 # Define once for both main plot and risk table
  p <- p + do.call(theme, theme_args)

  # Apply data-driven styling
  p <- apply_style_theme(
    plot = p,
    data = surv_data,
    x_var = time,
    y_var = "surv",
    flip = FALSE,
    base_size = base_size,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Custom theme adjustments
  if (!is.null(aspect.ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect.ratio)
  }

  # Add margin for risk table
  if (show_risk_table) {
    p <- p + ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5)
    )
  }

  # Build risk table if requested
  if (show_risk_table) {
    # Create risk table data
    risk_data <- data.frame()
    for (grp in group_levels) {
      grp_data <- surv_data[surv_data[[group_by]] == grp, ]
      n_at_risk <- sapply(x_breaks, function(t) {
        idx <- which(grp_data$time <= t)
        if (length(idx) > 0) {
          grp_data$n.risk[max(idx)]
        } else {
          0
        }
      })
      risk_data <- rbind(risk_data, data.frame(
        time = x_breaks,
        n.risk = n_at_risk,
        group = grp,
        stringsAsFactors = FALSE
      ))
    }
    names(risk_data)[3] <- group_by
    risk_data[[group_by]] <- factor(risk_data[[group_by]], levels = rev(group_levels))

    # Build risk table plot
    p_risk <- ggplot(risk_data, aes(x = .data$time, y = !!sym(group_by))) +
      ggplot2::geom_text(
        aes(label = .data$n.risk),
        size = risk_table_fontsize * text_size_scale,
        vjust = 0.5,
        color = "black"
      ) +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        limits = c(x_min, x_max),
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      ) +
      ggplot2::labs(
        x = xlab,
        y = "At risk",
        title = "Number at risk"
      ) +
      do.call(theme, theme_args)

    # Apply same styling as main plot
    p_risk <- apply_style_theme(
      plot = p_risk,
      data = risk_data,
      x_var = "time",
      y_var = group_by,
      flip = FALSE,
      base_size = base_size,
      legend.position = "none"
    )

    # Additional risk table specific styling
    p_risk <- p_risk + ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        color = "black",
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        face = "bold",
        hjust = 0
      ),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )

    # Remove x-axis title from main plot
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())

    # Combine plots
    p_combined <- p / p_risk +
      patchwork::plot_layout(
        ncol = 1,
        heights = c(1 - risk_table_height, risk_table_height)
      )

    return(p_combined)
  }

  # Add faceting if requested
  if (!is.null(facet_by)) {
    p <- add_facets(
      p, facet_by, facet_scales,
      facet_nrow, facet_ncol, facet_byrow
    )
  }

  return(p)
}

#' Kaplan-Meier Survival Plot
#'
#' @description
#' Create publication-ready Kaplan-Meier survival curves with optional risk tables,
#' confidence intervals, and statistical comparisons.
#'
#' This function provides a complete implementation of Kaplan-Meier survival analysis
#' visualization, supporting single or multiple groups, with automatic p-value
#' calculation using the log-rank test.
#'
#' @inheritParams parameters
#' @param time Column name for time variable (numeric).
#' @param status Column name for event status (1=event, 0=censored).
#' @param group_by Column(s) for grouping survival curves.
#' @param group_by_sep Separator for concatenating multiple group columns.
#' @param group_name Legend title for groups.
#' @param show_pval Show log-rank test p-value.
#' @param pval_method P-value calculation method ("logrank").
#' @param pval_digits Number of digits for p-value.
#' @param pval_size Text size for p-value.
#' @param pval_coord Position of p-value as c(x, y) where x is fraction of x-axis, y is absolute.
#' @param show_conf_int Show confidence interval ribbons.
#' @param conf_alpha Transparency for confidence interval ribbons.
#' @param show_median_line Show median survival lines: "none", "h", "v", "hv".
#' @param median_linetype Line type for median survival lines.
#' @param median_linewidth Line width for median survival lines.
#' @param line_width Width of survival curves.
#' @param show_risk_table Show risk table below plot.
#' @param risk_table_height Relative height of risk table (0-1).
#' @param risk_table_fontsize Font size for numbers in risk table.
#' @param show_censors Show censoring marks on curves.
#' @param censor_shape Shape for censoring marks.
#' @param censor_size Size for censoring marks.
#' @param censor_stroke Stroke width for censoring marks.
#' @param x_breaks Custom x-axis breaks (time points).
#' @param y_breaks Custom y-axis breaks.
#' @param x_min Minimum x-axis value.
#' @param x_max Maximum x-axis value.
#' @param y_min Minimum y-axis value (default: 0).
#' @param y_max Maximum y-axis value (default: 1).
#' @param ... Additional arguments passed to atomic plotting functions.
#'
#' @return A ggplot object, list of plots, or combined plots
# @keywords internal
#' @importFrom ggplot2 waiver
#' @examples
#' \donttest{
#' library(survival)
#'
#' # Basic Kaplan-Meier plot
#' KMPlot(data = lung, time = "time", status = "status")
#'
#' # Multiple groups with p-value
#' KMPlot(
#'   data = lung,
#'   time = "time",
#'   status = "status",
#'   group_by = "sex",
#'   show_pval = TRUE
#' )
#'
#' # With risk table
#' KMPlot(
#'   data = lung,
#'   time = "time",
#'   status = "status",
#'   group_by = "sex",
#'   show_risk_table = TRUE,
#'   show_pval = TRUE
#' )
#'
#' # With confidence intervals and median lines
#' KMPlot(
#'   data = lung,
#'   time = "time",
#'   status = "status",
#'   group_by = "sex",
#'   show_conf_int = TRUE,
#'   show_median_line = "hv",
#'   palette = "Set1"
#' )
#'
#' # Publication-ready plot
#' KMPlot(
#'   data = lung,
#'   time = "time",
#'   status = "status",
#'   group_by = "sex",
#'   show_risk_table = TRUE,
#'   show_pval = TRUE,
#'   show_conf_int = TRUE,
#'   show_median_line = "hv",
#'   palette = "jco",
#'   title = "Overall Survival by Sex",
#'   xlab = "Time (days)",
#'   ylab = "Survival Probability",
#'   theme_args = list(base_size = 14)
#' )
#' }
KMPlot <- function(
    data,
    time,
    status,
    group_by = NULL,
    group_by_sep = "_",
    group_name = NULL,
    split_by = NULL,
    split_by_sep = "_",
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    show_pval = TRUE,
    pval_method = "logrank",
    pval_digits = 4,
    pval_size = 4,
    pval_coord = c(0.05, 0.1),
    show_conf_int = FALSE,
    conf_alpha = 0.2,
    show_median_line = "none",
    median_linetype = 2,
    median_linewidth = 0.6,
    line_width = 1.3,
    show_risk_table = FALSE,
    risk_table_height = 0.25,
    risk_table_fontsize = 3.5,
    show_censors = TRUE,
    censor_shape = 3,
    censor_size = 4,
    censor_stroke = 0.5,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    aspect.ratio = NULL,
    x_breaks = NULL,
    y_breaks = waiver(),
    x_min = NULL,
    x_max = NULL,
    y_min = 0,
    y_max = 1,
    legend.position = "top",
    legend.direction = "horizontal",
    title = NULL,
    subtitle = NULL,
    xlab = "Time",
    ylab = "Survival Probability",
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    seed = 8525,
    axes = NULL,
    axis_titles = axes,
    guides = NULL,
    design = NULL,
    ...) {
  # Validate median line option
  show_median_line <- match.arg(show_median_line, c("none", "h", "v", "hv"))

  # Validate common arguments
  validate_common_args(
    seed = seed,
    facet_by = facet_by,
    split_by = split_by,
    theme = theme,
    palette = palette,
    alpha = alpha,
    aspect.ratio = aspect.ratio,
    legend.position = if (is.character(legend.position)) legend.position else "top",
    legend.direction = legend.direction
  )

  # Process theme
  theme <- process_theme(theme)

  # Validate split_by column
  split_by <- validate_columns(
    data, split_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = split_by_sep
  )

  # Collect all parameters for passing to atomic function
  params <- as.list(environment())
  params$data <- NULL # Remove data from params

  # Build plot using standard workflow
  build_plot(
    data = data,
    atomic_fn = KMPlotAtomic,
    params = params,
    split_by = split_by,
    facet_by = facet_by,
    combine = combine
  )
}

# --- Source: ggforge/R/coxplot.R ---

#' Cox Proportional-Hazards Model Visualization
#'
#' @description
#' Create publication-ready Cox proportional-hazards model visualizations including
#' hazard ratio curves, forest plots (simple and detailed), with support for
#' parallel computation and multiple visualization types.
#'
#' @inheritParams parameters
#' @param data A data frame containing survival time, event status, and variables for Cox analysis.
#' @param time Column name for time variable (numeric).
#' @param event Column name for event status (1=event, 0=censored).
#' @param vars Column name(s) for variables to include in Cox regression. If NULL, uses `var`.
#' @param var Column name for single variable (used for curve plot or when vars is NULL).
#' @param plot_type Type of Cox plot: "curve" (hazard ratio curve), "forest" (simple forest plot),
#'   or "forest2" (detailed forest plot with HR and p-value columns).
#' @param scale Logical. Whether to standardize variables to z-scores before Cox regression.
#' @param nonExpression_ratio Numeric. Threshold ratio for filtering non-expressed genes.
#'   If the proportion of zero values exceeds this ratio, the variable will be filtered out.
#' @param parallel Logical. Whether to perform parallel computation for multiple variables.
#' @param n_cores Integer. Number of cores to use for parallel computation.
#'   Default is detectCores() - 6, minimum 1.
#' @param ribbon_color Color for confidence interval ribbon (curve plot).
#' @param ribbon_alpha Alpha transparency for ribbon (curve plot).
#' @param line_color Color for hazard ratio line (curve plot).
#' @param line_type Line type for hazard ratio curve.
#' @param line_width Line width for hazard ratio curve.
#' @param show_cindex Logical. Whether to show concordance index (Cindex) on plot.
#' @param text_size Size of annotation text.
#' @param text_digit Number of significant digits for text display.
#' @param text_face Font face for annotation text.
#' @param point_colors Vector of colors for points in forest plot.
#'   Default: c("#ED6355", "#118ab2", "grey") for Risky, Protective, NoSig.
#' @param point_size Size of points in forest plot.
#' @param point_border_width Border width of points.
#' @param point_border_size Size of point borders.
#' @param line_colors Colors for error bars in forest plot.
#' @param cutoff_vline_type Line type for cutoff vertical line (HR = 1).
#' @param cutoff_vline_width Line width for cutoff vertical line.
#' @param cutoff_vline_color Color for cutoff vertical line.
#' @param x_log10_scale Logical. Whether to use log10 scale for x-axis in forest plot.
#' @param text_colors Text colors in forest2 plot.
#' @param digits Number of significant digits for numeric display.
#' @param rel_width Relative widths of columns in forest2 plot.
#'   Default: c(0.8, 1.4, 1.2, 0.6) for Variable, HR plot, HR text, P-value.
#' @param ... Additional arguments passed to atomic plotting functions.
#'
#' @return A ggplot object, patchwork object, or list of plots depending on settings
# @keywords internal
#' @importFrom survival coxph Surv
#' @importFrom stats predict as.formula quantile complete.cases
#' @importFrom parallel detectCores
#' @importFrom ggplot2 waiver
#'
#' @examples
#' \donttest{
#' # Prepare example data
#' library(survival)
#' data(lung)
#' lung$status <- lung$status - 1 # Convert to 0/1
#'
#' # Single variable hazard ratio curve
#' CoxPlot(
#'   data = lung,
#'   time = "time",
#'   event = "status",
#'   var = "age",
#'   plot_type = "curve"
#' )
#'
#' # Multiple variables forest plot (simple)
#' CoxPlot(
#'   data = lung,
#'   time = "time",
#'   event = "status",
#'   vars = c("age", "ph.ecog", "ph.karno", "pat.karno"),
#'   plot_type = "forest"
#' )
#'
#' # Multiple variables forest plot (detailed)
#' CoxPlot(
#'   data = lung,
#'   time = "time",
#'   event = "status",
#'   vars = c("age", "ph.ecog", "ph.karno", "pat.karno"),
#'   plot_type = "forest2"
#' )
#' }
CoxPlot <- function(
    data,
    time = "time",
    event = "event",
    vars = NULL,
    var = NULL,
    plot_type = c("curve", "forest", "forest2"),
    scale = FALSE,
    nonExpression_ratio = 1,
    parallel = FALSE,
    n_cores = max(1, parallel::detectCores() - 6),
    # Curve plot parameters
    ribbon_color = "#EFA63A",
    ribbon_alpha = 0.6,
    line_color = "#3a6ea5",
    line_type = 1,
    line_width = 1,
    show_cindex = TRUE,
    text_size = 3.5,
    text_digit = 3,
    text_face = "bold.italic",
    # Forest plot parameters
    point_colors = c("#ED6355", "#118ab2", "grey"),
    point_size = 2,
    point_border_width = 0.5,
    point_border_size = 3,
    line_colors = point_colors,
    cutoff_vline_type = 2,
    cutoff_vline_width = 0.5,
    cutoff_vline_color = "grey30",
    x_log10_scale = FALSE,
    # Forest2 plot parameters
    text_colors = "black",
    digits = 3,
    rel_width = c(0.8, 1.4, 1.2, 0.6),
    # Standard ggforge parameters
    split_by = NULL,
    split_by_sep = "_",
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "bottom",
    legend.direction = "horizontal",
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    seed = 8525,
    axes = NULL,
    axis_titles = axes,
    guides = NULL,
    design = NULL,
    ...) {
  # Validate plot type
  plot_type <- match.arg(plot_type)

  # Handle var/vars logic
  if (is.null(vars) && is.null(var)) {
    stop("Either 'vars' or 'var' must be specified", call. = FALSE)
  }
  if (is.null(vars)) {
    vars <- var
  }

  # Validate plot type and variables
  if (plot_type == "curve" && length(vars) > 1) {
    warning("Curve plot only supports single variable. Using first variable: ", vars[1], call. = FALSE)
    vars <- vars[1]
  }

  # Validate parameters
  .validate_cox_params(
    n_cores = n_cores,
    digits = digits,
    rel_width = rel_width,
    nonExpression_ratio = nonExpression_ratio
  )

  # Validate common arguments
  validate_common_args(
    seed = seed,
    facet_by = facet_by,
    split_by = split_by,
    theme = theme,
    palette = palette,
    alpha = alpha,
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Process theme
  theme <- process_theme(theme)

  # Validate columns
  time <- validate_columns(data, time)
  event <- validate_columns(data, event)
  vars <- validate_columns(data, vars, allow_multi = TRUE)
  split_by <- validate_columns(
    data, split_by,
    force_factor = TRUE,
    allow_multi = TRUE,
    concat_multi = TRUE,
    concat_sep = split_by_sep
  )

  # Collect all parameters
  params <- as.list(environment())
  params$data <- NULL

  # Build plot using standard workflow
  build_plot(
    data = data,
    atomic_fn = CoxPlotAtomic,
    params = params,
    split_by = split_by,
    facet_by = facet_by,
    combine = combine
  )
}

#' Cox Plot Atomic Function
#'
#' @description
#' Creates a single Cox plot without splitting
#'
#' @inheritParams CoxPlot
#' @keywords internal
#' @importFrom ggplot2 ggplot aes
#' @importFrom dplyr %>%
#' @importFrom rlang "%||%"
CoxPlotAtomic <- function(
    data,
    time,
    event,
    vars,
    var,
    plot_type,
    scale,
    nonExpression_ratio,
    parallel,
    n_cores,
    ribbon_color,
    ribbon_alpha,
    line_color,
    line_type,
    line_width,
    show_cindex,
    text_size,
    text_digit,
    text_face,
    point_colors,
    point_size,
    point_border_width,
    point_border_size,
    line_colors,
    cutoff_vline_type,
    cutoff_vline_width,
    cutoff_vline_color,
    x_log10_scale,
    text_colors,
    digits,
    rel_width,
    facet_by = NULL,
    facet_scales = "fixed",
    facet_ncol = NULL,
    facet_nrow = NULL,
    facet_byrow = TRUE,
    theme = "theme_ggforge",
    theme_args = list(),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "bottom",
    legend.direction = "horizontal",
    seed = 8525,
    ...) {
  # Get ggplot function
  ggplot <- get_ggplot()

  # Route to appropriate plotting function
  if (plot_type == "curve") {
    p <- .cox_curve_plot(
      data = data,
      var = vars[1],
      time = time,
      event = event,
      ribbon_color = ribbon_color,
      ribbon_alpha = ribbon_alpha,
      line_color = line_color,
      line_type = line_type,
      line_width = line_width,
      text_size = text_size,
      text_digit = text_digit,
      text_face = text_face,
      xlab = xlab %||% vars[1],
      ylab = ylab %||% "Relative HR",
      title = title,
      show_cindex = show_cindex,
      theme = theme,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio
    )
  } else {
    # 统一计算 Cox 表（forest 和 forest2 共享）
    cox_table <- .get_cox_table(
      data = data,
      vars = vars,
      time = time,
      event = event,
      scale = scale,
      nonExpression_ratio = nonExpression_ratio,
      parallel = parallel,
      n_cores = n_cores
    )

    # 根据类型创建森林图
    if (plot_type == "forest") {
      p <- .cox_forest_plot(
        data = cox_table,
        point_colors = point_colors,
        point_size = point_size,
        point_border_width = point_border_width,
        point_border_size = point_border_size,
        line_colors = line_colors,
        line_type = line_type,
        line_width = line_width,
        cutoff_vline_type = cutoff_vline_type,
        cutoff_vline_width = cutoff_vline_width,
        cutoff_vline_color = cutoff_vline_color,
        x_log10_scale = x_log10_scale,
        xlab = xlab %||% "HR (95%CI)",
        title = title,
        subtitle = subtitle,
        legend.position = legend.position,
        legend.direction = legend.direction,
        theme = theme,
        theme_args = theme_args,
        aspect.ratio = aspect.ratio
      )
    } else { # forest2
      p <- .cox_forest2_plot(
        data = cox_table,
        point_colors = point_colors,
        point_size = point_size,
        point_border_width = point_border_width,
        point_border_size = point_border_size,
        line_colors = line_colors,
        line_type = line_type,
        line_width = line_width,
        cutoff_vline_type = cutoff_vline_type,
        cutoff_vline_width = cutoff_vline_width,
        cutoff_vline_color = cutoff_vline_color,
        x_log10_scale = x_log10_scale,
        text_colors = text_colors,
        digits = digits,
        rel_width = rel_width,
        theme = theme,
        theme_args = theme_args
      )
    }
  }

  # Add faceting if requested (forest2 doesn't support facet well)
  if (!is.null(facet_by) && plot_type != "forest2") {
    p <- add_facets(
      p, facet_by, facet_scales,
      facet_nrow, facet_ncol, facet_byrow
    )
  }

  return(p)
}

# =============================================================================
# Parameter Validation
# =============================================================================

#' Validate Cox-specific parameters
#' @keywords internal
.validate_cox_params <- function(n_cores, digits, rel_width, nonExpression_ratio) {
  # Validate n_cores
  if (n_cores < 1) {
    n_cores <- max(1, parallel::detectCores() - 6)
    message("n_cores was < 1, reset to ", n_cores)
  }

  # Validate digits
  if (digits < 1) {
    stop("'digits' must be at least 1", call. = FALSE)
  }

  # Validate rel_width
  if (length(rel_width) != 4) {
    stop("'rel_width' must have exactly 4 elements (Variable, HR plot, HR text, P-value)", call. = FALSE)
  }

  # Validate nonExpression_ratio
  if (nonExpression_ratio < 0 || nonExpression_ratio > 1) {
    stop("'nonExpression_ratio' must be between 0 and 1", call. = FALSE)
  }

  invisible(NULL)
}

# =============================================================================
# Cox Table Computation
# =============================================================================

#' Compute Cox Table
#'
#' @description
#' Internal function to compute Cox proportional-hazards model summary table
#' with optional parallel computation
#'
#' @keywords internal
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @importFrom dplyr rename
.get_cox_table <- function(
    data,
    vars,
    time = "time",
    event = "event",
    scale = FALSE,
    nonExpression_ratio = 1,
    parallel = FALSE,
    n_cores = 1) {
  # Rename columns using dplyr
  data_copy <- data %>%
    dplyr::rename(time = !!time, event = !!event)

  # Filter variables by expression ratio (only for numeric variables)
  vardata <- data_copy[, vars, drop = FALSE]

  # Identify numeric columns for filtering and scaling
  numeric_vars <- sapply(vardata, is.numeric)

  if (any(numeric_vars)) {
    # Filter numeric variables by expression ratio
    numeric_data <- vardata[, numeric_vars, drop = FALSE]
    ind_numeric <- colMeans(numeric_data == 0, na.rm = TRUE) <= nonExpression_ratio
    vars_keep_numeric <- names(numeric_data)[ind_numeric]
  } else {
    vars_keep_numeric <- character(0)
  }

  # Keep all non-numeric variables (factors, characters)
  vars_keep_nonnumeric <- names(vardata)[!numeric_vars]

  # Combined filtered variables
  vars <- c(vars_keep_numeric, vars_keep_nonnumeric)

  if (length(vars) == 0) {
    stop("No variables passed filtering threshold (nonExpression_ratio = ",
      nonExpression_ratio, ")",
      call. = FALSE
    )
  }

  # Scale only numeric variables if requested
  if (scale && length(vars_keep_numeric) > 0) {
    data_copy[, vars_keep_numeric] <- scale(data_copy[, vars_keep_numeric, drop = FALSE])
  }

  # Define function to fit Cox model
  fit_cox_model <- function(var_name) {
    tryCatch(
      {
        formula <- as.formula(paste("survival::Surv(time, event) ~", var_name))
        fit <- survival::coxph(formula, data = data_copy)
        fit_summary <- summary(fit)

        # Check for multiple coefficients (categorical variable with >2 levels)
        n_coef <- nrow(fit_summary$coefficients)

        if (n_coef > 1) {
          level_names <- rownames(fit_summary$coefficients)
          clean_names <- gsub(paste0("^", var_name), "", level_names)

          result <- data.frame(
            Feature = paste0(var_name, "_", clean_names),
            HR = fit_summary$coefficients[, 2],
            HR.95H = fit_summary$conf.int[, 4],
            HR.95L = fit_summary$conf.int[, 3],
            COX_Pval = fit_summary$coefficients[, 5],
            Cindex = rep(fit_summary$concordance[1], n_coef),
            Cindex_se = rep(fit_summary$concordance[2], n_coef),
            row.names = NULL,
            stringsAsFactors = FALSE
          )
        } else {
          result <- data.frame(
            Feature = var_name,
            HR = fit_summary$coefficients[, 2],
            HR.95H = fit_summary$conf.int[, 4],
            HR.95L = fit_summary$conf.int[, 3],
            COX_Pval = fit_summary$coefficients[, 5],
            Cindex = fit_summary$concordance[1],
            Cindex_se = fit_summary$concordance[2],
            row.names = NULL,
            stringsAsFactors = FALSE
          )
        }

        return(result)
      },
      error = function(e) {
        warning("Failed to fit Cox model for variable: ", var_name,
          "\n  Error: ", e$message,
          call. = FALSE
        )
        return(NULL)
      }
    )
  }

  # Perform computation (parallel or sequential)
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("furrr", quietly = TRUE)) {
      warning("Packages 'future' and 'furrr' required for parallel computation. ",
        "Running sequentially.",
        call. = FALSE
      )
      parallel <- FALSE
    }
  }

  if (parallel) {
    future::plan(future::multisession, workers = n_cores)
    result <- furrr::future_map_dfr(vars, fit_cox_model)
  } else {
    result <- do.call(rbind, lapply(vars, fit_cox_model))
  }

  # Remove NULL results
  result <- result[!is.na(result$Feature), ]

  if (nrow(result) == 0) {
    stop("Failed to fit Cox models for any variables", call. = FALSE)
  }

  return(result)
}

# =============================================================================
# Shared Forest Plot Functions
# =============================================================================

#' Prepare forest plot data
#' @description Shared data preparation for both forest plot types
#' @keywords internal
.prepare_forest_data <- function(data, point_colors, digits = NULL) {
  # Order by HR
  data <- data[order(data$HR), ]
  data$Feature <- factor(data$Feature, levels = data$Feature)

  # Determine type based on P value and HR
  data$Type <- ifelse(data$COX_Pval > 0.05, "NoSig",
    ifelse(data$HR > 1, "Risky", "Protective")
  )
  data$Type <- factor(data$Type, levels = c("Risky", "Protective", "NoSig"))

  # Set color names
  names(point_colors) <- c("Risky", "Protective", "NoSig")

  # Add formatted text if digits provided (for forest2)
  if (!is.null(digits)) {
    data$HRS <- sprintf(
      "%.*f [%.*f-%.*f]",
      digits, data$HR,
      digits, data$HR.95L,
      digits, data$HR.95H
    )
    data$Plabel <- signif(data$COX_Pval, digits)
  }

  list(data = data, point_colors = point_colors)
}

#' Create base forest plot
#' @description Shared base plot creation for forest plots
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_vline geom_errorbar geom_point scale_color_manual
.create_forest_base <- function(
    data,
    point_colors,
    point_size,
    point_border_width,
    point_border_size,
    line_colors,
    line_type,
    line_width,
    cutoff_vline_type,
    cutoff_vline_width,
    cutoff_vline_color,
    x_log10_scale) {
  # Get x-axis scale
  x_scale <- if (x_log10_scale) {
    ggplot2::scale_x_log10()
  } else {
    ggplot2::scale_x_continuous()
  }

  # Create error bar geom
  error_bar <- if (all(point_colors == line_colors)) {
    ggplot2::geom_errorbar(
      ggplot2::aes(xmax = .data$HR.95H, xmin = .data$HR.95L, color = .data$Type),
      linetype = line_type,
      linewidth = line_width, width = 0
    )
  } else {
    ggplot2::geom_errorbar(
      ggplot2::aes(xmax = .data$HR.95H, xmin = .data$HR.95L),
      linetype = line_type,
      color = line_colors[1],
      linewidth = line_width, width = 0
    )
  }

  # Build base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$HR, y = .data$Feature)) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = cutoff_vline_type,
      linewidth = cutoff_vline_width,
      color = cutoff_vline_color
    ) +
    error_bar +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Type),
      shape = 15,
      size = point_size
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Type),
      shape = 22,
      size = point_border_size,
      stroke = point_border_width
    ) +
    ggplot2::scale_color_manual(values = point_colors) +
    x_scale

  list(plot = p, x_scale = x_scale, error_bar = error_bar)
}

#' Create text geom for forest2
#' @description Helper to create text geoms with optional coloring
#' @keywords internal
#' @importFrom ggplot2 geom_text aes
.create_text_geom <- function(mapping, text_size, point_colors, text_colors, ...) {
  if (all(point_colors == text_colors)) {
    ggplot2::geom_text(mapping, size = text_size, show.legend = FALSE, ...)
  } else {
    ggplot2::geom_text(mapping, color = text_colors[1], size = text_size, show.legend = FALSE, ...)
  }
}

# =============================================================================
# Plot Functions
# =============================================================================

#' Cox Curve Plot
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_text labs
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous expansion
#' @importFrom survival coxph Surv
#' @importFrom stats predict as.formula complete.cases
#' @importFrom dplyr arrange mutate %>%
#' @importFrom rlang .data "%||%"
.cox_curve_plot <- function(
    data,
    var,
    time,
    event,
    ribbon_color,
    ribbon_alpha,
    line_color,
    line_type,
    line_width,
    text_size,
    text_digit,
    text_face,
    xlab,
    ylab,
    title,
    show_cindex,
    theme,
    theme_args,
    aspect.ratio) {
  # Rename columns
  data_clean <- data %>%
    dplyr::rename(time = !!time, event = !!event)
  data_clean <- data_clean[complete.cases(data_clean[, c("time", "event", var)]), ]

  # Fit Cox model
  formula <- as.formula(paste("survival::Surv(time, event) ~", var))
  aCox <- survival::coxph(formula, data = data_clean)

  # Extract statistics
  fit_summary <- summary(aCox)
  p_val <- fit_summary$coefficients[, 5]
  coef <- fit_summary$coefficients[, 1]
  cindex <- fit_summary$concordance[1]

  # Predict risk
  aPred <- stats::predict(aCox, type = "risk", se.fit = TRUE)
  hr <- aPred$fit
  high <- hr + 2 * aPred$se.fit
  low <- hr - 2 * aPred$se.fit

  # Create result dataframe
  # Convert var to numeric for plotting (handle both numeric and factor)
  var_values <- data_clean[[var]]
  if (is.factor(var_values) || is.character(var_values)) {
    var_numeric <- as.numeric(as.factor(var_values))
    var_labels <- levels(as.factor(var_values))
  } else {
    var_numeric <- var_values
    var_labels <- NULL
  }

  dd <- data.frame(
    var_original = var_values,
    var = var_numeric,
    HR = hr,
    HR.95H = high,
    HR.95L = low
  ) %>%
    dplyr::arrange(.data$var) %>%
    dplyr::mutate(x = 1:length(.data$var))

  # Prepare label text (use first coefficient if multiple)
  coef_first <- coef[1]
  p_val_first <- p_val[1]

  lb <- if (show_cindex) {
    paste0(
      "HR = ", signif(exp(coef_first), text_digit), "\n",
      "P = ", signif(p_val_first, text_digit), "\n",
      "C = ", signif(cindex, text_digit)
    )
  } else {
    paste0(
      "HR = ", signif(exp(coef_first), text_digit), "\n",
      "P = ", signif(p_val_first, text_digit)
    )
  }

  # Determine text position (use first coefficient if multiple)
  coef_sign <- coef[1]
  text_data <- if (coef_sign > 0) {
    data.frame(var = min(dd$var), HR.95H = max(dd$HR.95H))
  } else {
    data.frame(var = max(dd$var), HR.95H = max(dd$HR.95H))
  }
  text_hjust <- if (coef_sign > 0) -0.1 else 1.1

  # Get base size
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")

  # Generate plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = dd,
      ggplot2::aes(x = .data$var, y = .data$HR, ymin = .data$HR.95L, ymax = .data$HR.95H),
      fill = ribbon_color, alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(
      data = dd,
      ggplot2::aes(x = .data$var, y = .data$HR),
      linewidth = line_width, color = line_color, linetype = line_type
    ) +
    ggplot2::geom_text(
      data = text_data,
      ggplot2::aes(x = .data$var, y = .data$HR.95H),
      label = lb, size = text_size,
      vjust = 1, hjust = text_hjust,
      fontface = text_face, color = "black"
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.01, 0.01))) +
    do.call(theme, theme_args)

  # Add appropriate x scale
  if (!is.null(var_labels)) {
    # Factor variable: use discrete scale with labels
    p <- p + ggplot2::scale_x_continuous(
      breaks = seq_along(var_labels),
      labels = var_labels,
      expand = ggplot2::expansion(c(0, 0))
    )
  } else {
    # Numeric variable: use continuous scale
    p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0)))
  }

  # Apply styling
  p <- apply_style_theme(
    plot = p,
    data = dd,
    x_var = "var",
    y_var = "HR",
    flip = FALSE,
    base_size = base_size,
    legend.position = "none"
  )

  # Add custom theme adjustments
  if (!is.null(aspect.ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect.ratio)
  }

  return(p)
}

#' Cox Forest Plot (Simple)
#' @keywords internal
#' @importFrom ggplot2 labs theme
#' @importFrom rlang "%||%" .data
.cox_forest_plot <- function(
    data,
    point_colors,
    point_size,
    point_border_width,
    point_border_size,
    line_colors,
    line_type,
    line_width,
    cutoff_vline_type,
    cutoff_vline_width,
    cutoff_vline_color,
    x_log10_scale,
    xlab,
    title,
    subtitle,
    legend.position,
    legend.direction,
    theme,
    theme_args,
    aspect.ratio) {
  # Prepare data (shared function)
  prepared <- .prepare_forest_data(data, point_colors)
  data <- prepared$data
  point_colors <- prepared$point_colors

  # Get base size
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")

  # Create base forest plot (shared function)
  base_plot <- .create_forest_base(
    data, point_colors, point_size, point_border_width, point_border_size,
    line_colors, line_type, line_width,
    cutoff_vline_type, cutoff_vline_width, cutoff_vline_color,
    x_log10_scale
  )

  # Add labels and theme
  p <- base_plot$plot +
    ggplot2::labs(x = xlab, y = NULL, color = NULL, title = title, subtitle = subtitle) +
    do.call(theme, theme_args)

  # Apply styling
  p <- apply_style_theme(
    plot = p,
    data = data,
    x_var = "HR",
    y_var = "Feature",
    flip = FALSE,
    base_size = base_size,
    legend.position = legend.position,
    legend.direction = legend.direction
  )

  # Add custom theme adjustments
  if (!is.null(aspect.ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect.ratio)
  }

  return(p)
}

#' Cox Forest Plot (Detailed)
#' @keywords internal
#' @importFrom ggplot2 ggplot aes scale_y_continuous coord_flip theme_void
#' @importFrom ggplot2 labs theme element_blank element_text margin
#' @importFrom patchwork plot_layout
#' @importFrom rlang "%||%" .data
.cox_forest2_plot <- function(
    data,
    point_colors,
    point_size,
    point_border_width,
    point_border_size,
    line_colors,
    line_type,
    line_width,
    cutoff_vline_type,
    cutoff_vline_width,
    cutoff_vline_color,
    x_log10_scale,
    text_colors,
    digits,
    rel_width,
    theme,
    theme_args) {
  # Prepare data (shared function, with formatted text)
  prepared <- .prepare_forest_data(data, point_colors, digits = digits)
  data <- prepared$data
  point_colors <- prepared$point_colors

  # Get base size
  base_size <- theme_args$base_size %||% ggforge_option("theme.base_size")
  text_size_scale <- base_size / 12
  text_size <- 4 * text_size_scale

  # Plot 1: Variable names
  p1 <- ggplot2::ggplot(data) +
    .create_text_geom(
      ggplot2::aes(x = .data$Feature, y = 0, label = .data$Feature, color = .data$Type),
      text_size, point_colors, text_colors,
      hjust = 0
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02)) +
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::labs(x = NULL, y = NULL, title = "Variable") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5),
      plot.title = build_element_text("font.title", base_size, hjust = 0, vjust = 2)
    )

  # Plot 2: HR and confidence intervals (use shared base function)
  base_plot <- .create_forest_base(
    data, point_colors, point_size, point_border_width, point_border_size,
    line_colors, line_type, line_width,
    cutoff_vline_type, cutoff_vline_width, cutoff_vline_color,
    x_log10_scale
  )

  p2 <- base_plot$plot +
    ggplot2::labs(x = NULL, y = NULL, color = NULL, title = "Hazard ratio") +
    do.call(theme, theme_args) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = build_element_text("font.axis_text.continuous", base_size),
      axis.title = build_element_text("font.axis_title", base_size),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 0, 5.5, 0),
      plot.title = build_element_text("font.title", base_size, hjust = 0.5, vjust = 0),
      legend.position = "none",
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )

  # Plot 3: HR text
  p3 <- ggplot2::ggplot(data) +
    .create_text_geom(
      ggplot2::aes(x = .data$Feature, y = 0.01, label = .data$HRS, color = .data$Type),
      text_size, point_colors, text_colors
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02)) +
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::labs(x = NULL, y = NULL, title = "HR (95%CI)") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 0, 5.5, 0),
      plot.title = build_element_text("font.title", base_size, hjust = 0.5, vjust = 2)
    )

  # Plot 4: P-value
  p4 <- ggplot2::ggplot(data) +
    .create_text_geom(
      ggplot2::aes(x = .data$Feature, y = 0.01, label = .data$Plabel, color = .data$Type),
      text_size, point_colors, text_colors
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02)) +
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::labs(x = NULL, y = NULL, title = "Pval") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 0),
      plot.title = build_element_text("font.title", base_size, hjust = 0.5, vjust = 2)
    )

  # Combine plots
  return(p1 + p2 + p3 + p4 + patchwork::plot_layout(nrow = 1, widths = rel_width))
}
