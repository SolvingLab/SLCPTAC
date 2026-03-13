# ==============================================================================
# Plotting Layer - Complete Implementation Following Original Style
# ==============================================================================
# All plotting functions return list(plot, width, height)
# Dynamic sizing based on y-axis text length and number of items
# ==============================================================================


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Format caption with automatic line breaks
#' @description
#' Breaks long captions at comma+space boundaries when exceeding max_len
#' @keywords internal
.format_caption <- function(caption_text, max_len = 80) {
  if (is.null(caption_text) || nchar(caption_text) <= max_len) {
    return(caption_text)
  }

  # Split by comma+space
  parts <- strsplit(caption_text, ", ")[[1]]

  if (length(parts) == 1) {
    return(caption_text) # No commas to split on
  }

  formatted_lines <- ""
  current_line <- parts[1]

  for (i in 2:length(parts)) {
    test_line <- paste0(current_line, ", ", parts[i])
    if (nchar(test_line) > max_len) {
      formatted_lines <- paste0(formatted_lines, current_line, ",\n")
      current_line <- parts[i]
    } else {
      current_line <- test_line
    }
  }
  formatted_lines <- paste0(formatted_lines, current_line)

  return(formatted_lines)
}


# ==============================================================================
# Scenario 1: 1 continuous vs 1 continuous
# ==============================================================================

#' Plot for Scenario 1 (Internal)
#' @keywords internal
.plot_scenario1 <- function(data, stats, var1_feat, var2_feat, scenario_info) {
  n_cancers <- scenario_info$n_cancers
  cancer_type <- unique(data$cancer_type)

  # Convert feature labels to column names
  var1_col <- .extract_colname_from_label(c(var1_feat), data)[1]
  var2_col <- .extract_colname_from_label(c(var2_feat), data)[1]

  if (n_cancers == 1) {
    # Single cancer: CorPlot
    plot_data <- data[complete.cases(data[, c(var1_col, var2_col)]), ]

    # 提取基因名和模态，去掉癌种信息
    # 格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"
    xlab <- gsub(",\\s*[^)]+\\)", ")", var1_feat)
    ylab <- gsub(",\\s*[^)]+\\)", ")", var2_feat)

    plot <- CorPlot(
      data = plot_data,
      x = var1_col,
      y = var2_col,
      title = paste0("CPTAC-", cancer_type[1]),
      xlab = xlab,
      ylab = ylab,
      anno_items = c("pearson", "p", "n")
    )

    width <- 4.5
    height <- 4
  } else {
    # Multi-cancer: LollipopPlot
    stats_plot <- stats[!is.na(stats$r), ]
    stats_plot$neg_log10_p <- -log10(pmax(stats_plot$p, .Machine$double.xmin))

    # 场景1：使用var1_feature作为y轴（显示每个癌种的var1特征）
    stats_plot$y_label <- factor(stats_plot$var1_feature, levels = rev(unique(stats_plot$var1_feature)))

    # 从var2_feat提取基因名和模态（去掉癌种信息）
    # 格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"
    var2_label <- gsub(",\\s*[^)]+\\)", ")", var2_feat)

    plot <- LollipopPlot(
      data = stats_plot,
      x = "r",
      y = "y_label",
      fill_by = "neg_log10_p",
      size_by = 3,
      fill_name = "-log10(p)",
      palette = "RdBu",
      title = "CPTAC-Database",
      xlab = paste0("Correlation with ", var2_label),
      ylab = ""
    )

    width <- 6
    height <- max(2 + n_cancers * 0.3, 3)
  }

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Scenario 2: 1 continuous vs multiple continuous
# ==============================================================================

#' Plot for Scenario 2 (Internal)
#' @keywords internal
.plot_scenario2 <- function(data, stats, var1_features, var2_features, scenario_info) {
  n_cancers <- scenario_info$n_cancers
  cancer_type <- unique(data$cancer_type)

  # Determine which is single vs multiple
  if (scenario_info$n_var1 == 1) {
    single_var <- var1_features[1]
    multi_vars <- var2_features
    feature_col <- "var2_feature"
  } else {
    single_var <- var2_features[1]
    multi_vars <- var1_features
    feature_col <- "var1_feature"
  }

  if (n_cancers == 1) {
    # Single cancer: LollipopPlot
    stats_plot <- stats[!is.na(stats$r), ]
    if (nrow(stats_plot) == 0) {
      stop("No valid correlations found for plotting", call. = FALSE)
    }
    p_col <- if ("p_value" %in% colnames(stats_plot)) "p_value" else "p"
    stats_plot$neg_log10_p <- -log10(pmax(stats_plot[[p_col]], .Machine$double.xmin))

    # 提取特征并去掉癌种信息（格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"）
    stats_plot$feature <- gsub(",\\s*[^)]+\\)", ")", stats_plot[[feature_col]])
    stats_plot$feature <- factor(stats_plot$feature, levels = rev(unique(stats_plot$feature)))

    # x轴标题：显示与单个变量的相关性
    single_var_label <- gsub(",\\s*[^)]+\\)", ")", single_var)

    plot <- LollipopPlot(
      data = stats_plot,
      x = "r",
      y = "feature",
      fill_by = "neg_log10_p",
      size_by = 3,
      fill_name = "-log10(p)",
      palette = "RdBu",
      title = paste0("CPTAC-", cancer_type[1]),
      xlab = paste0("Correlation with ", single_var_label),
      ylab = ""
    )

    n_vars <- length(multi_vars)
    width <- 7
    height <- max(2 + n_vars * 0.3, 3)
  } else {
    # Multi-cancer: DotPlot
    stats_plot <- stats[!is.na(stats$r), ]
    if (nrow(stats_plot) == 0) {
      stop("No valid correlations found for plotting", call. = FALSE)
    }
    p_col <- if ("p_value" %in% colnames(stats_plot)) "p_value" else "p"
    stats_plot$neg_log10_p <- -log10(pmax(stats_plot[[p_col]], .Machine$double.xmin))

    # 使用 var1_feature 和 var2_feature 展示完整的相关性矩阵
    plot <- DotPlot(
      data = stats_plot,
      x = "var1_feature",
      y = "var2_feature",
      size_by = "neg_log10_p",
      fill_by = "r",
      size_name = "-log10(p)",
      fill_name = "Correlation",
      palette = "RdBu",
      x_text_angle = 45,
      title = "CPTAC-Database",
      xlab = "",
      ylab = ""
    )

    # 计算图形尺寸
    n_var1_features <- length(unique(stats_plot$var1_feature))
    n_var2_features <- length(unique(stats_plot$var2_feature))
    width <- max(3 + n_var1_features * 0.45, 5)
    height <- max(3 + n_var2_features * 0.3, 4)
  }

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Scenario 3: Multiple continuous vs multiple continuous
# ==============================================================================

#' Plot for Scenario 3 (Internal)
#' @keywords internal
.plot_scenario3 <- function(data, stats, var1_features, var2_features, scenario_info) {
  n_cancers <- scenario_info$n_cancers
  cancer_type <- unique(data$cancer_type)

  stats_plot <- stats[!is.na(stats$r), ]

  # 过滤掉自己和自己的相关性（对角线）
  stats_plot <- stats_plot[stats_plot$var1_feature != stats_plot$var2_feature, ]

  stats_plot$neg_log10_p <- -log10(pmax(stats_plot$p, .Machine$double.xmin))

  if (n_cancers == 1) {
    # Single cancer: DotPlot
    # 去掉特征标签中的癌种信息（格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"）
    stats_plot$var1_label <- gsub(",\\s*[^)]+\\)", ")", stats_plot$var1_feature)
    stats_plot$var2_label <- gsub(",\\s*[^)]+\\)", ")", stats_plot$var2_feature)

    plot <- DotPlot(
      data = stats_plot,
      x = "var1_label",
      y = "var2_label",
      size_by = "neg_log10_p",
      fill_by = "r",
      size_name = "-log10(p)",
      fill_name = "Correlation",
      palette = "RdBu",
      x_text_angle = 45,
      title = paste0("CPTAC-", cancer_type[1]),
      xlab = "",
      ylab = ""
    )

    n_var1 <- length(var1_features)
    n_var2 <- length(var2_features)
    width <- max(4 + n_var1 * 0.45, 5)
    height <- max(3 + n_var2 * 0.3, 4)
  } else {
    # Multi-cancer: DotPlot with expanded y-axis
    # var2_feature已经包含癌种信息，无需再添加

    plot <- DotPlot(
      data = stats_plot,
      x = "var1_feature",
      y = "var2_feature",
      size_by = "neg_log10_p",
      fill_by = "r",
      size_name = "-log10(p)",
      fill_name = "Correlation",
      palette = "RdBu",
      x_text_angle = 45,
      title = "CPTAC-Database",
      xlab = "",
      ylab = ""
    )

    n_var1 <- length(var1_features)
    n_var2 <- length(var2_features)
    n_cancers_in_data <- length(unique(stats_plot$cancer_type))
    # 多癌种情况下，y轴特征数 = var2数量 × 癌种数量
    n_var2_total <- length(unique(stats_plot$var2_feature))
    width <- max(4 + n_var1 * 0.35, 5)
    height <- max(3 + n_var2_total * 0.35, 4)
  }

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Scenario 4: 1 continuous vs 1 categorical → BoxPlot
# ==============================================================================

#' Plot for Scenario 4 (Internal)
#' @keywords internal
.plot_scenario4 <- function(data, stats, var1_features, var2_features, scenario_info) {
  # Determine which is continuous vs categorical
  if (scenario_info$var1_class == "continuous") {
    con_var <- var1_features[1]
    cat_var <- var2_features[1]
  } else {
    con_var <- var2_features[1]
    cat_var <- var1_features[1]
  }

  # Extract column names
  con_col <- .extract_colname_from_label(c(con_var), data)[1]
  cat_col <- .extract_colname_from_label(c(cat_var), data)[1]

  n_cancers <- length(unique(data$cancer_type))
  cancer_type <- unique(data$cancer_type)

  if (n_cancers == 1) {
    # Single cancer: BoxPlot
    # .create_classic_boxplot expects (data, cat_var, con_var, title, xlab, ylab, stats)
    title <- paste0("CPTAC-", cancer_type[1])

    # X轴（分类变量）：只显示变量名
    xlab <- gsub("\\s*\\(.*", "", cat_var)
    # Y轴（连续变量）：保留模态信息，去掉癌种
    # 格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"
    ylab <- gsub(",\\s*[^)]+\\)", ")", con_var)

    plot <- .create_classic_boxplot(data, cat_col, con_col, title, xlab, ylab, stats)

    # Get number of groups for dynamic width
    n_groups <- length(unique(data[[cat_col]][!is.na(data[[cat_col]])]))
    width <- max(2.3 + n_groups * 0.35, 3)
    height <- 4.5
  } else {
    # Multi-cancer: multiple box plots
    plot_list <- list()

    # X轴（分类变量）：只显示变量名
    xlab <- gsub("\\s*\\(.*", "", cat_var)
    # Y轴（连续变量）：保留模态信息，去掉癌种
    ylab <- gsub(",\\s*[^)]+\\)", ")", con_var)

    for (ct in cancer_type) {
      ct_data <- data[data$cancer_type == ct, ]

      # 为每个癌种提取各自的列名（关键修复！）
      ct_con_cols <- grep(paste0("^", ct, "_.*_", gsub(".*_", "", con_col)), colnames(ct_data), value = TRUE)
      ct_cat_cols <- grep(paste0("^", ct, "_.*_", gsub(".*_", "", cat_col)), colnames(ct_data), value = TRUE)

      if (length(ct_con_cols) > 0 && length(ct_cat_cols) > 0 && nrow(ct_data) >= 3) {
        ct_con_col <- ct_con_cols[1]
        ct_cat_col <- ct_cat_cols[1]

        ct_stats <- stats[grepl(ct, stats$var1_feature) | grepl(ct, stats$var2_feature), ]
        p <- .create_classic_boxplot(ct_data, ct_cat_col, ct_con_col, paste0("CPTAC-", ct), xlab, ylab, ct_stats)
        plot_list[[ct]] <- p
      }
    }

    if (length(plot_list) == 0) {
      stop("No valid plots generated for any cancer type", call. = FALSE)
    }

    layout <- .calc_facet_layout(length(plot_list))
    plot <- patchwork::wrap_plots(plot_list, ncol = layout$ncol, nrow = layout$nrow)

    # Smaller width per plot for multiple boxplots
    width <- layout$ncol * 3
    height <- layout$nrow * 4.5
  }

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Scenario 5-6: Multiple BoxPlots
# S5: 1 continuous vs 多 categorical
# S6: 多 continuous vs 1 categorical
# ==============================================================================

#' Plot for Scenarios 5-6 (Internal)
#' @keywords internal
.plot_scenario5_6 <- function(data, stats, var1_features, var2_features, scenario_info) {
  n_cancers <- length(unique(data$cancer_type))
  cancer_type <- unique(data$cancer_type)

  # Determine which is single continuous vs multiple categorical (or vice versa)
  if (scenario_info$var1_class == "continuous" && length(var1_features) == 1) {
    # S5: 1 continuous vs 多 categorical
    con_vars <- var1_features
    cat_vars <- var2_features
  } else if (scenario_info$var2_class == "continuous" && length(var2_features) == 1) {
    # S5: 多 categorical vs 1 continuous
    con_vars <- var2_features
    cat_vars <- var1_features
  } else if (scenario_info$var1_class == "continuous") {
    # S6: 多 continuous vs 1 categorical
    con_vars <- var1_features
    cat_vars <- var2_features
  } else {
    # S6: 1 categorical vs 多 continuous
    con_vars <- var2_features
    cat_vars <- var1_features
  }

  # Create multiple boxplots
  plot_list <- list()

  for (con_var in con_vars) {
    for (cat_var in cat_vars) {
      con_col <- .extract_colname_from_label(c(con_var), data)[1]
      cat_col <- .extract_colname_from_label(c(cat_var), data)[1]

      if (con_col %in% colnames(data) && cat_col %in% colnames(data)) {
        # Get stats for this pair
        pair_stats <- stats[
          (grepl(gsub("\\s*\\(.*", "", con_var), stats$var1_feature) &
            grepl(gsub("\\s*\\(.*", "", cat_var), stats$var2_feature)) |
            (grepl(gsub("\\s*\\(.*", "", cat_var), stats$var1_feature) &
              grepl(gsub("\\s*\\(.*", "", con_var), stats$var2_feature)),
        ]

        # Extract cancer type from feature label
        cancer_match <- regmatches(con_var, regexpr("\\([^,]+,\\s*([^)]+)\\)", con_var))
        if (length(cancer_match) > 0) {
          ct <- gsub(".*,\\s*([^)]+)\\)", "\\1", cancer_match)
        } else {
          ct <- cancer_type[1]
        }

        # .create_classic_boxplot expects (data, cat_var, con_var, title, xlab, ylab, stats)
        title <- paste0("CPTAC-", ct)

        # X轴（分类变量）：只显示变量名
        xlab <- gsub("\\s*\\(.*", "", cat_var)
        # Y轴（连续变量）：保留模态信息，去掉癌种
        # 格式: "GENE (Modal, Cancer)" -> "GENE (Modal)"
        ylab <- gsub(",\\s*[^)]+\\)", ")", con_var)

        p <- .create_classic_boxplot(data, cat_col, con_col, title, xlab, ylab, pair_stats)
        plot_name <- paste(gsub("\\s*\\(.*", "", con_var), gsub("\\s*\\(.*", "", cat_var), sep = "_vs_")
        plot_list[[plot_name]] <- p
      }
    }
  }

  if (length(plot_list) == 0) {
    stop("No valid plots generated for any variable combination", call. = FALSE)
  }

  layout <- .calc_facet_layout(length(plot_list))
  plot <- patchwork::wrap_plots(plot_list, ncol = layout$ncol, nrow = layout$nrow)

  # Smaller width per plot for multiple boxplots
  width <- layout$ncol * 3
  height <- layout$nrow * 4.5

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Scenario 7: 多 continuous vs 多 categorical → DotPlot/Heatmap
# ==============================================================================

#' Plot for Scenario 7 (Internal)
#' @keywords internal
.plot_scenario7 <- function(data, stats, var1_features, var2_features, scenario_info) {
  n_cancers <- length(unique(data$cancer_type))
  cancer_type <- unique(data$cancer_type)

  # For multiple categorical vs multiple categorical OR mixed scenarios
  # Use categorical heatmap or DotPlot

  if (scenario_info$var1_class == "categorical" && scenario_info$var2_class == "categorical") {
    # Both categorical: use barplot for 1v1, 1v多, 多v1; use heatmap only for 多v多
    n_var1 <- length(var1_features)
    n_var2 <- length(var2_features)

    if (n_var1 == 1 || n_var2 == 1) {
      # 1 vs 1, 1 vs 多, or 多 vs 1: use percentage barplot(s)
      plot_list <- list()

      # Get all combinations
      for (v1_feat in var1_features) {
        for (v2_feat in var2_features) {
          var1_col <- .extract_colname_from_label(c(v1_feat), data)[1]
          var2_col <- .extract_colname_from_label(c(v2_feat), data)[1]

          # Get stats for this pair
          pair_stats <- stats[
            (stats$var1_feature == v1_feat & stats$var2_feature == v2_feat) |
              (stats$var1_feature == v2_feat & stats$var2_feature == v1_feat),
          ]

          # Extract cancer from feature label
          cancer_match <- regmatches(v1_feat, regexpr(",\\s*([^)]+)\\)", v1_feat))
          if (length(cancer_match) > 0) {
            ct <- gsub(",\\s*|\\)", "", cancer_match)
          } else {
            ct <- cancer_type[1]
          }

          p <- .create_percentage_barplot(
            data = data,
            var1 = var1_col,
            var2 = var2_col,
            cancer_label = ct,
            stats = pair_stats,
            var1_label = v1_feat,
            var2_label = v2_feat
          )

          plot_name <- paste(gsub("\\s*\\(.*", "", v1_feat),
            gsub("\\s*\\(.*", "", v2_feat),
            ct,
            sep = "_"
          )
          plot_list[[plot_name]] <- p
        }
      }

      if (length(plot_list) == 0) {
        stop("No valid plots generated", call. = FALSE)
      }

      if (length(plot_list) == 1) {
        # Single plot
        plot <- plot_list[[1]]
        var1_col <- .extract_colname_from_label(var1_features, data)[1]
        n_categories <- length(unique(data[[var1_col]][!is.na(data[[var1_col]])]))
        width <- 2.6 + n_categories * 0.35
        height <- 5
      } else {
        # Multiple plots: use patchwork
        layout <- .calc_facet_layout(length(plot_list))
        plot <- patchwork::wrap_plots(plot_list, ncol = layout$ncol, nrow = layout$nrow)

        width <- layout$ncol * 3.3
        height <- layout$nrow * 5
      }

      attr(plot, "width") <- width
      attr(plot, "height") <- height
      return(plot)
    } else {
      # 多 vs 多: use heatmap
      return(.plot_categorical_heatmap(data, stats, var1_features, var2_features))
    }
  } else {
    # Mixed: use DotPlot showing effect sizes
    # This is similar to S3 but for mixed variable types
    stats_plot <- stats[!is.na(stats$p_value), ]

    if (nrow(stats_plot) == 0) {
      stop("No valid statistics for plotting", call. = FALSE)
    }

    # Determine metric for plotting (correlation, effect size, etc.)
    if ("r" %in% colnames(stats_plot)) {
      metric <- "r"
      metric_name <- "Correlation"
    } else if ("effect_size" %in% colnames(stats_plot)) {
      metric <- "effect_size"
      metric_name <- "Effect Size"
    } else if ("cramers_v" %in% colnames(stats_plot)) {
      metric <- "cramers_v"
      metric_name <- "Cramer's V"
    } else {
      # Default: use p-value
      stats_plot$metric <- -log10(pmax(stats_plot$p_value, .Machine$double.xmin))
      metric <- "metric"
      metric_name <- "-log10(p)"
    }

    stats_plot$neg_log10_p <- -log10(pmax(stats_plot$p_value, .Machine$double.xmin))

    plot <- DotPlot(
      data = stats_plot,
      x = "var1_feature",
      y = "var2_feature",
      size_by = "neg_log10_p",
      fill_by = if (metric %in% colnames(stats_plot)) metric else "neg_log10_p",
      size_name = "-log10(p)",
      fill_name = metric_name,
      palette = "RdBu",
      x_text_angle = 45,
      title = if (n_cancers == 1) paste0("CPTAC-", cancer_type[1]) else "CPTAC Multi-Cancer",
      xlab = "",
      ylab = ""
    )

    n_vars1 <- length(var1_features)
    n_vars2 <- length(var2_features)
    width <- max(6, 4 + n_vars1 * 0.4)
    height <- max(5, 3 + n_vars2 * 0.4)

    attr(plot, "width") <- width
    attr(plot, "height") <- height
    return(plot)
  }

  if (n_cancers == 1) {
    # Single cancer: single boxplot
    pair_stats <- stats[stats$categorical_var == cat_var &
      stats$continuous_var == con_var, ]

    plot <- .create_classic_boxplot(
      data, cat_var, con_var,
      title = paste0("CPTAC-", cancer_type[1]),
      subtitle = paste0(con_var, " by ", cat_var),
      stats = pair_stats
    )

    n_groups <- length(unique(data[[cat_var]][!is.na(data[[cat_var]])]))
    width <- 2.3 + n_groups * 0.35
    height <- 4.5
  } else {
    # Multi-cancer: multiple boxplots
    plot_list <- list()

    for (ct in cancer_type) {
      ct_data <- data[data$cancer_type == ct, ]

      if (nrow(ct_data) < 3 || !cat_var %in% colnames(ct_data) ||
        !con_var %in% colnames(ct_data)) {
        next
      }

      valid_cat_data <- ct_data[[cat_var]][!is.na(ct_data[[cat_var]]) &
        !is.na(ct_data[[con_var]])]
      if (length(unique(valid_cat_data)) < 2) {
        next
      }

      ct_stats <- stats[stats$cancer_type == ct &
        stats$categorical_var == cat_var &
        stats$continuous_var == con_var, ]

      if (nrow(ct_stats) == 0) {
        next
      }

      p <- .create_classic_boxplot(
        ct_data, cat_var, con_var,
        title = paste0("CPTAC-", ct),
        stats = ct_stats
      )
      plot_list[[ct]] <- p
    }

    if (length(plot_list) == 0) {
      stop("No valid plots generated for any cancer type")
    }

    layout <- .calc_facet_layout(length(plot_list))
    plot <- patchwork::wrap_plots(plot_list, ncol = layout$ncol, nrow = layout$nrow)

    avg_groups <- mean(sapply(names(plot_list), function(ct) {
      if (ct %in% cancer_type) {
        ct_data <- data[data$cancer_type == ct, ]
        length(unique(ct_data[[cat_var]][!is.na(ct_data[[cat_var]])]))
      } else {
        3
      }
    }))

    single_width <- 2 + avg_groups * 0.35
    width <- layout$ncol * single_width
    height <- layout$nrow * 4.5
  }

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Create percentage barplot (Internal)
#' @keywords internal
.create_percentage_barplot <- function(data, var1, var2, cancer_label, stats = NULL,
                                       var1_label = NULL, var2_label = NULL) {
  plot_data <- data[complete.cases(data[, c(var1, var2)]), ]

  # 智能轴分配逻辑：
  # 优先级：Clinical > 其他变量 > Mutation
  # Clinical在x轴（分类），Mutation作为fill（二分类）
  should_swap <- FALSE

  is_var1_clinical <- grepl("_Clinical$", var1)
  is_var2_clinical <- grepl("_Clinical$", var2)
  is_var1_mutation <- grepl("_Mutation$", var1)
  is_var2_mutation <- grepl("_Mutation$", var2)

  if (is_var1_mutation && is_var2_clinical) {
    # var1是Mutation，var2是Clinical → 交换（Clinical在x轴）
    should_swap <- TRUE
  } else if (is_var1_mutation && !is_var2_mutation && !is_var2_clinical) {
    # var1是Mutation，var2是其他 → 交换（其他在x轴）
    should_swap <- TRUE
  } else if (is_var1_clinical && is_var2_clinical) {
    # 都是Clinical → 选择分类数多的在x轴
    n_cat1 <- length(unique(plot_data[[var1]][!is.na(plot_data[[var1]])]))
    n_cat2 <- length(unique(plot_data[[var2]][!is.na(plot_data[[var2]])]))
    if (n_cat2 > n_cat1) {
      should_swap <- TRUE
    }
  }

  if (should_swap) {
    temp <- var1
    var1 <- var2
    var2 <- temp
    # Also swap labels
    if (!is.null(var1_label) && !is.null(var2_label)) {
      temp_label <- var1_label
      var1_label <- var2_label
      var2_label <- temp_label
    }
  }

  contingency <- table(plot_data[[var1]], plot_data[[var2]])
  percentage_data <- data.frame()

  for (i in 1:nrow(contingency)) {
    for (j in 1:ncol(contingency)) {
      row_sum <- sum(contingency[i, ])
      if (row_sum > 0) {
        percentage_data <- rbind(percentage_data, data.frame(
          category = rownames(contingency)[i],
          group = colnames(contingency)[j],
          value = contingency[i, j] / row_sum * 100,
          count = contingency[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  percentage_data$category <- factor(percentage_data$category)
  percentage_data$group <- factor(percentage_data$group)

  caption_text <- NULL
  if (!is.null(stats) && nrow(stats) > 0) {
    # 构建caption，优先显示OR信息
    caption_parts <- paste0(stats$test_method[1], "\np = ", format.pval(stats$p_value[1], digits = 3))

    # 添加OR信息（如果有）
    if ("odds_ratio" %in% colnames(stats) && !is.na(stats$odds_ratio[1])) {
      caption_parts <- paste0(
        caption_parts,
        ",\nOR = ", round(stats$odds_ratio[1], 3),
        " (log2 = ", round(stats$log2_or[1], 2), ")"
      )
    } else if ("effect_size" %in% colnames(stats) && !is.na(stats$effect_size[1])) {
      caption_parts <- paste0(
        caption_parts,
        ", Effect size = ", round(stats$effect_size[1], 3)
      )
    }

    if (stats$significant[1]) {
      caption_parts <- paste0(caption_parts, " *")
    }

    caption_text <- caption_parts
  }

  # Dynamic colors
  group_levels <- levels(percentage_data$group)
  n_groups <- length(group_levels)

  if (all(group_levels %in% c("WildType", "Mutation"))) {
    fill_colors <- c("WildType" = "#41A98E", "Mutation" = "#ED6355")
  } else {
    color_palette <- c("#41A98E", "#ED6355", "#EFA63A", "#3a6ea5", "#9b59b6", "#e74c3c", "#2ecc71", "#34495e")
    fill_colors <- color_palette[1:n_groups]
    names(fill_colors) <- group_levels
  }

  # Use labels if provided, otherwise use column names
  x_axis_label <- if (!is.null(var1_label)) {
    gsub("\\s*\\(.*", "", var1_label) # Extract gene name only
  } else {
    var1
  }

  # Legend title: use var2 gene name
  legend_title <- if (!is.null(var2_label)) {
    gsub("\\s*\\(.*", "", var2_label) # Extract gene name only
  } else {
    var2
  }

  plot <- ggplot2::ggplot(percentage_data, ggplot2::aes(x = category, y = value, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(value, 1), "%\n(n=", count, ")")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3, color = "black", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.6, 0.6))) +
    ggplot2::labs(
      title = paste0("CPTAC-", cancer_label),
      x = x_axis_label,
      y = "Percentage (%)",
      fill = legend_title,
      caption = caption_text
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
      axis.title = ggplot2::element_text(size = 12, colour = "black", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12, colour = "black", face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9)
    )

  return(plot)
}


#' Create classic boxplot (Internal)
#' @keywords internal
.create_classic_boxplot <- function(data, cat_var, con_var, title, xlab, ylab, stats = NULL) {
  plot_data <- data[complete.cases(data[, c(cat_var, con_var)]), ]

  if (!is.factor(plot_data[[cat_var]])) {
    plot_data[[cat_var]] <- factor(plot_data[[cat_var]])
  }

  cat_levels <- levels(plot_data[[cat_var]])
  n_groups <- length(cat_levels)
  color_palette <- c("#41A98E", "#ED6355", "#EFA63A", "#3a6ea5", "#9b59b6", "#e74c3c", "#2ecc71", "#34495e")

  if (n_groups <= length(color_palette)) {
    colors <- color_palette[1:n_groups]
  } else {
    colors <- grDevices::colorRampPalette(color_palette)(n_groups)
  }
  names(colors) <- cat_levels

  caption_text <- NULL
  if (!is.null(stats) && nrow(stats) > 0) {
    group_sizes <- table(plot_data[[cat_var]])
    caption_text <- paste(
      paste(names(group_sizes), "=", group_sizes, collapse = ", "),
      paste0("\nTest: ", stats$test_method[1]),
      paste0("\nEffect size = ", round(stats$effect_size[1], 3)),
      paste0(", p = ", format.pval(stats$p_value[1], digits = 3)),
      ifelse(stats$significant[1], " *", ""),
      sep = ""
    )
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = !!rlang::sym(cat_var),
    y = !!rlang::sym(con_var),
    fill = !!rlang::sym(cat_var)
  )) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      caption = caption_text
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
      axis.title = ggplot2::element_text(size = 12, colour = "black", face = "bold"),
      axis.text = ggplot2::element_text(size = 11, colour = "black"),
      legend.position = "none",
      plot.caption = ggplot2::element_text(hjust = 0, size = 9)
    )

  return(p)
}


#' Calculate facet layout (Internal)
#' @keywords internal
.calc_facet_layout <- function(n_plots) {
  if (n_plots <= 4) {
    # 1-4个：一行排列
    list(ncol = n_plots, nrow = 1)
  } else if (n_plots == 5) {
    # 5个：3,2
    list(ncol = 3, nrow = 2)
  } else if (n_plots == 6) {
    # 6个：3,3
    list(ncol = 3, nrow = 2)
  } else if (n_plots == 7) {
    # 7个：4,3
    list(ncol = 4, nrow = 2)
  } else if (n_plots == 8) {
    # 8个：4,4
    list(ncol = 4, nrow = 2)
  } else if (n_plots == 9) {
    # 9个：3,3,3
    list(ncol = 3, nrow = 3)
  } else if (n_plots == 10) {
    # 10个：4,4,2
    list(ncol = 4, nrow = 3)
  } else if (n_plots == 11) {
    # 11个：4,4,3
    list(ncol = 4, nrow = 3)
  } else if (n_plots == 12) {
    # 12个：4,4,4
    list(ncol = 4, nrow = 3)
  } else {
    # >13个：保持4列
    list(ncol = 4, nrow = ceiling(n_plots / 4))
  }
}


#' Create categorical heatmap (Internal)
#' @keywords internal
.plot_categorical_heatmap <- function(data, stats, var1_features, var2_features) {
  all_cat_vars <- unique(c(var1_features, var2_features))
  n_vars <- length(all_cat_vars)

  or_matrix <- matrix(NA,
    nrow = n_vars, ncol = n_vars,
    dimnames = list(all_cat_vars, all_cat_vars)
  )
  pval_matrix <- matrix(NA,
    nrow = n_vars, ncol = n_vars,
    dimnames = list(all_cat_vars, all_cat_vars)
  )
  sig_matrix <- matrix("",
    nrow = n_vars, ncol = n_vars,
    dimnames = list(all_cat_vars, all_cat_vars)
  )

  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) {
        or_matrix[i, j] <- 0
        pval_matrix[i, j] <- 1
        next
      }

      var_i <- all_cat_vars[i]
      var_j <- all_cat_vars[j]

      # Use var1_feature and var2_feature for matching (feature labels, not column names)
      stat_row <- stats[
        (stats$var1_feature == var_i & stats$var2_feature == var_j) |
          (stats$var1_feature == var_j & stats$var2_feature == var_i),
      ]

      if (nrow(stat_row) > 0) {
        # 优先使用 log2(OR)，如果没有则使用 effect_size
        if ("log2_or" %in% colnames(stat_row) && !is.na(stat_row$log2_or[1])) {
          # 使用 log2(OR): 正值=共突变，负值=互斥
          effect_val <- stat_row$log2_or[1]
        } else if ("effect_size" %in% colnames(stat_row) && !is.na(stat_row$effect_size[1])) {
          # 备选：使用 effect_size
          effect_val <- stat_row$effect_size[1]
        } else {
          # 最后：使用 -log10(p)
          p_val <- stat_row$p_value[1]
          effect_val <- -log10(pmax(p_val, .Machine$double.xmin))
        }

        p_val <- stat_row$p_value[1]

        or_matrix[i, j] <- effect_val
        or_matrix[j, i] <- effect_val

        pval_matrix[i, j] <- p_val
        pval_matrix[j, i] <- p_val

        if (p_val < 0.001) {
          sig_matrix[i, j] <- "***"
          sig_matrix[j, i] <- "***"
        } else if (p_val < 0.01) {
          sig_matrix[i, j] <- "**"
          sig_matrix[j, i] <- "**"
        } else if (p_val < 0.05) {
          sig_matrix[i, j] <- "*"
          sig_matrix[j, i] <- "*"
        }
      }
    }
  }

  plot_matrix <- or_matrix
  plot_matrix[is.infinite(plot_matrix)] <- NA

  cancer_type <- unique(data$cancer_type)

  plot <- Heatmap(
    data = plot_matrix,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "right",
    column_names_side = "bottom",
    row_name_annotation = FALSE,
    column_name_annotation = FALSE,
    border = TRUE,
    na_col = "grey90",
    palcolor = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    cell_type = "label",
    label = function(x, i, j) {
      sig <- ComplexHeatmap::pindex(sig_matrix, i, j)
      val <- ComplexHeatmap::pindex(plot_matrix, i, j)
      ifelse(is.na(val), "",
        ifelse(val == 0, "",
          paste0(round(val, 2), sig)
        )
      )
    },
    label_size = 10,
    title = NULL,
    legend.position = "right",
    values_by = "log2(OR)"
  )

  base_size <- 4 + n_vars * 0.4
  width <- base_size + 0.7
  height <- base_size

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  attr(plot, "plot_type") <- "heatmap"
  return(plot)
}


# ==============================================================================
# ENRICHMENT PLOTS
# ==============================================================================

#' Network Plot (Internal)
#' @keywords internal
.plot_network <- function(stats, var1_name, edge_metric = "correlation",
                          query_omics = NULL, genome_omics = NULL, cancer_type = NULL,
                          analysis_type = "correlation", method = NULL) {
  # Determine target column name
  target_col <- if ("gene" %in% colnames(stats)) {
    "gene"
  } else if ("comparison_variable" %in% colnames(stats)) {
    "comparison_variable"
  } else {
    stop("Cannot find gene/variable column in stats")
  }

  # Split into positive and negative
  if (edge_metric == "logFC") {
    stats_pos <- stats[stats[[edge_metric]] > 0, ]
    stats_neg <- stats[stats[[edge_metric]] < 0, ]
  } else {
    stats_pos <- stats[stats[[edge_metric]] > 0, ]
    stats_neg <- stats[stats[[edge_metric]] < 0, ]
  }

  # Create two separate network plots
  create_single_network <- function(data_subset, color, direction_label, keep_negative = FALSE) {
    if (nrow(data_subset) == 0) {
      return(NULL)
    }

    network_data <- data.frame(
      source = var1_name,
      target = data_subset[[target_col]],
      weight = data_subset[[edge_metric]]
    )

    graph <- igraph::graph_from_data_frame(network_data, directed = FALSE)

    igraph::V(graph)$color <- ifelse(
      igraph::V(graph)$name == var1_name,
      "#00b4d8",
      color
    )

    igraph::V(graph)$size <- ifelse(igraph::V(graph)$name == var1_name, 15, 6)

    # 设置节点文本样式：中心节点加粗、更大
    igraph::V(graph)$text_size <- ifelse(igraph::V(graph)$name == var1_name, 3.5, 2.5)
    igraph::V(graph)$text_face <- ifelse(igraph::V(graph)$name == var1_name, "bold", "plain")

    # For layout, always use absolute weight (FR layout requires positive weights)
    igraph::E(graph)$abs_weight <- abs(igraph::E(graph)$weight)

    # For edge width: 左右图各自独立scale，用绝对值
    igraph::E(graph)$display_weight <- abs(igraph::E(graph)$weight)

    # Set weight for layout (must be positive)
    igraph::E(graph)$weight <- abs(igraph::E(graph)$weight)

    # For legend, show correct sign
    p <- ggraph::ggraph(graph, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(edge_width = display_weight), edge_color = "#BEBEBE99") +
      ggraph::geom_node_point(ggplot2::aes(size = size, color = color), show.legend = FALSE) +
      ggraph::geom_node_text(ggplot2::aes(label = name, size = text_size, fontface = text_face),
        color = "black", repel = FALSE, show.legend = FALSE
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_size_identity() +
      ggraph::scale_edge_width(
        range = c(0.3, 1.6),
        name = if (analysis_type == "DEA") "logFC" else "Correlation",
        labels = function(x) if (keep_negative) -abs(x) else abs(x) # 右图显示负值
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(title = direction_label) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", colour = "black")
      )

    return(p)
  }

  # Positive plot: use absolute values in legend
  # Negative plot: keep negative values in legend
  plot_pos <- create_single_network(stats_pos, "#e76f51",
    if (analysis_type == "DEA") "Up-regulated (logFC > 0)" else "Positive (r > 0)",
    keep_negative = FALSE
  )
  plot_neg <- create_single_network(stats_neg, "#41A98E",
    if (analysis_type == "DEA") "Down-regulated (logFC < 0)" else "Negative (r < 0)",
    keep_negative = TRUE
  )

  # Build main title and caption
  title_text <- if (!is.null(cancer_type)) {
    paste0("CPTAC-", cancer_type)
  } else {
    "CPTAC Network Plot"
  }

  subtitle_text <- if (!is.null(query_omics) && !is.null(genome_omics)) {
    paste0(
      "Top 50 up/down ", genome_omics, " genes ",
      if (analysis_type == "DEA") "associated with" else "correlated with",
      " ", var1_name, " (", query_omics, ")"
    )
  } else {
    NULL
  }

  caption_text <- if (analysis_type == "DEA") {
    "Method: Differential expression (limma), FDR < 0.05"
  } else if (!is.null(method)) {
    paste0("Method: ", tools::toTitleCase(method), " correlation, FDR < 0.05")
  } else {
    NULL
  }

  # Combine two plots
  if (!is.null(plot_pos) && !is.null(plot_neg)) {
    plot <- patchwork::wrap_plots(plot_pos, plot_neg, ncol = 2) +
      patchwork::plot_annotation(
        title = title_text,
        subtitle = subtitle_text,
        caption = caption_text,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
          plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray40")
        )
      )
  } else if (!is.null(plot_pos)) {
    plot <- plot_pos + ggplot2::labs(title = title_text, subtitle = subtitle_text, caption = caption_text)
  } else if (!is.null(plot_neg)) {
    plot <- plot_neg + ggplot2::labs(title = title_text, subtitle = subtitle_text, caption = caption_text)
  } else {
    stop("No data to plot")
  }

  # NetworkPlot: 使用属性附加width和height
  attr(plot, "width") <- 14
  attr(plot, "height") <- 7

  return(plot)
}


#' GSEA Paired DotPlot (Internal)
#' @keywords internal
.plot_gsea_paired <- function(gsea_stats, var_name, omics_type, cancer_types, enrich_type, GO_ont = NULL, method = NULL, top_n = 15) {
  library(dplyr)
  library(rlang)

  # 使用top_n参数控制显示的pathway数量
  pos_pathways <- gsea_stats %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice(1:min(top_n, n()))
  neg_pathways <- gsea_stats %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice(1:min(top_n, n()))

  # Build title
  cancer_str <- if (length(cancer_types) > 1) "Database" else cancer_types[1]
  main_title <- paste0(enrich_type, " Enrichment of ", var_name, " ", omics_type, " in CPTAC-", cancer_str)

  # Build caption
  if (omics_type == "Mutation") {
    caption_text <- paste0("Method: DEA (limma) → GSEA (", enrich_type)
  } else if (!is.null(method)) {
    caption_text <- paste0("Method: ", tools::toTitleCase(method), " correlation → GSEA (", enrich_type)
  } else {
    caption_text <- paste0("Method: GSEA (", enrich_type)
  }

  if (!is.null(GO_ont) && enrich_type == "GO") {
    caption_text <- paste0(caption_text, " ", GO_ont)
  }
  caption_text <- paste0(caption_text, ")")

  if (nrow(pos_pathways) > 0 && nrow(neg_pathways) > 0) {
    # 根据变量类型设置子图标题
    if (omics_type == "Mutation") {
      pos_title <- "Positive Terms (Mutation vs WildType)"
      neg_title <- "Negative Terms (Mutation vs WildType)"
    } else {
      pos_title <- "Positive Terms (Positive Correlation)"
      neg_title <- "Negative Terms (Negative Correlation)"
    }

    # Create separate plots for positive and negative
    plot_pos <- .custom_gsea_dot(
      res = pos_pathways,
      type = "Positive",
      show.term.num = nrow(pos_pathways),
      Select.P = "NP",
      cutoff.P = 1,
      title = pos_title,
      legend.position = "right"
    )

    # 修复：负NES使用BrBG色系的reverse，不用abs
    brbrg_colors <- c(
      "#003c30", "#01665e", "#35978f", "#80cdc1", "#c7eae5",
      "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005"
    )

    plot_neg <- .custom_gsea_dot(
      res = neg_pathways,
      type = "Negative",
      show.term.num = nrow(neg_pathways),
      Select.P = "NP",
      cutoff.P = 1,
      title = neg_title,
      legend.position = "right",
      colors = rev(brbrg_colors) # reverse色系：NES越负越深，legend显示负值
    )

    # Combine with patchwork
    plot <- patchwork::wrap_plots(plot_pos, plot_neg, ncol = 2) +
      patchwork::plot_annotation(
        title = main_title,
        caption = caption_text,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, colour = "black"),
          plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray40")
        )
      )

    n_pathways <- max(nrow(pos_pathways), nrow(neg_pathways))

    # 修复：固定基础宽度 + 基于最大ytext字符长度
    all_pathway_labels <- c(as.character(pos_pathways$Description), as.character(neg_pathways$Description))
    max_label_length <- if (length(all_pathway_labels) > 0) max(nchar(all_pathway_labels), na.rm = TRUE) else 50

    # 固定基础宽度 + 字符长度调整
    base_width <- 14 # 固定基础宽度
    width_adjust <- max_label_length * 0.10 # 基于最大字符长度
    width <- base_width + width_adjust
    width <- min(width, 22) # 设置上限

    height <- max(6, 3 + n_pathways * 0.25)
  } else if (nrow(pos_pathways) > 0) {
    # 只有positive pathways
    pos_title <- if (omics_type == "Mutation") {
      "Positive Terms (Mutation vs WildType)"
    } else {
      "Positive Terms (Positive Correlation)"
    }

    pos_pathways$neg_log10_p <- -log10(pmax(pos_pathways$pvalue, .Machine$double.xmin))
    plot <- DotPlot(
      data = pos_pathways,
      x = "NES",
      y = "Description",
      size_by = "neg_log10_p",
      fill_by = "NES",
      size_name = "-log10(p)",
      fill_name = "NES",
      x_text_angle = 45,
      title = pos_title,
      xlab = "NES",
      ylab = ""
    )
    # Calculate width based on pathway label length
    max_label_length <- max(nchar(as.character(pos_pathways$Description)), na.rm = TRUE)
    width <- max(8, 5 + max_label_length * 0.08)
    height <- max(6, 3 + nrow(pos_pathways) * 0.25)
  } else if (nrow(neg_pathways) > 0) {
    # 只有negative pathways
    neg_title <- if (omics_type == "Mutation") {
      "Negative Terms (Mutation vs WildType)"
    } else {
      "Negative Terms (Negative Correlation)"
    }

    neg_pathways$neg_log10_p <- -log10(pmax(neg_pathways$pvalue, .Machine$double.xmin))
    plot <- DotPlot(
      data = neg_pathways,
      x = "NES",
      y = "Description",
      size_by = "neg_log10_p",
      fill_by = "NES",
      size_name = "-log10(p)",
      fill_name = "NES",
      x_text_angle = 45,
      title = neg_title,
      xlab = "NES",
      ylab = ""
    )
    # Calculate width based on pathway label length
    max_label_length <- max(nchar(as.character(neg_pathways$Description)), na.rm = TRUE)
    width <- max(8, 5 + max_label_length * 0.08)
    height <- max(6, 3 + nrow(neg_pathways) * 0.25)
  } else {
    stop("No significant pathways found for GSEA plot")
  }

  # GSEA Paired: 使用属性附加width和height
  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


#' DotPlot Paired (Internal)
#' @keywords internal
.plot_dotplot_paired <- function(all_stats, analysis_type = "DEA", cancer_types,
                                 genome_omics = "Protein", is_mutation = FALSE,
                                 use_mean = FALSE, feature_list = NULL, method = "pearson", top_n = 30) {
  library(dplyr)
  library(rlang)

  if (analysis_type == "DEA") {
    metric <- "logFC"
    p_col <- "pvalue"
    fill_name <- "logFC"
  } else {
    # Correlation: check which column name exists
    if ("correlation" %in% colnames(all_stats)) {
      metric <- "correlation"
      p_col <- "p_value"
    } else if ("r" %in% colnames(all_stats)) {
      metric <- "r"
      p_col <- if ("p_value" %in% colnames(all_stats)) "p_value" else "pvalue"
    } else {
      stop("No correlation column found in stats")
    }
    fill_name <- "Correlation"
  }

  # Standardize gene column name
  if (!"Gene" %in% colnames(all_stats) && "gene" %in% colnames(all_stats)) {
    all_stats$Gene <- all_stats$gene
  }

  # 步骤1: 每个变量分别选top N上调和top N下调基因
  # 注意：这里获取的是基因名列表，分开存储positive和negative
  top_up_genes <- all_stats %>%
    filter(!!rlang::sym(metric) > 0) %>%
    group_by(var_name) %>%
    arrange(!!rlang::sym(p_col)) %>%
    slice(1:top_n) %>%
    pull(Gene)

  top_down_genes <- all_stats %>%
    filter(!!rlang::sym(metric) < 0) %>%
    group_by(var_name) %>%
    arrange(!!rlang::sym(p_col)) %>%
    slice(1:top_n) %>%
    pull(Gene)

  # 步骤2: 分别去重，得到positive基因列表和negative基因列表
  up_gene_list <- unique(top_up_genes)
  down_gene_list <- unique(top_down_genes)

  # 步骤3: 分别从完整统计结果中提取（不过滤正负值！）
  # 对于up_gene_list中的基因，提取它们在所有变量中的所有数据（包括正值和负值）
  up_data <- all_stats %>%
    filter(Gene %in% up_gene_list)

  # 对于down_gene_list中的基因，提取它们在所有变量中的所有数据（包括正值和负值）
  down_data <- all_stats %>%
    filter(Gene %in% down_gene_list)

  # Standardize metric column name for plotting
  if (analysis_type == "Correlation") {
    if ("r" %in% colnames(up_data) && !"correlation" %in% colnames(up_data)) {
      up_data$correlation <- up_data$r
      down_data$correlation <- down_data$r
    }
    metric <- "correlation"
  } else {
    metric <- "logFC"
  }

  # Calculate -log10(p) for size
  p_col_use <- if (analysis_type == "DEA") "pvalue" else p_col
  up_data$neg_log10_p <- -log10(pmax(up_data[[p_col_use]], .Machine$double.xmin))
  down_data$neg_log10_p <- -log10(pmax(down_data[[p_col_use]], .Machine$double.xmin))

  # 使用RdBu palette同时显示正负值
  # 左图：显示来自positive基因列表的基因（但包含所有变量的正负值）
  plot_up <- DotPlot(
    data = up_data,
    x = "var_name",
    y = "Gene",
    size_by = "neg_log10_p",
    fill_by = metric,
    size_name = "-log10(p)",
    fill_name = fill_name,
    palette = "RdBu", # 红蓝渐变：正值红，负值蓝
    x_text_angle = 45,
    title = if (analysis_type == "DEA") "Up-regulated Genes" else "Positive Correlation Genes",
    xlab = "",
    ylab = ""
  )

  # 右图：显示来自negative基因列表的基因（但包含所有变量的正负值）
  plot_down <- DotPlot(
    data = down_data,
    x = "var_name",
    y = "Gene",
    size_by = "neg_log10_p",
    fill_by = metric,
    size_name = "-log10(p)",
    fill_name = fill_name,
    palette = "RdBu", # 红蓝渐变：正值红，负值蓝
    x_text_angle = 45,
    title = if (analysis_type == "DEA") "Down-regulated Genes" else "Negative Correlation Genes",
    xlab = "",
    ylab = ""
  )

  # Build title/subtitle/caption
  n_cancers <- length(unique(cancer_types))
  cancer_str <- if (n_cancers > 1) "Database" else cancer_types[1]
  is_signature <- use_mean && !is.null(feature_list) && length(feature_list) > 1

  title_text <- paste0("CPTAC-", cancer_str)

  # 修复：移除subtitle（用户要求）
  # 修复：移除Features信息（用户要求：图中已能看出来）

  caption_text <- if (is_signature) {
    method_str <- tools::toTitleCase(method)
    paste0("Method: ", method_str, " correlation with mean expression")
  } else {
    n_features <- length(unique(all_stats$var_name))
    paste0(
      "Method: ",
      if (analysis_type == "DEA") "Differential expression (DEA)" else paste0(tools::toTitleCase(method), " correlation")
    )
  }

  plot <- patchwork::wrap_plots(plot_up, plot_down, ncol = 2) +
    patchwork::plot_annotation(
      title = title_text,
      subtitle = NULL, # 修复：移除subtitle
      caption = caption_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
        plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray40")
      )
    )

  n_vars <- length(unique(all_stats$var_name))
  max_genes_up <- length(unique(up_data$Gene))
  max_genes_down <- length(unique(down_data$Gene))
  max_genes <- max(max_genes_up, max_genes_down)

  # Calculate width based on y-axis label length (user requirement!)
  all_gene_labels <- c(as.character(up_data$Gene), as.character(down_data$Gene))
  max_label_length <- if (length(all_gene_labels) > 0) max(nchar(all_gene_labels), na.rm = TRUE) else 10

  # Base width + x-axis variables + y-axis label length
  width <- max(10, 5 + n_vars * 0.8 + max_label_length * 0.08)
  height <- max(6, 3 + max_genes * 0.25)

  # 设置合理上限（防止图片过大不实用）
  MAX_WIDTH <- 50
  MAX_HEIGHT <- 100

  if (width > MAX_WIDTH || height > MAX_HEIGHT) {
    message(sprintf(
      "\n[Note] Large plot detected: %.1f × %.1f inches (%d genes × %d variables)",
      width, height, max_genes, n_vars
    ))
    if (height > MAX_HEIGHT) {
      message(sprintf(
        "       Height limited to %d inches for readability (was %.1f inches)",
        MAX_HEIGHT, height
      ))
      height <- MAX_HEIGHT
    }
    if (width > MAX_WIDTH) {
      width <- MAX_WIDTH
    }
  }

  # DotPlot Paired: 使用属性附加width和height
  attr(plot, "width") <- width
  attr(plot, "height") <- height

  # 返回筛选后的数据（用于stats）
  filtered_stats <- rbind(
    up_data %>% select(-neg_log10_p),
    down_data %>% select(-neg_log10_p)
  )
  attr(plot, "filtered_stats") <- filtered_stats

  return(plot)
}


#' GSEA Matrix DotPlot (Internal)
#' @keywords internal
.plot_gsea_matrix <- function(all_gsea_stats, enrich_type, GO_ont = NULL,
                              cancer_types, method = NULL,
                              use_mean = FALSE, feature_list = NULL, top_n = 10, omics_type = NULL) {
  library(dplyr)
  library(rlang)

  # 步骤1: 每个变量分别选top N正向和top N负向pathways（基于p-value）
  # 注意：这里获取的是pathway名称列表，分开存储positive和negative
  top_pos_pathways <- all_gsea_stats %>%
    filter(NES > 0) %>%
    group_by(var_name) %>%
    arrange(pvalue) %>%
    slice(1:top_n) %>%
    pull(Description)

  top_neg_pathways <- all_gsea_stats %>%
    filter(NES < 0) %>%
    group_by(var_name) %>%
    arrange(pvalue) %>%
    slice(1:top_n) %>%
    pull(Description)

  # 步骤2: 分别去重，得到positive pathway列表和negative pathway列表
  pos_pathway_list <- unique(top_pos_pathways)
  neg_pathway_list <- unique(top_neg_pathways)

  # 步骤3: 分别从完整GSEA结果中提取（不过滤NES正负值！）
  # 对于pos_pathway_list中的pathways，提取它们在所有变量中的所有数据（包括正值和负值NES）
  pos_data <- all_gsea_stats %>%
    filter(Description %in% pos_pathway_list)

  # 对于neg_pathway_list中的pathways，提取它们在所有变量中的所有数据（包括正值和负值NES）
  neg_data <- all_gsea_stats %>%
    filter(Description %in% neg_pathway_list)

  # 创建pathway_label
  pos_data$pathway_label <- sapply(pos_data$Description, function(desc) {
    if (is.na(desc)) {
      return(desc)
    }
    if (nchar(desc) > 60) paste0(substr(desc, 1, 57), "...") else desc
  })

  neg_data$pathway_label <- sapply(neg_data$Description, function(desc) {
    if (is.na(desc)) {
      return(desc)
    }
    if (nchar(desc) > 60) paste0(substr(desc, 1, 57), "...") else desc
  })

  pos_data$neg_log10_p <- -log10(pmax(pos_data$qvalue, .Machine$double.xmin))
  neg_data$neg_log10_p <- -log10(pmax(neg_data$qvalue, .Machine$double.xmin))

  # 根据omics_type设置子图标题
  if (!is.null(omics_type) && omics_type == "Mutation") {
    pos_title <- "Positive Terms (Mutation vs WildType)"
    neg_title <- "Negative Terms (Mutation vs WildType)"
  } else {
    pos_title <- "Positive Terms (Positive Correlation)"
    neg_title <- "Negative Terms (Negative Correlation)"
  }

  # 使用RdBu palette同时显示正负值
  # 左图：显示来自positive pathway列表的pathways（但包含所有变量的正负NES）
  plot_pos <- DotPlot(
    data = pos_data,
    x = "var_name",
    y = "pathway_label",
    size_by = "neg_log10_p",
    fill_by = "NES",
    size_name = "-log10(qvalue)",
    fill_name = "NES",
    palette = "RdBu", # 红蓝渐变：正NES红，负NES蓝
    x_text_angle = 45,
    title = pos_title,
    xlab = "",
    ylab = ""
  )

  # 右图：显示来自negative pathway列表的pathways（但包含所有变量的正负NES）
  plot_neg <- DotPlot(
    data = neg_data,
    x = "var_name",
    y = "pathway_label",
    size_by = "neg_log10_p",
    fill_by = "NES",
    size_name = "-log10(qvalue)",
    fill_name = "NES",
    palette = "RdBu", # 红蓝渐变：正NES红，负NES蓝
    x_text_angle = 45,
    title = neg_title,
    xlab = "",
    ylab = ""
  )

  # Build title (不要caption)
  n_cancers <- length(unique(cancer_types))
  cancer_str <- if (n_cancers > 1) "Database" else cancer_types[1]

  title_text <- paste0(enrich_type, " Enrichment Matrix in CPTAC-", cancer_str)

  plot <- patchwork::wrap_plots(plot_pos, plot_neg, ncol = 2) +
    patchwork::plot_annotation(
      title = title_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black")
      )
    )

  n_vars <- length(unique(all_gsea_stats$var_name))
  max_pathways_pos <- if (nrow(pos_data) > 0) length(unique(pos_data$pathway_label)) else 0
  max_pathways_neg <- if (nrow(neg_data) > 0) length(unique(neg_data$pathway_label)) else 0
  max_pathways <- max(max_pathways_pos, max_pathways_neg, na.rm = TRUE)

  # Calculate width based on y-axis pathway label length (user requirement!)
  all_pathway_labels <- c(as.character(pos_data$pathway_label), as.character(neg_data$pathway_label))
  max_label_length <- if (length(all_pathway_labels) > 0) max(nchar(all_pathway_labels), na.rm = TRUE) else 50

  # Base width + x-axis variables + y-axis label length
  # Increased base and multipliers for better readability
  width <- max(14, 8 + n_vars * 1.5 + max_label_length * 0.10, na.rm = TRUE)
  height <- max(6, 3 + max_pathways * 0.25, na.rm = TRUE)

  # Ensure valid dimensions
  width <- ifelse(is.finite(width) && width > 0, width, 12)
  height <- ifelse(is.finite(height) && height > 0, height, 8)

  # 设置合理上限（防止图片过大不实用）
  MAX_WIDTH <- 50
  MAX_HEIGHT <- 100

  if (width > MAX_WIDTH || height > MAX_HEIGHT) {
    message(sprintf(
      "\n[Note] Large plot detected: %.1f × %.1f inches (%d pathways × %d variables)",
      width, height, max_pathways, n_vars
    ))
    if (height > MAX_HEIGHT) {
      message(sprintf(
        "       Height limited to %d inches for readability (was %.1f inches)",
        MAX_HEIGHT, height
      ))
      height <- MAX_HEIGHT
    }
    if (width > MAX_WIDTH) {
      width <- MAX_WIDTH
    }
  }

  # GSEA Matrix: 使用属性附加width和height
  attr(plot, "width") <- width
  attr(plot, "height") <- height

  # 返回筛选后的数据（用于stats）
  filtered_stats <- rbind(
    pos_data %>% select(-neg_log10_p, -pathway_label),
    neg_data %>% select(-neg_log10_p, -pathway_label)
  )
  attr(plot, "filtered_stats") <- filtered_stats

  return(plot)
}


# ==============================================================================
# SURVIVAL PLOTS
# ==============================================================================

#' KM + Cox Combined Plot (Internal)
#' @keywords internal
.plot_km_cox_combined <- function(km_fit, cox_model_stats, data, time_col, event_col, group_col, var_name, omics_type, cancer_type, surv_type, var_col = NULL) {
  # Prepare data for both KM and Cox plots
  plot_data <- data
  plot_data$time <- plot_data[[time_col]]
  plot_data$event <- plot_data[[event_col]]
  plot_data$Group <- plot_data[[group_col]]

  # Get the feature column for Cox curve
  if (is.null(var_col)) {
    feature_cols <- setdiff(colnames(plot_data), c("time", "event", "Group", "cancer_type", time_col, event_col, group_col))
    if (length(feature_cols) == 0) {
      stop("No feature column found for Cox curve. Available columns: ", paste(colnames(plot_data), collapse = ", "))
    }
    var_col <- feature_cols[1]
  }

  # 1. KM plot using KMPlot
  km_plot <- KMPlot(
    data = plot_data,
    time = "time",
    status = "event",
    group_by = "Group",
    group_name = paste0(var_name, " (", omics_type, ")"),
    show_pval = TRUE,
    show_risk_table = FALSE,
    palcolor = c("#41A98E", "#ED6355"),
    title = paste0("CPTAC-", cancer_type, " (", surv_type, ")"),
    xlab = "Time (years)",
    ylab = "Survival Probability",
    legend.position = "bottom"
  )

  # 2. Cox Curve using CoxPlot
  cox_plot <- CoxPlot(
    data = plot_data,
    time = "time",
    event = "event",
    var = var_col,
    plot_type = "curve",
    show_cindex = TRUE,
    xlab = paste0(var_name, " (", omics_type, ")"),
    title = paste0("CPTAC-", cancer_type, " (", surv_type, ")")
  )

  # 3. Combine with patchwork
  combined_plot <- patchwork::wrap_plots(km_plot, cox_plot, ncol = 2, widths = c(1, 1.3))

  width <- 10 # Reduced to 10 inches (KM + Cox curve)
  height <- 5

  attr(combined_plot, "width") <- width
  attr(combined_plot, "height") <- height
  return(combined_plot)
}


#' Forest Plot (Internal)
#' @keywords internal
.plot_forest <- function(cox_stats, surv_type, cancer_type) {
  # 直接使用cox_stats绘制forest plot
  forest_data <- cox_stats
  forest_data$variable <- factor(forest_data$variable, levels = rev(forest_data$variable))

  plot <- ggplot2::ggplot(forest_data, ggplot2::aes(y = variable, x = hr)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = hr_lower, xmax = hr_upper),
      width = 0.2,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = ifelse(p_value < 0.05, "Significant", "Not significant")),
      size = 4
    ) +
    ggplot2::scale_color_manual(
      values = c("Significant" = "#ED6355", "Not significant" = "#41A98E"),
      name = ""
    ) +
    ggplot2::labs(
      title = paste0("CPTAC-", cancer_type, " (", surv_type, ")"),
      caption = "Method: Cox proportional hazards model",
      x = "Hazard Ratio (95% CI)",
      y = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray40"),
      axis.title = ggplot2::element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, color = "black"),
      axis.text.y = ggplot2::element_text(hjust = 1, size = 12, color = "black"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10)
    )

  n_vars <- nrow(forest_data)
  max_label_length <- max(nchar(as.character(forest_data$variable)), na.rm = TRUE)
  width <- max(5.5, 3.5 + max_label_length * 0.08)
  height <- max(4, 1 + n_vars * 0.4)

  attr(plot, "width") <- width
  attr(plot, "height") <- height
  return(plot)
}


#' Custom GSEA Dot Plot
#'
#' @keywords internal
.custom_gsea_dot <- function(res, type, show.term.num = 15, Select.P = "NP", cutoff.P = 0.05,
                             colors = c(
                               "#003c30", "#01665e", "#35978f", "#80cdc1", "#c7eae5",
                               "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005"
                             ),
                             y.label.position = "right", title = NULL, legend.position = "right",
                             theme.plot = ggplot2::theme_bw(base_rect_size = 1.5)) {
  library(dplyr)
  library(ggplot2)

  # Filter and sort
  if (Select.P == "NP") {
    if (type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")) {
      r <- res %>%
        filter(NES > 0) %>%
        mutate(sig = -log10(pvalue)) %>%
        arrange(desc(sig))
    } else {
      r <- res %>%
        filter(NES < 0) %>%
        mutate(sig = -log10(pvalue)) %>%
        arrange(desc(sig))
    }
    x.lab <- bquote(~ -Log[10] ~ italic("P-value"))
  } else if (Select.P == "FDR") {
    if (type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")) {
      r <- res %>%
        filter(NES > 0) %>%
        mutate(sig = -log10(p.adjust)) %>%
        arrange(desc(sig))
    } else {
      r <- res %>%
        filter(NES < 0) %>%
        mutate(sig = -log10(p.adjust)) %>%
        arrange(desc(sig))
    }
    x.lab <- bquote(~ -Log[10] ~ "FDR")
  }

  # Limit number of terms
  show.term.num <- ifelse(nrow(r) >= show.term.num, show.term.num, nrow(r))
  if (show.term.num == 0) {
    stop(paste0("No significant ", tolower(type), " term was enriched in res!"))
  }
  r <- r[seq_len(show.term.num), ]
  r <- r %>% mutate(Description = factor(Description, rev(Description)))

  # 修复：不使用abs(NES)，直接用原始NES
  # 对于负值，配合reverse colors自然实现越负越深
  is_positive <- type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")

  if (is_positive) {
    color.title <- "NES"
    fill_var <- "NES" # 正值直接用NES
  } else {
    color.title <- "NES"
    fill_var <- "NES" # 修复：负值也直接用NES，不用abs
  }

  # 创建plot
  p <- ggplot(r, aes(sig, Description, fill = .data[[fill_var]])) +
    geom_point(shape = 21, color = "black", size = 5) +
    labs(fill = color.title, x = x.lab, y = NULL, title = title) + # 移除size（没有映射）
    scale_fill_gradientn(colours = colors) + # colors已经是reverse的
    scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) +
    theme.plot +
    theme(
      axis.text.y = element_text(size = 13, colour = "black"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
      plot.title = element_text(size = 14, colour = "black", face = "bold"),
      legend.position = legend.position,
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 13, colour = "black", face = "bold"),
      legend.text = element_text(size = 11, colour = "black")
    )

  return(p)
}

#' Custom GSEA Dot Plot
#'
#' 修复：Negative类型不使用abs(NES)，直接用原始NES + reverse colors
#'
#' @keywords internal
.vis_gsea_dot_custom <- function(res,
                                 type = "Positive",
                                 show.term.num = 15,
                                 Select.P = "FDR",
                                 cutoff.P = 0.05,
                                 colors = c(
                                   "#003c30", "#01665e", "#35978f", "#80cdc1", "#c7eae5",
                                   "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005"
                                 ),
                                 add.bar.border = TRUE,
                                 bar.width = 0.6,
                                 y.label.position = "right",
                                 title = NULL,
                                 legend.position = "right",
                                 theme.plot = ggplot2::theme_bw(base_rect_size = 1.5)) {
  library(dplyr)
  library(ggplot2)

  # Filter and prepare data
  if (Select.P == "NP") {
    if (type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")) {
      r <- res %>%
        dplyr::filter(NES > 0) %>%
        dplyr::mutate(sig = -log10(pvalue)) %>%
        dplyr::arrange(dplyr::desc(sig))
    } else {
      r <- res %>%
        dplyr::filter(NES < 0) %>%
        dplyr::mutate(sig = -log10(pvalue)) %>%
        dplyr::arrange(dplyr::desc(sig))
    }
    x.lab <- bquote(~ -Log[10] ~ italic("P-value"))
  } else if (Select.P == "FDR") {
    if (type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")) {
      r <- res %>%
        dplyr::filter(NES > 0) %>%
        dplyr::mutate(sig = -log10(p.adjust)) %>%
        dplyr::arrange(dplyr::desc(sig))
    } else {
      r <- res %>%
        dplyr::filter(NES < 0) %>%
        dplyr::mutate(sig = -log10(p.adjust)) %>%
        dplyr::arrange(dplyr::desc(sig))
    }
    x.lab <- bquote(~ -Log[10] ~ "FDR")
  }

  show.term.num <- ifelse(nrow(r) >= show.term.num, show.term.num, nrow(r))

  if (show.term.num == 0) {
    stop(paste0("No significant ", tolower(type), " term was enriched in res!"))
  }

  r <- r[seq_len(show.term.num), ]
  r <- r %>% dplyr::mutate(Description = factor(Description, rev(Description)))

  # 关键修复：不使用abs(NES)
  is_negative <- !type %in% c("Positive", "positive", "pos", "p", "po", "P", "Po", "Pos")

  if (is_negative) {
    # 负值类型：直接用NES（负值）+ reverse colors
    color.title <- "NES" # 修复：不用abs(NES)
    fill_var <- "NES" # 修复：不用abs
    colors_use <- rev(colors) # reverse colors让越负越深
  } else {
    # 正值类型：正常
    color.title <- "NES"
    fill_var <- "NES"
    colors_use <- colors
  }

  p <- ggplot2::ggplot(r, ggplot2::aes(x = sig, y = Description, fill = .data[[fill_var]])) +
    ggplot2::geom_point(shape = 21, color = "black", size = 5) +
    ggplot2::labs(fill = color.title, x = x.lab, y = NULL, title = title) +
    ggplot2::scale_fill_gradientn(colours = colors_use) +
    ggplot2::scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) +
    theme.plot +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 13, colour = "black"),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 13, colour = "black", face = "bold"),
      plot.title = ggplot2::element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
      legend.position = legend.position,
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 13, colour = "black", face = "bold"),
      legend.text = ggplot2::element_text(size = 11, colour = "black")
    )

  return(p)
}
