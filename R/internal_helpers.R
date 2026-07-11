# Internal helper functions for MultipletR.
# None of these are exported; they are called by find_adaptive_thresholds_3t()
# and, through it, by detect_multiplets().

#' Test a distribution for bimodality (Hartigan dip test)
#' @param x Numeric vector.
#' @param alpha Significance level. Default 0.05.
#' @param min_cells_for_test Minimum values needed to run the test. Default 30.
#' @return A list: is_bimodal, p_value, dip_statistic, reason.
#' @keywords internal
test_bimodality <- function(x, alpha = 0.05, min_cells_for_test = 30) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) < min_cells_for_test) {
    return(list(
      is_bimodal = FALSE, p_value = NA_real_, dip_statistic = NA_real_,
      reason = paste0("Not enough cells (", length(x), "<", min_cells_for_test, ")")
    ))
  }
  tryCatch(
    {
      dip_result <- diptest::dip.test(x)
      list(
        is_bimodal    = dip_result$p.value < alpha,
        p_value       = dip_result$p.value,
        dip_statistic = as.numeric(dip_result$statistic),
        reason        = if (dip_result$p.value < alpha) "Bimodal detected" else "Unimodal"
      )
    },
    error = function(e) {
      list(
        is_bimodal = FALSE, p_value = NA_real_, dip_statistic = NA_real_,
        reason = paste0("Error: ", e$message)
      )
    }
  )
}

#' Find the primary mode of a distribution via kernel density
#' @param x Numeric vector.
#' @return The x-position of the density peak, or NA.
#' @keywords internal
find_mode <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) < 10) {
    return(NA_real_)
  }
  tryCatch(
    {
      d <- stats::density(x, n = 512)
      mode_value <- d$x[which.max(d$y)]
      if (is.null(mode_value) || length(mode_value) == 0) {
        return(NA_real_)
      }
      as.numeric(mode_value[1])
    },
    error = function(e) NA_real_
  )
}

#' Distribution-quality metrics for a set of selected cells
#'
#' Computes the number of cells, the percent overlap between the human and
#' mouse log-read distributions, the area overlap, the difference between the
#' distribution modes, and bimodality of each distribution.
#'
#' @param data Data frame with human_10x and mouse_10x columns.
#' @param min_cells_for_bimodality Minimum cells to run the bimodality test.
#' @return A named list of metrics.
#' @keywords internal
calculate_adaptive_metrics <- function(data, min_cells_for_bimodality = 30) {
  default_return <- list(
    n_cells = 0, pct_overlap = NA_real_, area_overlap = NA_real_,
    diff_modes = NA_real_, human_bimodal = FALSE, mouse_bimodal = FALSE,
    human_dip_pval = NA_real_, mouse_dip_pval = NA_real_,
    human_mode = NA_real_, mouse_mode = NA_real_
  )
  if (is.null(data) || nrow(data) < 10) {
    if (!is.null(data)) default_return$n_cells <- nrow(data)
    return(default_return)
  }
  human_reads <- data$human_10x
  mouse_reads <- data$mouse_10x
  if (is.null(human_reads) || is.null(mouse_reads)) {
    default_return$n_cells <- nrow(data)
    return(default_return)
  }
  human_log <- log(as.numeric(human_reads) + 1)
  mouse_log <- log(as.numeric(mouse_reads) + 1)
  human_log <- human_log[is.finite(human_log)]
  mouse_log <- mouse_log[is.finite(mouse_log)]
  n_cells <- nrow(data)
  if (length(human_log) < 10 || length(mouse_log) < 10) {
    default_return$n_cells <- n_cells
    return(default_return)
  }

  bimodal_human <- test_bimodality(human_log, min_cells_for_test = min_cells_for_bimodality)
  bimodal_mouse <- test_bimodality(mouse_log, min_cells_for_test = min_cells_for_bimodality)

  pct_overlap <- NA_real_
  area_overlap <- NA_real_
  mode_human <- NA_real_
  mode_mouse <- NA_real_
  diff_modes <- NA_real_
  tryCatch(
    {
      x_all <- c(human_log, mouse_log)
      range_min <- min(x_all, na.rm = TRUE) - 0.1 * diff(range(x_all))
      range_max <- max(x_all, na.rm = TRUE) + 0.1 * diff(range(x_all))
      if (is.finite(range_min) && is.finite(range_max) && range_max > range_min) {
        d_human <- stats::density(human_log, from = range_min, to = range_max, n = 2000)
        d_mouse <- stats::density(mouse_log, from = range_min, to = range_max, n = 2000)
        dx <- diff(d_human$x)
        area_human <- sum((d_human$y[-1] + d_human$y[-length(d_human$y)]) / 2 * dx)
        area_mouse <- sum((d_mouse$y[-1] + d_mouse$y[-length(d_mouse$y)]) / 2 * dx)
        overlap_y <- pmin(d_human$y, d_mouse$y)
        overlap_area <- sum((overlap_y[-1] + overlap_y[-length(overlap_y)]) / 2 * dx)
        area_overlap <- overlap_area
        pct_human <- overlap_area / area_human
        pct_mouse <- overlap_area / area_mouse
        pct_overlap <- ((pct_human + pct_mouse) / 2) * 100

        # Mode from the SAME ranged density (matches the manuscript)
        mode_human <- d_human$x[which.max(d_human$y)]
        mode_mouse <- d_mouse$x[which.max(d_mouse$y)]
        diff_modes <- abs(mode_human - mode_mouse)
      }
    },
    error = function(e) {}
  )

  list(
    n_cells = n_cells, pct_overlap = pct_overlap, area_overlap = area_overlap,
    diff_modes = diff_modes, human_bimodal = bimodal_human$is_bimodal,
    mouse_bimodal = bimodal_mouse$is_bimodal,
    human_dip_pval = bimodal_human$p_value, mouse_dip_pval = bimodal_mouse$p_value,
    human_mode = mode_human, mouse_mode = mode_mouse
  )
}

#' Initialize the conservative starting thresholds
#' @param data Data frame with a total_reads column.
#' @param percent_lower Lower percent-mouse threshold (T2). Default 30.
#' @param percent_upper Upper percent-mouse threshold (T1). Default 70.
#' @param reads_percentile Lower total-reads threshold as a percentile. Default 25.
#' @return A list of the initial thresholds.
#' @keywords internal
initialize_thresholds <- function(data,
                                  percent_lower = 30,
                                  percent_upper = 70,
                                  reads_percentile = 25) {
  th <- list(
    t1_mouse_pct      = percent_upper,
    t2_human_pct      = percent_lower,
    t3_lower_reads    = as.numeric(stats::quantile(data$total_reads, reads_percentile / 100)),
    upper_reads_fixed = max(data$total_reads, na.rm = TRUE)
  )
  # Compatibility bridge (used by other internal code)
  th$t4_right_pct <- th$t1_mouse_pct
  th$t3_left_pct <- th$t2_human_pct
  th$t1_lower_reads <- th$t3_lower_reads
  th$t2_upper_reads <- th$upper_reads_fixed
  th
}

#' Select cells within the current threshold rectangle
#' @param data Data frame with total_reads and pct_mouse_10x columns.
#' @param thresholds A thresholds list from initialize_thresholds().
#' @return The subset of rows inside the rectangle.
#' @keywords internal
select_cells_in_rectangle <- function(data, thresholds) {
  data[
    data$total_reads >= thresholds$t3_lower_reads &
      data$total_reads <= thresholds$upper_reads_fixed &
      data$pct_mouse_10x >= thresholds$t2_human_pct &
      data$pct_mouse_10x <= thresholds$t1_mouse_pct,
  ]
}

#' Expand one threshold by a single step
#' @param thresholds Current thresholds list.
#' @param which_threshold 1 = T1 (mouse% up), 2 = T2 (human% down), 3 = T3 (reads down).
#' @param step_size Step magnitude.
#' @param min_reads Minimum total reads (floor for T3).
#' @param max_reads Maximum total reads (unused; upper boundary is fixed).
#' @return The updated thresholds list.
#' @keywords internal
expand_threshold <- function(thresholds, which_threshold, step_size,
                             min_reads, max_reads) {
  new_thresholds <- thresholds
  if (which_threshold == 1) {
    new_thresholds$t1_mouse_pct <- min(100, thresholds$t1_mouse_pct + step_size)
  } else if (which_threshold == 2) {
    new_thresholds$t2_human_pct <- max(0, thresholds$t2_human_pct - step_size)
  } else if (which_threshold == 3) {
    new_thresholds$t3_lower_reads <- max(min_reads, thresholds$t3_lower_reads - step_size)
  }
  new_thresholds
}

#' Decide whether a threshold should stop expanding
#' @param current_metrics Metrics for the candidate selection.
#' @param previous_metrics Metrics for the last accepted selection.
#' @param overlap_drop_threshold Stop if overlap drops more than this from the previous step.
#' @param mode_diff_increase_threshold Stop if the mode difference increases more than this.
#' @return A list: stop (logical) and reason (character).
#' @keywords internal
should_stop_expanding <- function(current_metrics, previous_metrics,
                                  overlap_drop_threshold = 10,
                                  mode_diff_increase_threshold = 0.9) {
  if (is.na(current_metrics$n_cells) || current_metrics$n_cells < 30) {
    return(list(stop = FALSE, reason = "OK - need more cells"))
  }
  if (isTRUE(current_metrics$human_bimodal)) {
    return(list(stop = TRUE, reason = "HUMAN became BIMODAL (dip p < 0.05)"))
  }
  if (isTRUE(current_metrics$mouse_bimodal)) {
    return(list(stop = TRUE, reason = "MOUSE became BIMODAL (dip p < 0.05)"))
  }
  if (!is.null(previous_metrics) &&
    !is.na(previous_metrics$pct_overlap) &&
    !is.na(current_metrics$pct_overlap)) {
    overlap_drop <- previous_metrics$pct_overlap - current_metrics$pct_overlap
    if (overlap_drop > overlap_drop_threshold) {
      return(list(stop = TRUE, reason = paste0(
        "Overlap dropped ", round(overlap_drop, 1),
        "% from previous"
      )))
    }
  }
  if (!is.null(previous_metrics) &&
    !is.na(previous_metrics$diff_modes) &&
    !is.na(current_metrics$diff_modes)) {
    mode_diff_increase <- current_metrics$diff_modes - previous_metrics$diff_modes
    if (mode_diff_increase > mode_diff_increase_threshold) {
      return(list(stop = TRUE, reason = paste0(
        "Mode diff increased ", round(mode_diff_increase, 2),
        " from previous"
      )))
    }
  }
  list(stop = FALSE, reason = "OK")
}
