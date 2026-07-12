#' Adaptive three-threshold multiplet detection (internal engine)
#'
#' Round-robin expansion of three thresholds (upper percent-mouse, lower
#' percent-mouse, lower total reads) from a conservative starting region until
#' the human and mouse read distributions begin to separate. The upper reads
#' boundary is fixed at the maximum. This is the internal engine called by
#' \code{\link{detect_multiplets}}.
#'
#' @param data Data frame with columns total_reads, pct_mouse_10x, human_10x,
#'   mouse_10x.
#' @param percent_lower Lower percent-mouse starting threshold (T2). Default 30.
#' @param percent_upper Upper percent-mouse starting threshold (T1). Default 70.
#' @param reads_percentile Lower total-reads starting threshold as a percentile.
#'   Default 25.
#' @param step_size_reads Expansion step for the reads threshold. Default 500.
#' @param step_size_pct Expansion step for the percent thresholds. Default 1.
#' @param overlap_drop_threshold Stop when overlap drops more than this percent
#'   below the running peak. Default 10.
#' @param mode_diff_increase_threshold Stop when the mode difference increases
#'   more than this. Default 0.9.
#' @param overlap_patience Consecutive below-peak steps tolerated before
#'   stopping. Default 3.
#' @param min_cells_for_bimodality Minimum cells to run the bimodality test.
#'   Default 30.
#' @param verbose Print progress. Default FALSE.
#'
#' @return A list with the final and initial thresholds, stop reasons,
#'   iteration counts, initial and final metrics, the detected multiplets, the
#'   expansion history, and the total number of rounds.
#' @keywords internal
find_adaptive_thresholds_3t <- function(data,
                                        percent_lower = 30,
                                        percent_upper = 70,
                                        reads_percentile = 25,
                                        step_size_reads = 500,
                                        step_size_pct = 1,
                                        overlap_drop_threshold = 10,
                                        mode_diff_increase_threshold = 0.9,
                                        overlap_patience = 3,
                                        min_cells_for_bimodality = 30,
                                        verbose = FALSE) {
  # verbose output helper: cat-style formatting, routed through message()
  .vcat <- function(..., sep = " ") {
    message(paste(..., sep = sep), appendLF = FALSE)
  }
  if (verbose) {
    .vcat("===================================================================\n")
    .vcat("ADAPTIVE 3-THRESHOLD ALGORITHM\n")
    .vcat("T1 (mouse %), T2 (human %), T3 (total reads)\n")
    .vcat("Upper reads boundary fixed at max\n")
    .vcat("Stopping: bimodality + per-step drop + peak drop (patience) + mode diff\n")
    .vcat("===================================================================\n\n")
  }
  max_reads <- max(data$total_reads, na.rm = TRUE)
  min_reads <- min(data$total_reads, na.rm = TRUE)
  if (verbose) {
    .vcat("PARAMETERS:\n")
    .vcat("  step_size_reads:", step_size_reads, "\n")
    .vcat("  step_size_pct:", step_size_pct, "%\n")
    .vcat("  overlap_drop_threshold:", overlap_drop_threshold, "%\n")
    .vcat("  overlap_patience:", overlap_patience, "steps\n")
    .vcat("  mode_diff_increase_threshold:", mode_diff_increase_threshold, "\n")
    .vcat("  min_cells_for_bimodality:", min_cells_for_bimodality, "\n\n")
    .vcat("DATA RANGE:\n")
    .vcat("  Min reads: ", min_reads, "\n")
    .vcat("  Max reads: ", max_reads, " (upper boundary fixed here)\n\n")
  }
  thresholds <- initialize_thresholds(data, percent_lower, percent_upper, reads_percentile)
  initial_thresholds_copy <- thresholds
  if (verbose) {
    .vcat("INITIAL THRESHOLDS:\n")
    .vcat("  T1 (mouse % upper):  ", thresholds$t1_mouse_pct, "% (TO BE OPTIMIZED)\n")
    .vcat("  T2 (human % lower):  ", thresholds$t2_human_pct, "% (TO BE OPTIMIZED)\n")
    .vcat("  T3 (total reads):    ", round(thresholds$t3_lower_reads), " (Q1 - TO BE OPTIMIZED)\n")
    .vcat("  Upper reads boundary: ", round(thresholds$upper_reads_fixed), " (FIXED)\n\n")
  }
  initial_cells <- select_cells_in_rectangle(data, thresholds)
  initial_metrics <- calculate_adaptive_metrics(initial_cells, min_cells_for_bimodality)
  if (verbose) {
    .vcat("INITIAL METRICS:\n")
    .vcat("  N Cells:       ", initial_metrics$n_cells, "\n")
    .vcat("  % Overlap:     ", round(initial_metrics$pct_overlap, 2), "%\n")
    .vcat("  Area Overlap:  ", round(initial_metrics$area_overlap, 4), "\n")
    .vcat("  Diff in Modes: ", round(initial_metrics$diff_modes, 4), "\n")
    .vcat("  Human Bimodal: ", initial_metrics$human_bimodal, "\n")
    .vcat("  Mouse Bimodal: ", initial_metrics$mouse_bimodal, "\n\n")
  }
  can_expand <- c(TRUE, TRUE, TRUE)
  stop_reasons <- c("Still expanding", "Still expanding", "Still expanding")
  threshold_names <- c("T1 up", "T2 down", "T3 down")
  iterations_count <- c(0, 0, 0)
  history <- data.frame()
  previous_metrics <- initial_metrics
  peak_overlap <- initial_metrics$pct_overlap
  below_peak_count <- 0
  if (verbose) {
    .vcat("----------------------------------------------------------------\n")
    .vcat("STARTING ROUND-ROBIN EXPANSION (T1, T2, T3)\n")
    .vcat("----------------------------------------------------------------\n\n")
  }
  round_num <- 0
  while (any(can_expand)) {
    round_num <- round_num + 1
    if (round_num > 1000) {
      if (verbose) .vcat("\n*** Safety limit reached (1000 rounds) ***\n")
      for (t in seq_len(3)) {
        if (can_expand[t]) stop_reasons[t] <- "Safety limit (1000 rounds)"
      }
      break
    }
    if (verbose && round_num %% 10 == 1) {
      .vcat("--- Round", round_num, "---\n")
    }
    for (t in seq_len(3)) {
      if (!can_expand[t]) next
      iterations_count[t] <- iterations_count[t] + 1
      step <- if (t == 3) step_size_reads else step_size_pct
      new_thresholds <- expand_threshold(thresholds, t, step, min_reads, max_reads)
      boundary_reached <- FALSE
      boundary_type <- ""
      if (t == 1 && new_thresholds$t1_mouse_pct >= 100) {
        boundary_reached <- TRUE
        boundary_type <- "data limit (100% mouse)"
      }
      if (t == 2 && new_thresholds$t2_human_pct <= 0) {
        boundary_reached <- TRUE
        boundary_type <- "data limit (0% mouse)"
      }
      if (t == 3 && new_thresholds$t3_lower_reads <= min_reads) {
        boundary_reached <- TRUE
        boundary_type <- paste0("data limit (min reads = ", round(min_reads), ")")
      }
      if (boundary_reached) {
        can_expand[t] <- FALSE
        b_cells <- select_cells_in_rectangle(data, new_thresholds)
        b_metrics <- calculate_adaptive_metrics(b_cells, min_cells_for_bimodality)
        if (isTRUE(b_metrics$human_bimodal) || isTRUE(b_metrics$mouse_bimodal)) {
          bm <- if (isTRUE(b_metrics$human_bimodal) && isTRUE(b_metrics$mouse_bimodal)) {
            "BOTH"
          } else if (isTRUE(b_metrics$human_bimodal)) "HUMAN" else "MOUSE"
          stop_reasons[t] <- paste0("Boundary rejected - would become bimodal (", bm, ")")
          if (verbose) {
            kept_val <- switch(t,
              thresholds$t1_mouse_pct,
              thresholds$t2_human_pct,
              thresholds$t3_lower_reads
            )
            .vcat("  ", threshold_names[t], " -> STOPPED at ", kept_val,
              " (boundary would be bimodal: ", bm, ", after ", iterations_count[t], " steps)\n",
              sep = ""
            )
          }
        } else {
          stop_reasons[t] <- paste0("Data boundary: ", boundary_type)
          thresholds <- new_thresholds
          if (!is.na(b_metrics$pct_overlap)) {
            peak_overlap <- max(peak_overlap, b_metrics$pct_overlap, na.rm = TRUE)
          }
          if (verbose) {
            final_val <- switch(t,
              new_thresholds$t1_mouse_pct,
              new_thresholds$t2_human_pct,
              new_thresholds$t3_lower_reads
            )
            .vcat("  ", threshold_names[t], " -> STOPPED at ", final_val,
              " (", boundary_type, " after ", iterations_count[t], " steps)\n",
              sep = ""
            )
          }
        }
        next
      }
      new_cells <- select_cells_in_rectangle(data, new_thresholds)
      new_metrics <- calculate_adaptive_metrics(new_cells, min_cells_for_bimodality)
      stop_check <- should_stop_expanding(
        new_metrics, previous_metrics,
        overlap_drop_threshold,
        mode_diff_increase_threshold
      )
      if (!stop_check$stop &&
        (isTRUE(new_metrics$human_bimodal) || isTRUE(new_metrics$mouse_bimodal))) {
        bm <- if (isTRUE(new_metrics$human_bimodal) && isTRUE(new_metrics$mouse_bimodal)) {
          "BOTH"
        } else if (isTRUE(new_metrics$human_bimodal)) "HUMAN" else "MOUSE"
        stop_check$stop <- TRUE
        stop_check$reason <- paste0("Combined selection became bimodal (", bm, ")")
      }
      if (!stop_check$stop && !is.na(new_metrics$pct_overlap) &&
        (peak_overlap - new_metrics$pct_overlap) > overlap_drop_threshold) {
        below_peak_count <- below_peak_count + 1
        if (below_peak_count >= overlap_patience) {
          stop_check$stop <- TRUE
          stop_check$reason <- paste0(
            "Overlap ",
            round(peak_overlap - new_metrics$pct_overlap, 1),
            "% below peak (", round(peak_overlap, 1), "%) for ",
            overlap_patience, " steps"
          )
        }
      } else if (!is.na(new_metrics$pct_overlap)) {
        below_peak_count <- 0
      }
      current_val <- switch(t,
        new_thresholds$t1_mouse_pct,
        new_thresholds$t2_human_pct,
        new_thresholds$t3_lower_reads
      )
      hist_row <- data.frame(
        round = round_num, threshold = threshold_names[t], threshold_num = t,
        value = current_val, n_cells = new_metrics$n_cells,
        pct_overlap = new_metrics$pct_overlap,
        area_overlap = new_metrics$area_overlap,
        diff_modes = new_metrics$diff_modes,
        human_bimodal = new_metrics$human_bimodal, mouse_bimodal = new_metrics$mouse_bimodal,
        human_dip_pval = new_metrics$human_dip_pval, mouse_dip_pval = new_metrics$mouse_dip_pval,
        stop_reason = stop_check$reason, stringsAsFactors = FALSE
      )
      history <- rbind(history, hist_row)
      if (stop_check$stop) {
        can_expand[t] <- FALSE
        stop_reasons[t] <- stop_check$reason
        last_good_val <- switch(t,
          thresholds$t1_mouse_pct,
          thresholds$t2_human_pct,
          thresholds$t3_lower_reads
        )
        if (verbose) {
          .vcat("  ", threshold_names[t], " -> STOPPED at ", last_good_val,
            " (", stop_check$reason, " after ", iterations_count[t], " steps)\n",
            sep = ""
          )
        }
      } else {
        thresholds <- new_thresholds
        previous_metrics <- new_metrics
        if (!is.na(new_metrics$pct_overlap)) {
          peak_overlap <- max(peak_overlap, new_metrics$pct_overlap, na.rm = TRUE)
        }
      }
    }
  }
  if (verbose) .vcat("\n*** All thresholds stopped at round", round_num, "***\n")
  final_cells <- select_cells_in_rectangle(data, thresholds)
  final_metrics <- calculate_adaptive_metrics(final_cells, min_cells_for_bimodality)
  if (verbose) {
    .vcat("\n")
    .vcat("===================================================================\n")
    .vcat("FINAL OPTIMAL THRESHOLDS (3-THRESHOLD)\n")
    .vcat("===================================================================\n\n")
    .vcat("THRESHOLD DETAILS:\n")
    .vcat(sprintf(
      "  T1 (mouse %%): %d%% -> %d%% (%d steps) [%s]\n",
      initial_thresholds_copy$t1_mouse_pct, thresholds$t1_mouse_pct,
      iterations_count[1], stop_reasons[1]
    ))
    .vcat(sprintf(
      "  T2 (human %%): %d%% -> %d%% (%d steps) [%s]\n",
      initial_thresholds_copy$t2_human_pct, thresholds$t2_human_pct,
      iterations_count[2], stop_reasons[2]
    ))
    .vcat(sprintf(
      "  T3 (reads): %d -> %d (%d steps) [%s]\n",
      round(initial_thresholds_copy$t3_lower_reads), round(thresholds$t3_lower_reads),
      iterations_count[3], stop_reasons[3]
    ))
    .vcat(sprintf(
      "  Upper reads boundary: %d (FIXED)\n",
      round(thresholds$upper_reads_fixed)
    ))
    .vcat("\nFINAL METRICS:\n")
    .vcat("  N Multiplets:  ", final_metrics$n_cells, "\n")
    .vcat("  % Overlap:     ", round(final_metrics$pct_overlap, 2), "%\n")
    .vcat("  Area Overlap:  ", round(final_metrics$area_overlap, 4), "\n")
    .vcat("  Diff in Modes: ", round(final_metrics$diff_modes, 4), "\n")
    .vcat("  Human Bimodal: ", final_metrics$human_bimodal, "\n")
    .vcat("  Mouse Bimodal: ", final_metrics$mouse_bimodal, "\n\n")
  }
  return(list(
    thresholds = thresholds, initial_thresholds = initial_thresholds_copy,
    stop_reasons = stop_reasons, iterations_count = iterations_count,
    initial_metrics = initial_metrics, final_metrics = final_metrics,
    detected_multiplets = final_cells, history = history, total_rounds = round_num
  ))
}
