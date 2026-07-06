utils::globalVariables(c("our_call", "total_reads", "pct_mouse_10x",
                         "call_10x", "mouse_10x", "human_10x"))
#' Percent plots: total reads vs percent mouse, 10x vs our classification
#'
#' Draws two side-by-side scatter plots of total reads (x) against percent
#' mouse content (y): one colored by the original 10x classification, one by
#' our multiplet classification.
#'
#' @param cc Data frame with columns total_reads, pct_mouse_10x, call_10x,
#'   and our_classification.
#' @return A patchwork object (two panels side by side).
#' @importFrom patchwork wrap_plots
#' @keywords internal
.plot_percent <- function(cc) {
  # left: 10x classification
  p_10x <- ggplot2::ggplot(cc,
      ggplot2::aes(x = total_reads, y = pct_mouse_10x, color = call_10x)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.0) +
    ggplot2::scale_color_manual(values = c(
      "Human"     = "blue",
      "Mouse"     = "green3",
      "Multiplet" = "purple")) +
    ggplot2::labs(x = "Total reads", y = "Percent mouse (%)",
                  color = "10X call", title = "10X classification") +
    ggplot2::coord_cartesian(ylim = c(-5, 105)) +
    ggplot2::theme_minimal(base_size = 12)

  # right: our classification (3-way, matching the main analysis figures)
  cc$our_call <- ifelse(cc$our_classification == "Multiplet", "Multiplet (Ours)",
                        ifelse(cc$pct_mouse_10x < 50, "Human (not multiplet)",
                               "Mouse (not multiplet)"))
  p_our <- ggplot2::ggplot(cc,
                           ggplot2::aes(x = total_reads, y = pct_mouse_10x, color = our_call)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.0) +
    ggplot2::scale_color_manual(values = c(
      "Multiplet (Ours)"      = "red",
      "Human (not multiplet)" = "lightblue",
      "Mouse (not multiplet)" = "lightgreen")) +
    ggplot2::labs(x = "Total reads", y = "Percent mouse (%)",
                  color = "Our call", title = "Our classification") +
    ggplot2::coord_cartesian(ylim = c(-5, 105)) +
    ggplot2::theme_minimal(base_size = 12)

  p_10x + p_our
}

#' Total-reads plots: mouse reads vs human reads, 10x vs our classification
#'
#' Draws two side-by-side scatter plots of mouse reads (x) against human reads
#' (y): one colored by the original 10x classification, one by our multiplet
#' classification.
#'
#' @param cc Data frame with columns human_10x, mouse_10x, call_10x,
#'   and our_classification.
#' @return A patchwork object (two panels side by side).
#' @keywords internal
.plot_total_reads <- function(cc) {
  # left: 10x classification
  p_10x <- ggplot2::ggplot(cc,
      ggplot2::aes(x = mouse_10x, y = human_10x, color = call_10x)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.0) +
    ggplot2::scale_color_manual(values = c(
      "Human"     = "blue",
      "Mouse"     = "green3",
      "Multiplet" = "purple")) +
    ggplot2::labs(x = "Mouse reads", y = "Human reads",
                  color = "10X call", title = "10X classification") +
    ggplot2::theme_minimal(base_size = 12)

  # right: our classification (3-way, matching the main analysis figures)
  cc$our_call <- ifelse(cc$our_classification == "Multiplet", "Multiplet (Ours)",
                        ifelse(cc$pct_mouse_10x < 50, "Human (not multiplet)",
                               "Mouse (not multiplet)"))
  p_our <- ggplot2::ggplot(cc,
                           ggplot2::aes(x = mouse_10x, y = human_10x, color = our_call)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.0) +
    ggplot2::scale_color_manual(values = c(
      "Multiplet (Ours)"      = "red",
      "Human (not multiplet)" = "lightblue",
      "Mouse (not multiplet)" = "lightgreen")) +
    ggplot2::labs(x = "Mouse reads", y = "Human reads",
                  color = "Our call", title = "Our classification") +
    ggplot2::theme_minimal(base_size = 12)
  patchwork::wrap_plots(p_10x, p_our)
}
