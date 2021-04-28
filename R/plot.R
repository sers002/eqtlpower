#' Plot the Estimated eQTL Power Results
#'
#' Generate a simple visualisation of the estimated statistical power from
#' either the \code{\link{analyse_coef_results}} or the analytic functions.
#'
#' @param power.dt Result dataset from either \code{analyse_coef_results} or \code{analytic_calculation_fdr}.
#' @param by A variable to colour lines by.
#' @param title Plot title
#'
#' @return A ggplot of the power curves
#' @export
#'
#' @examples
#' \donttest{
#' data("associations")
#' data("ld")
#'
#' sa_coefs <- generate_coefficients(associations.dt, ld.dt)
#' sa_results <- analyse_coef_results(associations.dt, sa_coefs)
#'
#' plot_power_results(sa_results)
#' }
#'
#' bh_results.dt <- analytic_calculation_fdr(
#'   associations.dt[qtl_tag == TRUE],
#'   by = c("cis_flag")
#' )
#'
#' plot_power_results(bh_results.dt, by = "cis_flag")
plot_power_results <- function(
  power.dt,
  by = NULL,
  title = "Estimated Statistical Power for eQTL Study"
) {
  sample_size <-
    power <-
    NULL

  if (is.null(by)) {
    p <- ggplot2::ggplot(power.dt, ggplot2::aes(x = sample_size, y = power)) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  } else {
    p <- ggplot2::ggplot(power.dt, ggplot2::aes_string(x = "sample_size", y = "power", colour = by)) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  }

  p +
    ggplot2::labs(
      x = "Sample Size",
      y = "Power",
      title = title
    )
}
