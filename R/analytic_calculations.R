#' Calculate Whether \eqn{p}-value Threshold Results in a Given FDR
#'
#' This function calculates the FDR for a given set of \eqn{p}-values if the
#' threshold is set at \code{pval_threshold}. It achieves this by solving
#' \eqn{(M \pi_0 \times pval_threshold)/(\sum(power for each eQTL))}. This
#' function is used within the \code{\link{analytic_calculation_fdr}} function.
#'
#' @param pval_threshold \eqn{p}-value threshold to check against
#' @param M Number of comparisons being made
#' @param q Desired level of FDR
#' @param ncp Non-centrality parameter
#' @param pi0 \eqn{pi_0} parameter for pFDR
#'
#' @return The value of the equation given in the description.
#' @export
#'
#' @examples
#' data("associations")
#' test_ncp <- 515 * 2 * associations.dt$maf * (1 - associations.dt$maf) * associations.dt$beta^2
#' fdr_p_threshold(0.05, ncp = test_ncp)
fdr_p_threshold <- function(pval_threshold, M = length(ncp), q = 0.05, ncp, pi0 = 1) {
  ## Chisquare critical value corresponding to p-value threshold
  chisq_threshold = stats::qchisq(1 - pval_threshold, 1)

  ## Power for each SNP using its NCP
  curr_power <- 1-stats::pchisq(chisq_threshold, 1, ncp = ncp)

  if (M != length(ncp)) {
    ((M * pi0 * pval_threshold) / (sum(curr_power) + M * pi0 * pval_threshold )) - q
  } else {
    ((M * pi0 * pval_threshold) / sum(curr_power)) - q
  }
}

#' Analytic Power Calculation for Benjamini-Hochberg FDR
#'
#' This function calculates the long-run statistical power to detect eQTLs when
#' using the Benjamini-Hochberg procedure to determine significance.
#'
#' @param associations A dataset of association data. It must contain the
#' columns described below (\code{pid}, \code{snp}, \code{beta}, \code{maf})
#' @param var_y The variance of the response
#' @param sample_sizes A set of sample sizes to test
#' @param n_tests The number of tests in the full dataset
#' @param q FDR threshold
#' @param by A set of variable names to aggregate power over
#' @param sample_values Should omitted variants (\code{n_tests - nrow(associations)})
#' be sampled from \code{associations}
#' @param sample_prop If \code{sample_values = TRUE}, the proportion of the
#' sampled values that should be sampled from \code{associations}. The
#' remaining proportion of sampled values (\code{1 - sample_prop}) will be null
#' values.
#' @param verbose If \code{TRUE}, messages will be output during the course of
#' the function.
#'
#' @return The estimated statistical power for each gene-SNP association at an
#' FDR of \code{q}.
#' @export
#'
#' @examples
#' data("associations")
#' bh_results.dt <- analytic_calculation_fdr(associations.dt)
analytic_calculation_fdr <- function(
  associations,
  var_y = 1,
  sample_sizes = c(50, 100, 250, 500),
  n_tests = nrow(associations),
  q = 0.05,
  by = NULL,
  sample_values = FALSE,
  sample_prop = 0,
  verbose = TRUE
) {
  maf <-
    var_error <-
    NCP <-
    BH <-
    sample_size <-
    power <-
    method <-
    NULL

  associations <- data.table::data.table(associations)
  associations[, var_error := var_y - (beta^2 * 2 * maf * (1 - maf))]
  associations <- associations[var_error > 0]

  analytic_results_fdr_n <- data.table::data.table()

  if (verbose) {
    message(" FDR threshold used (q):     ", q)
    message(" Range of error variance:    ", prettyNum(min(associations$var_error)), " - ", prettyNum(max(associations$var_error)))
    message(" Number of tests considered: ", prettyNum(n_tests, big.mark = ","))
    message("")
  }

  if (sample_values & n_tests == nrow(associations)) {
    stop(" sample_values should only be used when n_tests is more than nrow(associations)")
  }

  if (sample_values & n_tests != nrow(associations)) {
    if (!data.table::between(sample_prop, 0, 1, incbounds = FALSE)) {
      stop("Sample proportion should be between 0 and 1")
    }

    if (verbose) {
      message(" Sampling for omitted variants (", prettyNum(sample_prop * 100), "% of given associations)\n")
    }

    omitted_h2 <- stats::rbinom(n_tests - nrow(associations), size = 1, prob = sample_prop)

    omitted_h2[omitted_h2 == 1] <- sample(
      (2 * associations$maf * (1 - associations$maf) * (associations$beta^2))/associations$var_error,
      size = sum(omitted_h2 == 1),
      replace = TRUE
    )
  } else if (!sample_values & n_tests != nrow(associations)) {
    if (verbose) {
      message(" All omitted associations assumed to be null\n")
    }
  }

  for (N in sample_sizes) {
    if (verbose) {
      message(" N = ", N)
    }

    associations[, NCP := (2 * N * maf * (1 - maf) * (beta^2))/var_error]

    if (sample_values) {
      ncp_vals <- c(associations$NCP, omitted_h2 * N)
    } else {
      ncp_vals <- associations$NCP
    }

    ## Find p-value threshold that results in a pFDR of q
    threshold_optim <- stats::uniroot(
      fdr_p_threshold,
      interval  = c(.Machine$double.eps, 1),
      ncp       = ncp_vals,
      M         = n_tests,
      q         = q,
      extendInt = "yes"
    )

    pval_threshold = threshold_optim$root

    if (verbose) {
      message("  p-value threshold for q: ", prettyNum(pval_threshold))
    }

    critical_val = stats::qchisq(1 - pval_threshold, 1)

    associations[, BH := 1 - stats::pchisq(critical_val, 1, NCP) ]

    associations[, sample_size := N]

    analytic_results_fdr_n <- rbind(analytic_results_fdr_n, associations)
  }

  analytic_results_fdr_n <- data.table::melt(
    analytic_results_fdr_n,
    measure.vars = "BH",
    variable.name = "method",
    value.name = "power"
  )

  if (is.null(by)) {
    by <- "sample_size"
  } else {
    by <- c(by, "sample_size")
  }

  analytic_results_fdr_n <- analytic_results_fdr_n[
    ,
    list(
      power = mean(power),
      method = unique(method)
    ),
    by = by
  ]

  analytic_results_fdr_n[, n_tests := n_tests]

  analytic_results_fdr_n
}


#' Analytic Power Calculation for a Fixed p-value Threshold
#'
#' This function calculates the long-run statistical power to detect eQTLs when
#' using a fixed p-value threshold or the Bonferroni method to determine
#' significance.
#'
#' @param associations A dataset of association data. It must contain the
#' columns described below (\code{pid}, \code{snp}, \code{beta}, \code{maf})
#' @param var_y The variance of the response
#' @param sample_sizes A set of sample sizes to test
#' @param n_tests The number of tests in the full dataset
#' @param pval_threshold The level the FWER is controlled at by the Bonferroni
#' adjustment.
#' @param by A set of variable names to aggregate power over
#' @param verbose If \code{TRUE}, messages will be output during the course of
#' the function.
#'
#' @return The estimated statistical power for each gene-SNP association when
#' the \eqn{p}-value threshold is set at \code{pval_threshold}.
#' @export
#'
#' @examples
#' data("associations")
#' bonf_results.dt <- analytic_calculation_fixed(associations.dt)
analytic_calculation_fixed <- function(
  associations,
  var_y = 1,
  sample_sizes = c(50, 100, 250, 500),
  n_tests = nrow(associations),
  pval_threshold = 0.05 / n_tests,
  verbose = TRUE,
  by = NULL
) {
  var_error <-
    maf <-
    NCP <-
    Bonferroni <-
    sample_size <-
    power <-
    method <-
    NULL

  associations <- data.table::data.table(associations)
  associations[, var_error := var_y - (beta^2 * 2 * maf * (1 - maf))]
  associations <- associations[var_error > 0]

  analytic_results_fix_n <- data.table::data.table()

  if (verbose) {
    message(" p-value threshold used:     ", pval_threshold)
    message(" Range of error variance:    ", prettyNum(min(associations$var_error)), " - ", prettyNum(max(associations$var_error)))
    message(" Number of tests considered: ", prettyNum(n_tests, big.mark = ","))
    message("")
  }

  for (N in sample_sizes) {
    if (verbose) {
      message(" N = ", N)
    }

    associations[, NCP := (2 * N * maf * (1 - maf) * (beta^2))/var_error]

    critical_val = stats::qchisq(1 - pval_threshold, 1)

    associations[, Bonferroni := 1 - stats::pchisq(critical_val, 1, NCP) ]

    associations[, sample_size := N]

    analytic_results_fix_n <- rbind(analytic_results_fix_n, associations)
  }

  analytic_results_fix_n <- data.table::melt(
    analytic_results_fix_n,
    measure.vars = "Bonferroni",
    variable.name = "method",
    value.name = "power"
  )

  if (is.null(by)) {
    by <- "sample_size"
  } else {
    by <- c(by, "sample_size")
  }

  analytic_results_fix_n <- analytic_results_fix_n[
    ,
    list(
      power = mean(power),
      method = unique(method)
    ),
    by = by
  ]

  analytic_results_fix_n[, n_tests := n_tests]

  analytic_results_fix_n
}
