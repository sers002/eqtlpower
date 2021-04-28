#' Calculate statistical power to detect eQTLs based on adjusted p-values
#'
#' This function takes a set of coefficient results from
#' \code{\link{generate_coefficients}} and a adjusted p-value matrix generated
#' based on these values and calculates how many eQTLs were detected. Note this
#' function is not intended for typical usage as it is used within the
#'  \code{\link{analyse_coef_results}} function.
#'
#' @param coef_results A set of results for a given sample size from the
#' \code{generate_coefficients} function
#' @param adjusted_matrix A matrix of adjusted p-values
#' @param method_name A label for the method used to create
#' \code{adjusted_matrix}
#' @param by A set of variable names to aggregate power over
#' @param threshold A threshold for values in \code{adjusted_matrix} to be
#' considered significant
#' @param interval Confidence interval level
#'
#' @return The mean statistical power to detect eQTLs over the realisations for each group according to the \code{by} variables.
#' @export
sa_power_calc <- function(
  coef_results,
  adjusted_matrix,
  method_name,
  threshold = 0.05,
  interval  = 0.95,
  by        = NULL
) {
  value <-
    method <-
    NULL

  if (!is.null(by)) {
    coef_results <- coef_results[, by, with = FALSE]
  }
  ## Combine effect size/MAF bin and QTL metadata with discovery (adjusted p < 0.05),
  ## then calculate the percent detected for each group within the bins for a
  ## given realisation
  power_dt <- data.table(
    coef_results,
    adjusted_matrix < threshold
  )[,
    lapply(.SD, function(x) sum(x)/.N),
    by = by,
    .SDcols = paste0("V", 1:ncol(adjusted_matrix))
  ]

  alpha_lwr <- (1 - interval)/2
  alpha_upr <- 1 - (1 - interval)/2

  ## Calculate the average power over all realisations. Calculate the Monte Carlo
  ## interval based on the quantiles of power
  power_dt <- melt(
    power_dt,
    measure.vars = paste0("V", 1:ncol(adjusted_matrix))
  )[
    ,
    list(
      power = mean(value),
      power_lwr = stats::quantile(value, probs = alpha_lwr),
      power_upr = stats::quantile(value, probs = alpha_upr)
    ),
    by = by
  ]

  power_dt[, method := method_name]

  data.table::setkeyv(power_dt, cols = by)

  power_dt
}


#' Calculate statistical power to detect eGenes based on adjusted p-values
#'
#' This function takes a set of coefficient results from
#' \code{\link{generate_coefficients}} and a adjusted p-value matrix generated
#' based on these values and calculates how many eGenes were detected. Note this
#' function is not intended for typical usage as it is used within the
#'  \code{\link{analyse_coef_results}} function.
#'
#' @param coef_results A set of results for a given sample size from the
#' \code{generate_coefficients} function
#' @param adjusted_matrix A matrix of adjusted p-values
#' @param method_name A label for the method used to create
#' \code{adjusted_matrix}
#' @param by A set of variable names to aggregate power over
#' @param threshold A threshold for values in \code{adjusted_matrix} to be
#' considered significant
#' @param interval Confidence interval level
#'
#' @return The mean statistical power to detect eGenes over the realisations for each group according to the \code{by} variables.
#' @export
sa_power_calc_genes <- function(
  coef_results,
  adjusted_matrix,
  method_name,
  threshold = 0.05,
  interval  = 0.95,
  by        = NULL
) {
  value <-
    method <-
    NULL

  by <- NULL
  coef_results <- coef_results[, c("pid", by), with = FALSE]

  ## Combine effect size/MAF bin and QTL metadata with discovery (adjusted p < 0.05),
  ## then calculate the percent detected for each group within the bins for a
  ## given realisation
  power_dt <- data.table(
    coef_results,
    adjusted_matrix < threshold
  )[,
    lapply(.SD, any),
    .SDcols = paste0("V", 1:ncol(adjusted_matrix)),
    by = c("pid", by)
  ]

  alpha_lwr <- (1 - interval)/2
  alpha_upr <- 1 - (1 - interval)/2

  ## Calculate the average power over all realisations. Calculate the Monte Carlo
  ## interval based on the quantiles of power
  power_dt <- melt(
    power_dt,
    measure.vars = paste0("V", 1:ncol(adjusted_matrix))
  )[
    ,
    list(
      power = mean(value),
      power_lwr = stats::quantile(value, probs = alpha_lwr),
      power_upr = stats::quantile(value, probs = alpha_upr)
    ),
    by = by
  ]

  power_dt[, method := method_name]

  # data.table::setkeyv(power_dt, cols = by)

  power_dt
}

#' Analyse semi-analytic results
#'
#' This function takes the output of \code{\link{generate_coefficients}}
#' and calculates the average statistical power of either eQTLs
#' (\code{type = "eqtl"}) or eGenes (\code{type = "gene"}). The statistical
#' power is averaged within the groups defined by the \code{by} variable.
#'
#' @param associations A data frame containing the real associations. This
#' dataset must contain the columns described below (\code{pid}, \code{snp},
#'  \code{maf}, \code{beta}, and any columns defined in the \code{by} argument)
#' @param coef_results A set of results for a given sample size from the
#' \code{generate_coefficients} function
#' @param n_tests The number of tests in the full dataset
#' @param adjust_methods The methods to adjust \eqn{p}-values by to account for multiple comparisons. Options include: \code{BH} (Benjamini-Hochberg), \code{BY} (Benjamini–Yekutieli), \code{bonferroni} (Bonferroni).
#' @param by A set of variable names to aggregate power over
#' @param var_y The variance of the response
#' @param threshold Either a single value (if only one adjustment method is
#' used) or a named vector giving the thresholds used for significance.
#' @param interval Confidence interval level
#' @param type Type of power analysis, either \code{"eqtl"} or \code{"gene"}.
#' @param verbose If \code{TRUE}, messages will be output during the course of
#' the function.
#'
#' @return The average statistical power to detect either eQTLs or eGenes using the chosen adjustment methods for each sample size contained in \code{coef_results}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("associations")
#' data("ld")
#'
#' sa_coefs <- generate_coefficients(associations.dt, ld.dt)
#' sa_results <- analyse_coef_results(associations.dt, sa_coefs)
#' sa_results_genes <- analyse_coef_results(
#'   associations.dt,
#'   sa_coefs,
#'   type = "gene"
#' )
#' }
analyse_coef_results <- function(
  associations,
  coef_results,
  n_tests        = nrow(associations),
  adjust_methods = c("BH"),
  var_y          = 1,
  by             = NULL,
  threshold      = 0.05,
  interval       = 0.95,
  type           = "eqtl",
  verbose        = TRUE
) {
  ## This is required to prevent NOTEs when building package, as recommended by
  ## the data.table developers ("Importing data.table" vignette).
  var_error <-
    maf <-
    beta_hat_var <-
    sample_size <-
    NULL

  coef_results <- data.table::merge.data.table(
    associations,
    coef_results,
    by = c("pid", "snp")
  )

  ## Calculate the error variance
  coef_results[, var_error := var_y - (beta^2 * 2 * maf * (1 - maf))]
  coef_results[, beta_hat_var := var_error/(2 * sample_size * maf * (1 - maf))]

  # coef_results[var_error >]
  ## Datasets to store results of loop in
  power_results <- data.table::data.table()

  sample_sizes <- unique(coef_results$sample_size)

  calc_fun <- if (type == "eqtl") {
    if (verbose) {
      message("Estimating power of detection for eQTLs")
    }
    sa_power_calc
  } else if (type == "gene") {
    if (verbose) {
      message("Estimating power of detection for eGenes")
    }
    sa_power_calc_genes
  } else {
    stop("The type argument must be one of: 'eqtl', 'gene'")
  }

  if (verbose) {
    if (!is.null(by)) {
      message(" Aggregating by:     ", paste(by, collapse = ", "))
    }
    message(" Adjustment methods: ", paste(adjust_methods, collapse = ", "), "\n")

    message(" Range of error variance:    ", prettyNum(min(coef_results$var_error)), " - ", prettyNum(max(coef_results$var_error)))
    message(" Number of tests considered: ", prettyNum(n_tests, big.mark = ","))
    message(
      " Sample sizes:               ",
      paste(prettyNum(sample_sizes, big.mark = ","), collapse = "; ")
    )

    message("")
  }

  for (curr_n in sample_sizes) {
    if (verbose) {
      message(" N = ", curr_n)
    }

    coef_results_n <- coef_results[sample_size == curr_n]

    ## For each SNP, divide by its SE - realisations are the V1, V2, ... V100 columns
    z_mat <- as.matrix(coef_results_n[, grep("V[0-9]+", names(coef_results_n)), with = FALSE])
    z_mat <- sweep(z_mat, MARGIN = 1, sqrt(coef_results_n$beta_hat_var), FUN = "/")

    p_mat <- matrix(
      stats::pnorm(-abs(z_mat), 0, 1, lower = TRUE) + stats::pnorm(abs(z_mat), 0, 1, lower = FALSE),
      ncol = ncol(z_mat)
    )

    rm(z_mat)
    gc()

    coef_results_n <- coef_results_n[
      ,
      c(
        "pid",
        "snp",
        "beta",
        "maf",
        "beta_hat_var",
        "sample_size",
        by
      ),
      with = FALSE
    ]

    power_results_n <- data.table::data.table()

    if ("BH" %in% adjust_methods) {
      if (verbose) {
        message("  BH")
      }

      if ("BH" %in% names(threshold)) {
        bh_threshold <- threshold["BH"]
      } else if (length(threshold) == 1) {
        bh_threshold <- threshold
      } else {
        stop("If named vector is provided, it must have a threshold for each method")
      }

      f_mat_filt  <- apply(p_mat, 2, stats::p.adjust, method = "BH", n = n_tests)
      power_dt_filt_bh <- calc_fun(
        coef_results_n,
        f_mat_filt,
        method_name = "BH",
        by = by,
        interval = interval,
        threshold = bh_threshold
      )
      power_dt_filt_bh[, n_tests := n_tests]
      rm(f_mat_filt)

      power_results_n <- rbind(power_results_n, power_dt_filt_bh)
    }

    if ("BY" %in% adjust_methods) {
      if (verbose) {
        message("  BY")
      }

      if ("BY" %in% names(threshold)) {
        by_threshold <- threshold["BY"]
      } else if (length(threshold) == 1) {
        by_threshold <- threshold
      } else {
        stop("If named vector is provided, it must have a threshold for each method")
      }

      f_mat_filt  <- apply(p_mat, 2, stats::p.adjust, method = "BY", n = n_tests)
      power_dt_filt_bh <- calc_fun(coef_results_n, f_mat_filt, "BY", by = by, interval = interval, threshold = by_threshold)
      power_dt_filt_bh[, n_tests := n_tests]
      rm(f_mat_filt)

      power_results_n <- rbind(power_results_n, power_dt_filt_bh)
    }

    if ("bonferroni" %in% adjust_methods) {
      if (verbose) {
        message("  Bonferroni")
      }

      if ("bonferroni" %in% names(threshold)) {
        bon_threshold <- threshold["bonferroni"]
      } else if (length(threshold) == 1) {
        bon_threshold <- threshold
      } else {
        stop("If named vector is provided, it must have a threshold for each method")
      }

      f_mat_filt  <- apply(p_mat, 2, stats::p.adjust, method = "bonferroni", n = n_tests)
      power_dt_filt_bh <- calc_fun(coef_results_n, f_mat_filt, "bonferroni", by = by, interval = interval, threshold = bon_threshold)
      power_dt_filt_bh[, n_tests := n_tests]
      rm(f_mat_filt)

      power_results_n <- rbind(power_results_n, power_dt_filt_bh)
    }

    power_results_n[, sample_size := curr_n]

    power_results <- rbind(power_results, power_results_n)

    rm(p_mat)
    rm(power_results_n, coef_results_n)

    gc()
  }

  power_results
}
