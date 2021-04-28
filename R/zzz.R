#' @import data.table
#' @importFrom data.table ":="
NULL

#' eqtlpower: A package for estimating eQTL Study Power
#'
#' The eqtlpower package provides two methods to estimate statistical power
#' for eQTL studies: semi-analytic and analytic.
#'
#' @section Analytic Power Calculation:
#' The eqtlpower package provides two methods for calculating statistical power analytically depending on the multiple comparisons adjustment required. The \code{\link{analytic_calculation_fdr}} function provides analytic calculation for the Benjamini-Hochberg procedure and the \code{\link{analytic_calculation_fixed}} function can be used for either fixed p-value thresholds (e.g. p < 0.05) or the Bonferroni adjustment.
#'
#' @section Semi-Analytic Power Calculation:
#' The semi-analytic power calculation method generates realisations of SNP-gene marginal coefficients based on the linkage disequilibrium that exists between variants. This method requires the use of two functions sequentially. Firstly, the user can generate a number of realisations of the coefficients using the \code{\link{generate_coefficients}} function; this results in a dataset containing a column for each realisation of the coefficients. This dataset can either be analysed by the user with custom code for a novel methodology, or the \code{\link{analyse_coef_results}} function can be used to calculate the estimated statistical power based on the multiple testing methods contained in \code{\link[stats]{p.adjust}}. Note there is an option to choose whether to calculate the power for either variants or genes in the analysis function.
#'
#' @section Visualising the Estimated Power:
#' The estimated statistical power for a number of sample sizes can be visualised using the \code{\link{plot_power_results}} function. This provides a basic visualisation of the results as a ggplot2 object which can be used with any of the usual ggplot2 layers. Note that this is just a basic method of visualising results; advanced users can easily produce their own custom plots using the power results datasets directly.
#'
#' @docType package
#' @name eqtlpower
NULL



#' Example Gene-SNP Associations
#'
#' This dataset is an example of required data for the \code{associations}
#' argument in the \code{\link{generate_coefficients}} and analytic calculation
#' functions.
#'
#' @format A data.table with 1000 rows and 6 variables:
#' \describe{
#'   \item{pid}{gene ID}
#'   \item{snp}{SNP ID}
#'   \item{qtl_tag}{whether the SNP is considered a true eQTL}
#'   \item{beta}{true marginal effect between the SNP and the gene}
#'   \item{maf}{minor allele frequency of SNP}
#'   \item{cis_flag}{whether the SNP is within 1 megabases of the gene}
#' }
#'
#' @name associations.dt
#' @docType data
#' @usage data("associations")
NULL

#' Example Estimated Linkage Disequilibrium
#'
#' This dataset is an example of required data for the \code{ld}
#' argument in the \code{\link{generate_coefficients}}. Each row of this dataset
#' is a SNP pair and the estimated linkage disequilibrium between them. This
#' dataset was produced by the PLINK program.
#'
#' @format A data.table with 1064 rows and 9 variables:
#' \describe{
#'   \item{CHR_A}{chromosome on which SNP A is located}
#'   \item{BP_A}{location on the chromosome in base pairs where SNP A is located}
#'   \item{SNP_A}{ID of SNP A}
#'   \item{MAF_A}{minor allele frequency of SNP A}
#'   \item{CHR_B}{chromosome on which SNP B is located}
#'   \item{BP_B}{location on the chromosome in base pairs where SNP B is located}
#'   \item{SNP_B}{ID of SNP B}
#'   \item{MAF_B}{minor allele frequency of SNP B}
#'   \item{R}{estimated linkage disequilibrium between SNP A and SNP B}
#' }
#'
#' @name ld.dt
#' @docType data
#' @usage data("ld")
NULL
