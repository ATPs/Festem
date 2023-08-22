#' Festem: Feature Selection and Differential Expression Gene Analysis via EM-test
#'
#' Directly selecting differentially expressed genes for single-cell clustering analyses via EM-test. "Festem" also provides EM-test to test whether the number of component in mixtures of normal, Poisson or negative binomial distributions are larger than one. See Chen, Wang, et al (2023) <https://doi.org/10.1101/2023.07.26.550670> and Wang, Chen and Xi (2023) <https://doi.org/10.48550/arXiv.2306.12671> for more details.

#' @name Festem
#' @docType package
#' @useDynLib Festem
## usethis namespace: start
#' @import RcppEigen
#' @importFrom nloptr nl.grad
#' @importFrom nloptr nloptr
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnbinom
#' @importFrom stats lm
#' @importFrom stats p.adjust
#' @importFrom stats pchisq
#' @importFrom stats quantile
#' @importFrom stats rpois
#' @importFrom stats var
## usethis namespace: end
NULL
