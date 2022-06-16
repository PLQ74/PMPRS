#' Title
#'
#' @param b_exp a vector of SNP effects on the exposure
#' @param maf_x minor allele frequency
#'
#' @return population mean polygenic risk scores
#' @export
#'
#' @examples
pmprs_x <- function(b_exp, maf_x){
  sum(b_exp * ((2 * maf_x * (1 - maf_x)) + 2 * maf_x^2))
}
#' Title
#'
#' @param b_exp a vector of SNP effects on the exposure
#' @param maf_x minor allele frequency
#' @param se_exp a vector of stander error of b_exp
#'
#' @return the variance of pmprs_x
#' @export
#'
#' @examples
var_x <- function(b_exp, maf_x, se_exp){
  sum(((2 * maf_x * (1 - maf_x)) + 2 * maf_x^2)^2) * sum(se_exp^2)
}
#' Title
#'
#' @param b_out a vector of SNP effects on the outcome
#' @param maf_y minor allele frequency
#'
#' @return population mean polygenic risk scores
#' @export
#'
#' @examples
pmprs_y <- function(b_out, maf_y){
  sum(b_out * ((2 * maf_y * (1 - maf_y)) + 2 * maf_y^2))
}
#' Title
#'
#' @param b_out a vector of SNP effects on the outcome
#' @param maf_y minor allele frequency
#' @param se_out a vector of stander error of b_out
#'
#' @return the variance of pmprs_y
#' @export
#'
#' @examples
var_y <- function(b_out, maf_y, se_out){
  sum(((2 * maf_y * (1 - maf_y)) + 2 * maf_y^2)^2) * sum(se_out^2)
}
#' Title
#'
#' @param b_exp a vector of SNP effects on exposure
#' @param b_out a vector of SNP effects on the outcome
#' @param se_exp a vector of stander error of b_exp
#' @param se_out a vector of stander error of b_out
#' @param maf_x minor allele frequency
#' @param maf_y minor allele frequency
#'
#' @return risk score ratio
#' @export
#'
#' @examples
RSR <- function(b_exp, b_out, se_exp, se_out, maf_x, maf_y){
  sum((b_out * ((2 * maf_y * (1 - maf_y)) + 2 * maf_y^2))) /
    sum((b_exp * ((2 * maf_x * (1 - maf_x)) + 2 * maf_x^2)))
}
#' Title
#'
#' @param b_exp
#' @param b_out
#' @param se_exp
#' @param se_out
#' @param maf_x
#' @param maf_y
#'
#' @return risk score ratio and its 95% confidence interval
#' @export
#'
#' @examples
RSR.CI <- function(b_exp, b_out, se_exp, se_out, maf_x, maf_y){
  library(Rmisc)
  RSR = ((b_out * ((2 * maf_y * (1 - maf_y)) + 2 * maf_y^2))) /
    ((b_exp * ((2 * maf_x * (1 - maf_x)) + 2 * maf_x^2)))
  CI(RSR, ci = 0.95)
}
