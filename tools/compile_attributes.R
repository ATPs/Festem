#!/usr/bin/env Rscript

# Regenerate Rcpp exports
if (!require("Rcpp")) install.packages("Rcpp")
Rcpp::compileAttributes() 