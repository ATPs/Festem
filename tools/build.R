#!/usr/bin/env Rscript

# Install required packages
if (!require("devtools")) install.packages("devtools")

# Build and check package
devtools::document()
devtools::build() 