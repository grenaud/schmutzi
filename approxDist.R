#!/usr/bin/env Rscript-3.2.0
library(fitdistrplus)
library(MASS)


args=(commandArgs(TRUE))

data <- read.table(args[1]);


df<-fitdistr(data$V1, "lognormal")

print(df);
