#!/usr/bin/env Rscript
library(fitdistrplus)
library(MASS)


args=(commandArgs(TRUE))

data <- read.table(args[1]);

posdata<-data[data$V1>=1,]

df<-fitdistr(posdata, "lognormal")

print(df);
