#!/usr/bin/env Rscript


args=(commandArgs(TRUE))

data <- read.table(args[1]);

M <- as.table(rbind(c(data[1,]$V1, data[1,]$V2), c(data[1,]$V3,data[1,]$V4)))
dimnames(M) <- list(gender = c("Deam","noDeam"),
                       party = c("DeamE", "noDeamE"))
(Xsq <- chisq.test(M))  # Prints test summary

M <- as.table(rbind(c(data[2,]$V1, data[2,]$V2), c(data[2,]$V3,data[2,]$V4)))
dimnames(M) <- list(gender = c("Deam","noDeam"),
                       party = c("DeamE", "noDeamE"))
(Xsq <- chisq.test(M))  # Prints test summary





#plot(density(dataa$V1),xlab="xlab",ylab="ylab",main="main",col="red");
#legend("topright", title="sample", c("data"), lty=c(1), lwd=c(2.5),col=c("red")) 
