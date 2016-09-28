#!/usr/bin/env Rscript

#usage RSCRIPT [data log] [pdf out] <title> <expected contamination> 

args=(commandArgs(TRUE))

data <- read.table(args[1]);

pdffile <- args[2];
title <- "Posterior probability for contamination\nusing deamination patterns";
if(length(args)>2){ #title
title <- args[3];    
}


pdf(pdffile);
plot(data$V2,data$V3,xlab="Contamination",ylab="Posterior probability (log scale)",main=title,pch="*",col="darkred",lty=1);

if(length(args)>3){ #title
abline(v=args[4],lty=2);
}


                                        #legend("topright", title="sample", c("data"), lty=c(1), lwd=c(2.5),col=c("red")) 
dev.off();
