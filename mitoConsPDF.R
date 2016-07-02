#!/usr/bin/env Rscript-3.2.0


args=(commandArgs(TRUE))

data <- read.table(args[1],header=TRUE,stringsAsFactors=FALSE,colClasses=c("character","character","character",rep("numeric",8)));
#data <- read.table(args[1],header=TRUE,stringsAsFactors=FALSE);
data$pos<-as.numeric(data$pos);

data<- data[complete.cases(data),]


pdf(paste(args[1],".qual.pdf",sep=""));
plot(data$pos,data$qual,xlab="Position on the MT (bp)",ylab="Posterior probability (PHRED scale)",main="Posterior probability with respect to the position on the MT",col="red",pch=4,lty=1,cex=0.05);
dev.off();


pdf(paste(args[1],".cov.pdf",sep=""));
plot(data$pos,data$cov,xlab="Position on the MT (bp)",ylab="Coverage",main="Coverage with respect to the position on the MT",col="red",pch=4,lty=1,cex=0.05);
dev.off();


pdf(paste(args[1],".div.pdf",sep=""));

plot( density(data$pos[data$refBase!=data$base]),xlab="Position on the MT (bp)",ylab="Presence of mismatch",main="Divergence with respect to the position on the MT",col="red",pch=4,lty=1);
#plot(data$pos,data$refBase!=data$base,xlab="Position on the MT (bp)",ylab="Presence of mismatch",main="Divergence with respect to the position on the MT",col="red",pch=4,lty=1);
dev.off();

#legend("topright", title="sample", c("data"), lty=c(1), lwd=c(2.5),col=c("red")) 
