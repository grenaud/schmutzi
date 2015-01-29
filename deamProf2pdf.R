#!/usr/bin/env Rscript


#argument : [damage-patterns substitution as input] [output pdf] [title]
library(gplots)
library(ggplot2)
library(reshape)


args=(commandArgs(TRUE))


data<-read.table(args[1],
                 skip=1,
                 col.names=c("A->C",
                   "A->G",
                   "A->T",
                   "C->A",
                   "C->G",
                   "C->T",
                   "G->A",
                   "G->C",
                   "G->T",
                   "T->A",
                   "T->C",
                   "T->G"
                   ));
data<-cbind(seq(1,length(data[,1])),data)
colnames(data)[1]<-"pos"

meltdata<-melt.data.frame(data,id.vars="pos",measure.vars=seq(2,13,by=1));
meltdata$variable<-sub(".","",sub(".","->",meltdata$variable,fixed=TRUE),fixed=TRUE)              


pdf(args[2]);
ggplot(meltdata)  + geom_line(aes(x=pos, y=value,colour=variable)) +xlab("Position on fragment (bp)") + ggtitle(args[3])+ ylab("Probability of substitution") + scale_colour_discrete(name = "Type of error")
dev.off();

