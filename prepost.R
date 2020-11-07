setwd("~/Dropbox/WCMC/Phase 3/AOC/R")
library(esc)
library(metafor)
library(plotrix)

main <- read.table("prepost.csv",sep=",",header=TRUE)
#To analyze subset of data, change above file to be called "main" and use line below:
effect1 <- subset.data.frame(main,main$grp1sd!="NA")
#effect <- subset.data.frame(effect1,effect1$App!="MoodKit")
effect <- subset.data.frame(effect1,effect1$Outlier=="N")

#This code is for pre/post comparison
effect$First.Author..last.name.only. <- as.character(effect$First.Author..last.name.only.)
effect$First.Author..last.name.only. <- sort(effect$First.Author..last.name.only.)
effect$grp1m <- as.numeric(effect$grp1m)
effect$grp2m <- as.numeric(effect$grp2m)
effect$grp1sd <-as.numeric(effect$grp1sd)
effect$grp2sd <-as.numeric(effect$grp2sd)
effect$ri <- as.numeric(effect$ri)
meta <- escalc(measure="SMCC",m1i=effect$grp1m,m2i=effect$grp2m,sd1i=effect$grp1sd,
               sd2i=effect$grp2sd,ni=effect$grp1n,ri=effect$ri,slab=effect$First.Author..last.name.only.)
summary <- summary(meta) 
summary

metacont <- rma(yi,vi,data=meta,method="SJ",test="knha")
summary(metacont)

forest <- forest.rma(metacont)
forest

covar <- effect$total.fx
#change covar to whatever covar you want, and then rerun the code
metareg <- rma(yi,vi,mods=factor(covar),data=meta)
metareg

#check for outliers
regtest(metacont)
trimfill(metacont)

#Bubble plot - http://www.metafor-project.org/doku.php/plots:meta_analytic_scatterplot
preds <- predict(metareg,transf=exp)
size <- 1/sqrt(metareg$vi)
size <-size/max(size)
plot(NA, NA, xlim=c(0,5),ylim=c(0,4),xlab="Total PSD Categories",
     ylab="Standardized Mean Change",las=1,bty="l")

symbols(effect$total.cats,exp(meta$yi),circles=size,inches=FALSE,add=TRUE,bg="black")
