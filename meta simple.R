#Meta analysis code for RCTs only
setwd("~/Dropbox/WCMC/Phase 3/AOC/R")
library(esc)
library(meta)

effect <- read.table("Effectsize.csv",sep=",",header=TRUE)
#To analyze subset of data, chain above file to be called "main" and use line below:
#effect <- subset.data.frame(main,main$grp1sd!="NA")
#effect <- subset.data.frame(effect1,effect1$First.Author!="Moberg")

#Efficacy meta-analysis
metacont <- metacont(effect$grp1n,effect$grp1m,effect$grp1sd,effect$grp2n,
                     effect$grp2m,effect$grp2sd,data = effect, sm = "SMD",
                     method.smd="Hedges",studlab=effect$First.Author,hakn=TRUE,
                     method.tau = "SJ")

#Study completion rate effect size
metacont <- metabin(effect$grp1event,effect$grp1n,effect$grp2event,effect$grp2n,
                    data=effect,sm="RR",method="Inverse",studlab=effect$First.Author,
                    hakn=TRUE, method.tau = "SJ")
metacont
#metacont differs from metagen in the method for estimating tau and adjustment for random effects model
#Use random effects model

forest <- forest(metacont)
forest

ctrl.subgroup <- update.meta(metacont,byvar=effect$Overall.RoB,comb.random=TRUE, comb.fixed=FALSE)
ctrl.subgroup
forest.sub <- forest(ctrl.subgroup)
forest.sub

covar <- effect$Change.Cats
#change covar to whatever covar you want, and then rerun the code
metareg <- metareg(metacont,covar,intercept=FALSE)
metareg

bubble <- bubble(metareg,ylim=c(-1.25,0.75),xlim=c(-0.5,2.5),studlab=TRUE,
                 ylab="Treatment effect (log risk ratio)",
                 xlab="Difference in PSD features used in app vs. control")
bubble

#check for outliers
library(dmetar)
fo <- find.outliers(metacont)
forest(fo, col.predict = "blue")

inf.analysis <- InfluenceAnalysis(x = metacont, random = TRUE)
summary(inf.analysis)
plot(inf.analysis, "influence")
plot(inf.analysis, "baujat")
plot(inf.analysis, "es")
plot(inf.analysis, "i2")

#multiple metaregression - PSD FEATURES, EFFICACY
#Check for collinearity
library(PerformanceAnalytics)
library(metafor)
#check if the multiple covariates are correlated
#check which variables you want
effect[0,]
cor <- cor(effect[,16:17])
chart.Correlation(effect[,16:18])

model1 <- rma(yi = metacont$TE, 
              sei = metacont$seTE, 
              data = metacont, 
              method = "SJ", 
              mods = effect$Change.Cats,
              intercept=FALSE,
              test = "knha")
model1

model2 <- rma(yi = metacont$TE, 
              sei = metacont$seTE, 
              data = metacont, 
              method = "SJ", 
              mods = ~effect$Change.Fx + effect$Overall.RoB-1, 
              test="knha",
              intercept=FALSE)
model2

anova(model1, model2)

interaction.model <- rma(yi=metacont$TE,
                         sei=metacont$seTE, 
                         data=metacont, 
                         method = "SJ", 
                         mods = ~ effect$Overall.RoB*effect$Change.Cats-1, 
                         test="knha",
                         intercept=FALSE)
interaction.model

permutest(model2)

#funnel plot
funnel(metacont,xlab = "g",studlab = TRUE)

funnel(metacont, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("darkblue","blue","lightblue"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkblue","blue","lightblue"))

#need dmetar package
eggers.test(x = metacont)
metacont.trimfill<-trimfill(metacont)
metacont.trimfill
funnel(metacont.trimfill,xlab = "Hedges' g")

#Risk of bias plot
library(robvis)
data_rob2 <- read.table("Rob2.csv",sep=",",header=TRUE)
data_rob2
rob_summary <- rob_summary(data=data_rob2,tool="ROB2", overall=TRUE)
rob_summary
