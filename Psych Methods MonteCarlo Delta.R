##########################################################################################
# This R program runs a Monte Carlo study of meta-analysis methods.
#  It compares methods of estimating the lower bound of the credibility
#   and setting the confidence interval for that lower bound estimate.
#   Lower bound estimates by (a) Hedges' z, and Schmidt & Hunter with and without k correction.
#  The paramater set is... condition, delta, SD.delta, number of studies(k), 
#   sample sizes (average N).
#  We follow the simulation design used by 
#  Sanchez-Meca & Marin-Martinez (2008).  Confidence intervals for the overall effect size
#   in random-effects meta-analysis.  Psychological Methods, 13, 31-48.
##########################################################################################
#############################Clear Console##################################
rm(list=ls()) # cleans (empties)global environment
dev.off()     # cleans plots/graphics
cat("\014")   # cleans out the console
############################################################################
# If not installed, install package e1071
#install.packages('e1071')
# Run this line just once and then comment out.
# library(e1071)
# install.packages("beepr")

options(scipen=999) # this suppresses scientific notation (exponential notation)
install.packages("metafor")
install.packages("boot")
install.packages("sn")
install.packages("Hmisc")
library(metafor)              # need this for Morris estimates
library(MASS)                 # need this library for multivariate sampling
library(boot)                 # need this for bootstrapping
library(sn)
library(Hmisc)
###########################################################################
# set parameters for the simulation
###########################################################################
set.seed(1234)
#Shape of delta distribution, values at Normal, LeftSkew, and RightSkew
Condition = "LeftSkew"
# number of meta-analysis in the simulation run = counter for the main simulation loop for Monte Carlo
numbermetas <- 2000
# delta values at .5 and .8
delta <- .5                         # underlying mean of infinite-sample effect sizes
# tau-squared = .04, .08, .16, .32
# tau = .2, .283, .4, .566 (tau is SD.delta)
SD.delta <-.283                       # underlying SD of infinite-sample effect sizes
# sample sizes
avesamsize <- 30                    # pick one average sample size: 30, 50, 80, or 100
# number of independent samples in one meta-analysis (value of k)
# values are set to one of these:   10, 20, 40, 100 (skip 5)
numbersamples <- 10                 # set k the number of studies in 1 meta-analysis

###########GENERATING DISTIBUTIONS OF DELTA#################
if (Condition == 'Normal'){
  CR80L.pop <- delta-1.28*SD.delta                  # 80 percent credibility bounds
  CR80U.pop <- delta+1.28*SD.delta
}
# Find data for skewed distributions
#
if (Condition != 'Normal') {
  deltaskew = 10 #delta skew is set constant across all skewed distributions
  #xi is mean
  #omega is sd
  #alpha is skewness
  skew.parms <- c(1.6439, .691531, 1) # .5, .283
  ##Initial Skewed Distribution Generation
  distribution = rsn(n = 10000000, xi=delta,  omega=SD.delta, alpha=deltaskew)
  # selected for admissible bound of the correlation
  delta.distribution = sample(distribution, size = 10000000, replace = T)
  right.skew <- delta.distribution*skew.parms[1] - skew.parms[2]
  left.skew <- (right.skew*-1) + skew.parms[3]
  CR80L.pop <- delta-1.282*SD.delta  # population 80% credibility interval lower bound
  CR80U.pop <- delta+1.282*SD.delta  # upper bound
}
if (Condition == 'LeftSkew') {delta.distribution <- left.skew
true.lower.bound <- quantile(delta.distribution, .1)
delta = mean(delta.distribution)          # population correlation
SD.delta = sd(delta.distribution)          # population standard deviation
CR80L.pop <- quantile(delta.distribution,.1)  # population 80% credibility interval lower bound
CR80U.pop <- quantile(delta.distribution,.9)  # upper bound
}
if (Condition == 'RightSkew') {delta.distribution <- right.skew
true.lower.bound <- quantile(delta.distribution, .1)
delta = mean(delta.distribution)          # population correlation
SD.delta = sd(delta.distribution)          # population standard deviation
CR80L.pop <- quantile(delta.distribution,.1)  # population 80% credibility interval lower bound
CR80U.pop <- quantile(delta.distribution,.9)  # upper bound
}

ssizes <- 1:20
ssizes <- matrix(ssizes,5,5)
ssizes[1, ] <-c(12,16,18,20,84)    # average N = 30
ssizes[2, ] <-c(32,35,38,40,104)   # ave N = 50
ssizes[3, ] <-c(62,66,68,70,134)   # ave N = 80
ssizes[4, ] <-c(82,86,88,90,154)   # ave N = 100
# ssizes
#
if (avesamsize == 30) {ssizes1 <- ssizes[1, ]}
if (avesamsize == 50) {ssizes1 <- ssizes[2, ]}
if (avesamsize == 80) {ssizes1 <- ssizes[3, ]}
if (avesamsize == 100) {ssizes1 <- ssizes[4, ]}
#

#
Nreps <- numbersamples/5           # Nreps is the number of times each set of sample sizes (ssizes) is repeated
#

# collect the parameters
input.parms <- cbind(delta,SD.delta, avesamsize, numbersamples, CR80L.pop, CR80U.pop)
input.parms
##############################################################################
#  Find quantities for the Lawless Interval
#  Set p, the percentile of the underlying distribution for the
#  Credibility Interval.  If you want the conventional 80 percent
#  (2-tailed) interval, set p =.10 to cut off the bottom 10 percent.
p <- .10               # percentile of the normal for lower bound of CR
#  Set alphaLB, the error rate for the confidence interval for the LOWER BOUND.
#  This is the CI for the bottom estimate, not for the mean or the SD of effect sizes.
alphaLB <- .05
#
alpha1 <- alphaLB/2
alpha2 <- 1 - alpha1
alpha3 <- 1- alphaLB
#
true.lower.bound <- qnorm(p,mean=delta, sd=SD.delta)
#
N <- numbersamples                                # information needed for tolerance intervals
zp <- qnorm(1-p)                                  # for 80 CR z about 1.28
delta.NC <- sqrt(N)*zp                               # find noncentrality parameter
#
upper.bound <- (qt(alpha1,N-1,delta.NC))/sqrt(N)      # find t values for Hahn-Meeker tolerance intervals
middle.bound <- (qt(alpha3,N-1,delta.NC))/sqrt(N)
lower.bound <- (qt(alpha2,N-1,delta.NC))/sqrt(N)
############################################################################
# set placeholders for the simuation to hold the results
Hedz <- 1:numbermetas
Higt <- 1:numbermetas
HSCR <- 1:numbermetas
Hedzbot <- 1:numbermetas
Hedztop <- 1:numbermetas
Higtbot <- 1:numbermetas
Higttop <- 1:numbermetas
HSbot <- 1:numbermetas
HStop <- 1:numbermetas
HedM <- 1:numbermetas
SDMeans <- 1:numbermetas
Results <- matrix(-9,numbermetas,20) #placeholder for results - 
#  replications by the number of desired items 

############################################################################
# MAIN Loop - run as many times as you want replications
# set/reset placeholders for one meta
for (j in 1:numbermetas){
  out1 <- 1:numbersamples
  n1all<- 1:numbersamples
  n2all<- 1:numbersamples
  # Compute 1 meta-analysis Loop - compute study data and meta-analyze
  sam=0
  for (rep in 1:Nreps) {
    for (i in 1:5){  
      # set the local value of delta, sampled from the parameter set
      delta.i <-rnorm(1,delta,SD.delta)
      sam = sam+1
      n1 <-ssizes1[i]/2
      n2 <-n1
      sample1 <- rnorm(n1,delta.i,1) #experimental group data
      sample2 <- rnorm(n2,0,1) #control group data
      df1 <- n1-1
      df2 <- n2-1
      Ms1 = mean(sample1)
      Vs1 =var(sample1)
      Ms2 = mean(sample2)
      Vs2 = var(sample2)
      d = (Ms1-Ms2)/sqrt((Vs1*df1+Vs2*df2)/(df1+df2)) # compute local value of d
      #output d and sample size for each group
      out1[sam] <- d
      n1all[sam] <- n1
      n2all[sam] <- n2
    } # end meta-analysis loop (inner loop)
  } # end reps loop
  # we have sampled data for 'numbersamples' or k studies
  meta <- cbind(out1,n1all,n2all) # meta has d, n1 and n2 - what we need for a meta-analysis
  k <- length(out1)               # number of effect sizes (studies)
  bootdata <- data.frame(out1,n1all, n2all)
  bootdata
  ##############################################HEDGES###############################################
  # compute the Hedges meta-analysis in d (use g, unbiasd estimator of standarized mean difference)
  Vdi <- ((meta[,2] + meta[,3]) / (meta[,2] * meta[,3])) + (meta[,1]^2 / (2 * (meta[,2] + meta[,3])))
  gi <- (1 - (3 / (4 * (meta[,2]+meta[,3] - 3)-1))) * meta[,1]
  Vgi <- (1 - (3 / (4 * (meta[,2]+meta[,3] - 3)-1)))^2 * Vdi
  # run metafor for meta-analysis
  res1 <- rma(yi=gi, vi=Vgi,
              control=list(maxiter=1000, stepadj=.5), method="REML")               
  HedLB <- res1$b - 1.282*sqrt(res1$tau2)  # compute 80% lower bound
  # res1$bet is the mean
  # res1$tau2 is the random-effects variance component
  res1
  HedLB
  #
  #####################################################
  # Analytic function for Hedges
  #Analytic Approach For 95% Tolerance Interval Around 10 Percntile
  M <- res1$b
  sd <- sqrt(res1$tau2)
  
  Hedges.Percentile.Lower.CI <- M - lower.bound*sd
  Hedges.Percentile.Middle.CI <- M - middle.bound*sd
  Hedges.Percentile.Upper.CI <- M - upper.bound*sd
  Hedges.Percentile.CI.Width <- Hedges.Percentile.Upper.CI-Hedges.Percentile.Lower.CI
  ################################################Hunter & Schmidt####################################
  # compute the H&S meta-analysis in d (use g, unbiased estimator of standarized mean difference)
  #
  HSVdi <- ((meta[,2] + meta[,3]) / (meta[,2] * meta[,3])) + (meta[,1]^2 / (2 * (meta[,2] + meta[,3])))
  HSwi <- (meta[,2] + meta[,3])
  HSgi <- (1-(3/((4*(meta[,2]+meta[,3]-2))-1)))*meta[,1]
  HSVgi <- ((1-(3/((4*(meta[,2]+meta[,3]-2))-1)))^2)*HSVdi
  HSwigi <- HSwi*HSgi
  HSsumwigi <- sum(HSwigi)
  HSsumwi <- sum(HSwi)
  HSmeang <- HSsumwigi/HSsumwi
  HSdevg <- HSgi-HSmeang
  HSdevgsq <- HSdevg*HSdevg
  HSNdevgsq <- (meta[,2]+meta[,3])*HSdevgsq
  HSsumNdevgsq <- sum(HSNdevgsq)
  HSsumN <- sum(meta[,2]+meta[,3])
  HSVg <-HSsumNdevgsq/HSsumN
  HSVgk <- (k)/(k-1)*HSVg
  HSSEg <- (HSVg)^.5
  HSNVgi <- (meta[,2]+meta[,3])*HSVgi
  HSsumNVgi <- sum(HSNVgi)
  HSVser <-HSsumNVgi/HSsumN
  HStausq <- HSVg-HSVser
  HStausqk <- HSVgk-HSVser
  if(HStausq <= 0) {HStausq <- .0000001}
  if(HStausqk <= 0) {HStausqk <- .0000001}
  HStau <- HStausq^.5
  HStauk <- HStausqk^.5
  HSVgk <- HSVg/numbersamples #needed for variance of the mean calculation
  HSVgksqrt <- HSVgk^.5 
  HSl80 <- HSmeang - 1.282 * HStau  #HS lower bound
  HSl80k <- HSmeang - 1.282 * HStauk #HS-k lower bound
  #H&S Analytic Approach For 95% Tolerance Interval Around 10 Percentile
  HS.Percentile.Lower.CI <- HSmeang - lower.bound*HStau
  HS.Percentile.Upper.CI <- HSmeang - upper.bound*HStau
  HS.Percentile.CI.Width <- HS.Percentile.Upper.CI - HS.Percentile.Lower.CI
  #
  HSk.Percentile.Lower.CI <- HSmeang - lower.bound*HStauk
  HSk.Percentile.Upper.CI <- HSmeang - upper.bound*HStauk
  HSk.Percentile.CI.Width <- HSk.Percentile.Upper.CI - HSk.Percentile.Lower.CI
  ##########################################################
  # Bootstrap estimates for Hedges, HS, HS k-corrected method
  Boot.f <- function(d, i){ 
    d2 <- d[i,]
    boot.out1 <-d2$out1
    boot.n1 <- d2$n1all
    boot.n2 <- d2$n2all
    #Hedges Method
    b.Vdi <- ((boot.n1+boot.n2) / (boot.n1 * boot.n2)) + (boot.out1^2 / (2 * (boot.n1+boot.n2)))
    b.gi <- (1-(3/((4*(boot.n1+boot.n2-2))-1)))*boot.out1
    b.Vgi <- ((1-(3/((4*(boot.n1+boot.n2-2))-1)))^2)*b.Vdi
    res2 <- rma(yi=b.gi, vi=b.Vgi,
                control=list(maxiter=1000, stepadj=.5), method="REML")               
    HedLB2 <- res2$b - 1.282*sqrt(res2$tau2)
    #HS and HS k-corrected Method
    b.HSVdi <- ((boot.n1+boot.n2) / (boot.n1 * boot.n2)) + (boot.out1^2 / (2 * (boot.n1+boot.n2)))
    b.HSwi <- (boot.n1+boot.n2)
    b.HSgi <- (1-(3/((4*(boot.n1+boot.n2-2))-1)))*boot.out1
    b.HSVgi <- ((1-(3/((4*(boot.n1+boot.n2-2))-1)))^2)*b.HSVdi
    b.HSwigi <- b.HSwi*b.HSgi
    b.HSsumwigi <- sum(b.HSwigi)
    b.HSsumwi <- sum(b.HSwi)
    b.HSmeang <- b.HSsumwigi/b.HSsumwi
    b.HSdevg <- b.HSgi-b.HSmeang
    b.HSdevgsq <- b.HSdevg^2
    b.HSNdevgsq <- (boot.n1+boot.n2)*b.HSdevgsq
    b.HSsumNdevgsq <- sum(b.HSNdevgsq)
    b.HSsumN <- sum(boot.n1+boot.n2)
    b.HSVg <-b.HSsumNdevgsq/b.HSsumN
    b.HSVgk <- (k)/(k-1)*b.HSVg
    b.HSSEg <- (b.HSVg)^.5
    b.HSNVgi <- (boot.n1+boot.n2)*b.HSVgi
    b.HSsumNVgi <- sum(b.HSNVgi)
    b.HSVser <-b.HSsumNVgi/b.HSsumN
    b.HStausq <- b.HSVg-b.HSVser
    b.HStausqk <- b.HSVgk-b.HSVser
    if(b.HStausq <= 0) {b.HStausq <- .0000001}
    if(b.HStausqk <= 0) {b.HStausqk <- .0000001}
    b.HStau <- b.HStausq^.5
    b.HStauk <- b.HStausqk^.5
    b.HSl80 <- b.HSmeang - 1.282 * b.HStau  #HS lower bound
    b.HSl80k <- b.HSmeang - 1.282 * b.HStauk #HS-k lower bound
    boot.lowerbounds <- c(HedLB2,b.HSl80,b.HSl80k)
    return(boot.lowerbounds)
  } # end Boot.f function for bootstrapping
  # Run bootstrap
  Bootstrap <- boot(bootdata, Boot.f, R = 2000)
  #Hedges Bootstrapped Statistic
  Hed.bootstrap1 <- Bootstrap$t0[1]
  #HS Bootstrapped Statistic
  HS.bootstrap1 <- Bootstrap$t0[2]
  #HS k-corrected Bootstrapped Statistic
  HS.K.bootstrap1 <- Bootstrap$t0[3]
  #Hedges Bootstrap Confidence Interval
  Hed.bootstrap2 <- boot.ci(Bootstrap, type="bca",index=1)
  Hed.bootstrapLB <- Hed.bootstrap2$bca[4]
  Hed.bootstrapUB <- Hed.bootstrap2$bca[5]
  Hed.bootstrap.CI.Width <- Hed.bootstrapUB-Hed.bootstrapLB
  Hed.bootstrapLB
  Hed.bootstrapUB
  #HS Bootstrap Confidence Interval
  HS.bootstrap2 <- boot.ci(Bootstrap, type="bca",index=2)
  HS.bootstrapLB <- HS.bootstrap2$bca[4]
  HS.bootstrapUB <- HS.bootstrap2$bca[5]
  HS.bootstrap.CI.Width <- HS.bootstrapUB-HS.bootstrapLB
  HS.bootstrapLB
  HS.bootstrapUB
  #HS k-corrected Bootstrap Confidence Interval
  HS.K.bootstrap2 <- boot.ci(Bootstrap, type="bca",index=3)
  HS.K.bootstrapLB <- HS.K.bootstrap2$bca[4]
  HS.K.bootstrapUB <- HS.K.bootstrap2$bca[5]
  HS.K.bootstrap.CI.Width <- HS.K.bootstrapUB-HS.K.bootstrapLB
  HS.K.bootstrapLB
  HS.K.bootstrapUB
  
  ###Coverage Check###
  HS.Percentile.coverage <-0
  if(HS.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= HS.Percentile.Upper.CI) {HS.Percentile.coverage<-1}
  HS.bootstrap.coverage <-0
  if(HS.bootstrapLB<= true.lower.bound & true.lower.bound<=  HS.bootstrapUB) {HS.bootstrap.coverage<-1}
  HS.K.Percentile.coverage <-0
  if(HSk.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= HSk.Percentile.Upper.CI) {HS.K.Percentile.coverage<-1}
  HS.K.bootstrap.coverage <-0
  if(HS.K.bootstrapLB<= true.lower.bound & true.lower.bound<=  HS.K.bootstrapUB) {HS.K.bootstrap.coverage<-1}
  Hedges.Percentile.coverage <-0
  if(Hedges.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= Hedges.Percentile.Upper.CI) {Hedges.Percentile.coverage<-1}
  Hedges.bootstrap.coverage <-0
  if(Hed.bootstrapLB<= true.lower.bound & true.lower.bound<= Hed.bootstrapUB) {Hedges.bootstrap.coverage<-1}
  ##########################################
  # Output results
  ##########################################
  Results[j, ] <- cbind(delta, SD.delta, CR80L.pop, # 1, 2, 3
                        HSl80, HSl80k, HedLB[1], # 4, 5, 6
                        HS.Percentile.coverage, HS.bootstrap.coverage,         # 7, 8
                        HS.K.Percentile.coverage, HS.K.bootstrap.coverage,      # 9, 10
                        Hedges.Percentile.coverage, Hedges.bootstrap.coverage, # 11, 12
                        HS.Percentile.CI.Width, HSk.Percentile.CI.Width, Hedges.Percentile.CI.Width[1],  #13, 14, 15
                        HS.bootstrap.CI.Width, HS.K.bootstrap.CI.Width, Hed.bootstrap.CI.Width, #16, 17, 18
                        numbersamples,avesamsize)   #19, 20                                                            
} # END outer loop.
# Rename the columns of the Results

resultsnames = c("delta", "SD.delta", "true lower bound",                         #1, 2, 3
                 "HS Lower Bound", "HS.K Lower Bound",  "Hedges Lower Bound",   # 4 5, 6
                 "HS Percentile Coverage", "HS Bootstrap Coverage",             # 7, 8
                 "HS.K Percentile Coverage", "HS.K Bootstrap Coverage",         # 9, 10
                 "Hedges Percentile Coverage", "Hedges Bootstrap Coverage",     #11, 12
                 "HS Percentile CI Width", "HS k Percentile CI Width", "Hedges Percentile CI Width", #13, 14, 15
                 "HS bootstrap CI Width", "HS k bootstrap CI Width",  "Hedges bootstrap CI Width", #16, 17, 18
                 "Number of Samples", "Average Sample Size")        #19, 20                                                           
colnames(Results) = resultsnames

# Calculate the biases for the 3 estimators.
HS.Bias <- sum(Results[,4]-Results[,3])/numbermetas
HS.K.Bias <-sum(Results[,5]-Results[,3])/numbermetas
Hedges.Bias <-sum(Results[,6]-Results[,3])/numbermetas

# Calculate RMSE for each estimator of the lower bound
HS.RMSE <- sqrt(sum((Results[,4]-Results[,3])^2)/numbermetas)
HS.K.RMSE <- sqrt(sum((Results[,5]-Results[,3])^2)/numbermetas)
Hedges.RMSE <- sqrt(sum((Results[,6]-Results[,3])^2)/numbermetas)
#
#Calculate average Width of confidence interval estimate
HS.Percentile.width <- mean(Results[,13])
HS.Bootstrap.width <- mean(Results[,14])
HS.k.Percentile.width <- mean(Results[,15])
HS.k.Bootstrap.width <- mean(Results[,16])
Hedges.Percentile.width <- mean(Results[,17])
Hedges.Bootstrap.width <- mean(Results[,18])
#
#Calculate the coverage for the estimates
HS.Percentile.coverage <-mean(Results[,7])
HS.Bootstrap.coverage <-mean(Results[,8])
HS.K.Percentile.coverage <-mean(Results[,9])
HS.K.Bootstrap.coverage <-mean(Results[,10])
Hedges.Percentile.coverage <-mean(Results[,11])
Hedges.Bootstrap.coverage <-mean(Results[,12])

# Create database for output to collect across runs
Output <- data.frame(delta, SD.delta, avesamsize, numbersamples, numbermetas,
                     HS.Bias, HS.K.Bias, Hedges.Bias,
                     HS.RMSE, HS.K.RMSE, Hedges.RMSE,
                     HS.Percentile.coverage, HS.Bootstrap.coverage,
                     HS.K.Percentile.coverage, HS.K.Bootstrap.coverage,
                     Hedges.Percentile.coverage, Hedges.Bootstrap.coverage,
                     HS.Percentile.width,HS.Bootstrap.width,HS.k.Percentile.width,HS.k.Bootstrap.width,
                     Hedges.Percentile.width,Hedges.Bootstrap.width,Condition)

Output


# THIS WILL EXPORT RESULTS TO WORKING DIRECTORY 
# NEED TO CHANGE RESULTS TO CONDITION NUMBER; i.e condition saved as results1.csv 
write.table(Results, file = "results.csv", row.names=F, append=T, col.names=T, sep=",")

# THIS WILL EXPORT OUTPUT TO WORKING DIRECTORY 
write.table(Output, file = "output.csv", row.names=F, append=T, col.names=T, sep=",")
