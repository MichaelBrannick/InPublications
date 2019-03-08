## Monte Carlo Meta-analysis - CI for the lower bound with indirect range restriction ##
#Installing, Loading Necessary Packages
# install.packages("MASS", lib="/data/Rpackages/")
# install.packages("boot", lib="/data/Rpackages/")
# install.packages("beepr")
#install.packages("truncnorm")
##################################################
# Data generation for meta-analysis
#  with range restriction and unreliability
#  Based on on notation by Le & Schmidt 2006
# S = suitability - unobserved variable responsible for indirect range restriction
# T = True score for X (Test or Trait)
# P = True score for Y (Performance)
# X = predictor variable (observed X, observed test score)
# Y = criterion variable (observed Y, observed criterion score)
#  Path diagram for the unobserved variables S -> T -> P (Le & Schmidt, p 418)
#  Path diagram for the unobserved to observed variables T -> X and P -> Y; others are zero
#  We will assign the relations between S and T, T and P at the structural level
#  The correlation between T and P can vary (random effects)
#  The correlation between S and P is the product of the correlations ST, TP
#  We will assign values of the reliability of X and of Y (can vary)
#  The correlation between the latent variables and X and Y are found by tracing rule
#   (e.g., rSX = rSTrTX; rSY = rSTrTPrPY)
#  Range restriction is found by sampling, sorting on S, selecting and computing uT or uX.
#  Local reliability is found by computing rTX-squared and rPY-squared for local studies.
#  We need a correlation matrix 5x5
#  Variable 1 = S
#  Variable 2 = T
#  Variable 3 = P
#  Variable 4 = X
#  Variable 5 = Y
#  The objective of the meta-analysis is to find the average correlation
#  between T and P and the SD of the correlation between T and P
#  given rxy, rxxi, ryyi
#  Indirect range restriction happens with S; direct RR with X
#############################Clear Console##################################
rm(list=ls()) # cleans (empties)global environment
dev.off()     # cleans plots/graphics
cat("\014")   # cleans out the console
############################################################################
#              Open libraries
install.packages("metafor")
install.packages("boot")
install.packages("sn")
install.packages("Hmisc")
library(metafor)              # need this for Morris estimates
library(MASS)                 # need this library for multivariate sampling
library(boot)                 # need this for bootstrapping
library(sn)
library(Hmisc)
#library(truncnorm)           # need to find population value of range restriction, u
#library(psych)               # for describing distributions
##############################################################################
################## PART 1: Set the Parameters#################################
#     Fixed factors:
# Condition: Normal, LeftSkew, RightSkew
# means rho:  25, .50 (Paterson approx 15=.15, 50=.26, 85th=.42 percentiles)
# SD(rho):    .05, .10, .20 (.08, .13, .20 Paterson approx 15, 50, 85th percentiles)
# Rho ST (correlation between suitability and test): .5, .8 (Le & Schmidt)
# Degree of indirect range restriction (IRR) - selection ratio on S: .1, .5, .9 (1 in 10, 5 in 10, 9 in 10)
# Number of samples per meta: 8,25,50,100 (Brannick)
#     Random factors:
# Sample size per study: (SmallN skewed distribution)
# Reliability of X and Y (Le & Schmidt 2006)
# Replications is 2000, so numbermetas = 2000
############################################################################
# Fixed factors:
numbermetas = 2000     # number of metas, each based on numbersamples
numbersamples = 8      # k, the number of studies per meta-analysis
rho.mean = .25          # population correlation
sigmarho = .1          # population standard deviation
RhoST = 1              # population correlation between suitability true score and test true score, value of 1 results in no indirect range restriction
IRR = 1                # selection ratio for indirect range restriction, value of 1 results in no range restriction
Condition = "LeftSkew" #Shape of the underlying effect size distribution
# Random factors:
smallN2 = rgamma(1000000, shape = .64, scale = 165)
smallN2 = ceiling(smallN2)
smallN = sample(smallN2[smallN2 > 29], size = 1000000, replace = T)
# Simulate reliability in predictor and criterion, value of 1 results in perfect reliability
rxxp <- c(1)
#
ryyp <- c(1)
##################################################################################
# housekeeping required to run the simulation

rhobox <- matrix(-9,5,5)+10*diag(5) # placeholder for correlation matrix for sampling
mu <- rep(0,5)                  # multivariate means of zero
rhonames <- cbind("S", "T", "P", "X", "Y")
colnames(rhobox) <- rhonames
rownames(rhobox) <- rhonames

if (Condition == 'Normal'){
  distribution <- rnorm(10000000,rho.mean,sigmarho)  #10 million random numbers from normal
  rho.distribution <-                           #select within admissible bounds of correlation
    sample(distribution[distribution <= .99 & distribution >= -.99], size = 10000000, replace = T)
  popCR.L <- rho-1.28*sigmarho                  # 80 percent credibility bounds
  popCR.U <- rho+1.28*sigmarho
  pop.width.CR <- popCR.U-popCR.L               # width of credibility interval
  pop.CR <- cbind(popCR.L, popCR.U, pop.width.CR)
}
# Find data for skewed distributions
#
if (Condition != 'Normal') {
  rhoskew = 10 #rho skew is set constant across all skewed distributions
  #xi is mean
  #omega is sd
  #alpha is skewness
  skew.parms <- c(.99, .608, .2919, .5) # .25, .1
  ##Initial Skewed Distribution Generation
  distribution = rsn(n = 10000000, xi=rho.mean,  omega=sigmarho, alpha=rhoskew)
  summary(distribution)
  # selected for admissible bound of the correlation
  rho.distribution = sample(distribution[distribution < skew.parms[1]], size = 10000000, replace = T)
  right.skew <- rho.distribution/skew.parms[2] - skew.parms[3]
  left.skew <- (right.skew*-1) + skew.parms[4]
  p <- .10               # percentile of the normal for lower bound of CR
  true.lower.bound <- qnorm(p,mean=rho.mean, sd=sigmarho)
}
if (Condition == 'LeftSkew') {rho.distribution <- left.skew
true.lower.bound <- quantile(rho.distribution, .1)
rho.mean = mean(rho.distribution)          # population correlation
sigmarho = sd(rho.distribution)          # population standard deviation
}
if (Condition == 'RightSkew') {rho.distribution <- right.skew
true.lower.bound <- quantile(rho.distribution, .1)
rho.mean = mean(rho.distribution)          # population correlation
sigmarho = sd(rho.distribution)          # population standard deviation
}
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
#
N <- numbersamples                                # information needed for tolerance intervals
zp <- qnorm(1-p)                                  # for 80 CR z about 1.28
delta <- sqrt(N)*zp                               # find noncentrality parameter
#
upper.bound <- (qt(alpha1,N-1,delta))/sqrt(N)      # find t values for Hahn-Meeker tolerance intervals
middle.bound <- (qt(alpha3,N-1,delta))/sqrt(N)
lower.bound <- (qt(alpha2,N-1,delta))/sqrt(N)
# find underlying correlations between T & P (true score on tests and true score on performance)
pop.r <- rnorm(1000000,rho.mean,sigmarho)              # sample 1 million correlations from population
pop.rs <- sample(pop.r[pop.r > -.99 & pop.r <.99], size = 1000000, replace = T) # resample to be in bounds
#
out1 <- 1:numbersamples                            # placeholder for data for one meta (inner loop)
Results <- matrix(-9,numbermetas,38)                # placeholder for simulation results (outer loop)
##########################   START the SIMULATION    ###################################
# outer loop - do this to replicate meta-analyses
for (j in 1:numbermetas){             # outer loop - do this to replicate metas
  ri.unrestrict <- 1:numbersamples     # Placeholders for meta-analysis sampled ri
  rxxa.pop <- ri.unrestrict            # sampled from rxxp - population value
  ryya.pop <- ri.unrestrict            # sampled from ryyp - population value
  rTP.pop <- ri.unrestrict             # sampled from underlying rho, SDrho
  #
  ri.un <- ri.unrestrict               # sample values with sampling error, no RR
  rxx.un <- ri.unrestrict
  ryy.un <- ri.unrestrict
  #
  ri.ind <- ri.unrestrict              # sample values with indirect RR
  rxx.ind <- ri.unrestrict
  ryy.ind <- ri.unrestrict
  uX.ind <- ri.unrestrict
  num.meta <- ri.unrestrict
  
  # Set the population values for collection of studies in a meta
  rxxa <-  sample(rxxp, size=numbersamples, replace = T)     # sample from reliability of x
  ryya <-  sample(ryyp, size=numbersamples, replace = T)
  rho <- sample(pop.rs, size=numbersamples, replace = T)
  # Find sample sizes for the studies and the numbers needed for indirect selection
  numberpeople <- sample(smallN, size = numbersamples, replace = T)
  n.unrestrict <- round(numberpeople/IRR)
  # fill the population matrix to find the data for each study
  for (i in 1:numbersamples){   # inner loop: do the folowing for number of samples for a meta
    rST <- RhoST      #rhobox [2,1]-correlation between S and T
    rTP <- rho[i]        #[3,2]
    rTx <- sqrt(rxxa[i])  #[4,2]
    rPy <- sqrt(ryya[i])  #[5,3]
    rSP <- rST*rTP     #[3,1]
    rSx <- rST*rTx     #[4,1]
    rSy <- rST*rTP*rPy #[5,1]
    rTy <- rTP*rPy     #[5,2]
    rPx <- rTx*rTP     #[4,3]
    rxy <- rTP*rTx*rPy #[5,4]
    #
    rhobox[2,1] <- rST
    rhobox[3,1] <- rSP
    rhobox[4,1] <- rSx
    rhobox[5,1] <- rSy
    rhobox[3,2] <- rTP
    rhobox[4,2] <- rTx
    rhobox[5,2] <- rTy
    rhobox[4,3] <- rPx
    rhobox[5,3] <- rPy
    rhobox[5,4] <- rxy
    #
    rhobox[1,2] <- rhobox[2,1]
    rhobox[1,3] <- rhobox[3,1]
    rhobox[1,4] <- rhobox[4,1]
    rhobox[1,5] <- rhobox[5,1]
    rhobox[2,3] <- rhobox[3,2]
    rhobox[2,4] <- rhobox[4,2]
    rhobox[2,5] <- rhobox[5,2]
    rhobox[3,4] <- rhobox[4,3]
    rhobox[3,5] <- rhobox[5,3]
    rhobox[4,5] <- rhobox[5,4]
    # unrestricted sample
    rTP.pop[i] <- rho[i]
    #
    sample.rs <-mvrnorm(n.unrestrict[i],mu=mu,Sigma=rhobox) # sample ni data for obs r from attenuated rho
    cor1 <-cor(sample.rs)                           # compute unrestricted correlation from sampled data
    ri.un[i] <- cor1[5,4]                            # save the sampled correlation for output
    rxxa.pop[i] <- rhobox[4,2]                                 # save the pop reliability of x
    ryya.pop[i] <- rhobox[5,3]                                 # save the pop reliability of y
    rxx.un[i] <- cor1[4,2]                            # save sample sqrt(rel.x)
    ryy.un[i] <- cor1[5,3]                            # save sample sqrt(rel.y)
    sdX.un <- sd(sample.rs[,4])                     # save SD(x) unrestricted for RR calcuations
    #
    # Indirect range restriction - select on S
    #
    sample2 <- sample.rs[order(-sample.rs[,1]),]  # column 1 is S
    nsubset <- round(n.unrestrict[i]*IRR)    # find number of people to select top down
    sample3 <- sample2[1:nsubset, ]               # select top (prop.select) percent
    sdX.ind <- sd(sample3[,4])                    # find SD for X in selected group
    num.meta[i] <- nsubset                        # number of people in study for meta
    # compute correlation for indirect range restricted people
    corr.indir <- cor(sample3)
    colnames(corr.indir) <- rhonames
    rownames(corr.indir) <- rhonames
    ri.ind[i] <- corr.indir[5,4]                          # rXY with indirect RR
    rxx.ind[i] <- corr.indir[4,2]                         # rxx with indirect RR
    ryy.ind[i] <- corr.indir[5,3]                         # ryy with indirect RR
    uX.ind[i] <-sdX.ind/sdX.un                            # sample range restriction in X
    #
  } # end inner loop
  
  
  # we have 3 sets of numbers
  # 1.  the population values of ri, rxx, ryy, and RR
  # 2.  the unrestricted correlations and values of rxx, ryy, and u
  # 3.  correlations, rxx, ryy and u for indirect R
  ####################  Part 3: Compute a meta-analysis on simulated data #########################
  ri <-ri.ind
  rxxi <- rxx.ind^2
  ryyi <- ryy.ind^2
  ni <- numberpeople
  k <- length(ni)
  ui <- uX.ind
  bootdata <- data.frame(ri,rxxi, ryyi, ui, ni)
  
  ########################################
  #  Computations of H&S indirect RR
  #
  ########################################
  nR <- ri*ni                                 # weight r by N
  Nsum <- sum(ni)                             # find sum of N
  rbar <- sum(nR)/Nsum                        # find sample-weighted mean r
  V.obs <- sum(ni*(ri-rbar)^2)/(sum(ni))      # find weighted observed variance
  RXXa <-1-ui^2*(1-rxxi)                      # find the reliability of the unrestricted sample
  UT <- 1/sqrt((ui^2-(1-RXXa))/RXXa)              # find UT, disattenuation for indirect RR
  rdiss.rel <- ri/sqrt(rxxi*ryyi)                 # find correlation disattenuated for reliability in X and Y
  rC.2 <-
    (rdiss.rel*UT)/(sqrt((UT^2-1)*rdiss.rel^2+1)) # find corrected correlations
  A.compound.2 <- ri/rC.2                         # find compound correction factor
  wi.2 <- ni*A.compound.2^2                       # find the weights
  rbarC.2 <- sum(wi.2*rC.2)/sum(wi.2)             # find the mean corrected correlation
  V.rC.2 <- sum(wi.2*(rC.2-rbarC.2)^2)/sum(wi.2)  # find the corrected total variance
  V.rC.k <- V.rC.2*(k/(k-1))
  V.eo.2 <- (1-rbar^2)^2/(ni-1)                   # simple observed error variance for each study
  V.ec2 <- V.eo.2/A.compound.2^2                  # find the corrected error variance
  err.adj.2 <- 1/((UT^2-1)*ri^2+1)                # find the error adjustment for indirect range restriction
  V.ve.2 <- sum(wi.2*V.ec2*err.adj.2^2)/sum(wi.2) # find the error variance for the meta
  V.rho.2 <- V.rC.2-V.ve.2                        # find the random-effects variance component
  V.rho.k <- V.rC.k-V.ve.2
  if (V.rho.2 < 0) {V.rho.2 <- 0}                 # if REVC is less than zero, set to zero
  if (V.rho.k < 0) {V.rho.k <- 0}
  SD.rho.2 <- sqrt(V.rho.2)                       # find SD rho
  SD.rho.k <- sqrt(V.rho.k)
  HS.CR80.L <- rbarC.2 - 1.28*SD.rho.2             # find the lower bound of the credibility interval
  HS.K.CR80.L <- rbarC.2 - 1.28*SD.rho.k             # lower bound, k corrected
  #
  #H&S Analytic Approach For 95% Tolerance Interval Around 10 Percentile
  HS.Percentile.Lower.CI <- rbarC.2 - lower.bound*SD.rho.2
  HS.Percentile.Upper.CI <- rbarC.2 - upper.bound*SD.rho.2
  HS.Percentile.CI.Width <- HS.Percentile.Upper.CI - HS.Percentile.Lower.CI
  HSk.Percentile.Lower.CI <- rbarC.2 - lower.bound*SD.rho.k
  HSk.Percentile.Upper.CI <- rbarC.2 - upper.bound*SD.rho.k
  HSk.Percentile.CI.Width <- HSk.Percentile.Upper.CI - HSk.Percentile.Lower.CI
  
  ######Bootstrapping Approach Using Hunter-Schmidt Credibility Interval
  #Function that produces Hunter-Schmit Lower-bound of Credibility Interval
  HS.f <- function(d, i){
    d2 <- d[i,]
    boot.ri <- d2$ri
    k <- length(boot.ri)
    boot.ni <- d2$ni
    boot.rxxi <-d2$rxxi
    boot.ryyi <-d2$ryyi
    boot.ui <- d2$ui
    b.nR <- boot.ri*boot.ni                                    # weight r by N
    b.Nsum <- sum(boot.ni)                                     # find sum of N
    b.rbar <- sum(b.nR)/b.Nsum                                 # find sample-weighted mean r
    b.V.obs <- sum(boot.ni*(boot.ri-b.rbar)^2)/(sum(boot.ni))  # find weighted observed variance
    b.RXXa <-1-boot.ui^2*(1-boot.rxxi)                         # find the reliability of the unrestricted sample
    b.UT <- 1/sqrt((boot.ui^2-(1-b.RXXa))/b.RXXa)              # find UT, disattenuation for indirect RR
    b.rdiss.rel <- boot.ri/sqrt(boot.rxxi*boot.ryyi)           # find correlation disattenuated for reliability in X and Y
    b.rC.2 <-
      (b.rdiss.rel*b.UT)/(sqrt((b.UT^2-1)*b.rdiss.rel^2+1))    # find corrected correlations
    b.A.compound.2 <- boot.ri/b.rC.2                           # find compound correction factor
    b.wi.2 <- boot.ni*b.A.compound.2^2                         # find the weights
    b.rbarC.2 <- sum(b.wi.2*rC.2)/sum(b.wi.2)                  # find the mean corrected correlation
    b.V.rC.2 <- sum(b.wi.2*(b.rC.2-b.rbarC.2)^2)/sum(b.wi.2)   # find the corrected total variance
    b.V.eo.2 <- (1-b.rbar^2)^2/(boot.ni-1)                     # simple observed error variance for each study
    b.V.ec2 <- b.V.eo.2/b.A.compound.2^2                       # find the corrected error variance
    b.err.adj.2 <- 1/((b.UT^2-1)*boot.ri^2+1)                  # find the error adjustment for indirect range restriction
    b.V.ve.2 <- sum(b.wi.2*b.V.ec2*b.err.adj.2^2)/sum(b.wi.2)  # find the error variance for the meta
    b.V.rho.2 <- b.V.rC.2-b.V.ve.2                             # find the random-effects variance component
    if (b.V.rho.2 < 0) {b.V.rho.2 <- 0}                        # if REVC is less than zero, set to zero
    b.SD.rho.2 <- sqrt(b.V.rho.2)                              # find SD rho
    return(b.rbarC.2 - 1.28*b.SD.rho.2)                        # find the lower bound of the credibility interval
  } # end HS bootstrapped lower bound of credibility interval funtion
  
  #Bootstrapped Statistic
  HS.bootstrap <- boot(bootdata, HS.f, R = 2000)
  HS.bootstrap1 <- HS.bootstrap$t0
  #Confidence Interval
  HS.bootstrap2 <- boot.ci(HS.bootstrap, type="bca")
  HS.bootstrapLB <- HS.bootstrap2$bca[4]
  HS.bootstrapUB <- HS.bootstrap2$bca[5]
  HS.bootstrap.CI.Width <- HS.bootstrapUB-HS.bootstrapLB
  
  ######Bootstrapping Approach Using Hunter-Schmidt K-corrected Credibility Interval
  #Function that produces Hunter-Schmit K-corrected Lower-bound of Credibility Interval
  HS.K.f <- function(d, i){
    d2 <- d[i,]
    boot.ri <- d2$ri
    k <- length(boot.ri)
    boot.ni <- d2$ni
    boot.rxxi <-d2$rxxi
    boot.ryyi <-d2$ryyi
    boot.ui <- d2$ui
    b.nR <- boot.ri*boot.ni                                    # weight r by N
    b.Nsum <- sum(boot.ni)                                     # find sum of N
    b.rbar <- sum(b.nR)/b.Nsum                                 # find sample-weighted mean r
    b.V.obs <- sum(boot.ni*(boot.ri-b.rbar)^2)/(sum(boot.ni))  # find weighted observed variance
    b.RXXa <-1-boot.ui^2*(1-boot.rxxi)                         # find the reliability of the unrestricted sample
    b.UT <- 1/sqrt((boot.ui^2-(1-b.RXXa))/b.RXXa)              # find UT, disattenuation for indirect RR
    b.rdiss.rel <- boot.ri/sqrt(boot.rxxi*boot.ryyi)           # find correlation disattenuated for reliability in X and Y
    b.rC.2 <-
      (b.rdiss.rel*b.UT)/(sqrt((b.UT^2-1)*b.rdiss.rel^2+1))    # find corrected correlations
    b.A.compound.2 <- boot.ri/b.rC.2                           # find compound correction factor
    b.wi.2 <- boot.ni*b.A.compound.2^2                         # find the weights
    b.rbarC.2 <- sum(b.wi.2*rC.2)/sum(b.wi.2)                  # find the mean corrected correlation
    b.V.rC.2 <- sum(b.wi.2*(b.rC.2-b.rbarC.2)^2)/sum(b.wi.2)   # find the corrected total variance
    b.V.rC.2 <- b.V.rC.2*(k/(k-1))
    b.V.eo.2 <- (1-b.rbar^2)^2/(boot.ni-1)                     # simple observed error variance for each study
    b.V.ec2 <- b.V.eo.2/b.A.compound.2^2                       # find the corrected error variance
    b.err.adj.2 <- 1/((b.UT^2-1)*boot.ri^2+1)                  # find the error adjustment for indirect range restriction
    b.V.ve.2 <- sum(b.wi.2*b.V.ec2*b.err.adj.2^2)/sum(b.wi.2)  # find the error variance for the meta
    b.V.rho.2 <- b.V.rC.2-b.V.ve.2                             # find the random-effects variance component
    if (b.V.rho.2 < 0) {b.V.rho.2 <- 0}                        # if REVC is less than zero, set to zero
    b.SD.rho.2 <- sqrt(b.V.rho.2)                              # find SD rho
    return(b.rbarC.2 - 1.28*b.SD.rho.2)                        # find the lower bound of the credibility interval
  } # end HS bootstrapped lower bound of credibility interval funtion
  
  #Bootstrapped Statistic
  HS.K.bootstrap <- boot(bootdata, HS.K.f, R = 2000)
  HS.K.bootstrap1 <- HS.K.bootstrap$t0
  #Confidence Interval
  HS.K.bootstrap2 <- boot.ci(HS.K.bootstrap, type="bca")
  HS.K.bootstrapLB <- HS.K.bootstrap2$bca[4]
  HS.K.bootstrapUB <- HS.K.bootstrap2$bca[5]
  HS.K.bootstrap.CI.Width <- HS.K.bootstrapUB-HS.K.bootstrapLB
  #############  Morris ###########################################################
  morr.ri <- rC.2                                    # fully corrected correlations
  morr.V.eo <- V.ec2/err.adj.2^2                     # corrected error variances
  morris.dat <- data.frame(cbind(morr.ri,morr.V.eo)) # collect  estimates
  morris1 <- rma(yi=morr.ri,vi=morr.V.eo,data=morris.dat,
                 control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
  Morris.M.rho <- morris1$b                              # output the mean
  rownames(Morris.M.rho) <- c()                          # strips value of intercept label
  colnames(Morris.M.rho) <- c("Morris.M.rho")
  Morris.CI95.L <- morris1$ci.lb                         # output the lower CI bound
  Morris.CI95.U <- morris1$ci.ub                         # output the upper CI bound
  Morris.V.rho <- morris1$tau2                           # random-effects variance component
  Morris.SD.rho <- sqrt(morris1$tau2)                    # output for Morris RHO sd
  Morris.CR80.L <- (Morris.M.rho - 1.28*Morris.SD.rho)   # Lower CR Bound
  Morris.CR80.U <- (Morris.M.rho + 1.28*Morris.SD.rho)   # Upper CR Bound
  
  #Analytic Approach For 95% Tolerance Interval Around 10 Percntile
  M <- Morris.M.rho
  sd <- Morris.SD.rho
  
  Morris.Percentile.Lower.CI <- M - lower.bound*sd
  Morris.Percentile.Middle.CI <- M - middle.bound*sd
  Morris.Percentile.Upper.CI <- M - upper.bound*sd
  Morris.Percentile.CI.Width <- Morris.Percentile.Upper.CI-Morris.Percentile.Lower.CI
  
  ######Bootstrapping Approach Using Morris Credibility Interval
  #Function that produces Morris Lower-bound of Credibility Interval
  Morris.f <- function(d, i){
    d2 <- d[i,]
    boot.ri <- d2$ri
    k <- length(boot.ri)
    boot.ni <- d2$ni
    boot.rxxi <-d2$rxxi
    boot.ryyi <-d2$ryyi
    boot.ui <- d2$ui
    b.nR <- boot.ri*boot.ni                                    # weight r by N
    b.Nsum <- sum(boot.ni)                                     # find sum of N
    b.rbar <- sum(b.nR)/b.Nsum                                 # find sample-weighted mean r
    b.RXXa <-1-boot.ui^2*(1-boot.rxxi)                         # find the reliability of the unrestricted sample
    b.UT <- 1/sqrt((boot.ui^2-(1-b.RXXa))/b.RXXa)              # find UT, disattenuation for indirect RR
    b.rdiss.rel <- boot.ri/sqrt(boot.rxxi*boot.ryyi)           # find correlation disattenuated for reliability in X and Y
    b.rC.2 <-
      (b.rdiss.rel*b.UT)/(sqrt((b.UT^2-1)*b.rdiss.rel^2+1))    # find corrected correlations
    b.A.compound.2 <- boot.ri/b.rC.2                           # find compound correction facto
    b.V.eo.2 <- (1-b.rbar^2)^2/(boot.ni-1)                     # simple observed error variance for each study
    b.V.ec2 <- b.V.eo.2/b.A.compound.2^2                       # find the corrected error variance
    b.err.adj.2 <- 1/((b.UT^2-1)*boot.ri^2+1)                  # find the error adjustment for indirect range restriction
    b.var.i <- b.V.ec2/b.err.adj.2 #
    morris.dat <- data.frame(cbind(b.rC.2,b.var.i)) # collect  estimates
    morris1 <- rma(yi=b.rC.2,vi=b.var.i,data=morris.dat,
                   control=list(maxiter=1000, stepadj=.5), method="REML")        # run the random-effects meta with REML
    Morris.M.rho <- morris1$b                              # output the mean
    rownames(Morris.M.rho) <- c()                          # strips value of intercept label
    colnames(Morris.M.rho) <- c("Morris.M.rho")
    #    Morris.CI95.L <- morris1$ci.lb                         # output the lower CI bound
    #    Morris.CI95.U <- morris1$ci.ub                         # output the upper CI bound
    #    Morris.V.rho <- morris1$tau2                           # random-effects variance component
    Morris.SD.rho <- sqrt(morris1$tau2)
    return(Morris.M.rho - 1.28*Morris.SD.rho)
  }
  #Bootstrapped Statistic
  Morris.bootstrap <- boot(bootdata, Morris.f, R = 2000)
  Morris.bootstrap1 <- Morris.bootstrap$t0
  #Confidence Interval
  Morris.bootstrap2 <- boot.ci(Morris.bootstrap, type="bca")
  Morris.bootstrapLB <- Morris.bootstrap2$bca[4]
  Morris.bootstrapUB <- Morris.bootstrap2$bca[5]
  Morris.bootstrap.CI.Width <- Morris.bootstrapUB-Morris.bootstrapLB
  
  ###Coverage Check###
  
  HS.Percentile.coverage <-0
  if(HS.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= HS.Percentile.Upper.CI) {HS.Percentile.coverage<-1}
  HS.bootstrap.coverage <-0
  if(HS.bootstrapLB<= true.lower.bound & true.lower.bound<=  HS.bootstrapUB) {HS.bootstrap.coverage<-1}
  HS.K.Percentile.coverage <-0
  if(HSk.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= HSk.Percentile.Upper.CI) {HS.K.Percentile.coverage<-1}
  HS.K.bootstrap.coverage <-0
  if(HS.K.bootstrapLB<= true.lower.bound & true.lower.bound<=  HS.K.bootstrapUB) {HS.K.bootstrap.coverage<-1}
  Morris.Percentile.coverage <-0
  if(Morris.Percentile.Lower.CI<= true.lower.bound & true.lower.bound<= Morris.Percentile.Upper.CI) {Morris.Percentile.coverage<-1}
  Morris.bootstrap.coverage <-0
  if(Morris.bootstrapLB<= true.lower.bound & true.lower.bound<= Morris.bootstrapUB) {Morris.bootstrap.coverage<-1}
  
  ###Check if REVC Estimated is greater than 0###
  Non.Zero.HS.REVC <- 0
  if(V.rho.2>0) {Non.Zero.HS.REVC<-1}
  Non.Zero.HS.K.REVC <- 0
  if(V.rho.k>0) {Non.Zero.HS.K.REVC<-1}
  Non.Zero.Morris.REVC <- 0
  if(Morris.V.rho>0) {Non.Zero.Morris.REVC<-1}
  
  # Check if REVC Is Greater than 0 but Coverage is 0
  # HS.PercentileNonZeroREVC.NoCoverage <- 0
  # if(Non.Zero.HS.REVC <= 0) {HS.PercentileNonZeroREVC.NoCoverage<-NA}
  # if(Non.Zero.HS.REVC >= 1 & HS.Percentile.coverage<=0) {HS.PercentileNonZeroREVC.NoCoverage<-1}
  # HS.K.PercentileNonZeroREVC.NoCoverage <- 0
  # if(Non.Zero.HS.K.REVC <= 0) {HS.K.PercentileNonZeroREVC.NoCoverage<-NA}
  # if(Non.Zero.HS.K.REVC >= 1 & HS.K.Percentile.coverage<=0) {HS.K.PercentileNonZeroREVC.NoCoverage<-1}
  # Morris.PercentileNonZeroREVC.NoCoverage <- 0
  # if(Non.Zero.Morris.REVC <= 0) {Morris.PercentileNonZeroREVC.NoCoverage<-NA}
  # if(Non.Zero.Morris.REVC >= 1 & Morris.Percentile.coverage<=0) {Morris.PercentileNonZeroREVC.NoCoverage<-1}
  
  ##########################################
  # Output results
  ##########################################
  Results[j, ] <- cbind(rho.mean, sigmarho, true.lower.bound, # 1, 2, 3
                        HS.CR80.L, HS.K.CR80.L, Morris.CR80.L, # 4, 5, 6
                        HS.Percentile.coverage, HS.bootstrap.coverage,         # 7, 8
                        HS.K.Percentile.coverage,HS.K.bootstrap.coverage,      # 9, 10
                        Morris.Percentile.coverage, Morris.bootstrap.coverage, # 11, 12
                        HS.Percentile.Lower.CI, HS.Percentile.Upper.CI,         # 13, 14
                        HS.bootstrapLB, HS.bootstrapUB,                         # 15, 16
                        HSk.Percentile.Lower.CI, HSk.Percentile.Upper.CI,     # 17, 18
                        HS.K.bootstrapLB, HS.K.bootstrapUB,                     # 19, 20
                        Morris.Percentile.Lower.CI, Morris.Percentile.Upper.CI, # 21, 22
                        Morris.bootstrapLB, Morris.bootstrapUB,                 # 23, 24
                        rbarC.2, SD.rho.2, SD.rho.k,                            # 25, 26, 27
                        Morris.M.rho, Morris.SD.rho,                            # 28, 29
                        Non.Zero.HS.REVC, Non.Zero.HS.K.REVC, Non.Zero.Morris.REVC, #30, 31, 32
                        HS.Percentile.CI.Width,HS.bootstrap.CI.Width,HSk.Percentile.CI.Width, #33, 34, 35
                        HS.K.bootstrap.CI.Width, Morris.Percentile.CI.Width,Morris.bootstrap.CI.Width) #36,37,38
} # END outer loop.
# Rename the columns of the Results
resultsnames = c("rho", "sigmarho", "true lower bound",                         #1, 2, 3
                 "HS Lower Bound", "HS.K Lower Bound",  "Morris Lower Bound",   #4, 5, 6
                 "HS Percentile Coverage", "HS Bootstrap Coverage",             # 7, 8
                 "HS.K Percentile Coverage", "HS.K Bootstrap Coverage",         # 9, 10
                 "Morris Percentile Coverage", "Morris Bootstrap Coverage",     #11, 12
                 "HS Percentile LB", "HS Percentile UB",                        #13, 14
                 "HS Bootstrap LB", "HS Bootstrap UB",                          #15, 16
                 "HS K Percentile LB", "HS K Percentile UB",                    #17, 18
                 "HS K Bootstrap LB", "HS K Bootstrap UB",                      #19, 20
                 "Morris Percentile LB", "Morris Percentile UB",                #21, 22
                 "Morris Bootstrap LB", "Morris Bootstrap UB",                 #23, 24
                 "HSMean", "HSSDrho", "HSKSDrho",                              #25, 26, 27
                 "MorrisMean", "MorrisSDrho",                                  #28, 29
                 "HSproNonZero", "HSKpropNonZer", "MorrisproNonZero",          #30, 31, 32
                 "HS Percentile CI Width", "HS K Percentile CI Width", "HS bootstrap CI Width", #33, 34, 35
                 "HS K bootstrap CI Width", "Morris Percentile CI Width", "Morris bootstrap CI Width") #36, 37 ,38
colnames(Results) = resultsnames

# Calculate the biases for the 3 estimators.

HS.Bias <- sum(Results[,4]-Results[,3])/numbermetas
HS.K.Bias <-sum(Results[,5]-Results[,3])/numbermetas
Morris.Bias <-sum(Results[,6]-Results[,3])/numbermetas

# Calculate the RMSE for the 3 estimators.

HS.RMSE <- sqrt(sum((Results[,4]-Results[,3])^2)/numbermetas)
HS.K.RMSE <- sqrt(sum((Results[,5]-Results[,3])^2)/numbermetas)
Morris.RMSE <- sqrt(sum((Results[,6]-Results[,3])^2)/numbermetas)

#Calculate the average confidence interval estimate

HS.Percentile.AvgLB <- sum(Results[,13])/numbermetas
HS.Percentile.AvgUB <- sum(Results[,14])/numbermetas
HS.Bootstrap.AvgLB <- sum(Results[,15])/numbermetas
HS.Bootstrap.AvgUB <- sum(Results[,16])/numbermetas
HS.K.Percentile.AvgLB <- sum(Results[,17])/numbermetas
HS.K.Percentile.AvgUB <- sum(Results[,18])/numbermetas
HS.K.Bootstrap.AvgLB <- sum(Results[,19])/numbermetas
HS.K.Bootstrap.AvgUB <- sum(Results[,20])/numbermetas
Morris.Percentile.AvgLB <-sum(Results[,21])/numbermetas
Morris.Percentile.AvgUB <-sum(Results[,22])/numbermetas
Morris.Bootstrap.AvgLB <-sum(Results[,23])/numbermetas
Morris.Bootstrap.AvgUB <-sum(Results[,24])/numbermetas

#Calculate the coverage for the 8 types of estimates

HS.Percentile.coverage <-mean(Results[,7])
HS.Bootstrap.coverage <-mean(Results[,8])
HS.K.Percentile.coverage <-mean(Results[,9])
HS.K.Bootstrap.coverage <-mean(Results[,10])
Morris.Percentile.coverage <-mean(Results[,11])
Morris.Bootstrap.coverage <-mean(Results[,12])
HS.Mean <- mean(Results[,25])
Morris.Mean <- mean(Results[,28])
HS.SDrho <- mean(Results[,26])
HS.K.SDrho <- mean(Results[,27])
Morris.M.SDrho <- mean(Results[,29])
#Calculate Proportion of Times REVC Estimated to be Greater than 0
Nonzero.HS.REVC <- mean(Results[,30])
Nonzero.HS.K.REVC <- mean(Results[,31])
Nonzero.Morris.REVC <- mean(Results[,32])

#Calculate Proportion of Times Coverage is 0 When REVC Is greater than 0
#HS.Percentile.NonZeroREVC.NoCoverage <-mean(Results[,48], na.rm = TRUE)
#HS.K.Percentile.NonZeroREVC.NoCoverage <-mean(Results[,49], na.rm = TRUE)
#Morris.Percentile.NonZeroREVC.NoCoverage <-mean(Results[,50], na.rm = TRUE)


#Caclulate Proption of Times REVC Greater than 0 but Coverage is 0


Output <- data.frame(rho.mean, sigmarho, numbersamples, RhoST, IRR,
                     HS.Bias, HS.K.Bias, Morris.Bias,
                     HS.RMSE, HS.K.RMSE, Morris.RMSE,
                     HS.Percentile.coverage, HS.Bootstrap.coverage,
                     HS.K.Percentile.coverage, HS.K.Bootstrap.coverage,
                     Morris.Percentile.coverage, Morris.Bootstrap.coverage,
                     HS.Percentile.AvgLB, HS.Percentile.AvgUB, HS.Bootstrap.AvgLB, HS.Bootstrap.AvgUB,
                     HS.K.Percentile.AvgLB, HS.K.Percentile.AvgUB, HS.K.Bootstrap.AvgLB, HS.K.Bootstrap.AvgUB,
                     Morris.Percentile.AvgLB, Morris.Percentile.AvgUB,
                     Morris.Bootstrap.AvgLB, Morris.Bootstrap.AvgUB,HS.Mean,Morris.Mean, HS.SDrho,
                     HS.K.SDrho, Morris.M.SDrho, Nonzero.HS.REVC, Nonzero.HS.K.REVC, Nonzero.Morris.REVC)

Output

# THIS WILL EXPORT RESULTS TO WORKING DIRECTORY 
# NEED TO CHANGE RESULTS TO CONDITION NUMBER; i.e condition saved as results1.csv 
write.table(Results, file = "Results.csv", row.names=F, append=T, col.names=T, sep=",")

# THIS WILL EXPORT Output TO WORKING DIRECTORY 
write.table(Output, file = "Output.csv", row.names=F, append=T, col.names=T, sep=",")
