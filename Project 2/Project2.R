setwd("C:/Users/flode/OneDrive - Universite de Liege/Bureau/HDDA/Project 2 - 14")

## Generation of the data
mu1 = rep(0,12)
mu2 = c(rep(0,9), rep(4,3))
P9 = matrix(0.75, nrow=9, ncol=9); diag(P9) <- 1
P3 = matrix(0.75, nrow=3, ncol=3); diag(P3) <- 1
sigma = rbind(cbind(P9, matrix(0, nrow=9, ncol=3)), cbind(matrix(0, nrow=3, ncol=9), P3))

n = 150 # sample size
epsilon = c(0, 0.1, 0.2)
Nsim = 500

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Q 3.1 Outlier Detection
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

rob_PTP <- matrix(NA, nrow=Nsim, ncol=3)
cla_PTP <- matrix(NA, nrow=Nsim, ncol=3)
rob_PFP <- matrix(NA, nrow=Nsim, ncol=3)
cla_PFP <- matrix(NA, nrow=Nsim, ncol=3)

set.seed(42)
for (i in 1:Nsim)
{
  # For epsilon = 0 we just need to compute for x1 (normal distribution with mean 0)
  library(mvtnorm)
  x1 <- rmvnorm(n, mu1, sigma)
  d_mahalanobis <- mahalanobis(x1, colMeans(x1), cov(x1))
  
  library(MASS)
  h <- floor(n*0.75)
  rob_est <- cov.rob(x1, quantile.used=h, method="mcd", cor=FALSE)
  rob_dist <- mahalanobis(x1, rob_est$center, rob_est$cov)
  
  # The PTP is 0 since there is no sample in x2 
  # In this part, observations coming from x2, which is centered at mu2, are considered as being corrupted
  
  # We compute the PFP, prop of samples in x1 detected as outliers
  rob_PFP[i,1] <- 100*sum(rob_dist > qchisq(0.95,12))/n
  cla_PFP[i,1] <- 100*sum(d_mahalanobis > qchisq(0.95,12))/n
  
  # Now we compute everything for each of the other values of epsilon, but we have
  # 2 different groups x1 and x2.
  for (j in 2:3)
  {
    weights <- c(1-epsilon[j], epsilon[j])
    type <- sample(1:2, size=n, replace=TRUE, prob=weights)
    ni <- table(type)
    
    library(mvtnorm)
    x1 <- rmvnorm(ni[1], mu1, sigma)
    x2 <- rmvnorm(ni[2], mu2, sigma)
    F <- rbind(x1, x2)
    
    # Classical mahalanobis distance
    d_mahalanobis <- mahalanobis(F, colMeans(F), cov(F))
    
    # Robust distance with covariance matrix estimated by the mcd estimator
    library(MASS)
    h <- floor(n*0.75)
    rob_est <- cov.rob(F, quantile.used=h, method="mcd", cor=FALSE)
    rob_dist <- mahalanobis(F, rob_est$center, rob_est$cov)
    
    # We compute the PTP, prop of samples in x2 detected as outliers
    # NB : the chi-square cutoff relies on the normality assumption of the data
    rob_PTP[i,j] <- 100*sum(rob_dist[(ni[1]+1):n] > qchisq(0.95,12))/ni[2]
    cla_PTP[i,j] <- 100*sum(d_mahalanobis[(ni[1]+1):n] > qchisq(0.95,12))/ni[2]
    
    # We compute the PFP, prop of samples in x1 detected as outliers
    rob_PFP[i,j] <- 100*sum(rob_dist[1:ni[1]] > qchisq(0.95,12))/ni[1]
    cla_PFP[i,j] <- 100*sum(d_mahalanobis[1:ni[1]] > qchisq(0.95,12))/ni[1]
  }
}

avg_rob_PTP <- colMeans(rob_PTP, na.rm=TRUE)
avg_cla_PTP <- colMeans(cla_PTP, na.rm=TRUE)
avg_rob_PFP <- colMeans(rob_PFP, na.rm=TRUE)
avg_cla_PFP <- colMeans(cla_PFP, na.rm=TRUE)
sd_rob_PTP <- sapply(as.data.frame(rob_PTP),sd)
sd_cla_PTP <- sapply(as.data.frame(cla_PTP),sd)
sd_rob_PFP <- sapply(as.data.frame(rob_PFP),sd)
sd_cla_PFP <- sapply(as.data.frame(cla_PFP),sd)

par(mfrow=c(1,3))
for (i in 1:3)
{
  data <- data.frame(x1 = cla_PTP[,i], x2 = rob_PTP[,i], x3 = cla_PFP[,i], x4 = rob_PFP[,i])
  colnames(data) <- c("cla_PTP", "rob_PTP", "cla_PFP", "rob_PFP")
  boxplot(data, main = paste("boxplot for epsilon =", epsilon[i]),ylim=c(0,100))
}

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Q 3.2 Outlier detection after dimension reduction
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

n = 150
e = 0.2
weights <- c(1-e, e)

set.seed(42) ; type <- sample(1:2, size=n, replace=TRUE, prob=weights)
ni <- table(type)

set.seed(42)
x1 <- rmvnorm(ni[1], mu1, sigma)
x2 <- rmvnorm(ni[2], mu2, sigma)
colVec <- c(rep("#FF0000", ni[1]), rep("#0000FF", ni[2]))
F <- rbind(x1, x2)

# ----------
# Q.3.2.1 a) - tSNE

library(Rtsne)
set.seed(42) ; restSNE1_FF <- Rtsne(X=F, perplexity=1, pca=FALSE, normalize=FALSE)
set.seed(42) ; restSNE2_FF <- Rtsne(X=F, perplexity=2, pca=FALSE, normalize=FALSE)
set.seed(42) ; restSNE5_FF <- Rtsne(X=F, perplexity=5, pca=FALSE, normalize=FALSE)
set.seed(42) ; restSNE10_FF <- Rtsne(X=F, perplexity=10, pca=FALSE, normalize=FALSE)
set.seed(42) ; restSNE30_FF <- Rtsne(X=F, perplexity=30, pca=FALSE, normalize=FALSE)

par(mfrow=c(2,3))
plot(restSNE1_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=1")
plot(restSNE2_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=2")
plot(restSNE5_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=5")
plot(restSNE10_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=10")
plot(restSNE30_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=30")
par(mfrow=c(1,1))

# --------------------
set.seed(42) ; restSNE1_TF <- Rtsne(X=F, perplexity=1, pca=TRUE, normalize=FALSE)
set.seed(42) ; restSNE2_TF <- Rtsne(X=F, perplexity=2, pca=TRUE, normalize=FALSE)
set.seed(42) ; restSNE5_TF <- Rtsne(X=F, perplexity=5, pca=TRUE, normalize=FALSE)
set.seed(42) ; restSNE10_TF <- Rtsne(X=F, perplexity=10, pca=TRUE, normalize=FALSE)
set.seed(42) ; restSNE30_TF <- Rtsne(X=F, perplexity=30, pca=TRUE, normalize=FALSE)

par(mfrow=c(2,3))
plot(restSNE1_TF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=1")
plot(restSNE2_TF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=2")
plot(restSNE5_TF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=5")
plot(restSNE10_TF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=10")
plot(restSNE30_TF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=30")
par(mfrow=c(1,1))

# --------------------
set.seed(42) ; restSNE1_FT <- Rtsne(X=F, perplexity=1, pca=FALSE, normalize=TRUE)
set.seed(42) ; restSNE2_FT <- Rtsne(X=F, perplexity=2, pca=FALSE, normalize=TRUE)
set.seed(42) ; restSNE5_FT <- Rtsne(X=F, perplexity=5, pca=FALSE, normalize=TRUE)
set.seed(42) ; restSNE10_FT <- Rtsne(X=F, perplexity=10, pca=FALSE, normalize=TRUE)
set.seed(42) ; restSNE30_FT <- Rtsne(X=F, perplexity=30, pca=FALSE, normalize=TRUE)

par(mfrow=c(2,3))
plot(restSNE1_FT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=1")
plot(restSNE2_FT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=2")
plot(restSNE5_FT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=5")
plot(restSNE10_FT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=10")
plot(restSNE30_FT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=30")
par(mfrow=c(1,1))

# --------------------
set.seed(42) ; restSNE1_TT <- Rtsne(X=F, perplexity=1, pca=TRUE, normalize=TRUE)
set.seed(42) ; restSNE2_TT <- Rtsne(X=F, perplexity=2, pca=TRUE, normalize=TRUE)
set.seed(42) ; restSNE5_TT <- Rtsne(X=F, perplexity=5, pca=TRUE, normalize=TRUE)
set.seed(42) ; restSNE10_TT <- Rtsne(X=F, perplexity=10, pca=TRUE, normalize=TRUE)
set.seed(42) ; restSNE30_TT <- Rtsne(X=F, perplexity=30, pca=TRUE, normalize=TRUE)

par(mfrow=c(2,3))
plot(restSNE1_TT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=1")
plot(restSNE2_TT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=2")
plot(restSNE5_TT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=5")
plot(restSNE10_TT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=10")
plot(restSNE30_TT$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main="Perplexity=30")
par(mfrow=c(1,1))

# Plotting graph for epsilon = 0; 0.1; 0.2 with perplexity = 10 (pca=FALSE and normalize=FALSE)
n = 150
e = c(0,0.1,0.2)
par(mfrow=c(1,3))
for (i in 1:3)
{
  weights <- c(1-e[i], e[i])
  set.seed(42) ; type <- sample(1:2, size=n, replace=TRUE, prob=weights)
  ni <- table(type)
  set.seed(42)
  x1 <- rmvnorm(ni[1], mu1, sigma)
  if(e[i]!=0)
  {
    x2 <- rmvnorm(ni[2], mu2, sigma)
    F_dis <- rbind(x1, x2)
  }
  else
  {
    F_dis <- x1
    ni[2] <- 0
  }
  colVec <- c(rep("#FF0000", ni[1]), rep("#0000FF", ni[2]))
  set.seed(42) ; restSNE10_FF <- Rtsne(X=F_dis, perplexity=10, pca=FALSE, normalize=FALSE)
  plot(restSNE10_FF$Y, col=colVec, xlab="tSNE1", ylab="tSNE2", main=paste("Perplexity=10 | Epsilon=", e[i],sep=""))
}
par(mfrow=c(1,1))

# ----------
# Q.3.2.1 b) - PCA

cov(F)
cor(F)
sapply(as.data.frame(F),sd)
res_cov <- princomp(F)
summary(res_cov)
res_cov$loadings

par(mfrow=c(3,1))
for(i in 1:3)
{
  barplot(res_cov$loadings[,i], main=paste("comp",i))
}
par(mfrow=c(1,1))

plot(res_cov$scores[,1],res_cov$scores[,2], col=colVec, main="PCA on covariance matrix", xlab="PC1",ylab="PC2")

# ----------
# Q.3.2.2 - Simulation
# We selected tSNE method
# NB : Pay attention to use chi-square cutoff with 2 dof because of the dimension reduction

Nsim <- 100
epsilon = c(0,0.1,0.2)
rob_PTP <- matrix(NA, nrow=Nsim, ncol=3)
cla_PTP <- matrix(NA, nrow=Nsim, ncol=3)
rob_PFP <- matrix(NA, nrow=Nsim, ncol=3)
cla_PFP <- matrix(NA, nrow=Nsim, ncol=3)
set.seed(42) ;
for (i in 1:Nsim)
{
  # For epsilon = 0 we just need to compute for x1 (normal distribution with mean 0)
  library(mvtnorm)
  x1 <- rmvnorm(n, mu1, sigma)
  library(Rtsne)
  restSNE10_FF <- Rtsne(X=x1, perplexity=10, pca=FALSE, normalize=FALSE)
  F_reduc <- restSNE10_FF$Y
  
  # Classical mahalanobis distance
  d_mahalanobis <- mahalanobis(F_reduc, colMeans(F_reduc), cov(F_reduc))
  
  library(MASS)
  h <- floor(n*0.75)
  rob_est <- cov.rob(F_reduc, quantile.used=h, method="mcd", cor=FALSE)
  rob_dist <- mahalanobis(F_reduc, rob_est$center, rob_est$cov)
  
  # The PTP is 0 since there is no sample in x2 
  # In this part, observations coming from x2, which is centered at mu2, are considered as being corrupted
  
  # We compute the PFP, prop of samples in x1 detected as outliers
  rob_PFP[i,1] <- 100*sum(rob_dist > qchisq(0.95,2))/n
  cla_PFP[i,1] <- 100*sum(d_mahalanobis > qchisq(0.95,2))/n
  
  # Now we compute everything for each of the other values of epsilon, but we have
  # 2 different groups x1 and x2.
  for (j in 2:3)
  {
    weights <- c(1-epsilon[j], epsilon[j])
    type <- sample(1:2, size=n, replace=TRUE, prob=weights)
    ni <- table(type)
    
    library(mvtnorm)
    x1 <- rmvnorm(ni[1], mu1, sigma)
    x2 <- rmvnorm(ni[2], mu2, sigma)
    F <- rbind(x1, x2)
    restSNE10_FF <- Rtsne(X=F, perplexity=10, pca=FALSE, normalize=FALSE)
    F_reduc <- restSNE10_FF$Y
    
    # Classical mahalanobis distance
    d_mahalanobis <- mahalanobis(F_reduc, colMeans(F_reduc), cov(F_reduc))
    
    # Robust distance with covariance matrix estimated by the mcd estimator
    library(MASS)
    h <- floor(n*0.75)
    rob_est <- cov.rob(F_reduc, quantile.used=h, method="mcd", cor=FALSE)
    rob_dist <- mahalanobis(F_reduc, rob_est$center, rob_est$cov)
    
    # We compute the PTP, prop of samples in x2 detected as outliers
    # NB : the chi-square cutoff relies on the normality assumption of the data
    rob_PTP[i,j] <- 100*sum(rob_dist[(ni[1]+1):n] > qchisq(0.95,2))/ni[2]
    cla_PTP[i,j] <- 100*sum(d_mahalanobis[(ni[1]+1):n] > qchisq(0.95,2))/ni[2]
    
    # We compute the PFP, prop of samples in x1 detected as outliers
    rob_PFP[i,j] <- 100*sum(rob_dist[1:ni[1]] > qchisq(0.95,2))/ni[1]
    cla_PFP[i,j] <- 100*sum(d_mahalanobis[1:ni[1]] > qchisq(0.95,2))/ni[1]
  }
}

avg_rob_PTP <- colMeans(rob_PTP, na.rm=TRUE)
avg_cla_PTP <- colMeans(cla_PTP, na.rm=TRUE)
avg_rob_PFP <- colMeans(rob_PFP, na.rm=TRUE)
avg_cla_PFP <- colMeans(cla_PFP, na.rm=TRUE)
sd_rob_PTP <- sapply(as.data.frame(rob_PTP),sd)
sd_cla_PTP <- sapply(as.data.frame(cla_PTP),sd)
sd_rob_PFP <- sapply(as.data.frame(rob_PFP),sd)
sd_cla_PFP <- sapply(as.data.frame(cla_PFP),sd)

par(mfrow=c(1,3))
for (i in 1:3)
{
  data <- data.frame(x1 = cla_PTP[,i], x2 = rob_PTP[,i], x3 = cla_PFP[,i], x4 = rob_PFP[,i])
  colnames(data) <- c("cla_PTP", "rob_PTP", "cla_PFP", "rob_PFP")
  boxplot(data, main = paste("boxplot for epsilon =", epsilon[i]), ylim=c(0,100))
}
par(mfrow=c(1,1))

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Q 3.3 Supervised classification on the mixture
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# ---------------
# Q 3.3.b.ii)
# ---------------

n = 150
e = c(0.1, 0.2)

#For epsilon = 0.1

weights <- c(1-e[1], e[1])

set.seed(42) ; class <- sample(1:2, size=n, replace=TRUE, prob=weights)
ni <- table(class)

set.seed(42)
library(MASS)
library(mvtnorm)
x1 <- rmvnorm(ni[1], mu1, sigma)
x2 <- rmvnorm(ni[2], mu2, sigma)
colVec <- c(rep("#FF0000", ni[1]), rep("#0000FF", ni[2]))
F <- rbind(x1, x2)


F <- as.data.frame(F)
F$Class <- c(rep(0, ni[1]), rep(1, ni[2]))

ldafull <- lda(x=F[,1:12], grouping=F$Class)
g <- 2 ; n <- dim(F)[1]
l1 <- (g-1)*(ldafull$svd)^2 / (n-g)
( gamma1 <- l1/(1+l1) )

# First we are going to try to remove one variable, the one that allows to do
# the smallest loss of discriminant power
l <- gamma_0_1 <- NULL
for(i in 1:12)
{
  var <- (1:12)[-i]
  ldaTest <- lda(x=F[,var], grouping=F$Class)
  l[i] <- (g-1)*(ldaTest$svd)^2 / (n-g)
  gamma_0_1[i] <- l[i]/(1+l[i])
}
gamma_0_1

# Seeing the vector gamma of the different discriminant powers obtained when
# leaving variable out one-by-one, we are going to be able to delete some variables

# After having tried to remove different numbers of variables, we observed that
# deleting 7 variables is a good choice to minimize the loss of discriminant power

copy_F_0_1 <- F
left_out_0_1 <- rep(0,9)

for (j in 1:9)
{
  left_out_0_1[j] <- which.max(gamma_0_1)
  copy_F_0_1 <- subset(copy_F_0_1, select = -left_out_0_1[j])
  l <- gamma_0_1 <- NULL
  for(i in 1:(12-j))
  {
    var <- (1:(12-j))[-i]
    ldaTest <- lda(x=copy_F_0_1[,var], grouping=copy_F_0_1$Class)
    l[i] <- (g-1)*(ldaTest$svd)^2 / (n-g)
    gamma_0_1[i] <- l[i]/(1+l[i])
  }
  gamma_0_1
}

# To see the variables that we didn't delete, we just have to open the dataframe 
# copy_F_0_1 and to look to the variables that are still here 
# e.g. if there is still V7, V10, V12, we know that we deleted
# V1, V2, V3, V4, V5, V6, V8, V9, V11

ldafull <- lda(x=copy_F_0_1[,1:(12-j)], grouping=copy_F_0_1$Class)
g <- 2 ; n <- dim(copy_F_0_1)[1]
l1 <- (g-1)*(ldafull$svd)^2 / (n-g)
( gamma1 <- l1/(1+l1) )

#For epsilon = 0.2

weights <- c(1-e[2], e[2])

set.seed(42) ; class <- sample(1:2, size=n, replace=TRUE, prob=weights)
ni <- table(class)

set.seed(42)
library(MASS)
library(mvtnorm)
x1 <- rmvnorm(ni[1], mu1, sigma)
x2 <- rmvnorm(ni[2], mu2, sigma)
colVec <- c(rep("#FF0000", ni[1]), rep("#0000FF", ni[2]))
F <- rbind(x1, x2)


F <- as.data.frame(F)
F$Class <- c(rep(1, ni[1]), rep(2, ni[2]))

ldafull <- lda(x=F[,1:12], grouping=F$Class)
g <- 2 ; n <- dim(F)[1]
l1 <- (g-1)*(ldafull$svd)^2 / (n-g)
( gamma1 <- l1/(1+l1) )

# Variable selection:
l <- gamma_0_2 <- NULL
for(i in 1:12)
{
  var <- (1:12)[-i]
  ldaTest <- lda(x=F[,var], grouping=F$Class)
  l[i] <- (g-1)*(ldaTest$svd)^2 / (n-g)
  gamma_0_2[i] <- l[i]/(1+l[i])
}
gamma_0_2
copy_F_0_2 <- F
left_out_0_2 <- rep(0,9)
for (j in 1:9)
{
  left_out_0_2[j] <- which.max(gamma_0_2)
  copy_F_0_2 <- subset(copy_F_0_2, select = -left_out_0_2[j])
  l <- gamma_0_2 <- NULL
  for(i in 1:(12-j))
  {
    var <- (1:(12-j))[-i]
    ldaTest <- lda(x=copy_F_0_2[,var], grouping=copy_F_0_2$Class)
    l[i] <- (g-1)*(ldaTest$svd)^2 / (n-g)
    gamma_0_2[i] <- l[i]/(1+l[i])
  }
  gamma_0_2
}
ldafull <- lda(x=copy_F_0_2[,1:(12-j)], grouping=copy_F_0_2$Class)
g <- 2 ; n <- dim(copy_F_0_2)[1]
l1 <- (g-1)*(ldafull$svd)^2 / (n-g)
( gamma1 <- l1/(1+l1) )

# -----------
# Q 3.3.c
# -----------

Nsim = 500
PTP = matrix(0, nrow=Nsim, ncol=2)
PFP = matrix(0, nrow=Nsim, ncol=2)
e = c(0.1,0.2)
set.seed(42)
for (i in 1:Nsim)
{
  for (j in 1:2)
  {
    weights <- c(1-e[j], e[j])
    type <- sample(1:2, size=n, replace=TRUE, prob=weights)
    ni <- table(type)
    x1 <- rmvnorm(ni[1], mu1, sigma)
    x2 <- rmvnorm(ni[2], mu2, sigma)
    F <- rbind(x1, x2)
    F <- as.data.frame(F)
    F$Class <- c(rep(0, ni[1]), rep(1, ni[2]))
    colVec <- c(rep("#FF0000", ni[1]), rep("#0000FF", ni[2]))
    
    if (e[j]==0.1)
      lda <- lda(x=F[,c(7,10,12)], grouping=F$Class, CV=TRUE) 
    else 
      lda <- lda(x=F[,c(10,11,12)], grouping=F$Class, CV=TRUE) 
    
    postProb <- lda$posterior[,2]
    ( ConfMat_LDA <- table(F$Class, postProb >= 0.5) )
    PTP[i,j] = 100*ConfMat_LDA[2,2]/sum(ConfMat_LDA[2,])
    PFP[i,j] = 100*ConfMat_LDA[1,2]/sum(ConfMat_LDA[1,])
  }
}
(avg_PFP <- colMeans(PTP))
(avg_PFP <- colMeans(PFP))
sd_PTP <- sapply(as.data.frame(PTP),sd)
sd_PFP <- sapply(as.data.frame(PFP),sd)

par(mfrow=c(1,2))
for (i in 1:2)
{
  data <- data.frame(x1 = PTP[,i], x2 = PFP[,i])
  colnames(data) <- c("PTP", "PFP")
  boxplot(data, main = paste("boxplot for epsilon =", e[i]))
}
par(mfrow=c(1,1))