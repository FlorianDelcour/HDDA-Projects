setwd("C:/Users/flode/OneDrive - Universite de Liege/Bureau/HDDA/Project 1")
wine <- read.table(file = "wine.txt", header = TRUE, sep = ";")
attach(wine)


#Q1

# Number of miss
length(wine[wine == 'NA'])

# We don't have any missing data so no need to handle those

#Q2

library(psych)
multi.hist(wine, global = FALSE) 
p_before <- vector()
p_after <- vector()

for(i in 1:13)
{
  p_before[i] <- shapiro.test(wine[,i])$p.value
  if(shapiro.test(wine[,i])$p.value < 0.05)
    wine[,i] <- log(1 + wine[,i])
  p_after[i] <- shapiro.test(wine[,i])$p.value
}
multi.hist(wine, global = FALSE) 


#Q3
# First we compute the classical mahalanobis distance
d_mahalanobis <- mahalanobis(wine, colMeans(wine), cov(wine))


# Then we compute it but with the covariance matrix estimated by the mcd estimator
library(MASS)
set.seed(5)
n <- dim(wine)[1]
h <- floor(n*0.75)
rob_est <- cov.rob(wine, quantile.used=h, method = "mcd",cor = TRUE)
rob_dist <- mahalanobis(wine, rob_est$center, rob_est$cov)

# We do a DD-plot
plot(d_mahalanobis, rob_dist, main="DD-plot", xlab="Mahalanobis Distance", ylab="Robust Distance")
abline(h=qchisq(0.95,13), lty=2, col=2)
abline(v=qchisq(0.95,13), lty=2, col=2)

# We compute the number of outliers (with 95% cutoff)

rob_outliers <- sum(rob_dist > qchisq(0.95,13))
cla_outliers <- sum(d_mahalanobis > qchisq(0.95,13))


#Q4
cla_cor <- cor(wine)
rob_cor <- rob_est$cor
library(corrplot)
corrplot(cla_cor)
corrplot(rob_cor)
max_correlations <- vector()
max_correlations[1] <- 1
for(i in 2:4)
{
  max_correlations[i] <- max(abs(rob_cor[rob_cor < abs(max_correlations[i-1])-0.001]))
}

plot(wine[,6], wine[,7], xlab = colnames(wine)[6], 
     ylab = colnames(wine)[7], main = "Most correlated couple")

plot(wine[,7], wine[,12], xlab = colnames(wine)[7], 
     ylab = colnames(wine)[12], main = "Second most correlated couple")

plot(wine[,7], wine[,9], xlab = colnames(wine)[7], 
     ylab = colnames(wine)[9], main = "Third most correlated couple")

plot(wine[,8], wine[,12], xlab = colnames(wine)[8], 
     ylab = colnames(wine)[12], main = "Biggest negative correlation")


#Q5
source("plot_graph.R")
library(glasso)
wine2 <- scale(wine)
rob_est2 <- cov.rob(wine2, quantile.used=h, method = "mcd",cor = TRUE)
n <- dim(wine2)[1]
res <- glasso(s=rob_est2$cov, rho=0.6, nobs=h)
theta <- res$wi
plot_graph(invCov=theta,  varnames=names(wine))
#make plot of the evolution of number of edges with regard to rho used for glasso
nb_edges <- vector()
x <- vector()
for(i in 0:10)
{
  invCov <- glasso(s=rob_est2$cov, i/10, nobs=h)$wi
  nb_edges[i+1] <- length(invCov[abs(invCov)>=1e-16])
  x[i+1] = i/10
}
plot(x, (nb_edges-13)/2, xlab = "Penalization constant", ylab = "Number of edges")






