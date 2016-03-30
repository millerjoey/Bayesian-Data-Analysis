# Homework 5, BDA
setwd("/home/joey/Desktop/Bayesian Data Analysis/")
library(mvtnorm)
library(ggplot2)

alpha <- seq(from=-5, to=10, length.out = 500)
beta <- seq(from=-10, to=40, length.out = 500)
meanV <- c(0,10)
covM <- matrix(data = c(4, 10, 10, 100), ncol = 2)
y <- c(0,1,3,5)
x <- c(-0.86, -0.3, -0.05, 0.73)
n <- c(5,5,5,5)

lgpost <- function(alpha, beta) {
  z <- -(1/2)*alpha^2+alpha*(beta-10)/10-(1/50)*(beta-10)^2 + alpha*(sum(y)) + beta*(sum(x*y)) - 5*log(1+exp(alpha+beta*x[1])) - 5*log(1+exp(alpha+beta*x[2])) - 5*log(1+exp(alpha+beta*x[3])) - 5*log(1+exp(alpha+beta*x[4]))
  return(z)
}

# Compute vectors of all combinations of alpha, beta in grid
alphaVec <- rep(alpha, times=500)
betaVec <- rep(beta, each=500)

df <- data.frame("alpha"=alphaVec, "beta"=betaVec)

lz <- lgpost(alphaVec, betaVec)
z <- exp(lz)
df <- cbind.data.frame(df, z)
names(df)[3] <- "prob"

df$prob <- df$prob/sum(df$prob)

View(df)

ggplot(data=df, aes(x=alpha, y=beta, z=prob)) +
  geom_contour(bins=10) +
  xlim(c(-4,10)) + 
  ylim(c(-10,40)) +
  labs(title="Posterior distribution of alpha and beta")

dev.print(pdf, "PostAlphaBeta")

# Sample from posterior, using conditional sampling.
# First, create matrix of densities:
pMat <- t(matrix(df$prob, nrow=500))


marginalBeta <- apply(pMat, MARGIN = 1, sum)
a <- vector(length=1000)
b <- vector(length=1000)
for (i in 1:1000) {
  b[i] <- sample(size = 1, x = beta, prob = marginalBeta)
  a[i] <- sample(size = 1, x = alpha, prob = pMat[which(beta==b[i]),])
}
# Adding an runif so the distribution is continuous.
a <- a + runif(n=1000, min = -0.01503006, max=0.01503006)
b <- b + runif(n=1000, min= -0.0501002, max=0.0501002)

dfSamples <- data.frame(a,b)

ggplot(dfSamples, aes(x=a, y=b)) +
  geom_point(size=0.1) +
  labs(x="alpha", y="beta", title="One Thousand Samples from Posterior") +
  xlim(c(-5,10)) + ylim(c(-10,40))

dev.print(pdf, "PostSamples")

# Now, plotting -alpha/beta for the LD50:
ggplot() + 
  geom_histogram(aes(x=-(a/b)), color="black", bins = 40) +
  labs(x="LD 50", title="Samples from Posterior of LD50")
  
dev.print(pdf, "LD50Samples")

# Contour plot for Prior:

lgprior <- function(alpha, beta) {
  z <- -(1/2)*alpha^2+alpha*(beta-10)/10-(1/50)*(beta-10)^2
  return(z)
}

alphaVec <- rep(alpha, times=500)
betaVec <- rep(beta, each=500)

df <- data.frame("alpha"=alphaVec, "beta"=betaVec)

lz <- lgprior(alphaVec, betaVec)
z <- exp(lz)
df <- cbind.data.frame(df, z)
names(df)[3] <- "prob"

df$prob <- df$prob/sum(df$prob)

View(df)

ggplot(data=df, aes(x=alpha, y=beta, z=prob)) +
  geom_contour(bins=10) +
  xlim(c(-4,10)) + 
  ylim(c(-10,40)) +
  labs(title="Prior distribution of alpha and beta")

dev.print(pdf, "PriorAlphaBeta")

# Contour plot for likelihood:

lgliklihood <- function(alpha, beta) {
  z <- alpha*(sum(y)) + beta*(sum(x*y)) - 5*log(1+exp(alpha+beta*x[1])) - 5*log(1+exp(alpha+beta*x[2])) - 5*log(1+exp(alpha+beta*x[3])) - 5*log(1+exp(alpha+beta*x[4]))
  return(z)
}

alphaVec <- rep(alpha, times=500)
betaVec <- rep(beta, each=500)

df <- data.frame("alpha"=alphaVec, "beta"=betaVec)

lz <- lgliklihood(alphaVec, betaVec)
z <- exp(lz)
df <- cbind.data.frame(df, z)
names(df)[3] <- "prob"

df$prob <- df$prob/sum(df$prob)

View(df)

ggplot(data=df, aes(x=alpha, y=beta, z=prob)) +
  geom_contour(bins=10) +
  xlim(c(-4,10)) + 
  ylim(c(-10,40)) +
  labs(title="Likelihood for alpha and beta")

dev.print(pdf, "LikelihoodAlphaBeta")

# Approximation to the posterior in problem 4.1
theta <- seq(from = -3.313, to=4.687, length.out = 1000)
dens <- dnorm(theta, mean=0.687, sd = sqrt(.706))
df <- cbind.data.frame(theta, dens)

ggplot(data=df, aes(x=theta, y=dens)) +
  geom_line() +
  labs(title="N(0.687, .707) distribution", y="density")

dev.print(pdf, "DensNormalApprox.pdf")
