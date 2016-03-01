library("ggplot2")

# From BDA edition 3.

# Question 1a
y <- seq(-10, 10, length.out=1000)
dens <- dnorm(y, mean = 1, sd = 2) + dnorm(y, mean = 2, sd = 2)
qplot(x = y, y = dens, geom = "line", main = "Density", xlim = c(-5, 8))

# Question 8c
# n=10
mu <- (180/40^2 + 150*10/20^2)/(1/40^2+10/20^2)
tau <- (1/40^2 + 10/20^2)^-1
upper <- qnorm(0.975, mean = mu, sd = sqrt(tau))
lower <- qnorm(0.025, mean = mu, sd = sqrt(tau))


mu <- (180/40^2 + 150*10/20^2)/(1/40^2+10/20^2)
tau <- (1/40^2 + 10/20^2)^-1 + 20^2
upper <- qnorm(0.975, mean = mu, sd = sqrt(tau))
lower <- qnorm(0.025, mean = mu, sd = sqrt(tau))

# Question 8d
# n=100
mu <- (180/40^2 + 150*100/20^2)/(1/40^2+100/20^2)
tau <- (1/40^2 + 100/20^2)^-1
upper <- qnorm(0.975, mean = mu, sd = sqrt(tau))
lower <- qnorm(0.025, mean = mu, sd = sqrt(tau))


mu <- (180/40^2 + 150*100/20^2)/(1/40^2+100/20^2)
tau <- (1/40^2 + 100/20^2)^-1 + 20^2
upper <- qnorm(0.975, mean = mu, sd = sqrt(tau))
lower <- qnorm(0.025, mean = mu, sd = sqrt(tau))

# Question 9a
x <- seq(0, 1, length.out=1000)
dens <- dbeta(x, shape1 = 1, shape2 = 2/3)
qplot(x = x, y = dens, geom = "line", main = "Density", ylim=c(0,4))

# Question 9b
x <- seq(0, 1, length.out=1000)
dens <- dbeta(x, shape1 = 651, shape2 = (350+2/3))
qplot(x = x, y = dens, geom = "line", main = "Density", ylim=c(0,30))

# Question 11a
theta <- seq(0, 100, length.out = 10000)
y <- c(43, 44, 45, 46.5, 47.5)
likelihood <- function(y, theta) {
  lprod <- 0
  for(i in 1:length(y)) {
    lprod <- log(1/(1+(y[i]-theta)^2)) + lprod
  }
  return(exp(lprod))
}

post <- likelihood(y, theta)
post <- post/sum(post)

qplot(x = theta, y = post, geom = "line", main = "PMF", xlim=c(40, 50))

# Question 11b (assumes workspace of Q11a)
postSamples <- sample(x = theta, prob = post, size = 1000)
qplot(postSamples, geom = "histogram", main = "Sampled Values from Posterior", bins=30, xlim=c(35, 55))

# Question 11c (assumes workspace of Q11a and Q11b)
y <- seq(0, 110, length.out = 1000)
yNew <- rep(NA, times=1000)
for (j in 1:length(postSamples)) {
  prob <- lapply(X = y, FUN = likelihood, theta=postSamples[j])
  yNew[j] <- sample(x=y, prob = as.numeric(prob), size = 1)
}

qplot(yNew, geom = "histogram", main = "Posterior Predictive Distribution", bins=40)
