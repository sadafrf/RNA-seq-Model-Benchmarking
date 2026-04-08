library(ggplot2)
library(dplyr)
library(tidyr)
library(rjags)
library(coda)


counts <- matrix(nrow=1000,ncol=100)
for (i in nrow(counts)){
  rownames(counts) <- paste0("Gene", i) 
}

#Initialize alpha dispersion parameter with gamma prior. This prior is uninformative.
#The names for the dispersion list are set to be the rownames for the count matrix
alpha <- rgamma(n=1000, shape=2,rate=0.5)
names(alpha) <- rownames(counts)

#Initalize bernoulli (1 or 0) if a gene is DE expressed with probability 0.1. 
#The names of the probability list are set to be the rownames for the count matrix
de_p <- rbinom(n=1000, size=1, prob=0.1)
names(de_p) <- rownames(counts)

#Initialize delta as a normal distribution with mean 0 and variance 2
delta <- rnorm(n=length(de_p == 1), mu = 1, sd = sqrt(2))
names(delta) <- rownames(counts[de_p == 1, ])

#Initialize gene average (mu) as a poisson prior with sampling rate = 10 
mu <- rpois(n=1000, lambda=10)
names(mu) <- rownames(counts)

colnames(counts) <- "Case" 
colnames(counts)[1:50] <- "Control"

for (g in nrow(counts)) {
  for (i in c("Case", "Control")) {
    if ((de_p[g] != 0) && (i == "Case")) {
      counts[g,i] <- rnbinom(n=50, mu= mu[g] + delta[g], size = alpha[g])
    }
    else {
      counts[g,i] <- rnbinom(n=50, mu = mu[g], size = alpha[g])
    }
  }
}

print("Data is successfully simulated")




