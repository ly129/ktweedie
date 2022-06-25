## code to prepare `dat` dat goes here
library(tweedie)

rand_tweedie<- function(mu,...) {
  Y <- rtweedie(1, mu = mu,...)
  Y
}

phi <- 0.5
rho <- 1.5
P <- 50
N <- 200
beta.true <- c(6, -4, 3, 2, -2, rep(0, P-5))

x <- matrix(rnorm(N * P), N, P)
Fx <- sin(x) %*% beta.true
mu <- exp(Fx)
y = sapply(mu, rand_tweedie, xi = rho, phi = phi)

dat <- list(x, y)
names(dat) <- c("x", "y")

usethis::use_data(dat, overwrite = TRUE)
