# Gibbs Sampling algorithm for a triangle without and with bonding 

library(dplyr) 
library(ggplot2) 

#%% Gibbs Sampling algoritm for the Ising model without U ####### 

gsi <- function (n) {
  path <- matrix(rep(-1, (n * 3 + 1) * 3), ncol = 3) 
  i = 1 
  for (iteration in 1:n) { 
    path[i + 1, 1:3] = path[i, 1:3] 
    if (path[i, 2] == path[i, 3]) {
      if (runif(1) < 16 / 17) path[i + 1, 1] = path[i, 2] 
      else path[i + 1, 1] = -path[i, 2] 
    }
    else {
      if (runif(1) < 0.5) path[i + 1, 1] = -1 
      else path[i + 1, 1] = 1 
    }
    i <- i + 1 
    path[i + 1, 1:3] = path[i, 1:3] 
    if (path[i, 3] == path[i, 1]) {
      if (runif(1) < 16 / 17) path[i + 1, 2] = path[i, 3] 
      else path[i + 1, 2] = -path[i, 3] 
    }
    else {
      if (runif(1) < 0.5) path[i + 1, 2] = -1 
      else path[i + 1, 2] = 1 
    }
    i <- i + 1 
    path[i + 1, 1:3] = path[i, 1:3] 
    if (path[i, 1] == path[i, 2]) {
      if (runif(1) < 16 / 17) path[i + 1, 3] = path[i, 1] 
      else path[i + 1, 3] = -path[i, 1] 
    }
    else {
      if (runif(1) < 0.5) path[i + 1, 3] = -1 
      else path[i + 1, 3] = 1 
    }
    i <- i + 1 
  } 
  path 
}

#%% Gibbs Sampling algoritm for the Ising model with U ####### 
# We have separate rows only for X-iterations 
# Column 4 connects 1 and 2, Column 5 connects 2 and 3, Column 6 connects 3 and 1 

gsib <- function (n) {
  p <- 0.75 
  path <- matrix(rep(c( -1, -1, -1, 0, 0, 0 ), n * 3 + 1), ncol = 6) 
  i = 1 
  for (iteration in 1:n) {
    path[i + 1, ] = path[i, ] 
    if (!path[i, 4] && !path[i, 6]) { 
      if (runif(1) < 0.5) path[i + 1, 1] = -path[i + 1, 1] 
    } 
    i <- i + 1 
    path[i + 1, ] = path[i, ] 
    if (!path[i, 5] && !path[i, 4]) { 
      if (runif(1) < 0.5) path[i + 1, 2] = -path[i + 1, 2] 
    } 
    i <- i + 1 
    path[i + 1, ] = path[i, ] 
    if (!path[i, 6] && !path[i, 5]) { 
      if (runif(1) < 0.5) path[i + 1, 3] = -path[i + 1, 3] 
    } 
    i <- i + 1 
    if (path[i, 4] || path[i, 1] == path[i, 2]) {
      path[i, 4] <- rbinom(1, 1, p) 
    }
    if (path[i, 5] || path[i, 2] == path[i, 3]) {
      path[i, 5] <- rbinom(1, 1, p) 
    }
    if (path[i, 6] || path[i, 3] == path[i, 1]) {
      path[i, 6] <- rbinom(1, 1, p) 
    }
  }
  path 
}

#%% tests 

time <- Sys.time() 
stdev <- matrix(rep(0, 6), ncol = 1) 
for (i in 1:100) {
  path <- data.frame(gsi(1000)) 
  colnames(path) <- c("x1", "x2", "x3") 
  s <- path %>% select(x1, x2, x3) %>% group_by(x1, x2, x3) %>% summarise(mean = n() / nrow(path)) 
  stdev <- stdev + ( s["mean"] - c( 8/19, 1/38, 1/38, 1/38, 1/38, 1/38, 1/38, 8/19 ) ) ** 2 
}
colnames(stdev) <- c( "stdev" ) 
bind_cols(s[, 1:3], sqrt( stdev / 99 )) 
print(Sys.time() - time) 

time <- Sys.time() 
stdevb <- matrix(rep(0, 6), ncol = 1) 
for (i in 1:100) {
  pathb <- data.frame(gsib(1000)) 
  colnames(pathb) <- c("x1", "x2", "x3", "u12", "u23", "u31") 
  s <- pathb %>% select(x1, x2, x3) %>% group_by(x1, x2, x3) %>% summarise(mean = n() / nrow(pathb)) 
  stdevb <- stdevb + ( s["mean"] - c( 8/19, 1/38, 1/38, 1/38, 1/38, 1/38, 1/38, 8/19 ) ) ** 2 
}
colnames(stdevb) <- c( "stdev" ) 
bind_cols(s[, 1:3], sqrt( stdevb / 99 )) 
print(Sys.time() - time) 
