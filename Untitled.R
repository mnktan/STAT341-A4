happy <- read.csv("world_happiness.csv", header = TRUE) 
srsSampIndex <- read.table("srsSampIndex.txt")$V1 
srsSamp <- happy[srsSampIndex, ]

### Problem 1 ###

## automated functions ##

popSize <- function(pop) {
  nrow(as.data.frame(pop))
}
sampSize <- function(samp) {
  popSize(samp)
}

createInclusionProbFn <- function(pop, sampSize) {
  N <- popSize(pop)
  n <- sampSize
  function(u) {
    n/N
  }
}

createJointInclusionProbFn <- function(pop, sampSize) {
  N <- popSize(pop)
  n <- sampSize
  function(u, v) {
    ## Note that the answer depends on whether u and v are the same or different
    if (u == v) {
      n/N
    } else {
      (n * (n - 1))/(N * (N - 1))
    }
  }
}

createHTestimator <- function(pi_u_fn) {
  function(samp, variateFn) {
    Reduce(`+`, Map(function(u) {
      variateFn(u)/pi_u_fn(u)
    }, samp), init = 0)
  }
}


createHTVarianceEstimator <- function(pop, pi_u_fn, pi_uv_fn) {
  function(samp, variateFn) {
    Reduce(`+`, Map(function(u) {
      pi_u <- pi_u_fn(u)
      y_u <- variateFn(u)
      Reduce(`+`, Map(function(v) {
        pi_v <- pi_u_fn(v)
        pi_uv <- pi_uv_fn(u, v)
        y_v <- variateFn(v)
        Delta_uv <- pi_uv - pi_u * pi_v
        
        result <- (Delta_uv * y_u * y_v)
        result <- result/(pi_uv * pi_u * pi_v)
        result
      }, samp), init = 0)
    }, samp), init = 0)
    
  }
}

## automated functions ##

# ai
createvariateFnN <- function(popData, variate, N = 1) {
  function(u) {
    popData[u, variate]/N
  }
}

N = 156
n = 50

inclusionProb <- createInclusionProbFn(1:N, sampSize=n)
happyHTestimator <- createHTestimator(inclusionProb)

happyAvg <- createvariateFnN(happy, "Score", N=N)
happyAvg_HT_srswor <- happyHTestimator(srsSampIndex, happyAvg) # value is 5.3863
print(happyAvg_HT_srswor)

#aii

estVarHT <- function(y_u, pi_u, pi_uv) {
  ## y_u = an n element array containing the variate values for the sample
  ## pi_u = an n element array containing the (marginal) inclusion ## probabilities for the sample
  ## pi_uv = an nxn matrix containing the joint inclusion probabilities ## for the sample
  delta <- pi_uv - outer(pi_u, pi_u)
  estimateVar <- sum((delta/pi_uv) * outer(y_u/pi_u, y_u/pi_u)) 
  return(abs(estimateVar))
}

y_u <- srsSamp$Score/N
pi_u <- rep(n/N, n)
pi_uv <- matrix((n * (n - 1)) / (N * (N - 1)),
                nrow = n, ncol = n)
diag(pi_uv) <- pi_u

se_HT_srswor <- sqrt(estVarHT(y_u, pi_u, pi_uv))
print(se_HT_srswor) # 0.1341426

## alt method

inclusionJointProb <- createJointInclusionProbFn(1:N, sampSize=n)
happyVarianceEstimator <- createHTVarianceEstimator(1:N, 
                                                    pi_u_fn=inclusionProb, 
                                                    pi_uv_fn=inclusionJointProb)

happyAvg_se_srswor <- sqrt(happyVarianceEstimator(srsSampIndex, happyAvg))
print(happyAvg_se_srswor)

## alt method

#aiii
# interval is (5.118015 5.654585)
happy_srswor <- happyAvg_HT_srswor + 2 * c(-1, 1) * happyAvg_se_srswor
print(happy_srswor)


# b

happy_sizes <- seq(5, 100, by=5)
happy_biases <- c(NA);
happy_var <- c(NA);
happy_mse <- c(NA);
happy_cov <- c(NA);

est <- rep(0, 10000)
ci <- matrix(0, nrow = 10000, ncol = 2)
happyMean <- mean(happy$Score)

for (val in happy_sizes) {
  pi_u = rep(val/N, val)
  pi_uv = matrix((val * (val - 1)) / (N * (N - 1)),
                  nrow = val, ncol = val)
  diag(pi_uv) = pi_u
  
  for (i in 1:10000) {
    samp = sample(happy$Score, size=val, replace=FALSE)
    y_u = samp/N
    est[i] = sum(y_u/pi_u)
    se = sqrt(estVarHT(y_u, pi_u, pi_uv))
    ci[i, ] = sum(y_u/pi_u) + 2 * c(-1, 1) * se
  }
  
  # bias, variance, MSE, and coverage
  bias_srswor = mean(est - happyMean)
  variance_srswor = var(est)
  MSE_srswor = mean((est - happyMean)^2)
  coverage = apply(X=ci, MARGIN = 1, FUN = function(u) {
    happyMean >= u[1] & happyMean <= u[2]
  })
  
  # add values to the columns
  happy_biases[val/5] = bias_srswor;
  happy_var[val/5] = variance_srswor;
  happy_mse[val/5] = MSE_srswor;
  happy_cov[val/5] = mean(coverage);
}

# plots of bias, variance, mse, coverage vs n
par(mfrow = c(1, 4))
plot(happy_sizes, happy_biases, main="Bias vs n", 
     xlab="n", ylab="Bias", ylim=c(-0.1, 0.1))
abline(h=0, lty=3, col="red")
plot(happy_sizes, happy_var, main="Variance vs n", 
     xlab="n", ylab="Variance", ylim=c(0, 0.25))
plot(happy_sizes, happy_biases, main="MSE vs n", 
     xlab="n", ylab="MSE", ylim=c(0, 0.25))
plot(happy_sizes, happy_cov, main="Coverage vs n", 
     xlab="n", ylab="Coverage", ylim=c(0, 1))
abline(h=0.95, lty=3, col="red")


### QUESTION 2 ###

# c)

stratLabel <- seq(1:156)
hapcont <- as.character(happy$Continent)

for (i in 1:156) {
  j = hapcont[i]
  if (j == "Africa")
    stratLabel[i] <- 1
  else if (j == "Asia")
    stratLabel[i] <- 2
  else if (j == "Europe")
    stratLabel[i] <- 3
  else if (j == "North America")
    stratLabel[i] <- 4
  else if (j == "Oceania")
    stratLabel[i] <- 5
  else
    stratLabel[i] <- 6
}

happy <- data.frame(happy, stratLabel)
table(happy$stratLabel)


### functions for II ### 
getInclusionProbStrat <- function(stratLabel, stratSampSize) {
  H <- length(stratSampSize)
  stratSize <- as.numeric(table(stratLabel))
  N <- sum(stratSize)
  pi_u <- rep(0, N)
  for (h in 1:H) {
    pi_u[which(stratLabel == h)] <- stratSampSize[h]/stratSize[h]
  }
  return(pi_u)
}

getJointInclusionProbStrat <- function(stratLabel, stratSampSize) {
  H <- length(stratSampSize)
  stratSize <- as.numeric(table(stratLabel))
  N <- sum(stratSize)
  pi_uv <- matrix(0, N, N)
  for (u in 1:N) {
    for (v in 1:N) {
      if (u == v) {
        pi_uv[u, v] <- stratSampSize[stratLabel[u]]/stratSize[stratLabel[u]]
      } else {
        if (stratLabel[u] == stratLabel[v]) {
          pi_uv[u, v] <- (stratSampSize[stratLabel[u]] * (stratSampSize[stratLabel[u]] - 1))/(stratSize[stratLabel[u]] * (stratSize[stratLabel[u]] - 1))
        } else {
          pi_uv[u, v] <- (stratSampSize[stratLabel[u]]/stratSize[stratLabel[u]]) * (stratSampSize[stratLabel[v]]/stratSize[stratLabel[v]])
        }
      }
    }
  }
  return(pi_uv)
}

stratRS <- function(stratLabel, stratSampSize){
  H <- length(stratSampSize)
  sampIndex <- list()
  for(h in 1:H){
    sampIndex[[h]] <- sample(which(stratLabel == h), size = stratSampSize[h], replace = FALSE)
  }
  return(unlist(sampIndex))
}
### functions for II ### 


# d)

# i)
stratSampIndex <- read.table("stratSampIndex.txt")$V1 
stratSamp <- happy[stratSampIndex, ]

N = 156
n = 50

stratIncProb <- getInclusionProbStrat(stratLabel, c(14, 13, 15, 4, 1, 3))  # marginal inc. probs for units in sample
strat_y_u <- stratSamp$Score/N  # variate values being summed in the sample

stratAvg_HT_srswor <- sum(strat_y_u/stratIncProb[stratSampIndex])
print(stratAvg_HT_srswor)  # 5.415015

# ii)

stratJointIncProb <- getJointInclusionProbStrat(stratLabel, c(14, 13, 15, 4, 1, 3))
diag(stratJointIncProb) <- stratIncProb
strat_se_HT_srswor <- sqrt(estVarHT(strat_y_u, stratIncProb[stratSampIndex],
                                    stratJointIncProb[stratSampIndex,][,stratSampIndex])) # 0.128552

# iii)
# interval is (5.157911 5.672119)
strat_srswor <- stratAvg_HT_srswor + 2 * c(-1, 1) * strat_se_HT_srswor
print(strat_srswor)

# e)

strat_est <- rep(0, 10000)
strat_ci <- matrix(0, nrow=10000, ncol=2)

for (i in 1:10000) {
  samp <- stratRS(stratLabel, c(14, 13, 15, 4, 1, 3))
  y_u <- happy$Score[samp]/N
  strat_est[i] <- sum(y_u/stratIncProb[samp])
  se <- sqrt(estVarHT(y_u, stratIncProb[samp],
                     stratJointIncProb[samp,][,samp]))
  strat_ci[i, ] <- sum(y_u/stratIncProb[samp]) + 2 * c(-1, 1) * se
}

par(mfrow = c(1,2))
hist(strat_est, col="lightgrey", main="HT Estimates, SRSWOR (n=50)", xlab="Average Happiness Score")
abline(v=happyMean, col="red", lwd=2)
strat_coverage <- apply(X=strat_ci, MARGIN=1, function(u) {
  happyMean >= u[1] & happyMean <= u[2]
})
plot(0, type="n", ylim=c(0,100), xlim=c(min(strat_ci[,1]),max(strat_ci[,2])),
     xlab="95% Confidence Intervals", ylab="Sample Number", 
     main=paste0("Coverage Prob. = ", round(100*mean(strat_coverage),2), "%"))
for (i in 1:10000) {
  segments(x0=strat_ci[i, 1], y0=i, x1=ci[i,2], y1=i, col=adjustcolor("gray",alpha=0.3))
}
abline(v=happyMean, col="red", lwd=2)

# f)

strat_bias_srswor = mean(strat_est - happyMean)
strat_variance_srswor = var(strat_est)
# old: strat_variance_srswor + strat_bias_srswor^2
# 0.009427777
strat_MSE_srswor = mean((strat_est - happyMean)^2)
strat_cov = mean(strat_coverage)

# results for bias, variance, MSE, coverage prob.
cat("Q2e results:\n ","bias = ", strat_bias_srswor,
             " variance = ", strat_variance_srswor,
             " MSE = ", strat_MSE_srswor,
             " coverage prob. = ", strat_cov, "\n")


# results from question 1b (n=50)
cat("Q1b results:\n", "bias = ",happy_biases[10],
             " variance = ", happy_var[10],
             " MSE = ", happy_mse[10],
             " coverage prob. = ", happy_cov[10], "\n")


### QUESTION 3 ###

# a)
pop <- list(pop1 = happy[order(happy$Health, decreasing=TRUE), ][1:78, ], 
            pop2 = happy[order(happy$Health, decreasing=TRUE), ][79:156, ])

par(mfrow = c(1, 2))
hist(pop[[1]]$Score, col=adjustcolor("firebrick", 0.7), freq = FALSE,
     xlab = "Happiness Score", main = "Healthy Countries", xlim = c(2, 8))
abline(v = mean(pop[[1]]$Score), col = "darkgreen", lwd = 2)
hist(pop[[2]]$Score, col=adjustcolor("blue", 0.7), freq = FALSE,
     xlab = "Happiness Score", main = "Unhealthy Countries", xlim = c(2, 8))
abline(v = mean(pop[[2]]$Score), col = "darkgreen", lwd = 2)

# c)

par(mfrow = c(1,1))

D <- function(pop) {
  abs(mean(pop[[1]]$Score) - mean(pop[[2]]$Score))
}
d_obs <- D(pop)
print(d_obs)

# ii)

# function from class
mixRandomly <- function(pop) { 
  pop1 <- pop$pop1
  n_pop1 <- nrow(pop1)
  
  pop2 <- pop$pop2 
  n_pop2 <- nrow(pop2)
  
  mix <- rbind(pop1, pop2)
  select4pop1 <- sample(1:(n_pop1 + n_pop2), n_pop1, replace = FALSE)
  
  new_pop1 <- mix[select4pop1, ] 
  new_pop2 <- mix[-select4pop1, ] 
  list(pop1 = new_pop1, pop2 = new_pop2)
}

diffPops <- sapply(1:5000, FUN = function(...) {
  D(mixRandomly(pop))
})

hist(diffPops, breaks = 20, main = "Randomly Mixed Populations",
     xlab = "Discrepancy D", col = adjustcolor("darkorchid", 0.4),
     xlim=c(0, 1.6))
abline(v = D(pop), col = "darkgreen", lwd = 2)

# iii)

mean(diffPops >= D(pop))

# d)

D1 <- function(pop) {
  abs(sd(pop[[1]]$Score) / sd(pop[[2]]$Score) - 1)
}
d1_obs <- D1(pop)
print(d1_obs)

# ii)

diffPops1 <- sapply(1:5000, FUN = function(...) {
  D1(mixRandomly(pop))
})

hist(diffPops1, breaks = 20, main = "Randomly Mixed Populations",
     xlab = "Discrepancy D", col = adjustcolor("darkorchid", 0.4))
abline(v = D1(pop), col = "darkgreen", lwd = 2)

# iii)

mean(diffPops1 >= D1(pop))






