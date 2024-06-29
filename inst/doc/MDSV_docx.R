## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

## ---- echo = TRUE, eval = TRUE------------------------------------------------
library(MDSV)

## ---- echo = TRUE, eval = TRUE------------------------------------------------
args(MDSVfit)

## ---- echo = TRUE, eval = TRUE------------------------------------------------
start.date <- as.Date("2000-01-03") 
end.date   <- as.Date("2019-08-31")
data(sp500)
sp500     <- sp500[rownames(sp500)>=start.date & rownames(sp500)<=end.date,]
N         <- 2
K         <- 3
LEVIER    <- TRUE

# Model estimation : univariate log-returns (ModelType = 0)
out       <- MDSVfit(K = K, N = N, data = sp500, ModelType = 0, LEVIER = LEVIER)
# Summary
summary(out)
# Plot
plot(out,c("nic"))

# Model estimation : univariate realized variances (ModelType = 1) without leverage
out       <- MDSVfit(K = K, N = N, data = sp500, ModelType = 1, LEVIER = FALSE)
# Summary
summary(out)
# Plot
plot(out,c("dis"))

# Model estimation : joint model (ModelType = 2)
out       <- MDSVfit(K = K, N = N, data = sp500, ModelType = 2, LEVIER = LEVIER)
# Summary
summary(out)
# Plot
plot(out,c("dis","nic"))


## ---- echo = TRUE, eval = TRUE------------------------------------------------
N         <- 3
K         <- 3
LEVIER    <- TRUE
para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72,
               l = 0.78, theta = 0.876)

# Model filtering : univariate log-returns (ModelType = 0)
out       <- MDSVfilter(K = K, N = N,data = sp500, para = para, ModelType = 0, LEVIER = LEVIER, calculate.VaR = TRUE, VaR.alpha = c(0.05))
# Summary
summary(out)
# Plot
plot(out)

para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72, shape = 2.10)

# Model filtering : univariate realized variances (ModelType = 1) without leverage
out       <- MDSVfilter(K = K, N = N, data = sp500, para = para, ModelType = 1, LEVIER = FALSE, calculate.VaR = TRUE, VaR.alpha = c(0.05))
# Summary
summary(out)
# Plot
plot(out)

para      <- c(omega = 0.52, a = 0.99, b = 2.77, sigma = 1.95, v0 = 0.72, 
              xi = -0.5, varphi = 0.93, delta1 = 0.93, delta2 = 0.04, shape = 2.10,
              l = 0.78, theta = 0.876)

# Model filtering : joint model (ModelType = 2)
out       <- MDSVfilter(K = K, N = N,data = sp500, para = para, ModelType = 2, LEVIER = LEVIER, calculate.VaR = TRUE, VaR.alpha = c(0.05))
# Summary
summary(out)
# Plot
plot(out)

## ---- echo = TRUE, eval = TRUE------------------------------------------------
N         <- 3
K         <- 3
LEVIER    <- TRUE

# Model forecasting : univariate log-returns (ModelType = 0) without leverage
out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = 0, LEVIER = FALSE)
out       <- MDSVboot(fit = out_fit, n.ahead = 100, n.bootpred = 10000, rseed = 349)
# Summary
summary(out)

# Model forecasting : univariate realized variances (ModelType = 1) without leverage
out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = 1, LEVIER = FALSE)
out       <- MDSVboot(fit = out_fit, n.ahead = 100, n.bootpred = 10000, rseed = 349)
# Summary
summary(out)

# Model bootstrap forecasting : joint model (ModelType = 2) with leverage
out_fit   <- MDSVfit(K = K, N = N, data = sp500, ModelType = 2, LEVIER = LEVIER)
out       <- MDSVboot(fit = out_fit, n.ahead = 100, n.bootpred = 10000, rseed = 349)
# Summary
summary(out)


## ---- echo = TRUE, eval = TRUE------------------------------------------------
args(MDSVsim)

## ---- echo = TRUE, eval = TRUE------------------------------------------------
N                <- 2
K                <- 3
ModelType        <- 2
LEVIER           <- FALSE
n.ahead          <- 100
forecast.length  <- 756
refit.every      <- 63
refit.window     <- "recursive"
calculate.VaR    <- TRUE
VaR.alpha        <- c(0.01, 0.05, 0.1)
cluster          <- parallel::makeCluster(7)
rseed            <- 125

# rolling forecasts
out<-MDSVroll(N=N, K=K, data=sp500, ModelType=ModelType, LEVIER=LEVIER, n.ahead = n.ahead,
              forecast.length = forecast.length, refit.every = refit.every, 
              refit.window = refit.window, window.size=NULL,calculate.VaR = calculate.VaR,
              VaR.alpha = VaR.alpha, cluster = cluster, rseed = rseed)

parallel::stopCluster(cluster)
# Summary
summary(out, VaR.test=TRUE, Loss.horizon = c(1,5,10,25,50,75,100), Loss.window = 756)
# plot
plot(out, plot.type=c("VaR","sigma","dens"))


