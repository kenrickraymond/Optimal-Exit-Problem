# Set options and load required libraries
options(repr.plot.width = 20, repr.plot.height = 10, warn = 1)
install.packages("pacman")
pacman::p_load(readxl, fastGHQuad, stats4, pracma, rootSolve)

# Load the dataset
df <- read_excel("AAPLNVDA5minNormalized.xlsx")
spread <- df[, 2][[1]]

# Plot the spread
plot(spread, type = "l", main = "Spread Plot")
vline <- as.numeric(seq(0, length(spread), by = 78))
abline(v = vline, col = "gray", lty = 2)

# Jump Detection & Estimation
m <- 78
delta <- 1 / m
noObs <- ceiling(seq_along(spread) / m)
spreadInterval <- split(spread, noObs)

rtInterval <- array()
noInterval <- length(spreadInterval)

for (i in 1:noInterval) {
  rtInterval[i] <- list(diff(spreadInterval[[i]]))
}

# Realized Variance (RV)
realizedVarInterval <- array()
for (i in 1:noInterval) {
  productTemp <- array()
  for (k in 1:length(rtInterval[[i]])) {
    productTemp[k] <- (rtInterval[[i]][k])^2
  }
  realizedVarInterval[i] <- sum(na.omit(productTemp))
}
realizedVarInterval <- na.omit(realizedVarInterval)

# Realized Bipower Variation (BV)
sumInterval <- array()
bipowerVarInterval <- array()
for (i in 1:noInterval) {
  productTemp <- array()
  for (k in 2:length(rtInterval[[i]])) {
    productTemp[k] <- abs(rtInterval[[i]][k]) * abs(rtInterval[[i]][k - 1])
  }
  sumInterval[i] <- sum(na.omit(productTemp))
  bipowerVarInterval[i] <- (pi / 2) * m / (m - 1) * sumInterval[[i]]
}
bipowerVarInterval <- na.omit(bipowerVarInterval)

# Jump Magnitude
RJV <- realizedVarInterval - bipowerVarInterval
RJ <- RJV / realizedVarInterval

# Tripower Quarticity (TP)
tripowerInterval <- array()
mu <- (2^(2 / 3)) * gamma(7 / 6) / gamma(1 / 2)
for (i in 1:length(rtInterval)) {
  productTemp <- array()
  for (k in 3:length(rtInterval[[i]])) {
    productTemp[k] <- abs(rtInterval[[i]][k - 2])^(4 / 3) * abs(rtInterval[[i]][k - 1])^(4 / 3) * abs(rtInterval[[i]][k])^(4 / 3)
  }
  sumInterval[i] <- sum(na.omit(productTemp))
  tripowerInterval[i] <- m * (mu)^(-3) * m / (m - 2) * sumInterval[[i]]
}

# Test Statistic ZJ
ZJ <- array()
for (i in 1:length(RJ)) {
  ZJ[i] <- RJ[i] / (((pi / 2)^2 + pi - 5) * delta * max(1, tripowerInterval[i] / bipowerVarInterval[i]))^(1 / 2)
}

# Jump Magnitude Jt
retCurrent <- array()
retPrevious <- array()
ret <- array()
intervalOne <- sign(spreadInterval[[1]][length(spreadInterval[[1]])] - spreadInterval[[1]][1])
ocSpread <- array(intervalOne)
for (i in 2:noInterval) {
  retCurrent[i] <- spreadInterval[[i]][length(spreadInterval[[i]])] - spreadInterval[[i]][1]
  retPrevious[i] <- spreadInterval[[i - 1]][length(spreadInterval[[i - 1]])] - spreadInterval[[i - 1]][1]
  ret[i] <- retCurrent[i] - retPrevious[i]
  rt <- append(ocSpread, ret)
}
sign <- sign(na.omit(rt))

alpha <- 0.95
x <- rnorm(100000, 0, 1)
pdf <- density(x)
f <- approxfun(pdf$x, pdf$y, yleft = 0, yright = 0)
cdf <- function(x) {
  integrate(f, -Inf, x)$value
}
invcdf <- function(q) {
  uniroot(function(x) {
    cdf(x) - q
  }, range(x))$root
}
threshold <- as.numeric(invcdf(alpha))
for (i in 1:length(ZJ)) {
  if (ZJ[i] > threshold) {
    ZJ[i] <- 1
  } else {
    ZJ[i] <- 0
  }
}
JT <- sign * (RJV * ZJ)^(1 / 2)

# Remove Jumps
jumpInterval <- which(JT != 0, arr.ind = T)
for (i in 1:length(spreadInterval)) {
  if (i %in% jumpInterval) {
    spreadInterval[[i]] <- replace(spreadInterval[[i]], 1, spreadInterval[[i]][1] - JT[i])
  }
}
lambda <- length(jumpInterval) / noInterval
k <- mean(JT)
removedJumpProcess <- unlist(spreadInterval)
plot(removedJumpProcess - spread, type = "l", main = "Removed Jump Process")

# Parameter Estimation
checkOU <- function(theta) {
  if (theta[1] <= 0) message("\nthe process is not stationary\n")
  if (theta[2] <= 0) stop("variance must be positive")
}
mydsOU <- function(x, theta, log = FALSE) {
  checkOU(theta)
  dnorm(x, mean = 0, sd = theta[2] / sqrt(2 * theta[1]), log = log)
}
ou.lik <- function(x) {
  function(theta1, theta2) {
    n <- length(x)
    dt <- deltat(x)
    -sum(mydsOU(x = x[2:n], theta = c(theta1, theta2), log = TRUE))
  }
}
ou.fit <- mle(ou.lik(spread), start = list(theta1 = 0.2, theta2 = 0.7), method = "L-BFGS-B")
ou.coe <- coef(ou.fit)
mu <- ou.coe[[1]]
sigma <- ou.coe[[2]]

# Numerical Solution
xmin <- min(spread)
xmax <- max(spread)
Nx <- 1170
dx <- (xmax - xmin) / Nx
x <- spread
f <- -mu * x
alpha <- -2
alphaSpread <- spread
alphaSpread[alphaSpread < alpha] <- 0
alphaIndex <- which(alphaSpread == 0)[1]
alphaIndex[is.na(alphaIndex)] <- 0
if (alphaIndex != 0) {
  plot(spread, type = "l", main = "Alpha Hitting Time")
  abline(v = alphaIndex)
  spread <- spread[1:alphaIndex]
  alpha <- spread[alphaIndex][1]
}

sigma2 <- sigma^2
dx2 <- dx^2
A <- matrix(0, nrow = length(spread), ncol = length(spread))
for (i in 1:(length(spread) - 1)) {
  a <- -mu * x[i] / dx - 1 / 2 * sigma2 / dx2
  c <- -1 / 2 * sigma2 / dx2
  A[i, i + 1] <- c
  A[i + 1, i] <- a
}
b <- mu * x / dx + sigma2 / dx2
diag(A) <- b[1:length(spread)]

rule <- gaussHermiteData(3000)
jumpIndex <- which(removedJumpProcess[1:length(spread)] - spread != 0)
integralSolve <- function(previousV) {
  sequenceX <- seq(previousV)
  sequenceXPlusY <- seq(previousV)
  for (i in 1:length(jumpIndex)) {
    sequenceXPlusY[jumpIndex[i]] <- sequenceXPlusY[jumpIndex[i]] + JT[i]
  }
  splinePoly <- splinefun(sequenceX, previousV, method = "fmm")
  splinePolyPlusy <- splinefun(sequenceXPlusY, previousV, method = "fmm")
  vOfX <- ghQuad(splinePoly, rule)
  vOfXPlusY <- ghQuad(splinePolyPlusy, rule)
  integral <- vOfXPlusY - vOfX
  return(integral)
}

epsilon <- 0.1
threshold <- 1e-5
initialIter <- 0.5
vZero <- array(initialIter, length(spread))
integral <- integralSolve(vZero)
fk <- f + lambda * integral
newV <- solve(A, fk)
while (epsilon > threshold) {
  oldV <- newV
  integral <- integralSolve(newV)
  fk <- fk + lambda * integral
  newV <- solve(A, fk)
  epsilon <- max(newV - oldV)
}
v <- newV

# Plot v
plot(v, type = "l", col = "blue", lwd = 2, main = "Plot of v", xlab = "Index", ylab = "v(x)")
