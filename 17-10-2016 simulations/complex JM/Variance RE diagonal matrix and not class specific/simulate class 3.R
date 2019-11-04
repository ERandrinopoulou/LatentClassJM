
if (ClassSim == 3) n = 900


age <- rnorm(n, 45, 15.69578) #age <- mat1$age


##############################
#n <- 100 # number of subjects
##############################
K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model 1
betas <- c("(Intercept)" = 0.0317, "Group1" = -1.9568, "Time1" =  -0.00578)
sigma.y <- 0.6920  # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" = 2.8526, "Age" = -0.1150) # coefficients for baseline covariates
alpha <-    0.5833#-0.5 # association parameter - value
Dalpha <- 0#2.9 # association parameter - slope

phi <- 1.8 #1.6458 # shape for the Weibull baseline hazard
mean.Cens <- 10 #12 # mean of the exponential distribution for the censoring mechanism

D <- diag(c(0.5337, 0.5560)^2)

################################################

Bkn <- c(0, 19.5)
kn <- c(2.1, 5.5)

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- sample(c(0,1), n, replace = TRUE) #rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
####################################
#age <- rnorm(n, 46.39673, 2.69578)
####################################
DF <- data.frame(year = times, drug = factor(rep(group, each = K)), age = rep(age, each = K))


#X <- model.matrix(~ 0 + drug  + ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
#Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
X <- model.matrix(~ drug + year, data = DF)
Z <- model.matrix(~ year, data = DF)


# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Age" = age)

################################################

#simulate random effects
library(MASS)

b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)


# simulate event times
eta.t <- as.vector(as.matrix(W) %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group0 <- 1 - group[i]
    group1 <- group[i]
    ages <- age[i]
    #NS <- ns(s, knots = kn, Boundary.knots = Bkn)
    #DNS <- dns(s, knots = kn, Boundary.knots = Bkn)
    #XX <- cbind(group0, group1, NS[, 1], NS[, 2], NS[, 3])
    #ZZ <- cbind(1, NS)
    
    #XXd <- cbind(DNS[, 1], DNS[, 2], DNS[, 3])
    #ZZd <- DNS
    
    #NS <- ns(s, 1)/((Sp/N)[2])
    NS <- s
    XX <- cbind(1, group1, (NS))
    ZZ <- cbind(1, NS)
    
    fval <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
    #fsl <- as.vector(XXd %*% betas[3:5] + rowSums(ZZd * b[rep(i, nrow(ZZd)), 2:4]))
    
    
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + fval * alpha) # + fsl * Dalpha)
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
  Up <- 5000
  tries <- 5
  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up + 5000
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- runif(n, 0, 2 * mean.Cens)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id <- data.frame(IDnr = unique(id), years = Time, status = event, age = W[, 2])
names(dat) <- c("echotime", "group", "age", "IDnr", "y0", "years", "status")

#############################
#ind2 <- unique(id)[Time > 19.5]
#'%!in%' <- function(x,y)!('%in%'(x,y))
#dat$years[dat$IDnr %in% ind2] <- 19.5
#dat$status[dat$IDnr %in% ind2] <- 0
#dat.id <- dat[!duplicated(dat$IDnr), ]
#############################

data <- dat
data.id <- dat.id

