# Longitudinal submodel
fm1 <- lme(y0 ~ group + echotime, data = data,
           na.action = na.omit,
           random =~ echotime | IDnr)
lmeObject <- fm1

timeVar <- "echotime"
lag <- 0
survMod <- "spline-PH"

id <- data$IDnr 
offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

Time <- data.id$years

# Survival submodel
W <- model.matrix(~ -1 + age, data.id)

# Design matrices
formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

#
data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)


mfX.id <- model.frame(TermsX, data = data.id)  
mfZ.id <- model.frame(TermsZ, data = data.id)  
Xtime <- model.matrix(formYx, mfX.id)
Ztime <- model.matrix(formYz, mfZ.id)

###
####################################################
eventD <- data.id$status
nT <- length(Time)
zeros <- numeric(nT)

y.long <- model.response(mfX, "numeric")
y <- list(y = y.long, offset = offset, logT = log(Time),
          eventD = eventD, zeros = zeros, lag = lag)

# Time and data for approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)

P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))

# Design matrices for approximation
mfX <- model.frame(TermsX, data = data.id2)   
mfZ <- model.frame(TermsZ, data = data.id2)   
Xs <- model.matrix(formYx, mfX)
Zs <- model.matrix(formYz, mfZ)

####################


# MCMC details
con <- list(program = "JAGS", n.chains = 1, n.iter = 50000,
            n.burnin = 25000, n.thin = 10, n.adapt = 500, K = 1000,
            C = 5000, working.directory = getwd(), bugs.directory = "C:/Program Files/WinBUGS14/",
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 3,
            bugs.seed = 1, quiet = FALSE)


x <- list(X = X, Z = Z,  W = if (survMod == "weibull-PH") {
  if (is.null(W)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                            W)
} else {
  if (is.null(W)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(W) == 1) cbind(W, rep(0, nT)) else W
  }
}
) 


# Baseline hazard
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- if (con$ObsTimes.knots) {
    Time
  } else {  Time[event == 1]    }
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr

WBH <- splineDesign(rr, Time, ord = con$ordSpline)
if (any(colSums(WBH) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")

# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation
WBHs <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(WBH = WBH, WBHs = WBHs))


#####################################################

#################################
ncX <- ncol(X)
ncZ <- ncol(Z)
ncW <- ncol(x$W)
ncWBH <- ncol(x$WBH)
ncWBHs <- ncol(x$WBHs)


C <- con$C

nb <- ncZ 

T <- cbind(1, data.id$age)

# priors/hyperpriors
beta <- "string"
for (jj in 1:Class) {
  beta[jj] <- paste0("betas", jj, " <- rep(0, ncX)", sep = "", collapse = "")
}


betaF <- function(){
  paste0(beta, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "beta.txt")

write.table(betaF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


betaVar <- "string"
for (jj in 1:Class) {
  betaVar[jj] <- paste0("var.betas", jj, " <- rep(con$K, ncX)", sep = "", collapse = "")
}


betaVarF <- function(){
  paste0(betaVar, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "betaVar.txt")

write.table(betaVarF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


al <- "string"
for (jj in 1:Class) {
  al[jj] <- paste0("alphas", jj, " <- 0", sep = "", collapse = "")
}


alF <- function(){
  paste0(al, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "alphas.txt")

write.table(alF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


alVar <- "string"
for (jj in 1:Class) {
  alVar[jj] <- paste0("var.alphas", jj, " <- 100", sep = "", collapse = "")
}


alVarF <- function(){
  paste0(alVar, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "alphasVar.txt")

write.table(alVarF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source("Adjusting code\\alphasVar.txt")

if (fixedGammas == TRUE) {
  gam <- paste0("gammas <- rep(0,(ncW))", sep = "", collapse = "")
} else {
 gam <- "string"
 for (jj in 1:Class) {
   gam[jj] <- paste0("gammas", jj, " <- rep(0,(ncW))", sep = "", collapse = "")
 }
}


gamF <- function(){
  paste0(gam, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "gammas.txt")

write.table(gamF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)

if (fixedGammas == TRUE) {
gamVar <- paste0("var.gammas <- rep(con$K, (ncW))", sep = "", collapse = "")
} else {
 gamVar <- "string"
 for (jj in 1:Class) {
   gamVar[jj] <- paste0("var.gammas", jj, " <- rep(con$K, (ncW))", sep = "", collapse = "")
 }
}

gamVarF <- function(){
  paste0(gamVar, "\n", sep = "", collapse = "")
}

filename <- file.path("Adjusting code", "gammasVar.txt")

write.table(gamVarF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


if (fixedBsgammas == TRUE) {
 Bs.gam <- paste0("Bs.gammas <- rep(0, (ncWBH))", sep = "", collapse = "")
} else {
 Bs.gam <- "string"
 for (jj in 1:Class) {
   Bs.gam[jj] <- paste0("Bs.gammas", jj, " <- rep(0, (ncWBH))", sep = "", collapse = "")
 }
}

Bs.gamF <- function(){
  paste0(Bs.gam, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "Bs.gammas.txt")

write.table(Bs.gamF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


if (fixedBsgammas == TRUE) {
 Bs.gamVar <- paste0("var.Bs.gammas <- rep(con$K, (ncWBH))", sep = "", collapse = "")
} else {
 Bs.gamVar <- "string"
 for (jj in 1:Class) {
   Bs.gamVar[jj] <- paste0("var.Bs.gammas", jj, " <- rep(con$K, (ncWBH))", sep = "", collapse = "")
 }
}


Bs.gamVarF <- function(){
  paste0(Bs.gamVar, "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "Bs.gammasVar.txt")

write.table(Bs.gamVarF(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)



mu0. <- "string"
for (jj in 1:Class) {
  mu0.[jj] <- paste0("mu0", jj, " <- rep(0,(ncZ))", sep = "", collapse = "")
}


mu0F <- function(){
  paste0(mu0., "\n", sep = "", collapse = "")
}


filename <- file.path("Adjusting code", "mu0.txt")

write.table(mu0F(), filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)

b <- cbind(data.matrix(ranef(lmeObject)))


nY <- nrow(b)

if (family == "gaussian") {
  sigma2 <- lmeObject$sigma^2
}

#####################################################################################################

#Data for jags
dataMCMCi <- paste0("N = nY, K = K, offset = offset, X = X, Xtime = Xtime, ncX = ncol(X), ncZ = ncol(Z), 
                    y = y$y, 
                    Xs = Xs, 
                    Z = Z, Ztime = Ztime,  
                    Zs = Zs, 
                    event = eventD, zeros = zeros,  
                    C = C, P = P,
                    wk = wk, nb = nb, 
                    W = x$W,
                    ncW = ncol(x$W), 
                    WBH = WBH, 
                    ncWBH = ncol(x$WBH),
                    WBHs = WBHs, 
                    ncWBHs = ncol(x$WBH),",
                    if (family == "gaussian") {
                      paste0("priorA.tau = 0.01,
                    priorB.tau = 0.01, ")
                    },
                    "prior.cl = rep(", a, ",", jj, ")", sep = "", collapse = "")

datamu<- "string"
for (jj in 1:Class) {
  datamu[jj] <- paste0("mu0", jj, " = mu0", jj, sep = "", collapse = "")
}


datainvD<- "string"
for (jj in 1:Class) {
  datainvD[jj] <- paste0("priorR.D", jj, " = diag(0.01, ncZ), priorK.D", jj, " = (ncZ)", sep = "", collapse = "")
}


databet<- "string"
for (jj in 1:Class) {
  databet[jj] <- paste0("priorMean.betas", jj, " = betas", jj, sep = "", collapse = "")
}

databetVar <- "string"
for (jj in 1:Class) {
  databetVar[jj] <- paste0("priorTau.betas", jj, " = diag(1/var.betas", jj, ")", sep = "", collapse = "")
}

if (fixedGammas == TRUE){
datagam <- paste0("priorMean.gammas = gammas", sep = "", collapse = "")
} else {
 datagam<- "string"
 for (jj in 1:Class) {
   datagam[jj] <- paste0("priorMean.gammas", jj, " = gammas", jj, sep = "", collapse = "")
 }
}

if (fixedGammas == TRUE){
  datagamVar <- paste0("priorTau.gammas = diag(1/var.gammas)", sep = "", collapse = "")
} else {
 datagamVar <- "string"
 for (jj in 1:Class) {
   datagamVar[jj] <- paste0("priorTau.gammas", jj, " = diag(1/var.gammas", jj, ")", sep = "", collapse = "")
 }
}

dataal <- "string"
for (jj in 1:Class) {
  dataal[jj] <- paste0("priorMean.alphas", jj, " = alphas", jj, sep = "", collapse = "")
}

dataalVar <- "string"
for (jj in 1:Class) {
  dataalVar[jj] <- paste0("priorTau.alphas", jj, " = 1/var.alphas", jj, sep = "", collapse = "")
}


if (fixedBsgammas == TRUE){
  databsgam <- paste0("priorMean.Bs.gammas = Bs.gammas", sep = "", collapse = "")
} else {
  databsgam<- "string"
  for (jj in 1:Class) {
    databsgam[jj] <- paste0("priorMean.Bs.gammas", jj, " = Bs.gammas", jj, sep = "", collapse = "")
  }
}

if (fixedBsgammas == TRUE){
  databsgamVar <- paste0("priorTau.Bs.gammas  = diag(1/var.Bs.gammas)", sep = "", collapse = "")
} else {
  databsgamVar <- "string"
  for (jj in 1:Class) {
    databsgamVar[jj] <- paste0("priorTau.Bs.gammas", jj, " = diag(1/var.Bs.gammas", jj, ")", sep = "", collapse = "")
  }
}


datamuF <- function(){
  paste0(datamu, sep = "", collapse = ", ")
}
datainvDF <- function(){
  paste0(datainvD, sep = "", collapse = ", ")
}
databetF <- function(){
  paste0(databet, sep = "", collapse = ", ")
}
databetVarF <- function(){
  paste0(databetVar, sep = "", collapse = ", ")
}
datagamF <- function(){
  paste0(datagam, sep = "", collapse = ", ")
}
datagamVarF <- function(){
  paste0(datagamVar, sep = "", collapse = ", ")
}
dataalF <- function(){
  paste0(dataal, sep = "", collapse = ", ")
}
dataalVarF <- function(){
  paste0(dataalVar, sep = "", collapse = ", ")
}
dataBsgamF <- function(){
  paste0(databsgam, sep = "", collapse = ", ")
}
dataBsgamVarF <- function(){
  paste0(databsgamVar, sep = "", collapse = ", ")
}


dataMCMC <- paste0("Data <- list(", dataMCMCi, ",\n", datamuF(), ",\n", datainvDF(), ",\n", databetF(), ",\n", databetVarF(),
                   ",\n", datagamF(), ",\n", datagamVarF(), ",\n", if (method == "JM") { paste0(dataalF(), ",\n", dataalVarF(), ",\n")  },
                   dataBsgamF(), ",\n", dataBsgamVarF(), ")", sep = "", collapse = "" )


filename <- file.path("Adjusting code", "MCMCdata.txt")

write.table(dataMCMC, filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)


# Paramteras to monitor
parbet <- "string"
for (jj in 1:Class) {
  parbet[jj] <- paste0("\"betas", jj, "\"", sep = "", collapse = "")
}

# parinvD <- "string"
# for (jj in 1:Class) {
#   parinvD [jj] <- paste0("\"inv.D", jj, "\"", sep = "", collapse = "")
# }
parb <- "string"
for (jj in 1:Class) {
  parb[jj] <- paste0("\"b", jj, "\"", sep = "", collapse = "")
}


if (fixedGammas == TRUE) {
   pargam <- paste0("\"gammas\"", sep = "", collapse = "")
} else {
 pargam <- "string"
 for (jj in 1:Class) {
   pargam[jj] <- paste0("\"gammas", jj, "\"", sep = "", collapse = "")
 }
}


paral <- "string"
for (jj in 1:Class) {
  paral[jj] <- paste0("\"alphas", jj, "\"", sep = "", collapse = "")
}

if (fixedBsgammas == TRUE) {
   parbsgam <- paste0("\"Bs.gammas\"", sep = "", collapse = "")
} else {
 parbsgam <- "string"
 for (jj in 1:Class) {
   parbsgam[jj] <- paste0("\"Bs.gammas", jj, "\"", sep = "", collapse = "")
 }
}



parbetF <- function(){
  paste0(parbet, sep = "", collapse = ", ")
}

# parinvDF <- function(){
#   paste0(parinvD, sep = "", collapse = ", ")
# }
parbF <- function(){
  paste0(parb, sep = "", collapse = ", ")
}
paralF <- function(){
  paste0(paral, sep = "", collapse = ", ")
}
pargamF <- function(){
  paste0(pargam, sep = "", collapse = ", ")
}
parbsgamF <- function(){
  paste0(parbsgam, sep = "", collapse = ", ")
}



parms <- paste0("parms <- c(", parbetF(), ",", if (family == "gaussian") {paste0("\"tau\"",",")}, {paste0("\"inv.D\"",",")}, 
                pargamF(), ", ", parbF(), ", ", 
                if (method == "JM") {paste0(paralF(), ", ")}, parbsgamF(), ", ",
                "\"pr\"", ", ", "\"v\"", " )")



filename <- file.path("Adjusting code", "params.txt")

write.table(parms, filename, append = FALSE, quote = FALSE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))


source(filename)

#####################
JAGSmodel(Class, family, method, hc, fixedGammas, fixedBsgammas, fixedInvD, RM_method)
#####################

Data$inv.D <- matrix(0, nb, nb)
diag(Data$inv.D) <- NA


# Model in JAGS
model.fit <- jags.model(file = "JMsuper.txt", data = Data, n.chains = con$n.chains, 
                        n.adapt = con$n.adapt, quiet = con$quiet)


update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
codaFit <- as.mcmc.list(res)


