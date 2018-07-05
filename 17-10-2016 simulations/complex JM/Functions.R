'%!in%' <- function(x,y)!('%in%'(x,y))


tpower <- function(x, t, p)
  # Function for truncated p-th power function
  (x - t) ^ p * (x > t)


bbase <- function(x, xl, xr, ndx, deg){
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B 
}


knots <- function(xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  seq(xl - deg * dx, xr + deg * dx, by = dx)
}


JAGSmodel <- function(jj, family, method, hc, fixedGammas, fixedBsgammas, RM_method){
  
  myt <- function(times = 1) {
    tb <- "    "
    paste(rep(tb, times), collapse = "")
  } 
  
  longbetas <- "string" 
  for (i in 1:jj) {
    longbetas[i] <- paste0("(equals(v[i], ", i, ") * betas", i, "[1:ncX])", sep = "", collapse = "")  
  }
  
  longb <- "string" 
  for (i in 1:jj) { 
    longb[i] <- paste0("(equals(v[i], ", i, ") * b", i, "[i, 1:ncZ])", sep = "", collapse = "")  
  } 
  
  longu <- "string" 
  for (i in 1:jj) { 
    longu[i] <- paste0("(equals(v[i], ", i, ") * u", i, "[i, 1:ncZ])", sep = "", collapse = "")  
  }
  
  longCont <- function() {
    long1 <- paste0(myt(), "for (j in offset[i]:(offset[i+1] - 1)) {\n")
    long2 <- paste0(myt(2), "muy[j] <- inprod(", paste0(longbetas, collapse = " + "), ", X[j, 1:ncX]) +\n", myt(2)," inprod(", paste0(longb, sep = "", collapse = "+ "), ", Z[j, 1:ncZ] )\n")                                                       
    
    long2hier <- paste0(myt(2), "muy[j] <- inprod(", paste0(longu, sep = "", collapse = " + "), ", Z[j, 1:ncZ])\n")                                                       
    
    long3 <- paste0(myt(2), "y[j] ~ dnorm(muy[j], tau)\n")
    long4 <- paste0(myt(), "}\n")  
    paste0(long1, if (hc == "TRUE") {paste0(long2hier)} else {paste0(long2)}, long3, long4)
  }
  
  longBin <- function() {
    long1 <- paste0(myt(), "for (j in offset[i]:(offset[i+1] - 1)) {\n" )
    long2 <- paste0(myt(2), "muy[j] <- inprod(", paste0(longbetas, collapse = " + "), ", X[j, 1:ncX]) +\n", myt(2)," inprod(", paste0(longb, sep = "", collapse = "+ "), ", Z[j, 1:ncZ] )\n")                                                       
    
    long2hier <- paste0(myt(2), "muy[j] <- inprod(", paste0(longu, sep = "", collapse = "+ "), ", Z[j, 1:ncZ] )\n")                                                       
    
    long3 <- paste0(myt(2), "Pr[j] <- max(1.00000E-05, min(0.99999, (exp(muy[j])/(1 + exp(muy[j])))))\n")
    long4 <- paste0(myt(2), "y[j] ~ dbin(Pr[j], 1)\n")
    long5 <- paste0(myt(), "}\n", collapse = NULL)  
    paste0(long1, if (hc == "TRUE") {paste0(long2hier)} else {paste0(long2)}, long3, long4, long5)
  }
  
  if (family == "gaussian") {
    long <- longCont
  } 
  if (family == "binomial") {
    long <- longBin
  } 
  
  if (fixedGammas == TRUE) {
    gammas <- paste0("gammas[1:ncW]", sep = "", collapse = "")  
  } else {
   gammas <- "string" 
   for (i in 1:jj) {
     gammas[i] <- paste0("(equals(v[i], ", i, ") * gammas", i, "[1:ncW])", sep = "", collapse = "")  
   }
  }
  
  if (fixedBsgammas == TRUE) {
    Bs.gammas <- paste0("Bs.gammas[1:ncWBH]", sep = "", collapse = "")  
  } else {
   Bs.gammas <- "string"
   for (i in 1:jj) {
     Bs.gammas[i] <-  paste0("(equals(v[i], ", i, ") * Bs.gammas", i, "[1:ncWBH])", sep = "", collapse = "")
   }
  }
  
  alphas <- "string"
  for (i in 1:jj) {
    alphas[i] <-  paste0("(equals(v[i], ", i, ") * alphas", i, ")", sep = "", collapse = "")
  }
  
  
  hazJM <- function() {
    fT <- paste0(myt(), "f.T[i] <- inprod(", paste0(longbetas, sep = "", collapse = " + "),", Xtime[i, 1:ncX]) + inprod(", paste0(longb, sep = "", collapse = "+ "), ", Ztime[i, 1:ncZ])\n", sep = "", collapse = NULL)
    etaBas <- paste0(myt(), "etaBaseline[i] <- inprod(", paste0(gammas, sep = "", collapse = " + "), ", W[i, 1:ncW])\n", collapse = "")
    log.h0 <- paste0(myt(), "log.h0.T[i] <- inprod(", paste0(Bs.gammas, sep = "", collapse = " + "), ", WBH[i, 1:ncWBH])\n", collapse = "")
    log.haz <- paste0(myt(), "log.hazard[i] <- log.h0.T[i] + etaBaseline[i] + (", paste0(alphas, sep = "", collapse = " + "), ") * f.T[i]\n", collapse = "")
    paste0(etaBas, log.h0, fT, log.haz)
  }
  
  hazLCJM <- function() {
    etaBas <- paste0(myt(), "etaBaseline[i] <- inprod(", paste0(gammas, sep = "", collapse = " + "), ", W[i, 1:ncW])\n", collapse = "")
    log.h0 <- paste0(myt(), "log.h0.T[i] <- inprod(", paste0(Bs.gammas, sep = "", collapse = " + "), ", WBH[i, 1:ncWBH])\n", collapse = "")
    log.haz <- paste0(myt(), "log.hazard[i] <- log.h0.T[i] + etaBaseline[i]\n", collapse = "")
    paste0(etaBas, log.h0, log.haz)
  }
  
  
  if (method == "JM") {
    haz <- hazJM
  } 
  if (method == "LCJM") {
    haz <- hazLCJM
  } 
  
  if (RM_method == TRUE) {
     c <- paste0(myt(), "v[i] ~ dcat(pr[]) \n")  
  } else {
     c <- paste0(myt(), "v[i] ~ dcat(pr[i, ]) \n")  
  }
  
  
  prior.b <- "string"
  for (i in 1:jj) {
    prior.b[i] <- paste0(myt(), "b", i, "[i, 1:nb] ~ dmnorm(mu0", i, "[], inv.D", i, "[, ])\n")  
  }
  
  
  #### D
  if (hc == TRUE) {
    
    hierlong <- "string"
    
    for (j in seq_along(colmns_HC)) {
      ii <- colmns_HC[[j]]
      j.incr <- j + 0
      longbetasN <- "string" 
      
      if (length(ii) > 1) {
        ind.cl <- paste0("c(", paste(ii, collapse = ", "), ")")
        for (i in 1:jj) {
          longbetasN[i] <- paste0("equals(v[i], ", i, ") * betas", i, "[", ind.cl, "]", sep = "", collapse = "")  
        }
        
        hierlong[j] <- paste0(myt(1), "mu.u[i, ", j.incr, "] <- inprod((",  paste0(longbetasN, collapse = " + "),
                              "), Xhc[i, ", ind.cl, "])\n", sep = "", collapse = "")
      } else {
        ind.cl <- ii
        for (i in 1:jj) {
          longbetasN[i] <- paste0("equals(v[i], ", i, ") * betas", i, "[", ind.cl, "]", sep = "", collapse = "")  
        }
        
        hierlong[j] <- paste0(paste0(myt(1), "mu.u[i, ", j.incr, "] <- ", paste0(longbetasN, collapse = " + "), "\n", 
                                     sep = "", collapse = ""))
      }
    }
    
    
    uf <- function(j){
      u <- "" 
      for (i in 1:jj) {
        u[i] <- paste0("equals(v[i], ", i, ") * u", i, "[i, ", j, "]" )
      }
      paste0(u)
    }
    
    prior.b.hier <- "" 
    for (i in 1:jj) {
      for (j in seq_along(colmns_HC)) {
        prior.b.hier <- paste0(prior.b.hier, myt(), "b", i, "[i, ", j, "] <- (",  paste0(uf(j), collapse = " + "), ") - mu.u[i, ", j, "]\n", sep = "", collapse = "")  
      }
    }
    
    
    prior.b.hier2 <- "string"
    for (i in 1:jj) {
      prior.b.hier2[i] <- paste0(myt(), "u", i, "[i, 1:nb] ~ dmnorm(mu.u[i, ], inv.D", i, "[, ]) \n") 
    }
    
  }
  ####
  
  survJM <- function() {
    surv1 <- paste0(myt(), "for (k in 1:K) {\n")
    fTs <- paste0(myt(2), "f.s[i, k] <- inprod(", paste0(longbetas, sep = "", collapse = " + "),", Xs[K * (i-1) + k, 1:ncX]) + inprod(", paste0(longb, sep = "", collapse = " + "), ", Zs[K * (i-1) + k, 1:ncZ])\n", sep = "", collapse = NULL)
    surv2 <- paste0(myt(), "}\n")
    phi <- paste0(myt(), "phi[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]\n", collapse = "")
    surv4 <- paste0(myt(), "zeros[i] ~ dpois(phi[i])\n")
    log.survival <- paste0(myt(), "log.survival[i] <- -exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])\n")
    c <- paste0(c, sep = "", collapse = " ")
    b <- paste0(prior.b, collapse = "")
    
    if (hc == "TRUE") {
      bhier <- paste0(prior.b.hier, collapse = "")  
      hierlong.mu <- paste(hierlong, collapse = "")
      prior.b.hier2 <- paste0(prior.b.hier2, collapse = "")  
    }
    
    log.h0 <- paste0(myt(2), "log.h0.s[i, k] <- inprod(", paste0(Bs.gammas, sep = "", collapse = " + "), ", WBHs[K * (i-1) + k, 1:ncWBHs])\n", collapse = "")
    survLong <- paste0(myt(2), "SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k] + (", paste0(alphas, sep = "", collapse = " + "), ") * f.s[i, k])\n", collapse = "")
    log.surv <- paste0(log.survival, sep = "", collapse = " ")
    paste0(surv1, fTs, log.h0, survLong, surv2, log.survival, phi, surv4, c, if (hc == "TRUE") {paste0(hierlong.mu)}, if (hc == "TRUE") {paste0(bhier)} else {paste0(b)},
           if (hc == "TRUE") {paste0(prior.b.hier2)})
  }
  
  survLCJM <- function() {
    surv1 <- paste0(myt(), "for (k in 1:K) {\n")
    surv2 <- paste0(myt(), "}\n")
    phi <- paste0(myt(), "phi[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]\n", collapse = "")
    surv4 <- paste0(myt(), "zeros[i] ~ dpois(phi[i])\n")
    log.survival <- paste0(myt(), "log.survival[i] <- -exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])\n")
    c <- paste0(c, sep = "", collapse = " ")
    b <- paste0(prior.b, collapse = "")
    
    if (hc == "TRUE") {
    bhier <- paste0(prior.b.hier, collapse = "")
    prior.b.hier2 <- paste0(prior.b.hier2, collapse = "") 
    }
    
    log.h0 <- paste0(myt(2), "log.h0.s[i, k] <- inprod(", paste0(Bs.gammas, sep = "", collapse = " + "), ", WBHs[K * (i-1) + k, 1:ncWBHs])\n", collapse = "")
    survLong <- paste0(myt(2), "SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k])\n", collapse = "")
    log.surv <- paste0(log.survival, sep = "", collapse = " ")
    paste0(surv1, log.h0, survLong, surv2, log.survival, phi, surv4, c, if (hc == "TRUE") {paste0(bhier)} else {paste0(b)})
  }
  
  
  if (method == "JM") {
    surv <- survJM
  } 
  if (method == "LCJM") {
    surv <- survLCJM
  } 
  
  

  if (RM_method == TRUE) {
       class1 <- paste0(myt(), "pr[1:", jj, "] ~ ddirch(prior.cl[]) \n", sep = "", collapse = "")
  } else {
       class1 <- paste0(myt(), "pr[i, 1:", jj, "] ~ ddirch(prior.cl[]) \n", sep = "", collapse = "")
  }


  prior.betas <- "string"
  for (i in 1:jj) {
    prior.betas[i] <- paste0(myt(), "betas", i, "[1:(ncX)] ~ dmnorm(priorMean.betas", i, "[], priorTau.betas", i, "[, ])\n")  
  }
  
  if (fixedGammas == TRUE) {
    prior.gammas <- paste0(myt(), "gammas[1:(ncW)] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])\n")  
  } else {
   prior.gammas <- "string"
   for (i in 1:jj) {
     prior.gammas[i] <- paste0(myt(), "gammas", i, "[1:(ncW)] ~ dmnorm(priorMean.gammas", i, "[], priorTau.gammas", i, "[, ])\n")  
   }
  }
  
  prior.alphas <- "string"
  for (i in 1:jj) {
    prior.alphas[i] <- paste0(myt(), "alphas", i, " ~ dnorm(priorMean.alphas", i, ", priorTau.alphas", i, ")\n")  
  }
  
  if (fixedBsgammas == TRUE) {
    prior.Bs.gammas <- paste0(myt(), "Bs.gammas[1:(ncWBH)] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, ])\n")  
  } else {
   prior.Bs.gammas <- "string"
   for (i in 1:jj) {
     prior.Bs.gammas[i] <- paste0(myt(), "Bs.gammas", i, "[1:(ncWBH)] ~ dmnorm(priorMean.Bs.gammas", i, "[], priorTau.Bs.gammas", i, "[, ])\n")  
   }
  }
  
  inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
  
  prior.invD <- "string"
  for (i in 1:jj) {
    prior.invD[i] <- paste0(myt(), "inv.D", i, "[1:nb, 1:nb] ~ dwish(priorR.D", i, "[, ], priorK.D", i, ")\n")  
  }
  
  
  
  priors <- function() {
    prior.betas <- paste0(prior.betas, collapse = "")
    prior.tau <- paste0(myt(), "tau ~ dgamma(priorA.tau, priorB.tau)\n", collapse = "")
    prior.invD <-  paste0(prior.invD, collapse = "")
    prior.gam <- paste0(prior.gammas, collapse = "")
    prior.als <- paste0(prior.alphas, collapse = "")
    prior.Bs.gam <- paste0(prior.Bs.gammas, collapse = "")
    if (RM_method == TRUE) {
      classtem <- paste0(class1, sep = "", collapse = " ")
      paste0(classtem, prior.betas, if (family == "gaussian") {paste0(prior.tau)}, prior.invD, prior.gam, if (method == "JM") {paste0(prior.als)}, prior.Bs.gam)
    } else {
      paste0(prior.betas, if (family == "gaussian") {paste0(prior.tau)}, prior.invD, prior.gam, if (method == "JM") {paste0(prior.als)}, prior.Bs.gam)
    }
  }
  
  
  model <- function() {
    beg <- paste0("model {\n for (i in 1:N) {\n",  sep = "", collapse = "")
    end <- paste0("}\n", sep = "", collapse = "")
    if (RM_method == TRUE) {
      paste0(beg, long(), haz(), surv(), " }\n", priors(), end)
    } else {
      classtem <- paste0(class1, sep = "", collapse = " ")
      paste0(beg, long(), haz(), surv(), classtem, " }\n", priors(), end) 
   }
  }
  
  
  model <- model()
  
  
  filename <- file.path("JMsuper.txt")
  
  write.table(model, filename, append = FALSE, quote = FALSE, sep = " ", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = FALSE, qmethod = c("escape", "double"))
  
  
}

#############################################################################
extractFrames <- function (formula, data) {
  Terms <- terms(formula)
  term_labels <- attr(Terms, "term.labels")
  which_RE <- grep("|", term_labels, fixed = TRUE)
  namesVars <- all.vars(formula)
  respVar <- as.character(formula)[2L]
  # Fixed Effects
  formYx <- paste(term_labels[-which_RE], collapse = " + ")
  formYx <- as.formula(paste(respVar, "~", formYx))
  TermsX <- terms(formYx, data = data)
  mfX <- model.frame(TermsX, data)
  X <- model.matrix(TermsX, data)
  # Random Effects
  spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
  idVar <- spl[2L]
  data <- data[complete.cases(data[namesVars]), ]
  id <- data[[idVar]]
  id <- match(id, unique(id))
  formYz <- paste(spl[1], collapse = " + ")
  formYz <- as.formula(paste(respVar, "~", formYz))
  TermsZ <- terms(formYz, data = data)
  Z <- model.matrix(TermsZ, data)
  # response variable
  y <- model.response(mfX)
  if (is.factor(y))
    y <- as.vector(unclass(y) - 1)
  # hierarchical centering
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  has_interceptX <- attr(TermsX, "intercept")
  has_interceptZ <- attr(TermsZ, "intercept")
  performHC <- has_interceptX && (has_interceptX == has_interceptZ)
  if (performHC) {
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z)) grep(terms.labs_Z[1L], colnames(X), fixed = TRUE)
    which_td <- unname(which(apply(X, 2, check_td, id = id)))
    all_TDterms <- unique(c(timeTerms, which_td))
    baseline <- seq_len(ncol(X))[-all_TDterms]
    ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions,
                                           nams2 = colnames(X)))
    ind_colmns2 <- seq_len(ncol(X))
    ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    data.id <- data[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
      data.id[[terms.labs_Z[1L]]] <- 1
      mfHC <- model.frame(TermsX, data = data.id)
      which.timeVar <- grep(terms.labs_Z[1L], names(mfHC), fixed = TRUE)
      mfHC[which.timeVar] <- lapply(mfHC[which.timeVar],
                                    function (x) { x[] <- 1; x })
      model.matrix(formYx, mfHC)
    } else {
      model.matrix(formYx, model.frame(TermsX, data = data.id))
    }
  }
  # extract results
  list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
       id = id, y = y, X = X, Z = Z, TermsX = TermsX,
       TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
       Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
       ncx = ncol(X), ncz = ncol(Z))
}




