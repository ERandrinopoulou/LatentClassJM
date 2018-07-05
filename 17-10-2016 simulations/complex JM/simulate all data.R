rm(list=ls(all=TRUE))

setwd("G:/2014-10 PostDoc/super JM/Latent class JM/17-10-2016 simulations/complex JM")

library(JMbayes)
library(rjags)
library(xtable)

ClassSim <- 2
Class <- 2
family <- "gaussian"
method <- "JM"
hc <- FALSE
fixedGammas <- FALSE
fixedBsgammas <- FALSE
RM_method <- TRUE


a <- 6.9

source("Functions.R")

dir.create(paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a ))
#################################################
#################################################
M <- 25

betas1mean <- matrix(0, M, Class*3)
tau1mean <- numeric(M)
gammasmean <- matrix(0, M, Class*2)
alphas1mean <- matrix(0, M, Class)
class <- list(M) 
weight <- list(M)
listRes<- list(M)
conv <- list(M)

for (l in 1:M) {
  
  print(l)
  set.seed(l)
  
  ########################
  
  source("simulate class 1.R")
  
  data_1 <- data
  data.id_1 <- data.id
  
  data <- rbind(data_1)
  data.id <- rbind(data.id_1)
  
  #######################
  if (ClassSim == 2 | ClassSim == 3) {
    source("simulate class 2.R")
    
    data_2 <- data
    data.id_2 <- data.id
    data_2$IDnr <- data_2$IDnr + tail(data.id_1$IDnr, n=1)
    data.id_2$IDnr <- data.id_2$IDnr + tail(data.id_1$IDnr, n=1)
    
    data <- rbind(data_1, data_2)
    data.id <- rbind(data.id_1, data.id_2)
  }
  
  
  #######################
  if (ClassSim == 3) {
    source("simulate class 3.R")
    
    data_3 <- data
    data.id_3 <- data.id
    data_3$IDnr <- data_3$IDnr + tail(data.id_2$IDnr, n=1)
    data.id_3$IDnr <- data.id_3$IDnr + tail(data.id_2$IDnr, n=1)
    
    
    data <- rbind(data_1, data_2, data_3)
    data.id <- rbind(data.id_1, data.id_2, data.id_3)
  }
  
  ########################
  try(source("model run.R"))
  ########################
  
  g <- geweke.diag(codaFit)
  
  parms1 <- parms[parms %!in% c("b1", "b2", "v", "pr")]
  
  bss <- g[[1]]$z
  n.sims <- nrow(bss)
  Zsc <- vector("list", length(parms1))
  names(Zsc) <- parms1
  for (p in seq_along(parms1)) {
    ii <- grep(paste("^", parms1[p], sep = ""), names(bss))
    Zsc[[p]] <- bss[ii]
  }
  
  conv[[l]] <- 2*pnorm(-abs(unlist(Zcs))) < 0.05

  
  
  bss <- do.call(rbind,codaFit)
  colnames(bss)
  n.sims <- nrow(bss)
  sims.list <- vector("list", length(parms))
  names(sims.list) <- parms
  for (p in seq_along(parms)) {
    ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
    sims.list[[p]] <- bss[, ii]
  }
  
  ### betas
  betas <- "string"
  for (jj in 1:Class) {
    betas[jj] <- paste0("sims.list$betas", jj, sep = "", collapse = "")
  }
  betasF <- function(){
    paste0(betas, sep = "", collapse = ", ")
  }
  
  sims.list$betas <- paste0("sims.list$betas <- cbind(",betasF(), ")")
  filename <- file.path("Adjusting code", "betasRes.txt")
  
  write.table(sims.list$betas, filename, append = FALSE, quote = FALSE, sep = " ", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = FALSE, qmethod = c("escape", "double"))
  
  source("Adjusting code\\betasRes.txt")
  
  if (Class %in% 1:9) betas1mean[l,] <- apply(sims.list$betas,2,mean)  
  if (Class == 10)  betas1mean[l,] <- apply(sims.list$betas,2,mean)[c(1:3,7:33)]
  
  ### tau
  tau1mean[l] <- mean(sims.list$tau)          
  
  
  ### gammas
  
  if (fixedGammas == TRUE) {
    gammas <- paste0("sims.list$gammas", sep = "", collapse = "")
  }  else {
    gammas <- "string"
    for (jj in 1:Class) {
      gammas[jj] <- paste0("sims.list$gammas", jj, sep = "", collapse = "")
    }
  }
  
  gammasF <- function(){
    paste0(gammas, sep = "", collapse = ", ")
  }
  
  sims.list$gammas <- paste0("sims.list$gammas <- cbind(",gammasF(), ")")
  filename <- file.path("Adjusting code", "gammasRes.txt")
  
  write.table(sims.list$gammas, filename, append = FALSE, quote = FALSE, sep = " ", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
              col.names = FALSE, qmethod = c("escape", "double"))
  
  source("Adjusting code\\gammasRes.txt")
  
  if (Class %in% 1:9) gammasmean[l,] <- apply(sims.list$gammas,2,mean)
  if (Class == 10) gammasmean[l,] <- apply(sims.list$gammas,2,mean)[c(1:2,5:22)]  
  
  if (method == "JM") {  
    ### alphas
    alphas <- "string"
    for (jj in 1:Class) {
      alphas[jj] <- paste0("sims.list$alphas", jj, sep = "", collapse = "")
    }
    
    alphasF <- function(){
      paste0(alphas, sep = "", collapse = ", ")
    }
    
    sims.list$alphas <- paste0("sims.list$alphas <- cbind(",alphasF(), ")")
    filename <- file.path("Adjusting code", "alphasRes.txt")
    
    write.table(sims.list$alphas, filename, append = FALSE, quote = FALSE, sep = " ", 
                eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
                col.names = FALSE, qmethod = c("escape", "double"))
    
    
    source("Adjusting code\\alphasRes.txt")
    
    if (Class %in% 1:9) alphas1mean[l,] <- apply(sims.list$alphas,2,mean)   
    if (Class == 10) alphas1mean[l,] <- apply(sims.list$alphas,2,mean)[c(1,3:11)]
  }  
  
  ### v
  class[[l]] <- sims.list$v  
  
  weight[[l]] <- sims.list$pr
  
  listRes[[l]] <- list(betas1mean[l,], tau1mean[l], gammasmean[l,], alphas1mean[l,], class[[l]], weight[[l]], conv[[l]])
  ress <- listRes[[l]]
  
  fileout <- paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a ,"/results",l,".RData")
  save(ress, file = fileout)
}





###############################################################################
############################################################################### 
######################## Obtain and present results ###########################
############################# NOT READY YET ###################################
###############################################################################
M <- 6

betas <- matrix(, M, Class*3)
tau <- numeric(M)
gammas <- matrix(, M, Class*2)
alphas <- matrix(, M, Class)
cl <- matrix(, M, Class)
conver <- list(M)

for (i in 1:M) {
  fileout <- paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a ,"/results",i,".RData")
  load(fileout)
  
  betas[i, ] <- as.numeric(ress[[1]])
  tau[i] <- as.numeric(ress[[2]])
  gammas[i, ] <- as.numeric(ress[[3]])
  alphas[i, ] <- as.numeric(ress[[4]])
  conver[[i]] <- ress[[7]]
  
  #cl[i, ] <- as.numeric(res[[5]])

  #if (Class == 2) {
  #  if (cl[i, 1] < 0.7) { ##### label switching needs to be fixed per iteration
  #    betas[i,] <- betas[i, c(4, 5, 6, 1, 2, 3)]
  #    gammas[i,] <- gammas[i, c(3, 4, 1, 2)]
  #    alphas[i,] <- alphas[i, c(2,1)]
  #  }
  #}
  
}

array(c(betas, gammas, alphas, c(,,3)))
array(c(as.matrix(betas), as.matrix(gammas), as.matrix(alphas)))

aic(betas, 1)

apply(betas, 2, mean)
apply(gammas, 2, mean)
apply(alphas, 2, mean)
apply(cl, 2, mean)


betas <- betas[-c(3,5,7,8,21,22,25,27,28, 32,35,37,40,45), ]
alphas <- alphas[-c(3,5,7,8,21,22,25,27,28, 32,35,37,40,45), ]

apply(betas, 2, mean)
apply(alphas, 2, mean)

par(mfrow= c(2,3))
  boxplot(betas[,1], ylim = c(7.5,8.5))
  abline(h = 8.0317, col = "red")

  boxplot(betas[,2], ylim = c(-4.3,-6.5))
  abline(h = -5.8568, col = "red")
  
  boxplot(betas[,3], ylim = c(-0.19,0.10))
  abline(h = -0.1578, col = "red")
  
  boxplot(betas[,4], ylim = c(-7.5,-8.5))
  abline(h = -8.0317, col = "red")
  
  boxplot(betas[,5], ylim = c(9,15))
  abline(h = 12.2066, col = "red")
  
  boxplot(betas[,6], ylim = c(-0.19,1))
  abline(h = 0.4578, col = "red")
  
  
par(mfrow= c(2,3))
  boxplot(betas[,1])
  abline(h = 8.0317, col = "red")
  
  boxplot(betas[,2])
  abline(h = -5.8568, col = "red")
  
  boxplot(betas[,3])
  abline(h = -0.1578, col = "red")
  
  boxplot(betas[,4])
  abline(h = -8.0317, col = "red")
  
  boxplot(betas[,5])
  abline(h = 12.2066, col = "red")
  
  boxplot(betas[,6])
  abline(h = 0.4578, col = "red")
  
par(mfrow= c(1,2))
  boxplot(alphas[,1])
  abline(h = 0.3833, col = "red")
  
  boxplot(alphas[,2])
  abline(h = 0.0833, col = "red")
  
  

boxplot(betas, ylim = c(-10, 15))

# create data for segments
n <- ncol(betas)
# width of each boxplot is 0.8
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
# these are the y-coordinates for the horizontal lines
# that you need to set to the desired values.
y0s <- c(8.0317, -5.8568, -0.1578, -8.0317, 12.2066, 0.4578)

# add segments
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red")
