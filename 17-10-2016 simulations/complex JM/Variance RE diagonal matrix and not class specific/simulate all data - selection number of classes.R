rm(list=ls(all=TRUE))

library(JMbayes)
library(rjags)
library(xtable)
library(splines)

ClassSim <- 3
Class <- 6
family <- "gaussian"
method <- "JM"
hc <- FALSE
fixedGammas <- FALSE
fixedBsgammas <- FALSE
fixedInvD <- TRUE
RM_method <- TRUE


a <- 6.9

source("Functions.R")

dir.create(paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a ))
#################################################
#################################################
M <- 150

betas1mean <- matrix(0, M, Class*3)
tau1mean <- numeric(M)
gammasmean <- matrix(0, M, Class*2)
alphas1mean <- matrix(0, M, Class)
class <- list(M) 
weight <- list(M)
conv <- list(M)
listRes<- list(M)

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
  
  parms1 <- parms[parms %!in% c("b1", "b2", "b3", "b4", "b5", "b6", "v", "pr")]
  
  bss <- g[[1]]$z
  n.sims <- nrow(bss)
  Zsc <- vector("list", length(parms1))
  names(Zsc) <- parms1
  for (p in seq_along(parms1)) {
    ii <- grep(paste("^", parms1[p], sep = ""), names(bss))
    Zsc[[p]] <- bss[ii]
  }
  
  
  
  conv[[l]] <- 2*pnorm(-abs(unlist(Zsc))) < 0.05
  
  
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





##################################################################################
##################################################################################
#################################### Obtain results ##############################
##################################################################################
##################################################################################
Class <- 6
ClassSim <- 1
M <- 150
a <- 6.9
method <-  "JM"

betas <- matrix(, M, Class*3)
tau <- numeric(M)
gammas <- matrix(, M, Class*2)
alphas <- matrix(, M, Class)
cl_0.01 <- numeric(M)
cl_0.02 <- numeric(M)
cl_0.05 <- numeric(M)
cl_0.08 <- numeric(M)
cl_0.10 <- numeric(M)
cl_0.12 <- numeric(M)
cl_0.15 <- numeric(M)
cl_0.18 <- numeric(M)
Mweights <- matrix(, M, Class)


vecCHECK <- matrix(, M, Class)


Mod <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for (i in 1:M) {
  fileout <- paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a, "/results",i,".RData")
  load(fileout)
  
  betas[i, ] <- as.numeric(ress[[1]])
  tau[i] <- as.numeric(ress[[2]])
  gammas[i, ] <- as.numeric(ress[[3]])
  alphas[i, ] <- as.numeric(ress[[4]])
  W <- as.matrix(ress[[5]])
  N <- dim(W)[2]
  Mweights[i, ] <- apply(ress[[6]], 2, mean)
  
  EPTC_0.01 <- apply(W, 1, function(x){sum(sum(x==1)<=N/100, sum(x==2)<=N/100, sum(x==3)<=N/100, sum(x==4)<=N/100,
                                           sum(x==5)<=N/100, sum(x==6)<=N/100)})
  
  EPTC_0.02 <- apply(W, 1, function(x){sum(sum(x==1)<=2*N/100, sum(x==2)<=2*N/100, sum(x==3)<=2*N/100, sum(x==4)<=2*N/100,
                                           sum(x==5)<=2*N/100, sum(x==6)<=2*N/100)})
  
  EPTC_0.05 <- apply(W, 1, function(x){sum(sum(x==1)<=5*N/100, sum(x==2)<=5*N/100, sum(x==3)<=5*N/100, sum(x==4)<=5*N/100,
                                           sum(x==5)<=5*N/100, sum(x==6)<=5*N/100)})
  
  EPTC_0.08 <- apply(W, 1, function(x){sum(sum(x==1)<=8*N/100, sum(x==2)<=8*N/100, sum(x==3)<=8*N/100, sum(x==4)<=8*N/100,
                                           sum(x==5)<=8*N/100, sum(x==6)<=8*N/100)})
  
  EPTC_0.10 <- apply(W, 1, function(x){sum(sum(x==1)<=10*N/100, sum(x==2)<=10*N/100, sum(x==3)<=10*N/100, sum(x==4)<=10*N/100,
                                           sum(x==5)<=10*N/100, sum(x==6)<=10*N/100)})
  
  EPTC_0.12 <- apply(W, 1, function(x){sum(sum(x==1)<=12*N/100, sum(x==2)<=12*N/100, sum(x==3)<=12*N/100, sum(x==4)<=12*N/100,
                                           sum(x==5)<=12*N/100, sum(x==6)<=12*N/100)})
  
  EPTC_0.15 <- apply(W, 1, function(x){sum(sum(x==1)<=15*N/100, sum(x==2)<=15*N/100, sum(x==3)<=15*N/100, sum(x==4)<=15*N/100,
                                           sum(x==5)<=15*N/100, sum(x==6)<=15*N/100)})
  
  EPTC_0.18 <- apply(W, 1, function(x){sum(sum(x==1)<=18*N/100, sum(x==2)<=18*N/100, sum(x==3)<=18*N/100, sum(x==4)<=18*N/100,
                                           sum(x==5)<=18*N/100, sum(x==6)<=18*N/100)})
  
  
  cl_0.01[i] <- Class - Mod(EPTC_0.01)
  cl_0.02[i] <- Class - Mod(EPTC_0.02)
  cl_0.05[i] <- Class - Mod(EPTC_0.05)
  cl_0.08[i] <- Class - Mod(EPTC_0.08)
  cl_0.10[i] <- Class - Mod(EPTC_0.10)
  cl_0.12[i] <- Class - Mod(EPTC_0.12)
  cl_0.15[i] <- Class - Mod(EPTC_0.15)
  cl_0.18[i] <- Class - Mod(EPTC_0.18)
  
  par(mfrow = c(2,3))
  for (j in 1: Class) {
    plot(ress[[6]][,j], type = "l", main = paste0(i))
    #if (head(ress[[6]][,j], n = 1) - tail(ress[[6]][,j], n = 1) > 0.1) { vecCHECK[i,j] <- i}
  }
    
}


#vecNonConv <- sort(unique(vecCHECK[!is.na(vecCHECK)]))
#vec1 <- 1:81
#vecConv <- vec1[! vec1 %in% vecNonConv]

#M <- length(vecConv)

vec <- which(cl_0.10 == 1)
vecCONVER <- sort(c(vec, c(8,10,12,15,23,24,25,30,32,35,38,39,42,44,45,46,52,56,59,61,69,74,79)))

vecCon1sim3 <- c(2,7,8,9,12,14,15,16,17,18,20,26,27,28,30,31,32,33,34,35,36,38,39,40,41,42,43,44,45,46,51,54,61,
            62,63,64,65,66,67,68,69,70,71,74,78,80,81,82,85,86,87,89,90,93,94,95,96,97,98,100,105,105,107,
            108,114,117,118,119,121,125,127,129,134,137,138,140,143,145,148,149)

vecCon2sim3 <- c(1:150)[c(1:150) %!in% c(3,5,19,21,23,25,37,55,58,60,72,88,91,103,109,113,116,122,123,124,126,136,141,144)]

vecCon2sim2 <- c(1:150)[c(1:150) %!in% c(3,8,9,11,14,15,16,17,18,22,23,24,29,30,35,36,37,40,41,43,44,54,55,57,59,60,70,71,
                                         75,78,84,92,95,96,98,100,102,106,108,110,111,120,142,144,145,146,148,150)]

vecCon2sim1 <- c(1:150)[c(1:150) %!in% c(3,4,5,6,7,9,14,15,16,17,18,19,22,25,28,29,30,31,33,34,36,40,41,47,48,49,51,57,58,
                                         60,62,66,67,68,70,71,76,77,78,80,81,87,88,89,93,95,96,98,104,110,111,112,114,115,116,
                                         119,120,122,125,132,133,139,140,143,144,145,156,147,150)]

sum(cl_0.01[1:M] == ClassSim)/M
sum(cl_0.02[1:M] == ClassSim)/M
sum(cl_0.05[1:M] == ClassSim)/M
sum(cl_0.08[1:M] == ClassSim)/M
sum(cl_0.10[1:M] == ClassSim)/M
sum(cl_0.12[1:M] == ClassSim)/M
sum(cl_0.15[1:M] == ClassSim)/M

Mod(cl_0.01[1:M])
Mod(cl_0.02[1:M])
Mod(cl_0.05[1:M])
Mod(cl_0.08[1:M])
Mod(cl_0.10[1:M])
Mod(cl_0.12[1:M])
Mod(cl_0.15[1:M])


sum(cl_0.01[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.02[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.05[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.08[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.10[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.12[vecCon2] == ClassSim)/length(vecCon2)
sum(cl_0.15[vecCon2] == ClassSim)/length(vecCon2)

Mod(cl_0.01[vecCon2])
Mod(cl_0.02[vecCon2])
Mod(cl_0.05[vecCon2])
Mod(cl_0.08[vecCon2])
Mod(cl_0.10[vecCon2])
Mod(cl_0.12[vecCon2])
Mod(cl_0.15[vecCon2])

par(mfrow = c(2,3))
pie(table(cl_0.01), main = "cutoff 1%", xlab = paste("mode = ", Mod(cl_0.01)), mgp = c(0, 0, 0))
pie(table(cl_0.02), main = "cutoff 2%", xlab = paste("mode = ", Mod(cl_0.02)), mgp = c(0, 0, 0))
pie(table(cl_0.05), main = "cutoff 5%", xlab = paste("mode = ", Mod(cl_0.05)), mgp = c(0, 0, 0))
pie(table(cl_0.08), main = "cutoff 8%", xlab = paste("mode = ", Mod(cl_0.08)), mgp = c(0, 0, 0))
pie(table(cl_0.10), main = "cutoff 10%", xlab = paste("mode = ", Mod(cl_0.10)), mgp = c(0, 0, 0))
pie(table(cl_0.12), main = "cutoff 15%", xlab = paste("mode = ", Mod(cl_0.12)), mgp = c(0, 0, 0))



# scenario I
prob_a0.01 <- c(0.22, 0.25, 0.4, 0.34, 0.27, 0.18) 
prob_a0.1 <- c(0.21, 0.26, 0.27, 0.25, 0.22, 0.11) 
prob_a6.9 <- c(0.11, 0.15, 0.48, 0.45, 0.43, 0.43) 
prob_a6.99 <- c(0.06, 0.07, 0.42, 0.51, 0.56, 0.55) 

mode_a0.01 <- c(5,4,3,2,2,2) 
mode_a0.1 <- c(4,3,2,2,2,1) 
mode_a6.9 <- c(6,4,3,3,3,1) 
mode_a6.99 <- c(6,5,3,3,3,3) 

names(prob_a0.01) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a0.1) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.9) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.99) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")

par(mfrow = c(2,2))
bb <- barplot(prob_a0.01, main = "% correct classes, alpha = 0.01 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.01,")"))
bb <- barplot(prob_a0.1, main = "% correct classes, alpha = 0.1 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.1,")"))
bb <- barplot(prob_a6.9, main = "% correct classes, alpha = 6.9 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.9,")"))
bb <- barplot(prob_a6.99, main = "% correct classes, alpha = 6.99 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.99,")"))
title("ScI", line = -1, outer=TRUE)

# scenario II
prob_a0.01 <- c(0.05, 0.09, 0.22, 0.32, 0.28, 0.33) 
prob_a0.1 <- c(0.01, 0.13, 0.37, 0.42, 0.45, 0.41) 
prob_a6.9 <- c(0.12, 0.21, 0.54, 0.47, 0.47, 0.37) 
prob_a6.99 <- c(0.26, 0.32, 0.44, 0.44, 0.4, 0.31) 

mode_a0.01 <- c(4,4,3,3,3,1) 
mode_a0.1 <- c(4,3,2,2,2,1) 
mode_a6.9 <- c(3,3,2,2,1,1) 
mode_a6.99 <- c(4,2,2,1,1,1) 

names(prob_a0.01) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a0.1) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.9) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.99) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")

par(mfrow = c(2,2))
bb <- barplot(prob_a0.01, main = "% correct classes, alpha = 0.01 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.01,")"))
bb <- barplot(prob_a0.1, main = "% correct classes, alpha = 0.1 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.1,")"))
bb <- barplot(prob_a6.9, main = "% correct classes, alpha = 6.9 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.9,")"))
bb <- barplot(prob_a6.99, main = "% correct classes, alpha = 6.99 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.99,")"))
title("ScII", line = -1, outer=TRUE)

# scenario III
prob_a0.01 <- c(0.01, 0.02, 0.06, 0.12, 0.15, 0.23) 
prob_a0.1 <- c(0.05, 0.08, 0.15, 0.17, 0.18, 0.28) 
prob_a6.9 <- c(0, 0.04, 0.09, 0.15, 0.23, 0.46) 
prob_a6.99 <- c(0, 0.05, 0.13, 0.19, 0.27, 0.49) 

mode_a0.01 <- c(5,4,3,3,2,2) 
mode_a0.1 <- c(4,3,3,3,2,2) 
mode_a6.9 <- c(4,4,3,2,2,1) 
mode_a6.99 <- c(4,4,3,2,2,1) 

names(prob_a0.01) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a0.1) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.9) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")
names(prob_a6.99) <- c("cutoff 1%", "cutoff 2%", "cutoff 5%", "cutoff 8%", "cutoff 10%", "cutoff 15%")

par(mfrow = c(2,2))
bb <- barplot(prob_a0.01, main = "% correct classes, alpha = 0.01 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.01,")"))
bb <- barplot(prob_a0.1, main = "% correct classes, alpha = 0.1, \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a0.1,")"))
bb <- barplot(prob_a6.9, main = "% correct classes, alpha = 6.9 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.9,")"))
bb <- barplot(prob_a6.99, main = "% correct classes, alpha = 6.99 \n(mode correct classes)", ylab = "%", ylim = c(0,0.55))
text(bb, 0.03, paste("(",mode_a6.99,")"))
title("ScIII", line = -1, outer=TRUE)


##################################################################################
##################################################################################
#################################### Obtain results ##############################
##################################################################################
##################################################################################
Class <- 8
ClassSim <- 3
M <- 5
a <- 6.9
method <-  "JM"

betas <- matrix(, M, Class*3)
tau <- numeric(M)
gammas <- matrix(, M, Class*2)
alphas <- matrix(, M, Class)
cl_0.01 <- numeric(M)
cl_0.02 <- numeric(M)
cl_0.05 <- numeric(M)
cl_0.08 <- numeric(M)
cl_0.10 <- numeric(M)
cl_0.12 <- numeric(M)
weight <- list()

Mod <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for (i in 1:M) {
  fileout <- paste0("Results", ClassSim, "ClassSelection", Class, noquote(method), "alpha", a, "/results",i,".RData")
  load(fileout)
  
  betas[i, ] <- as.numeric(ress[[1]])
  tau[i] <- as.numeric(ress[[2]])
  gammas[i, ] <- as.numeric(ress[[3]])
  alphas[i, ] <- as.numeric(ress[[4]])
  W <- as.matrix(ress[[5]])
  N <- dim(W)[2]
  #  weight[[i]] <- res[[6]] 
  
  EPTC_0.01 <- apply(W, 1, function(x){sum(sum(x==1)<=N/100, sum(x==2)<=N/100, sum(x==3)<=N/100, sum(x==4)<=N/100,
                                           sum(x==5)<=N/100, sum(x==6)<=N/100, sum(x==7)<=N/100, sum(x==8)<=N/100)})
  
  EPTC_0.02 <- apply(W, 1, function(x){sum(sum(x==1)<=2*N/100, sum(x==2)<=2*N/100, sum(x==3)<=2*N/100, sum(x==4)<=2*N/100,
                                           sum(x==5)<=2*N/100, sum(x==6)<=2*N/100, sum(x==7)<=2*N/100, sum(x==8)<=2*N/100)})
  
  EPTC_0.05 <- apply(W, 1, function(x){sum(sum(x==1)<=5*N/100, sum(x==2)<=5*N/100, sum(x==3)<=5*N/100, sum(x==4)<=5*N/100,
                                           sum(x==5)<=5*N/100, sum(x==6)<=5*N/100, sum(x==7)<=5*N/100, sum(x==8)<=5*N/100)})
  
  EPTC_0.08 <- apply(W, 1, function(x){sum(sum(x==1)<=8*N/100, sum(x==2)<=8*N/100, sum(x==3)<=8*N/100, sum(x==4)<=8*N/100,
                                           sum(x==5)<=8*N/100, sum(x==6)<=8*N/100, sum(x==7)<=8*N/100, sum(x==8)<=8*N/100)})
  
  EPTC_0.10 <- apply(W, 1, function(x){sum(sum(x==1)<=10*N/100, sum(x==2)<=10*N/100, sum(x==3)<=10*N/100, sum(x==4)<=10*N/100,
                                           sum(x==5)<=10*N/100, sum(x==6)<=10*N/100, sum(x==7)<=10*N/100, sum(x==8)<=10*N/100)})
  
  EPTC_0.12 <- apply(W, 1, function(x){sum(sum(x==1)<=12*N/100, sum(x==2)<=12*N/100, sum(x==3)<=12*N/100, sum(x==4)<=12*N/100,
                                           sum(x==5)<=12*N/100, sum(x==6)<=12*N/100, sum(x==7)<=12*N/100, sum(x==8)<=12*N/100)})
  
  
  
  cl_0.01[i] <- Class - Mod(EPTC_0.01)
  cl_0.02[i] <- Class - Mod(EPTC_0.02)
  cl_0.05[i] <- Class - Mod(EPTC_0.05)
  cl_0.08[i] <- Class - Mod(EPTC_0.08)
  cl_0.10[i] <- Class - Mod(EPTC_0.10)
  cl_0.12[i] <- Class - Mod(EPTC_0.12)
  
  par(mfrow = c(2,3))
  for (i in 1: Class) {
    plot(ress[[6]][,i], type = "l")
  }
  
}

sum(cl_0.01 == ClassSim)/M
sum(cl_0.02 == ClassSim)/M
sum(cl_0.05 == ClassSim)/M
sum(cl_0.08 == ClassSim)/M
sum(cl_0.10 == ClassSim)/M
sum(cl_0.12 == ClassSim)/M

Mod(cl_0.01)
Mod(cl_0.02)
Mod(cl_0.05)
Mod(cl_0.08)
Mod(cl_0.10)
Mod(cl_0.12)

