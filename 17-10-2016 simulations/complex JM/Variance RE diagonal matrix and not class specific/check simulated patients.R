ClassSim = 3

vec_len_test1 <- numeric(150)
vec_len_test2 <- numeric(150)
vec_len_test3 <- numeric(150)


vec_len_test11 <- numeric(150)
vec_len_test22 <- numeric(150)
vec_len_test33 <- numeric(150)

for (l in 1:150){
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

  vec_len_test1[l] <- c(dim(data.id_1)[1])
  vec_len_test2[l] <- c(dim(data.id_2)[1])
  vec_len_test3[l] <- c(dim(data.id_3)[1])
  
  
  vec_len_test11[l] <- c(dim(data_1)[1])
  vec_len_test22[l] <- c(dim(data_2)[1])
  vec_len_test33[l] <- c(dim(data_3)[1])
}
