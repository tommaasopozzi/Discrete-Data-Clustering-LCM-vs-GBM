devtools::load_all(".") 

library(TPthesis)

library(logging)
library(poLCA)
library(label.switching)
library(mcclust.ext)
library(mcclust)
library(Rcpp)
library(RcppArmadillo)
library(poLCA)
library(aricode)
library(cluster)
library(Rmixmod)
library(MBCbook)
library(FactoMineR)
library(cluster)
library(Rmixmod)
library(MBCbook)
library(FactoMineR)
library(proxy)
library(label.switching)
library(abind)
library(klaR)
library(future.apply)
library(coda)

setwd("exec")
Rcpp::sourceCpp("tp_thesis_cpp_function.cpp")
Rcpp::sourceCpp("tp_gbpp_thesis.cpp")
source("tp_thesis_R_function.R")


simulazione_1.1 <- list()
n_replicate <- 2

for(i in 1:n_replicate){
  simulazione_1.1[[i]] <- simulation_1_function(
    n = 100,
    tau = rep(1/3, 3),
    delta = 0.25,
    m = c(4,4,4,3,3,3)
  )
  cat("Replica simulazione", i, "\n")
}


  
simulazione_1.1 <- do.call(rbind, lapply(simulazione_1.1, as.data.frame))
save(simulazione_1.1, file = "simulazione_1.1.RData")

