rm(list=ls());graphics.off()
setwd("G:/My Drive/Università/Magistrale/Secondo anno M/Bayesian Statistical Modelling/Progetto")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("G:/My Drive/Università/Magistrale/Secondo anno M/Bayesian Statistical Modelling/Progetto/funzioni_multinomiale_armadillo_finale.cpp")
library(poLCA)
library(mcclust.ext)
library(coda)
library(mclust)
library(salso)
library(mclust)
library(ggplot2)
library(aricode)
epsilon = 0.0000000001


# Simulazione 2 -----------------------------------------------------------

#Considero k = 4 con proporzioni simili e pochi livelli per variabile. 
#Effettuo 4 diverse simulazioni con 3 livelli di distanza dei cluster
#Delta = 0.2, 0.5 e 0.7 ed N = 200, 400, 600, 800
#Ricordarsi di salvare i tempi 

#Variabili con livelli simili 2 cluster molto distanti tra loro 

res = matrix(0, nrow = 12, ncol = 5)
res_lcm =  matrix(0, nrow = 12, ncol = 5)

set.seed(1)
m = c(3,15,10,3,13,4)
d = length(m)
tau = c(0.2, 0.5, 0.3)
K = length(tau)
delta = c(0.25, 0.4, 0.6)
prob = prob_poLCA(K = K, d = d, m_vec = m, delta = delta[1])
sim.data_1 <- poLCA.simdata(N = 200, probs = prob, P = tau)
data_1 = sim.data_1$dat; head(data_1)

#Gibbs-Sampling-----------
start.time <- Sys.time()
gibbs_1 = kmultinomial_gibbs(x = as.matrix(data_1), k = K, lambda = 1, R = 10000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_1$loss))
n <- 200
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_1$G[i, ])
}

#Diagnostica--------------
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

#Clustering---------------

G_VI_1 <- minVI(comp.psm(gibbs_1$G))$cl
ARI_1_VI <- adjustedRandIndex(G_VI_1, sim.data_1$trueclass)
G_binder_1 <- point_Binder((gibbs_1$G))
ARI_1_Binder <- adjustedRandIndex(G_binder_1, sim.data_1$trueclass)
NMI(sim.data_1$trueclass, G_VI_1)

#LCM----------------------
start.time_lcm <- Sys.time()
gibbs_lcm_1 <- gibbs_lcm_standard(Y = data_1, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 200
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_1$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_1 <- minVI(comp.psm(gibbs_lcm_1$cluster))$cl
ARI_1_lcm_VI <- adjustedRandIndex(G_VI_lcm_1, sim.data_1$trueclass)
G_binder_lcm_1 <- point_Binder((gibbs_lcm_1$cluster))
ARI_1_lcm_Binder <- adjustedRandIndex(G_binder_lcm_1, sim.data_1$trueclass)


#Salvo i risultati--------

res[1, 1] <- 200
res[1, 2] <- ARI_1_VI
res[1, 3] <- ARI_1_Binder
res[1, 4] <- delta[1]
res[1, 5] <- time_in_sec

res_lcm[1, 1] <- 200
res_lcm[1, 2] <- ARI_1_lcm_VI
res_lcm[1, 3] <- ARI_1_lcm_Binder
res_lcm[1, 4] <- delta[1]
res_lcm[1, 5] <- time_in_sec_lcm


#Aumento n---------------

sim.data_2 <- poLCA.simdata(N = 400, probs = prob, P = tau)
data_2 = sim.data_2$dat; head(data_2)

start.time <- Sys.time()
gibbs_2 = kmultinomial_gibbs(x = as.matrix(data_2), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)



coda::traceplot(coda::as.mcmc(gibbs_2$loss))
n <- 400
h <- c()
for (i in 1:12000){
  h[i] <-  H(gibbs_2$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_2 <- minVI(comp.psm(gibbs_2$G))$cl
ARI_2_VI <- adjustedRandIndex(G_VI_2, sim.data_2$trueclass)
G_binder_2 <- point_Binder((gibbs_2$G))
ARI_2_Binder <- adjustedRandIndex(G_binder_2, sim.data_2$trueclass)
res[2, 1] <- 400
res[2, 2] <- ARI_2_VI
res[2, 3] <- ARI_2_Binder
res[2, 4] <- delta[1]
res[2, 5] <- time_in_sec




#LCM----------------------

start.time_lcm <- Sys.time()
gibbs_lcm_2 <- gibbs_lcm_standard(Y = data_2, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)


n <- 400
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_2$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_2 <- minVI(comp.psm(gibbs_lcm_2$cluster))$cl
ARI_2_lcm_VI <- adjustedRandIndex(G_VI_lcm_2, sim.data_2$trueclass)
G_binder_lcm_2 <- point_Binder((gibbs_lcm_2$cluster))
ARI_2_lcm_Binder <- adjustedRandIndex(G_binder_lcm_2, sim.data_2$trueclass)

res_lcm[2, 1] <- 400
res_lcm[2, 2] <- ARI_2_lcm_VI
res_lcm[2, 3] <- ARI_2_lcm_Binder
res_lcm[2, 4] <- delta[1]
res_lcm[2, 5] <- time_in_sec_lcm



#N = 600----------------
sim.data_3 <- poLCA.simdata(N = 600, probs = prob, P = tau)
data_3 = sim.data_3$dat; head(data_3)

start.time <- Sys.time()
gibbs_3 = kmultinomial_gibbs(x = as.matrix(data_3), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

coda::traceplot(coda::as.mcmc(gibbs_3$loss))
n <- 600
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_3$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_3 <- minVI(comp.psm(gibbs_3$G))$cl
ARI_3_VI <- adjustedRandIndex(G_VI_3, sim.data_3$trueclass)
G_binder_3 <- point_Binder((gibbs_3$G))
ARI_3_Binder <- adjustedRandIndex(G_binder_3, sim.data_3$trueclass)
res[3, 1] <- 600
res[3, 2] <- ARI_3_VI
res[3, 3] <- ARI_3_Binder
res[3, 4] <- delta[1]
res[3, 5] <- time_in_sec

start.time_lcm <- Sys.time()
gibbs_lcm_3 <- gibbs_lcm_standard(Y = data_3, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 600
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_3$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_3 <- minVI(comp.psm(gibbs_lcm_3$cluster))$cl
ARI_3_lcm_VI <- adjustedRandIndex(G_VI_lcm_3, sim.data_3$trueclass)
G_binder_lcm_3 <- point_Binder((gibbs_lcm_3$cluster))
ARI_3_lcm_Binder <- adjustedRandIndex(G_binder_lcm_3, sim.data_3$trueclass)

res_lcm[3, 1] <- 600
res_lcm[3, 2] <- ARI_3_lcm_VI
res_lcm[3, 3] <- ARI_3_lcm_Binder
res_lcm[3, 4] <- delta[1]
res_lcm[3, 5] <- time_in_sec_lcm


#N = 800----------------


sim.data_4 <- poLCA.simdata(N = 800, probs = prob, P = tau)
data_4 = sim.data_4$dat; head(data_4)

start.time <- Sys.time()
gibbs_4 = kmultinomial_gibbs(x = as.matrix(data_4), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_4$loss))
n <- 800
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_4$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_4 <- minVI(comp.psm(gibbs_4$G))$cl
ARI_4_VI <- adjustedRandIndex(G_VI_4, sim.data_4$trueclass)
G_binder_4 <- point_Binder((gibbs_4$G))
ARI_4_Binder <- adjustedRandIndex(G_binder_4, sim.data_4$trueclass)

res[4, 1] <- 800
res[4, 2] <- ARI_4_VI
res[4, 3] <- ARI_4_Binder
res[4, 4] <- delta[1]
res[4, 5] <- time_in_sec


start.time_lcm <- Sys.time()
gibbs_lcm_4 <- gibbs_lcm_standard(Y = data_4, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)


n <- 800
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_4$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_4 <- minVI(comp.psm(gibbs_lcm_4$cluster))$cl
ARI_4_lcm_VI <- adjustedRandIndex(G_VI_lcm_4, sim.data_4$trueclass)
G_binder_lcm_4 <- point_Binder((gibbs_lcm_4$cluster))
ARI_4_lcm_Binder <- adjustedRandIndex(G_binder_lcm_4, sim.data_4$trueclass)


res_lcm[4, 1] <- 800
res_lcm[4, 2] <- ARI_4_lcm_VI
res_lcm[4, 3] <- ARI_4_lcm_Binder
res_lcm[4, 4] <- delta[1]
res_lcm[4, 5] <- time_in_sec_lcm



# Delta[2] ----------------------------------------------------------------

prob = prob_poLCA(K = K, d = d, m_vec = m, delta = delta[2])
sim.data_1 <- poLCA.simdata(N = 200, probs = prob, P = tau)
data_1 = sim.data_1$dat; head(data_1)

#Gibbs-Sampling-----------

start.time <- Sys.time()
gibbs_1 = kmultinomial_gibbs(x = as.matrix(data_1), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_1$loss))
n <- 200
h <- c()
for (i in 1:12000){
  h[i] <-  H(gibbs_1$G[i, ])
}

#Diagnostica--------------
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

#Clustering---------------

G_VI_1 <- minVI(comp.psm(gibbs_1$G))$cl
ARI_1_VI <- adjustedRandIndex(G_VI_1, sim.data_1$trueclass)
G_binder_1 <- point_Binder((gibbs_1$G))
ARI_1_Binder <- adjustedRandIndex(G_binder_1, sim.data_1$trueclass)
NMI(sim.data_1$trueclass, G_VI_1)

#LCM----------------------
start.time_lcm <- Sys.time()
gibbs_lcm_1 <- gibbs_lcm_standard(Y = data_1, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 200
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_1$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_1 <- minVI(comp.psm(gibbs_lcm_1$cluster))$cl
ARI_1_lcm_VI <- adjustedRandIndex(G_VI_lcm_1, sim.data_1$trueclass)
G_binder_lcm_1 <- point_Binder((gibbs_lcm_1$cluster))
ARI_1_lcm_Binder <- adjustedRandIndex(G_binder_lcm_1, sim.data_1$trueclass)


#Salvo i risultati--------

res[5, 1] <- 200
res[5, 2] <- ARI_1_VI
res[5, 3] <- ARI_1_Binder
res[5, 4] <- delta[2]
res[5, 5] <- time_in_sec

res_lcm[5, 1] <- 200
res_lcm[5, 2] <- ARI_1_lcm_VI
res_lcm[5, 3] <- ARI_1_lcm_Binder
res_lcm[5, 4] <- delta[2]
res_lcm[5, 5] <- time_in_sec_lcm

#Aumento n---------------

sim.data_2 <- poLCA.simdata(N = 400, probs = prob, P = tau)
data_2 = sim.data_2$dat; head(data_2)

start.time <- Sys.time()
gibbs_2 = kmultinomial_gibbs(x = as.matrix(data_2), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_2$loss))
n <- 400
h <- c()
for (i in 1:12000){
  h[i] <-  H(gibbs_2$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_2 <- minVI(comp.psm(gibbs_2$G))$cl
ARI_2_VI <- adjustedRandIndex(G_VI_2, sim.data_2$trueclass)
G_binder_2 <- point_Binder((gibbs_2$G))
ARI_2_Binder <- adjustedRandIndex(G_binder_2, sim.data_2$trueclass)
res[6, 1] <- 400
res[6, 2] <- ARI_2_VI
res[6, 3] <- ARI_2_Binder
res[6, 4] <- delta[2]
res[6, 5] <- time_in_sec



#LCM----------------------

start.time_lcm <- Sys.time()
gibbs_lcm_2 <- gibbs_lcm_standard(Y = data_2, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 400
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_2$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_2 <- minVI(comp.psm(gibbs_lcm_2$cluster))$cl
ARI_2_lcm_VI <- adjustedRandIndex(G_VI_lcm_2, sim.data_2$trueclass)
G_binder_lcm_2 <- point_Binder((gibbs_lcm_2$cluster))
ARI_2_lcm_Binder <- adjustedRandIndex(G_binder_lcm_2, sim.data_2$trueclass)

res_lcm[6, 1] <- 400
res_lcm[6, 2] <- ARI_2_lcm_VI
res_lcm[6, 3] <- ARI_2_lcm_Binder
res_lcm[6, 4] <- delta[2]
res_lcm[6, 5] <- time_in_sec_lcm



#N = 600----------------
sim.data_3 <- poLCA.simdata(N = 600, probs = prob, P = tau)
data_3 = sim.data_3$dat; head(data_3)

start.time <- Sys.time()
gibbs_3 = kmultinomial_gibbs(x = as.matrix(data_3), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_3$loss))
n <- 600
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_3$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_3 <- minVI(comp.psm(gibbs_3$G))$cl
ARI_3_VI <- adjustedRandIndex(G_VI_3, sim.data_3$trueclass)
G_binder_3 <- point_Binder((gibbs_3$G))
ARI_3_Binder <- adjustedRandIndex(G_binder_3, sim.data_3$trueclass)
res[7, 1] <- 600
res[7, 2] <- ARI_3_VI
res[7, 3] <- ARI_3_Binder
res[7, 4] <- delta[2]
res[7, 5] <- time_in_sec

start.time_lcm <- Sys.time()
gibbs_lcm_3 <- gibbs_lcm_standard(Y = data_3, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 600
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_3$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_3 <- minVI(comp.psm(gibbs_lcm_3$cluster))$cl
ARI_3_lcm_VI <- adjustedRandIndex(G_VI_lcm_3, sim.data_3$trueclass)
G_binder_lcm_3 <- point_Binder((gibbs_lcm_3$cluster))
ARI_3_lcm_Binder <- adjustedRandIndex(G_binder_lcm_3, sim.data_3$trueclass)

res_lcm[7, 1] <- 600
res_lcm[7, 2] <- ARI_3_lcm_VI
res_lcm[7, 3] <- ARI_3_lcm_Binder
res_lcm[7, 4] <- delta[2]
res_lcm[7, 5] <- time_in_sec_lcm


#N = 800----------------
sim.data_4 <- poLCA.simdata(N = 800, probs = prob, P = tau)
data_4 = sim.data_4$dat; head(data_4)

start.time <- Sys.time()
gibbs_4 = kmultinomial_gibbs(x = as.matrix(data_4), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
coda::traceplot(coda::as.mcmc(gibbs_4$loss))
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

n <- 800
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_4$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_4 <- minVI(comp.psm(gibbs_4$G))$cl
ARI_4_VI <- adjustedRandIndex(G_VI_4, sim.data_4$trueclass)
G_binder_4 <- point_Binder((gibbs_4$G))
ARI_4_Binder <- adjustedRandIndex(G_binder_4, sim.data_4$trueclass)
res[8, 1] <- 800
res[8, 2] <- ARI_4_VI
res[8, 3] <- ARI_4_Binder
res[8, 4] <- delta[2]
res[8, 5] <- time_in_sec

start.time_lcm <- Sys.time()
gibbs_lcm_4 <- gibbs_lcm_standard(Y = data_4, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)


n <- 800
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_4$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_4 <- minVI(comp.psm(gibbs_lcm_4$cluster))$cl
ARI_4_lcm_VI <- adjustedRandIndex(G_VI_lcm_4, sim.data_4$trueclass)
G_binder_lcm_4 <- point_Binder((gibbs_lcm_4$cluster))
ARI_4_lcm_Binder <- adjustedRandIndex(G_binder_lcm_4, sim.data_4$trueclass)
library(aricode)
NMI(c1 = sim.data_4$trueclass, c2 = G_VI_lcm_4)
# 0.9169864


res_lcm[8, 1] <- 800
res_lcm[8, 2] <- ARI_4_lcm_VI
res_lcm[8, 3] <- ARI_4_lcm_Binder
res_lcm[8, 4] <- delta[2]
res_lcm[8, 5] <- time_in_sec_lcm


# Delta[3] ----------------------------------------------------------------

prob = prob_poLCA(K = K, d = d, m_vec = m, delta = delta[3])
sim.data_1 <- poLCA.simdata(N = 200, probs = prob, P = tau)
data_1 = sim.data_1$dat; head(data_1)

#Gibbs-Sampling-----------

start.time <- Sys.time()
gibbs_1 = kmultinomial_gibbs(x = as.matrix(data_1), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_1$loss))
n <- 200
h <- c()
for (i in 1:12000){
  h[i] <-  H(gibbs_1$G[i, ])
}

#Diagnostica--------------
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

#Clustering---------------

G_VI_1 <- minVI(comp.psm(gibbs_1$G))$cl
ARI_1_VI <- adjustedRandIndex(G_VI_1, sim.data_1$trueclass)
G_binder_1 <- point_Binder((gibbs_1$G))
ARI_1_Binder <- adjustedRandIndex(G_binder_1, sim.data_1$trueclass)
NMI(sim.data_1$trueclass, G_VI_1)

#LCM----------------------
start.time_lcm <- Sys.time()
gibbs_lcm_1 <- gibbs_lcm_standard(Y = data_1, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)


gibbs_lcm_2 <- gibbs_lcm_standard_cpp(Y = as.matrix(data_1), G = K, niter = 12000, nburn = 2000, b = 4, c = 4)

n <- 200
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_1$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_1 <- minVI(comp.psm(gibbs_lcm_1$cluster))$cl
ARI_1_lcm_VI <- adjustedRandIndex(G_VI_lcm_1, sim.data_1$trueclass)
G_binder_lcm_1 <- point_Binder((gibbs_lcm_1$cluster))
ARI_1_lcm_Binder <- adjustedRandIndex(G_binder_lcm_1, sim.data_1$trueclass)


#Salvo i risultati--------

res[9, 1] <- 200
res[9, 2] <- ARI_1_VI
res[9, 3] <- ARI_1_Binder
res[9, 4] <- delta[3]
res[9, 5] <- time_in_sec

res_lcm[9, 1] <- 200
res_lcm[9, 2] <- ARI_1_lcm_VI
res_lcm[9, 3] <- ARI_1_lcm_Binder
res_lcm[9, 4] <- delta[3]
res_lcm[9, 5] <- time_in_sec_lcm


#Aumento n---------------

sim.data_2 <- poLCA.simdata(N = 400, probs = prob, P = tau)
data_2 = sim.data_2$dat; head(data_2)

start.time <- Sys.time()
gibbs_2 = kmultinomial_gibbs(x = as.matrix(data_2), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_2$loss))


n <- 400
h <- c()
for (i in 1:12000){
  h[i] <-  H(gibbs_2$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_2 <- minVI(comp.psm(gibbs_2$G))$cl
ARI_2_VI <- adjustedRandIndex(G_VI_2, sim.data_2$trueclass)
G_binder_2 <- point_Binder((gibbs_2$G))
ARI_2_Binder <- adjustedRandIndex(G_binder_2, sim.data_2$trueclass)

res[10, 1] <- 400
res[10, 2] <- ARI_2_VI
res[10, 3] <- ARI_2_Binder
res[10, 4] <- delta[3]
res[10, 5] <- time_in_sec



#LCM----------------------

start.time_lcm <- Sys.time()
gibbs_lcm_2 <- gibbs_lcm_standard(Y = data_2, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)


n <- 400
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_2$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_2 <- minVI(comp.psm(gibbs_lcm_2$cluster))$cl
ARI_2_lcm_VI <- adjustedRandIndex(G_VI_lcm_2, sim.data_2$trueclass)
G_binder_lcm_2 <- point_Binder((gibbs_lcm_2$cluster))
ARI_2_lcm_Binder <- adjustedRandIndex(G_binder_lcm_2, sim.data_2$trueclass)

res_lcm[10, 1] <- 400
res_lcm[10, 2] <- ARI_2_lcm_VI
res_lcm[10, 3] <- ARI_2_lcm_Binder
res_lcm[10, 4] <- delta[3]
res_lcm[10, 5] <- time_in_sec_lcm



#N = 600----------------
sim.data_3 <- poLCA.simdata(N = 600, probs = prob, P = tau)
data_3 = sim.data_3$dat; head(data_3)

start.time <- Sys.time()
gibbs_3 = kmultinomial_gibbs(x = as.matrix(data_3), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_3$loss))
n <- 600
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_3$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_3 <- minVI(comp.psm(gibbs_3$G))$cl
ARI_3_VI <- adjustedRandIndex(G_VI_3, sim.data_3$trueclass)
G_binder_3 <- point_Binder((gibbs_3$G))
ARI_3_Binder <- adjustedRandIndex(G_binder_3, sim.data_3$trueclass)
res[11, 1] <- 600
res[11, 2] <- ARI_3_VI
res[11, 3] <- ARI_3_Binder
res[11, 4] <- delta[3]
res[11, 5] <- time_in_sec

start.time_lcm <- Sys.time()
gibbs_lcm_3 <- gibbs_lcm_standard(Y = data_3, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 600
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_3$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_3 <- minVI(comp.psm(gibbs_lcm_3$cluster))$cl
ARI_3_lcm_VI <- adjustedRandIndex(G_VI_lcm_3, sim.data_3$trueclass)
G_binder_lcm_3 <- point_Binder((gibbs_lcm_3$cluster))
ARI_3_lcm_Binder <- adjustedRandIndex(G_binder_lcm_3, sim.data_3$trueclass)

res_lcm[11, 1] <- 600
res_lcm[11, 2] <- ARI_3_lcm_VI
res_lcm[11, 3] <- ARI_3_lcm_Binder
res_lcm[11, 4] <- delta[3]
res_lcm[11, 5] <- time_in_sec_lcm


#N = 800----------------
sim.data_4 <- poLCA.simdata(N = 800, probs = prob, P = tau)
data_4 = sim.data_4$dat; head(data_4)

start.time <- Sys.time()
gibbs_4 = kmultinomial_gibbs(x = as.matrix(data_4), k = K, lambda = 1, R = 12000, burn_in = 2000, nstart = 100, trace = T)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

coda::traceplot(coda::as.mcmc(gibbs_4$loss))
n <- 800
h <- c()
for (i in 1:8000){
  h[i] <-  H(gibbs_4$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))
G_VI_4 <- minVI(comp.psm(gibbs_4$G))$cl
ARI_4_VI <- adjustedRandIndex(G_VI_4, sim.data_4$trueclass)
G_binder_4 <- point_Binder((gibbs_4$G))
ARI_4_Binder <- adjustedRandIndex(G_binder_4, sim.data_4$trueclass)
res[12, 1] <- 800
res[12, 2] <- ARI_4_VI
res[12, 3] <- ARI_4_Binder
res[12, 4] <- delta[3]
res[12, 5] <- time_in_sec

start.time_lcm <- Sys.time()
gibbs_lcm_4 <- gibbs_lcm_standard(Y = data_4, G = K, niter = 12000, nburn = 2000, b = 4, c = 4, progress_interval = 50)
end.time_lcm <- Sys.time()
time_in_sec_lcm <- as.numeric(end.time_lcm - start.time_lcm)

n <- 800
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_lcm_4$cluster[i ,])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_lcm_4 <- minVI(comp.psm(gibbs_lcm_4$cluster))$cl
ARI_4_lcm_VI <- adjustedRandIndex(G_VI_lcm_4, sim.data_4$trueclass)
G_binder_lcm_4 <- point_Binder((gibbs_lcm_4$cluster))
ARI_4_lcm_Binder <- adjustedRandIndex(G_binder_lcm_4, sim.data_4$trueclass)
library(aricode)
NMI(c1 = sim.data_4$trueclass, c2 = G_VI_lcm_4)
# 0.9169864


res_lcm[12, 1] <- 800
res_lcm[12, 2] <- ARI_4_lcm_VI
res_lcm[12, 3] <- ARI_4_lcm_Binder
res_lcm[12, 4] <- delta[3]
res_lcm[12, 5] <- time_in_sec_lcm

save(res, file = "simulazione2_bpm.RData")
save(res_lcm, file = "simulazione2_lcm.RData")


######################Funzioni R########################################


alpha_k_jh = function(k, h, m_j, delta){
  h_check <- ((k-1) %% m_j) + 1
  if(h_check == h){
    return(1 / m_j + (1 - delta)*(m_j-1)/m_j)
  }else{
    return((1-1/m_j-(1-delta)*(m_j-1)/m_j)/(m_j-1))
  }
}

prob <- function(K, d, m_vec, delta = 0.1){
  prob_list = list()
  for(k in 1:K){
    prob_matrix = matrix(0, nrow = max(m_vec), ncol = d)
    for(j in 1:d){
      for(h in 1:m[j]){
        prob_matrix[h, j] <- alpha_k_jh(k = k, h = h, m_j = m[j], delta = delta)
      }
    }
    prob_list[[k]] <- prob_matrix
  }
  return(prob_list)
}


prob_poLCA <- function(K, d, m_vec, delta = 0.1) {
  prob_list <- vector("list", d)
  
  for (j in 1:d) {
    m_j <- m_vec[j]
    mat <- matrix(0, nrow = K, ncol = m_j)
    for (k in 1:K) {
      for (h in 1:m_j) {
        mat[k, h] <- alpha_k_jh(k = k, h = h, m_j = m_j, delta = delta)
      }
    }
    prob_list[[j]] <- mat
  }
  
  return(prob_list)
}




#Kmultinomial_Gibbs
kmultinomial_gibbs <- function(x, k, lambda = 1, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE) {
  
  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)
  
  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  if (trace) {
    cat("Initialization of the algorithm\n")
  }
  
  G <- sample(1:k, n, replace = T)
  fit_map <- kmulticlass_R_function(X = x, G = G, freq = as.numeric(table(G)), k = k)
  G_map <- fit_map$Cluster
  freq_map <- as.numeric(table(G_map))
  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }
  
  fit <- Gibbs_kmulticlass_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = T)
  
  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]
  
  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$loss
  fit
}


kmulticlass_R_function <- function(X, G, freq, k = 3) {
  n = nrow(X)
  d = ncol(X)
  
  losses = numeric(k)
  convergence = F
  
  total_loss_old = -1
  C = max(apply(X, 2, max))
  x_tilde =  array(0, dim = c(C, d, k))
  
  total_loss = 0
  for (j in 1:k) {
    idx_cluster = which(G == j)
    x_tilde[, , j] = multiclass_centroid_C(x_k = X[idx_cluster, ], x = X, alpha = 0.5)$matrix_prob
    total_loss = total_loss + loss_multiclass_C(X[idx_cluster, ], x = X, alpha = 0.5)
    # cat("Loss_function : " total_loss, "\n")
  }
  
  while (!convergence) {
    for (i in 1:n) {
      for (g in 1:k) {
        prob = numeric(d)
        for (j in 1:d) {
          c <- as.integer(X[i, j])
          prob[j] <- max(x_tilde[c, j, g], epsilon)
        }
        losses[g] <- -sum(log(prob))
      }
      
      freq[G[i]] = freq[G[i]] - 1
      
      if (freq[G[i]] > 0) {
        G[i] = which.min(losses)
      }
      freq[G[i]] = freq[G[i]] + 1
    }
    total_loss = 0
    x_tilde =  array(0, dim = c(C, d, k))
    for (j in 1:k) {
      idx_cluster = which(G == j)
      x_tilde[, , j] = multiclass_centroid_C(x_k = X[idx_cluster, ], x = X, alpha = 0.3)$matrix_prob
      total_loss = total_loss + loss_multiclass_C(X[idx_cluster, ], x = X, alpha = 0.3)
    }
    if ((total_loss_old - total_loss) == 0) {
      convergence = T
    } else{
      total_loss_old = total_loss
      cat("Total_loss", total_loss, "\n")
    }
  }
  return(list("Cluster" = G,
              "Centers" = x_tilde,
              "Loss" = total_loss))
}


kmultinomial_select <- function(x, k_max, nstart = 1) {
  n <- nrow(x)
  p <- ncol(x)
  x <- matrix(x, n, p)
  loss <- numeric(k_max)
  ncluster <- 1:k_max
  
  for (k in ncluster) {
    G <- sample(1:k, n, replace = T)
    fit <- kmulticlass_R_function(X = x, G = G, k = k, freq = as.numeric(table(G)))
    loss[k] <- fit$Loss
  }
  p <- ggplot(data = data.frame(ncluster = ncluster, loss = loss), aes(x = ncluster, y = loss)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("Number of clusters") +
    ylab("Loss function")+
    scale_x_continuous(breaks = 1:k_max)
  p
}

point_Binder <- function(partition_matrix){
  
  n <- ncol(partition_matrix)
  PSM <- matrix(0, ncol = n, nrow = n)
  dist_Binder <- c()
  
  for(i in 1:n){
    for(j in 1:i){
      PSM[i,j] <- PSM[j,i] <- mean(partition_matrix[,i] == partition_matrix[,j])
    }
  }
  for(i in 1:nrow(partition_matrix)){
    temp_mat <- matrix(as.numeric(sapply(partition_matrix[i,], 
                                         function(z) z == partition_matrix[i,])), ncol = n)
    dist_Binder[i] <- sum((temp_mat - PSM)^2)
  }
  point_estimate <- partition_matrix[which.min(dist_Binder),]
  
  return(point_estimate)
}

H <- function(partizione){
  N <- n
  kn <- max(partizione)
  nj <- c()
  for (j in 1:kn){
    nj[j] <- sum(partizione == j)
  }
  somma <- sum(nj*log(nj))
  log(N, base = 2) - 1/N*somma
}



#LCM function 

gibbs_lcm_standard <- function(Y, G, niter, nburn, b, c, progress_interval = 50) {
  
  start_time <- Sys.time()
  
  n <- nrow(Y)
  d <- ncol(Y)
  mj <- apply(Y, 2, max)
  cat(sprintf("Dataset: %d observation, %d variable\n", n, d))
  out_cluster <- matrix(0, nrow = niter, ncol = n)
  
  alpha <- array(0, dim = c(d, max(mj), G)) 
  for(g in 1:G) {
    for(j in 1:d) {
      alpha[j, 1:mj[j], g] <- rdirichlet(1, rep(1, mj[j])) #Inizio simulando con valori casuali dalla dirichlet
    }
  }
  
  tau <- rdirichlet(1, rep(1, G)) #Inizializzo con valori casuali dalla dirichlet 
  z <- sample(1:G, size = n, replace = TRUE, prob = tau) #Inizializzo le allocazioni in maniera casuale 
  cat("Start Gibbs sampling...\n")
  cat(sprintf("Parameters: G=%d gruppi, niter=%d iterations, hyperparameters b=%g, c=%g\n", 
              G, niter, b, c))
  pb <- txtProgressBar(min = 1, max = niter, style = 3)
  
  for(m in 1:niter) {
    
    # Mostra info aggiuntive a intervalli
    if (m %% progress_interval == 0 || niter == 1) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      estimated_total <- elapsed * (niter / m)
      remaining <- estimated_total - elapsed
      
      cat("\nIterazione", m, "di", niter, 
          sprintf("(%.1f%%)", m/niter*100),
          sprintf("- Tempo trascorso: %.2f min, Stimato rimanente: %.2f min\n", 
                  as.numeric(elapsed), as.numeric(remaining)))
    }
    
    
    for(i in 1:n) {
      temp_prob <- rep(0, G)
      
      for(g in 1:G) {
        log_prob <- log(tau[g])
        for(j in 1:d) {
          if(!is.na(Y[i,j])) {
            log_prob <- log_prob + log(alpha[j, Y[i,j], g])
          }
        }
        temp_prob[g] <- exp(log_prob)
      }
      
      temp_prob <- temp_prob / sum(temp_prob)
      z[i] <- sample(1:G, size = 1, prob = temp_prob)
    }
    for(g in 1:G) {
      for(j in 1:d) {
        counts <- rep(0, mj[j])
        for(h in 1:mj[j]) {
          counts[h] <- sum(Y[z == g, j] == h, na.rm = TRUE)
        }
        alpha[j, 1:mj[j], g] <- rdirichlet(1, c + counts)
      }
    }
    
    n_g <- tabulate(z, nbins = G)
    tau <- rdirichlet(1, b + n_g)
    
    out_cluster[m, ] <- z
    
    setTxtProgressBar(pb, m)
  }
  close(pb)
  
  return(list(
    cluster = out_cluster[-c(1:nburn),],
    alpha = alpha,
    tau = tau
  ))
}

# Funzione di supporto per il campionamento dalla distribuzione Dirichlet
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  result <- matrix(0, nrow = n, ncol = k)
  
  for(i in 1:n) {
    gamma_samples <- rgamma(k, shape = alpha, scale = 1)
    result[i,] <- gamma_samples / sum(gamma_samples)
  }
  
  if(n == 1) {
    return(as.vector(result))
  } else {
    return(result)
  }
}

