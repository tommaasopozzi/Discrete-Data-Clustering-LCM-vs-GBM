library(Rcpp)
library(RcppArmadillo)
sourceCpp("funzioni_multinomiale_armadillo_finale.cpp")

library(coda)
library(reshape2) # per melt()


load("alcohol_&_grades.RData")
load("Funzioni_Fond.RData") # Funzioni alla fine del file Finale.R


epsilon = 0.0000000001
data 

# set.seed(123)
p_select <- kmultinomial_select(as.matrix(data), k_max = 8, nstart = 250) 
p_select 

# Point estimate ---------------------------

n <- dim(data)[1]
k <- 4
set.seed(123)
G <- sample(1:k, n, replace = T)
fit <- kmulticlass_R_function(as.matrix(data), G, freq = as.numeric(table(G)))
G_map <- fit$Cluster
fit$Loss

# Distance Based -------------------------------------------------------

fit_gibbs <- kmultinomial_gibbs(x = as.matrix(data), k = 4, lambda = 1, R = 10000, burn_in = 5000, nstart = 100, trace = T)

traceplot(as.mcmc(fit_gibbs$loss))
h <- c()
for (i in 1:10000){
  h[i] <-  H(fit_gibbs$G[i, ])
}

#Diagnostica--------------
traceplot(as.mcmc(h))
plot(as.mcmc(h))

# Alternative point estimates -----------------

library(mcclust.ext)
G_VI <- minVI(mcclust::comp.psm(fit_gibbs$G))$cl
G_binder_1 <- point_Binder((fit_gibbs$G))


# Model Based ---------------------------------------------------------

# Da rivedere
gibbs_lcm <- gibbs_lcm_standard(Y = data, G = 4, niter = 10000, nburn = 5000, b = 4, c = 4, progress_interval = 50)

h <- c()
for (i in 1:500){
  h[i] <-  H(gibbs_lcm$cluster[i ,])
}

#Diagnostica--------------
traceplot(as.mcmc(h))
plot(as.mcmc(h))

G_VI_lcm_1 <- minVI(comp.psm(gibbs_lcm$cluster))$cl
G_binder_lcm_1 <- point_Binder((gibbs_lcm$cluster))

















