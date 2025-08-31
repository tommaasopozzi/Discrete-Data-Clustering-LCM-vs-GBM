#Simulazione Prior
rm(list=ls());graphics.off()
setwd("G:/My Drive/Università/Magistrale/Secondo anno M/Bayesian Statistical Modelling/Progetto")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("G:/My Drive/Università/Magistrale/Secondo anno M/Bayesian Statistical Modelling/Progetto/funzioni_multinomiale_armadillo_finale2.cpp")
library(poLCA)
library(mcclust.ext)
library(coda)
library(mclust)
library(salso)
library(mclust)
library(ggplot2)
library(aricode)
epsilon = 0.0000000001


set.seed(1)
m = c(3,3,3,4,4,4)
d = length(m)
tau = c(0.3, 0.7)
K = length(tau)
delta = c(0.25, 0.4, 0.65)
prob = prob_poLCA(K = K, d = d, m_vec = m, delta = delta[3])
sim.data_1 <- poLCA.simdata(N = 200, probs = prob, P = tau)
data_1 = sim.data_1$dat; head(data_1)


res_prior = matrix(0, ncol = 3, nrow = 8)



#Come strutturare la simulazione? 

#Sono interessato a vedere gli effetti delle prior quando faticavo a clusterizzare 

#Primo approccio beta = (1,1)
#Secondo approccio beta = (10, 10)
#Terzo Approccio beta = (5, 15)
#Quarto approccio beta = (15, 5) #missspecificazione della prior 
#Quinto aprroccio beta = (3, 7)
#Sesto approccio beta = ceiling del metodo frequentista 

rownames(res_prior) = c("Prior Uniforme", "beta = (1,1)", "beta = (10, 10)", "beta = (5, 15)",
                  "beta = (100, 1)", "beta = (3, 7)", "Empirical Bayes", "Lcm")



beta_1 = c(1,1)
beta_2 = c(10,10)
beta_3 = c(5, 15)
beta_4 = c(1,100)
beta_5 = c(3,7)






#No prior 

start.time <- Sys.time()
gibbs_noprior <- kmultinomial_gibbs(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

G_VI_noprior <- minVI(comp.psm(gibbs_noprior$G))$cl
(ARI_VI_noprior <- adjustedRandIndex(G_VI_noprior, sim.data_1$trueclass))
G_binder_noprior <- point_Binder((gibbs_noprior$G))
ARI_noprior_Binder <- adjustedRandIndex(G_binder_noprior, sim.data_1$trueclass)

res_prior[1, 1] <- ARI_VI_noprior
res_prior[1, 2] <- ARI_noprior_Binder
res_prior[1, 3] <- time_in_sec



#Primo approccio beta = (1,1)

start.time <- Sys.time()
gibbs_1 <- kmultinomial_gibbs_prior(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, beta = c(0.3,0.7))
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


coda::traceplot(coda::as.mcmc(gibbs_1$loss))
n <- 200
h <- c()
for (i in 1:10000){
  h[i] <-  H(gibbs_1$G[i, ])
}
traceplot(coda::as.mcmc(h))
plot(coda::as.mcmc(h))

G_VI_1 <- minVI(comp.psm(gibbs_1$G))$cl
(ARI_VI_1 <- adjustedRandIndex(G_VI_1, sim.data_1$trueclass))
G_binder_1 <- point_Binder((gibbs_1$G))
ARI_1_Binder <- adjustedRandIndex(G_binder_1, sim.data_1$trueclass)

res_prior[2, 1] <- ARI_VI_1
res_prior[2, 2] <- ARI_1_Binder
res_prior[2, 3] <- time_in_sec




# Secondo approccio beta = (10, 10)

start.time <- Sys.time()
gibbs_2 <- kmultinomial_gibbs_prior(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, beta = beta_2)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

G_VI_2 <- minVI(comp.psm(gibbs_2$G))$cl
(ARI_VI_2 <- adjustedRandIndex(G_VI_2, sim.data_1$trueclass))
G_binder_2 <- point_Binder((gibbs_2$G))
ARI_2_Binder <- adjustedRandIndex(G_binder_2, sim.data_1$trueclass)

res_prior[3, 1] <- ARI_VI_2
res_prior[3, 2] <- ARI_2_Binder
res_prior[3, 3] <- time_in_sec



start.time <- Sys.time()
gibbs_3 <- kmultinomial_gibbs_prior(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, beta = beta_3)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


G_VI_3 <- minVI(comp.psm(gibbs_3$G))$cl
(ARI_VI_3 <- adjustedRandIndex(G_VI_3, sim.data_1$trueclass))
G_binder_3 <- point_Binder((gibbs_3$G))
ARI_3_Binder <- adjustedRandIndex(G_binder_3, sim.data_1$trueclass)


res_prior[4, 1] <- ARI_VI_3
res_prior[4, 2] <- ARI_3_Binder
res_prior[4, 3] <- time_in_sec

start.time <- Sys.time()
gibbs_4 <- kmultinomial_gibbs_prior(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, beta = beta_4)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


G_VI_4 <- minVI(comp.psm(gibbs_4$G))$cl
(ARI_VI_4 <- adjustedRandIndex(G_VI_4, sim.data_1$trueclass))
G_binder_4 <- point_Binder((gibbs_4$G))
ARI_4_Binder <- adjustedRandIndex(G_binder_4, sim.data_1$trueclass)

res_prior[5, 1] <- ARI_VI_4
res_prior[5, 2] <- ARI_4_Binder
res_prior[5, 3] <- time_in_sec


start.time <- Sys.time()
gibbs_5 <- kmultinomial_gibbs_prior(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, beta = beta_5)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)


G_VI_5<- minVI(comp.psm(gibbs_5$G))$cl
(ARI_VI_5 <- adjustedRandIndex(G_VI_5, sim.data_1$trueclass))
G_binder_5 <- point_Binder((gibbs_5$G))
ARI_5_Binder <- adjustedRandIndex(G_binder_5, sim.data_1$trueclass)

res_prior[6, 1] <- ARI_VI_5
res_prior[6, 2] <- ARI_5_Binder
res_prior[6, 3] <- time_in_sec


start.time <- Sys.time()
gibbs_6 <- kmultinomial_gibbs_prior_eb(as.matrix(data_1), k = 2, R = 10000, burn_in = 2000, c = 10)
end.time <- Sys.time()
time_in_sec <- as.numeric(end.time - start.time)

G_VI_6 <- minVI(comp.psm(gibbs_6$G))$cl
(ARI_VI_6 <- adjustedRandIndex(G_VI_6, sim.data_1$trueclass))
G_binder_6 <- point_Binder((gibbs_6$G))
ARI_6_Binder <- adjustedRandIndex(G_binder_6, sim.data_1$trueclass)



res_prior[7, 1] <- ARI_VI_6
res_prior[7, 2] <- ARI_6_Binder
res_prior[7, 3] <- time_in_sec











#LCM
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


res_prior[8, 1] <- ARI_1_lcm_VI
res_prior[8, 2] <- ARI_1_lcm_Binder
res_prior[8, 3] <- time_in_sec_lcm

G = sample(1:2, 200, T)

start.time_k <- Sys.time()
kmulticlass = kmulticlass_R_function(X = as.matrix(data_1), G = G, freq = as.numeric(table(G)), k = 2)
end.time_k <- Sys.time()
time_in_sec_k <- as.numeric(end.time_lcm - start.time_lcm)


ARI_kmulticlass <- adjustedRandIndex(kmulticlass$Cluster, sim.data_1$trueclass)


res_prior = rbind(res_prior, c(ARI_kmulticlass, ARI_kmulticlass, time_in_sec_k))


colnames(res_prior) = c("ARI_VI", "ARI_binder", "Time (min)")
rownames(res_prior) = c("Prior Uniforme", "beta = (1,1)", "beta = (10, 10)", "beta = (5, 15)",
                        "beta = (100, 1)", "beta = (3, 7)", "Empirical Bayes", "Lcm", "K-Bregman")



xtable(res_prior, digits = c(0,4,4,4))

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


kmultinomial_gibbs_prior <- function(x, k, lambda = 1, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE, beta) {
  
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
  
  fit <- Gibbs_kmulticlass_prior_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = T, beta)
  
  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]
  
  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$loss
  fit
}



kmultinomial_gibbs_prior_eb <- function(x, k, lambda = 1, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE, c) {
  
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
  
  beta = ceiling(freq_map/n*c)
  fit <- Gibbs_kmulticlass_prior_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = T, beta)
  
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

