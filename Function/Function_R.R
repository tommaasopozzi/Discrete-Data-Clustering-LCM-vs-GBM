# Funzioni Utili 


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


multiclass_centroid = function(x_k){
  n_k = nrow(x_k)
  d = ncol(x_k)
  #Calcolo delle P(x_{kj} = c)
  
  C = apply(x_k, 2, max)
  p = matrix(0, nrow = max(C), ncol = d)
  for(j in 1:d){
    for(c in 1:C[j]){
      p[c, j] = sum(x_k[,j] == c)/n_k
    }
  }
  return(list(matrix_prob = p, classi = C))
}

loss_multiclass = function(x_k) {
  x_tilde = multiclass_centroid(x_k)$matrix_prob
  C = multiclass_centroid(x_k)$classi
  n_k = nrow(x_k)
  d = ncol(x_k)
  loss = 0
  for (i in 1:n_k) {
    for (j in 1:d) {
      for (c in 1:C[j]) {
        if (x_k[i, j] == c) {
          # print(loss)
          # print(log(x_tilde[c, j]))
          # print(j)
          # print(c)
          loss = loss - log(x_tilde[c, j])  # Conta solo se l'elemento appartiene alla classe
        }
      }
    }
  }
  return(loss)
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
      x_tilde[, , j] = multiclass_centroid_C(x_k = X[idx_cluster, ], x = X, alpha = 0.5)$matrix_prob
      total_loss = total_loss + loss_multiclass_C(X[idx_cluster, ], x = X, alpha = 0.5)
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


kmulticlass <- function(x, k, nstart = 1, trace = FALSE) {
  
  # Integrity checks
  n <- nrow(x)
  
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  p <- ncol(x)
  # x <- matrix(x, n, p)
  
  # Random allocation with equal sizes
  
  if (k == 1) {
    # Initialize randomly with Gibbs-sampling
    G <- rep(1, n)
    freq <- as.numeric(table(G))
    best_fit <- kmulticlass_R_function(X = x, G = G, freq = freq)
    return(best_fit)
  }
  
  # Random allocation with equal sizes
  G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
  freq <- as.numeric(table(G))
  
  best_fit <- kmulticlass_R_function(X = x, G = G, freq = freq)
  if (nstart >= 2) {
    for (r in 2:nstart) {
      G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
      freq <- as.numeric(table(G))
      fit <- kmulticlass_R_function(X = x, G = G, freq = freq)
      if (fit$Loss < best_fit$Loss) best_fit <- fit
    }
  }
  return(best_fit)
}


kmultinomial_select <- function(x, k_max, nstart = 1) {
  n <- nrow(x)
  p <- ncol(x)
  x <- matrix(x, n, p)
  loss <- numeric(k_max)
  ncluster <- 1:k_max
  
  for (k in ncluster) {
    fit <- kmulticlass(x = x,  nstart = nstart, k = k, trace = FALSE)
    loss[k] <- fit$loss
  }
  p <- ggplot(data = data.frame(ncluster = ncluster, loss = loss), aes(x = ncluster, y = loss)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("Number of clusters") +
    ylab("Loss function")
  p
}


kmulticlass_gibbs <- function(x, k, lambda = 1, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE) {
  
  n <- nrow(x)
  d <- ncol(x)
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  if (trace) {
    cat("Initialization of the algorithm\n")
  }
  
  
  fit_map <- sample(1:k, n, T)
  G_map <- fit_map
  freq_map <- as.numeric(table(G_map))
  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }
  
  fit <- Gibbs_kmulticlass_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = trace)
  
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]
  
  fit$G_map <- G_map
  fit
}



#######################

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
  
  
  fit_map <- kmulticlass(x = x, k = k, trace = F, nstart = nstart)
  G_map <- fit_map$Cluster
  freq_map <- as.numeric(table(G_map))
  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }
  
  fit <- Gibbs_kmulticlass_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = T)
  
  cat("Controllo struttura fit:\n")
  print(str(fit))
  
  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  #fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]
  
  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$Loss
  fit
}


save(H, point_Binder, multiclass_centroid, loss_multiclass, kmulticlass_R_function, kmulticlass,
     kmultinomial_select, kmulticlass_gibbs, kmultinomial_gibbs, file = "funzioni.Rdata")
