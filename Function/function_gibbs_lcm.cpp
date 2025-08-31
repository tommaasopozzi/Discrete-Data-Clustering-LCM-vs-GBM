#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Funzione per campionare dalla distribuzione Dirichlet
vec rdirichlet(const vec& alpha) {
  vec sample(alpha.n_elem);
  for (unsigned int i = 0; i < alpha.n_elem; ++i) {
    sample(i) = R::rgamma(alpha(i), 1.0);
  }
  return sample / sum(sample);
} 

// Campionamento discreto da vettore di probabilità
int sample_index(const vec& prob) {
  double u = R::runif(0.0, 1.0);
  double cum = 0.0;
  for (unsigned int i = 0; i < prob.n_elem; ++i) {
    cum += prob(i);
    if (u <= cum) return i;
  } 
  return prob.n_elem - 1; // fallback
} 

// [[Rcpp::export]]
List gibbs_lcm_standard_cpp(const arma::mat& Y, int G, int niter, int nburn,
                            double b, double c, int progress_interval = 50) {
  
  int n = Y.n_rows;
  int d = Y.n_cols;
  ivec mj(d);
  
  for (int j = 0; j < d; ++j) {
    mj(j) = Y.col(j).max(); // valori interi 1, ..., mj
  } 
  
  Rcout << "Dataset: " << n << " observation, " << d << " variable\n";
  
  imat out_cluster(niter, n, fill::zeros);
  int mj_max = mj.max();
  cube alpha(d, mj_max, G, fill::zeros);
  
  // Inizializzazione di alpha
  for (int g = 0; g < G; ++g) {
    for (int j = 0; j < d; ++j) {
      vec temp = rdirichlet(vec(mj(j), fill::ones));
      for (int h = 0; h < mj(j); ++h) {
        alpha(j, h, g) = temp(h);
      }
    } 
  }
  
  // Inizializzazione di tau
  vec tau = rdirichlet(vec(G, fill::ones));
  
  // Inizializzazione delle allocazioni z
  ivec z(n);
  for (int i = 0; i < n; ++i) {
    z(i) = sample_index(tau) + 1;
  } 
  
  Rcout << "Start Gibbs sampling...\n";
  Rcout << "Parameters: G=" << G << " gruppi, niter=" << niter << " iterations, hyperparameters b=" << b << ", c=" << c << "\n";
  
  for (int m = 0; m < niter; ++m) {
    if ((m + 1) % progress_interval == 0 || niter == 1) {
      Rcout << "\nIterazione " << (m + 1) << " di " << niter << " (" << ((double)(m + 1) / niter * 100) << "%)\n";
    } 
    
    // Aggiorna z
    for (int i = 0; i < n; ++i) {
      vec log_prob(G, fill::zeros);
      for (int g = 0; g < G; ++g) {
        log_prob(g) = std::log(tau(g));
        for (int j = 0; j < d; ++j) {
          double yij = Y(i, j);
          if (!std::isnan(yij)) {
            int val = static_cast<int>(yij) - 1;
            if (val >= 0 && val < mj(j)) {
              log_prob(g) += std::log(alpha(j, val, g));
            }
          }
        } 
      }
      vec prob = exp(log_prob - max(log_prob)); // stabilità numerica
      prob = prob / sum(prob);
      z(i) = sample_index(prob) + 1;
    } 
    
    // Aggiorna alpha
    for (int g = 0; g < G; ++g) {
      for (int j = 0; j < d; ++j) {
        vec counts(mj(j), fill::zeros);
        for (int i = 0; i < n; ++i) {
          if (z(i) == (g + 1) && !std::isnan(Y(i, j))) {
            int val = static_cast<int>(Y(i, j)) - 1;
            if (val >= 0 && val < mj(j)) {
              counts(val) += 1;
            }
          } 
        }
        vec new_alpha = rdirichlet(c + counts);
        for (int h = 0; h < mj(j); ++h) {
          alpha(j, h, g) = new_alpha(h);
        }
      } 
    }
    
    // Aggiorna tau
    vec n_g(G, fill::zeros);
    for (int i = 0; i < n; ++i) {
      n_g(z(i) - 1) += 1;
    }
    tau = rdirichlet(b + n_g);
    
    for (int i = 0; i < n; ++i) {
      out_cluster(m, i) = z(i);
    }
  } 
  
  imat cluster_out = out_cluster.rows(nburn, niter - 1);
  
  return List::create(
    Named("cluster") = cluster_out,
    Named("alpha") = alpha,
    Named("tau") = tau
  );
} 


