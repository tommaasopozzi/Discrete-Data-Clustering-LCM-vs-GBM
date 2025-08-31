#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List multiclass_centroid_C(const arma::mat& x_k,  const arma::mat& x, double alpha = 0.5) {
  int n_k = x_k.n_rows;
  int d = x_k.n_cols;
  
  // Troviamo il numero massimo di classi per ogni colonna
  arma::ivec C(d);
  for (int j = 0; j < d; j++) {
    C[j] = static_cast<int>(max(x.col(j)));  
  }
  
  // Creiamo la matrice delle probabilità
  int max_C = max(C);
  arma::mat p(max_C, d, arma::fill::zeros);
  
  // Calcoliamo P(x_{kj} = c)
  for (int j = 0; j < d; j++) {
    for (int c = 1; c <= C[j]; c++) {
      double n_jc = accu(x_k.col(j) == c);
      p(c - 1, j) = (n_jc + alpha) / (n_k + C[j] * alpha);
    }
  } 
  
  return List::create(Named("matrix_prob") = p, Named("classi") = C);
} 

// [[Rcpp::export]]
double loss_multiclass_C(const arma::mat& x_k, const arma::mat& x,double alpha = 0.5) {
  List result = multiclass_centroid_C(x_k, x, alpha);
  arma::mat x_tilde = as<arma::mat>(result["matrix_prob"]);
  arma::ivec C = as<arma::ivec>(result["classi"]);
  
  int n_k = x_k.n_rows;
  int d = x_k.n_cols;
  double loss = 0.0;
  
  double epsilon = 1e-10;
  
  for (int i = 0; i < n_k; i++) {
    for (int j = 0; j < d; j++) {
      int c = static_cast<int>(x_k(i, j));//statistic_cast trasforma il valore della matrice in intero 
      double prob = std::max(x_tilde(c - 1, j), epsilon); 
      loss -= log(prob);
    }
  } 
  
  return loss;
} 


// [[Rcpp::export]]
List Gibbs_kmulticlass_C(int R, const arma::mat& X, arma::vec G, arma::vec freq, double lambda, bool trace){
  
  // Inizializzo:
  //n: numero di osservazioni 
  //K: numero di cluster 
  //R: numero di iterazioni campionamento MCMC
  //prob: vettore in cui inserire P(S_i = j|-)
  //clusters: vettore di possibili etichette del cluster 
  
  int n = X.n_rows;
  int K = freq.n_elem;

  int R_show = floor(R/5);
  
  arma::vec prob(K);  
  arma::vec lprob(K);
  IntegerVector clusters = Range(1,K);
  bool skip;
  
  // Inizializzo le matrici di output
  arma::mat G_out(R,n);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;
  
  //Calcolo della loss function per i k clusters 
  double total_loss=0;
  for(int j=0; j < K; j++){
    arma::uvec idx_cluster; idx_cluster = find(G == j + 1); 
    total_loss = total_loss + loss_multiclass_C(X.rows(idx_cluster), X);
  } 
  
  
  // Inizio Gibbs Sampling
  for(int r = 0; r < R; r++){
    
    //Devo aggiornare le probabilità di allocazione per ogni osservazione
    for(int i = 0; i < n; i++) {
      
      // Rimuovo l'i-esima osservazione dal suo cluster
      freq(G(i)-1) = freq(G(i)-1) - 1; 
      
      // Controllo se n_k(i) = 1
      if(freq(G(i)-1) == 0){
        skip = true;
      } else { 
        skip = false;
      } 
      
      if(!skip){
        // Alloco l'i-esima osservazione in un nuovo cluster o nel suo precedente 
        G(i) = arma::datum::nan;
        for(int j = 0; j < K; j++){
          idx_cluster = find(G == j + 1); 
          arma::uvec idx_cluster_i(freq(j)+1); idx_cluster_i.head(freq(j)) = idx_cluster; idx_cluster_i.tail(1) = i;
          lprob(j) = - lambda*(loss_multiclass_C(X.rows(idx_cluster_i), X) - loss_multiclass_C(X.rows(idx_cluster), X));
        } 
        
        lprob = lprob - max(lprob); prob  = exp(lprob); prob  = prob/sum(prob);
        G(i) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0]; //Nuovo cluster di appartenenza
      } 
      freq(G(i)-1) = freq(G(i)-1) + 1;
    } 

    //Aggiorno la loss
    double total_loss=0;
    for(int j=0; j < K; j++){
      arma::uvec idx_cluster; idx_cluster = find(G == j + 1); 
      total_loss = total_loss + loss_multiclass_C(X.rows(idx_cluster), X);
    } 
    
    if(trace) {if((r+1)%R_show==0) {Rprintf("Iteration: %i \n", r+1);}}
    
    G_out.row(r) = trans(G);
    total_loss_out(r) = total_loss;
  } 
  return(List::create(Named("G") = G_out, Named("loss") = total_loss_out));
} 



// [[Rcpp::export]]
List kmulticlass_C(const arma::mat& x, arma::vec G, arma::vec freq, int k, bool trace){
  
  // Initialization
  int n = x.n_rows;
  int d = x.n_cols;
  
  // Output is a huge matrix
  //arma::vec G(n);
  arma::uvec idx_cluster;
  arma::vec losses(k);
  bool convergence = false;
  double total_loss_old = -1; 
  double total_loss;
  //std::vector<arma::mat> m(k);
  
  // Initialization
  total_loss = 0;
  for(int j = 0; j < k; j++){
    // Identify the elements within the cluster
    idx_cluster = find(G == j + 1); 
    //m[j] = multiclass_centroid_C(x.rows(idx_cluster), 0.5)["matrix_prob"];
    //Rprintf("n_k: %f \n", m(j,1));
    total_loss = total_loss + loss_multiclass_C(x.rows(idx_cluster), x);
  } 
  
  // Iterations
  while(!convergence){
    // Cluster allocation
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++){
        losses(j) = -loss_multiclass_C(x.rows(idx_cluster), x);
      } 
      // Compute the probability to move into another cluster
      freq(G(i)-1) += -1;
      
      // Allocate the element only if the new frequency is positive.
      if(freq(G(i)-1) > 0){
        G(i) = index_min(losses) + 1;
      } 
      freq(G(i)-1) += 1;
    } 
    
    total_loss = 0;
    for(int j = 0; j < k; j++){
      // Identify the elements within the cluster
      idx_cluster = find(G == j + 1); 
      //m[j] = multiclass_centroid_C(x.rows(idx_cluster), 0.5)["matrix_prob"];
      //Rprintf("n_k: %f \n", m(j,1));
      total_loss = total_loss + loss_multiclass_C(x.rows(idx_cluster), x);
    } 
    
    
    if((total_loss_old - total_loss) == 0){convergence=true;} else {
      if(trace) {Rprintf("Loss function: %f \n", total_loss);}
      total_loss_old = total_loss;
    } 
    
  } 
  return(List::create(Named("cluster") = G, Named("loss") = total_loss));
} 


// [[Rcpp::export]]
double loss_binary_C(const arma::mat& x_cluster){
  int nk = x_cluster.n_rows;
  double nkd = double(nk);
  
  arma::mat x_tilde = nkd /(nkd+1)*arma::mean(x_cluster) + 0.5/(nkd+1);
  arma::mat a = log(x_tilde / (1 - x_tilde));
  arma::mat b = log(1 - x_tilde);
  double out = 0;
  for(int i=0; i < nk; i++){out += - sum(x_cluster.row(i) % a + b);}
  return(out);
}


// [[Rcpp::export]]
List Gibbs_kbinary_C(int R, const arma::mat& X, arma::vec G, arma::vec freq, double lambda, bool trace){
  
  // Initialization
  int n = X.n_rows;
  int K = freq.n_elem;
  
  // How many times the output is displayed
  int R_show = floor(R/5);
  
  // Internal vectors and quantities
  arma::vec prob(K);  // Probabilities of allocation
  arma::vec lprob(K); // Log-probabilities of allocation
  IntegerVector clusters = Range(1,K); // Possible clusters
  bool skip;
  
  // Output is a huge matrix 
  arma::mat G_out(R,n);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;
  
  // Update the total_loss
  double total_loss=0;
  for(int j=0; j < K; j++){
    arma::uvec idx_cluster; idx_cluster = find(G == j + 1); 
    total_loss = total_loss + loss_binary_C(X.rows(idx_cluster));
  }
  
  
  // Cycle of the Gibbs Sampling
  for(int r = 0; r < R; r++){
    
    // Cycle of the Observations
    for(int i = 0; i < n; i++) {
      
      // Which is the old frequency?
      freq(G(i)-1) = freq(G(i)-1) - 1; // Reduces the associated frequencies
      
      // Is the cluster of G_i a singleton?
      if(freq(G(i)-1) == 0){
        skip = true;
      } else {
        skip = false;
      }
      
      if(!skip){
        // The ith element is not allocated
        G(i) = arma::datum::nan;
        // Compute the probability to move into another cluster
        for(int j = 0; j < K; j++){
          // Identify the elements within the cluster
          idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
          arma::uvec idx_cluster_i(freq(j)+1); idx_cluster_i.head(freq(j)) = idx_cluster; idx_cluster_i.tail(1) = i;
          lprob(j) = - lambda*(loss_binary_C(X.rows(idx_cluster_i)) - loss_binary_C(X.rows(idx_cluster)));
        }
        
        // Log-sum-exp trick
        lprob = lprob - max(lprob); prob  = exp(lprob); prob  = prob/sum(prob);
        // Sample the new value
        G(i) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
      }
      // Update the frequency
      freq(G(i)-1) = freq(G(i)-1) + 1;
    }
    
    // Update the total_loss and the w parameter
    double total_loss=0;
    for(int j=0; j < K; j++){
      // Update the total loss
      arma::uvec idx_cluster; idx_cluster = find(G == j + 1); // Now G_i has been
      total_loss = total_loss + loss_binary_C(X.rows(idx_cluster));
    }
    
    if(trace) {if((r+1)%R_show==0) {Rprintf("Iteration: %i \n", r+1);}}
    
    G_out.row(r) = trans(G);
    total_loss_out(r) = total_loss;
  }
  return(List::create(Named("G") = G_out, Named("loss") = total_loss_out));
}





// [[Rcpp::export]]
List Gibbs_kmulticlass_prior_C(int R, const arma::mat& X, arma::vec G, arma::vec freq, double lambda, arma::vec beta, bool trace){
  
  // Inizializzo:
  //n: numero di osservazioni 
  //K: numero di cluster 
  //R: numero di iterazioni campionamento MCMC
  //prob: vettore in cui inserire P(S_i = j|-)
  //clusters: vettore di possibili etichette del cluster 
  
  int n = X.n_rows;
  int K = freq.n_elem;
  
  int R_show = floor(R/5);
  
  arma::vec prob(K);  
  arma::vec lprob(K);
  IntegerVector clusters = Range(1,K);
  bool skip;
  
  // Inizializzo le matrici di output
  arma::mat G_out(R,n);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;
  
  //Calcolo della loss function per i k clusters 
  double total_loss=0;
  for(int j=0; j < K; j++){
    arma::uvec idx_cluster; idx_cluster = find(G == j + 1); 
    total_loss = total_loss + loss_multiclass_C(X.rows(idx_cluster), X);
  }  
  
  
  // Inizio Gibbs Sampling
  for(int r = 0; r < R; r++){
    
    //Devo aggiornare le probabilità di allocazione per ogni osservazione
    for(int i = 0; i < n; i++) {
      
      // Rimuovo l'i-esima osservazione dal suo cluster
      freq(G(i)-1) = freq(G(i)-1) - 1; 
      
      // Controllo se n_k(i) = 1
      if(freq(G(i)-1) == 0){
        skip = true;
      } else {  
        skip = false;
      }  
      
      if(!skip){
        // Alloco l'i-esima osservazione in un nuovo cluster o nel suo precedente 
        G(i) = arma::datum::nan;
        for(int j = 0; j < K; j++){
          idx_cluster = find(G == j + 1); 
          int n_j = idx_cluster.n_elem;
          double diff = lgamma(beta(j) + n_j - 1) - lgamma(beta(j) + n_j);;
          arma::uvec idx_cluster_i(freq(j)+1); idx_cluster_i.head(freq(j)) = idx_cluster; idx_cluster_i.tail(1) = i;
          lprob(j) = - lambda*(diff + loss_multiclass_C(X.rows(idx_cluster_i), X) - loss_multiclass_C(X.rows(idx_cluster), X));
        }  
        

        
        lprob = lprob - max(lprob); prob  = exp(lprob); prob  = prob/sum(prob);
        
        if (any(prob != prob)) {
          Rcpp::Rcout << "NaN in prob at iteration " << r << ", obs " << i << std::endl;
          Rcpp::Rcout << "lprob: " << lprob.t() << std::endl;
          Rcpp::Rcout << "beta: " << beta.t() << std::endl;
          Rcpp::Rcout << "freq: " << freq.t() << std::endl;
        }
        G(i) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0]; //Nuovo cluster di appartenenza
      }  
      freq(G(i)-1) = freq(G(i)-1) + 1;
    } 
     
    //Aggiorno la loss
    double total_loss=0;
    for(int j=0; j < K; j++){
      arma::uvec idx_cluster; idx_cluster = find(G == j + 1); 
      total_loss = total_loss + loss_multiclass_C(X.rows(idx_cluster), X);
    }  
    
    if(trace) {if((r+1)%R_show==0) {Rprintf("Iteration: %i \n", r+1);}}
    
    G_out.row(r) = trans(G);
    total_loss_out(r) = total_loss;
  }  
  return(List::create(Named("G") = G_out, Named("loss") = total_loss_out));
}  