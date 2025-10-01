#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Define the T_function using Armadillo
// [[Rcpp::export]]
double T_function(double x) {
  double Tx = - std::pow(2 * M_PI, -0.5) * exp(-x * x / 2);
  return Tx;
}

// Define the delta_T function using Armadillo
// [[Rcpp::export]]
double delta_T(double xi_n, double xi_n_minus, double mu_k, double sigma2_k) {
  double result = T_function((xi_n - mu_k) / sqrt(sigma2_k)) - T_function((xi_n_minus - mu_k) / sqrt(sigma2_k));
  return result;
}


// [[Rcpp::export]]
double delta_F(double xi_n, double xi_n_minus, double mu_k, double sigma2_k) {
  double f1 = R::pnorm(xi_n, mu_k, sqrt(sigma2_k), true, false); // true for lower.tail
  double f2 = R::pnorm(xi_n_minus, mu_k, sqrt(sigma2_k), true, false); // true for lower.tail
  double result = f1 - f2;
  // Rcout<<"delta_F"<<std::endl;
  // Rcout<<"f1"<<f1<<std::endl;
  // Rcout<<"f2"<<f2<<std::endl;
  return result;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double F(double x, arma::vec w, arma::vec u, arma::vec s) {
  int n = w.n_elem;
  double p_sum = 0;
  for(int i=0; i<n; i++){
    // double mean_i = u[i];
    double p = R::pnorm(x, u[i], s[i],true, false);
    p_sum = p_sum + w[i] * p;
  }
  return p_sum;
}

double G_func(double x, List param){
  double p = param["p"];
  double G = F(x, param["w"], param["u"], param["s"]) - p;
  return G;
}

// Define the uniroot-like function in C++
// [[Rcpp::export]]
double myUniroot(List func_param, double lower, double upper, double tol = 1.0e-9) {
  double f_lower = G_func(lower,func_param);
  // double f_upper = G_func(upper,func_param);
  while (std::abs(upper - lower) > tol) {
    double mid = (lower + upper) / 2.0;
    double f_mid = G_func(mid,func_param);
    if (f_mid == 0.0) {
      // Rcout<<"case 1"<<std::endl;
      return mid;
    } else if (f_lower * f_mid <= 0) {
      // Rcout<<"case 2"<<std::endl;
      upper = mid;
      // f_upper = f_mid;
    } else {
      // Rcout<<"-----case 3-----"<<std::endl;
      // Rcout<<"f_lower"<<f_lower<<std::endl;
      // Rcout<<"f_mid"<<f_mid<<std::endl;
      lower = mid;
      f_lower = f_mid;
    }
    // Rcout<<"mid"<<mid<<std::endl;
  }
  return (lower + upper) / 2.0;
}


// Define the F_inv function using myUniroot for root finding
// [[Rcpp::export]]
double F_inv(double p, arma::vec w, arma::vec u, arma::vec s) {
  List param_list = List::create(Named("p") = p , Named("w") = w, Named("u") = u, Named("s") = s);
  // Function G_func("G_func");
  double root = myUniroot(param_list, -2000, 2000);
  return root;
}



// Define the W2 function using Armadillo
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double W2_old(const arma::vec data, const arma::vec c_alloc, List prior_list) {
  double alpha = as<double>(prior_list["alpha"]);
  double kappa0 = as<double>(prior_list["kappa0"]);
  double theta = as<double>(prior_list["theta"]);
  double a = as<double>(prior_list["a"]);
  double b = as<double>(prior_list["b"]);
  
  arma::vec cluster = arma::unique(c_alloc);
  int K = cluster.size();
  int N = data.size();
  arma::vec data_ordered = arma::sort(data);
  
  arma::vec alpha_post(K, arma::fill::zeros);
  arma::vec mu_post(K, arma::fill::zeros);
  arma::vec sigma2_post(K, arma::fill::zeros);
  
  // posterior mean
  for (int k = 0; k < K; ++k) {
    arma::uvec I = arma::find(c_alloc == cluster(k));
    int nk = I.size();
    arma::vec yk = data.elem(I);
    alpha_post(k) = alpha + nk;
    double kappa_post = kappa0 + (double)nk;
    mu_post(k) = (kappa0 * theta + arma::sum(yk)) / kappa_post;
    double a_post = a + ((double)nk + 1) * 0.5;
    double b_post = b + 0.5 * arma::sum(arma::pow(yk - arma::mean(yk), 2)) + 0.5 * kappa0 * (double)nk * std::pow(arma::mean(yk) - theta, 2) / kappa_post;
    sigma2_post(k) = b_post / (a_post - 1);
  }
  arma::vec weights = alpha_post / arma::sum(alpha_post);
  
  
  arma::mat F_mat(N, K, arma::fill::zeros);
  arma::mat T_mat(N, K, arma::fill::zeros);
  arma::vec part2(K, arma::fill::zeros);
  arma::vec part3(K, arma::fill::zeros);
  
  for (int i_n = 0; i_n < N; ++i_n) {
    double i_double = (double)i_n;
    
    // Rcout<<"input4"<<sqrt(sigma2_post)<<std::endl;
    
    double xi_n = F_inv((i_double + 1.0) / N, weights, mu_post, sqrt(sigma2_post));
    double xi_n_minus = F_inv(i_double / N, weights, mu_post, sqrt(sigma2_post));
    
    // Rcout<<"xi_n"<<xi_n<<std::endl;
    // Rcout<<"xi_n_minus"<<xi_n_minus<<std::endl;
    
    for (int k = 0; k < K; ++k) {
      F_mat(i_n, k) = delta_F(xi_n, xi_n_minus, mu_post(k), sigma2_post(k));
      T_mat(i_n, k) = delta_T(xi_n, xi_n_minus, mu_post(k), sigma2_post(k));
    }
  }
  
  
  
  arma::vec sum_F_k(K, arma::fill::zeros);
  arma::vec sum_T_k(K, arma::fill::zeros);
  
  for (int k = 0; k < K; ++k) {
    for (int i_n = 0; i_n < N; ++i_n) {
      sum_F_k(k) += data_ordered(i_n) * F_mat(i_n, k);
      sum_T_k(k) += data_ordered(i_n) * T_mat(i_n, k);
    }
    part2(k) = std::pow(mu_post(k), 2) + sigma2_post(k);
    part3(k) = mu_post(k) * sum_F_k(k) + sqrt(sigma2_post[k]) * sum_T_k(k);
    // Rcout<<"mu_post(k) * sum_F_k(k)"<<(mu_post(k) * sum_F_k(k))<<std::endl;
    // Rcout<<"std::pow(sigma2_post(k), 1/2) * sum_T_k(k)"<< (std::pow(sigma2_post(k), 1/2) * sum_T_k(k))<<std::endl;
    // Rcout<<"sqrt(sigma2_post[k]) * sum_T_k[k]"<<sqrt(sigma2_post[k]) * sum_T_k[k]<<std::endl;
    // Rcout<<"part3(k)"<<part3(k)<<std::endl;
  }
  
  arma::mat F_sub = F_mat.submat(0, 0, 2, K-1);
  arma::mat T_sub = T_mat.submat(0, 0, 2, K-1);
  
  // Rcout<<"F_mat-Rcpp"<<F_mat<<std::endl;
  // Rcout<<"T_mat-Rcpp"<<T_sub<<std::endl;
  
  // Rcout<<"part3"<<part3<<std::endl;
  
  
  // double result = - 2 * arma::sum(weights % part3);
  double result = (1/(double)N) * arma::dot(data, data) + arma::sum(weights % part2) - 2 * arma::sum(weights % part3);
  return result;
  
}

// Define the W2 function using Armadillo, now using the posterior parameters already as input
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double W2(const arma::vec data, const arma::vec weights_post, arma::vec mu_post, arma::vec sigma2_post) {
  
  int N = data.size();
  int K = weights_post.size();
  arma::vec data_ordered = arma::sort(data);
  
  arma::mat F_mat(N, K, arma::fill::zeros);
  arma::mat T_mat(N, K, arma::fill::zeros);
  arma::vec part2(K, arma::fill::zeros);
  arma::vec part3(K, arma::fill::zeros);
  
  for (int i_n = 0; i_n < N; ++i_n) {
    double i_double = (double)i_n;
    
    // Rcout<<"input4"<<sqrt(sigma2_post)<<std::endl;
    
    double xi_n = F_inv((i_double + 1.0) / N, weights_post, mu_post, sqrt(sigma2_post));
    double xi_n_minus = F_inv(i_double / N, weights_post, mu_post, sqrt(sigma2_post));
    
    // Rcout<<"xi_n"<<xi_n<<std::endl;
    // Rcout<<"xi_n_minus"<<xi_n_minus<<std::endl;
    
    for (int k = 0; k < K; ++k) {
      F_mat(i_n, k) = delta_F(xi_n, xi_n_minus, mu_post(k), sigma2_post(k));
      T_mat(i_n, k) = delta_T(xi_n, xi_n_minus, mu_post(k), sigma2_post(k));
    }
  }
  
  arma::vec sum_F_k(K, arma::fill::zeros);
  arma::vec sum_T_k(K, arma::fill::zeros);
  
  for (int k = 0; k < K; ++k) {
    for (int i_n = 0; i_n < N; ++i_n) {
      sum_F_k(k) += data_ordered(i_n) * F_mat(i_n, k);
      sum_T_k(k) += data_ordered(i_n) * T_mat(i_n, k);
    }
    part2(k) = std::pow(mu_post(k), 2) + sigma2_post(k);
    part3(k) = mu_post(k) * sum_F_k(k) + sqrt(sigma2_post[k]) * sum_T_k(k);
    // Rcout<<"mu_post(k) * sum_F_k(k)"<<(mu_post(k) * sum_F_k(k))<<std::endl;
    // Rcout<<"std::pow(sigma2_post(k), 1/2) * sum_T_k(k)"<< (std::pow(sigma2_post(k), 1/2) * sum_T_k(k))<<std::endl;
    // Rcout<<"sqrt(sigma2_post[k]) * sum_T_k[k]"<<sqrt(sigma2_post[k]) * sum_T_k[k]<<std::endl;
    // Rcout<<"part3(k)"<<part3(k)<<std::endl;
  }
  
  arma::mat F_sub = F_mat.submat(0, 0, 2, K-1);
  arma::mat T_sub = T_mat.submat(0, 0, 2, K-1);
  
  // Rcout<<"F_mat-Rcpp"<<F_mat<<std::endl;
  // Rcout<<"T_mat-Rcpp"<<T_sub<<std::endl;
  
  // Rcout<<"part3"<<part3<<std::endl;
  
  
  // double result = - 2 * arma::sum(weights % part3);
  double result = (1/(double)N) * arma::dot(data, data) + arma::sum(weights_post % part2) - 2 * arma::sum(weights_post % part3);
  return result;
  
}


//// [[Rcpp::depends(RcppArmadillo)]]
//// [[Rcpp::export]]
// arma::vec min_W2(const arma::mat c_samples, const arma::vec data, List prior_list) {
//   int M = c_samples.n_rows;
//   arma::vec W2_vec(M);
//   
//   for (int i = 0; i < M; ++i) {
//     if (i % 100 == 0) {
//       Rcpp::Rcout << i << std::endl;
//     }
//     W2_vec(i) = W2(data, c_samples.row(i).t(), prior_list);
//   }
//   // Rcout<<"W2"<<W2_vec<<std::endl;
//   
//   int min_m = arma::index_min(W2_vec);
//   return c_samples.row(min_m).t();
// }