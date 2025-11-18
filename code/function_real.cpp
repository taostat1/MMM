#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <trng/yarn2.hpp>
#include <trng/uniform_dist.hpp>
#include <trng/normal_dist.hpp>
#include <trng/lcg64.hpp>
#include <trng/discrete_dist.hpp>
#include <trng/gamma_dist.hpp>
#include <omp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(rTRNG)]]
// TRNG >= 4.22 requires C++11
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


double U_theta(const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  vec U = (theta_miss - mu).t() * prec * (theta_miss - mu) / 2.0 -
    sum(log(normcdf(lambda0_miss + lambda1_miss % theta_miss)));
  return U(0);
}


void U_theta_grad(vec &U_grad, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  U_grad = prec * (theta_miss - mu) - 
    (lambda1_miss % normpdf(lambda0_miss + lambda1_miss % theta_miss)) / normcdf(lambda0_miss + lambda1_miss % theta_miss);
}


void HMC_theta(vec &theta_new, double &log_r, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss,
               const vec &lambda1_miss, const vec &p_rnorm, const double &epsilon, const int &num_step) {
  //The kinetic energy has the simplest form sum(p^2)/2
  vec U_grad;
  
  theta_new = theta_miss;
  vec p_new = p_rnorm;
  
  // a half step for momentum at the beginning
  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon *  U_grad / 2.0;
  
  // full steps for position and momentum
  for (int i = 0; i < (num_step-1); i++) {
    theta_new += epsilon * p_new;
    
    U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
    p_new -= epsilon * U_grad;
  }
  theta_new += epsilon * p_new;
  
  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon * U_grad / 2.0;
  
  p_new = -p_new;
  log_r = U_theta(theta_miss, mu, prec, lambda0_miss, lambda1_miss) - U_theta(theta_new, mu, prec, lambda0_miss, lambda1_miss) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
}


void update_theta(mat &theta_t, const mat &ind_zero, const mat &mu_t, const cube &invcov_t, const cube &cov_t, const vec &lambda0_t, 
                  const vec &lambda1_t, const vec &group_t, const int &N, const unsigned long &seed, const double &epsilon = 0.01, 
                  const int &num_step = 50) {
  int i;
  
#pragma omp parallel shared(theta_t) private(i) num_threads(24)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::uniform_dist<> r_unif(0,1);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (i = 0; i < N; i++) {
    
    vec ind_zero_i = ind_zero.row(i).t();
    uvec ind_0 = find(ind_zero_i == true);
    
    if (ind_0.n_elem > 0) {
      vec theta_i = theta_t.row(i).t();
      uvec ind_obs = find(ind_zero_i == false);
      
      vec theta_i_obs = theta_i(ind_obs);
      
      vec theta_i_0 = theta_i(ind_0);
      
      vec mu_i = mu_t.col(group_t(i));
      mat cov_i = cov_t.slice(group_t(i));
      mat invcov_i = invcov_t.slice(group_t(i));
      
      //calculate the acceptance probability
      vec mu_obs = mu_i(ind_obs);
      vec mu_0 = mu_i(ind_0);
      
      mat invcov_i_21 = invcov_i.submat(ind_obs, ind_0);
      
      mat prec_cond = invcov_i.submat(ind_0, ind_0);
      
      mat cov_obs_inv = invcov_i.submat(ind_obs, ind_obs) - invcov_i_21 * inv(prec_cond) * invcov_i_21.t();
      
      mat cov_21 = cov_i.submat(ind_0, ind_obs);
      mat cov_0 = cov_i.submat(ind_0, ind_0);
      
      vec mu_cond = mu_0 + cov_21 * cov_obs_inv * (theta_i_obs - mu_obs);
      
      vec p_rnorm(ind_0.n_elem);
      for (int t = 0; t < ind_0.n_elem; t++) {
        p_rnorm(t) = r_norm(rx);
      }
      double tmp_unif = r_unif(rx);
      
      vec lambda0_miss = lambda0_t(ind_0);
      vec lambda1_miss = lambda1_t(ind_0);
      
      vec theta_star_0;
      double log_r;
      HMC_theta(theta_star_0, log_r, theta_i_0, mu_cond, prec_cond, lambda0_miss, lambda1_miss, p_rnorm, epsilon, num_step);
      
      if (tmp_unif < exp(log_r)) {
        theta_i(ind_0) = theta_star_0;
        theta_t.row(i) = theta_i.t();
      }
    }
  }
}
}


void update_mu(mat &mu_t, const mat &theta_t, const cube &invcov_t, const vec &group_t, const int &G, const int &K, const double &eta_mu = 0, const double &tau_sq_mu = 1) {
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    vec tmp1 = sum(theta_t.rows(ind_k), 0).t();
    
    mat invcov_k = invcov_t.slice(k);
    
    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);
    
    vec mu_k = mvnrnd(MU, COV, 1);
    mu_t.col(k) = mu_k;
  }
}


void update_invcov(cube &invcov_t, cube &cov_t, const cube &edge_t, const mat &theta_t, const mat &mu_t, const vec &group_t, 
                   const double &ssp_v0, const double &ssp_v1, const double &ssp_l, const int &G, const int &K, const unsigned long &seed) {
  int k;
#pragma omp parallel shared(invcov_t, cov_t) private(k) num_threads(K)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    mat theta_k = theta_t.rows(ind_k).t();
    
    mat theta_mu = theta_k.each_col() - mu_t.col(k);
    mat S = theta_mu * theta_mu.t();
    
    mat edge_k = edge_t.slice(k);
    mat V = edge_k * ssp_v1 * ssp_v1;
    uvec ind_n_v1 = find(edge_k == 0);
    V(ind_n_v1).fill(ssp_v0 * ssp_v0); 
    V.diag().fill(0);
    
    mat invcov_k = invcov_t.slice(k);
    mat cov_k = cov_t.slice(k);
    
    vec G_vec = regspace(0, G-1);
    
    for (int g = 0; g < G; g++) {
      uvec ind_ng = find(G_vec != g);
      vec v12_g = V.col(g);
      mat v12_inv = diagmat(1/v12_g(ind_ng));
      vec cov_k_g = cov_k.col(g);
      mat w11_inv = cov_k.submat(ind_ng, ind_ng) - cov_k_g(ind_ng) * cov_k_g(ind_ng).t() / cov_k_g(g);
      mat C = (S(g,g) + ssp_l) * w11_inv + v12_inv;
      
      mat C_chol_inv = inv(chol(C));
      
      vec s12_g = S.col(g);
      vec mu_w12 = -C_chol_inv * C_chol_inv.t() * s12_g(ind_ng);
      
      vec rnorm_vec(G-1);
      for (int t = 0; t < (G-1); t++) {
        rnorm_vec(t) = r_norm(rx);
      }
      vec w12 = mu_w12 + C_chol_inv * rnorm_vec;
      
      trng::gamma_dist<> r_gamma(ind_k.n_elem/2.0 + 1.0, 2.0/(S(g,g)+ssp_l));
      
      double w_v = r_gamma(rx);
      vec w11_inv_w12 = w11_inv * w12;
      vec w22 = w_v + w12.t() * w11_inv_w12;
      
      vec w_update(G);
      w_update(ind_ng) = w12;
      w_update(g) = w22(0);
      
      invcov_k.col(g) = w_update;
      invcov_k.row(g) = w_update.t();
      
      cov_k.submat(ind_ng, ind_ng) = w11_inv + w11_inv_w12 * w11_inv_w12.t() / w_v;
      vec cov_k_ng = - w11_inv_w12 / w_v;
      
      vec cov_update(G);
      cov_update(ind_ng) = cov_k_ng;
      cov_update(g) = 1 / w_v;
      
      cov_k.col(g) = cov_update;
      cov_k.row(g) = cov_update.t();
    }
    cov_t.slice(k) = cov_k;
    invcov_t.slice(k) = invcov_k;
  }
}
}


void update_edge(cube &edge_t, const cube &invcov_t, const double &ssp_v0, const double &ssp_v1, const double &ssp_pi, const int &G, const int &K, const unsigned long &seed) {
  int g;
#pragma omp parallel shared(edge_t) private(g) num_threads(24)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::uniform_dist<> r_unif(0,1);
#pragma omp for schedule(auto)
  for (g = 0; g < G; g++) {
    for (int k = 0; k < K; k++) {
      for (int i = g+1; i < G; i++) {
        double log_p = invcov_t(g, i, k) * invcov_t(g, i, k) * (1.0 / (2.0*ssp_v1*ssp_v1) - 1.0 / (2.0*ssp_v0*ssp_v0)) +
          log(ssp_v1) - log(ssp_v0) + log(1.0 - ssp_pi) - log(ssp_pi);
        double r = 1.0/(1.0 + exp(log_p));
        double tmp = r_unif(rx);
        if (tmp < r) {
          edge_t(g, i, k) = 1;
          edge_t(i, g, k) = 1;
        } else {
          edge_t(g, i, k) = 0;
          edge_t(i, g, k) = 0;
        }
      } 
    }
  }
}
}

double U_lam(const vec &lambda, const vec &theta_miss, const vec &theta_obs, 
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec U = - sum(log(1.0 - normcdf(lambda(0) + lambda(1) * theta_obs))) - sum(log(normcdf(lambda(0) + lambda(1) * theta_miss))) +
    (lambda(0) - lam0_0) * (lambda(0) - lam0_0) / (2.0 * sigma2_lam0) + (lambda(1) - lam1_0) * (lambda(1) - lam1_0) / (2.0 * sigma2_lam1);
  return U(0);
}


void U_lam_grad(vec &U_grad, const vec &lambda, const vec &theta_miss, const vec &theta_obs, 
                const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec part1 = normpdf(lambda(0) + lambda(1) * theta_obs) / (1.0 - normcdf(lambda(0) + lambda(1) * theta_obs));
  vec part2 = - normpdf(lambda(0) + lambda(1) * theta_miss) / normcdf(lambda(0) + lambda(1) * theta_miss);
  
  U_grad(0) = sum(part1) + sum(part2) + (lambda(0) - lam0_0) / sigma2_lam0;
  U_grad(1) = sum(part1 % theta_obs) + sum(part2 % theta_miss) + (lambda(1) - lam1_0) / sigma2_lam1;
}


void HMC_lam(vec &lambda_new, double &log_r, const vec &lambda, const vec &theta_miss, const vec &theta_obs, const vec &p_rnorm,
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1,
             const double &epsilon, const int &num_step) {
  //The kinetic energy has the simplest form sum(p^2)/2
  vec U_grad(2);
  
  lambda_new = lambda;
  vec p_new = p_rnorm;
  
  // a half step for momentum at the beginning
  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon *  U_grad / 2.0;
  // full steps for position and momentum
  for (int i = 0; i < (num_step-1); i++) {
    lambda_new += epsilon * p_new;
    
    U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
    p_new -= epsilon * U_grad;
  }
  lambda_new += epsilon * p_new;
  
  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon * U_grad / 2.0;
  
  p_new = -p_new;
  log_r = U_lam(lambda, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) - 
    U_lam(lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
  
}



void update_lambda(vec &lambda0_t, vec &lambda1_t, const mat &theta_t, const mat &ind_zero, const int &G, const unsigned long &seed,
                   const double &lam0_0=2, const double &lam1_0=-2, const double &sigma2_lam0=0.25, const double &sigma2_lam1=0.25,
                   const double &epsilon = 0.01, const int &num_step = 50) {
  int g;
  
#pragma omp parallel shared(lambda0_t, lambda1_t) private(g) num_threads(24)
{
  trng::yarn2 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  trng::uniform_dist<> r_unif(0,1);
  trng::normal_dist<> r_norm(0,1);
#pragma omp for schedule(auto)
  for (g = 0; g < G; g++) {
    vec lambda_g(2);
    lambda_g(0) = lambda0_t(g);
    lambda_g(1) = lambda1_t(g);
    
    vec theta_g = theta_t.col(g);
    uvec ind_0= find(ind_zero.col(g) == true);
    uvec ind_obs = find(ind_zero.col(g) == false);
    
    vec theta_g_obs = theta_g(ind_obs);
    
    vec theta_g_miss = theta_g(ind_0);
    
    vec p_rnorm(2);
    p_rnorm(0) = r_norm(rx);
    p_rnorm(1) = r_norm(rx);
    
    double tmp_unif = r_unif(rx);
    
    vec lambda_g_new;
    double log_r;
    HMC_lam(lambda_g_new, log_r, lambda_g, theta_g_miss, theta_g_obs, p_rnorm, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon, num_step);
    
    if (tmp_unif < exp(log_r) && lambda_g_new(1) < 0) {
      lambda0_t(g) = lambda_g_new(0);
      lambda1_t(g) = lambda_g_new(1);
    }
  }
}
}

vec rDirichlet(const vec &alpha_vec) {
  vec tmp(alpha_vec.n_elem, fill::zeros);
  for (unsigned int i = 0; i < alpha_vec.n_elem; i++) {
    tmp(i) = log(randg(distr_param(alpha_vec(i), 1.0)));
  }
  tmp = tmp - max(tmp);
  tmp = exp(tmp);
  tmp = tmp/sum(tmp);
  return tmp;
}


void update_pi(vec &pi_t, const vec &group_t, const vec &gam, const int &K) {
  vec tmp0(K, fill::zeros);
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  pi_t = rDirichlet(tmp0+gam);
}


void update_group(vec &group_t, const cube &invcov_t, const mat &theta_t, const mat &mu_t, const vec &pi_t, const int &G, const int &N, const int &K, const unsigned long &seed) {
  uvec pi_n0 = find(pi_t > 0);
  int i;
  vec invcov_logdet(K);
  for (int k = 0; k < K; k++) {
    invcov_logdet(k) = log(det(invcov_t.slice(k)))/2.0;
  }
#pragma omp parallel shared(group_t, pi_n0, invcov_logdet) private(i) num_threads(24)
{
  trng::lcg64 rx;
  rx.seed(seed);
  int size = omp_get_num_threads();
  int rank = omp_get_thread_num();
  rx.split(size, rank);
  
#pragma omp for schedule(auto) 
  for (i = 0; i < N; i++) {
    vec theta_i = theta_t.row(i).t();
    vec tmp(K);
    tmp.fill(- datum::inf);
    for (unsigned int k_pi = 0; k_pi < pi_n0.n_elem; k_pi++) {
      vec tmp_k_pi = invcov_logdet(pi_n0(k_pi)) - 
        (theta_i - mu_t.col(pi_n0(k_pi))).t() * invcov_t.slice(pi_n0(k_pi)) * (theta_i - mu_t.col(pi_n0(k_pi)))/2.0;
      tmp(pi_n0(k_pi)) = tmp_k_pi(0);
    }
    tmp.replace(datum::inf, pow(10, 308));
    vec tmp_new = tmp - max(tmp);
    tmp_new = exp(tmp_new);
    vec prob = tmp_new % pi_t;
    prob = prob / sum(prob);
    
    trng::discrete_dist dist_C(prob.begin(), prob.end());
    group_t(i) = dist_C(rx);
    
  }
}
}


// [[Rcpp::export]]
List MCMC_full(const int num_iter, const int num_save, mat theta_t, mat ind_zero, mat mu_t,
               cube invcov_t, cube cov_t, cube edge_t, vec group_t, vec lambda0_t, vec lambda1_t, 
               vec pi_t, vec gam, const int G, const int N, const int K, 
               double ssp_v0, const double ssp_v1, const double ssp_l, double ssp_pi,
               double epsilon_theta = 0.2, int num_step_theta = 20,
               double eta_mu = 0, double tau_sq_mu = 1,
               double lam0_0=2, double lam1_0=-2, double sigma2_lam0=0.25, 
               double sigma2_lam1=0.25, double epsilon_lam = 0.01, int num_step_lam = 10) {
  group_t = group_t - 1;
  
  int save_start = num_iter - num_save;
  mat group_save(K, N, fill::zeros);
  mat mu_save(G, K, fill::zeros);
  cube invcov_save(G, G, K, fill::zeros);
  // cube cov_save(G, G, K, fill::zeros);
  cube edge_save(G, G, K, fill::zeros);
  
  mat theta_save(N, G, fill::zeros);
  
  mat group_iter(N, num_save);
  cube mu_iter(G, K, num_save);
  
  arma::field<cube> invcov_iter(num_save);
  //save the results in the iterations for the first column of each matrix
  cube edge_iter(G, K, num_save);
  mat theta_iter(N, num_save);
  ////////////////////////////////
  
  mat pi_iter(K, num_save);
  mat lam0_iter(G, num_save);
  mat lam1_iter(G, num_save);
  
  long seed = randi(distr_param(1,1000000));
  
  for (int t_iter = 0; t_iter < num_iter; t_iter++) {
    seed = randi(distr_param(1,1000000));
    update_theta(theta_t, ind_zero, mu_t, invcov_t, cov_t, lambda0_t, lambda1_t, group_t, N, seed, epsilon_theta, num_step_theta);
    
    seed = randi(distr_param(1,1000000));
    update_lambda(lambda0_t, lambda1_t, theta_t, ind_zero, G, seed, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon_lam, num_step_lam);
    
    update_mu(mu_t, theta_t, invcov_t, group_t, G, K, eta_mu, tau_sq_mu);
    
    seed = randi(distr_param(1,1000000)); 
    update_invcov(invcov_t, cov_t, edge_t, theta_t, mu_t, group_t, ssp_v0, ssp_v1, ssp_l, G, K, seed);
    
    seed = randi(distr_param(1,1000000)); 
    update_edge(edge_t, invcov_t, ssp_v0, ssp_v1, ssp_pi, G, K, seed);
    
    seed = randi(distr_param(1,1000000)); 
    update_group(group_t, invcov_t, theta_t, mu_t, pi_t, G, N, K, seed);
    
    update_pi(pi_t, group_t, gam, K);
    
    if (t_iter >= save_start) {
      int save_i = t_iter - save_start;
      for (int C_i = 0; C_i < N; C_i++) {
        group_save(group_t(C_i),C_i) += 1;
      }
      
      mu_save += mu_t;
      invcov_save += invcov_t;
      edge_save += edge_t;
      theta_save += theta_t;
      
      // cov_save += cov_t;
      // pi_save += pi_t;
        
      group_iter.col(save_i) = group_t;
      mu_iter.slice(save_i) = mu_t;
      invcov_iter(save_i) = invcov_t;
      edge_iter.slice(save_i) = edge_t.col(0);
      theta_iter.col(save_i) = theta_t.col(0);
      
      lam0_iter.col(save_i) = lambda0_t;
      lam1_iter.col(save_i) = lambda1_t;
      
      pi_iter.col(save_i) = pi_t;
    }
  }
  
  mu_save /= num_save;
  invcov_save /= num_save;
  edge_save /= num_save;
  theta_save /= num_save;
  // pi_save /= num_save;
  // cov_save /= num_save;
  
  vec group_est(N);
  for (int C_i = 0; C_i < N; C_i++) {
    group_est(C_i) = group_save.col(C_i).index_max() + 1;
  }
  
  return List::create(Named("group_est")=group_est, Named("mu_est")=mu_save, 
                      Named("invcov_est")=invcov_save, Named("edge_est")=edge_save, 
                      Named("theta_est")=theta_save, Named("pi_iter")=pi_iter,
                      // Named("cov_est")=cov_save,
                      Named("mu_iter")=mu_iter, 
                      Named("invcov_iter")=invcov_iter, Named("lam0_iter")=lam0_iter, 
                      Named("lam1_iter")=lam1_iter,
                      Named("group_iter")=group_iter, Named("edge_iter")=edge_iter,
                      Named("theta_iter")=theta_iter);
}





// [[Rcpp::export]]
mat update_mu_R(mat theta_t, cube invcov_t, vec group_t, int G, int K, double eta_mu = 0, double tau_sq_mu = 1) {
  mat mu_updated(G, K);
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    vec tmp1 = sum(theta_t.rows(ind_k), 0).t();
    
    mat invcov_k = invcov_t.slice(k);
    
    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);
    
    vec mu_k = mvnrnd(MU, COV, 1);
    mu_updated.col(k) = mu_k;
  }
  return mu_updated;
}

// [[Rcpp::export]]
vec update_pi_R(vec group_t, vec gam, int K) {
  vec pi_new(K);
  vec tmp0(K, fill::zeros);
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  pi_new = rDirichlet(tmp0+gam);
  return pi_new;
}
