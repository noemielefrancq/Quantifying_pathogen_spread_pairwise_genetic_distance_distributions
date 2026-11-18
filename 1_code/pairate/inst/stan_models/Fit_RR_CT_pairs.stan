functions {
  matrix transition_matrix(matrix Q, real t) {
    matrix[rows(Q),rows(Q)] res;
    if (rows(Q) > 2) {
      res = matrix_exp(Q*t);
    } else {
      real a = Q[1,2];
      real b = Q[2,1];
      real e = exp(-(a+b)*t);
      res[1,1] = (b+a*e)/(a+b);
      res[1,2] = a*(1-e)/(a+b);
      res[2,1] = b*(1-e)/(a+b);
      res[2,2] = (a+b*e)/(a+b);
    }
    return res;
  }

  real pair_loglik(array[] int y, vector D_ref, vector cr, real n_diff, real alpha) {
    array[num_elements(D_ref)] int keep;
    int n = 0;
    for (bin in 1:num_elements(D_ref)) {
      if (D_ref[bin] > 0) {
        n += 1;
        keep[n] = bin;
      }
    }
    vector[n] P = D_ref[keep[1:n]].*cr[keep[1:n]];
    real phi = (1-alpha)/alpha;
    if (is_nan(sum(P)) || is_inf(sum(P)) || min(P) == 0 || phi == 0 || is_inf(phi)) return negative_infinity();
    return dirichlet_multinomial_lpmf(y[keep[1:n]] | P/sum(P)*phi)+neg_binomial_2_lpmf(sum(y) | n_diff*sum(P), phi);
  }
}

data {
  int<lower=2> K;
  int<lower=1> N_pairs;
  array[N_pairs] int<lower=1, upper=K> pair_a;
  array[N_pairs] int<lower=1, upper=K> pair_b;
  int<lower=1> N_rate_params;
  array[K,K] int<lower=0, upper=N_rate_params> rate_map;

  int<lower=1> max_N_bins;
  array[N_pairs] int<lower=1, upper=max_N_bins> N_bins_pair;
  vector[max_N_bins] bin_break_min;
  vector[max_N_bins] bin_break_max;

  int<lower=1> N_grid;
  int<lower=1> N_boot;
  real<lower=0> rate_dispersion_prior;

  array[N_pairs,max_N_bins,N_boot] int<lower=0> hist_diff;
  array[N_pairs,max_N_bins,N_boot] int<lower=0> hist_within_a;
  array[N_pairs,max_N_bins,N_boot] int<lower=0> hist_within_b;
  array[N_pairs,N_boot] real<lower=1> n_diff;
  array[N_pairs,N_boot] real<lower=1> n_within_a;
  array[N_pairs,N_boot] real<lower=1> n_within_b;

  int<lower=1> N_sets;
  array[N_sets] int<lower=1, upper=N_pairs> set_pair;
  array[N_sets] vector<lower=0, upper=1>[3] ref_weight;
  vector<lower=0>[N_sets] set_weight;
}

transformed data {
  array[N_sets,N_boot] vector[max_N_bins] D_ref;
  for (s in 1:N_sets) {
    int p = set_pair[s];
    for (boot in 1:N_boot) {
      D_ref[s,boot] = ref_weight[s,1]*to_vector(hist_within_a[p,:,boot])/n_within_a[p,boot]+ref_weight[s,2]*to_vector(hist_within_b[p,:,boot])/n_within_b[p,boot]+ref_weight[s,3]*to_vector(hist_diff[p,:,boot])/n_diff[p,boot];
    }
  }
}

parameters {
  real<lower=-15, upper=15> mu;
  sum_to_zero_vector[N_rate_params] epsilon;
  real<lower=0, upper=1> alpha;
}

transformed parameters {
  vector[N_rate_params] rates_log = mu+epsilon;
  vector[N_rate_params] rates = exp(rates_log);

  matrix[K,K] Q = rep_matrix(0.0, K, K);
  for (i in 1:K) {
    for (j in 1:K) {
      if (rate_map[i,j] > 0) Q[i,j] = rates[rate_map[i,j]];
    }
  }
  for (i in 1:K) {
    Q[i,i] = -sum(Q[i,1:K]);
  }

  array[N_sets,max_N_bins] real cr = rep_array(0.0, N_sets, max_N_bins);
  for (bin in 1:max_N_bins) {
    real width = bin_break_max[bin]-bin_break_min[bin];
    for (g in 1:N_grid) {
      matrix[K,K] T = transition_matrix(Q, (bin_break_min[bin]+(g-0.5)/N_grid*width)/2);
      for (s in 1:N_sets) {
        int p = set_pair[s];
        real p_ab = dot_product(T[pair_a[p],], T[pair_b[p],]);
        real p_ref = ref_weight[s,1]*dot_product(T[pair_a[p],], T[pair_a[p],])+ref_weight[s,2]*dot_product(T[pair_b[p],], T[pair_b[p],])+ref_weight[s,3]*p_ab;
        cr[s,bin] += p_ab/p_ref/N_grid;
      }
    }
  }

  vector[N_sets*N_boot] log_lik;
  for (s in 1:N_sets) {
    int p = set_pair[s];
    int nb = N_bins_pair[p];
    for (boot in 1:N_boot) {
      log_lik[(boot-1)*N_sets+s] = set_weight[s]*pair_loglik(hist_diff[p,1:nb,boot], D_ref[s,boot,1:nb], to_vector(cr[s,1:nb]), n_diff[p,boot], alpha)/N_boot;
    }
  }
}

model {
  mu ~ cauchy(-1, 1);
  epsilon ~ normal(0, rate_dispersion_prior);
  alpha ~ beta(1, 99);
  target += sum(log_lik);
}

generated quantities {
  vector[N_rate_params] rates_log_prior;
  {
    real mu_prior = cauchy_rng(-1, 1);
    while (mu_prior > 15 || mu_prior < -15) {
      mu_prior = cauchy_rng(-1, 1);
    }
    vector[N_rate_params] z = to_vector(normal_rng(rep_vector(0, N_rate_params), rate_dispersion_prior));
    rates_log_prior = mu_prior+z-mean(z);
  }
}
