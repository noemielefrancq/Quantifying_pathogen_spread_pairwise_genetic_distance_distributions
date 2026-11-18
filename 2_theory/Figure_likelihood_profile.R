#########################################################################################################################################
## Likelihood profiles on simulated trees
## The log-likelihood below is the same as in Fit_RR_CT_pairs.stan
#########################################################################################################################################
source('1_code/load_pairate.R', chdir = T)

#########################################################################################################################################
## Log-likelihood, as in Fit_RR_CT_pairs.stan
#########################################################################################################################################
transition_matrix = function(Q, t) {
  res = matrix(0, nrow(Q), nrow(Q))
  if(nrow(Q) > 2) {
    res = as.matrix(Matrix::expm(Q*t))
  }else{
    a = Q[1,2]
    b = Q[2,1]
    e = exp(-(a+b)*t)
    res[1,1] = (b+a*e)/(a+b)
    res[1,2] = a*(1-e)/(a+b)
    res[2,1] = b*(1-e)/(a+b)
    res[2,2] = (a+b*e)/(a+b)
  }
  return(res)
}

ldirmult = function(y, alpha) {
  return(lgamma(sum(alpha))-lgamma(sum(alpha)+sum(y))+sum(lgamma(alpha+y)-lgamma(alpha))+lgamma(sum(y)+1)-sum(lgamma(y+1)))
}

pair_loglik = function(y, D_ref, cr, n_diff, alpha) {
  keep = which(D_ref > 0)
  P = D_ref[keep]*cr[keep]
  return(ldirmult(y[keep], P/sum(P)*(1-alpha)/alpha)+dnbinom(sum(y), size = (1-alpha)/alpha, mu = n_diff*sum(P), log = T))
}

## Sum of log_lik in Fit_RR_CT_pairs.stan, for given rates and alpha
loglik = function(rates, alpha, stan_data) {
  K = stan_data$K
  Q = matrix(0, K, K)
  for(i in 1:K) {
    for(j in 1:K) {
      if(stan_data$rate_map[i,j] > 0) Q[i,j] = rates[stan_data$rate_map[i,j]]
    }
  }
  for(i in 1:K) {
    Q[i,i] = -sum(Q[i,1:K])
  }

  cr = matrix(0, stan_data$N_sets, stan_data$max_N_bins)
  for(bin in 1:stan_data$max_N_bins) {
    width = stan_data$bin_break_max[bin]-stan_data$bin_break_min[bin]
    for(g in 1:stan_data$N_grid) {
      T_mat = transition_matrix(Q, (stan_data$bin_break_min[bin]+(g-0.5)/stan_data$N_grid*width)/2)
      for(s in 1:stan_data$N_sets) {
        p = stan_data$set_pair[s]
        a = stan_data$pair_a[p]
        b = stan_data$pair_b[p]
        p_ab = sum(T_mat[a,]*T_mat[b,])
        p_ref = stan_data$ref_weight[s,1]*sum(T_mat[a,]*T_mat[a,])+stan_data$ref_weight[s,2]*sum(T_mat[b,]*T_mat[b,])+stan_data$ref_weight[s,3]*p_ab
        cr[s,bin] = cr[s,bin]+p_ab/p_ref/stan_data$N_grid
      }
    }
  }

  log_lik = rep(0, stan_data$N_sets*stan_data$N_boot)
  for(s in 1:stan_data$N_sets) {
    p = stan_data$set_pair[s]
    nb = stan_data$N_bins_pair[p]
    for(boot in 1:stan_data$N_boot) {
      D_ref = stan_data$ref_weight[s,1]*stan_data$hist_within_a[p,,boot]/stan_data$n_within_a[p,boot]+stan_data$ref_weight[s,2]*stan_data$hist_within_b[p,,boot]/stan_data$n_within_b[p,boot]+stan_data$ref_weight[s,3]*stan_data$hist_diff[p,,boot]/stan_data$n_diff[p,boot]
      log_lik[(boot-1)*stan_data$N_sets+s] = stan_data$set_weight[s]*pair_loglik(stan_data$hist_diff[p,1:nb,boot], D_ref[1:nb], cr[s,1:nb], stan_data$n_diff[p,boot], alpha)/stan_data$N_boot
    }
  }
  return(sum(log_lik))
}

## Profile log-likelihood of the one-rate (ER) model: for each rate, maximum over alpha
profile_likelihood_ER = function(hists, rate_grid) {
  stan_data = make_stan_data(hists, hists$n_bins_to_fit, build_rate_index(hists$demes_in_tree, 'ER')$rate_map)
  log_lik = rep(NA, length(rate_grid))
  alpha_hat = rep(NA, length(rate_grid))
  for(i in 1:length(rate_grid)) {
    opt = optimize(function(alpha) loglik(rate_grid[i], alpha, stan_data), c(1e-10, 1-1e-10), maximum = T)
    log_lik[i] = opt$objective
    alpha_hat[i] = opt$maximum
  }
  return(data.frame(rate = rate_grid,
                    log_lik = log_lik,
                    alpha_hat = alpha_hat))
}
#########################################################################################################################################

#########################################################################################################################################
## Profiles on simulated trees
#########################################################################################################################################
migrationrates = c(0.01, 0.1, 1, 10)
res_migrationrates = list()
for(m in 1:length(migrationrates)) {
  tree_data = treeio::read.beast(paste0('3_simulations_ReMaster/example_trees/Structured_Coalescent_AB_exponential_pop/AB_base/Structured_Coalescent_AB_m_',
                                        migrationrates[m],
                                        '_NeB_1000_rB_0.1_Nb_1000.tree'))
  res_list = list()
  for(i in 1:length(tree_data)) {
    print(i)
    tree_tmp = tree_data[[i]]
    tip_data = as.data.frame(tree_tmp@data)
    tip_data = tip_data[match(1:(tree_tmp@phylo$Nnode+1), tip_data$node),]
    hists = compute_histograms(tree = tree_tmp@phylo,
                               traits = as.character(tip_data$type),
                               times = as.numeric(tip_data$time),
                               hist_dist_min_max = c(0, 10000),
                               hist_break_method = 'Manual',
                               first_bin_npairs = 10,
                               first_bin_pmin = 5E-5,
                               max_n_bins = 200,
                               collection_time_diff_max = NULL,
                               bootstrap = 0)
    res_list[[i]] = profile_likelihood_ER(hists, rate_grid = 10^seq(-4, 3, length.out = 200))
  }
  res_migrationrates[[m]] = res_list
}

# pdf('2_theory/Figure_method/Likelihood_profiles_15102026.pdf', width = 10, height = 3)
par(mfrow = c(1, 4))
for(m in 1:length(migrationrates)) {
  results = res_migrationrates[[m]]
  plot(NULL, xlim = c(0.001, 100), ylim = c(-100, 0), log = 'x', xlab = 'Migration rate', ylab = 'llk', bty = 'n', yaxt = 'n', main = paste0('m = ', migrationrates[m]))
  axis(2, las = 2)
  for(i in 1:length(results)) {
    lines(results[[i]]$rate, results[[i]]$log_lik-max(results[[i]]$log_lik))
    idx = which(results[[i]]$log_lik-max(results[[i]]$log_lik) == 0)
    abline(v = results[[i]]$rate[idx], lty = 3, col = 'grey25')
  }
  abline(v = migrationrates[m], lty = 2, col = 'red', lwd = 1.5)
}
# dev.off()
#########################################################################################################################################
