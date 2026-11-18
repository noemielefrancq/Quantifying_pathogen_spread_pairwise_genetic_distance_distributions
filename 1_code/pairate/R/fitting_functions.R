#' Build Stan data from pair histograms
#'
#' For each pair of demes (a, b): between histogram, within-a and within-b histograms
#' and the numbers of pairs at any distance N(ab), N(aa), N(bb), for the observed data
#' (`bootstrap = 0`) or the first `bootstrap` replicates. The fitted sets are one per pair
#' (`by_category = FALSE`, reference D_ref(d) = N(d|aa) / (4 N(aa)) + N(d|bb) / (4 N(bb)) + N(d|ab) / (2 N(ab))) or two per pair
#' (`by_category = TRUE`, a->b with reference N(d|aa) / (3 N(aa)) + 2 N(d|ab) / (3 N(ab)) and b->a with reference
#' N(d|bb) / (3 N(bb)) + 2 N(d|ab) / (3 N(ab)), each with weight 1/2 in the likelihood).
#'
#' @param hists Output of `compute_histograms()` or `compute_histograms_from_dist()`.
#' @param n_bins_to_fit Number of bins to fit (scalar or one value per pair).
#' @param rate_map Integer K x K matrix of rate indices (see `build_rate_index()`).
#' @param by_category Logical; fit by deme (TRUE) or by pair (FALSE).
#' @param N_grid Integration grid size per bin.
#' @param bootstrap Number of bootstrap replicates to use (0 for the observed data).
#' @param rate_dispersion_prior SD of the normal prior on the log-rate deviations.
#' @return List formatted for `Fit_RR_CT_pairs.stan`.
#' @keywords internal
make_stan_data = function(hists, n_bins_to_fit, rate_map, by_category = F, N_grid = DEFAULT_N_GRID, bootstrap = 0, rate_dispersion_prior = 1) {
  if(bootstrap > 0) {
    if(is.null(hists$pair_hists_boot) || dim(hists$pair_hists_boot$between)[3] < bootstrap) {
      stop('Not enough bootstrap replicates in the histograms: compute them with bootstrap >= ', bootstrap, '.')
    }
    pair_hists = hists$pair_hists_boot
  }else{
    pair_hists = hists$pair_hists
  }
  n_boot = max(bootstrap, 1)
  N_pairs = nrow(hists$pair_mat)
  sets = make_sets(hists, by_category)
  n_bins = length(hists$common_breaks)-1
  n_bins_pair = pmin(rep(n_bins_to_fit, length.out = N_pairs), n_bins)
  max_n_bins_fit = max(n_bins_pair)
  hist_diff = array(pair_hists$between, c(N_pairs, n_bins, n_boot))[,1:max_n_bins_fit,,drop = F]
  hist_within_a = array(pair_hists$within_a, c(N_pairs, n_bins, n_boot))[,1:max_n_bins_fit,,drop = F]
  hist_within_b = array(pair_hists$within_b, c(N_pairs, n_bins, n_boot))[,1:max_n_bins_fit,,drop = F]
  return(list(
    K = length(hists$demes_in_tree),
    N_pairs = N_pairs,
    pair_a = as.array(match(hists$pair_mat[,1], hists$demes_in_tree)),
    pair_b = as.array(match(hists$pair_mat[,2], hists$demes_in_tree)),
    N_rate_params = max(rate_map),
    rate_map = rate_map,
    max_N_bins = max_n_bins_fit,
    N_bins_pair = as.array(n_bins_pair),
    bin_break_min = hists$common_breaks[1:max_n_bins_fit],
    bin_break_max = hists$common_breaks[2:(max_n_bins_fit+1)],
    N_grid = N_grid,
    N_boot = n_boot,
    rate_dispersion_prior = rate_dispersion_prior,
    hist_diff = hist_diff,
    hist_within_a = hist_within_a,
    hist_within_b = hist_within_b,
    n_diff = array(pair_hists$n_between, c(N_pairs, n_boot)),
    n_within_a = array(pair_hists$n_within_a, c(N_pairs, n_boot)),
    n_within_b = array(pair_hists$n_within_b, c(N_pairs, n_boot)),
    N_sets = nrow(sets),
    set_pair = as.array(sets$pair),
    ref_weight = cbind(sets$pi_a, sets$pi_b, sets$pi_ab),
    set_weight = as.array(sets$weight)
  ))
}

#' Sets of histograms fitted: one per pair (by pair) or two per pair (by deme)
#'
#' @param hists Output of `compute_histograms()`.
#' @param by_category Logical; fit by deme (TRUE) or by pair (FALSE).
#' @return data.frame with the pair, the weights of the within-a, within-b and between histograms
#'   in the reference (`pi_a`, `pi_b`, `pi_ab`), the likelihood weight and a label of each set.
#' @keywords internal
make_sets = function(hists, by_category) {
  pair_mat = hists$pair_mat
  N_pairs = nrow(pair_mat)
  if(!by_category) return(data.frame(pair = 1:N_pairs,
                                      pi_a = 1/4,
                                      pi_b = 1/4,
                                      pi_ab = 1/2,
                                      weight = 1,
                                      label = hists$pair_labels))
  return(data.frame(pair = rep(1:N_pairs, each = 2),
                    pi_a = rep(c(1/3, 0), N_pairs),
                    pi_b = rep(c(0, 1/3), N_pairs),
                    pi_ab = 2/3,
                    weight = 0.5,
                    label = as.vector(rbind(paste(pair_mat[,1], pair_mat[,2], sep = '->'), paste(pair_mat[,2], pair_mat[,1], sep = '->')))))
}

#' Fit migration rates from pairwise distance histograms
#'
#' By pair (default), for each pair of demes (a, b), the between-deme histogram of pairwise distances
#' is fitted against a reference distribution P(d): the distances of random pairs in a population where demes
#' a and b are equally represented, i.e. within-a, within-b and between pairs in proportions 1/4, 1/4 and 1/2
#' (each histogram divided by its number of pairs at any distance, so it does not depend on sampling):
#' D_ref(d) = N(d|aa) / (4 N(aa)) + N(d|bb) / (4 N(bb)) + N(d|ab) / (2 N(ab)). The model distribution of
#' between pairs is P(d|ab) = D_ref(d) cr_ab(d), with cr_ab(d) = P(A_d|ab) / P(A_d) the colocation
#' ratio averaged over each bin: the probability that the ancestors of a between pair are in the same deme at their MRCA (d/2),
#' relative to P(A_d) = P(A_d|aa) / 4 + P(A_d|bb) / 4 + P(A_d|ab) / 2 with the same weights, from the
#' continuous-time Markov chain with rate matrix Q. By deme (`by_category = TRUE`), each pair gives two
#' sets with the same population restricted to the set's pair types: a->b with reference
#' N(d|aa) / (3 N(aa)) + 2 N(d|ab) / (3 N(ab)) and P(A_d) = P(A_d|aa) / 3 + 2 P(A_d|ab) / 3,
#' and b->a likewise, each weighted 1/2 in the likelihood. The between histogram is modelled with a
#' Dirichlet-multinomial on P(d|ab) normalised over the fitted bins (shape) and a negative binomial with
#' mean N(ab) times the sum of P(d|ab) over the fitted bins (total). Bins with an empty reference are not
#' fitted. Bootstrap replicates are combined by averaging their log-likelihoods. Log-rates are mu + epsilon,
#' with epsilon ~ N(0, rate_dispersion_prior) summing to 0.
#'
#' @param tree A `phylo` object (not needed if `precomputed_hists` is given).
#' @param traits Vector/factor of tip traits (same order as `tree$tip.label`).
#' @param times Optional numeric vector of tip collection times.
#' @param demes Optional vector of deme labels to keep.
#' @param model Rate structure: 'ER' (one rate), 'SYM' (one rate per pair of demes), 'ARD' (one rate per
#'   ordered pair), or a K x K integer matrix as the `model` argument of `ace()` in `ape` (equal integers
#'   mean equal rates, 0 no transition).
#' @param by_category Logical; fit by deme (TRUE, sets a->b and b->a) or by pair (FALSE). Fitting by deme is
#'   biased when demes differ in demography (population size, growth), as each set's reference holds only one deme.
#' @param rate_dispersion_prior SD of the normal prior on the log-rate deviations from their mean
#'   (e.g. 100 for an essentially flat prior).
#' @param hist_dist_min_max Optional numeric length-2 vector for distance bounds.
#' @param hist_break_method Break method: 'Sturges', 'Scott', 'FD', or 'Manual'.
#' @param first_bin_npairs For 'Manual', number of pairs in first bin or a vector of breaks. Bins are closed on the left, [a, b) (the last one also includes its upper break).
#' @param first_bin_pmin For 'Manual', minimum probability in first bin.
#' @param max_n_bins Maximum number of bins.
#' @param collection_time_diff_max Max allowed time difference between tips (only used if `times` are given).
#' @param bootstrap Integer; number of bootstrap replicates to fit (0 for the observed data only). By default
#'   (NULL), none when histograms are computed from the tree, and all the replicates of `precomputed_hists`.
#'   With `precomputed_hists`, a supplied value replaces theirs and cannot exceed the replicates they contain.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of parallel chains.
#' @param iter_warmup Warmup iterations.
#' @param iter_sampling Sampling iterations.
#' @param adapt_delta Stan adapt_delta.
#' @param max_treedepth Stan max_treedepth.
#' @param refresh Stan refresh rate.
#' @param verbose Logical; print progress messages.
#' @param N_grid Integration grid size per bin.
#' @param precomputed_hists Optional output of `compute_histograms()` or `compute_histograms_from_dist()`.
#' @param keep_raw_data Logical; if FALSE, drop the raw data (distance matrix) from results.
#' @param keep_stan_data Logical; if FALSE, omit Stan data from results (needed for plots).
#' @param compiled_model Optional precompiled cmdstanr model.
#' @param force_recompile Logical; force recompiling the Stan model.
#' @param stan_data_override Optional named list replacing elements of the Stan data.
#' @param n_bins_to_fit Optional number of bins to fit (scalar or one value per pair); see `compute_histograms()`.
#' @return A list with the fit object, rate estimates (`estimates`: median and 95\% credible interval of each rate),
#'   MCMC diagnostics (`diagnostics`: posterior summary with R-hat and ESS of `mu`, `rates` and `alpha`, and
#'   sampler diagnostics; kept in the result so they remain available when the fit is saved and read elsewhere),
#'   parsed chains (`rates`, `rates_log`, `cr` (draws x sets x bins), `alpha`, ...),
#'   rate labels and matrices, Stan data and metadata.
#' @export
fit_hist_model = function(tree = NULL,
                          traits = NULL,
                          times = NULL,
                          demes = NULL,
                          model = 'SYM',
                          by_category = F,
                          rate_dispersion_prior = 1,
                          hist_dist_min_max = DEFAULT_HIST_DIST_MIN_MAX,
                          hist_break_method = DEFAULT_HIST_BREAK_METHOD,
                          first_bin_npairs = DEFAULT_FIRST_BIN_NPAIRS,
                          first_bin_pmin = DEFAULT_FIRST_BIN_PMIN,
                          max_n_bins = DEFAULT_MAX_N_BINS,
                          collection_time_diff_max = DEFAULT_COLLECTION_TIME_DIFF_MAX,
                          bootstrap = NULL,
                          chains = DEFAULT_CHAINS,
                          parallel_chains = DEFAULT_PARALLEL_CHAINS,
                          iter_warmup = DEFAULT_ITER_WARMUP,
                          iter_sampling = DEFAULT_ITER_SAMPLING,
                          adapt_delta = DEFAULT_ADAPT_DELTA,
                          max_treedepth = DEFAULT_MAX_TREEDEPTH,
                          refresh = DEFAULT_REFRESH,
                          verbose = T,
                          N_grid = DEFAULT_N_GRID,
                          precomputed_hists = NULL,
                          keep_raw_data = T,
                          keep_stan_data = T,
                          compiled_model = NULL,
                          force_recompile = F,
                          stan_data_override = NULL,
                          n_bins_to_fit = NULL) {
  if(is.null(precomputed_hists)) {
    if(is.null(tree) || is.null(traits)) stop('Either supply tree + traits, or precomputed_hists.')
    if(is.null(bootstrap)) bootstrap = 0
    if(verbose) cat('Computing histograms...\n')
    precomputed_hists = compute_histograms(tree = tree, traits = traits, times = times, demes = demes,
                                           hist_dist_min_max = hist_dist_min_max, hist_break_method = hist_break_method,
                                           first_bin_npairs = first_bin_npairs, first_bin_pmin = first_bin_pmin,
                                           max_n_bins = max_n_bins, collection_time_diff_max = collection_time_diff_max,
                                           bootstrap = bootstrap, n_bins_to_fit = n_bins_to_fit,
                                           verbose = verbose)
  }
  hists = precomputed_hists
  if(is.null(bootstrap)) bootstrap = hists$bootstrap
  if(is.null(n_bins_to_fit)) n_bins_to_fit = hists$n_bins_to_fit
  if(length(unique(n_bins_to_fit)) == 1) n_bins_to_fit = n_bins_to_fit[1]

  rate_info = build_rate_index(hists$demes_in_tree, model)
  stan_data = make_stan_data(hists, n_bins_to_fit, rate_info$rate_map, by_category, N_grid, bootstrap, rate_dispersion_prior)
  for(nm in names(stan_data_override)) stan_data[[nm]] = stan_data_override[[nm]]

  stan_model = compiled_model
  if(is.null(stan_model)) stan_model = compile_stan_model(stan_file_path(DEFAULT_STAN_FILE), force_recompile = force_recompile)
  init = lapply(1:chains, function(x) list(mu = rnorm(1, -3, 1),
                                            epsilon = array(0, stan_data$N_rate_params),
                                            alpha = 0.01))
  fit = stan_model$sample(data = stan_data, chains = chains, parallel_chains = parallel_chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                          adapt_delta = adapt_delta, max_treedepth = max_treedepth, refresh = refresh, init = init)

  result = list(
    all = list(hists = hists$hists_all,
               data = hists$data_all),
    fit = fit,
    chains = tryCatch(read_chains(fit), error = function(e) NULL),
    pair_labels = hists$pair_labels,
    set_labels = make_sets(hists, by_category)$label,
    N_bins_to_fit = n_bins_to_fit,
    rate_labels = rate_info$rate_labels,
    rate_index = rate_info$rate_index,
    rate_map = rate_info$rate_map,
    stan_data = stan_data,
    demes = hists$demes_in_tree,
    model = rate_info$model,
    bootstrap_n = bootstrap,
    by_category = by_category
  )
  if(!keep_raw_data) result$all$data = NULL
  if(!keep_stan_data) result$stan_data = NULL
  if(!is.null(result$chains)) {
    result$estimates = rate_estimates(result)
    result$diagnostics = list(summary = as.data.frame(fit$summary(c('mu', 'rates', 'alpha'))),
                              sampler = fit$diagnostic_summary(quiet = T))
    result$waic = compute_waic(result)
    result$loo = compute_loo(result)
  }
  if(verbose) {
    if(is.null(result$chains)) {
      cat('No chains finished successfully, see result$fit$output(1).\n')
    }else{
      cat('Computing the diagnostics..\n')
      print_diagnostics(result)
      cat('\nEstimated rates:\n')
      print_estimates(result)
      cat('\n')
    }
  }
  return(result)
}

#' Rate estimates with credible intervals
#'
#' @param result List returned by `fit_hist_model()`.
#' @return data.frame with the rate index, demes, posterior median and 95\% credible interval of each rate.
#' @keywords internal
rate_estimates = function(result) {
  rates_mci = t(apply(as.matrix(result$chains$rates), 2, median_and_ci))
  return(data.frame(rate_index = 1:nrow(rates_mci),
                    demes = result$rate_labels,
                    median = rates_mci[,1],
                    lo_ci = rates_mci[,2],
                    hi_ci = rates_mci[,3]))
}

#' Print rate estimates with credible intervals
#'
#' @param result List returned by `fit_hist_model()`.
#' @return Invisibly returns a data.frame of estimates.
#' @export
print_estimates = function(result) {
  if(is.null(result$chains) || is.null(result$chains$rates)) {
    cat('No rates found in result.\n')
    return(invisible(NULL))
  }
  out = rate_estimates(result)
  cat('Rate index matrix:\n')
  print(result$rate_index, quote = F)
  cat('\nParameter estimates:\n')
  print(out, row.names = F)
  return(invisible(out))
}

#' Print basic MCMC diagnostics
#'
#' Uses the diagnostics stored in the result, so it also works on a fit saved and read on another machine.
#'
#' @param result List returned by `fit_hist_model()`.
#' @param rhat_max Maximum acceptable R-hat.
#' @param ess_min Minimum acceptable ESS.
#' @return Invisibly returns the posterior summary.
#' @export
print_diagnostics = function(result, rhat_max = 1.01, ess_min = 100) {
  fit_summary = result$diagnostics$summary
  if(is.null(fit_summary)) fit_summary = as.data.frame(result$fit$summary(c('mu', 'rates', 'alpha')))
  max_rhat = max(fit_summary$rhat, na.rm = T)
  min_ess_bulk = min(fit_summary$ess_bulk, na.rm = T)
  min_ess_tail = min(fit_summary$ess_tail, na.rm = T)
  cat(sprintf('Diagnostics: max R-hat = %.3f; min ESS bulk = %.1f; min ESS tail = %.1f', max_rhat, min_ess_bulk, min_ess_tail))
  if(!is.null(result$diagnostics$sampler)) cat('; divergent transitions =', sum(result$diagnostics$sampler$num_divergent))
  cat('\n')

  failed = NULL
  if(max_rhat > rhat_max) failed = c(failed, 'R-hat')
  if(min_ess_bulk < ess_min) failed = c(failed, 'ESS bulk')
  if(min_ess_tail < ess_min) failed = c(failed, 'ESS tail')
  if(length(failed) > 0) message('Some diagnostics failed: ', paste(failed, collapse = ', '), '\nPlease increase chain length and/or double check parameters')
  return(invisible(fit_summary))
}
