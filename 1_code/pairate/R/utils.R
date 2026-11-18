## Defaults
DEFAULT_STAN_FILE = 'Fit_RR_CT_pairs.stan'

DEFAULT_HIST_DIST_MIN_MAX = c(0, 10000)
DEFAULT_HIST_BREAK_METHOD = 'Manual'
DEFAULT_FIRST_BIN_NPAIRS = 10
DEFAULT_FIRST_BIN_PMIN = 5E-5
DEFAULT_MAX_N_BINS = 200
DEFAULT_COLLECTION_TIME_DIFF_MAX = 1

DEFAULT_CHAINS = 3
DEFAULT_PARALLEL_CHAINS = 3
DEFAULT_ITER_WARMUP = 1000
DEFAULT_ITER_SAMPLING = 1000
DEFAULT_ADAPT_DELTA = 0.95
DEFAULT_MAX_TREEDEPTH = 15
DEFAULT_REFRESH = 500
DEFAULT_N_GRID = 5

#' Resolve a Stan file path within the installed package
#'
#' @param filename Stan file name (e.g., `Fit_RR_CT_pairs.stan`) or a path.
#' @return Absolute path to the Stan file.
#' @export
stan_file_path = function(filename) {
  if(file.exists(filename)) return(filename)
  return(system.file('stan_models', filename, package = 'pairate'))
}

#' Compile the Stan model
#'
#' The executable is stored next to the Stan file, and rebuilt when missing or older than the Stan file.
#'
#' @param stan_file Path to the `.stan` source file.
#' @param force_recompile Logical; recompile even if an up-to-date binary exists.
#' @return A compiled cmdstanr model object.
#' @export
compile_stan_model = function(stan_file, force_recompile = F) {
  exe_path = file.path(dirname(stan_file), tools::file_path_sans_ext(basename(stan_file)))
  if(!file.exists(exe_path) || file.mtime(stan_file) > file.mtime(exe_path)) force_recompile = T
  return(cmdstanr::cmdstan_model(stan_file, exe_file = exe_path, force_recompile = force_recompile))
}

#' Parse cmdstanr CSV table or fit into posterior draws
#'
#' @param table_or_fit Numeric matrix from cmdstanr CSV (no comments) or a cmdstanr fit.
#' @return Named list of arrays/vectors of posterior draws.
#' @export
read_chains = function(table_or_fit){
  table = table_or_fit
  if(inherits(table_or_fit, 'CmdStanFit')) table = as.matrix(do.call(rbind, lapply(table_or_fit$output_files(), read.csv, comment.char = '#')))
  if(typeof(table) != 'double') {
    message('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_names = gsub('\\[|\\]|,', '.', colnames(table))
  col_names = sub('\\.$', '', gsub('\\.+', '.', col_names))
  col_variables = sapply(col_names, function(x) strsplit(x, '.', fixed = T)[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(sum(!is.na(a)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      tmp = col_names[which(!is.na(a))]
      ndims = 0
      dims = NULL
      empty = F
      while(!empty && ndims < 10){
        l = sapply(tmp, function(x) strsplit(x, ".", fixed = T)[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(!all(is.na(l))) {
          ndims = ndims+1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      a = which(!is.na(match(col_variables, variable_names[i])))
      Chains[[i]] = array(table[,a], dim = c(nchains, dims))
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}

#' Median and 95\% credible interval helper
#'
#' @param v Numeric vector of samples.
#' @return Numeric vector of length 3: median, 2.5\% quantile, 97.5\% quantile.
#' @export
median_and_ci = function(v) {
  return(c(as.numeric(median(v, na.rm = T)),
           as.numeric(quantile(v, probs = 0.025, na.rm = T)),
           as.numeric(quantile(v, probs = 0.975, na.rm = T))))
}

#' Choose a default number of histogram bins to fit
#'
#' Ratio R(h) of the within-deme distribution D(h|same) = (D(h|aa) + D(h|bb)) / 2 to the between-deme
#' distribution D(h|ab), summed over pairs of demes, each histogram divided by its number of pairs at any
#' distance. Fits the bins before the first one where R(h) drops below 1 (from there on the between
#' histogram is no longer depleted and the bins mostly add genealogical noise), or all bins if it never
#' does, with a minimum of `min_bins` bins. The normalisation makes the choice independent of the number
#' of sequences per deme.
#'
#' @param pair_hists Pair histograms (`pair_hists` from `compute_histograms()`).
#' @param min_bins Minimum number of bins.
#' @return Integer number of bins to fit.
#' @keywords internal
choose_bins_to_fit = function(pair_hists, min_bins = 5) {
  same = (pair_hists$within_a/pair_hists$n_within_a+pair_hists$within_b/pair_hists$n_within_b)/2
  ratio = colSums(same)/colSums(pair_hists$between/pair_hists$n_between)
  below = which(is.finite(ratio) & ratio < 1)
  n_bins_fit = length(ratio)
  if(length(below) > 0) n_bins_fit = below[1]-1
  return(min(length(ratio), max(n_bins_fit, min_bins)))
}

#' Build rate labels, index matrix and rate map
#'
#' @param demes Character vector of deme names.
#' @param model_type Either a character shortcut ("ER", "SYM", "ARD") or a
#'   numeric matrix specifying the rate structure (like `ace()`'s model argument):
#'   square, positive integers off-diagonal giving the rate parameter of each
#'   transition (equal integers = equal rates, 0 = no transition).
#' @return List with `model` ('ER', 'SYM', 'ARD' or 'CUSTOM'), `rate_labels`, `rate_index` (character matrix
#'   for printing), `rate_map` (integer matrix used by the Stan model) and `n_rate_params`.
#' @keywords internal
build_rate_index = function(demes, model_type) {
  K = length(demes)
  model = 'CUSTOM'
  if(!is.matrix(model_type)) model = toupper(model_type)
  if(model == 'CUSTOM') {
    rate_map = model_type
    if(nrow(rate_map) != K || ncol(rate_map) != K) stop('model matrix must be ', K, 'x', K, ' (one row/column per deme: ', paste(demes, collapse = ', '), ').')
  }else{
    if(!(model %in% c('ER', 'SYM', 'ARD'))) stop('model must be "ER", "SYM", "ARD", or a numeric matrix.')
    rate_map = matrix(0, K, K)
    if(model == 'ER') rate_map[] = 1
    if(model == 'SYM') {
      deme_pairs = t(combn(K, 2))
      for(p in 1:nrow(deme_pairs)) rate_map[deme_pairs[p,1],deme_pairs[p,2]] = rate_map[deme_pairs[p,2],deme_pairs[p,1]] = p
    }
    if(model == 'ARD') {
      r = 0
      for(i in 1:K) for(j in 1:K) if(i != j) {
        r = r+1
        rate_map[i,j] = r
      }
    }
  }
  diag(rate_map) = 0
  n_rate_params = as.integer(max(rate_map))
  rate_index = matrix('.', K, K, dimnames = list(demes, demes))
  rate_index[rate_map > 0] = rate_map[rate_map > 0]
  sep = '->'
  if(model == 'SYM') sep = '-'
  rate_labels = sapply(1:n_rate_params, function(r) {
    idx = which(rate_map == r, arr.ind = T)
    if(model == 'SYM') idx = idx[idx[,1] < idx[,2],,drop = F]
    return(paste(paste(demes[idx[,1]], demes[idx[,2]], sep = sep), collapse = '/'))
  })
  if(model == 'ER') rate_labels = 'all'
  return(list(model = model,
              rate_labels = rate_labels,
              rate_index = rate_index,
              rate_map = rate_map,
              n_rate_params = n_rate_params))
}

#' Compute WAIC for a fitted model
#'
#' Works from parsed chains (no fit object needed).
#'
#' @param result List returned by `fit_hist_model()`.
#' @return A `loo::waic` object.
#' @export
compute_waic = function(result) {
  return(loo::waic(as.matrix(result$chains$log_lik)))
}

#' Compute LOO-CV for a fitted model
#'
#' Works from parsed chains (no fit object needed).
#'
#' @param result List returned by `fit_hist_model()`.
#' @return A `loo::loo` object.
#' @export
compute_loo = function(result) {
  return(loo::loo(as.matrix(result$chains$log_lik)))
}
