# Remove histogram bins with 0 counts by recomputing with pruned breaks
remove_empty_bins = function(hist_obj, data_values) {
  if(any(hist_obj$counts == 0)) {
    idx = which(hist_obj$counts > 0)
    new_breaks = hist_obj$breaks[sort(unique(c(idx, idx+1)))]
    hist_obj = hist(data_values, breaks = new_breaks, right = F, plot = F)
  }
  return(hist_obj)
}

#' Pre-process a phylogenetic tree for histogram analysis
#'
#' Builds a tip metadata data.frame and a sparse pairwise distance matrix
#' (upper triangle only) from a phylogeny. Optionally filters by collection
#' time window, subsamples tips, and removes pairs outside a distance range.
#'
#' @param phy A `phylo` object.
#' @param traits Vector/factor of tip traits (same order as `phy$tip.label`).
#' @param collection_times Optional vector of tip collection times.
#' @param hist_dist_min_max Optional numeric length-2 vector for distance bounds.
#' @param subsampling_size Optional integer; number of tips to sample.
#' @param collection_times_to_sample Optional numeric length-2 vector of time window.
#' @param collection_time_diff_max Optional numeric; max allowed time difference between tips.
#' @return A list with:
#' 1) data.frame of selected tips (id, traits, time if provided)
#' 2) sparse upper-triangular matrix of pairwise distances.
#' @export
pre_process_tree = function(phy,
                            traits,
                            collection_times = NULL,

                            hist_dist_min_max = NULL,
                            subsampling_size = NULL,
                            collection_times_to_sample = NULL,
                            collection_time_diff_max = NULL) {
  tips = phy$tip.label
  ntips = length(tips)
  if(length(traits) != ntips) stop("Found length(traits) != length(phy$tip.label), but there should be one trait per tip. ")

  ###########################################
  ## Step 1: build a df for selected isolates
  ###########################################
  df = data.frame(id = tips,
                  traits = traits)
  if(!is.null(collection_times)){
    if(length(collection_times) != ntips) stop("Found length(collection_times) != length(phy$tip.label), but there should be one trait per tip. ")
    df$time = collection_times
  }

  # Sampling
  weights = rep(1, ntips)
  # If specific time window
  if(!is.null(collection_times_to_sample)){
    weights[which(df$time < collection_times_to_sample[1] | df$time > collection_times_to_sample[2])] = 0
  }
  # If subsampling size
  if(!is.null(subsampling_size)){
    if(length(which(weights != 0)) > subsampling_size){
      sample_ids = sample(x = df$id, size = subsampling_size, prob = weights, replace = F)
    }else{
      sample_ids = df$id
      sample_ids = sample_ids[which(weights != 0)]
    }
  }else{ # If no subsampling size
    sample_ids = df$id
    sample_ids = sample_ids[which(weights != 0)]
  }
  df = df[df$id %in% sample_ids,]

  ###########################################
  ## Step 2: for the selected isolates, build a distance matrix
  ###########################################
  # Distance matrix from phylogeny
  phy = ape::drop.tip(phy, tip = setdiff(phy$tip.label, df$id))
  dist_mat = ape::cophenetic.phylo(phy)
  dist_mat = dist_mat[df$id,df$id,drop = F]
  dist_mat[which(dist_mat == 0)] = 1E-20
  if(!is.null(hist_dist_min_max)) {
    idx = which(dist_mat < min(hist_dist_min_max) | dist_mat > max(hist_dist_min_max), arr.ind = T)
    dist_mat[idx] = 0
    if(length(which(dist_mat == 0))/(ncol(dist_mat)^2) > 0.99) warning('With your hist_dist_min_max bounds, more than 99% of pairs are removed. Continuing computation anyways.')
  }
  dist_mat[lower.tri(dist_mat)] = 0
  diag(dist_mat) = 0
  dist_mat = Matrix::Matrix(dist_mat, sparse = T)

  # Collection time matrix - if collection_times_to_sample exist
  if(!is.null(collection_time_diff_max)){
    collection_time_mat = Matrix::Matrix(abs(outer(df$time[match(dist_mat@Dimnames[[1]], df$id)], df$time[match(dist_mat@Dimnames[[2]], df$id)], FUN = '-')), sparse = T)
    keep = Matrix::Matrix(collection_time_mat <= collection_time_diff_max)
    dist_mat = Matrix::Matrix(dist_mat*keep, sparse = T)
    if(length(dist_mat@x) == 0) stop('With your collection_times_to_sample bounds, there are no pairs left. Check your bounds.')
  }

  # Store pre processed data
  pre_processed_data = list()
  pre_processed_data[[1]] = df
  pre_processed_data[[2]] = dist_mat

  return(pre_processed_data)
}

#' Compute histogram breaks
#'
#' @param distances Pairwise distances (all pairs).
#' @param distances_between Pairwise distances of between-deme pairs, used by 'Manual' with a number of pairs.
#' @param hist_break_method Break method: 'Sturges', 'Scott', 'FD', or 'Manual'.
#' @param first_bin_npairs For 'Manual', number of pairs in first bin or a vector of breaks. Bins are closed on the left, [a, b) (the last one also includes its upper break).
#' @param first_bin_pmin For 'Manual', minimum probability in first bin.
#' @param max_n_bins Maximum number of bins.
#' @param verbose Logical; print progress messages.
#' @return Numeric vector of breaks.
#' @keywords internal
compute_breaks = function(distances, distances_between, hist_break_method, first_bin_npairs, first_bin_pmin, max_n_bins, verbose = T) {
  if(hist_break_method %in% c('Sturges', 'Scott', 'FD')) {
    if(verbose) message('Using ', hist_break_method, ' histogram breaks.')
    breaks = hist(distances, breaks = hist_break_method, plot = F)$breaks
  }else if(hist_break_method == 'Manual' && is.numeric(first_bin_npairs) && length(first_bin_npairs) > 1) {
    if(verbose) message('Using first_bin_npairs as histogram breaks.')
    breaks = as.numeric(first_bin_npairs)
  }else if(hist_break_method == 'Manual' && is.numeric(first_bin_npairs)) {
    if(verbose) message('Default histogram breaks: equal-width bins, with ', first_bin_npairs, ' between-deme pairs in the first bin (first_bin_npairs).')
    breaks = NULL
    distances_between = sort(distances_between)
    if(length(distances_between) > first_bin_npairs) {
      if(first_bin_npairs/length(distances_between) < first_bin_pmin) {
        bin_width = distances_between[ceiling(length(distances_between)*first_bin_pmin)]
      }else{
        bin_width = distances_between[first_bin_npairs]
      }
      breaks = remove_empty_bins(hist(distances, breaks = seq(floor(min(distances)), max(distances)+bin_width, bin_width), right = F, plot = F), distances)$breaks
    }
    if(length(breaks) <= 2) {
      if(verbose) message('Not enough pairs for first_bin_npairs = ', first_bin_npairs, ', using Scott histogram breaks instead.')
      breaks = remove_empty_bins(hist(distances, breaks = 'Scott', right = F, plot = F), distances)$breaks
    }
  }else{
    stop('hist_break_method not recognised, or Manual breaks not in the right format.')
  }
  if(!is.null(max_n_bins) && length(breaks) > max_n_bins+1) breaks = breaks[1:(max_n_bins+1)]
  return(breaks)
}

#' Count pairs of sequences within and between demes
#'
#' Counts all pairs (at any distance) sampled within `collection_time_diff_max` of each
#' other, or all pairs when no times are given. Sequences can be counted several times
#' (bootstrap multiplicities `m`).
#'
#' @param deme Integer deme index of each sequence.
#' @param times Collection times, or NULL.
#' @param mult Multiplicity of each sequence.
#' @param K Number of demes.
#' @param deme_pairs Two-column matrix of deme pairs.
#' @param collection_time_diff_max Max allowed time difference between sequences.
#' @return List with `within` (one value per deme) and `between` (one value per pair).
#' @keywords internal
count_pairs = function(deme, times, mult, K, deme_pairs, collection_time_diff_max) {
  if(is.null(times)) {
    n_seq = sapply(1:K, function(k) sum(mult[deme == k]))
    n_seq2 = sapply(1:K, function(k) sum(mult[deme == k]^2))
    return(list(within = (n_seq^2-n_seq2)/2,
                between = n_seq[deme_pairs[,1]]*n_seq[deme_pairs[,2]]))
  }
  ## Number of pairs (one time in ta, one in tb, both sorted) sampled within collection_time_diff_max of each other
  n_close = function(ta, tb) return(sum(findInterval(ta+collection_time_diff_max, tb)-findInterval(ta-collection_time_diff_max, tb, left.open = T)))
  times_by_deme = lapply(1:K, function(k) sort(rep(times[deme == k], mult[deme == k])))
  within = sapply(1:K, function(k) (n_close(times_by_deme[[k]], times_by_deme[[k]])-sum(mult[deme == k]^2))/2)
  between = sapply(1:nrow(deme_pairs), function(p) n_close(times_by_deme[[deme_pairs[p,1]]], times_by_deme[[deme_pairs[p,2]]]))
  return(list(within = within,
              between = between))
}

# Bootstrap multiplicities: each deme resampled with replacement at its own size
resample_demes = function(deme, K) {
  mult = numeric(length(deme))
  for(k in 1:K) {
    idx = which(deme == k)
    mult[idx] = tabulate(sample(length(idx), replace = T), length(idx))
  }
  return(mult)
}

#' Compute within- and between-deme histograms for each pair of demes
#'
#' Shared by `compute_histograms()` and `compute_histograms_from_dist()`. For each pair of
#' demes (a, b), returns the within-a, within-b and between-ab histograms, and the number of
#' pairs of each type (all distances). Bootstrap replicates resample the sequences of each deme
#' at its own size.
#'
#' @param data_all List with a data.frame of sequences (id, traits, time) and a sparse
#'   upper-triangular matrix of pairwise distances.
#' @param demes Character vector of demes.
#' @param hist_break_method,first_bin_npairs,first_bin_pmin,max_n_bins Histogram breaks options.
#' @param collection_time_diff_max Max allowed time difference between sequences (NULL if no times).
#' @param bootstrap Number of bootstrap replicates (0 for none).
#' @param n_bins_to_fit Optional number of bins to fit (scalar or one value per pair). By default the bins before
#'   the first one where the within/between ratio drops below 1 (at least 5; see `choose_bins_to_fit()`).
#' @param verbose Logical; print progress messages.
#' @return Histogram object (see `compute_histograms()`).
#' @keywords internal
compute_histograms_helper = function(data_all, demes,
                                     hist_break_method, first_bin_npairs, first_bin_pmin, max_n_bins,
                                     collection_time_diff_max, bootstrap,
                                     n_bins_to_fit = NULL, verbose = T) {
  df = data_all[[1]]
  K = length(demes)
  if(K < 2) stop('Need at least 2 demes.')

  ###########################################
  ## Step 1: histograms for each pair
  ###########################################
  deme_pairs = t(combn(K, 2))
  N_pairs = nrow(deme_pairs)
  cat_index = diag(1:K)
  for(p in 1:N_pairs) cat_index[deme_pairs[p,1],deme_pairs[p,2]] = cat_index[deme_pairs[p,2],deme_pairs[p,1]] = K+p
  deme = match(as.character(df$traits), demes)
  times = df$time
  if(is.null(collection_time_diff_max)) times = NULL

  pairs = Matrix::summary(data_all[[2]])
  if(nrow(pairs) == 0) stop('No pairs of sequences left.')
  pairs$cat = cat_index[cbind(deme[pairs$i], deme[pairs$j])]
  breaks = compute_breaks(pairs$x, pairs$x[pairs$cat > K], hist_break_method, first_bin_npairs, first_bin_pmin, max_n_bins, verbose)
  n_bins = length(breaks)-1
  pairs$bin = cut(pairs$x, breaks = breaks, include.lowest = T, right = F, labels = F)
  pairs = pairs[!is.na(pairs$bin),]
  pairs$cat = factor(pairs$cat, levels = 1:(K+N_pairs))
  pairs$bin = factor(pairs$bin, levels = 1:n_bins)

  ## Histograms (categories x bins) and numbers of pairs, for sequence multiplicities m, in pair format
  count_hists = function(sub, mult) {
    hist_counts = tapply(mult[sub$i]*mult[sub$j], sub[,c('cat', 'bin')], sum, default = 0)
    N = count_pairs(deme, times, mult, K, deme_pairs, collection_time_diff_max)
    return(list(within_a = hist_counts[deme_pairs[,1],,drop = F],
                within_b = hist_counts[deme_pairs[,2],,drop = F],
                between = hist_counts[K+1:N_pairs,,drop = F],
                n_within_a = N$within[deme_pairs[,1]],
                n_within_b = N$within[deme_pairs[,2]],
                n_between = N$between))
  }
  pair_hists = count_hists(pairs, rep(1, nrow(df)))

  counts_all = unclass(table(pairs$cat, pairs$bin))
  hist_data = data.frame(distance = (breaks[-1]+breaks[-length(breaks)])/2,
                         same = colSums(counts_all[1:K,,drop = F]),
                         diff = colSums(counts_all[K+1:N_pairs,,drop = F]))
  hist_data$total = hist_data$same+hist_data$diff
  rownames(hist_data) = NULL
  hists_all = list(hist_data = hist_data)

  if(is.null(n_bins_to_fit)) n_bins_to_fit = choose_bins_to_fit(pair_hists)
  if(length(unique(n_bins_to_fit)) == 1) n_bins_to_fit = n_bins_to_fit[1]
  if(!(length(n_bins_to_fit) %in% c(1, N_pairs))) stop('n_bins_to_fit must be a scalar or a vector of length N_pairs (', N_pairs, ').')
  n_bins_to_fit = pmin(as.integer(n_bins_to_fit), n_bins)

  ###########################################
  ## Step 2: bootstrap 
  ###########################################
  pair_hists_boot = NULL
  if(bootstrap > 0) {
    if(verbose) message('Computing ', bootstrap, ' bootstrap replicates.')
    pair_hists_boot = list(within_a = array(0, c(N_pairs, n_bins, bootstrap)),
                           within_b = array(0, c(N_pairs, n_bins, bootstrap)),
                           between = array(0, c(N_pairs, n_bins, bootstrap)),
                           n_within_a = matrix(0, N_pairs, bootstrap),
                           n_within_b = matrix(0, N_pairs, bootstrap),
                           n_between = matrix(0, N_pairs, bootstrap))
    for(r in 1:bootstrap) {
      boot_hists = count_hists(pairs, resample_demes(deme, K))
      pair_hists_boot$within_a[,,r] = boot_hists$within_a
      pair_hists_boot$within_b[,,r] = boot_hists$within_b
      pair_hists_boot$between[,,r] = boot_hists$between
      pair_hists_boot$n_within_a[,r] = boot_hists$n_within_a
      pair_hists_boot$n_within_b[,r] = boot_hists$n_within_b
      pair_hists_boot$n_between[,r] = boot_hists$n_between
    }
  }

  return(list(data_all = data_all,
              hists_all = hists_all,
              common_breaks = breaks,
              n_bins_to_fit = n_bins_to_fit,
              demes_in_tree = demes,
              pair_mat = matrix(demes[deme_pairs], ncol = 2),
              pair_labels = paste(demes[deme_pairs[,1]], demes[deme_pairs[,2]], sep = '-'),
              bootstrap = bootstrap,
              pair_hists = pair_hists,
              pair_hists_boot = pair_hists_boot))
}

#' Compute pairwise distance histograms from a phylogeny
#'
#' Computes, for each pair of demes (a, b), the histograms of pairwise distances within a,
#' within b and between a and b, together with the number of pairs of each type. These are
#' the data fitted by `fit_hist_model()`. Optional bootstrap replicates resample the sequences
#' of each deme at its own size.
#'
#' @param tree A `phylo` object.
#' @param traits Vector/factor of tip traits (same order as `tree$tip.label`).
#' @param times Optional numeric vector of tip collection times.
#' @param demes Optional vector of deme labels to keep.
#' @param hist_dist_min_max Optional numeric length-2 vector for distance bounds.
#' @param hist_break_method Break method: 'Sturges', 'Scott', 'FD', or 'Manual'.
#' @param first_bin_npairs For 'Manual', number of pairs in first bin or a vector of breaks. Bins are closed on the left, [a, b) (the last one also includes its upper break).
#' @param first_bin_pmin For 'Manual', minimum probability in first bin.
#' @param max_n_bins Maximum number of bins.
#' @param collection_time_diff_max Max allowed time difference between tips (only used if `times` are given).
#' @param bootstrap Integer; number of bootstrap replicates (0 for none).
#' @param n_bins_to_fit Optional number of bins to fit (scalar or one value per pair). By default the bins before
#'   the first one where the within/between ratio drops below 1 (at least 5; see `choose_bins_to_fit()`).
#' @param plot_ratio Logical; plot the observed ratios.
#' @param plot_ratio_by_deme Logical; if TRUE, plot ratios per deme pair.
#' @param verbose Logical; print progress messages.
#' @return A list with `data_all`, `hists_all`, `common_breaks`, `n_bins_to_fit`, `demes_in_tree`,
#'   `pair_mat`, `pair_labels`, `bootstrap`, `pair_hists` (within_a, within_b,
#'   between: pairs x bins; n_within_a, n_within_b, n_between: number of pairs) and
#'   `pair_hists_boot` (same with a last dimension for replicates, NULL without bootstrap).
#' @export
compute_histograms = function(tree,
                              traits,
                              times = NULL,
                              demes = NULL,
                              hist_dist_min_max = DEFAULT_HIST_DIST_MIN_MAX,
                              hist_break_method = DEFAULT_HIST_BREAK_METHOD,
                              first_bin_npairs = DEFAULT_FIRST_BIN_NPAIRS,
                              first_bin_pmin = DEFAULT_FIRST_BIN_PMIN,
                              max_n_bins = DEFAULT_MAX_N_BINS,
                              collection_time_diff_max = DEFAULT_COLLECTION_TIME_DIFF_MAX,
                              bootstrap = 0,
                              n_bins_to_fit = NULL,
                              plot_ratio = F,
                              plot_ratio_by_deme = F,
                              verbose = T) {
  traits = as.character(traits)
  if(!is.null(demes)) {
    keep = which(traits %in% demes)
    tree = ape::drop.tip(tree, tree$tip.label[-keep])
    traits = traits[keep]
    if(!is.null(times)) times = times[keep]
  }
  demes_in_tree = sort(unique(traits))

  if(is.null(times)) collection_time_diff_max = NULL
  data_all = pre_process_tree(
    phy = tree,
    traits = as.factor(traits),
    collection_times = times,
    hist_dist_min_max = hist_dist_min_max,
    collection_time_diff_max = collection_time_diff_max
  )
  hists = compute_histograms_helper(data_all, demes_in_tree, hist_break_method, first_bin_npairs, first_bin_pmin, max_n_bins,
                                    collection_time_diff_max, bootstrap, n_bins_to_fit, verbose)
  if(plot_ratio) plot_ratios_data(hists, by_deme = plot_ratio_by_deme)
  return(hists)
}

#' Compute pairwise distance histograms from a precomputed distance matrix
#'
#' Same as `compute_histograms()`, from a matrix of pairwise distances instead of a tree (e.g. for
#' very large datasets, see `pairsnp`). The output can be passed to `fit_hist_model()`.
#'
#' @param dist_mat Square numeric matrix (or sparse `Matrix`) of pairwise distances, with sequence IDs
#'   as row and column names. Only the upper triangle is used.
#' @param traits Named character/factor vector of deme labels, names matching the row names of `dist_mat`.
#' @param times Optional named numeric vector of collection times, names matching the row names of `dist_mat`.
#' @inheritParams compute_histograms
#' @return Same as `compute_histograms()`.
#' @export
compute_histograms_from_dist = function(
    dist_mat,
    traits,
    times = NULL,
    demes = NULL,
    hist_dist_min_max = DEFAULT_HIST_DIST_MIN_MAX,
    hist_break_method = DEFAULT_HIST_BREAK_METHOD,
    first_bin_npairs = DEFAULT_FIRST_BIN_NPAIRS,
    first_bin_pmin = DEFAULT_FIRST_BIN_PMIN,
    max_n_bins = DEFAULT_MAX_N_BINS,
    collection_time_diff_max = DEFAULT_COLLECTION_TIME_DIFF_MAX,
    bootstrap = 0,
    n_bins_to_fit = NULL,
    plot_ratio = F,
    plot_ratio_by_deme = F,
    verbose = T
) {
  ids = rownames(dist_mat)
  traits = traits[ids]
  if(!is.null(times)) times = times[ids]

  if(!is.null(demes)) {
    ids = ids[as.character(traits[ids]) %in% demes]
    traits = traits[ids]
    if(!is.null(times)) times = times[ids]
    dist_mat = dist_mat[ids,ids,drop = F]
  }

  ## Upper triangle only, zeros (identical sequences) kept as 1E-20
  if(!inherits(dist_mat, 'sparseMatrix')) {
    dist_mat = as.matrix(dist_mat)
    dist_mat[dist_mat == 0] = 1E-20
    dist_mat = Matrix::Matrix(dist_mat, sparse = T)
  }
  dist_mat = Matrix::triu(dist_mat, k = 1)
  dist_mat@x[dist_mat@x == 0] = 1E-20

  if(!is.null(hist_dist_min_max)) {
    dist_mat@x[dist_mat@x < min(hist_dist_min_max) | dist_mat@x > max(hist_dist_min_max)] = 0
    dist_mat = Matrix::drop0(dist_mat)
  }

  df = data.frame(id = ids,
                  traits = as.factor(as.character(traits[ids])))
  if(!is.null(times)) df$time = as.numeric(times[ids])
  if(is.null(times)) collection_time_diff_max = NULL
  if(!is.null(collection_time_diff_max)) {
    pairs = Matrix::summary(dist_mat)
    too_far = abs(df$time[pairs$i]-df$time[pairs$j]) > collection_time_diff_max
    dist_mat@x[too_far] = 0
    dist_mat = Matrix::drop0(dist_mat)
  }
  data_all = list(df,
                  dist_mat)

  demes_in_tree = sort(unique(as.character(traits)))
  hists = compute_histograms_helper(data_all, demes_in_tree, hist_break_method, first_bin_npairs, first_bin_pmin, max_n_bins,
                                    collection_time_diff_max, bootstrap, n_bins_to_fit, verbose)
  if(plot_ratio) plot_ratios_data(hists, by_deme = plot_ratio_by_deme)
  return(hists)
}
