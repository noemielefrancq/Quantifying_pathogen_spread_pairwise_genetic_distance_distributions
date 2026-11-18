## Internal plotting helpers

# Plot one observed-vs-predicted scatter panel (log-log) with optional CI and R^2
plot_obs_vs_pred_panel = function(obs, pred, pred_lo, pred_hi, plot_ci, main_lab) {
  obs[obs <= 0] = NA
  pred[pred <= 0] = NA
  lim = range(c(obs, pred), na.rm = T)
  plot(obs, pred, pch = 16, col = 'black',
       xlab = 'Observed proportion', ylab = 'Predicted proportion',
       xlim = lim, ylim = lim,
       log = 'xy', bty = 'n', main = main_lab, las = 2, xaxt = 'n')
  axis(1, las = 1)
  abline(0, 1, col = 'red', lty = 2)
  if(plot_ci) arrows(obs, pred_lo, obs, pred_hi, length = 0.05, angle = 90, code = 0, col = 'black')
  r2 = summary(lm(log10(pred) ~ log10(obs)))$r.squared
  text(10^par('usr')[1], 10^par('usr')[4], labels = paste0('R^2 = ', sprintf('%.2f', r2)), adj = c(0, 1))
}

# Plot one ratio panel (distance vs ratio, log y) with CI polygon
plot_ratio_panel = function(x, y_obs, y_pred, y_lo, y_hi, obs_lo, obs_hi, main_lab, ylim = NULL) {
  if(is.null(ylim)) {
    y = c(y_obs, y_pred, y_lo, y_hi, obs_lo, obs_hi)
    ylim = range(y[is.finite(y) & y > 0])
  }
  plot(x, y_obs, pch = 16, col = 'black', xlab = 'Distance',
       ylab = 'Ratio', main = main_lab,
       ylim = ylim, las = 2, log = 'y', xaxt = 'n', xlim = c(0, max(x, na.rm = T)))
  axis(1, las = 1)
  ok = is.finite(y_lo) & is.finite(y_hi)
  polygon(c(x[ok], rev(x[ok])), c(y_lo[ok], rev(y_hi[ok])), col = rgb(1, 0.6, 0, 0.25), border = NA)
  lines(x, y_pred, col = 'orange', lwd = 2)
  abline(h = 1, lty = 2, col = 'gray40')
  arrows(x, obs_lo, x, obs_hi, length = 0.05, angle = 90, code = 0, col = 'black')
  legend('topright',
         legend = c('Observed', 'Predicted', '95% CI'),
         col = c('black', 'orange', rgb(1, 0.6, 0, 0.25)),
         pch = c(16, NA, 15),
         lty = c(NA, 1, NA),
         lwd = c(NA, 2, NA),
         pt.cex = c(1, NA, 2),
         bty = 'n')
}

# Compute R^2 between observed and predicted proportion vectors
compute_r2 = function(obs, pred) {
  ok = is.finite(obs) & is.finite(pred)
  return(1-sum((obs[ok]-pred[ok])^2)/sum((obs[ok]-mean(obs[ok]))^2))
}

# Median and 95% interval of a 3D array over dimension `over`
summarise_dim = function(x, over) {
  x[!is.finite(x)] = NA
  mci = apply(x, setdiff(1:3, over), median_and_ci)
  return(list(median = array(mci[1,,], dim(mci)[-1]),
              lo = array(mci[2,,], dim(mci)[-1]),
              hi = array(mci[3,,], dim(mci)[-1])))
}

# Rows normalised to sum to 1 over the fitted bins (NA beyond)
normalise_fitted = function(m, fitted) {
  m[!fitted] = NA
  return(m/rowSums(m, na.rm = T))
}

#' Observed and predicted distributions for each pair of demes
#'
#' From a fit, computes for each fitted set (pair of demes, or directed set a->b when fitted by
#' deme) and bin, with each histogram divided by its number of pairs at any distance: the observed
#' within-deme distribution D(d|same) (D(d|aa) and D(d|bb) averaged with the set's reference weights,
#' (D(d|aa) + D(d|bb)) / 2 by pair) and between-deme distribution D(d|ab) (median and 95% over bootstrap
#' replicates), and the predicted ones (draws x sets x bins), P(d|ab) = D_ref(d) cr(d) and
#' D_ref(d) (1 - w cr(d)) / (1 - w), with w the weight of between pairs in the reference (1/2 by pair).
#' Distributions are normalised over the fitted bins of each set; ratios are D(d|same) / D(d|ab), predicted
#' (1 - w cr(d)) / ((1 - w) cr(d)).
#'
#' @param result List returned by `fit_hist_model()` (with `keep_stan_data = TRUE`).
#' @param draw Optional draw index; if NULL, all draws are used.
#' @return List of arrays (sets x bins, or sets x bins x replicates), and `fitted` (logical, sets x bins).
#' @keywords internal
pair_fit_summary = function(result, draw = NULL) {
  stan_data = result$stan_data
  if(is.null(stan_data) || is.null(result$chains$cr)) stop('Need result$stan_data and result$chains$cr (fit with keep_stan_data = TRUE).')
  N_sets = stan_data$N_sets
  set_pair = as.vector(stan_data$set_pair)
  n_bins = stan_data$max_N_bins
  R = stan_data$N_boot
  ref_weight = stan_data$ref_weight
  n_bins_set = as.vector(stan_data$N_bins_pair)[set_pair]
  fitted = matrix(F, N_sets, n_bins)
  for(i in 1:N_sets) fitted[i,1:n_bins_set[i]] = T

  ## Observed distributions, each histogram divided by its number of pairs
  same = array(0, c(N_sets, n_bins, R))
  between = array(0, c(N_sets, n_bins, R))
  obs_same = array(0, c(N_sets, n_bins, R))
  obs_diff = array(0, c(N_sets, n_bins, R))
  for(r in 1:R) {
    D_within_a = matrix(stan_data$hist_within_a[,,r], stan_data$N_pairs)[set_pair,,drop = F]/stan_data$n_within_a[set_pair,r]
    D_within_b = matrix(stan_data$hist_within_b[,,r], stan_data$N_pairs)[set_pair,,drop = F]/stan_data$n_within_b[set_pair,r]
    between[,,r] = matrix(stan_data$hist_diff[,,r], stan_data$N_pairs)[set_pair,,drop = F]/stan_data$n_diff[set_pair,r]
    same[,,r] = (ref_weight[,1]*D_within_a+ref_weight[,2]*D_within_b)/(ref_weight[,1]+ref_weight[,2])
    obs_same[,,r] = normalise_fitted(matrix(same[,,r], N_sets), fitted)
    obs_diff[,,r] = normalise_fitted(matrix(between[,,r], N_sets), fitted)
  }
  D_ref = apply((1-ref_weight[,3])*same+ref_weight[,3]*between, c(1, 2), mean)

  ## Predicted distributions and ratio, for each draw
  cr = array(result$chains$cr, c(dim(result$chains$cr)[1], N_sets, n_bins))
  if(!is.null(draw)) cr = cr[draw,,,drop = F]
  pred_same = array(0, dim(cr))
  pred_diff = array(0, dim(cr))
  pred_ratio = array(0, dim(cr))
  for(k in 1:dim(cr)[1]) {
    cr_draw = matrix(cr[k,,], N_sets)
    pred_same[k,,] = normalise_fitted(D_ref*(1-ref_weight[,3]*cr_draw), fitted)
    pred_diff[k,,] = normalise_fitted(D_ref*cr_draw, fitted)
    pred_ratio[k,,] = (1-ref_weight[,3]*cr_draw)/((1-ref_weight[,3])*cr_draw)
  }

  obs_ratio = summarise_dim(same/between, 3)
  pred_ratio = summarise_dim(pred_ratio, 1)
  pred_same = summarise_dim(pred_same, 1)
  pred_diff = summarise_dim(pred_diff, 1)
  return(list(distance = (stan_data$bin_break_min+stan_data$bin_break_max)/2,
              N_bins_pair = n_bins_set,
              fitted = fitted,
              pair_labels = result$set_labels,
              observed_same = apply(obs_same, c(1, 2), median),
              observed_diff = apply(obs_diff, c(1, 2), median),
              observed_same_boot = obs_same,
              observed_diff_boot = obs_diff,
              predicted_same = pred_same$median,
              predicted_same_lo = pred_same$lo,
              predicted_same_hi = pred_same$hi,
              predicted_diff = pred_diff$median,
              predicted_diff_lo = pred_diff$lo,
              predicted_diff_hi = pred_diff$hi,
              ratio_observed = obs_ratio$median,
              ratio_observed_lo = obs_ratio$lo,
              ratio_observed_hi = obs_ratio$hi,
              ratio_predicted = pred_ratio$median,
              ratio_predicted_lo = pred_ratio$lo,
              ratio_predicted_hi = pred_ratio$hi))
}

## Exported plotting functions

#' Compute variance explained (R-squared) for fitted histograms
#'
#' R-squared between observed and predicted same/diff proportions over all pairs and bins,
#' and for a null (uniform) prediction.
#'
#' @param result List returned by `fit_hist_model()`.
#' @return A list with `r2` and `null_r2`.
#' @export
compute_variance_explained = function(result) {
  fit_summary = pair_fit_summary(result)
  obs = c(fit_summary$observed_same[fit_summary$fitted], fit_summary$observed_diff[fit_summary$fitted])
  pred = c(fit_summary$predicted_same[fit_summary$fitted], fit_summary$predicted_diff[fit_summary$fitted])
  null = rep(1/matrix(fit_summary$N_bins_pair, nrow(fit_summary$fitted), ncol(fit_summary$fitted))[fit_summary$fitted], 2)
  return(list(r2 = compute_r2(obs, pred),
              null_r2 = compute_r2(obs, null)))
}

#' Plot observed vs predicted histogram proportions
#'
#' For each pair of demes, plots observed vs predicted proportions of the within-deme ("same") and
#' between-deme ("different") distributions over the fitted bins.
#'
#' @param result List returned by `fit_hist_model()`.
#' @param plot_ci Logical; include 95\% credible intervals.
#' @param plot Logical; if FALSE, only return the values.
#' @return Observed and predicted proportions and ratios (see `pair_fit_summary()`), invisibly when plotting.
#' @export
plot_proportion_fits = function(result, plot_ci = F, plot = T) {
  fit_summary = pair_fit_summary(result)
  if(!plot) return(fit_summary)
  par(mfrow = c(length(fit_summary$pair_labels), 2), las = 2, mgp = c(3, 0.7, 0))
  for(p in seq_along(fit_summary$pair_labels)) {
    n_bins = fit_summary$N_bins_pair[p]
    plot_obs_vs_pred_panel(obs = fit_summary$observed_same[p,1:n_bins], pred = fit_summary$predicted_same[p,1:n_bins],
                           pred_lo = fit_summary$predicted_same_lo[p,1:n_bins], pred_hi = fit_summary$predicted_same_hi[p,1:n_bins],
                           plot_ci = plot_ci, main_lab = paste('Pair', fit_summary$pair_labels[p], '- same'))
    plot_obs_vs_pred_panel(obs = fit_summary$observed_diff[p,1:n_bins], pred = fit_summary$predicted_diff[p,1:n_bins],
                           pred_lo = fit_summary$predicted_diff_lo[p,1:n_bins], pred_hi = fit_summary$predicted_diff_hi[p,1:n_bins],
                           plot_ci = plot_ci, main_lab = paste('Pair', fit_summary$pair_labels[p], '- different'))
  }
  return(invisible(fit_summary))
}

#' Plot ratio fits (observed vs predicted)
#'
#' Plots, for each pair of demes, the ratio R(d) of the within-deme ("same") to the between-deme
#' ("different") distribution on a log scale: observed (median and 95\% over bootstrap replicates)
#' and predicted (posterior median and 95\%, or a single draw).
#'
#' @param result List returned by `fit_hist_model()`.
#' @param draw Optional draw index (ignored if `use_median = TRUE`).
#' @param use_median Logical; plot the posterior median prediction.
#' @param same_y_lim Logical; use the same y-axis range for all pairs.
#' @param plot Logical; if FALSE, only return the values.
#' @return Observed and predicted ratios and proportions (see `pair_fit_summary()`), invisibly when plotting.
#' @export
plot_ratios_fit = function(result, draw = NULL, use_median = T, same_y_lim = F, plot = T) {
  if(!is.null(draw) && use_median) stop('Choose either use_median = TRUE or a draw index, not both.')
  if(use_median) draw = NULL
  fit_summary = pair_fit_summary(result, draw)
  if(!plot) return(fit_summary)
  n_pairs = length(fit_summary$pair_labels)
  n_col = ceiling(sqrt(n_pairs))
  par(mfrow = c(ceiling(n_pairs/n_col), n_col), mar = c(4, 5, 2, 1), bty = 'n')
  ylim = NULL
  if(same_y_lim) {
    y = c(fit_summary$ratio_observed[fit_summary$fitted], fit_summary$ratio_observed_lo[fit_summary$fitted], fit_summary$ratio_observed_hi[fit_summary$fitted],
          fit_summary$ratio_predicted[fit_summary$fitted], fit_summary$ratio_predicted_lo[fit_summary$fitted], fit_summary$ratio_predicted_hi[fit_summary$fitted])
    ylim = range(y[is.finite(y) & y > 0])
  }
  for(p in 1:n_pairs) {
    n_bins = fit_summary$N_bins_pair[p]
    plot_ratio_panel(x = fit_summary$distance[1:n_bins], y_obs = fit_summary$ratio_observed[p,1:n_bins], y_pred = fit_summary$ratio_predicted[p,1:n_bins],
                     y_lo = fit_summary$ratio_predicted_lo[p,1:n_bins], y_hi = fit_summary$ratio_predicted_hi[p,1:n_bins],
                     obs_lo = fit_summary$ratio_observed_lo[p,1:n_bins], obs_hi = fit_summary$ratio_observed_hi[p,1:n_bins],
                     main_lab = fit_summary$pair_labels[p], ylim = ylim)
  }
  return(invisible(fit_summary))
}

#' Plot ratios from precomputed histograms
#'
#' Plots the observed ratio R(d) of the within-deme distribution ("same", average of the within-a and
#' within-b distributions) to the between-deme distribution ("different"), each histogram divided by its
#' number of pairs at any distance, summed over pairs of demes, or for each pair (`by_deme = TRUE`).
#' Intervals are 95\% over bootstrap replicates, if available.
#'
#' @param pre_h Output from `compute_histograms()` or `compute_histograms_from_dist()`.
#' @param by_deme Logical; plot per-deme pairs instead of overall.
#' @return Invisibly returns plotted ratios.
#' @export
plot_ratios_data = function(pre_h, by_deme = F) {
  distance = pre_h$hists_all$hist_data$distance
  n_bins = length(distance)
  N_pairs = nrow(pre_h$pair_mat)
  if(is.null(pre_h$pair_hists_boot)) {
    pair_hists = pre_h$pair_hists
    n_boot = 1
  }else{
    pair_hists = pre_h$pair_hists_boot
    n_boot = dim(pair_hists$between)[3]
  }
  within_a = array(pair_hists$within_a, c(N_pairs, n_bins, n_boot))
  within_b = array(pair_hists$within_b, c(N_pairs, n_bins, n_boot))
  between = array(pair_hists$between, c(N_pairs, n_bins, n_boot))
  n_within_a = matrix(pair_hists$n_within_a, N_pairs, n_boot)
  n_within_b = matrix(pair_hists$n_within_b, N_pairs, n_boot)
  n_between = matrix(pair_hists$n_between, N_pairs, n_boot)
  ratios = array(NA, c(N_pairs+1, n_bins, n_boot))
  for(r in 1:n_boot) {
    same = (matrix(within_a[,,r], N_pairs)/n_within_a[,r]+matrix(within_b[,,r], N_pairs)/n_within_b[,r])/2
    diff = matrix(between[,,r], N_pairs)/n_between[,r]
    ratios[1:N_pairs,,r] = same/diff
    ratios[N_pairs+1,,r] = colSums(same)/colSums(diff)
  }
  ratios[!is.finite(ratios) | ratios <= 0] = NA
  mci = summarise_dim(ratios, 3)

  par(mfrow = c(1, 1), las = 2, bty = 'n')
  if(by_deme){
    rows = 1:N_pairs
    cols = grDevices::rainbow(N_pairs)
  }else{
    rows = N_pairs+1
    cols = 'black'
  }
  ylim = range(c(mci$median[rows,], mci$lo[rows,], mci$hi[rows,]), finite = T)
  plot(NA, xlim = range(distance), ylim = ylim, log = 'y', xaxt = 'n', xlab = 'Distance', ylab = 'Ratio (same/diff)',
       main = 'Ratio of histograms across demes')
  axis(1, las = 1)
  abline(h = 1, lty = 2, col = 'gray40')
  for(k in seq_along(rows)) {
    points(distance, mci$median[rows[k],], pch = 16, col = cols[k])
    if(n_boot > 1) arrows(distance, mci$lo[rows[k],], distance, mci$hi[rows[k],], length = 0.05, angle = 90, code = 0, col = cols[k])
  }
  if(by_deme) legend('topright', legend = pre_h$pair_labels, col = cols, pch = 16, bty = 'n')
  if(by_deme) return(invisible(list(ratio_by_pair = mci$median[1:N_pairs,,drop = F])))
  return(invisible(list(ratio_shared = mci$median[N_pairs+1,])))
}
