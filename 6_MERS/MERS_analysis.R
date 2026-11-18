## MERS analysis - transmission between hosts
library(ape)
library(stringr)
library(lubridate)
source('1_code/load_pairate.R', chdir = T)

#########################################################################################################################################
## Load data
#########################################################################################################################################
Months = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec')
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}

######################################
## MERS
######################################
metadata_MERS = read.csv('6_MERS/Data/MERS_CoV_274_CDS.dates.txt', sep = '\t', header = F)
tree_MERS = treeio::read.beast('6_MERS/Data/summary_relaxed_part2.tree')
tree_MERS = ladderize(tree_MERS@phylo)
alignment_MERS = seqinr::read.fasta('6_MERS/Data/MERS_CoV_274.fasta')
genome_length_MERS = NA
substitution_rate_MERS = NA

colnames(metadata_MERS) = c('seqname', 'date')
metadata_MERS$host = sapply(metadata_MERS$seqname, function(x)str_split(x, '\\|')[[1]][3])

metadata_MERS = metadata_MERS[which(metadata_MERS$date > max(min(metadata_MERS$date[which(metadata_MERS$host == 'human')]), min(metadata_MERS$date[which(metadata_MERS$host == 'camel')])) &
                                    metadata_MERS$date < min(max(metadata_MERS$date[which(metadata_MERS$host == 'human')]), max(metadata_MERS$date[which(metadata_MERS$host == 'camel')]))),]
tree_MERS = drop.tip(tree_MERS, tree_MERS$tip.label[which(is.na(match(tree_MERS$tip.label, metadata_MERS$seqname)))])
metadata_MERS = metadata_MERS[match(tree_MERS$tip.label, metadata_MERS$seqname),]

write.csv(metadata_MERS, '6_MERS/Data/Table_S7_List_isolates_MERS_CoV.csv')
#########################################################################################################################################

#########################################################################################################################################
## Compute ratios
#########################################################################################################################################
traits_MERS = metadata_MERS$host[match(tree_MERS$tip.label, metadata_MERS$seqname)]
times_MERS = metadata_MERS$date[match(tree_MERS$tip.label, metadata_MERS$seqname)]
demes_MERS = sort(unique(traits_MERS))
hists_MERS = compute_histograms(tree = tree_MERS,
                           traits = traits_MERS,
                           times = times_MERS,
                           demes = demes_MERS, 
                           
                           hist_dist_min_max = c(0,10),
                           hist_break_method = c('Manual'),
                           # first_bin_npairs = 10,
                           first_bin_npairs = c(seq(0, 1, 0.1), seq(1.5, 10, 0.5)),
                           first_bin_pmin = 5E-7,
                           max_n_bins = 1000,
                           
                           collection_time_diff_max = 0.1,

                           bootstrap = 400)

plot_ratios_data(hists_MERS)
#########################################################################################################################################

#########################################################################################################################################
## Fit ratios
#########################################################################################################################################
fit_MERS_ER = fit_hist_model(
  tree = tree_MERS,
  traits = traits_MERS,
  times = times_MERS,
  demes = demes_MERS,
  model = matrix(c(0,1,1,0), 2, byrow = T),
  parallel_chains = 3, chains = 3, 
  n_bins_to_fit = 16,
  iter_warmup = 2000, iter_sampling = 2000,
  
  precomputed_hists = hists_MERS, 
  
  bootstrap = 0
)
res_rr = plot_ratios_fit(fit_MERS_ER)
res_prop = plot_proportion_fits(fit_MERS_ER)
#########################################################################################################################################


#########################################################################################################################################
## Figure 
#########################################################################################################################################
# pdf(file = '6_MERS/Figures/FigureSXX_MERS_raw.pdf', 
    # width = 7.5, height = 5)
par(mfrow = c(2,3), oma = c(0,0,0,0), 
    mar = c(2.5,2.5,3,2.5))
root_height = max(metadata_MERS$date)-max(node.depth.edgelength(tree_MERS))
traits_MERS_color = factor(traits_MERS)
levels(traits_MERS_color) = c('orange2', 'royalblue')
traits_MERS_color = as.character(traits_MERS_color)
plot(tree_MERS, show.tip.label = F)
tiplabels(pch = 16, col = traits_MERS_color, cex = 1)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(2012, 2015, 1)-(root_height),
             lab_axis = seq(2012, 2015, 1), lwd = 1, tck=-0.025, mgp = c(0,0.5,0))
legend('bottomleft', pch = 16, legend = c('Camel', 'Human'), col = c('orange2', 'royalblue'), bty ='n')

## Estimated rate
plot(NULL, xlim = c(0,3), ylim = c(0.001, 10), log = 'y', bty = 'n', yaxt = 'n',
     xlab = '', ylab = 'Estimated rate', main = 'Estimated rate', xaxt = 'n')
axis(1, at = c(1,2), labels = c('Prior', 'Estimate'), padj = 1)
axis(2, las = 2, at = 10^seq(-3,1,1), labels = format(10^seq(-3,1,1), scientific = F))
prior = 2*exp(fit_MERS_ER$chains$rates_log_prior)
prior_mci = median_and_ci(prior)
points(1, prior_mci[1], pch = 16, col = 'grey', cex = 1.25)
arrows(x0 = 1, x1 = 1, y0 = prior_mci[2], y1 = prior_mci[3], 
       length = 0.05, angle = 90, code = 0, col = 'grey', lwd = 2)
rate = 2*fit_MERS_ER$chains$rates
rate_mci = median_and_ci(rate)
points(2, rate_mci[1], pch = 16, col = 'darkred', cex = 1.25)
arrows(x0 = 2, x1 = 2, y0 = rate_mci[2], y1 = rate_mci[3], 
       length = 0.05, angle = 90, code = 0, col = 'darkred', lwd = 2)

## Compare to BDMM-prime
log_BDMM = beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.10.log'))
log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.9.log')))
log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.8.log')))
log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.7.log')))
log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.6.log')))
d_c = log_BDMM[,'becomeUninfectiousRateSPEpi.camel']
d_h = log_BDMM[,'becomeUninfectiousRateSPEpi.human']
g_c = d_c*(log_BDMM[,'ReSPEpi.camel']-1)
g_h = d_h*(log_BDMM[,'ReSPEpi.human']-1)
growth = (g_c+g_h+sqrt((g_c-g_h)^2+4*d_c*log_BDMM[,'ReAmongDemesSPEpi.camel_to_human']*d_h*log_BDMM[,'ReAmongDemesSPEpi.human_to_camel']))/2
rate_BDMM = growth-g_h
ymax = max(hist(rate_BDMM, breaks = 10^seq(-3,1.5,0.1), plot = F)$density, hist(rate, breaks = 10^seq(-2,1.5,0.05), plot = F)$density)
hist(rate_BDMM, border = NA, yaxt = 'n', ylim = c(0, 1.3*ymax), xlim = c(0.01, 10^1.5),
     freq = F, log = 'x', breaks = 10^seq(-3,1.5,0.1), col = adjustcolor('#1B9E77', alpha.f = 0.5),
     main = 'Camel -> Human', ylab = 'Density', xlab = 'Rate (1/year)')
hist(rate, add = T, border = NA,
     freq = F, breaks = 10^seq(-2,1.5,0.05), col = adjustcolor('darkred', alpha.f = 0.5))
axis(2, las = 2)
## Median and 95% CrI of each model
t = median_and_ci(rate_BDMM)
points(t[1], 1.25*ymax, pch = 16, col = '#1B9E77')
arrows(x0 = t[2], x1 = t[3], y0 = 1.25*ymax, y1 = 1.25*ymax, length = 0.05, angle = 90, code = 0, col = '#1B9E77', lwd = 2)
t = median_and_ci(rate)
points(t[1], 1.15*ymax, pch = 16, col = 'darkred')
arrows(x0 = t[2], x1 = t[3], y0 = 1.15*ymax, y1 = 1.15*ymax, length = 0.05, angle = 90, code = 0, col = 'darkred', lwd = 2)
# legend('topright', legend = c('BDMM-prime', 'This model'), pch = 15, col = c('royalblue', 'darkred'))

## Ratio overall: observed ((D(d|camel) + D(d|human)) / 2) / D(d|camel-human) (median and 95% over bootstrap replicates), and predicted by the fit (median and 95% CrI)
b = hists_MERS$pair_hists_boot
ratio_obs = ((t(t(b$within_a[1,,])/colSums(b$within_a[1,,]))+t(t(b$within_b[1,,])/colSums(b$within_b[1,,])))/2)/t(t(b$between[1,,])/colSums(b$between[1,,]))
ratio_obs[!is.finite(ratio_obs) | ratio_obs <= 0] = NA
ratio_obs = apply(ratio_obs, 1, median_and_ci)
distance = hists_MERS$hists_all$hist_data$distance
n_fit = res_rr$N_bins_pair[1]
plot(distance, ratio_obs[1,], log = 'y', pch = 16, bty = 'n', xaxt = 'n', yaxt = 'n',
     ylab = 'Ratio', xlab = 'Distance', ylim = c(0.05, 500), xlim = c(0,4), main = 'Ratio overall')
arrows(x0 = distance, x1 = distance, y0 = ratio_obs[2,], y1 = ratio_obs[3,], length = 0.05, angle = 90, code = 0)
polygon(x = c(res_rr$distance[1:n_fit], rev(res_rr$distance[1:n_fit])), y = c(res_rr$ratio_predicted_lo[1,1:n_fit], rev(res_rr$ratio_predicted_hi[1,1:n_fit])),
        col = adjustcolor('darkred', alpha.f = 0.5), border = F)
lines(res_rr$distance[1:n_fit], res_rr$ratio_predicted[1,1:n_fit], col = 'darkred')
axis(1, las = 1, at = seq(0,4,1))
axis(2, las = 2, at = 10^seq(-1,2), labels = format(10^seq(-1,2), scientific = F))
abline(h=1, lty = 2)
legend('topright', pch = c(16, NA), lty = c(NA, 1), col = c('black', 'darkred'), legend = c('Data', 'Model fit'), bty = 'n')

## Pred vs Observed proportions, same host and different hosts
plot(c(res_rr$observed_same[1,1:n_fit], res_rr$observed_diff[1,1:n_fit]), c(res_rr$predicted_same[1,1:n_fit], res_rr$predicted_diff[1,1:n_fit]), pch = 16, log = 'xy', xaxt = 'n', yaxt = 'n', xlim = c(0.001, 1), ylim = c(0.001, 1),
     bty ='n', xlab = 'Observed', ylab = 'Predicted', main = 'Model fit', col = rep(c('darkred', 'darkred'), each = n_fit))
axis(2, las = 2, at = c(0.001, 0.01, 0.1, 1))
axis(1, las = 1, at = c(0.001, 0.01, 0.1, 1))
abline(a=0, b=1, lty = 2)
# dev.off()
#########################################################################################################################################


#########################################################################################################################################
# Compare to BDMM-prime result
#########################################################################################################################################
# log_BDMM = beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.10.log'))
# log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.9.log')))
# log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.8.log')))
# log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.7.log')))
# log_BDMM = rbind(log_BDMM, beastio::readLog(filenames = c('6_MERS/Data/BEAST_log/mers_relaxed_part2.6.log')))
# colnames(log_BDMM)
# becomeUninfectiousRateSPEpicamel = log_BDMM[,27]
# ReAmongDemesSPEpicamel_to_human  = log_BDMM[,40]
# 
# becomeUninfectiousRateSPEpihuman = log_BDMM[,28]
# ReAmongDemesSPEpihuman_to_camel  = log_BDMM[,41]
# ##    ReAmongDeme * N generation per year = 0.5558129 (95%CI: 0.2543609, 1.0476226)
# 
# par(mfrow = c(1,2), mar = c(4,4,3,0))
# hist(ReAmongDemesSPEpicamel_to_human*becomeUninfectiousRateSPEpicamel, border = NA, yaxt = 'n',
#      freq = F, log = 'x', breaks = 10^seq(-2,1.5,0.1), col = adjustcolor('royalblue', alpha.f = 0.5),
#      main = 'Camel -> Human', ylab = 'Density', xlab = 'Rate (1/year)')
# hist(fit_MERS_ER$chains$rates, add = T, border = NA,
#      freq = F, breaks = 10^seq(-2,1.5,0.05), col = adjustcolor('darkred', alpha.f = 0.5))
# axis(2, las = 2)
# legend('topleft', legend = c('BDMM-prime', 'This model'), pch = 15, col = c('royalblue', 'darkred'))
# 
# hist(ReAmongDemesSPEpihuman_to_camel*becomeUninfectiousRateSPEpihuman, border = NA, yaxt = 'n',
#      freq = F, log = 'x', breaks = 10^seq(-5.5,1.5,0.1), col = adjustcolor('royalblue', alpha.f = 0.5),
#      main = 'Human -> Camel', ylab = 'Density', xlab = 'Rate (1/year)')
# axis(2, las = 2)
# legend('topright', legend = c('BDMM-prime', 'This model'), pch = 15, col = c('royalblue', 'darkred'))
# 
# hist((ReAmongDemesSPEpihuman_to_camel*becomeUninfectiousRateSPEpihuman+ReAmongDemesSPEpicamel_to_human*becomeUninfectiousRateSPEpicamel)/2, border = NA, yaxt = 'n',
#      freq = F, log = 'x', breaks = 10^seq(-2,1.5,0.1), col = adjustcolor('royalblue', alpha.f = 0.5),
#      main = 'Human -> Camel', ylab = 'Density', xlab = 'Rate (1/year)')
# axis(2, las = 2)
# 
# 
# par(mfrow = c(1,2), mar = c(4,4,3,0))
# hist(ReAmongDemesSPEpicamel_to_human*becomeUninfectiousRateSPEpicamel, border = NA, yaxt = 'n',
#      freq = F, log = 'x', breaks = 10^seq(-7,1.5,0.05), col = adjustcolor('royalblue', alpha.f = 0.5),
#      main = 'Camel -> Human', ylab = 'Density', xlab = 'Rate (1/year)')
# hist(fit_MERS_ARD$chains$trans_rate[,2], add = T, border = NA,
#      freq = F, breaks = 10^seq(-7,1.5,0.05), col = adjustcolor('darkred', alpha.f = 0.5))
# axis(2, las = 2)
# legend('topleft', legend = c('BDMM-prime', 'This model'), pch = 15, col = c('royalblue', 'darkred'))
# 
# hist(ReAmongDemesSPEpihuman_to_camel*becomeUninfectiousRateSPEpihuman, border = NA, yaxt = 'n',
#      freq = F, log = 'x', breaks = 10^seq(-7,1.5,0.05), col = adjustcolor('royalblue', alpha.f = 0.5),
#      main = 'Human -> Camel', ylab = 'Density', xlab = 'Rate (1/year)')
# hist(fit_MERS_ARD$chains$trans_rate[,1], add = T, border = NA,
#      freq = F, breaks = 10^seq(-7,1.5,0.05), col = adjustcolor('darkred', alpha.f = 0.5))
# axis(2, las = 2)
# legend('topright', legend = c('BDMM-prime', 'This model'), pch = 15, col = c('royalblue', 'darkred'))
#########################################################################################################################################



