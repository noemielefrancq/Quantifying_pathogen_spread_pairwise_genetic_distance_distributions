#########################################################################################################################################
## Spread across countries worldwide - SC2 NCBI data - Figures
#########################################################################################################################################

#########################################################################################################################################
## Libraries
#########################################################################################################################################
library(ape)
library(stringr)
library(lubridate)
library(Matrix)
library(MetBrewer)

source('1_code/load_pairate.R', chdir = T)

median.and.ci = median_and_ci = function(v){ return( c(median(v, na.rm = T), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
#########################################################################################################################################

#########################################################################################################################################
## Parameters SC2
#########################################################################################################################################
genome_length_SC2 = 29903
substitution_rate_SC2 = 5*10^-4 ## Nextstrain
#########################################################################################################################################

#########################################################################################################################################
## Load full data, distances from tree (iqtree)
#########################################################################################################################################
## Read file precomputed on cluster
dist_mat = readRDS('4_SC2_geography/Data_countries/Seqs/mat_sparse_7_countries.rds')
seq_names = colnames(dist_mat)
n = length(seq_names)
colnames(dist_mat) = rownames(dist_mat) = 1:ncol(dist_mat)
dist_mat@x = dist_mat@x*genome_length_SC2 ## Rescale branch lengths, as iqtree default is N subs per site
dist_mat@x[dist_mat@x > 40] <- 0
dist_mat <- drop0(dist_mat)

# Metadata from NCBI
metadata_raw = read.csv('4_SC2_geography/Data_countries/Seqs/Table_S1_List_isolates_SC2_geo.csv')
metadata_raw = metadata_raw[match(sapply(seq_names,function(x)stringr::str_split(x, '\\.1')[[1]][1]), 
                                  metadata_raw$Accession),]
metadata_raw$id = 1:n

#########################################################################################################################################

#########################################################################################################################################
## Subset data for specific countries
#########################################################################################################################################
countries = c("Denmark", "Japan",  "Slovakia",
              "Switzerland", "Thailand",  "United Kingdom", "USA")
countries_iso2 = c('DK', 'JP', 'SK', 'CH', 'TH', 'UK', 'USA')

countries_colors = met.brewer(name = 'Java', n = length(countries_iso2))

keep = which(!is.na(match(metadata_raw$Country, countries)))
metadata_raw = metadata_raw[keep,]
dist_mat = dist_mat[keep,]
dist_mat = dist_mat[,keep]

traits = setNames(metadata_raw$Country, metadata_raw$id)
times  = setNames(metadata_raw$Collection_Date, metadata_raw$id)
bootstrap = 200
#########################################################################################################################################

######################################
## Read previous results
######################################
fit = readRDS('4_SC2_geography/Data_countries/Seqs/Fit_7_countries_0_15_2000_2000_relrate.rds')
hists = readRDS('4_SC2_geography/Data_countries/Seqs/Hists_7_countries_relative.rds')
countries_pairs_colors = met.brewer(name = 'Signac', n = length(hists$pair_labels))
######################################


#########################################################################################################################################
## Figure 3: SARS-CoV-2 spread between countries
#########################################################################################################################################
####################### A-B Data counts ####################### 
# pdf(file = '4_SC2_geography/Figures/Figure3/Figure3AB_Data_counts.pdf', width = 8, height = 2.5)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 2))
par(mar = c(2,2,2,2))

## Number of sequences through time, and countries
time_bins = seq(2024.49,2025.51,1/52)
N_seqs = matrix(NA, ncol = length(countries), nrow = length(time_bins)-1)
for(i in 1:length(countries)){
  N_seqs[,i] = table(cut(metadata_raw$Collection_Date[which(metadata_raw$Country == countries[i])], breaks = time_bins, include.lowest = T, ))
}

plot(colSums(N_seqs), type = 'h', ylim = c(0,5000), xlim = c(0.5, 7.5), lwd = 15, col = countries_colors,
     main = 'Number of sequences (2024-2025)', xlab = 'Country', ylab = 'Number', bty = 'n', yaxt = 'n', xaxt = 'n', lend = "butt")
axis(2, las = 2)
axis(1, at = 1:length(countries_iso2), labels = countries_iso2)

plot(NULL, xlim = c(2024.5,2025.5), ylim = c(1, max(N_seqs)), 
     yaxt = 'n', xaxt = 'n', bty = 'n',
     xlab = 'Time', ylab = 'Count', main = 'N sequences through time')
for(i in 1:length(countries)){
  lines((time_bins[-1] + time_bins[-length(time_bins)]) / 2, N_seqs[,i], col = countries_colors[i], lwd = 2) 
}
axis(2, las = 2)
axis(1, at = decimal_date(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))),
     paste0(month(as_date(c('2024/07/01', '2024/09/01', '2024/11/01', '2025/01/01', '2025/03/01', '2025/05/01', '2025/07/01'))), ' - ', c(24,24,24,25,25,25,25)))
# legend('topright', pch = 16, legend = countries_iso2, col = countries_colors)
# dev.off()
##################################################################### 

####################### C-D Ratios ####################### 
# pdf(file = '4_SC2_geography/Figures/Figure3/Figure3CD_Ratios_data_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))

res_within = res_between = matrix(0, nrow = dim(hists$pair_hists_boot$within_a)[2], ncol = 200)
for(i in 1:dim(hists$pair_hists_boot$within_a)[1]){
  for(b in 1:200){
    res_within[,b] = res_within[,b] + (hists$pair_hists_boot$within_a[i,,b]/sum(hists$pair_hists_boot$within_a[i,,b])) + (hists$pair_hists_boot$within_b[i,,b]/sum(hists$pair_hists_boot$within_b[i,,b]))
    res_between[,b] = res_between[,b] + (hists$pair_hists_boot$between[i,,b]/sum(hists$pair_hists_boot$between[i,,b]))
  }
}
res_within = apply(res_within, MARGIN = 2, function(x)x/sum(x))
res_between = apply(res_between, MARGIN = 2, function(x)x/sum(x))
tmp = sapply(1:ncol(res_within), function(x)res_within[,x]/res_between[,x])
tmp[which(is.infinite(tmp))] = NA
tmp = apply(tmp, MAR = 1, median_and_ci)
tmp[2,which(tmp[2,] == 0)] = 0.001
tmp[3,which(is.infinite(tmp[3,]))] = 1000
plot(((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2), tmp[1,], log = 'y', pch = 16, bty = 'n', 
     xaxt = 'n', yaxt = 'n',
     xlim = c(0,15), ylim = c(0.5, 500), ylab = 'Ratio', xlab = 'Genetic distance', main = 'Across all countries')
axis(1, at = c(0,5,10,15))
axis(2, at = c(0.5, 1, 10, 100, 500), labels = c('<0.5', 1, 10, 100, '>500'), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100, 1000))
arrows(x0 = ((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2), y0 = tmp[2,],
       x1 = ((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2), y1 = tmp[3,],
       length = 0.05, angle = 90, code = 0)
abline(h = 1, lty = 2)

plot(NULL, log = 'y', pch = 16, bty = 'n',  xaxt = 'n',  yaxt = 'n',
     xlim = c(0,15), ylim = c(0.5, 500), ylab = 'Ratio', xlab = 'Genetic distance', main = 'Per country pair')
axis(1, at = c(0,5,10,15))
axis(2, at = c(0.5, 1, 10, 100, 500), labels = c('<0.5', 1, 10, 100, '>500'), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100, 1000))
abline(h = 1, lty = 2)
for(p in 1:length(hists$pair_labels)){
  res_within = res_between = matrix(0, nrow = dim(hists$pair_hists_boot$within_a)[2], ncol = 200)
  for(b in 1:200){
    res_within[,b] = res_within[,b] + (hists$pair_hists_boot$within_a[p,,b]/sum(hists$pair_hists_boot$within_a[p,,b])) + (hists$pair_hists_boot$within_b[p,,b]/sum(hists$pair_hists_boot$within_b[p,,b]))
    res_between[,b] = res_between[,b] + (hists$pair_hists_boot$between[p,,b]/sum(hists$pair_hists_boot$between[p,,b]))
  }
  res_within = apply(res_within, MARGIN = 2, function(x)x/sum(x))
  res_between = apply(res_between, MARGIN = 2, function(x)x/sum(x))
  tmp = sapply(1:ncol(res_within), function(x)res_within[,x]/res_between[,x])
  tmp[which(is.infinite(tmp))] = NA
  tmp = apply(tmp, MAR = 1, median_and_ci)
  tmp[2,which(tmp[2,] == 0)] = 0.001
  tmp[3,which(is.infinite(tmp[3,]))] = 1000
  print(p)
  print(tmp[,1])
  points(((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2), tmp[1,], 
         pch = 16, col = countries_pairs_colors[p])
  lines(((hists$common_breaks[-1]+hists$common_breaks[-length(hists$common_breaks)])/2), tmp[1,], 
        col = countries_pairs_colors[p], lwd = 0.5)
}

plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(0,1), ylim = c(1,1))
legend('center', col = countries_pairs_colors, pch = 15, legend = hists$pair_labels, bty = 'n', cex = 0.5)
# dev.off()
#####################################################################

####################### E Model fits - proportions and RR #############################
# pdf(file = '4_SC2_geography/Figures/Figure3/Figure3E_Fits_proportions_ratio_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))

mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
p_names = which(lower.tri(mobility_flow_mat_sym)==T, arr.ind = T)
p_names = paste0(p_names[,2], p_names[,1])

fit_summary = plot_proportion_fits(fit, plot = F)
pair_colors = c()
for(p in 1:fit$stan_data$N_sets){
  pr = fit$stan_data$set_pair[p]
  p_id = paste0(sort(c(fit$stan_data$pair_a[pr], fit$stan_data$pair_b[pr])), collapse = '')
  p_id = match(p_id, p_names)
  pair_colors[p] = countries_pairs_colors[p_id]
}

plot(NULL, xlim = c(10^-(4),0.5), ylim = c(10^-(4),0.5), bty = 'n', xaxt = 'n', yaxt = 'n', log = 'xy', main = 'P same', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(1*10^(seq(-4,-1, 1)), 0.5), labels = c(1*10^(seq(-4,-1, 1)), 0.5))
axis(2, las = 2, at = c(1*10^(seq(-4,-1, 1)), 0.5), labels = c(1*10^(seq(-4,-1, 1)), 0.5))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  bins = which(fit_summary$fitted[p,])
  points(fit_summary$observed_same[p,bins], fit_summary$predicted_same[p,bins], pch = 16, col = pair_colors[p])
}
plot(NULL, xlim = c(10^-(4),0.5), ylim = c(10^-(4),0.5), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = 'P diff', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(1*10^(seq(-4,-1, 1)), 0.5), labels = c(1*10^(seq(-4,-1, 1)), 0.5))
axis(2, las = 2, at = c(1*10^(seq(-4,-1, 1)), 0.5), labels = c(1*10^(seq(-4,-1, 1)), 0.5))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  bins = which(fit_summary$fitted[p,])
  points(fit_summary$observed_diff[p,bins], fit_summary$predicted_diff[p,bins], pch = 16, col = pair_colors[p])
}
plot(NULL, xlim = c(0.5,500), ylim = c(0.5,500), bty = 'n', xaxt = 'n',  yaxt = 'n', log = 'xy', main = 'Ratio', xlab = 'Observed', ylab = 'Predicted (median)')
axis(1, at = c(0.1, 1, 10, 100, 500), labels = c('<0.1', 1, 10, 100, '500'), las = 1, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
axis(2, at = c(0.1, 1, 10, 100, 500), labels = c('<0.1', 1, 10, 100, '500'), las = 2, tick = c(0.1, 0.5, 1, 5, 10, 50, 100))
abline(a=0, b=1, lty= 2)
for(p in 1:fit$stan_data$N_sets){
  pr = fit$stan_data$set_pair[p]
  w = fit$stan_data$ref_weight[p,]
  bins = which(fit_summary$fitted[p,])
  ref = 0
  for(i in 1:fit$stan_data$N_boot){
    same = (w[1]*fit$stan_data$hist_within_a[pr,bins,i]/fit$stan_data$n_within_a[pr,i] + w[2]*fit$stan_data$hist_within_b[pr,bins,i]/fit$stan_data$n_within_b[pr,i])/(w[1]+w[2])
    ref = ref + ((1-w[3])*same + w[3]*fit$stan_data$hist_diff[pr,bins,i]/fit$stan_data$n_diff[pr,i])/fit$stan_data$N_boot
  }
  pred_same = t(apply(1 - w[3]*fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
  pred_diff = t(apply(fit$chains$cr[,p,bins], MAR = 1, function(x)(x*ref)/sum(x*ref)))
  data = fit_summary$observed_same[p,bins]/fit_summary$observed_diff[p,bins]
  points(data, apply(pred_same/pred_diff, MARGIN = 2, FUN = median), pch = 16, col = pair_colors[p])
}
# dev.off()
###########################################

####################### F Matrices ####################### 
# pdf(file = '4_SC2_geography/Figures/Figure3/Figure3F_Matrices_estimates_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
rates = fit$chains$rates
ts =365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median.and.ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median.and.ci(rates[,i,])
}
estimated_mat = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
estimated_mat[lower.tri(estimated_mat)] = res_ci[1,]/ts
estimated_mat = estimated_mat+t(estimated_mat)
diag(estimated_mat) = NA
colnames(estimated_mat) = hists$demes_in_tree
rownames(estimated_mat) = hists$demes_in_tree
estimated_mat = log10(estimated_mat)

zlim = range( log10(mobility_flow_mat_sym), na.rm = T)
cols = hcl.colors(50, "YlOrRd", rev = TRUE)
mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
image(1:nrow(mobility_flow_mat_sym), 1:ncol(mobility_flow_mat_sym), 
      log10(mobility_flow_mat_sym),
      axes = FALSE, xlab = "", ylab = "",
      col = cols, zlim = zlim,
      main = names(estimated_mat))
axis(1, at = 1:nrow(mobility_flow_mat_sym), labels = rownames(mobility_flow_mat_sym), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(mobility_flow_mat_sym), labels = colnames(mobility_flow_mat_sym), las = 2, cex.axis = 0.7)

zlim = range(estimated_mat, na.rm = T)
cols = hcl.colors(50, "YlOrRd", rev = TRUE)
colnames(estimated_mat) = rownames(mobility_flow_mat_sym)
rownames(estimated_mat) = rownames(mobility_flow_mat_sym)
mat = estimated_mat
image(1:nrow(mat), 1:ncol(mat), mat,
      axes = FALSE, xlab = "", ylab = "",
      col = cols, zlim = zlim,
      main = names(estimated_mat))
axis(1, at = 1:nrow(mat), labels = rownames(mat), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(mat), labels = colnames(mat), las = 2, cex.axis = 0.7)

plot.new()
par(usr = c(0, 1, zlim[1], zlim[2]))
rasterImage(as.raster(matrix(rev(cols), ncol = 1)), 0.2, zlim[1], 0.4, zlim[2])
axis(4, at = pretty(zlim), labels = 10^pretty(zlim), las = 1, cex.axis = 0.8)
# dev.off()
################################################################# 

####################### H Plot against mobility data ####################### 
# pdf(file = '4_SC2_geography/Figures/Figure3/Figure3H_Rates_vs_data_13102026.pdf', width = 8, height = 2.5)
par(mfrow = c(1,3), mar = c(2,2,1,2), oma = c(1,1,1,1))
plot.new()
## Raw rates
idx_steps = c(1:nrow(fit$chains$rates))
rates = fit$chains$rates[idx_steps,]
ts = 365/(substitution_rate_SC2*genome_length_SC2)
# ts = 1
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median_and_ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median_and_ci(rates[,i,])
}

## Plot against mobility data
mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
mobility_flow_mat_sym[upper.tri(mobility_flow_mat_sym)] = NA
plot(mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))],
     res_ci[1,]/ts, pch = 16, 
     ylim = c(1e-4, 1e-1),  xlim = c(1e-4, 0.5), bty = 'n', yaxt = 'n', log = 'xy', col = countries_pairs_colors, xaxt = 'n',
     xlab = 'Normalized number of trips between countries per year', ylab = 'Rates from genetic distances (per day)')
axis(1, at = 5*10^seq(-4,-1,1), labels = format(5*10^seq(-4,-1,1), scientific = T))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
arrows(x0 = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))],
       y0 = res_ci[2,]/ts, 
       x1 = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))],
       y1 = res_ci[3,]/ts, 
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)
## Compute model fit
slopes = NULL
intercept = NULL
o = order(mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))])
x = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))][o]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)[o]
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-10, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the countries
mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
mat = mobility_flow_mat_sym
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value

plot.new()
a = summary(m)
legend('topleft', title = paste0('p.value = ', round(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt), digits = 5)), legend = '', bty = 'n')
# dev.off()
############################################################################ 
####################### end ####################### 

#########################################################################################################################################
## Figure SXX: Plot rates and Rates vs covariates (N seqs and geo distance)
#########################################################################################################################################
# pdf(file = '4_SC2_geography/Figures/FigureSXX_Rates_and_covariates-geo-Nseqs_13102026.pdf', width = 6, height = 6)
layout(matrix(c(1,1,2,3), nrow = 2, byrow = T), widths = c(1, 1))
par(mar = c(2,2,4,2), oma = c(1,1,1,1))

## Raw rates
idx_steps = c(1:nrow(fit$chains$rates))
rates = fit$chains$rates[idx_steps,]
ts = 365/(substitution_rate_SC2*genome_length_SC2)
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median_and_ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median_and_ci(rates[,i,])
}
plot(1:length(fit$rate_labels), res_ci[1,]/ts, pch = 16, 
     ylim = c(1e-4, 1e-1), bty = 'n', yaxt = 'n', xaxt = 'n',
     log = 'y', col = countries_pairs_colors, 
     xlab = 'Country pair', ylab = 'Rates from genetic distances (per day)')
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
axis(1, 1:length(fit$rate_labels), fit$rate_labels)
arrows(x0 = 1:length(fit$rate_labels), y0 = res_ci[2,]/ts, 
       x1 = 1:length(fit$rate_labels), y1 = res_ci[3,]/ts, code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)


N_seqs_countries = matrix(0, ncol = length(hists$demes_in_tree), nrow = length(hists$demes_in_tree))
colnames(N_seqs_countries) = hists$demes_in_tree
rownames(N_seqs_countries) = hists$demes_in_tree
for(i in 1:length(hists$demes_in_tree)){
  for(j in 1:length(hists$demes_in_tree)){
    N = length(which(traits == hists$demes_in_tree[i])) + length(which(traits == hists$demes_in_tree[j]))
    N_seqs_countries[i,j] = N
  }
}
diag(N_seqs_countries) = NA
N_seqs_countries_chosen = N_seqs_countries
N_seqs_countries_chosen[upper.tri(N_seqs_countries_chosen)] = NA
plot(N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))],  # [lower.tri(N_seqs_countries_chosen)], 
     res_ci[1,]/ts, pch = 16, xlim = c(1e2, 1e4), xaxt = 'n',
     ylim = c(1e-4, 1e-1), bty = 'n', yaxt = 'n', log = 'xy',
     col = countries_pairs_colors,
     xlab = 'Nb. sequences per pair', ylab = 'Rates from genetic distances (per day)')
axis(1, at = 10^seq(2,4,1), labels = format(10^seq(2,4,1), scientific = T))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
arrows(x0 = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))], #[lower.tri(N_seqs_countries_chosen)], 
       y0 = res_ci[2,]/ts,
       x1 = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))], #[lower.tri(N_seqs_countries_chosen)], 
       y1 = res_ci[3,]/ts,
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)

## Compute model fit
slopes = NULL
intercept = NULL
o = order(N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))])
x = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))][o]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)[o]
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 10, 0.1))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the countries
N_seqs_countries_chosen = N_seqs_countries
mat = N_seqs_countries_chosen
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value


## Plot against geographical distance (between capital cities)
compute_country_distances = function(countries) {                                                      
  if (!requireNamespace("maps", quietly = TRUE)) stop('Package "maps" is required. Install with 
  install.packages("maps").')                                                                            
  
  db = maps::map("world", plot = FALSE, fill = TRUE)                                                   
  
  haversine = function(lat1, lon1, lat2, lon2) {
    R = 6371
    dlat = (lat2 - lat1) * pi / 180
    dlon = (lon2 - lon1) * pi / 180
    a = sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
    R * 2 * atan2(sqrt(a), sqrt(1 - a))
  }
  
  # Compute centroid for a country by averaging polygon coordinates
  get_centroid = function(country) {
    idx = grep(paste0("^", country, "(:|$)"), db$names)
    if (length(idx) == 0) stop('Country not found: ', country,
                               '. Check maps::map("world")$names for valid names.')
    x = db$x[unlist(lapply(idx, function(i) {
      start = ifelse(i == 1, 1, sum(is.na(db$x[1:(i-1)])) + i)
      # simpler: just use range matching
      NULL
    }))]
    # Extract polygon coords for matching regions
    sub = maps::map("world", regions = db$names[idx], plot = FALSE, fill = TRUE)
    lon = sub$x[!is.na(sub$x)]
    lat = sub$y[!is.na(sub$y)]
    c(mean(lat), mean(lon))
  }
  
  n = length(countries)
  coords = t(sapply(countries, get_centroid))
  
  d = matrix(0, nrow = n, ncol = n, dimnames = list(countries, countries))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d[i, j] = haversine(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2])
      d[j, i] = d[i, j]
    }
  }
  d
}
geo_distance_countries = compute_country_distances(c("Denmark", "Japan", "Slovakia", "Switzerland", "Thailand", "UK", "USA"))
diag(geo_distance_countries) = NA

geo_distance_countries_chosen = geo_distance_countries
geo_distance_countries_chosen[upper.tri(geo_distance_countries_chosen)] = NA
plot(geo_distance_countries_chosen[which(!is.na(geo_distance_countries_chosen))], #[lower.tri(geo_distance_countries_chosen)], 
     res_ci[1,]/ts, pch = 16, xaxt = 'n', 
     xlim = c(500, 15000),
     ylim = c(1e-4, 1e-1), bty = 'n', yaxt = 'n', log = 'xy',
     col = countries_pairs_colors,
     xlab = 'Geographical distance (km)', ylab = 'Rates from genetic distances (per day)')
axis(1, at = c(500, 1000, 5000, 10000))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
arrows(x0 = geo_distance_countries_chosen[which(!is.na(geo_distance_countries_chosen))], #[lower.tri(geo_distance_countries_chosen)], 
       y0 = res_ci[2,]/ts,
       x1 = geo_distance_countries_chosen[which(!is.na(geo_distance_countries_chosen))], #[lower.tri(geo_distance_countries_chosen)], 
       y1 = res_ci[3,]/ts,
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)
## Compute model fit
slopes = NULL
intercept = NULL
o = order(geo_distance_countries_chosen[lower.tri(geo_distance_countries_chosen)])
x = geo_distance_countries_chosen[lower.tri(geo_distance_countries_chosen)][o]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)[o]
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-10, 10, 0.1))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the countries
geo_distance_countries_chosen = geo_distance_countries
mat = geo_distance_countries_chosen
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt<=median(slopes[draws])))/length(slopes_permt)) ## p value

# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: Ratio fits
#########################################################################################################################################
val = plot_ratios_fit(fit, plot = F)

# pdf(file = '4_SC2_geography/Figures/FigureSXX_Fits_ratios_per_pair_13102026.pdf', width = 5, height = 11)
obs_m = val$ratio_observed
obs_cimin = val$ratio_observed_lo
obs_cimax = val$ratio_observed_hi
pred_m = val$ratio_predicted
pred_cimin = val$ratio_predicted_lo
pred_cimax = val$ratio_predicted_hi

par(mfrow = c(7,3), cex.main = 0.7,
    mar = c(3, 3, 1, 1), oma = c(0, 0, 0, 0))
ylim = c(0.1, 500)
xlim = c(0,15)
xaxis = (fit$stan_data$bin_break_min+fit$stan_data$bin_break_max)/2
for(p in 1:fit$stan_data$N_pairs){
  plot(xaxis, obs_m[p,], log = 'y', ylim = ylim, bty = 'n', yaxt = 'n', xaxt = 'n', pch = 16, col = 'grey40',
       xlab = '', ylab = 'n', xlim = xlim,
       main = fit$set_labels[p])
  arrows(x0 = xaxis, x1 = xaxis, y0 = obs_cimin[p,], y1 = obs_cimax[p,],
         length = 0.05, angle = 90, code = 0, col = 'grey40')
  abline(h=1, lty = 2)
  axis(1, at = seq(0,15,5))
  axis(2, at = c(0.1, 1, 10, 100, 500), labels = c('<0.1', 1, 10, 100, '>500'), las = 2)
  
  p_id = paste0(sort(c(fit$stan_data$pair_a[p], fit$stan_data$pair_b[p])), collapse = '')
  p_id = match(p_id, p_names)
  
  lines(xaxis, pred_m[p,], col = countries_pairs_colors[p_id])
  polygon(x = c(xaxis, rev(xaxis)),
          y = c(pred_cimin[p,], rev(pred_cimax[p,])), border = F, col = adjustcolor(countries_pairs_colors[p_id], alpha.f = 0.2))
}
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
##Figure SXX: tree colored by country 
#########################################################################################################################################
tree = read.tree('4_SC2_geography/Data_countries/Seqs/SC2_geo_ncbi_7_countries.treefile')
tree = phangorn::midpoint(tree) ## to test
tree = ladderize(tree)
metadata = read.csv('4_SC2_geography/Data_countries/Seqs/Table_S1_List_isolates_SC2_geo.csv')
metadata = metadata[match(sub('\\..*', '', tree$tip.label), metadata$Accession),]
countries_iso2 = c('DK', 'JP', 'SK', 'CH', 'TH', 'UK', 'USA')
countries_colors = met.brewer(name = 'Java', n = length(countries_iso2))

countries = metadata$Country
countries = factor(countries)
levels(countries) = countries_colors
countries = as.character(countries)

# pdf(file = '4_SC2_geography/Figures/FigureSXX_Tree_SC2_countries.pdf', width = 19/2.56, height = 25/2.56)
par(mfrow = c(1,1))
plot.phylo(tree, show.tip.label = F, edge.width = 0.25)
tiplabels(pch = 16, col = countries, cex = 0.5)
add.scale.bar('topright', length = 2.5e-4)
# legend('topright', pch = 16, col = countries_colors, legend = c("Denmark", "Japan",  "Slovakia",
                                                                # "Switzerland", "Thailand",  "United Kingdom", "USA"), bty = 'n', title = 'Country')
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: Sensitivity to sampling
#########################################################################################################################################
fit_sensitivity = readRDS('4_SC2_geography/Data_countries/Seqs/Fit_7_countries_0_15_2000_2000_relrate_sensitivity.rds')
hist_sensitivity = readRDS('4_SC2_geography/Data_countries/Seqs/Hists_7_countries_relative_sensitivity.rds')
metadata_raw_sensitivity = readRDS('4_SC2_geography/Data_countries/Seqs/Metadata_raw_7_countries_relative_sensitivity.rds')

# pdf(file = '4_SC2_geography/Figures/FigureSXX_Subsampled_results_13102026.pdf', width = 6, height = 5.5)
par(mfrow = c(2,2), mar = c(2.5,3.5,1,2), oma = c(1,1,1,1))

## Number of sequences through time, and countries - FULL
time_bins = seq(2024.49,2025.51,1/52)
countries = c("Denmark", "Japan",  "Slovakia",
              "Switzerland", "Thailand",  "United Kingdom", "USA")
N_seqs = matrix(NA, ncol = length(countries), nrow = length(time_bins)-1)
for(i in 1:length(countries)){
  N_seqs[,i] = table(cut(metadata_raw$Collection_Date[which(metadata_raw$Country == countries[i])], breaks = time_bins, include.lowest = T, ))
}

## Number of sequences through time, and countries - subsampled
time_bins = seq(2024.49,2025.51,1/52)
N_seqs_subsampled = matrix(NA, ncol = length(countries), nrow = length(time_bins)-1)
for(i in 1:length(countries)){
  N_seqs_subsampled[,i] = table(cut(metadata_raw_sensitivity$Collection_Date[which(metadata_raw_sensitivity$Country == countries[i])], breaks = time_bins, include.lowest = T, ))
}

plot(colSums(N_seqs), type = 'h', ylim = c(0,5000), xlim = c(0.5, 7.5), 
     lwd = 15, col = adjustcolor(countries_colors, alpha.f = 0.15),
     main = 'Number of sequences (2024-2025)', xlab = 'Country', ylab = 'Number', 
     bty = 'n', yaxt = 'n', xaxt = 'n', lend = "butt")
axis(2, las = 2)
axis(1, at = 1:length(countries_iso2), labels = countries_iso2)
lines(colSums(N_seqs_subsampled), type = 'h', ylim = c(0,5000), xlim = c(0.5, 7.5), 
     lwd = 15, col = adjustcolor(countries_colors, alpha.f = 1),
     main = 'Number of sequences (2024-2025)', xlab = 'Country', ylab = 'Number', 
     bty = 'n', yaxt = 'n', xaxt = 'n', lend = "butt")
abline(h = 200, lty = 2)


## Rates full sample vs sub sample
## Raw rates 
idx_steps = c(1:nrow(fit_sensitivity$chains$rates))
rates = fit_sensitivity$chains$rates[idx_steps,]
ts = 365/(substitution_rate_SC2*genome_length_SC2)
# ts = 1
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median_and_ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median_and_ci(rates[,i,])
}

idx_steps = c(1:nrow(fit_sensitivity$chains$rates))
rates_full = fit$chains$rates[idx_steps,]
ts = 365/(substitution_rate_SC2*genome_length_SC2)
# ts = 1
rates_full = apply(rates_full, c(1, 2), median)
res_ci_full = matrix(0, nrow = 3, ncol = dim(rates_full)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates_full)) == 2) res_ci_full[,i] = median_and_ci(rates_full[,i])
  if(length(dim(rates_full)) == 3) res_ci_full[,i] = median_and_ci(rates_full[,i,])
}

plot(res_ci_full[1,]/ts, #[lower.tri(mobility_flow_mat_sym)], 
     res_ci[1,]/ts, pch = 16, 
     ylim = c(1e-4, 1e-1),  xlim = c(1e-4, 1e-1), bty = 'n', yaxt = 'n', log = 'xy', col = countries_pairs_colors, xaxt = 'n',
     xlab = 'Normalized number of trips between countries per year', ylab = 'Rates from genetic distances (per day)')
axis(1, las = 1, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
abline(a = 0, b = 1, lty = 2)
arrows(x0 = res_ci_full[1,]/ts,
       y0 = res_ci[2,]/ts, 
       x1 = res_ci_full[1,]/ts,
       y1 = res_ci[3,]/ts, 
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)
arrows(x0 = res_ci_full[2,]/ts,
       y0 = res_ci[1,]/ts, 
       x1 = res_ci_full[3,]/ts,
       y1 = res_ci[1,]/ts, 
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)


## Rates sub sample vs mobility
## Raw rates
idx_steps = c(1:nrow(fit_sensitivity$chains$rates))
rates = fit_sensitivity$chains$rates[idx_steps,]
ts = 365/(substitution_rate_SC2*genome_length_SC2)
# ts = 1
rates = apply(rates, c(1, 2), median)
res_ci = matrix(0, nrow = 3, ncol = dim(rates)[2])
for(i in 1:dim(rates)[2]){
  if(length(dim(rates)) == 2) res_ci[,i] = median_and_ci(rates[,i])
  if(length(dim(rates)) == 3) res_ci[,i] = median_and_ci(rates[,i,])
}

## Plot against mobility data
mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
mobility_flow_mat_sym[upper.tri(mobility_flow_mat_sym)] = NA
plot(mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))], #[lower.tri(mobility_flow_mat_sym)], 
     res_ci[1,]/ts, pch = 16, 
     ylim = c(1e-4, 1e-1),  xlim = c(1e-4, 0.5), bty = 'n', yaxt = 'n', log = 'xy', col = countries_pairs_colors, xaxt = 'n',
     xlab = 'Normalized number of trips between countries per year', ylab = 'Rates from genetic distances (per day)')
axis(1, at = 5*10^seq(-4,-1,1), labels = format(5*10^seq(-4,-1,1), scientific = T))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
arrows(x0 = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))], #[lower.tri(mobility_flow_mat_sym)], 
       y0 = res_ci[2,]/ts, 
       x1 = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))], #[lower.tri(mobility_flow_mat_sym)], 
       y1 = res_ci[3,]/ts, 
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)
## Compute model fit
slopes = NULL
intercept = NULL
o = order(mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))])
x = mobility_flow_mat_sym[which(!is.na(mobility_flow_mat_sym))][o]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)[o]
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(-10, 2, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the countries
mobility_flow_mat_sym = readRDS('4_SC2_geography/Data_countries/Mobility_between_countries/Mobility_data_between_countries_symmetric_ncbi.rds')
mat = mobility_flow_mat_sym
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value

traits = setNames(metadata_raw_sensitivity$Country, metadata_raw_sensitivity$id)
N_seqs_countries = matrix(0, ncol = length(hist_sensitivity$demes_in_tree), nrow = length(hist_sensitivity$demes_in_tree))
colnames(N_seqs_countries) = hist_sensitivity$demes_in_tree
rownames(N_seqs_countries) = hist_sensitivity$demes_in_tree
for(i in 1:length(hist_sensitivity$demes_in_tree)){
  for(j in 1:length(hist_sensitivity$demes_in_tree)){
    N = length(which(traits == hist_sensitivity$demes_in_tree[i])) + length(which(traits == hist_sensitivity$demes_in_tree[j]))
    N_seqs_countries[i,j] = N
  }
}
diag(N_seqs_countries) = NA
N_seqs_countries_chosen = N_seqs_countries
N_seqs_countries_chosen[upper.tri(N_seqs_countries_chosen)] = NA
plot(N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))],  # [lower.tri(N_seqs_countries_chosen)], 
     res_ci[1,]/ts, pch = 16, xlim = c(300, 400), xaxt = 'n',
     ylim = c(1e-4, 1e-1), bty = 'n', yaxt = 'n', log = 'y',
     col = countries_pairs_colors,
     xlab = 'Nb. sequences per pair', ylab = 'Rates from genetic distances (per day)')
axis(1, at = seq(300, 400, 25))
axis(2, las = 2, at = 10^seq(-4,-1,1), labels = format(10^seq(-4,-1,1), scientific = T))
arrows(x0 = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))], #[lower.tri(N_seqs_countries_chosen)], 
       y0 = res_ci[2,]/ts,
       x1 = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))], #[lower.tri(N_seqs_countries_chosen)], 
       y1 = res_ci[3,]/ts,
       code = 0, length = 0.05, angle = 90, col = countries_pairs_colors)
## Compute model fit
slopes = NULL
intercept = NULL
x = N_seqs_countries_chosen[which(!is.na(N_seqs_countries_chosen))]
for(ii in 1:nrow(rates)){
  y =  (rates[ii,]/ts)
  m <- lm(log(y) ~ log(x))
  slopes = c(slopes, m$coefficients[2])
  intercept = c(intercept, m$coefficients[1])
}
x_new = exp(seq(0, 7, length.out = 1000))
pred_post <- sapply(x_new, function(x) {exp(intercept) * x^slopes})
pred_med <- apply(pred_post, 2, median)
pred_lwr <- apply(pred_post, 2, quantile, probs = 0.025)
pred_upr <- apply(pred_post, 2, quantile, probs = 0.975)
polygon(x = c(x_new, rev(x_new)), y = c(pred_lwr, rev(pred_upr)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
lines(x_new, pred_med, lty = 2)

## Compute p-value: median slope over draws, for every relabelling of the countries
mat = N_seqs_countries
n = nrow(mat)
perms = as.matrix(expand.grid(rep(list(1:n), n)))
perms = perms[apply(perms, 1, function(r) length(unique(r)) == n),]
set.seed(1)
draws = sample(1:nrow(rates), size = 500)
slopes_permt = NULL
Y = log(rates[draws,]/ts)
for(k in 1:nrow(perms)){
  s = perms[k,]
  mat_tmp = mat[s, s]
  mat_tmp[upper.tri(mat_tmp)] = NA
  x = mat_tmp[which(!is.na(mat_tmp))]
  lx = log(x)-mean(log(x))
  slopes_k = as.vector(Y%*%lx)/sum(lx^2)
  slopes_permt = c(slopes_permt, median(slopes_k))
}
print(length(which(slopes_permt>=median(slopes[draws])))/length(slopes_permt)) ## p value

# dev.off()
#########################################################################################################################################
















