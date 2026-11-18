#########################################################################################################################################
## Exact structured coalescent versus co-location approximation
## m is the pair-level migration rate (each lineage migrates at rate m/2)
## For two contemporaneous sequences, d = 2t
#########################################################################################################################################
library(deSolve)
library(MetBrewer)

colours = met.brewer(name = 'Java', n = 4)

#########################################################################################################################################
## Functions
#########################################################################################################################################
## Effective population size as a function of time: a constant, or a function with its parameters
Ne_function = function(N, params = list()) {
  if(is.function(N)) return(function(t) do.call(N, c(list(t = t), params)))
  return(function(t) rep(N, length(t)))
}

## Exponentially varying effective population size, t backwards from sampling
Nexp = function(t, N0, r) {
  return(N0*exp(r*t))
}

## Normalise a density over the distances where it is shown (trapezoid rule)
normalise_density = function(d, y) {
  ok = !is.na(y)
  area = sum(diff(d[ok])*(head(y[ok], -1)+tail(y[ok], -1))/2)
  return(y/area)
}

## Exact structured coalescent: probabilities that the two lineages have not coalesced and are both in deme 1 (q11),
## in different demes (qD) or both in deme 2 (q22), with coalescence rates c1 = 1/(2 N1(t)) and c2 = 1/(2 N2(t)):
## dq11/dt = -(m+c1) q11 + (m/2) qD
## dqD/dt = m q11 - m qD + m q22
## dq22/dt = -(m+c2) q22 + (m/2) qD
solve_two_deme_SC = function(times, N1, N2, m, start, N1_params = list(), N2_params = list()) {
  N1_fun = Ne_function(N1, N1_params)
  N2_fun = Ne_function(N2, N2_params)
  if(start == '11') y0 = c(q11 = 1, qD = 0, q22 = 0)
  if(start == '22') y0 = c(q11 = 0, qD = 0, q22 = 1)
  if(start == 'D') y0 = c(q11 = 0, qD = 1, q22 = 0)
  derivs = function(t, y, pars) {
    c1 = 1/(2*N1_fun(t))
    c2 = 1/(2*N2_fun(t))
    dq11 = -(m+c1)*y['q11']+(m/2)*y['qD']
    dqD = m*y['q11']-m*y['qD']+m*y['q22']
    dq22 = -(m+c2)*y['q22']+(m/2)*y['qD']
    return(list(c(q11 = dq11,
                  qD = dqD,
                  q22 = dq22)))
  }
  times = sort(unique(c(0, times)))
  out = deSolve::ode(y = y0, times = times, func = derivs, parms = NULL, method = 'lsoda', rtol = 1e-10, atol = 1e-12)
  return(as.data.frame(out))
}

## Distance density of a pair sampled in 11, 22 or different demes (D): g(d) = 1/2 f(d/2), with f(t) = c1(t) q11(t) + c2(t) q22(t)
distance_density_SC = function(d, N1, N2, m, start, N1_params = list(), N2_params = list()) {
  N1_fun = Ne_function(N1, N1_params)
  N2_fun = Ne_function(N2, N2_params)
  t = d/2
  sol = solve_two_deme_SC(t, N1, N2, m, start, N1_params, N2_params)
  q11 = approx(sol$time, sol$q11, xout = t, rule = 2)$y
  q22 = approx(sol$time, sol$q22, xout = t, rule = 2)$y
  c1 = 1/(2*N1_fun(t))
  c2 = 1/(2*N2_fun(t))
  f_t = c1*q11+c2*q22
  return(0.5*f_t)
}

## Exact distance densities, with equal reference probabilities: P(11) = P(22) among same, and P(same) = P(diff) = 1/2
get_exact_distance_densities = function(d, N1, N2, m, N1_params = list(), N2_params = list()) {
  g11 = distance_density_SC(d, N1, N2, m, '11', N1_params, N2_params)
  g22 = distance_density_SC(d, N1, N2, m, '22', N1_params, N2_params)
  g_diff = distance_density_SC(d, N1, N2, m, 'D', N1_params, N2_params)
  g_same = 0.5*g11+0.5*g22
  return(list(g_same = g_same,
              g_diff = g_diff,
              g_overall = 0.5*g_same+0.5*g_diff))
}

## Co-location approximation: P(A_d | same, m) = 1/2 (1+exp(-m d)) and P(A_d | diff, m) = 1/2 (1-exp(-m d))
P_colocated_given_same = function(d, m) {
  return(0.5*(1+exp(-m*d)))
}

P_colocated_given_diff = function(d, m) {
  return(0.5*(1-exp(-m*d)))
}

## Exact and approximate P(d | same) and P(d | diff), the approximation being D(d) P(A_d | same or diff, m)/P(A_d | m), with P(A_d | m) = 1/2
## Not normalised: the exact densities integrate to 1 over all distances, the approximations to 1+delta (same) and 1-delta (diff)
compare_pd_given_same = function(d, N1, N2, m, N1_params = list(), N2_params = list()) {
  exact = get_exact_distance_densities(d, N1, N2, m, N1_params, N2_params)
  return(list(exact = exact$g_same,
              approx = exact$g_overall*(P_colocated_given_same(d, m)/0.5)))
}

compare_pd_given_diff = function(d, N1, N2, m, N1_params = list(), N2_params = list()) {
  exact = get_exact_distance_densities(d, N1, N2, m, N1_params, N2_params)
  return(list(exact = exact$g_diff,
              approx = exact$g_overall*(P_colocated_given_diff(d, m)/0.5)))
}

## Exact and approximate ratio P(d | same)/P(d | diff), the approximation being P(A_d | same, m)/P(A_d | diff, m)
compare_ratio = function(d, N1, N2, m, N1_params = list(), N2_params = list()) {
  exact = get_exact_distance_densities(d, N1, N2, m, N1_params, N2_params)
  return(list(exact = exact$g_same/exact$g_diff,
              approx = P_colocated_given_same(d, m)/P_colocated_given_diff(d, m)))
}

## Exact and approximate P(same | d), with P(same) = P(diff) = 1/2
compare_Psame_given_d = function(d, N1, N2, m, N1_params = list(), N2_params = list()) {
  exact = get_exact_distance_densities(d, N1, N2, m, N1_params, N2_params)
  return(list(exact = exact$g_same/(exact$g_same+exact$g_diff),
              approx = P_colocated_given_same(d, m)))
}

## Exact and approximate P(d | same) (quantity = 'pd_same'), P(d | diff) ('pd_diff'), their ratio ('ratio') or P(same | d) ('Psame'),
## not shown where Ne < 1 in both demes, P(d | same) and P(d | diff) normalised over the distances shown
compare_exact_approx = function(d, N1, N2, m, N1_params, N2_params, quantity) {
  if(quantity == 'pd_same') res = compare_pd_given_same(d, N1, N2, m, N1_params, N2_params)
  if(quantity == 'pd_diff') res = compare_pd_given_diff(d, N1, N2, m, N1_params, N2_params)
  if(quantity == 'ratio') res = compare_ratio(d, N1, N2, m, N1_params, N2_params)
  if(quantity == 'Psame') res = compare_Psame_given_d(d, N1, N2, m, N1_params, N2_params)
  N_max = pmax(Ne_function(N1, N1_params)(d/2), Ne_function(N2, N2_params)(d/2))
  res$exact[N_max < 0.00001] = NA
  res$approx[N_max < 0.00001] = NA
  if(quantity %in% c('pd_same', 'pd_diff')) {
    res$exact = normalise_density(d, res$exact)
    res$approx = normalise_density(d, res$approx)
  }
  return(res)
}

## Plot the exact (colours) and approximate (grey, dashed) curves of each scenario
plot_exact_vs_approx = function(d, dmax, results, quantity) {
  if(quantity == 'Psame') {
    plot(NULL, bty = 'n', yaxt = 'n', xlab = 'Pairwise distance, d', ylab = expression(P(same ~ '|' ~ d)), xlim = c(0, dmax), ylim = c(0, 1))
    axis(2, las = 2, at = seq(0, 1, 0.2))
  }
  if(quantity == 'ratio') {
    y = unlist(lapply(results, function(res) c(res$exact, res$approx)))
    plot(NULL, bty = 'n', yaxt = 'n', log = 'y', xlab = 'Pairwise distance, d', ylab = expression(P(d ~ '|' ~ same)/P(d ~ '|' ~ diff)), xlim = c(0, dmax), ylim = range(c(1, y[is.finite(y) & y > 0])))
    axis(2, las = 2)
    abline(h = 1, lty = 2, col = 'grey50')
  }
  if(quantity %in% c('pd_same', 'pd_diff')) {
    ymax = max(unlist(lapply(results, function(res) c(res$exact, res$approx))), na.rm = T)
    ylab = expression(P(d ~ '|' ~ diff))
    if(quantity == 'pd_same') ylab = expression(P(d ~ '|' ~ same))
    plot(NULL, bty = 'n', yaxt = 'n', xlab = 'Pairwise distance, d', ylab = ylab, xlim = c(0, dmax), ylim = c(0, 1.05*ymax))
    axis(2, las = 2)
  }
  for(i in 1:length(results)) {
    lines(d, results[[i]]$exact, lwd = 2, lty = 1, col = colours[i])
    lines(d, results[[i]]$approx, lwd = 2, lty = 2, col = adjustcolor('grey', alpha.f = 0.75))
  }
  legend('topright', legend = c('Exact structured coalescent', 'Co-location approximation'), lwd = 2, lty = c(1, 2), col = c('grey30', 'grey60'), bty = 'n', cex = 0.65)
}

## Plot the maximum absolute error (relative error for the ratio) between the exact and approximate curves of each scenario
plot_errors = function(x, results, xlab, quantity, xat = NULL) {
  errors = sapply(results, function(res) max(abs(res$exact-res$approx), na.rm = T))
  ylab = 'Max absolute error'
  if(quantity == 'ratio') {
    errors = sapply(results, function(res) max(abs(res$exact-res$approx)/res$approx, na.rm = T))
    ylab = 'Max relative error'
  }
  ylim = c(0, 0.25)
  if(quantity != 'Psame') ylim = c(0, max(0.05, 1.1*max(errors)))
  plot(NULL, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab, xlim = range(x), log = 'x', ylim = ylim)
  if(is.null(xat)) {
    axis(1, las = 1)
  }else{
    axis(1, las = 1, at = xat, labels = format(xat, scientific = F))
  }
  if(quantity == 'Psame') {
    axis(2, las = 2, at = seq(0, 0.25, 0.05))
  }else{
    axis(2, las = 2)
  }
  for(i in 1:length(x)) points(x[i], errors[i], pch = 16, col = colours[i])
  abline(h = 0, lty = 2, col = 'grey50')
}

## Panel A: constant and equal Ne
plot_panel_A = function(N10, m, dmax, quantity) {
  d = seq(0, dmax, length.out = 1000)
  t = d/2
  plot(NULL, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = 'Time before sampling, t', ylab = 'Ne', log = 'y', xlim = range(t), ylim = c(1, max(N10)))
  axis(1, las = 1)
  axis(2, las = 2, at = 10^seq(0, 5, 1), labels = format(10^seq(0, 5, 1), scientific = F))
  for(i in 1:length(N10)) abline(h = N10[i], lty = 1, lwd = 2, col = colours[i])
  results = lapply(1:length(N10), function(i) compare_exact_approx(d, N10[i], N10[i], m, list(), list(), quantity))
  plot_exact_vs_approx(d, dmax, results, quantity)
  plot_errors(N10, results, 'Effective population size', quantity, xat = 10^seq(1, 5, 1))
}

## Panel B: equal demes, exponentially varying Ne
plot_panel_B = function(N10, r10, m, dmax, quantity) {
  d = seq(0, dmax, length.out = 1000)
  t = d/2
  N_curves = lapply(r10, function(r) Nexp(t, N0 = N10, r = r))
  plot(NULL, bty = 'n', xlab = 'Time before sampling, t', ylab = 'Ne', xlim = range(t), ylim = range(unlist(N_curves), finite = T))
  for(i in 1:length(r10)) lines(t, N_curves[[i]], lwd = 2, col = colours[i])
  results = lapply(1:length(r10), function(i) {
    N_params = list(N0 = N10,
                    r = r10[i])
    return(compare_exact_approx(d, Nexp, Nexp, m, N_params, N_params, quantity))
  })
  plot_exact_vs_approx(d, dmax, results, quantity)
  plot_errors(abs(r10), results, 'Rate of change in effective population size', quantity)
}

## Panel C: deme 1 constant, deme 2 changing
plot_panel_C = function(N10, N20, r20, m, dmax, quantity) {
  d = seq(0, dmax, length.out = 1000)
  t = d/2
  N1_curve = rep(N10, length(t))
  N2_curves = lapply(r20, function(r) Nexp(t, N0 = N20, r = r))
  plot(NULL, bty = 'n', xlab = 'Time before sampling, t', ylab = 'Ne', xlim = range(t), ylim = range(c(N1_curve, unlist(N2_curves)), finite = T))
  lines(t, N1_curve, lwd = 2, col = 'black')
  for(i in 1:length(r20)) lines(t, N2_curves[[i]], lwd = 2, col = colours[i])
  results = lapply(1:length(r20), function(i) {
    N2_params = list(N0 = N20,
                     r = r20[i])
    return(compare_exact_approx(d, N10, Nexp, m, list(), N2_params, quantity))
  })
  plot_exact_vs_approx(d, dmax, results, quantity)
  plot_errors(abs(r20), results, 'Rate of change in deme 2 Ne', quantity)
}

## Panel D: demes changing in opposite directions
plot_panel_D = function(N10, N20, r10, r20, m, dmax, quantity) {
  d = seq(0, dmax, length.out = 1000)
  t = d/2
  N1_curves = lapply(r10, function(r) Nexp(t, N0 = N10, r = r))
  N2_curves = lapply(r20, function(r) Nexp(t, N0 = N20, r = r))
  plot(NULL, bty = 'n', xlab = 'Time before sampling, t', ylab = 'Ne', xlim = range(t), ylim = range(c(unlist(N1_curves), unlist(N2_curves)), finite = T))
  lines(t, N1_curves[[1]], lwd = 2, col = 'black')
  for(i in 1:length(r20)) lines(t, N2_curves[[i]], lwd = 2, col = colours[i])
  results = lapply(1:length(r20), function(i) {
    N1_params = list(N0 = N10,
                     r = r10[i])
    N2_params = list(N0 = N20,
                     r = r20[i])
    return(compare_exact_approx(d, Nexp, Nexp, m, N1_params, N2_params, quantity))
  })
  plot_exact_vs_approx(d, dmax, results, quantity)
  plot_errors(abs(r20), results, 'Rate of change in deme 2 Ne', quantity)
}
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: comparison of P(d | diff)
#########################################################################################################################################
# pdf(file = '2_theory/Figure_method/Figure_approximation_pddiff_structured_coalescent_15102026.pdf', width = 6, height = 8)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
plot_panel_A(N10 = c(1e2, 1e3, 1e4, 1e5), m = 0.1, dmax = 25, quantity = 'pd_diff')
plot_panel_B(N10 = 1e4, r10 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_diff')
plot_panel_C(N10 = 1e4, N20 = 1e4, r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_diff')
plot_panel_D(N10 = 1e4, N20 = 1e4, r10 = rep(0.1, 5), r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_diff')
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: comparison of P(d | same)
#########################################################################################################################################
# pdf(file = '2_theory/Figure_method/Figure_approximation_pdsame_structured_coalescent_15102026.pdf', width = 6, height = 8)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
plot_panel_A(N10 = c(1e2, 1e3, 1e4, 1e5), m = 0.1, dmax = 25, quantity = 'pd_same')
plot_panel_B(N10 = 1e4, r10 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_same')
plot_panel_C(N10 = 1e4, N20 = 1e4, r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_same')
plot_panel_D(N10 = 1e4, N20 = 1e4, r10 = rep(0.1, 5), r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'pd_same')
# dev.off()
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: comparison of the ratio P(d | same)/P(d | diff)
#########################################################################################################################################
# pdf(file = '2_theory/Figure_method/Figure_approximation_ratio_structured_coalescent_15102026.pdf', width = 6, height = 8)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
plot_panel_A(N10 = c(1e2, 1e3, 1e4, 1e5), m = 0.1, dmax = 25, quantity = 'ratio')
plot_panel_B(N10 = 1e4, r10 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'ratio')
plot_panel_C(N10 = 1e4, N20 = 1e4, r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'ratio')
plot_panel_D(N10 = 1e4, N20 = 1e4, r10 = rep(0.1, 5), r20 = -c(0.01, 0.1, 0.5, 1), m = 0.1, dmax = 25, quantity = 'ratio')
# dev.off()
#########################################################################################################################################
