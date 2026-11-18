#########################################################################################################################################
## Exact structured coalescent with asymmetric migration versus co-location approximation with the average rate (a + b)/2
## a: rate at which a lineage in deme 1 moves to deme 2, b: from deme 2 to deme 1 (per lineage, backwards in time), a + b fixed
## For two contemporaneous sequences, d = 2t
#########################################################################################################################################
library(deSolve)
library(MetBrewer)

colours = met.brewer(name = 'Java', n = 4)
asym = c(0.5, 0.75, 0.875, 1)
ab_sum = 0.1

#########################################################################################################################################
## Functions
#########################################################################################################################################
## Exponentially varying effective population size, t backwards from sampling
Nexp = function(t, N0, r) {
  return(N0*exp(r*t))
}

## Normalise a density over the distances where it is shown (trapezoid rule)
normalise_density = function(d, y) {
  area = sum(diff(d)*(head(y, -1)+tail(y, -1))/2)
  return(y/area)
}

## Exact structured coalescent: probabilities that the two lineages have not coalesced and are both in deme 1 (q11),
## in different demes (qD) or both in deme 2 (q22), with coalescence rates c1 = 1/(2 N1(t)) and c2 = 1/(2 N2(t)):
## dq11/dt = -(2a+c1) q11 + b qD
## dqD/dt = 2a q11 - (a+b) qD + 2b q22
## dq22/dt = -(2b+c2) q22 + a qD
## Distance density of a pair starting in 11, 22 or different demes (D): g(d) = 1/2 f(d/2), with f(t) = c1(t) q11(t) + c2(t) q22(t)
distance_density_SC = function(d, N1, N2, a, b, start) {
  derivs = function(t, y, pars) {
    c1 = 1/(2*N1(t))
    c2 = 1/(2*N2(t))
    return(list(c(-(2*a+c1)*y[1]+b*y[2],
                  2*a*y[1]-(a+b)*y[2]+2*b*y[3],
                  -(2*b+c2)*y[3]+a*y[2])))
  }
  y0 = c(q11 = 0, qD = 0, q22 = 0)
  y0[start] = 1
  t = d/2
  sol = as.data.frame(deSolve::ode(y = y0, times = sort(unique(c(0, t))), func = derivs, parms = NULL, method = 'lsoda', rtol = 1e-10, atol = 1e-12))
  q11 = approx(sol$time, sol$q11, xout = t, rule = 2)$y
  q22 = approx(sol$time, sol$q22, xout = t, rule = 2)$y
  return(0.5*(q11/(2*N1(t))+q22/(2*N2(t))))
}

## Co-location approximation: P(A_d | diff) / P(A_d), with P(A_d) = P(A_d | 11)/4 + P(A_d | 22)/4 + P(A_d | diff)/2, for rates a and b
cr_colocation = function(d, a, b) {
  e = exp(-(a+b)*d/2)
  P11 = (b+a*e)/(a+b)
  P12 = a*(1-e)/(a+b)
  P21 = b*(1-e)/(a+b)
  P22 = (a+b*e)/(a+b)
  p11 = P11^2+P12^2
  p22 = P21^2+P22^2
  pD = P11*P21+P12*P22
  return(pD/(p11/4+p22/4+pD/2))
}

## Exact P(d | diff) and its approximation D(d) cr(d), with D(d) the exact distance density of the reference (pairs 11, 22 and D weighted 1/4, 1/4, 1/2)
## and cr computed with the average rate (a + b)/2 or with the true rates, all normalised over the distances shown
compare_pd_diff = function(d, N1, N2, a, b) {
  g11 = distance_density_SC(d, N1, N2, a, b, 'q11')
  g22 = distance_density_SC(d, N1, N2, a, b, 'q22')
  gD = distance_density_SC(d, N1, N2, a, b, 'qD')
  D = g11/4+g22/4+gD/2
  return(list(exact = normalise_density(d, gD),
              approx_average = normalise_density(d, D*cr_colocation(d, (a+b)/2, (a+b)/2)),
              approx_true = normalise_density(d, D*cr_colocation(d, a, b))))
}

## One row: Ne of both demes through time, exact (colours) and approximate P(d | diff) with the average rate (grey, dashed), and max absolute errors
plot_row = function(N1, N2, dmax, legend_asym = F) {
  d = seq(0, dmax, length.out = 1000)
  t = d/2
  plot(NULL, bty = 'n', xlab = 'Time before sampling, t', ylab = 'Ne', xlim = range(t), ylim = range(c(0, N1(t), N2(t))))
  lines(t, N1(t), lwd = 2, col = 'black')
  lines(t, N2(t), lwd = 2, lty = 2, col = 'grey40')
  if(legend_asym) legend('bottomright', legend = c('Deme 1', 'Deme 2'), lwd = 2, lty = c(1, 2), col = c('black', 'grey40'), bty = 'n', cex = 0.65)

  results = lapply(asym, function(f) compare_pd_diff(d, N1, N2, f*ab_sum, (1-f)*ab_sum))
  ymax = max(unlist(results))
  plot(NULL, bty = 'n', yaxt = 'n', xlab = 'Pairwise distance, d', ylab = expression(P(d ~ '|' ~ diff)), xlim = c(0, dmax), ylim = c(0, 1.05*ymax))
  axis(2, las = 2)
  for(i in 1:length(asym)) {
    lines(d, results[[i]]$exact, lwd = 2, col = colours[i])
    lines(d, results[[i]]$approx_average, lwd = 2, lty = 2, col = adjustcolor('grey', alpha.f = 0.75))
  }
  if(legend_asym) legend('topleft', legend = c('Exact structured coalescent', 'Co-location, average rate'), lwd = 2, lty = c(1, 2), col = c('grey30', 'grey60'), bty = 'n', cex = 0.65)
  if(legend_asym) legend('bottomright', legend = paste0('a/(a+b) = ', asym), lwd = 2, col = colours, bty = 'n', cex = 0.65)

  err_average = sapply(results, function(res) max(abs(res$exact-res$approx_average)))
  err_true = sapply(results, function(res) max(abs(res$exact-res$approx_true)))
  plot(NULL, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = 'Asymmetry, a/(a+b)', ylab = 'Max absolute error', xlim = c(0.45, 1.05), ylim = c(0, max(0.05, 1.1*max(c(err_average, err_true)))))
  axis(1, at = asym, labels = asym)
  axis(2, las = 2)
  points(asym, err_average, pch = 16, col = colours)
  points(asym, err_true, pch = 1, col = colours)
  abline(h = 0, lty = 2, col = 'grey50')
  if(legend_asym) legend('topleft', legend = c('Average rate', 'True rates'), pch = c(16, 1), col = 'grey30', bty = 'n', cex = 0.65)
}
#########################################################################################################################################

#########################################################################################################################################
## Figure SXX: comparison of P(d | diff) with asymmetric migration
#########################################################################################################################################
# pdf(file = '2_theory/Figure_method/Figure_approximation_pddiff_asymmetry_102026.pdf', width = 6, height = 8)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
plot_row(function(t) Nexp(t, 1e4, 0), function(t) Nexp(t, 1e4, 0), dmax = 25, legend_asym = T)
plot_row(function(t) Nexp(t, 1e4, -0.5), function(t) Nexp(t, 1e4, -0.5), dmax = 25)
plot_row(function(t) Nexp(t, 1e4, 0), function(t) Nexp(t, 1e4, -0.5), dmax = 25)
plot_row(function(t) Nexp(t, 1e4, 0.1), function(t) Nexp(t, 1e4, -0.5), dmax = 25)
# dev.off()
#########################################################################################################################################
