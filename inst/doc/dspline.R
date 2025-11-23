## ----include=FALSE------------------------------------------------------------
knitr::opts_knit$set(global.par = TRUE)
knitr::opts_chunk$set(dev = "png", dpi = 300)

## ----include=FALSE------------------------------------------------------------
par(mar = c(4.5, 4.5, 0.5, 0.5))

## -----------------------------------------------------------------------------
library(dspline)

## -----------------------------------------------------------------------------
expr = expression(sin(x * (1 + x) * 7 * pi / 2) + (1 - 4 * (x - 0.5)^2) * 4)
f = function(x) eval(expr)
curve(f, n = 1000)

## -----------------------------------------------------------------------------
n = 100
xd = 1:n / (n+1)
x = seq(0, 1, length = 1000)

# First derivatives
fd = D(expr, "x")
u1 = as.numeric(eval(fd))
vhat1 = d_mat_mult(f(xd), 1, xd)

plot(x, u1, type = "l", ylab = "First derivative")
points(xd[2:n], vhat1, col = 2, type = "o")
legend("bottomleft", lty = 1, pch = c(NA, 21), col = 1:2,
       legend = c("Exact", "Discrete"))
rug(xd)

# Second derivatives
fd2 = D(fd, "x")
u2 = as.numeric(eval(fd2))
v2 = d_mat_mult(f(xd), 2, xd)

plot(x, u2, type = "l", ylab = "Second derivative")
points(xd[3:n], v2, col = 2, type = "o")
legend("bottomleft", legend = c("Exact", "Discrete"),
       lty = 1, pch = c(NA, 21), col = 1:2)
rug(xd)

## -----------------------------------------------------------------------------
y = f(xd) + rnorm(n, sd = 0.5)

knot_idx = 1:9 * 10 
res = dspline_solve(y, 3, xd, knot_idx)
yhat = res$fit # Fitted values from the regression of y onto discrete splines
N = res$mat # Discrete B-spline basis, for visualization purposes only!

plot(xd, y, xlab = "x")
points(xd, yhat, col = 2, pch = 19)
matplot(xd, N * 2, type = "l", lty = 1, add = TRUE)
rug(xd)

## -----------------------------------------------------------------------------
x = seq(0, 1, length = 1000)
fhat = dspline_interp(yhat, 3, xd, x, implicit = FALSE)

plot(xd, y, xlab = "x")
lines(x, fhat, col = 2, lwd = 2)
matplot(xd, N * 2, type = "l", lty = 1, add = TRUE)
abline(v = xd[knot_idx], col = 8)
rug(xd)

## -----------------------------------------------------------------------------
dd1 = d_mat_mult(yhat, 1, xd)
dd2 = d_mat_mult(yhat, 2, xd)

plot(xd[2:n], dd1, xlab = "x", ylab = "First discrete derivative")
abline(v = xd[knot_idx], col = 8)
rug(xd)

plot(xd[3:n], dd2, xlab = "x", ylab = "Second discrete derivative")
abline(v = xd[knot_idx], col = 8)
rug(xd)

## -----------------------------------------------------------------------------
inds = x > 0.05 # Exclude points near left boundary
x = x[inds]
fhat = fhat[inds]
dd3 = discrete_deriv(c(yhat, fhat), 3, xd, x)

plot(x, dd3, type = "l", ylab = "Third discrete derivative")
abline(v = xd[knot_idx], col = 8)
rug(xd[xd > min(x)])

## -----------------------------------------------------------------------------
res = dspline_solve(y, 3, xd, knot_idx, basis = "H")
yhat = res$fit # Fitted values from the regression of y onto discrete splines
H = res$mat # Falling factorial basis, for visualization purposes only!
sol = res$sol # Falling factorial basis coefficients, for expansion later

# Sanity check: the fitted values from H (instead of N) should look as before 
plot(xd, y, xlab = "x")
points(xd, yhat, col = 2, pch = 19)
matplot(xd, H * 80, type = "l", lty = 1, add = TRUE)
rug(xd)

# Now build analytic expansion in terms of falling factorial bases functions.
# Unfortunately we need to separate out the terms involving inequalities, since
# D() can't properly (symbolically) differentiate through inequalities, so this
# ends up being more complicated than it should be ...
poly_terms = sprintf("%f", sol[1])
ineq_terms = "1"
for (j in 2:length(sol)) {
  if (j <= 4) {
    x_prod = paste(
      sprintf("(x - %f)", xd[1:(j-1)]), 
      collapse = " * ")
    poly_terms = c(
      poly_terms, 
      sprintf("%f * %s / factorial(%i)", sol[j], x_prod, j-1))
    ineq_terms = c(ineq_terms, "1")
  }
  else {
    x_prod = paste(
      sprintf("(x - %f)", xd[knot_idx[j-4] - 2:0]), 
      collapse = " * ")
    poly_terms = c(
      poly_terms, 
      sprintf("%f * %s / factorial(3)", sol[j], x_prod))
    ineq_terms = c(
      ineq_terms, 
      sprintf("(x > %f)", xd[knot_idx[j-4]]))
  }
}

# Sanity check: the interpolant from this expression should look as before 
combined_terms = paste(
  paste(poly_terms, ineq_terms, sep = " * "),
  collapse = " + ")
fhat = eval(str2expression(combined_terms))
plot(xd, y, xlab = "x")
lines(x, fhat, col = 2, lwd = 2)
matplot(xd, N * 2, type = "l", lty = 1, add = TRUE)
abline(v = xd[knot_idx], col = 8)
rug(xd)

# Higher-order symbolic derivative function (from help file for D())
DD = function(expr, name, order = 1) {
   if(order == 1) D(expr, name)
   else DD(D(expr, name), name, order - 1)
}

# Finally, compute third derivatives of falling factorial basis expansion
d3 = numeric(length(x))
for (i in 1:length(x)) {
  if (x[i] <= xd[knot_idx[1]]) {
    expr = str2expression(paste(poly_terms[1:4], collapse = " + "))
  }
  else {
    j = max(which(x[i] > xd[knot_idx])) + 4
    expr = str2expression(paste(poly_terms[1:j], collapse = " + "))
  }
  d3[i] = eval(DD(expr, "x", 3), list(x = x[i]))
}

plot(x, d3, type = "l", ylab = "Third derivative")
lines(x, dd3, col = 2, lty = 2, lwd = 2)
abline(v = xd[knot_idx], col = 8)
legend("bottomleft", lty = 1:2, col = 1:2, 
       legend = c("Exact", "Discrete"))
rug(xd[xd > min(x)])

