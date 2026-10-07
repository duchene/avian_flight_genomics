library(ape)
t <- read.tree("inputs/tree_codes.nwk"); d <- read.delim("inputs/pheno_loading.tsv")
x <- setNames(d$trait, d$species)[t$tip.label]; n <- length(x); C <- vcv(t)
ll <- function(lam) {           # REML-free ML for an intercept-only Gaussian model
  V <- C * lam; diag(V) <- diag(C); Vi <- solve(V); o <- rep(1, n)
  mu <- sum(Vi %*% x) / sum(Vi); r <- x - mu; s2 <- drop(t(r) %*% Vi %*% r) / n
  -0.5 * (n * log(2 * pi * s2) + determinant(V)$modulus + n)
}
op <- optimize(ll, c(0, 1), maximum = TRUE)
lr0 <- 2 * (op$objective - ll(0)); lr1 <- 2 * (op$objective - ll(1))
cat(sprintf("lambda %.3f | logL %.2f | vs lambda=0: LR %.2f p %.2g | vs lambda=1 (BM): LR %.2f p %.2g\n",
  op$maximum, op$objective, lr0, pchisq(lr0, 1, lower = FALSE), lr1, pchisq(lr1, 1, lower = FALSE)))
# How much do ancestral values move if BM is replaced by the lambda-scaled tree?
tl <- t; int <- t$edge[, 2] > n
H <- node.depth.edgelength(t); tl$edge.length[!int] <- t$edge.length[!int] + 0  # placeholder
a1 <- ace(x, t, type = "continuous", method = "REML")$ace
lt <- t; lt$edge.length <- t$edge.length * op$maximum
tipd <- max(H) - H[t$edge[!int, 1]]; lt$edge.length[!int] <- lt$edge.length[!int] + (1 - op$maximum) * max(H) * 0 + t$edge.length[!int] * 0
cat("BM ancestral values: range", sprintf("%.2f", range(a1)), "\n")
