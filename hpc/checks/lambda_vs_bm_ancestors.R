library(ape)
t <- read.tree("inputs/tree_codes.nwk"); d <- read.delim("inputs/pheno_loading.tsv")
x <- setNames(d$trait, d$species)[t$tip.label]; n <- length(x)
lam <- 0.849
H <- max(node.depth.edgelength(t)); term <- t$edge[, 2] <= n
lt <- t; lt$edge.length <- t$edge.length * lam; lt$edge.length[term] <- lt$edge.length[term] + (1 - lam) * H
a_bm <- ace(x, t, type = "continuous", method = "REML")$ace
a_l  <- ace(x, lt, type = "continuous", method = "REML")$ace
cat(sprintf("ancestral BM vs lambda: r %.4f | max |diff| %.3f | mean |diff| %.3f | sd BM %.3f sd lambda %.3f\n",
  cor(a_bm, a_l), max(abs(a_bm - a_l)), mean(abs(a_bm - a_l)), sd(a_bm), sd(a_l)))
# branch values under --branch-trait average, all 166 branches
nv_bm <- c(x, a_bm); nv_l <- c(x, a_l)
b_bm <- (nv_bm[t$edge[,1]] + nv_bm[t$edge[,2]]) / 2; b_l <- (nv_l[t$edge[,1]] + nv_l[t$edge[,2]]) / 2
cat(sprintf("branch values BM vs lambda: r %.4f | max |diff| %.3f\n", cor(b_bm, b_l), max(abs(b_bm - b_l))))
