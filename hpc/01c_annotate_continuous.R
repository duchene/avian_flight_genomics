# 01c_annotate_continuous.R : annotate the tree with the continuous loading index.
#
# `pelican scan continuous` needs a trait value at EVERY node, because it assigns each
# branch a value from its two ends (--branch-trait). PELICAN's own tree-annotation only
# reconstructs discrete states, so internal nodes are estimated here instead.
#
# Ancestral values are the maximum likelihood estimates under Brownian motion (ape::ace,
# REML), fitted on the DATED tree, because Brownian variance accrues with time rather than
# with substitutions. The values are then written onto the substitution-scaled tree,
# which has the same topology and node order, because that is what PELICAN's model reads.
#
# Usage:  Rscript 01c_annotate_continuous.R [phenotype]     # default loading

library(ape)
ph <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(ph)) ph <- "loading"

dated <- read.tree("inputs/tree_codes.nwk")
subst <- read.tree("inputs/tree_subst.nwk")
stopifnot(identical(dated$tip.label, subst$tip.label),
          identical(dated$edge, subst$edge))

d <- read.delim(sprintf("inputs/pheno_%s.tsv", ph), stringsAsFactors = FALSE)
x <- setNames(as.numeric(d$trait), d$species)[dated$tip.label]
stopifnot(!anyNA(x))

fit <- ace(x, dated, type = "continuous", method = "REML")
node_val <- c(x, fit$ace)          # ape numbering: tips 1..n, then internal nodes
cat(sprintf("%s: %d tips, %d internal nodes, sigma2 %.4g\n",
            ph, Ntip(dated), dated$Nnode, fit$sigsq))
cat("tip range", sprintf("%.3f", range(x)), "| ancestral range",
    sprintf("%.3f", range(fit$ace)), "| root", sprintf("%.3f", fit$ace[1]), "\n")

# NHX writer. Every node carries [&&NHX:Trait=value], and every branch a length.
nhx <- function(node) {
  kids <- subst$edge[subst$edge[, 1] == node, 2]
  lab  <- if (node <= Ntip(subst)) subst$tip.label[node] else
    paste0("(", paste(vapply(kids, nhx, ""), collapse = ","), ")")
  e <- which(subst$edge[, 2] == node)
  len <- if (length(e)) sprintf(":%.9g", subst$edge.length[e]) else ""
  sprintf("%s%s[&&NHX:Trait=%.6g]", lab, len, node_val[node])
}
out <- sprintf("work/tree_%s.nhx", ph)
dir.create("work", showWarnings = FALSE)
writeLines(paste0(nhx(Ntip(subst) + 1), ";"), out)
cat("written", out, "\n")
