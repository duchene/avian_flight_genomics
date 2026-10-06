# 01b_prune_trees.R : one tree per gene, pruned to that gene's taxa.
#
# Only needed in PERGENE mode, which 00_setup.sh selects when PELICAN refuses a tree
# carrying taxa the alignment lacks. 7,087 of the 9,890 eligible genes have a taxon set no
# other gene shares, so there is nothing to be gained by grouping them.
#
# Needs ape:  conda install -c conda-forge r-ape

library(ape)
phy <- read.tree("inputs/tree_codes.nwk")
genes <- readLines("inputs/eligible_genes.txt")
dir.create("inputs/pruned_trees", showWarnings = FALSE, recursive = TRUE)

cat("pruning", length(genes), "trees from a", Ntip(phy), "tip tree\n")
n_ok <- 0
for (g in genes) {
  base <- sub("\\..*$", "", g)
  f <- file.path("inputs/pruned_trees", paste0(base, ".nwk"))
  if (file.exists(f)) { n_ok <- n_ok + 1; next }
  l <- readLines(file.path("alignments", g), warn = FALSE)
  tax <- sub("^>", "", l[startsWith(l, ">")])
  drop <- setdiff(phy$tip.label, tax)
  # Fewer than four tips leaves nothing for a phylogenetic test to work with.
  if (Ntip(phy) - length(drop) < 4) next
  write.tree(if (length(drop)) drop.tip(phy, drop) else phy, f)
  n_ok <- n_ok + 1
  if (n_ok %% 1000 == 0) cat("  ", n_ok, "\n")
}
cat("written", n_ok, "pruned trees into inputs/pruned_trees/\n")
