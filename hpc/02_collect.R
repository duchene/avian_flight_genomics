# 02_collect.R : gather the array output, control the false discovery rate, report.
#
# PELICAN writes one all_sites.tsv per run with these columns:
#
#   alignment site nseq naa  <prefilter>_pval  aagtr_pval  reduced_loglk full_loglk  tag
#
# aagtr_pval is the test that matters. It compares a model in which amino acid preferences
# at the site differ between conditions against one in which they do not. The prefilter
# column is multinomial_pval for the discrete test and anova_pval for the continuous one,
# and it exists only to avoid fitting the full model everywhere.
#
# TWO LEVELS OF CORRECTION, AND THEY ANSWER DIFFERENT QUESTIONS.
#   per site  Benjamini-Hochberg across every site in the screen. Answers "which sites".
#   per gene  the smallest p-value in a gene, Sidak-corrected for its length, then BH
#             across genes. Answers "which genes", which is what the gene-set work needs.
#
# Usage:  Rscript 02_collect.R [phenotype]     # default flight

ph <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(ph)) ph <- "flight"
root <- file.path("results", ph)
stopifnot(dir.exists(root))

f <- list.files(root, pattern = "^all_sites\\.tsv$", recursive = TRUE, full.names = TRUE)
cat("phenotype:", ph, "| all_sites.tsv files found:", length(f), "\n")
stopifnot(length(f) > 0)

d <- do.call(rbind, lapply(f, function(x)
  tryCatch(read.delim(x, stringsAsFactors = FALSE), error = function(e) NULL)))
cat("rows:", nrow(d), "| genes:", length(unique(d$alignment)), "\n")

pre <- intersect(c("multinomial_pval", "anova_pval"), names(d))[1]
cat("pre-filter column:", pre, "\n")
if ("tag" %in% names(d)) { cat("\nsite status tags:\n"); print(table(d$tag, useNA = "ifany")) }

# Sites the pre-filter never passed carry no aagtr_pval. They are not missing data, they
# are sites the screen declined to test, so they leave the denominator rather than being
# imputed. The count is reported because it defines what the correction is over.
d$aagtr_pval <- suppressWarnings(as.numeric(d$aagtr_pval))
tested <- d[is.finite(d$aagtr_pval), ]
cat("\nsites tested in full:", nrow(tested), "of", nrow(d),
    sprintf(" (%.1f%%)\n", 100 * nrow(tested) / nrow(d)))

# ---- per-site ----
tested$q_site <- p.adjust(tested$aagtr_pval, method = "BH")
cat("\nPER SITE, BH across all tested sites\n")
for (a in c(0.01, 0.05, 0.1))
  cat(sprintf("   q < %-5.2f : %6d sites in %4d genes\n", a,
              sum(tested$q_site < a),
              length(unique(tested$alignment[tested$q_site < a]))))

# ---- per-gene ----
# Sidak over the sites actually tested in that gene, which is the correct multiplicity for
# taking a minimum, then BH across genes.
g <- do.call(rbind, lapply(split(tested, tested$alignment), function(x) {
  k <- nrow(x); pmin_ <- min(x$aagtr_pval)
  data.frame(alignment = x$alignment[1], n_sites_tested = k,
             best_site = x$site[which.min(x$aagtr_pval)],
             p_best = pmin_, p_sidak = 1 - (1 - pmin_)^k, stringsAsFactors = FALSE)
}))
g$q_gene <- p.adjust(g$p_sidak, method = "BH")
g <- g[order(g$q_gene, g$p_sidak), ]
cat("\nPER GENE, Sidak within gene then BH across genes\n")
for (a in c(0.01, 0.05, 0.1))
  cat(sprintf("   q < %-5.2f : %4d genes\n", a, sum(g$q_gene < a)))

dir.create("results/summary", showWarnings = FALSE, recursive = TRUE)
write.csv(g, sprintf("results/summary/pelican_%s_per_gene.csv", ph), row.names = FALSE)
write.csv(tested[tested$q_site < 0.1, ],
          sprintf("results/summary/pelican_%s_sites_q10.csv", ph), row.names = FALSE)

cat("\nTop 20 genes:\n\n")
print(format(head(g, 20), digits = 3), row.names = FALSE)

# Calibration. A screen with no signal should give a uniform p-value distribution, so a
# median far from 0.5 means the test is miscalibrated rather than that biology is present.
cat(sprintf("\nCalibration: median aagtr_pval %.3f (0.5 expected under the null), ",
            median(tested$aagtr_pval)))
cat(sprintf("fraction below 0.05 %.4f (0.05 expected)\n",
            mean(tested$aagtr_pval < 0.05)))
cat("\nwritten results/summary/\n")
