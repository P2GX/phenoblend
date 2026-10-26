#!/usr/bin/env Rscript
# QQ and PP plots of the shared-term enrichment p-values from phenoblend-batch.
#
# usage: Rscript qq_plot.R results.tsv [shuffles_directory] [output.pdf|output.png]
#
#  results.tsv          table written by phenoblend-batch (rows with mode "observed")
#  shuffles_directory   optional: tables from runs with different --shuffle-seed; each gives one
#                       grey curve, together an empirical band for "no relation between an
#                       individual and the genes"
#  output               default qq_plot.pdf
#
# Only individuals with exactly two genes and one two-gene result are used. The p-values are
# randomized, p = P(X > k) + U * P(X = k) with U uniform, which makes them exactly uniform
# under the null even though the exact test is discrete.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript qq_plot.R results.tsv [shuffles_directory] [output]")
results_file <- args[1]
shuffle_dir  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NA
out_file     <- if (length(args) >= 3) args[3] else "qq_plot.pdf"

read_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]            # the provenance comments
  read.delim(text = lines, quote = "", comment.char = "", na.strings = "", fill = TRUE,
             stringsAsFactors = FALSE)
}

# one row per phenopacket: its two-gene result
two_gene_rows <- function(df, mode) {
  d <- df[df$mode == mode & df$status == "ok" & df$n_genes == 2, ]
  n_rows <- table(d$phenopacket_id)
  d[d$phenopacket_id %in% names(n_rows)[n_rows == 1], ]
}

randomized_p <- function(d, seed) {
  set.seed(seed)
  d$p_greater + runif(nrow(d)) * d$p_equal
}

real <- two_gene_rows(read_table(results_file), "observed")
u_real <- sort(randomized_p(real, 1))
n <- length(u_real)
if (n < 2) stop("no individuals with exactly two genes and mode 'observed' were found in ", results_file)

shuffled <- list()
if (!is.na(shuffle_dir)) {
  files <- Sys.glob(file.path(shuffle_dir, "*.tsv"))
  for (i in seq_along(files)) {
    d <- two_gene_rows(read_table(files[i]), "shuffled")
    if (nrow(d) > 0) shuffled[[length(shuffled) + 1]] <- sort(randomized_p(d, 100 + i))
  }
} else {
  d <- two_gene_rows(read_table(results_file), "shuffled")
  if (nrow(d) > 0) shuffled[[1]] <- sort(randomized_p(d, 100))
}

ks <- ks.test(u_real, "punif", alternative = "greater")   # more small p-values than chance
cat(sprintf("individuals with two genes: %d\n", n))
cat(sprintf("exact p < 0.05: %d (expected at most %.1f)\n",
            sum(real$p_value < 0.05), 0.05 * n))
cat(sprintf("exact p < 0.01: %d (expected at most %.1f)\n",
            sum(real$p_value < 0.01), 0.01 * n))
cat(sprintf("one-sided KS test of the randomized p-values: D+ = %.3f, p = %.2g\n",
            unname(ks$statistic), ks$p.value))
cat(sprintf("shuffled controls drawn: %d\n", length(shuffled)))

purple <- "#8D10CB"; grey <- adjustcolor("grey40", alpha.f = 0.35)
expected <- function(n) (seq_len(n) - 0.5) / n

if (grepl("\\.png$", out_file)) png(out_file, width = 2200, height = 1100, res = 220) else pdf(out_file, width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.1)

# ---- QQ plot on the -log10 scale: the small p-values, where enrichment shows ----
ex <- expected(n)
lower <- qbeta(0.025, seq_len(n), n - seq_len(n) + 1)
upper <- qbeta(0.975, seq_len(n), n - seq_len(n) + 1)
top <- max(-log10(u_real), -log10(ex)) * 1.05
plot(NA, xlim = c(0, max(-log10(ex))), ylim = c(0, top),
     xlab = expression("expected " * -log[10](p)), ylab = expression("observed " * -log[10](p)),
     main = "QQ plot", las = 1)
polygon(c(-log10(ex), rev(-log10(ex))), c(-log10(upper), rev(-log10(lower))),
        col = adjustcolor("grey80", alpha.f = 0.6), border = NA)
for (u in shuffled) lines(-log10(expected(length(u))), -log10(u), col = grey, lwd = 1)
abline(0, 1, col = "grey30", lty = 2)
points(-log10(ex), -log10(u_real), pch = 19, col = purple, cex = 0.8)
legend("topleft", bty = "n", cex = 0.85,
       legend = c(sprintf("observed (n = %d)", n),
                  if (length(shuffled)) sprintf("%d shuffled controls", length(shuffled)),
                  "95% band under the null"),
       col = c(purple, if (length(shuffled)) grey, "grey80"),
       pch = c(19, if (length(shuffled)) NA, 15), lty = c(0, if (length(shuffled)) 1, 0),
       pt.cex = c(0.8, if (length(shuffled)) NA, 1.4))

# ---- PP plot: cumulative fraction of individuals against the p-value ----
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "p-value (randomized)",
     ylab = "cumulative fraction of individuals", main = "PP plot", las = 1)
for (u in shuffled) lines(u, expected(length(u)), col = grey, lwd = 1)
abline(0, 1, col = "grey30", lty = 2)
lines(u_real, ex, col = purple, lwd = 2)
legend("bottomright", bty = "n", cex = 0.85,
       legend = sprintf("one-sided KS: D+ = %.2f, p = %.2g", unname(ks$statistic), ks$p.value))
invisible(dev.off())
cat("written:", out_file, "\n")