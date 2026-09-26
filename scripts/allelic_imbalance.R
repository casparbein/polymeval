library(VGAM)
library(tidyverse)
library(data.table)
library(RColorBrewer)

safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

query_files  <- snakemake@input[["query"]]
depth_files  <- snakemake@input[["depth"]]
sites_file   <- snakemake@input[["sites"]]
gc_file      <- snakemake@input[["gc"]]
sample_names <- unlist(strsplit(as.character(snakemake@params[["sample_names"]]), ","))
in_colors    <- snakemake@params[["colors"]]
min_dp       <- as.integer(snakemake@params[["min_dp"]])
min_gq       <- as.integer(snakemake@params[["min_gq"]])
out_vaf      <- snakemake@output[["vaf"]]
out_dropout  <- snakemake@output[["dropout"]]
out_sig      <- snakemake@output[["sig"]]
out_plots    <- snakemake@output[["plots"]]

if (length(in_colors) == 1 && nzchar(in_colors) && file.exists(in_colors)) {
  col_dict <- read_delim(in_colors, col_names = FALSE, show_col_types = FALSE)
  custom_colors <- setNames(col_dict$X2, col_dict$X1)
} else {
  palette_colors <- if (length(sample_names) > 12)
    colorRampPalette(brewer.pal(8, "Set2"))(length(sample_names)) else safe
  custom_colors <- setNames(palette_colors[seq_along(sample_names)], sort(sample_names))
}

HET     <- c("0/1", "1/0", "0|1", "1|0")
HOM_REF <- c("0/0", "0|0")
HOM_ALT <- c("1/1", "1|1")
NOCALL  <- c("./.", ".|.")
STATUS  <- c("het_pass", "het_lowdp", "het_lowgq", "het_filtered",
              "hom_ref", "hom_ref_filtered", "hom_alt", "alt_mismatch", "alt_nocall",
              "no_call", "missing")
CLASSES <- c("SNV", "INS", "DEL", "MNP")

## ---------------------------------------------------------------- statistics

## Exact two-sided binomial p-value against p = 0.5. The null is symmetric, so doubling
## the lower tail is exact rather than approximate, which lets it vectorise --
## binom.test() per row is unusable at a few million sites.
binom_p_half <- function(k, n) pmin(1, 2 * pbinom(pmin(k, n - k), n, 0.5))

## Williams (1982) moment estimator of beta-binomial overdispersion. Closed form, O(n).
## raw = TRUE returns the unfloored value: a negative estimate means the data are *under*
## dispersed relative to binomial, which is worth seeing rather than clamping to zero.
rho_moment <- function(k, n, raw = FALSE) {
  if (length(n) < 2) return(NA_real_)
  p   <- sum(k) / sum(n)
  x2  <- sum((k - n * p)^2 / (n * p * (1 - p)))
  est <- (x2 - length(n)) / sum(n - 1)
  if (raw) est else max(0, est)
}

## Maximum-likelihood rho on a capped random subsample. vglm routinely fails to converge
## when rho sits on the zero boundary; that is not an error and the moment estimate is the
## right fallback.
rho_betabinom <- function(k, n, max_n = 2e5) {
  if (!requireNamespace("VGAM", quietly = TRUE) || length(n) < 100)
    return(list(rho = NA_real_, method = "not attempted"))
  idx <- if (length(n) > max_n) sample.int(length(n), max_n) else seq_along(n)
  fit <- try(VGAM::vglm(cbind(k[idx], n[idx] - k[idx]) ~ 1, VGAM::betabinomial),
             silent = TRUE)
  if (inherits(fit, "try-error"))
    return(list(rho = NA_real_, method = "vglm failed (rho likely at boundary)"))
  list(rho = unname(VGAM::Coef(fit)[["rho"]]),
       method = sprintf("vglm (n = %d)", length(idx)))
}

## Exact two-sided beta-binomial p-value by the minimum-likelihood definition: the total
## probability of every outcome no more likely than the one observed. Correct for
## asymmetric mu too, where doubling a tail would not be. Cost is set by the number of
## distinct depths, not the number of sites.
bb_p_exact <- function(k, n, rho, mu = 0.5) {
  if (!is.finite(rho) || rho <= 0) return(binom_p_half(k, n))
  s <- (1 - rho) / rho
  a <- mu * s
  b <- (1 - mu) * s
  out <- numeric(length(k))
  for (nn in unique(n)) {
    dens <- VGAM::dbetabinom.ab(0:nn, nn, a, b)
    ## signif() before ranking so outcomes equal by symmetry but differing in the last
    ## ulp are treated as ties and get the same p-value
    key  <- signif(dens, 12)
    cs   <- cumsum(dens[order(key)])
    tab  <- cs[rank(key, ties.method = "max")]
    idx  <- which(n == nn)
    out[idx] <- pmin(1, tab[k[idx] + 1L])
  }
  out
}

## Distribution-free CI for the median from order statistics: exact, and one sort rather
## than a thousand resamples of a multi-million-row vector.
median_ci <- function(x, conf = 0.95) {
  x <- sort(x); n <- length(x); a <- (1 - conf) / 2
  c(median = median(x),
    lo = x[max(1, qbinom(a, n, 0.5))],
    hi = x[min(n, qbinom(1 - a, n, 0.5) + 1)])
}

## ---------------------------------------------------------------- the site universe

## The denominator is every true het, not every row the caller produced. That is the whole
## point of the dropout measurement, so the truth list is the left side of every join.
sites <- fread(sites_file, header = FALSE, col.names = c("chrom", "pos", "ref", "alt"))
sites[, `:=`(key4 = paste(chrom, pos, ref, alt, sep = "_"),
             key2 = paste0(chrom, "_", pos))]
sites[, var_class := fcase(
        nchar(ref) == 1L & nchar(alt) == 1L, "SNV",
        nchar(alt)  >  nchar(ref),           "INS",
        nchar(ref)  >  nchar(alt),           "DEL",
        default =                            "MNP")]

gc_tbl <- unique(fread(gc_file, header = FALSE, col.names = c("key2", "gc")), by = "key2")
sites  <- merge(sites, gc_tbl, by = "key2", all.x = TRUE)

classify <- function(query_path, depth_path, nm) {
  q <- fread(query_path, na.strings = c(".", "", "NA"),
             dec = ".", colClasses = c(ad = "character"))
  if (!any(grepl(",", q$ad, fixed = TRUE)))
    stop("no comma-separated AD values in ", query_path,
         " (class ", class(q$ad), ") -- fread mis-parsed the column")
  ad_ref <- as.integer(sub(",.*$",    "", q$ad))
  ad_alt <- as.integer(sub("^[^,]*,", "", q$ad))
  q[, `:=`(ad_ref = ad_ref,
           ad_alt = ad_alt,
           dp     = as.integer(dp),
           gq     = as.integer(gq),
           key4   = paste(chrom, pos, ref, alt, sep = "_"),
           key2   = paste0(chrom, "_", pos))]
  q <- unique(q, by = "key4")

  dep <- fread(depth_path, header = FALSE, col.names = c("chrom", "pos", "cov"))
  dep[, key2 := paste0(chrom, "_", pos)]

  d <- merge(sites, q[, .(key4, filter, gt, dp, ad_ref, ad_alt, gq)],
             by = "key4", all.x = TRUE)
  d <- merge(d, dep[, .(key2, cov)], by = "key2", all.x = TRUE)

  ## a truth site with no exact allele match, but where the caller did emit something,
  ## is a different failure from one it never mentioned
  d[, any_at_pos := key2 %in% q$key2]
  d[, ad_tot := fifelse(is.na(ad_ref) | is.na(ad_alt), NA_integer_, ad_ref + ad_alt)]

  ## first match wins, and NA is never TRUE, so every comparison on a possibly-missing
  ## field is guarded explicitly
  d[, status := fcase(
        is.na(gt) &  any_at_pos,                         "alt_mismatch",
        is.na(gt),                                       "missing",
        gt %in% NOCALL,                                  "no_call",
        gt %in% HOM_REF & !is.na(filter) & filter != "PASS", "hom_ref_filtered",
        gt %in% HOM_REF,                                 "hom_ref",
        gt %in% HOM_ALT,                                 "hom_alt",
        gt %in% HET & !is.na(filter) & filter != "PASS", "het_filtered",
        gt %in% HET & (is.na(ad_tot) | ad_tot < min_dp), "het_lowdp",
        gt %in% HET & (is.na(gq) | gq < min_gq),         "het_lowgq",
        gt %in% HET,                                     "het_pass",
        default =                                        "missing")]
  d[, `:=`(sample = nm, cov = fifelse(is.na(cov), 0L, as.integer(cov)))]
  d[]
}

cls <- rbindlist(lapply(seq_along(sample_names), function(i)
  classify(query_files[i], depth_files[i], sample_names[i])), use.names = TRUE)
cls[, `:=`(sample    = factor(sample,    levels = sample_names),
           status    = factor(status,    levels = STATUS),
           var_class = factor(var_class, levels = CLASSES))]

## ---------------------------------------------------------------- dropout table

overall <- cls[, .N, by = .(sample, var_class, status)][
  , frac := N / sum(N), by = .(sample, var_class)]

dropout_tbl <- cls[, .(n_truth_het = .N,
                       recovered   = sum(status == "het_pass"),
                       hom_ref     = sum(status == "hom_ref"),
                       hom_alt     = sum(status == "hom_alt"),
                       alt_mismatch = sum(status == "alt_mismatch"),
                       missing     = sum(status == "missing"),
                       median_cov  = as.numeric(median(cov))),
                   by = .(sample, var_class)]
dropout_tbl[, `:=`(recovery_rate = recovered / n_truth_het,
                   fn_rate       = 1 - recovered / n_truth_het,
                   dropout_ratio_strict = hom_ref / pmax(1, hom_alt),          # confident calls only
                   dropout_ratio_all    = (hom_ref + hom_ref_filtered) / pmax(1, hom_alt)]
write_tsv(dropout_tbl[order(sample, var_class)], out_dropout)

## ---------------------------------------------------------------- VAF and dispersion

dat <- cls[status == "het_pass"]
dat[, vaf := ad_alt / ad_tot]
setnames(dat, "ad_tot", "n")

rho_tbl <- dat[, {
  bb <- rho_betabinom(ad_alt, n)
  .(rho_moment     = rho_moment(ad_alt, n),
    rho_moment_raw = rho_moment(ad_alt, n, raw = TRUE),
    rho_ml         = bb$rho,
    rho_method     = bb$method)
}, by = .(sample, var_class)]
rho_tbl[, rho_used := fifelse(is.na(rho_ml), rho_moment, rho_ml)]

dat <- merge(dat, rho_tbl[, .(sample, var_class, rho_used)],
             by = c("sample", "var_class"))

## Both tests are kept. The binomial assumes independent reads at p = 0.5 and is
## anti-conservative whenever rho > 0; the gap between the two FDR rates is the size of
## that effect and is more informative than either number alone.
dat[, p_binom := binom_p_half(ad_alt, n)]
dat[, p_bb    := bb_p_exact(ad_alt, n, rho_used[1]), by = .(sample, var_class)]
dat[, `:=`(fdr_binom = p.adjust(p_binom, method = "BH"),
           fdr_bb    = p.adjust(p_bb,    method = "BH")), by = .(sample, var_class)]

vaf_tbl <- dat[, {
  ci <- median_ci(vaf)
  .(n_sites          = .N,
    median_depth     = as.numeric(median(n)),
    median_vaf       = ci[["median"]],
    vaf_ci_lo        = ci[["lo"]],
    vaf_ci_hi        = ci[["hi"]],
    frac_fdr05_binom = mean(fdr_binom < 0.05),
    frac_fdr05_bb    = mean(fdr_bb    < 0.05),
    n_fdr05_bb       = sum(fdr_bb < 0.05))
}, by = .(sample, var_class)]
vaf_tbl <- merge(vaf_tbl, rho_tbl, by = c("sample", "var_class"))
write_tsv(vaf_tbl[order(sample, var_class)], out_vaf)

fwrite(dat[fdr_bb < 0.05][order(fdr_bb),
           .(sample, var_class, chrom, pos, ref, alt,
             ad_ref, ad_alt, n, vaf, gc, p_bb, fdr_bb)],
       out_sig, sep = "\t")

## ---------------------------------------------------------------- plots

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
pct <- scales::percent

present <- CLASSES[CLASSES %in% unique(as.character(cls$var_class))]
core    <- intersect(c("SNV", "INS", "DEL"), present)

p_stack <- ggplot(overall, aes(sample, frac, fill = status)) +
  geom_col(width = 0.65) +
  facet_wrap(~ var_class, nrow = 1) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  scale_y_continuous(labels = pct) +
  labs(title = "Fate of every true heterozygous site",
       x = NULL, y = "fraction of truth-het sites", fill = NULL) +
  base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

## recovery against coverage is what separates "this polymerase loses hets" from
## "this sample was sequenced less deeply"
cls[, cov_bin := cut(cov, breaks = c(-1, 0, 5, 10, 15, 20, 30, 40, 60, Inf),
                     labels = c("0", "1-5", "6-10", "11-15", "16-20",
                                "21-30", "31-40", "41-60", ">60"))]

p_cov <- cls[, .(rate = mean(status == "het_pass"), n = .N),
             by = .(sample, var_class, cov_bin)][n >= 100] %>%
  ggplot(aes(cov_bin, rate, colour = sample, group = sample)) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  scale_y_continuous(labels = pct) +
  labs(title = "Het recovery against coverage",
       subtitle = "coverage from samtools depth, independent of the caller",
       x = "reads at site", y = "recovered as confident het", colour = NULL) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_gcdrop <- cls[!is.na(gc)][
  , gc_bin := cut(gc, breaks = seq(0, 1, by = 0.05))][
  , .(rate = mean(status == "het_pass"),
      drop = mean(status %in% c("hom_ref", "hom_alt")), n = .N),
    by = .(sample, var_class, gc_bin)][n >= 100] %>%
  melt(id.vars = c("sample", "var_class", "gc_bin", "n"),
       measure.vars = c("rate", "drop"), variable.name = "metric") %>%
  mutate(metric = recode(metric, rate = "recovered as het",
                                 drop = "called homozygous")) %>%
  ggplot(aes(gc_bin, value, colour = sample, group = sample)) +
  geom_line() + geom_point(size = 0.8) +
  facet_grid(metric ~ var_class, scales = "free_y") +
  scale_colour_manual(values = custom_colors) +
  scale_y_continuous(labels = pct) +
  labs(title = "Het recovery and dropout by local GC content",
       x = "GC fraction of the flanking window", y = NULL, colour = NULL) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_asym <- ggplot(dropout_tbl, aes(factor(sample, levels = sample_names), dropout_ratio_all,
                                  fill = factor(sample, levels = sample_names))) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Dropout asymmetry",
       subtitle = "hom_ref / hom_alt; 1 = symmetric, >1 = the alt allele is lost more often",
       x = NULL, y = "ratio", fill = NULL) +
  base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_density <- ggplot(dat, aes(vaf, colour = sample)) +
  geom_density(linewidth = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "VAF at recovered heterozygous sites",
       subtitle = "dashed line = unbiased expectation (0.5)",
       x = "alt / (ref + alt)", y = "density", colour = NULL) +
  base_theme

## Fixed-depth fit check. Conditioning on a single total depth is what makes this
## readable: pooled over depths the binomial becomes a mixture and looks wide for the
## wrong reason.
n_fit <- as.integer(names(which.max(table(dat$n))))
obs <- dat[n == n_fit & var_class %in% core,
           .(count = .N), by = .(sample, var_class, k = ad_alt)]
obs[, prop := count / sum(count), by = .(sample, var_class)]

curves <- rho_tbl[var_class %in% core, {
  bb <- if (is.finite(rho_used) && rho_used > 0) {
    s <- (1 - rho_used) / rho_used
    VGAM::dbetabinom.ab(0:n_fit, n_fit, s / 2, s / 2)
  } else dbinom(0:n_fit, n_fit, 0.5)   # rho == 0 collapses onto the binomial
  .(k = 0:n_fit, Binomial = dbinom(0:n_fit, n_fit, 0.5), `Beta-binomial` = bb)
}, by = .(sample, var_class)]
curves <- melt(curves, id.vars = c("sample", "var_class", "k"),
               variable.name = "model", value.name = "density")

p_fit <- ggplot(obs, aes(k, prop)) +
  geom_col(fill = "#8FA9D8", width = 1) +
  geom_vline(xintercept = n_fit / 2, colour = "grey30", linewidth = 0.3) +
  geom_line(data = curves, aes(k, density, colour = model, linetype = model),
            linewidth = 0.6) +
  facet_grid(var_class ~ sample) +
  scale_colour_manual(values = c(Binomial = "grey20", `Beta-binomial` = "#CC6677")) +
  scale_linetype_manual(values = c(Binomial = "dashed", `Beta-binomial` = "solid")) +
  labs(title = sprintf("Allele counts at fixed depth (ref + alt = %d)", n_fit),
       subtitle = "bars observed; the gap between the curves is the overdispersion",
       x = "alt reads", y = "density", colour = NULL, linetype = NULL) +
  base_theme

## Calibration: a correct null is uniform on (0,1), so its ECDF is the diagonal, and one
## that is too thin bows above it. Same statement as p_fit, in the units the FDR consumes.
p_cal <- melt(dat[, .(sample, var_class, Binomial = p_binom, `Beta-binomial` = p_bb)],
              id.vars = c("sample", "var_class"),
              variable.name = "test", value.name = "p") %>%
  ggplot(aes(p, colour = test)) +
  stat_ecdf(linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  facet_grid(var_class ~ sample) +
  scale_colour_manual(values = c(Binomial = "grey20", `Beta-binomial` = "#CC6677")) +
  labs(title = "p-value calibration",
       subtitle = "diagonal = correctly calibrated; bowing above it = anti-conservative",
       x = "p", y = "ECDF", colour = NULL) +
  base_theme

p_gcvaf <- dat[!is.na(gc)][
  , gc_bin := cut(gc, breaks = seq(0, 1, by = 0.05))][
  , .(med = median(vaf), lo = quantile(vaf, 0.25),
      hi = quantile(vaf, 0.75), n = .N), by = .(sample, var_class, gc_bin)][n >= 100] %>%
  ggplot(aes(gc_bin, med, colour = sample, group = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "VAF by local GC content",
       subtitle = "median with interquartile ribbon; bins with <100 sites dropped",
       x = "GC fraction of the flanking window", y = "VAF",
       colour = NULL, fill = NULL) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(out_plots, width = 11, height = 6)
print(p_stack); print(p_cov); print(p_gcdrop); print(p_asym)
print(p_density); print(p_fit); print(p_cal); print(p_gcvaf)
dev.off()