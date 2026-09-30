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
min_gq_hom   <- as.integer(snakemake@params[["min_gq_hom"]])
vaf_min_dp   <- as.integer(snakemake@params[["vaf_min_dp"]])
out_vaf      <- snakemake@output[["vaf"]]
out_dropout  <- snakemake@output[["dropout"]]
out_sig      <- snakemake@output[["sig"]]
out_plots    <- snakemake@output[["plots"]]
out_assign <- snakemake@output[["assign"]]
out_shift <- snakemake@output[["shift"]]


## New version Sunday 27 Sep 2026:
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
             "hom_ref", "hom_ref_filtered", "hom_alt",
             "alt_outcompeted", "alt_mismatch", "alt_nocall", "no_call", "missing")
CLASSES <- c("SNV", "INS", "DEL", "MNP")
HOMREF_OK <- c("PASS", "RefCall", ".")

## ---------------------------------------------------------------- statistics

## Exact two-sided binomial p-value against p = 0.5
binom_p_half <- function(k, n) pmin(1, 2 * pbinom(pmin(k, n - k), n, 0.5))

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
  ## AD arrives as "ref,alt" since bcftools norm -m -any leaves every record biallelic..
  q <- fread(query_path, na.strings = c(".", "", "NA"),
             dec = ".", colClasses = c(ad = "character", vaf_caller = "character"))
  ## fread can silently parse "," as a decimal sep, so stop when this happens
  if (!any(grepl(",", q$ad, fixed = TRUE)))
    stop("no comma-separated AD values in ", query_path,
         " (class ", class(q$ad), ") -- fread mis-parsed the column")
  ## Ref is second record
  ad_ref <- as.integer(sub(",.*$",    "", q$ad))
  
  ## Alt is first record
  ad_alt <- as.integer(sub("^[^,]*,", "", q$ad))
  
  ## Dp is caller-depth at that position
  q[, `:=`(ad_ref = ad_ref,
           ad_alt = ad_alt,
           dp     = as.integer(dp),
           gq     = as.integer(gq),
           vaf_caller = as.numeric(vaf_caller),
           key4   = paste(chrom, pos, ref, alt, sep = "_"),
           key2   = paste0(chrom, "_", pos))]
  q <- unique(q, by = "key4")
  
  ## Samtools depth - might differ from caller depth and ad_ref + ad_alt because of filtering
  dep <- fread(depth_path, header = FALSE, col.names = c("chrom", "pos", "cov"))
  dep[, key2 := paste0(chrom, "_", pos)]
  
  q[, matched := TRUE]
  d <- merge(sites, q[, .(key4, matched, filter, gt, dp, ad_ref, ad_alt, gq)],
             by = "key4", all.x = TRUE)
  d <- merge(d, dep[, .(key2, cov)], by = "key2", all.x = TRUE)
  
  ## What the caller actually genotyped at this position, carried onto the truth row.
  win <- q[gt %in% c(HET, HOM_ALT)][order(key2, -ad_alt)][, .(
    called_alt    = alt[1],
    called_gt     = gt[1],
    called_ad_alt = ad_alt[1],
    called_gq     = gq[1],
    n_alt_at_pos  = .N), by = key2]
  d <- merge(d, win, by = "key2", all.x = TRUE)
  
  ## a truth site with no exact allele match, but where the caller did emit something,
  ## can be ref_hom (mostly because of multiallelic records), or a mismatch:
  ## "saw something and declined to call it"  vs  "called a different allele":
  qpos <- q[order(key2, -ad_alt)][, .SD[1L], by = key2]
  d <- merge(d, qpos[, .(key2, gt_pos = gt)], by = "key2", all.x = TRUE)
  d[, any_at_pos  := key2 %in% q$key2]
  d[, other_called := !is.na(called_alt)]
  d[, ad_tot := fifelse(is.na(ad_ref) | is.na(ad_alt), NA_integer_, ad_ref + ad_alt)]
  
  ## fraction of the reads at the locus that the caller assigned to either allele. A site
  ## with DP 31 and AD 12,5 assigned 55% -- the rest support neither allele cleanly, which
  ## is what alignment ambiguity at an indel looks like from the VCF alone.
  d[, assigned_frac := fifelse(is.na(dp) | dp == 0L | is.na(ad_tot),
                               NA_real_, ad_tot / dp)]
  ## negative means the caller's allele was shorter than the truth's.
  d[, called_len_diff := fifelse(is.na(called_alt), NA_integer_,
                                 nchar(called_alt) - nchar(alt))]
  
  ## first match wins, and NA is never TRUE, so every comparison on a possibly-missing
  ## field is guarded explicitly
  d[, status := fcase(
    is.na(matched) & gt_pos %in% NOCALL,             "alt_nocall",
    is.na(matched) & any_at_pos,                     "alt_mismatch",
    is.na(matched),                                  "missing",
    is.na(gt) | gt %in% NOCALL,                      "no_call",
    ## After norm -m -any a 0/0 on a split record means "not this allele", but looks like hom_ref
    gt %in% HOM_REF & other_called &
      !is.na(ad_alt) & ad_alt > 0L,                "alt_outcompeted",
    gt %in% HOM_REF & other_called,                  "alt_mismatch",
    ## a RefCall, or a 0/0 the caller had no confidence in, is an abstention rather
    ## than a confident homozygous-reference call. Both count as not-recovered, but
    ## only the confident one belongs in dropout_ratio.
    gt %in% HOM_REF & (is.na(filter) | filter != "PASS" | !(filter %in% HOMREF_OK) |
                         is.na(gq) | gq < min_gq_hom), "hom_ref_filtered",
    gt %in% HOM_REF,                                 "hom_ref",
    gt %in% HOM_ALT,                                 "hom_alt",
    gt %in% HET & !is.na(filter) & filter != "PASS", "het_filtered",
    gt %in% HET & (is.na(ad_tot) | ad_tot < min_dp), "het_lowdp",
    gt %in% HET & (is.na(gq) | gq < min_gq),         "het_lowgq",
    gt %in% HET,                                     "het_pass",
    default =                                        "missing")]
  d[, lost_despite_support := status == "alt_outcompeted" &
      !is.na(called_ad_alt) & ad_alt > called_ad_alt]
  d[, `:=`(sample = nm, cov = as.integer(cov))]
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
                       hom_ref_filtered = sum(status == "hom_ref_filtered"),
                       hom_alt     = sum(status == "hom_alt"),
                       alt_outcompeted = sum(status == "alt_outcompeted"),
                       frac_lost_despite_support =
                         mean(lost_despite_support[status == "alt_outcompeted"]),
                       median_assigned = median(assigned_frac, na.rm = TRUE),
                       alt_mismatch = sum(status == "alt_mismatch"),
                       alt_nocall  = sum(status == "alt_nocall"),
                       missing     = sum(status == "missing"),
                       median_cov  = as.numeric(median(cov, na.rm = TRUE))),
                   by = .(sample, var_class)]
dropout_tbl[, `:=`(
  recovery_rate = recovered / n_truth_het,
  fn_rate       = 1 - recovered / n_truth_het)]

write_tsv(dropout_tbl[order(sample, var_class)], out_dropout)

## ---------------------------------------------------------------- VAF and dispersion

## Table for het_pass variants (those that would have been included in truth benchmarks)
dat <- cls[status == "het_pass" & ad_tot >= vaf_min_dp]
dat[, vaf := ad_alt / ad_tot]

setnames(dat, "ad_tot", "n")

## Binomial test for each site, and a FDR test for each site as well
dat[, p_binom   := binom_p_half(ad_alt, n)]
dat[, fdr_binom := p.adjust(p_binom, method = "BH"), by = .(sample, var_class)]

## Median depth as reported by the caller here, not the samtools depth cov
vaf_tbl <- dat[,
  .(n_sites          = .N,
    median_depth     = as.numeric(median(n)),
    median_assigned  = median(assigned_frac, na.rm = TRUE),
    median_vaf_caller = median(vaf_caller, na.rm = TRUE),
    median_vaf = median(vaf),
    frac_fdr05       = mean(fdr_binom < 0.05),
    n_fdr05          = sum(fdr_binom < 0.05)), by = .(sample, var_class)]

write_tsv(vaf_tbl[order(sample, var_class)], out_vaf)

## Variants passing here likely CNVs
fwrite(dat[fdr_binom < 0.05][order(fdr_binom),
                             .(sample, var_class, chrom, pos, ref, alt,
                               ad_ref, ad_alt, n, vaf, gc, p_binom, fdr_binom)],
       out_sig, sep = "\t")


## ---------------------------------------------------------------- plots

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
pct <- scales::percent

present <- CLASSES[CLASSES %in% unique(as.character(cls$var_class))]
core    <- intersect(c("SNV", "INS", "DEL"), present)

## Get level order and sample order right
potential_levels <- c("het_pass","hom_alt","hom_ref","alt_outcompeted",
                               "alt_mismatch","missing","hom_ref_filtered","het_filtered",
                               "het_lowdp","het_lowgq", "no_call")

overall$status <- factor(overall$status, levels = rev(potential_levels[potential_levels %in% unique(overall$status)]), 
                         ordered = TRUE)

overall_sort_order <- overall %>%
  group_by(sample) %>%
  mutate(
    "correct" = sum(N[status== "het_pass"])
  ) %>%
  arrange(correct) %>% 
  filter(status == "het_pass" & var_class == "SNV") %>%
  mutate(sample = factor(sample))
  
overall_mut <- overall%>% 
  mutate(sample = factor(sample, levels = overall_sort_order$sample, ordered = TRUE))

## From Spectral: Colors for stack plot
stack_colors <- c("alt_nocall" = "#9E0142",
                  "alt_mismatch" = "#F46D43",
                  "missing" = "#D53E4F",
                  "no_call" = "#FDAE61",
                  "hom_ref_filtered" = "#FEE08B",
                  "hom_ref" =  "#FFFFBF",
                  "hom_alt" = "#E6F598",
                  "het_filtered" = "#ABDDA4",
                  "het_lowdp" =  "#66C2A5",
                  "het_lowgq" = "#5E4FA2",
                  "het_pass" = "#3288BD")

p_stack <- ggplot(overall_mut, aes(frac, sample , fill = status)) +
  geom_col(width = 0.65) +
  facet_wrap(~ var_class, nrow = 1) +
  scale_fill_manual(values = stack_colors) +
  scale_x_continuous(labels = pct) +
  labs(title = "Fate of every true heterozygous site",
       x = NULL, y = "fraction of truth-het sites", fill = NULL) +
  base_theme

## recovery against coverage is what separates "this polymerase loses hets" from
## "this sample was sequenced less deeply"
## A confident het call at zero coverage is impossible, so if this fires either the
## depth track is incomplete or samtools depth is filtering more strictly than the caller
#stopifnot(cls[!is.na(cov) & cov == 0L & status == "het_pass", .N] == 0L)

## How many of mapped reads were assigned by caller
reads_assigned <- cls[!is.na(dp) & dp > 0, .(med_assigned = median(ad_tot / cov, na.rm = TRUE), n = .N),
    by = .(var_class, status)][order(var_class, status)]

write_tsv(reads_assigned, out_assign)

## Read-depth for correct calls
cap <- 2L * as.integer(median(cls$cov, na.rm = TRUE))
cls[, cov_plot := pmin(cov, cap + 1L)]     # cap + 1 is the ">cap" pile

p_cov <- cls[!is.na(cov_plot), .(rate = mean(status == "het_pass"), n = .N),
             by = .(sample, var_class, cov_plot)][n >= 20] %>%
  ggplot(aes(cov_plot, rate, colour = sample)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = cap + 0.5, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_x_continuous(breaks = scales::breaks_width(10)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = custom_colors) +
  labs(title = "Het recovery against coverage",
       subtitle = sprintf("one point per read depth; everything above %d pooled at %d", cap, cap + 1L),
       x = "reads at site", y = "recovered as confident het", colour = NULL) +
  base_theme

## For labelling the GC axis
#cats <- length(unique(cut(cls$gc, breaks = seq(0, 1, by = 0.05),
#    include.lowest = TRUE)))
#label_cats <- cats * 5

gc_breaks <- seq(0, 1, by = 0.05)
gc_labels <- as.character(head(gc_breaks, -1) * 100)   # 20 labels, unconditionally

p_gcdrop <- cls[!is.na(gc)][
  , gc_bin := cut(gc, breaks = gc_breaks, 
                  labels = gc_labels,
                  include.lowest = TRUE)][
    , .(rate = mean(status == "het_pass"),
        drop = mean(status %in% c("hom_alt")),
        ref_hom = mean(status %in% c("missing", "hom_ref")),
        miss = mean(status %in% c("alt_outcompeted", "alt_mismatch")), n = .N),
    by = .(sample, var_class, gc_bin)][n >= 50] %>%
  melt(id.vars = c("sample", "var_class", "gc_bin", "n"),
       measure.vars = c("rate", "drop", "ref_hom", "miss"), variable.name = "metric") %>%
  mutate(metric = recode(metric, rate = "recovered as het",
                         drop = "alt homozygous", 
                         ref_hom = "ref homozygous",
                         miss = "allele call mismatch")) %>%
  ggplot(aes(gc_bin, value, colour = sample, group = sample)) +
  geom_line() + 
  geom_point(size = 0.8) +
  facet_grid(metric ~ var_class, scales = "free_y") +
  scale_colour_manual(values = custom_colors) +
  scale_y_continuous(labels = pct) +
  labs(title = "Het recovery and dropout by local GC content",
       x = "GC fraction of the flanking window", y = NULL, colour = NULL) +
  base_theme

## VAF density plots
p_density <- ggplot(cls, aes(ad_alt/(ad_alt + ad_ref), colour = sample)) +
  geom_density(linewidth = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "VAF at recovered heterozygous sites",
       subtitle = "dashed line = unbiased expectation (0.5)",
       x = "alt / (ref + alt)", y = "density", colour = NULL) +
  base_theme

## The Caller's density
p_density_caller <- ggplot(cls, aes(vaf_caller, colour = sample)) +
  geom_density(linewidth = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "VAF at recovered heterozygous sites",
       subtitle = "dashed line = unbiased expectation (0.5) if all reads were assigned",
       x = "alt / depth", y = "density", colour = NULL) +
  base_theme


## Fixed-depth fit per sample, adding reads at ref+alt = fit assigned to alt, as well as
## reads at that site according to samtools depth
depth_tbl <- dat[var_class %in% core,
                 .(depth_fit = as.integer(names(which.max(table(n))))), by = sample]
setorder(depth_tbl, sample)     # sample is a factor, so this is sample_names order
lab <- sprintf("%s  (ref + alt = %d)", as.character(depth_tbl$sample), depth_tbl$depth_fit)
depth_tbl[, panel := factor(lab, levels = lab)]

dat[depth_tbl, depth_fit := i.depth_fit, on = "sample"]

obs <- dat[n == depth_fit & var_class %in% core,
           .(count = .N), by = .(sample, depth_fit, var_class, k = ad_alt, ar = cov - ad_alt)]
obs[, prop := count / sum(count), by = .(sample, var_class)]
obs <- merge(obs, depth_tbl[, .(sample, panel)], by = "sample")

## Long Table for plotting
obs_long <- obs %>%
  pivot_longer(cols = c(ar, k), names_to = "alt_read_class", values_to = "read_count") %>%
  group_by(var_class, sample, alt_read_class, read_count)%>%
  mutate(overall_prop = sum(prop)) %>%
  select(sample, var_class, alt_read_class, read_count, overall_prop) %>%
  distinct(.keep_all = TRUE)

## Read mapping mismatch for alt alleles  
p_fit <- ggplot(obs_long, aes(read_count, overall_prop)) +
  geom_col(aes(fill = sample, alpha = alt_read_class),
           position = "identity", width = 1, colour = NA) +
  geom_step(aes(linetype = alt_read_class), direction = "mid",
            colour = "grey20", linewidth = 0.4) +
  geom_vline(data = depth_tbl, aes(xintercept = depth_fit / 2),
             colour = "grey20", linewidth = 0.3) +
  facet_grid(var_class ~ sample) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_alpha_manual(values = c(ar = 0.35, k = 1),
                     labels = c(ar = "cov - alt", k = "alt")) +
  scale_linetype_manual(values = c(ar = "dotted", k = "solid"),
                        labels = c(ar = "cov - alt", k = "alt")) +
  labs(x = "reads (at alt + ref = median(depth))", y = "density", alpha = NULL, linetype = NULL) +
  theme_bw() + 
  coord_cartesian(xlim = c(0, 1.5 * max(depth_tbl$depth_fit))) +
  theme(legend.position = "bottom")


## calculate the ratio of alt to all mapped reads (according to samtools)
dat[, depth := n]
depth_grid <- dat[var_class %in% core, .N, by = .(sample, var_class, depth)][N >= 50]
setorder(depth_grid, sample, var_class, -N)
depth_grid <- depth_grid[, head(.SD, 2L * max(depth_tbl$depth_fit)), by = .(sample, var_class)]

obs_all <- dat[depth_grid[, .(sample, var_class, depth)],
               on = .(sample, var_class, depth)][
                 , .(count = .N), by = .(sample, var_class, depth, k = ad_alt, ar = cov - ad_alt)]

shift_read_table <- obs_all %>%
  group_by(sample, var_class) %>%
  summarise(count_k = sum(count*k),
            count_ar = sum(count*ar)) %>%
  mutate(alt_all_ratio = count_k/count_ar)

write_tsv(shift_read_table, out_shift)

## For all obersvations with more than 50 Vars
obs_all_long <- obs_all %>%
  pivot_longer(c(ar, k), names_to = "alt_read_class", values_to = "read_count") %>%
  group_by(sample, var_class, alt_read_class, read_count) %>%
  summarise(n_sites = sum(count), .groups = "drop") %>%
  group_by(sample, var_class, alt_read_class) %>%
  mutate(main_prop = n_sites / sum(n_sites)) %>%
  ungroup()

## plot obs_all
p_fit_all <- ggplot(obs_all_long, aes(read_count, main_prop)) +
  geom_col(aes(fill = sample, alpha = alt_read_class),
           position = "identity", width = 1, colour = NA) +
  geom_step(aes(linetype = alt_read_class), direction = "mid",
            colour = "grey20", linewidth = 0.4) +
  geom_vline(data = depth_tbl, aes(xintercept = depth_fit / 2),
             colour = "grey20", linewidth = 0.3) +
  facet_grid(var_class ~ sample)+
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_alpha_manual(values = c(ar = 0.35, k = 1),
                     labels = c(ar = "cov - alt", k = "alt")) +
  scale_linetype_manual(values = c(ar = "dotted", k = "solid"),
                        labels = c(ar = "cov - alt", k = "alt")) +
  labs(x = "reads (at alt + ref = median(depth))", y = "density", alpha = NULL, linetype = NULL) +
  coord_cartesian(xlim = c(0, 1.5 * max(depth_tbl$depth_fit))) +
  #xlim(c(0, depth_tbl$depth_fit*2)) + 
  theme_bw() + theme(legend.position = "bottom")

## Calibration: a correct null is uniform on (0,1), so its ECDF is the diagonal, and one
p_gcvaf <- dat[!is.na(gc)][
  , gc_bin := cut(gc, breaks = gc_breaks, 
                  labels = gc_labels,
                  include.lowest = TRUE)][
    , .(med = median(vaf), lo = quantile(vaf, 0.25),
        hi = quantile(vaf, 0.75), n = .N), by = .(sample, var_class, gc_bin)][n >= 100] %>%
  ggplot(aes(gc_bin, med, colour = sample, group = sample, fill = sample)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = NA) +
  geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ var_class, nrow = 1) +
  scale_colour_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "VAF by local GC content",
       subtitle = "median with interquartile ribbon; bins with <100 sites dropped",
       x = "GC % of the flanking window (100 bp)", y = "VAF",
       colour = NULL, fill = NULL) +
  base_theme

pdf(out_plots, width = 11, height = 6)
print(p_stack); print(p_cov); print(p_gcdrop);
print(p_density); print(p_density_caller), print(p_fit); print(p_fit_all); print(p_gcvaf)
dev.off()