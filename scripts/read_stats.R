library(tidyverse)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(scales)
library(patchwork)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input for plotting/sumstats
input_path_all = snakemake@params[["stats_path"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output (N(x) + compleasm tables, final plot)
gc_hist_hist_out = snakemake@output[["gc_hist_hist"]]
qv_plot_out = snakemake@output[["qv_plot"]]
standard_hist_plot_readeval_out = snakemake@output[["standard_hist_plot_readeval"]]
rdeval_like_out = snakemake@output[["rdeval_like"]]


## set colors for names
input_names <- c(input_names)
input_names <- unlist(strsplit(input_names, split = ","))
if (!is.null(input_names)){
  labels <- sort(input_names)
  if (!is.null(in_colors)) {
  col_dict <- read_delim(in_colors, col_names = FALSE)
  custom_colors <- setNames(col_dict$X2, col_dict$X1)
  } else {
  palette_colors <- safe
  custom_colors <- setNames(palette_colors, labels)
}
}

## find breakpoints
get_breaks <- function(min, max, steps = c(5,10)) {
  logstart <- floor(log10(min))
  logend <- ceiling(log10(max))
  
  breaks <- c()
  for (e in logstart:logend) {
  breaks <- c(breaks, steps*10^e)
  }
  
  return(breaks)
}


## Function to read and manipulate a rdeval dumped sequence report tsv file
read_rdeval <- function(path, name){
  
  rdeval_in <- fread(path, header = T)
  
  rdeval_in_named  <- rdeval_in %>%
    mutate(polymerase = name)
  
  gc_bins_hist <- rdeval_in_named %>%
    mutate(GC_bin = round(GC,1)) %>%
    group_by(Length, GC_bin) %>%
    summarize(num_reads = n()) %>%
    ungroup() %>%
    group_by(GC_bin) %>%
    mutate(relative_contribution = num_reads/sum(num_reads),
           mean_rl = sum(Length*num_reads)/sum(num_reads),
           median_rl = median(rep(Length,num_reads)),
           polymerase = name) %>%
    filter(sum(num_reads) >= 100)
  
  gc_bins_dens <- rdeval_in_named %>%
    mutate(GC_bin = round(GC,1)) %>%
    group_by(GC_bin) %>%
    mutate(mean_rl = mean(Length),
           median_rl = median(Length),
           polymerase = name) %>%
    filter(n() >= 100)
  
  m <- rdeval_in_named %>%
  count(Length) %>%
  summarise(max_value = max(n, na.rm = TRUE))
  
  scatter_breaks <- get_breaks(1, m$max_value)
  
  scatter <- ggplot(rdeval_in_named, aes(Length, `Average Quality`)) +
  geom_bin_2d(bins=100, binwidth = c(1, 1)) +
  theme_bw() +
  scale_fill_gradient(
    trans = "log10",
    low = "lightblue",
    high = "purple",
    name = "Count (log10)",
    breaks = scatter_breaks,
    label = comma) + 
  theme(legend.position = "left") +
  xlim(c(-1, mean(rdeval_in_named$Length) + 4*sd(rdeval_in_named$Length))) 
  
  print(max(rdeval_in_named$`Average Quality`))
  
  read_length_histogram <- ggplot(rdeval_in_named, aes(Length)) +
    geom_histogram(binwidth = 1, aes(fill = polymerase)) +
    scale_fill_manual(values = custom_colors) + 
    theme_void() +
    xlim(c(-1, mean(rdeval_in_named$Length) + 4*sd(rdeval_in_named$Length)))
    
  
  qc_histogram <- ggplot(rdeval_in_named, aes(`Average Quality`)) +
    geom_histogram(bins = 50, aes(fill = polymerase)) +
    scale_fill_manual(values = custom_colors) +
    theme_void() +
    coord_flip()
    #xlim(c(0,50)) 

  
  combined <- read_length_histogram + plot_spacer() + scatter + qc_histogram + 
    plot_layout(
      ncol = 2, 
      nrow = 2, 
      widths = c(4, 1),
      heights = c(1, 4),
      guides = "collect"
    )
  
  out_list <- list()
  out_list[[1]] <- rdeval_in_named
  out_list[[2]] <- gc_bins_hist
  out_list[[3]] <- gc_bins_dens
  out_list[[4]] <- combined
    
  return(out_list)
}

## create rdeval plots
rdeval_plots <- function(path) {
  
  ## set up variables
  tmp_path <- paste(path, "*rdeval_dump.tsv", sep="/")
  in_files <- Sys.glob(tmp_path)

  rdeval_plots <- list()
  rdeval_gc_bins <- list()
  rdeval_gc_dens <- list()
  rdeval_standard <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    if (!is.na(pattern)) {
    tmp_frame <- read_rdeval(in_files[i], name=pattern)

    rdeval_plots[[i]] <- tmp_frame[[4]]
    rdeval_gc_bins[[i]] <- tmp_frame[[2]]
    rdeval_gc_dens[[i]] <- tmp_frame[[3]]
    rdeval_standard[[i]] <- tmp_frame[[1]]
  }
  }
  
  gc_bin_frame <- rbindlist(rdeval_gc_bins)
  gc_dens_frame <- rbindlist(rdeval_gc_dens)
  standard_frame <- rbindlist(rdeval_standard)
  
  ## plot
  
  dens_plot <- ggplot(gc_dens_frame, aes(Length)) +
    geom_density(aes(fill = polymerase), alpha = 0.6) + 
    geom_vline(aes(xintercept = mean_rl), colour = 'black', size = 0.75) +
    geom_vline(aes(xintercept = median_rl), colour = 'grey', size = 0.75, linetype = "dashed") +
    facet_grid(GC_bin ~ polymerase) +
    scale_fill_manual(values = custom_colors) +
    theme_bw()
    
  standard_frame <- standard_frame %>%
  group_by(polymerase) %>%
  mutate(mean_length = mean(Length),
            median_length = median(Length))

  standard_hist_plot <- ggplot(standard_frame, aes(Length)) +
    geom_histogram(binwidth = 1, aes(fill = polymerase)) +
    scale_fill_manual(values = custom_colors) +
    geom_vline(aes(xintercept = mean_length), colour = 'black', size = 0.75) +
    geom_vline(aes(xintercept = median_length), colour = 'grey', size = 0.75, linetype = "dashed") +
    facet_grid(rows = vars(polymerase)) +
    theme_bw() +
    xlim(c(-1, mean(standard_frame$Length) + 4*sd(standard_frame$Length)))
    
  #rdeval_like <- patchwork::wrap_plots(rdeval_plots, 
  #                                nrow = 2, ncol = 5)
  
  out_list <- list()
  out_list[[1]] <- dens_plot
  out_list[[2]] <- standard_hist_plot
  out_list[[3]] <- rdeval_plots
  
  return(out_list)
}

## create qv plots
qv_plots <- function(path){
  
  read_qv <- function(file_path, name) {
    qv_file <- read_delim(file_path, col_names = T)
    qv_file <- qv_file %>%
      filter(fraction1 > 0) %>%
      mutate(polymerase = name)
  }
  
  ## set up variables
  tmp_path <- paste(path, "*qchist.txt", sep="/")
  in_files <- Sys.glob(tmp_path)
  qv_list <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    if (!is.na(pattern)) {
    tmp_frame <- read_qv(in_files[i], name=pattern)
    qv_list[[i]] <- tmp_frame
  }
  }
  
  qv_frame <- rbindlist(qv_list)
  
  ## plot
  
  qv_plot <- ggplot(qv_frame, aes(polymerase, fraction1)) +
    geom_bar(aes(fill = factor(`#Quality`)), position="stack", stat="identity") +
    scale_fill_viridis_d(name = "QV",direction = -1) +
    theme_bw() +
    ylab("percentage of reads")
  
  ## output
  qv_list <- list()
  qv_list[[1]] <- qv_plot
  qv_list[[2]] <- qv_frame
  
  return(qv_list)

}

## create output plots

## qv plot
qv_out <- qv_plots(input_path_all)
qv_plot <- qv_out[[1]]

## rdeval
rdeval_out <- rdeval_plots(input_path_all)

## Standard hist (rdeval)
standard_hist_plot_readeval <- rdeval_out[[2]]

## GC hist (rdeval)
gc_hist_hist <- rdeval_out[[1]]

## Rdeval-like
rdeval_like <- rdeval_out[[3]]
  
## write output
ggsave(qv_plot_out, plot = qv_plot, device = 'pdf', width = 12, height = 14, units = "in")
ggsave(standard_hist_plot_readeval_out, plot = standard_hist_plot_readeval, device = 'pdf', width = 12, height = 14, units = "in")
ggsave(gc_hist_hist_out, plot = gc_hist_hist, device = 'pdf', width = 12, height = 14, units = "in")


# Open PDF device
pdf(rdeval_like_out, width = 7, height = 7)

# Print each plot on a new page
for (p in rdeval_like) {
  print(p)
}

# Close the device to write file
dev.off()
