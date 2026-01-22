library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)
library(RColorBrewer)
library(scales)
library(R.utils)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input
input_path_pandepth = snakemake@params[["pandepth_path"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output 
output_boxplot = snakemake@output[["boxplot"]]
output_line_mean = snakemake@output[["mean_line"]]
output_line_median = snakemake@output[["median_line"]]
output_line_iqr = snakemake@output[["iqr_line"]]
output_bar_lci = snakemake@output[["lci_bar"]]

## set colors for names
## set colors for names
input_names <- c(input_names)
input_names <- unlist(strsplit(input_names, split = ","))
if (!is.null(input_names)){
  labels <- sort(input_names)
  if (in_colors != "") {
  col_dict <- read_delim(in_colors, col_names = FALSE)
  custom_colors <- setNames(col_dict$X2, col_dict$X1)
  } else {
  palette_colors <- safe
  custom_colors <- setNames(palette_colors, labels)
}
}

## Function to read pandepth infiles
## Including calculation of LCI/ genome breadth
read_pandepth <- function(pandepth_file, asm = "") {
  
  read_last <- function(file) {
    tail(readLines(file), n=1)
  }
  
  mead_depth <- sub(".*MeanDepth: *([0-9.]+).*", "\\1", read_last(pandepth_file))
  
  in_pandepth <- fread(pandepth_file, header = T)
  out_frame <- in_pandepth %>%
    mutate(asm = asm,
           NormDepth = MeanDepth / as.numeric(mead_depth),
           roundGC = round((as.numeric(`GC(%)`)/100), 2))

  out_frame_length <- sum(as.numeric(out_frame$Length), na.rm = T)
  
  out_frame_lci <- out_frame %>%
    filter(NormDepth < 0.5) %>%
    summarise(LCI_length = sum(as.numeric(Length), na.rm = T))
  
  out_genome_breadth <- out_frame %>%
    summarise(genome_breadth = sum(as.numeric(CoveredSite), na.rm = T))
  
  out_frame_final <- data.frame(asm = asm,
                          LCI = out_frame_lci/out_frame_length,
                          genome_breadth = out_genome_breadth/out_frame_length,
                          SDnorm = sd(out_frame$NormDepth, na.rm = T))

  out_list <- list()
  out_list[[1]] <- out_frame
  out_list[[2]] <- out_frame_final
  
  return(out_list)
}

## function to create GC vs NormDepth Boxplot:
gc_plot_detailed <- function(in_data)
{
  gc_box_plot <- ggplot(in_data, 
                          aes(factor(roundGC), NormDepth)) +
    coord_cartesian(ylim = c(0, 3.0), expand = F) + 
    geom_boxplot(aes(fill = asm), lwd = 0.05, outlier.shape = NA, alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "none") + 
    xlab('GC content') +
    ylab('Norm. Depth') + 
    scale_x_discrete(breaks= c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)) +
    scale_color_manual(values = custom_colors,
                           aesthetics = c("fill"))
  
  densNorm <- ggplot(in_data,
                      aes(NormDepth)) + 
    geom_density(aes(fill = asm), color = "black",  adjust = 3) +
    theme_void() +
    theme(
         axis.line.x = element_line(color = "black", linewidth = 0.5),
         axis.text.x = element_text(size = 8, color = "black", vjust = 0.3),
         axis.title.x = element_text(size = 8, angle = 90, vjust = 0.5),
         axis.ticks.x = element_line(color = "black", linewidth = 0.5),
         legend.position = "none") +
    #theme(legend.position = "none") + 
    geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T),
               linetype = "dashed") +
    geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T) + 
                 round(sd(in_data$NormDepth, na.rm = T),3), 
               linetype = "dashed") +
    geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T) - 
                 round(sd(in_data$NormDepth, na.rm = T),3), 
               linetype = "dashed") +
    coord_flip(xlim = c(0, 3.0), expand = F) +
    geom_area(
      aes(x = stage(NormDepth, after_stat = oob_censor(x, c(mean(in_data$NormDepth, na.rm = T) - 
                                                              round(sd(in_data$NormDepth, na.rm = T),3),
                                                            mean(in_data$NormDepth, na.rm = T) + 
                                                              round(sd(in_data$NormDepth, na.rm = T),3))))),
      stat = "density",
      adjust = 3,
      color = "grey", alpha = 0.2
    ) +
    scale_color_manual(values = custom_colors,
                           aesthetics = c("fill"))
   
  return(gc_box_plot)
}

dens_plot_detailed <- function(in_data, max_density)
{

  # gc_box_plot <- ggplot(in_data, 
  #                         aes(factor(roundGC), NormDepth)) +
  #   coord_cartesian(ylim = c(0, 3.0), expand = F) + 
  #   geom_boxplot(aes(fill = asm), lwd = 0.05, outlier.shape = NA, alpha = 0.5) +
  #   theme_bw() +
  #   theme(legend.position = "none") + 
  #   xlab('GC content') +
  #   ylab('Norm. Depth') + 
  #   scale_x_discrete(breaks= c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)) +
  #   scale_color_manual(values = custom_colors,
  #                          aesthetics = c("fill"))
  
  # densNorm <- ggplot(in_data,
  #                     aes(NormDepth)) + 
  #   geom_density(aes(fill = asm), color = "black",  adjust = 3) +
  #   #theme_void() +
  #   theme(
  # axis.text.y = element_blank(),
  # axis.title.y = element_blank(),
  # axis.ticks.y = element_blank(),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  # panel.background = element_blank(),
  # legend.position = "none",
  # ) + 
  #   geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T),
  #              linetype = "dashed") +
  #   geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T) + 
  #                round(sd(in_data$NormDepth, na.rm = T),3), 
  #              linetype = "dashed") +
  #   geom_vline(xintercept = mean(in_data$NormDepth, na.rm = T) - 
  #                round(sd(in_data$NormDepth, na.rm = T),3), 
  #              linetype = "dashed") +
  #   coord_flip(xlim = c(0, 3.0), ylim = c(0, 1.0), expand = F) +
  #   xlab("") +
  #   geom_area(
  #     aes(x = stage(NormDepth, after_stat = oob_censor(x, c(mean(in_data$NormDepth, na.rm = T) - 
  #                                                             round(sd(in_data$NormDepth, na.rm = T),3),
  #                                                           mean(in_data$NormDepth, na.rm = T) + 
  #                                                             round(sd(in_data$NormDepth, na.rm = T),3))))),
  #     stat = "density",
  #     adjust = 3,
  #     color = "grey", alpha = 0.2
  #   ) +
  #   scale_color_manual(values = custom_colors,
  #                          aesthetics = c("fill"))
                           
                           
  densNorm_q25 <- ggplot(in_data %>%
                        filter(NormDepth < 10),
                      aes(NormDepth)) + 
    geom_density(aes(fill = asm), color = "black", adjust = 3) +
            theme(
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          ) + 
    geom_vline(xintercept = median(in_data$NormDepth, na.rm = T),
               linetype = "dashed") +
    geom_vline(xintercept = quantile(in_data$NormDepth, 0.25, na.rm = T), 
               linetype = "dashed") +
    geom_vline(xintercept = quantile(in_data$NormDepth, 0.75, na.rm = T), 
               linetype = "dashed") +
    coord_flip(xlim = c(0, 3.0), ylim = c(0, max_density*1.1), expand = F) +
    xlab("") +
    geom_area(
      aes(x = stage(NormDepth, after_stat = oob_censor(x, 
                                                      quantile(in_data$NormDepth, 0.25, na.rm = T), 
                                                      quantile(in_data$NormDepth, 0.75, na.rm = T)))),
      stat = "density",
      adjust = 3,
      color = "grey", alpha = 0.1
    ) +
    scale_color_manual(values = custom_colors,
                           aesthetics = c("fill"))
   
  return(densNorm_q25)
}

## load all pandepth dataframes of a species
load_all_pandepth <- function(path) {
  print(path)
  
  ## set up variables
  tmp_path <- paste(path, "*_all.win.stat.gz", sep="/")
  in_files <- Sys.glob(tmp_path)
  pandepth_list <- list()
  lci_list <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    print(in_files[i])
    tmp_res <- read_pandepth(in_files[i], asm = pattern)
    pandepth_list[[i]] <- tmp_res[[1]]
    lci_list[[i]] <- tmp_res[[2]]
  }
  
  pandepth_list_long <- rbindlist(pandepth_list)
  LCI_breadth <- rbindlist(lci_list)
  
  ## rename columns
  LCI_breadth <- LCI_breadth %>%
    rename(polymerase=asm,
           LCI = LCI_length,
           "Genome Breadth" = genome_breadth)
  
  ## LCI/Genome Breadth plots
  LCI_plot <- ggplot(LCI_breadth, aes(polymerase, LCI)) +
    geom_col(aes(fill = polymerase)) + 
    theme_bw() +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_blank()) +
    ylim(c(0,0.5))
  
  genome_breadth_plot <- ggplot(LCI_breadth, aes(polymerase, `Genome Breadth`)) +
    geom_col(aes(fill = polymerase)) + 
    theme_bw() +
    coord_cartesian(ylim = c(0.85,1 )) +
    scale_fill_manual(values = custom_colors) +
    theme(axis.text.x = element_blank())
  
  LCI_breadth_plot <- (LCI_plot /
              genome_breadth_plot) +
              plot_layout(guides = "collect", axes = "collect")
  
  densGc <- ggplot(pandepth_list_long, aes(`GC(%)`)) + 
    geom_density(fill = "grey", color = "black") +
    theme_void() +
    geom_vline(xintercept = mean(pandepth_list_long$`GC(%)`, na.rm = T), 
               linetype = "dashed")
  
  iqr_data <- pandepth_list_long %>%
    arrange(`#Chr`, Start) %>% 
    group_by(`#Chr`, asm) %>%
    mutate(window = ceiling(row_number() / 500)) %>%
    ungroup() %>%
    group_by(`#Chr`, window, asm) %>%
    summarise(
      mean_depth = mean(NormDepth, na.rm = TRUE),
      median_depth = median(NormDepth,na.rm = TRUE), 
      iqr_depth = IQR(NormDepth,na.rm = TRUE),
      q25NormDepth = quantile(NormDepth, 0.25, na.rm = TRUE),
      q75NormDepth = quantile(NormDepth, 0.75, na.rm = TRUE)) %>%
    ungroup()
  
  chr_length <- pandepth_list_long %>%
    arrange(`#Chr`, Start) %>%
    group_by(`#Chr`) %>%
    summarise(chr_length  = max(Start)) %>%
    ungroup()
  
  longest_chr <- chr_length %>%
    filter(chr_length == max(chr_length))
  
  print("longest chromosome")
  print(longest_chr)

  #window_factor <- max(10^round(log10(longest_chr$chr_length)) / 10000, 1)
  ## window size = 200, number of windows desired = 500
  window_factor1 <- max(10^round(log10(longest_chr$chr_length)) / 200 / 200, 1)
  window_factor2 <- max(10^round(log10(longest_chr$chr_length)) / 200 / 40, 1)
  
  iqr_data_line <- pandepth_list_long %>%
    filter(`#Chr` == longest_chr$`#Chr`) %>%
    arrange(`#Chr`, Start) %>% 
    group_by(`#Chr`, asm) %>%
    mutate(window = ceiling(row_number() / window_factor1)) %>%
    ungroup() %>%
    group_by(`#Chr`, window, asm) %>%
    summarise(
      mean_depth = mean(NormDepth, na.rm = TRUE),
      median_depth = median(NormDepth,na.rm = TRUE), 
      iqr_depth = IQR(NormDepth,na.rm = TRUE),
      q25NormDepth = quantile(NormDepth, 0.25, na.rm = TRUE),
      q75NormDepth = quantile(NormDepth, 0.75, na.rm = TRUE)) %>%
    ungroup()   
    
  iqr_data_line2 <- pandepth_list_long %>%
    filter(`#Chr` == longest_chr$`#Chr`) %>%
    arrange(`#Chr`, Start) %>% 
    group_by(`#Chr`, asm) %>%
    mutate(window = ceiling(row_number() / window_factor2)) %>%
    ungroup() %>%
    group_by(`#Chr`, window, asm) %>%
    summarise(
      mean_depth = mean(NormDepth, na.rm = TRUE),
      median_depth = median(NormDepth,na.rm = TRUE), 
      iqr_depth = IQR(NormDepth,na.rm = TRUE),
      q25NormDepth = quantile(NormDepth, 0.25, na.rm = TRUE),
      q75NormDepth = quantile(NormDepth, 0.75, na.rm = TRUE)) %>%
    ungroup()
  
  iqr_plot1 <- ggplot(iqr_data_line, aes(window, median_depth)) +
    geom_line(aes(color = asm), alpha = 0.5) +
    theme_bw() + 
    scale_color_manual(values = custom_colors,
                       aesthetics = c("color")) +
    scale_y_continuous(limits = c(0, 2.0)) +
    facet_wrap(~`#Chr`)
  
  iqr_plot2 <- ggplot(iqr_data, aes(iqr_depth)) +
    geom_density(aes(fill = asm, color = asm), alpha = 0.05, adjust = 3) +
    theme_bw() + 
    scale_color_manual(values = custom_colors,
                       aesthetics = c("fill", "color")) +
    scale_x_continuous(limits = c(0, 1.5))
  
  iqr_plot3 <- ggplot(iqr_data, aes(median_depth)) +
    geom_density(aes(fill = asm, color = asm), alpha = 0.05, adjust = 3) +
    theme_bw() + 
    scale_color_manual(values = custom_colors,
                       aesthetics = c("fill", "color")) +
    scale_x_continuous(limits = c(0, 3.0))
  
  
  iqr_plot4 <- ggplot(iqr_data_line2,
                      aes(window, median_depth)) +
    geom_line(aes(color = asm), alpha = 0.5, linewidth = 0.1) +
    geom_ribbon(aes(ymin = ifelse(q25NormDepth < 0, 0, q25NormDepth), 
                    ymax = ifelse(q75NormDepth > 2.0, 2, q75NormDepth), 
                    color = asm, fill = asm),
                linewidth = 0.001, 
                alpha  = 0.15) +
    theme_bw() + 
    scale_color_manual(values = custom_colors,
                       aesthetics = c("color", "color", "fill")) +
    scale_y_continuous(limits = c(0, 2.0)) +
    facet_grid(`#Chr` ~ asm)
  
  iqr_plot <- ((iqr_plot2 + iqr_plot3) / 
                 (iqr_plot1 / iqr_plot4)) +
    plot_layout(heights = c(1,5), guides = "collect")
  
  ## Like this paper: https://www.nature.com/articles/s41467-024-51577-2#MOESM1
  
  # Source - https://stackoverflow.com/a
  # Posted by Fernando
  # Retrieved 2026-01-15, License - CC BY-SA 3.0
  
  RMSE = function(m, o){
    sqrt(mean((m - o)^2, na.rm = T))
  }
  
  cov_hetero <- iqr_data %>%
    group_by(asm) %>%
    summarise(meanNormDepth = mean(mean_depth, na.rm = T),
              coverage_heterogeneity = median(RMSE(meanNormDepth, mean_depth)))
  
  
  line_plot_data <- pandepth_list_long %>%
    group_by(roundGC, asm) %>%
    summarise(meanNormDepth = mean(NormDepth,na.rm = T),
              SDNormDepth = round(sd(NormDepth, na.rm = T),3),
              medianNormDepth = median(NormDepth, na.rm = T),
              q25NormDepth = quantile(NormDepth, 0.25, na.rm = TRUE),
              q75NormDepth = quantile(NormDepth, 0.75, na.rm = TRUE)) %>%
    ungroup()
  
  line_plot1 <- ggplot(line_plot_data, 
                       aes(roundGC, meanNormDepth)) +
    
    geom_ribbon(aes(ymin = ifelse(meanNormDepth - SDNormDepth < 0, 0, meanNormDepth - SDNormDepth), 
                    ymax = ifelse(meanNormDepth + SDNormDepth > 2.0, 2, meanNormDepth + SDNormDepth), 
                    color = asm, fill = asm), 
                linewidth = 0.001, 
                alpha  = 0.15) +
    geom_line(aes(color = asm), alpha = 0.75, size = 1) +
    scale_color_manual(values = custom_colors,
                       aesthetics = c("color", "color", "fill")) +
    scale_x_continuous(breaks= c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)) +
    scale_y_continuous(limits = c(0, 2.0)) +
    coord_cartesian(xlim = c(0, 0.95), expand = F) +
    theme_bw()
  
  line_plot2 <- ggplot(line_plot_data, 
                       aes(roundGC, medianNormDepth)) +
    geom_ribbon(aes(ymin = q25NormDepth, 
                    ymax = q75NormDepth, color = asm, fill = asm),
                linewidth = 0.001,  
                alpha  = 0.15) +
    geom_line(aes(color = asm), alpha = 0.75, size = 1) +
    scale_color_manual(values = custom_colors,
                       aesthetics = c("color", "color", "fill")) +
    scale_x_continuous(breaks= c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)) +
    scale_y_continuous(limits = c(0, 2.0)) +
    coord_cartesian(xlim = c(0, 0.95), expand = F) +
    theme_bw()
  
  sample_names <- unique(pandepth_list_long$asm)
  print(sample_names)
  str(pandepth_list_long)

  ## get maximum density for density plot display:
  max_dens_df <- pandepth_list_long %>%
    filter(NormDepth < 10)
  max_density <- max(sapply(sample_names, function(g) max(density(max_dens_df$NormDepth[max_dens_df$asm == g])$y)))

  gc_plots <- map(sample_names, ~{
    df <- pandepth_list_long %>%
      filter(asm == .x)
    out_plot <- gc_plot_detailed(df)
  })
  
  dens_plots <- map(sample_names, ~{
    df <- pandepth_list_long %>%
      filter(asm == .x)
    out_plot <- dens_plot_detailed(df, max_density)
  })
  
  all_out1 <- wrap_plots(gc_plots, nrow = length(sample_names), axes = "collect", guides = "collect")
  all_out2 <- wrap_plots(dens_plots, nrow = length(sample_names), axes = "collect", guides = "collect")
  
  all_out3 <- (all_out1 | all_out2) + plot_layout(widths = c(4,1))
  
  all_out4 <- ((densGc + plot_spacer() +
                  plot_layout(ncol = 2, width = c(4,1))) /
                 all_out3) + plot_layout(heights = c(1,10), guides = "collect")
  
  out_list <- list()
  out_list[[1]] <- all_out4
  out_list[[2]] <- line_plot1
  out_list[[3]] <- line_plot2
  out_list[[4]] <- iqr_plot
  out_list[[5]] <- LCI_breadth_plot
  
  return(out_list)
  
}

## Run and write output
out <- load_all_pandepth(input_path_pandepth)
ggsave(output_boxplot, plot = out[[1]], device = "pdf", width = 7, height = 16, units = "in" )
ggsave(output_line_mean, plot = out[[2]], device = "pdf", width = 12, height = 8, units = "in" )
ggsave(output_line_median, plot = out[[3]], device = "pdf", width = 12, height = 8, units = "in" )
ggsave(output_line_iqr, plot = out[[4]], device = "pdf", width = 12, height = 8, units = "in" )
ggsave(output_bar_lci, plot = out[[5]], device = "pdf", width = 12, height = 8, units = "in" )