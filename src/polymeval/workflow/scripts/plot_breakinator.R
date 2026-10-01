## proportion of chimeric alignments
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)
library(data.table)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input
input_path_chimeras = snakemake@params[["chimera_path"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output
output_plot = snakemake@output[["plot"]]

## set colors for names
input_names <- c(input_names)
input_names <- unlist(strsplit(input_names, split = ","))
if (!is.null(input_names)){
  labels <- sort(input_names)
  if (!is.null(in_colors)) {
  col_dict <- read_delim(in_colors, col_names = FALSE)
  custom_colors <- setNames(col_dict$X2, col_dict$X1)
  } else {
  if (length(input_names) > 12) {
  palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(labels))
  } else {
  palette_colors <- safe
  }
  custom_colors <- setNames(palette_colors, labels)
}
}


## Single breakinator file reading
breakinator_file <- function(path, lib = "") {
    tmp <- read_delim(path, col_names = T)
    out <- tmp %>%
    mutate(sample = lib,
            fold_reads_p = Fold_reads/Uniq_artifact_reads,
            fold_breaks_p = Fold_breaks/all_break,
            chim_reads_p = Chim_reads/Uniq_artifact_reads,
            chim_breaks_p = Chim_breaks/all_break) %>%
    pivot_longer(cols = c("#Reads_passed", 
                            "all_break", 
                            "Uniq_artifact_reads",
                            "Fold_reads", 
                            "Fold_breaks",
                            "Chim_reads",
                            "Chim_breaks",
                            "fold_reads_p",
                            "fold_breaks_p",
                            "chim_reads_p",
                            "chim_breaks_p"), 
                names_to = "read_feature", 
                values_to = "values") %>%
    select(sample, read_feature, values)

    return(out)
}

## Load and plot breakinator
breakinator_plotting <- function(path) {
  
  ## set up variables
  tmp_path <- paste(path, "*_breakinator_summary.txt", sep="/")
  in_files <- Sys.glob(tmp_path)
  breakinator_list <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    breakinator_list[[i]] <- breakinator_file(in_files[i], lib=pattern)
  }
  
  out_df <- rbindlist(breakinator_list)
  
  abs_plot <- ggplot(out_df %>% filter(read_feature %in% 
                         c("#Reads_passed","Uniq_artifact_reads","Fold_reads","Chim_reads")), 
       aes(sample, values))+
  geom_col(aes(fill = sample), position = "dodge2") +
  facet_wrap(~factor(read_feature, levels = 
                       c("#Reads_passed","Uniq_artifact_reads","Fold_reads","Chim_reads")),
             nrow = 1)+
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() + 
    theme(axis.text.x = element_blank())

  rel_plot <- ggplot(out_df %>% filter(read_feature %in% 
                         c("fold_reads_p",
                           "fold_breaks_p",
                           "chim_reads_p",
                           "chim_breaks_p")), 
       aes(sample, values))+
  geom_col(aes(fill = sample), position = "dodge2") +
  facet_wrap(~factor(read_feature, levels = 
                       c("fold_reads_p",
                         "fold_breaks_p",
                         "chim_reads_p",
                         "chim_breaks_p")),
             nrow = 1)+
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() + 
  theme(axis.text.x = element_blank()) 
  
  out_plot <- (abs_plot / rel_plot) +
    plot_layout(axes = "collect", guides = "collect")

  return(out_plot)
}

breakinator_plot <- breakinator_plotting(input_path_chimeras)

ggsave(output_plot, plot = breakinator_plot, device = "pdf", width = 12, height = 14, units = "in")