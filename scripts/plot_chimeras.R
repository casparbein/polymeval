## proportion of chimeric alignments
library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)
library(RColorBrewer)
library(scales)
library(ggtext)

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
  if (in_colors != "") {
  col_dict <- read_delim(in_colors, col_names = FALSE)
  custom_colors <- setNames(col_dict$X2, col_dict$X1)
  } else {
  palette_colors <- safe
  custom_colors <- setNames(palette_colors, labels)
}
}


### updated mapping type plot
decode_sam_flag <- function(flag) {
  meanings <- c(
    "0" = "primary alignment",
    "4" = "read unmapped",
    "16" = "primary alignment",
    "256" = "non-primary alignment",
    "272" = "non-primary alignment",
    "2048" = "suppl. (chimeric) alignment",
    "2064" = "suppl. (chimeric) alignment"
  )
  
  return(meanings[as.character(flag)])
}

## read in mapping types from chimera operations
map_file <- function(path, lib = "") {
  map_types <- read_delim(path, col_names = F)
  
  map_types <- map_types %>%
  mutate(description = sapply(X1, decode_sam_flag)) %>%
  select(-X1) %>%
  group_by(description) %>%
  summarize(
    total_count = sum(X2)
  ) %>%
  ungroup() %>%
  mutate(percentage = total_count/sum(total_count),
         library = lib)
  
  return(map_types)
}

## MappingFlags
read_mappings <- function(path) {
  
  ## set up variables
  tmp_path <- paste(path, "*_mappingType.txt", sep="/")
  in_files <- Sys.glob(tmp_path)
  map_list <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    map_list[[i]] <- map_file(in_files[i], lib=pattern)
  }
  
  out_df <- rbindlist(map_list)
  
  out_plot <- ggplot(out_df, aes(description, percentage)) +
    geom_col(aes(fill = library), position = "dodge2") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0.1, 1, by = 0.1), labels = comma) + 
    theme_bw() +
    guides(color = guide_legend(byrow = TRUE, ncol = 2)) +
    xlab("Alignment type") +
    ylab("fraction of alignments") +
    theme(legend.text = element_markdown(size = 10),
          legend.title = element_markdown(size = 12),
          axis.text= element_text(size = 12),
          axis.title= element_text(size = 15),
          legend.spacing.y = unit(0.2, 'cm'),
          legend.background = element_rect(color = "black"),
          plot.title = element_text(size = 15)) +
    scale_fill_manual(values = custom_colors)
  
  return(out_plot)
}

## Plot histograms plus mapped read lengths:
load_hist <- function(in_path, name, read_type) {
    hist <- read_delim(in_path, col_names = F)
    hist <- hist %>% 
      mutate(polymerase = name,
             read_type = read_type)
    
    return(hist)
}


## create plots for all of them:
chimera_histograms <- function(path) {
  
  ## set up variables
  tmp_path1 <- paste(path, "*_prim.hist.txt", sep="/")
  tmp_path2 <- paste(path, "*_suppl.hist.txt", sep="/")
  tmp_path3 <- paste(path, "*_rl.hist.txt", sep="/")
  in_files1 <- Sys.glob(tmp_path1)
  in_files2 <- Sys.glob(tmp_path2)
  in_files3 <- Sys.glob(tmp_path3)
  chimera_list <- list()
  
  ## read files
  for (i in seq_along(in_files1)) {
    pattern <- str_extract(in_files1[i], paste(labels, collapse="|"))
    raw <- load_hist(in_files3[i], pattern, "raw")
    prim <- load_hist(in_files1[i], pattern, "prim")
    supp <- load_hist(in_files2[i], pattern, "supp")
    all <- rbind(raw, prim, supp)
    
    chimera_list[[i]] <- all
  }
  
  out_df <- rbindlist(chimera_list)
  out_plot <- ggplot() +
    geom_col(out_df %>%
               filter(read_type == "raw"), 
             mapping = aes(X1, X2, fill= polymerase) ) +
    geom_line(out_df %>%
                filter(read_type == "supp" | read_type == "prim"),  
              mapping = aes(X1, X2, color= read_type), 
              alpha = 0.5) +  # Line
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(name = "alignment type", values = c("red", "blue"), 
                       labels = c("primary", "supplementary (chimeric)")) +
    xlim(c(1,20000)) +
    ylim(c(0,max(out_df$X2 + 200))) +
    theme_bw() +
    guides(fill = "none") +
    theme(legend.position = "none") + 
    xlab("read/alignment length") +
    ylab("count") +
    facet_wrap(~polymerase, ncol = 1)
    
    return(out_plot)
}
  

## Get output
maps <- read_mappings(input_path_chimeras)
hists <-chimera_histograms(input_path_chimeras)
combined_plot <- (maps / hists) +
  plot_layout(heights = c(1,3), guides = "collect", axes = "collect")
  
ggsave(output_plot, plot = combined_plot, device = "pdf", width = 12, height = 14, units = "in")
  
