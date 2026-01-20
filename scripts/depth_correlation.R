## Figure code for Figure 1 of polymerase benchmark paper
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(ggpubr)
library(ggpmisc)
library(colorspace)
library(purrr)
library(R.utils)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input
input_path_pandepth = snakemake@params[["pandepth_path"]]
contig_breaks = snakemake@input[["contig_breaks"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output
output_overall_matrix = snakemake@output[["overall_matrix"]]
output_breakpoint_matrix = snakemake@output[["breakpoint_matrix"]]
output_by_input = snakemake@output[["by_input"]]

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

## Custom sort function in order to retain the order that we use in all other plots
custom_sort <- function(x, order_vec) {
  
  ## Not necessary in cases where there is only one string per sample, but included for 
  ## strings separated by ,
  n_parts <- lengths(regmatches(x, gregexpr(",", x))) + 1
  sort_key <- n_parts
  
  ## Put Ranks to letters (see what happens with PolA and PolB)
  first_letters <- sub(",.*", "", x)
  first_rank <- match(first_letters, order_vec)
  
  second_letters <- ifelse(grepl(",", x), sub(".*,", "", x), NA_character_)
  second_rank <- match(second_letters, order_vec)

  order(sort_key, first_rank, second_rank, x)
}


## Function to read pandepth infiles
read_pandepth <- function(pandepth_file, asm = "") {
  
  read_last <- function(file) {
    tail(readLines(pandepth_file), n=1)
  }
  
  mead_depth <- sub(".*MeanDepth: *([0-9.]+).*", "\\1", read_last(pandepth_file))
  
  in_pandepth <- fread(pandepth_file, header = T)
  out_pandepth <- in_pandepth %>%
    mutate(polymerase = asm,
           MeanDepth = MeanDepth/as.numeric(mead_depth))
  
  return(out_pandepth)
}

## Function to create pairwise correlation plots
## prop will be set to 1 by default as computation should be possible,
## but if not, it can be set to lower values 
## (if only a part of the sample should be evaluated)
## limits for 2D Histogram are for now hard coded
merge_pairs <- function(i, j, prop = 1.0, corr_list = "", window_size = 200) {
    poly_names <- c(i,j)#c(names(corr_list)[i],names(corr_list)[j])
  
    df_i <- corr_list[[i]] %>%
      select( `#Chr`, Start, End, MeanDepth, polymerase) #%>%
      #rename(polymerase = MeanDepth)
    
    df_j <- corr_list[[j]] %>%
      select( `#Chr`, Start, End, MeanDepth, polymerase) #%>%
      #rename(polymerase = MeanDepth)
    
    ## Merge pandepth dataframes
    merged_ij <- full_join(df_i, df_j, by = c("#Chr", "Start", "End"))
    
    ## obsolete: overall R2 of by lm
    #model <- lm(MeanDepth.y ~ MeanDepth.x, data = merged_ij)
    #print(summary(model)$adj.r.squared)
    
    ## R2 for selected chromsome only
    #model2 <- lm(MeanDepth.y ~ MeanDepth.x, data = merged_ij %>% 
    #               filter(`#Chr` == chrom))
    #print(summary(model2)$adj.r.squared)
    
    ## function to find out color
    tile_color <- custom_colors[[which(names(custom_colors) == poly_names[1])]]
    
    corr_plot <- ggplot(merged_ij %>%
             slice_sample(prop = prop), 
           aes(MeanDepth.x, MeanDepth.y)) +
      geom_bin_2d(bins=400) +
      theme_bw() +
      scale_fill_gradient2(
        low = "grey", 
        mid = lighten(tile_color, 0.5),
        high = darken(tile_color,0.5),
        name = "Number of windows (200nt)",
        midpoint = 1000, ## for now hard-coded
        transform = "log10",
        label = comma,
        limits = c(1, 200000) ## for now hard-coded
      ) + 
      #geom_smooth(method = "lm", se=T, color="black", formula = y ~ x) +
      coord_cartesian(xlim = c(0,3.0),
                      ylim = c(0,3.0)) +
      ## R2 is not that informative I think, but we can turn this on in theory
      #stat_poly_eq(use_label(c("R2", "p"))) + ## 'eq'
      xlab(paste0("Norm. Depth (" , poly_names[1],")")) +
      ylab(paste0("Norm. Depth (" , poly_names[2],")")) +
      geom_abline(intercept= 0, slope = 1, linetype = "dashed") 
    
    out_list <- list()
    out_list[[1]] <- corr_plot
    out_list[[2]] <- merged_ij
    
    ## For now, only outplots are needed here
    return(corr_plot)
}

## do the N50 plot and summary
## parameters 'chrom' and 'window' were removed as the function can work on the whole pandepth file
## and limits are hard-coded
correlation_summary <- function(path) {
  
  ## set up variables
  tmp_path <- paste(path, "*.win.stat.gz", sep="/")
  in_files <- Sys.glob(tmp_path)
  corr_list <- list()
  
  
  ## read files (sort was added)
  for (i in seq_along(sort(in_files))) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    print(pattern)
    corr_list[[i]] <- read_pandepth(in_files[i], asm=pattern)
    names(corr_list)[i] <- pattern
  }
  
  ## get all combinations
  name_order <- names(corr_list)[custom_sort(names(corr_list), labels)]

  combs <- combn(name_order, 2, simplify = FALSE)
  
  ## Apply plotting function
  plots <- map(
    combs,
    ~ {
      i <- .x[1]
      j <- .x[2]
      ## window_size was removed here
      merge_pairs(i, j, prop = 1.0, corr_list = corr_list)
    }
  )
  
  ## Create a 'correlation matrix'
  layout_matrix <- matrix(list(plot_spacer()), 
                          nrow = length(corr_list)-1, 
                          ncol = length(corr_list)-1)
  
  max_j <- length(corr_list)-1
  max_i <- length(corr_list)-1
  plot_counter <- 1
  
  for (j in 1:max_j) {
    i <- 1
    
    while (i <= max_i) {
      layout_matrix[i, j][[1]] <- plots[[plot_counter]]
      i <- i + 1
      plot_counter <- plot_counter + 1
      }
    max_i <- max_i -1
    }
  
  # Assemble patchwork from matrix
  patchwork <- wrap_plots(layout_matrix, guides = 'collect')
  
  return(patchwork)
  
}

contig_break_matrix <- function(path, contig_break_path) {

  ## set up variables (like in previous function)
  tmp_path <- paste(path, "*.win.stat.gz", sep="/")
  in_files <- Sys.glob(tmp_path)
  corr_list <- list()
  
  
  ## read files (sort was added)
  for (i in seq_along(sort(in_files))) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    print(pattern)
    corr_list[[i]] <- read_pandepth(in_files[i], asm=pattern)
    names(corr_list)[i] <- pattern
  }
  
  ## Join all dataframes:
  corr_list_frame <- imap(corr_list, ~ {
      prefix <- unique(.x$polymerase)[1]
      
      .x %>%
        rename_with(~ paste0(prefix, "_", .x), .cols = -all_of(c("#Chr", "Start", "End", "GC(%)")))
    }) 

  ## Change into Dataframe
  corr_list_merged <- corr_list_frame %>%
    purrr::reduce(full_join, by = c("#Chr", "Start", "End","GC(%)"))
    
  ## Load contig break bed file
  contig_breaks <- read_delim(contig_break_path, col_names = TRUE)
  
  
  ####################
  ##  Break Matrix  ##
  ####################
  
  ## Create a small and a large version of windows (one for 'correlation like' matrix, one for summarized windows)
  ## Plot the same with contig break windows
  ## Zoom into clear break windows
  ## all_merged comes from calculate_contig_break
  contig_breaks_small <- contig_breaks %>%
    mutate(start_c = start + 750,
           end_c = ifelse(end - 750 <= start_c, end, end - 750)) %>%
    filter(end_c - start_c > 0  & end_c - start_c < 10000) %>%
    mutate(window_id = row_number())
  
  ## Give windows of +- 10000 to contig breaks
  contig_breaks_big <- contig_breaks %>%
    mutate(start_c = start - 10000,
           end_c = end + 10000) %>%
    filter(end_c - start_c > 0  & end_c - start_c < 100000) %>%
    mutate(window_id = row_number())

  ## Merge/Intersect with the data.table functionality against pandepth (corr_list_merged)
  setDT(contig_breaks_small)
  setDT(contig_breaks_big)
  setDT(corr_list_merged)
  
  ## merging datatables while retaining column names 
  ## see here: https://stackoverflow.com/questions/41637736/non-equi-joins-adding-all-columns-of-range-table-in-data-table-in-one-step
  ## Filtering here can be more stringent but shouldnt change drastically: chr = chr, Start <= end_c,  Start >= start_c, End >= start_c, End <= end_c
  merged_breaks_big <- corr_list_merged[contig_breaks_big, 
                               on = .(`#Chr` = chr, Start <= end_c, End >= start_c), 
  mget(c(paste0("x.", names(corr_list_merged)), paste0("i.", names(contig_breaks_big)))) %>% 
    set_names(c(names(corr_list_merged), names(contig_breaks_big)))]
    
    
  merged_breaks_small <- corr_list_merged[contig_breaks_small, 
                               on = .(`#Chr` = chr, Start <= end_c,  Start >= start_c,
                                      End >= start_c, End <= end_c), 
  mget(c(paste0("x.", names(corr_list_merged)), paste0("i.", names(contig_breaks_small)))) %>% 
    set_names(c(names(corr_list_merged), names(contig_breaks_small)))]
  

  ## Get polymerase names
  #str(contig_breaks_small)
  dropouts_flt <- merged_breaks_small %>%
    select(chr, Start, End, ends_with("MeanDepth"),)
  
  dropout_names <- dropouts_flt %>% 
    select(ends_with("MeanDepth")) %>%names()
    
  dropout_names_order <- dropout_names[custom_sort(dropout_names, labels)]
  
  combs <- combn(dropout_names_order, 2, simplify = FALSE)
  
  
  plots <- map(combs, ~ {
    col1 <- sym(.x[1])
    col2 <- sym(.x[2])
    tile_color <- custom_colors[which(names(custom_colors) == str_split_i(.x[1], "_", 1))]
    other_color <- custom_colors[which(names(custom_colors) == str_split_i(.x[2], "_", 1))]
    (merged_breaks_small %>% 
        select(chr, Start, End, !!col1, !!col2) %>%
        filter(!!col2 < 0.1 & !!col1 > 0.1) %>%
        select(!!col1) %>%
        ggplot(aes(!!col1)) +
        geom_density(fill = tile_color, alpha = 0.3) +
        theme_void() +
        coord_cartesian(xlim = c(0,2))) + 
      plot_spacer() + 
    (merged_breaks_small %>% 
      select(chr, Start, End, !!col1, !!col2) %>%
      ggplot(aes(!!col1, !!col2)) +
      geom_bin_2d(bins = 75) +
      scale_fill_gradient2(
        low = "grey", #lighten(tile_color, 0.9),
        mid = lighten(tile_color, 0.5),
        high = darken(tile_color,0.5),
        name = "Number of windows (200nt)",
        midpoint = 10,
        transform = "log10",
        label = comma,
        limits = c(1, 10000)
      ) +
      #stat_poly_eq(use_label(c("R2", "p"))) + 
      coord_cartesian(xlim = c(0,2), ylim = c(0,2)) +
      theme_bw() +
      geom_abline(intercept= 0, slope = 1, linetype = "dashed") +
      geom_vline(xintercept = 0.1, color = "black", linetype = "dotted") +
      geom_hline(yintercept = 0.1, color = "black", linetype = "dotted")) +
    (merged_breaks_small %>% 
       select(chr, Start, End, !!col1, !!col2) %>%
       filter(!!col1 < 0.1 & !!col2 > 0.1) %>%
       select(!!col2) %>%
       ggplot(aes(!!col2)) +
       geom_density(fill = other_color, alpha = 0.3) +
       theme_void()+
       coord_flip(xlim = c(0,2))) +
      plot_layout(
        ncol = 2, 
        nrow = 2, 
        widths = c(4, 1),
        heights = c(1, 4)
      )
  }
  )
  
  layout_matrix <- matrix(list(plot_spacer()), 
                          nrow = length(dropout_names_order)-1, 
                          ncol = length(dropout_names_order)-1)
  
  max_j <- length(dropout_names_order)-1
  max_i <- length(dropout_names_order)-1
  plot_counter <- 1
  
  for (j in 1:max_j) {
    i <- 1
    
    while (i <= max_i) {
      layout_matrix[i, j][[1]] <- plots[[plot_counter]]
      i <- i + 1
      plot_counter <- plot_counter + 1
    }
    max_i <- max_i -1
  }
  
  # Assemble patchwork from matrix
  patchwork <- wrap_plots(layout_matrix, guides = 'collect')
  
  ###################
  ## Window plots  ##
  ###################
  
  merged_breaks_big <- as_tibble(merged_breaks_big)
  
  ## Which windows are displayed is still hard-coded, but will probably be changed
  merged_breaks_big_long <- merged_breaks_big %>%
    group_by(window_id) %>%
    mutate(norm_start = Start - min(Start),
           norm_end = norm_start + 200) %>%
    pivot_longer(cols = ends_with("MeanDepth"), 
                  names_to = "Depth_polymerase",
                  values_to = "MeanDepth") %>%
    ungroup() %>%
    group_by(polymerase, norm_start, Depth_polymerase) %>%
    mutate(MeanDepth_mean = mean(MeanDepth),
           MeanDepth_sd = sd(MeanDepth),
           Mean_GC = mean(`GC(%)`),
           SD_GC = sd(`GC(%)`),
           centered = norm_start - 12000) %>%
    ungroup() %>%
    group_by(polymerase) %>%
    filter(n() > 2500, norm_start < 24000)
  
  ## only plot single drop-outs for now
  facets <- input_names
  
  ## GC window plots (all windows stacked on top of each other)
  GC_windows <- map(facets, ~ {
    current_break <- .x[1]
    
    ## For display on top of plot
    number_of_breaks <- contig_breaks %>%
      group_by(polymerase) %>%
      summarise(dropouts = n()) %>%
      filter(polymerase == current_break)
  
    gc_df <- merged_breaks_big_long %>%
      filter(polymerase == current_break) %>%
      mutate(polymerase_lib = str_split_i(Depth_polymerase, "_", 1)) %>%
      group_by(polymerase_lib, centered) %>%
      summarise(centered = max(centered),
                Mean_GC = median(Mean_GC, na.rm =T),
                SD_GC = median(SD_GC, na.rm = T))
    
    ggplot(gc_df, aes(centered,Mean_GC)) +
      geom_line(color = "black", linewidth = 0.75) +
      geom_ribbon(aes(ymin = ifelse(Mean_GC - SD_GC < 10, 10, Mean_GC - SD_GC), 
                     ymax = ifelse(Mean_GC + SD_GC > 70, 70, Mean_GC + SD_GC)), color = "grey",
                  alpha  = 0.1) +
      theme_void() +
      theme(
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(size = 8, color = "black", vjust = 0.3),
        axis.title.y = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.y = unit(0.1, "cm"),
        plot.title = element_text(size = 11)
      ) +
      scale_y_continuous(limits = c(5,75)) +
      ylab("Mean GC(%)") +
      ggtitle(paste0(current_break, "\n(n = ", number_of_breaks$dropouts , ")"))
  }
  )

  ## Coverage window plots (all windows stacked on top of each other)
  coverage_windows <- map(facets, ~ {
    current_break <- .x[1]
    
    
    coverage_df <- merged_breaks_big_long %>%
      filter(polymerase == current_break) %>%
      mutate(polymerase_lib = str_split_i(Depth_polymerase, "_", 1)) %>%
      group_by(polymerase_lib, centered) %>%
      summarise(centered = max(centered),
                MeanDepth_mean = median(MeanDepth_mean,na.rm = T),
                SD_mean = median(MeanDepth_sd, na.rm = T))
  
    ggplot(coverage_df, 
                    aes(centered,MeanDepth_mean)) +
      geom_line(aes(color = polymerase_lib), linewidth = 0.75) +
      geom_ribbon(aes(ymin = ifelse(MeanDepth_mean - SD_mean <0, 0, MeanDepth_mean - SD_mean), 
                      ymax = MeanDepth_mean + SD_mean, fill = polymerase_lib), 
                  alpha  = 0.1) +
      scale_color_manual(values = custom_colors,
                         aesthetics = c("color", "fill")) +
      theme_bw() +
      ylab("Mean Depth") +
      xlab("Position from contig break") +
      geom_vline(xintercept = 0) +
      scale_y_continuous(limits = c(0,2))
  }
  )
  
  ## Create faceted plot 'semi-manually'
  n_plots <- length(coverage_windows)
  chunk_size <- 4
  n_chunks <- ceiling(n_plots / chunk_size)
  
  ## Generate chunks automatically in groups of 4 (chunk size must be determined by lenght of labels in the future)
  chunks <- split(1:n_plots, ceiling(seq_along(1:n_plots) / chunk_size))

  ## Create pairs of plots
  row_pairs <- map(chunks, \(i) {
    top_row <- (wrap_plots(GC_windows[i], nrow = 1) + plot_layout(axes = "collect"))
    bot_row <- (wrap_plots(coverage_windows[i], nrow = 1) + plot_layout(axes = "collect"))
    (top_row / bot_row) + plot_layout(heights = c(1,1))
  })
  
  ## Final faceted plot
  facet_plot <- wrap_plots(row_pairs, nrow = n_chunks, guides = "collect")
  
  
  out_list <- list()
  out_list[[1]] <- patchwork
  out_list[[2]] <- facet_plot
  
  return(out_list)
  
}

## Create plots
## Big matrix plot
big_matrix <- correlation_summary(input_path_pandepth)

## Contig Break matrix plot and Contig Break windows
contig_break_plots <- contig_break_matrix(input_path_pandepth, contig_breaks)

## write output
ggsave(output_overall_matrix, plot = big_matrix, device = "pdf", width = 14, height = 14, units = "in")
ggsave(output_breakpoint_matrix, plot = contig_break_plots[[1]], device = "pdf", width = 14, height = 14, units = "in")
ggsave(output_by_input, plot = contig_break_plots[[2]], device = "pdf", width = 10, height = 10, units = "in")
