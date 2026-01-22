## Figure code for Figure 1 of polymerase benchmark paper
library(tidyverse)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(scales)
library(patchwork)

safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
## Input for plotting/sumstats
input_path_faidx = snakemake@params[["faidx_path"]]
input_path_compleasm = snakemake@params[["compleasm_path"]]
compleasm_database_name = snakemake@params[["compleasm_database"]]
input_path_merqury = snakemake@params[["merqury_path"]]
input_path_hifieval = snakemake@params[["hifieval_path"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output (N(x) + compleasm tables, final plot)
output_ng_table = snakemake@output[["ng_table"]]
output_compleasm_table = snakemake@output[["compleasm_table"]]
output_final_figure = snakemake@output[["final_plot"]]
out_hifieval_table = snakemake@output[["hifieval_table"]]
out_merqury_table = snakemake@output[["merqury_table"]]

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


## format labels for N50 plot
format_labels <- function(x) {
  sapply(x, function(val) {
    if (val >= 1e9){
      paste0(val / 1e9, "Gb")
    }  
      else if (val >= 1e6) {
      paste0(val / 1e6, "Mb")
    } else if (val >= 1e3) {
      paste0(val / 1e3, "Kb")
    } else {
      paste0(val, "B")
    }
  })
}

## find breakpoints
get_breaks <- function(min, max, steps = c(1,2.5,5)) {
  logstart <- floor(log10(min))
  logend <- ceiling(log10(max))
  
  breaks <- c()
  for (e in logstart:logend) {
  breaks <- c(breaks, steps*10^e)
  }
  
  labels <- format_labels(breaks)
  
  break_list <- list()
  break_list[[1]] <- breaks
  break_list[[2]] <- labels
  
  return(break_list)
}

## contig N50 plots for all polymerases
## load chrom size plots
prepare_n50 <- function(chrom_file_path = "", asm = "") {
  chrom_file <- read_delim(chrom_file_path, col_names = F)
  
  chrom_file_ord <- chrom_file %>%
    select(X1,X2) %>%
    arrange(desc(X2)) %>%
    rename(length = X2,
           chr_ID = X1) %>%
    mutate(cumul_size = cumsum(length),
           n_val = round((cumul_size/sum(length) *100), 3),
           order_val = as.integer((row_number()) )) %>%
    rename(setNames("cumul_size", paste("cumul_size", asm, sep = ".")),
           setNames("length", paste("length", asm, sep = ".")),
           setNames("n_val", paste("n_val", asm, sep = ".")),
           setNames("chr_ID", paste("chr_ID", asm, sep = ".")))
  
  return(chrom_file_ord)
}

## do the N50 plot and summary
n50_summary <- function(path) {
  
  ## set up variables
  tmp_path <- paste(path, "*.fa.fai", sep="/")
  in_files <- Sys.glob(tmp_path)
  chrom_list <- list()
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    
    if (!is.na(pattern)) {
    chrom_list[[i]] <- prepare_n50(in_files[i], asm=pattern)
  } 
  }
   
  chrom_list <- chrom_list[lengths(chrom_list) > 0]
  
  ## convert to dataframe
  chrom_df <- chrom_list %>% 
    reduce(full_join, by='order_val') %>%
    pivot_longer(-order_val, 
                 names_to = c(".value", "Assembly"), 
                 names_sep='\\.')
                 
                 
  ## get contig N50 values   
  n50_info <- chrom_df %>%
    mutate(n50_proxy = n_val - 50) %>%
    group_by(Assembly) %>%
    filter(n50_proxy > 0) %>%
    filter(n50_proxy  == min(n50_proxy, na.rm = T)) %>%
    ungroup()
    
  ## get top 10 assemblies
  top_10 <- n50_info %>%
    slice_max(order_by = length, n = 10)
  
  print(top_10)
  
  chrom_df <- chrom_df %>%
  filter(Assembly %in% top_10$Assembly)          
                                
                 
  ## get break points for plot
  break_list <- get_breaks(min(n50_info$length), max(n50_info$length))
  
  ## create N50 plot
  n50_graph_plot <- ggplot(chrom_df, aes(n_val, length, color = Assembly)) + #,size = Assembly 
    geom_step(size = 0.4) +
    theme_bw() +
    theme(axis.text= element_text(size = 15),
          axis.title= element_text(size = 18),
          legend.position = "none") + 
    scale_y_log10(
      breaks = break_list[[1]],
      label = break_list[[2]],
      limits = c(1000,max(chrom_df$length))) +
    annotation_logticks(sides = "l") +  
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) ) +
    ylab("Contig length (bp)") +
    xlab("N(x) value") +
    geom_vline(xintercept = 50, color = "red",linetype = "dashed") +
    scale_color_manual(values = custom_colors)
  
  
  ## write output list 
  n50_out_list <- list()
  n50_out_list[[1]] <- n50_info
  n50_out_list[[2]] <- n50_graph_plot
  n50_out_list[[3]] <- top_10
  return(n50_out_list)
  
}

## compleasm summary
read_compl <- function(path, name) {
  compl <- read_delim(path, col_names = F)
  compl <- compl %>%
    rename(BUSCO_class = X1) %>%
    mutate(polymerase = name)
  return (compl)
}

## create summary stats and output plot
compleasm_summary <- function(path, compleasm_label, top_10) {
  
  ## set up variables
  tmp_path <- paste(path, "*rf.txt", sep="/")
  in_files <- Sys.glob(tmp_path)
  compl_list <- list()
  
  ## create data table
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    if (!is.na(pattern)) {
    compl_list[[i]] <- read_compl(in_files[i], name=pattern)
  }
  }
  
  
  compl_list <- compl_list[lengths(compl_list) > 0]
  
  #compl_list <- compl_list[lengths(compl_list) > 0]
  compl_all <- rbindlist(compl_list)
  compl_all <- compl_all %>%
  filter(polymerase %in% top_10$Assembly)
  
  compl_all$category <- factor(compl_all$BUSCO_class, levels = c("M","I","F","D","S"), ordered = TRUE)
  
  ## color scheme
  custom_colors_s <- c("S" = "#274655",
                       "D" = "#2A9D8E",
                       "F" = "#E9C46B",
                       "I" = "#E9C48D",
                       "M" = "#E56F52")
                    
                    
  compl_all <- compl_all %>%
  mutate(wrapped_polymerase = str_replace_all(polymerase, "PLUS", "\\+\n")) 
  
  sel_order_compl <- 
    compl_all%>% 
    group_by(wrapped_polymerase) %>%
    mutate(
      "SD" = sum(X3[category %in% c("S", "D")])
    ) %>%
    arrange(SD) %>% 
    filter(category == "S") %>% 
    mutate(wrapped_polymerase = factor(wrapped_polymerase))
                       
  compl_all_mut <- compl_all%>% 
    mutate(wrapped_polymerase = factor(wrapped_polymerase, levels = sel_order_compl$wrapped_polymerase, ordered = TRUE))
  
  compleasm_plot <- ggplot(compl_all_mut, aes(X3, wrapped_polymerase, fill = category)) +
    geom_col() +
    theme_bw() +
    theme(axis.text= element_text(size = 6),
          axis.title= element_text(size = 12))+
    scale_fill_manual(name = "compleasm status", values = custom_colors_s, 
                      labels =  c("missing", "fragmented(split)", "fragmented(partial)", "complete (duplicated)", "complete (single)")) +
    xlab(compleasm_label) +
    ylab("")
  
  ## Compleasm Summary
  compl_all_out_table <- compl_all_mut %>%
    group_by(polymerase) %>%
    pivot_wider(names_from=BUSCO_class,values_from = c(X3,X2))%>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
    mutate(sum_intact = X3_S + X3_D,
           perc_intact = X2_S + X2_D,
           sum_all = X3_S + X3_D + X3_F + X3_I + X3_M) %>%
    select(-category, -wrapped_polymerase)

  compleasm_out_list <- list()
  
  compleasm_out_list[[1]] <- compl_all_out_table
  compleasm_out_list[[2]] <- compleasm_plot
  
  ## return outlist
  return(compleasm_out_list)
  
}


output_hifieval_readstats <- function(path, top_10)
{
  ## set up variables
  tmp_path <- paste(path, "*.summary.tsv", sep="/")
  in_files <- Sys.glob(tmp_path)
  hifieval_list <- list()
  
  read_hifieval <- function(path, name) {
    hifieval_reads <- fread(path, header = T)
    hifieval_reads <- hifieval_reads %>%
      mutate(polymerase = name)
    return(hifieval_reads)
  }
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
    if (!is.na(pattern)){
    hifieval_list[[i]] <- read_hifieval(in_files[i], name=pattern)
    }
  }
  
  ## convert to dataframe
  hifieval_list <- hifieval_list[lengths(hifieval_list) > 0]
  hifieval_df <- rbindlist(hifieval_list)


  hifieval_df <- hifieval_df %>%
  filter(polymerase %in% top_10$Assembly)
  
  ## collapse df
  hifieval_df_sum <- hifieval_df %>%
    group_by(polymerase) %>%
    summarise(corrected_bases = sum(num_cc)/sum(raw_end - raw_start),
              undercorrected_bases = sum(num_uc)/sum(raw_end - raw_start),
              overcorrected_bases = sum(num_oc)/sum(raw_end - raw_start)) %>%
    pivot_longer(cols = c(corrected_bases,
                          undercorrected_bases,
                          overcorrected_bases), 
                 names_to = "correction_class",
                 values_to = "fraction") 
  
  ## Adjust labels
  hifieval_df_sum <- hifieval_df_sum %>%
    mutate(wrapped_polymerase = str_replace_all(polymerase, "PLUS", "\\+\n"))
  
  ## plot
  hifieval_plot <- ggplot(hifieval_df_sum %>%
                              filter(correction_class %in% c("corrected_bases", "undercorrected_bases")),
                              aes(polymerase, fraction, fill = polymerase)) +
    geom_col() +
    #geom_col(position = "dodge2") +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = custom_colors,
    name = "polymerase",
    labels = unique(hifieval_df_sum$wrapped_polymerase))  +
    theme_bw() +
    theme(axis.text.x = element_blank()) + 
    ylab("% of corrected bases")
   
  hifieval_out <- list()
  
  ## Write output
  hifieval_out[[1]] <- hifieval_plot 
  hifieval_out[[2]] <- hifieval_df_sum
  
  return(hifieval_out)
}


merqury_asm_sum <- function(path, top_10)
{


  read_merqury_oneline <- function(path, name) {
    merqury <- read_delim(path, col_names = F)
    merqury <- merqury %>%
      mutate(polymerase = name)
    return(merqury)
  }
  
  ## set up variables
  tmp_path_qv <- paste(path, "*/*slf.qv", sep="/")
  in_files_qv <- Sys.glob(tmp_path_qv)
  tmp_path_com <- paste(path, "*/*completeness.stats", sep="/")
  in_files_com <- Sys.glob(tmp_path_com)
  merqury_list_qv <- list()
  merqury_list_com <- list()
  
  ## create data table
  for (i in seq_along(in_files_qv)) {
    pattern <- str_extract(in_files_qv[i], paste(labels, collapse="|"))
    if (!is.na(pattern)) {
    merqury_list_qv[[i]] <- read_merqury_oneline(in_files_qv[i], name=pattern)
  }
  }
  
  for (i in seq_along(in_files_com)) {
    pattern <- str_extract(in_files_com[i], paste(labels, collapse="|"))
    if (!is.na(pattern)) {
    merqury_list_com[[i]] <- read_merqury_oneline(in_files_com[i], name=pattern)
  }
  }
  
  merqury_list_qv <- merqury_list_qv[lengths(merqury_list_qv) > 0]
  mq_qv <- rbindlist(merqury_list_qv)
  
  merqury_list_com <- merqury_list_com[lengths(merqury_list_com) > 0]
  mq_com <- rbindlist(merqury_list_com)
  
  merqury_all <- mq_qv %>%
  left_join(mq_com, by = c("X1" = "X1")) %>%
  rename(polymerase = X1)
  
  merqury_all <- merqury_all %>%
  filter(polymerase %in% top_10$Assembly)
    
  er <- ggplot(merqury_all, aes(polymerase, X5.x)) +
  geom_col(aes(fill = polymerase)) +
  #ggtitle("Merqury assembly error rate") +
  ylab("error rate") + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
    theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 10),
        axis.title= element_text(size = 11))

  qv <- ggplot(merqury_all, aes(polymerase, X4.x)) +
  geom_col(aes(fill = polymerase)) +
  coord_cartesian(ylim = c(min(merqury_all$X4.x)-5, 60)) +
  #ggtitle("Merqury assembly quality value") +
  ylab("QV") + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
    theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 10),
        axis.title= element_text(size = 11))

  comp <- ggplot(merqury_all, aes(polymerase, X5.y)) +
  geom_col(aes(fill = polymerase)) +
  coord_cartesian(ylim = c(min(merqury_all$X5.y)-2, 100)) +
  #ggtitle("Merqury assembly completeness") +
  ylab("Completeness (%)") +
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 10),
        axis.title= element_text(size = 11))
  
  out_plot1 <- (er / qv /comp) +
  plot_layout(guides = "collect",
              axes = "collect") & theme(legend.position = "none")
              
  out_plot2 <- (er / qv /comp) +
  plot_layout(guides = "collect",
              axes = "collect")
              
  out_list <- list()
  out_list[[1]] <- out_plot1
  out_list[[2]] <- merqury_all
  out_list[[3]] <- out_plot2
  
  return(out_list)
}

## final output

## N(x) plot (and table)
n50_plot <- n50_summary(input_path_faidx)
top_10 <- n50_plot[[3]]

## compleasm plot (and table)
compleasm_plot <- compleasm_summary(input_path_compleasm, 
                                          compleasm_database_name, top_10)

## hifieval error stats
if (!is.null(input_path_hifieval)) {
hifieval_out <- output_hifieval_readstats(input_path_hifieval, top_10)
}

## Merqury output plot
merqury_out <- merqury_asm_sum(input_path_merqury, top_10)

## Final Plot
if (!is.null(input_path_hifieval)) {
final_plot <- (
  n50_plot[[2]] /
  compleasm_plot[[2]] /
  merqury_out[[1]] /
   hifieval_out[[1]]) +
   plot_layout(guides = "collect", heights = c(3,2,2,2)) + 
  plot_annotation(tag_levels = 'A') 
  } else { ## Final Plot without hifieval
final_plot <- (
  n50_plot[[2]] /
  compleasm_plot[[2]] /
  merqury_out[[3]]) +
   plot_layout(guides = "collect", heights = c(3,2,2)) + 
  plot_annotation(tag_levels = 'A')
}

## write output
ggsave(output_final_figure, plot = final_plot, device = "pdf", width = 12, height = 14, units = "in")
write_delim(n50_plot[[1]], output_ng_table, delim = "\t")
write_delim(compleasm_plot[[1]], output_compleasm_table, delim = "\t")
write_delim(merqury_out[[2]], out_merqury_table, delim = "\t")
if (!is.null(input_path_hifieval)) {
  write_delim(hifieval_out[[2]], out_hifieval_table, delim = "\t")
  }
