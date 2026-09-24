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
input_path_seqkit_stats = snakemake@input[["seqkit"]]
input_path_merqury = snakemake@params[["merqury_path"]]
input_path_hifieval = snakemake@params[["hifieval_path"]]
in_colors = snakemake@params[["colors"]]
input_names = snakemake@params[["asm_ids"]]
sample_names = snakemake@params[["sample_names"]]

## output (N(x) + compleasm tables, final plot)
output_ng_table = snakemake@output[["ng_table"]]
output_compleasm_table = snakemake@output[["compleasm_table"]]
output_final_figure = snakemake@output[["final_plot"]]
out_hifieval_table = snakemake@output[["hifieval_table"]]
out_merqury_table = snakemake@output[["merqury_table"]]

## set colors for names
input_names <- unlist(strsplit(input_names, split = ","))
sample_names <- unlist(strsplit(c(sample_names), split = ","))

## one row per assembly: which sample it came from, and which assembler built it
id_table <- tibble(asm_id = input_names) %>%
  mutate(sample    = if_else(str_detect(asm_id, "__"),
                             str_remove(asm_id, "__[^_]+$"), asm_id),
         assembler = if_else(str_detect(asm_id, "__"),
                             str_extract(asm_id, "[^_]+$"), "hifiasm"))
assemblers <- sort(unique(id_table$assembler))

if (!is.null(input_names)){
   if (!is.null(in_colors)) {
   col_dict <- read_delim(in_colors, col_names = FALSE)
   custom_colors <- setNames(col_dict$X2, col_dict$X1)
   } else {
    if (length(sample_names) > 12) {
    palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(sample_names))
    } else {
   palette_colors <- safe
   }
   custom_colors <- setNames(palette_colors[seq_along(sample_names)], sort(sample_names))
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

## find linear breakpoints
get_breaks_lin <- function(min, max, steps = 10) {
  breaks <- pretty(c(min, max), n = steps)
  labels <- format_labels(breaks)
  
  
  break_list <- list()
  break_list[[1]] <- breaks
  break_list[[2]] <- labels
  
  return(break_list)
}


## contig N50 plots for all polymerases
## load chrom size plots
prepare_n50 <- function(chrom_file_path = "", asm = "") {
  chrom_file <- chrom_file <- read_delim(chrom_file_path, col_names = FALSE,
+                           col_types = cols(X1 = col_character(), .default = col_double()))
  
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
n50_summary <- function(path, ids, scales) {
  
  ## set up variables
  tmp_path <- paste(path, "*.fa.fai", sep="/")
  in_files <- Sys.glob(tmp_path)
  chrom_list <- list()
  alt <- paste(ids[order(nchar(ids), decreasing = TRUE)], collapse = "|")
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], alt)
    
    if (!is.na(pattern)) {
    chrom_list[[i]] <- prepare_n50(in_files[i], asm=pattern)
  } 
  }
  
  print(chrom_list)
  
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
                 
                 
  ## get break points for plot
  break_list <- scales$n50_breaks
  
  ## create N50 plot
  chrom_df <- chrom_df %>% left_join(id_table, by = c("Assembly" = "asm_id"))
  n50_graph_plot <- ggplot(chrom_df, aes(n_val, length, color = sample)) +
    geom_step(size = 0.5) +
    theme_bw() +
    theme(axis.text= element_text(size = 15),
          axis.title= element_text(size = 18)) + 
    scale_y_log10(
      breaks = break_list[[1]],
      label = break_list[[2]],
      limits = scales$n50_limits) +
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
compleasm_summary <- function(path, compleasm_label, ids, scales) {
  
  ## set up variables
  tmp_path <- paste(path, "*rf.txt", sep="/")
  in_files <- Sys.glob(tmp_path)
  compl_list <- list()
  alt <- paste(ids[order(nchar(ids), decreasing = TRUE)], collapse = "|")
  
  ## create data table
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], alt)
    if (!is.na(pattern)) {
    compl_list[[i]] <- read_compl(in_files[i], name=pattern)
  }
  }
  
  
  compl_list <- compl_list[lengths(compl_list) > 0]
  
  #compl_list <- compl_list[lengths(compl_list) > 0]
  compl_all <- rbindlist(compl_list)
  compl_all$category <- factor(compl_all$BUSCO_class, levels = c("M","I","F","D","S"), ordered = TRUE)
  
  ## color scheme
  custom_colors_s <- c("S" = "#274655",
                       "D" = "#2A9D8E",
                       "F" = "#E9C46B",
                       "I" = "#E9C48D",
                       "M" = "#E56F52")
  
  sel_order_compl <- 
    compl_all%>% 
    group_by(polymerase) %>%
    mutate(
      "SD" = sum(X3[category %in% c("S", "D")])
    ) %>%
    arrange(SD) %>% 
    filter(category == "S") %>% 
    mutate(polymerase = factor(polymerase))
  
  #sel_order_compl <- 
  #  compl_all%>% 
  #  filter(category == "S") %>% 
  #  arrange(X3) %>% 
  #  mutate(polymerase = factor(polymerase))
  
  compl_all_mut <- compl_all%>% 
    mutate(polymerase = factor(polymerase, levels = sel_order_compl$polymerase, ordered = TRUE))
  
  compl_all_mut <- compl_all_mut %>%
    left_join(id_table, by = c("polymerase" = "asm_id"))
  
  compleasm_plot <- ggplot(compl_all_mut, aes(X3, sample, fill = category)) +
    geom_col() +
    theme_bw() +
    theme(axis.text= element_text(size = 12),
          axis.title= element_text(size = 18))+
    scale_fill_manual(name = "compleasm status", values = custom_colors_s, 
                      labels =  c("missing", "fragmented(split)", "fragmented(partial)", "complete (duplicated)", "complete (single)")) +
    xlab(compleasm_label)
  
  ## Compleasm Summary
  compl_all_out_table <- compl_all_mut %>%
    group_by(polymerase) %>%
    pivot_wider(names_from=BUSCO_class,values_from = c(X3,X2))%>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
    mutate(sum_intact = X3_S + X3_D,
           perc_intact = X2_S + X2_D,
           sum_all = X3_S + X3_D + X3_F + X3_I + X3_M) %>%
    select(-category)
  
  compleasm_out_list <- list()
  
  compleasm_out_list[[1]] <- compl_all_out_table
  compleasm_out_list[[2]] <- compleasm_plot
  
  ## return outlist
  return(compleasm_out_list)
  
}

## Function to read seqkit stats file
read_seqkit <- function(path, file_pattern)
{
seqkit <- read_delim(path, col_names = T)
seqkit_renamed <- seqkit %>%
  mutate(polymerase = str_split_i(str_replace(file, file_pattern, ""), "/", 2))
  
print(seqkit_renamed)  
seqkit_breaks <- get_breaks_lin(0, max(seqkit_renamed$sum_len)*1.1)

seqkit_sequenced_gigas <- ggplot(seqkit_renamed, aes(polymerase, sum_len)) +
  geom_col(aes(fill = polymerase)) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(labels = seqkit_breaks[[2]], 
  limits = c(0, max(seqkit_renamed$sum_len) + 10000),
  breaks = seqkit_breaks[[1]]) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 15),
        axis.title= element_text(size = 18))+
  ylab("Sequenced nts")

return(seqkit_sequenced_gigas)

}

output_hifieval_readstats <- function(path, ids, scales)
{
  ## set up variables
  tmp_path <- paste(path, "*.summary.tsv", sep="/")
  in_files <- Sys.glob(tmp_path)
  hifieval_list <- list()
  alt <- paste(ids[order(nchar(ids), decreasing = TRUE)], collapse = "|")
  
  read_hifieval <- function(path, name) {
    hifieval_reads <- fread(path, header = T)
    hifieval_reads <- hifieval_reads %>%
      mutate(polymerase = name)
    return(hifieval_reads)
  }
  
  ## read files
  for (i in seq_along(in_files)) {
    pattern <- str_extract(in_files[i], alt)
    if (!is.na(pattern)){
    hifieval_list[[i]] <- read_hifieval(in_files[i], name=pattern)
    }
  }
  
  ## convert to dataframe
  hifieval_list <- hifieval_list[lengths(hifieval_list) > 0]
  hifieval_df <- rbindlist(hifieval_list)
  print(hifieval_df)
  str(hifieval_df)
  
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
  
  ## plot
  ## Only includes corrected bases (as proxy for error rate in reads)
  hifieval_plot <- ggplot(hifieval_df_sum %>%
                              filter(correction_class %in% c("corrected_bases")),
                              aes(polymerase, fraction, fill = polymerase)) +
    #geom_col(position = "dodge2") +
    geom_col() +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(axis.text.x = element_blank()) + 
    ylab("% of corr. bases")
  
  hifieval_out <- list()
  
  ## Write output
  hifieval_out[[1]] <- hifieval_plot 
  hifieval_out[[2]] <- hifieval_df_sum
  
  return(hifieval_out)
}


merqury_asm_sum <- function(path, ids, scales)
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
  alt <- paste(ids[order(nchar(ids), decreasing = TRUE)], collapse = "|")
  
  ## create data table
  for (i in seq_along(in_files_qv)) {
    pattern <- str_extract(in_files_qv[i], alt)
    if (!is.na(pattern)) {
    merqury_list_qv[[i]] <- read_merqury_oneline(in_files_qv[i], name=pattern)
  }
  }
  
  for (i in seq_along(in_files_com)) {
    pattern <- str_extract(in_files_com[i], alt)
    if (!is.na(pattern)) {
    merqury_list_com[[i]] <- read_merqury_oneline(in_files_com[i], name=pattern)
  }
  }
  
  merqury_list_qv <- merqury_list_qv[lengths(merqury_list_qv) > 0]
  mq_qv <- rbindlist(merqury_list_qv)
  mq_qv <- mq_qv %>%
    rename(error_rate = X5, qv = X4)
  
  merqury_list_com <- merqury_list_com[lengths(merqury_list_com) > 0]
  mq_com <- rbindlist(merqury_list_com)
  mq_com <- mq_com %>%
    rename(completeness = X5)

  merqury_all <- mq_qv %>%
    left_join(mq_com, by = c("X1" = "X1")) %>%
    rename(polymerase = X1)
  
  print(merqury_all)
    
  er <- ggplot(merqury_all, aes(polymerase, error_rate)) +
  geom_col(aes(fill = polymerase)) +
  #ggtitle("Merqury assembly error rate") +
  ylab("error rate") + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
    theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 10),
        axis.title= element_text(size = 11))

  qv <- ggplot(merqury_all, aes(polymerase, qv)) +
  geom_col(aes(fill = polymerase)) +
  coord_cartesian(ylim = c(min(merqury_all$qv)-5, 60)) +
  #ggtitle("Merqury assembly quality value") +
  ylab("QV") + 
  theme_bw() +
  scale_fill_manual(values = custom_colors) +
    theme(axis.text.x= element_blank(),
        axis.text.y= element_text(size = 10),
        axis.title= element_text(size = 11))

  comp <- ggplot(merqury_all, aes(polymerase, completeness)) +
  geom_col(aes(fill = polymerase)) +
  coord_cartesian(ylim = c(min(merqury_all$completeness)-2, 100)) +
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

## seqkit plot
seqkit_giga <- read_seqkit(input_path_seqkit_stats, ".fastq.gz")

all_len <- map_dfr(Sys.glob(paste(input_path_faidx, "*.fa.fai", sep = "/")),
                   ~ read_delim(.x, col_names = FALSE) %>% select(len = X2))
scales <- list(
  n50_breaks = get_breaks(10000, 1.5 * max(all_len$len)),
  n50_limits = c(1000, max(all_len$len))
)

## One page per assembler, same layout on each
build_page <- function(asm) {
  ids  <- id_table$asm_id[id_table$assembler == asm]
  n50  <- n50_summary(input_path_faidx, ids, scales)
  comp <- compleasm_summary(input_path_compleasm, compleasm_database_name, ids, scales)
  merq <- merqury_asm_sum(input_path_merqury, ids, scales)
  ttl  <- if (length(assemblers) > 1) paste("Assembler:", asm) else NULL

  if (!is.null(input_path_hifieval) && asm != "flye") {
    hife <- output_hifieval_readstats(input_path_hifieval, ids, scales)
    page <- seqkit_giga / n50[[2]] / comp[[2]] / merq[[1]] / hife[[1]] +
            plot_layout(guides = "collect", heights = c(2,3,2,2,2))
  } else {
    page <- seqkit_giga / n50[[2]] / comp[[2]] / merq[[3]] +
            plot_layout(guides = "collect", heights = c(2,3,2,2))
  }
  page + plot_annotation(tag_levels = "A", title = ttl)
}


## Final Plot
pdf(output_final_figure, width = 12, height = 14)
for (asm in assemblers) print(build_page(asm))
dev.off()

## tables stay single files, with an assembler column
per_asm <- function(f) map_dfr(assemblers, function(asm) {
  ids <- id_table$asm_id[id_table$assembler == asm]
  f(ids) %>% mutate(assembler = asm)
})

write_delim(per_asm(function(i) n50_summary(input_path_faidx, i, scales)[[1]]),
            output_ng_table, delim = "\t")
write_delim(per_asm(function(i) compleasm_summary(input_path_compleasm,
                                                  compleasm_database_name, i, scales)[[1]]),
            output_compleasm_table, delim = "\t")
write_delim(per_asm(function(i) merqury_asm_sum(input_path_merqury, i, scales)[[2]]),
            out_merqury_table, delim = "\t")
if (!is.null(input_path_hifieval)) {
  hife_asm <- assemblers[assemblers != "flye"]
  write_delim(map_dfr(hife_asm, function(asm) {
    ids <- id_table$asm_id[id_table$assembler == asm]
    output_hifieval_readstats(input_path_hifieval, ids, scales)[[2]] %>%
      mutate(assembler = asm)
  }), out_hifieval_table, delim = "\t")
}
