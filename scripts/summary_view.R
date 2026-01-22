## Load libraries from environment
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(colorspace)
library(ggpubr)
library(S4Vectors)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input
input_path_paf = snakemake@params[["paf_path"]]
input_path_pandepth = snakemake@params[["pandepth_path"]]
input_bed_path = snakemake@input[["breakpoint_bed"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output
output_summary_plot = snakemake@output[["summary_plot"]]

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

## Read pandepth files and long join them
## Since SVbyEye cannot be installed form bioconductor or cran, I copied the readPaf function from it here (will cite it in the reference manual)
## Needed:
startTimedMessage <- function(...) {
  
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  ptm <- proc.time()
  return(ptm)
  
}


stopTimedMessage <- function(ptm) {
  
  time <- proc.time() - ptm
  message(" ", round(time[3],2), "s")
  
}

#' Read PAF from an input file
#'
#' This function takes an PAF output file from minimap2 and loads the file along
#' with user defined set of additional alignment tags (see PAF specification).
#'
#' @param paf.file A path to a PAF file containing alignments to be loaded.
#' @param include.paf.tags Set to \code{TRUE} if all additional PAF alignment tags should be included in the output.
#' @param restrict.paf.tags Define a set of PAF tag ids (e.g. NM, cg) to be reported in the output.
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols
#' @importFrom S4Vectors lapply
#' @return A \code{tibble} of loaded PAF alignments
#' @author David Porubsky
#' @export
#' @examples
#' ## Get PAF to read in ##
#' paf.file <- system.file("extdata", "test1.paf", package = "SVbyEye")
#' ## Read in PAF
#' paf.table <- readPaf(paf.file = paf.file)
#' ## Read in PAF including all tags
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE)
#' ## Read in PAF including CIGAR string only
#' paf.table <- readPaf(paf.file = paf.file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
#'
readPaf <- function(paf.file = NULL, include.paf.tags = TRUE, restrict.paf.tags = c("NM", "cg")) {
  ## Check user input ##
  if (is.null(paf.file)) {
    stop("Path to a PAF file to load is not defined!!!")
  }
  
  ## Check if file exists and is TAB delimited
  if (file.exists(paf.file)) {
    con <- file(paf.file, "r")
    first.line <- readLines(con, n = 1)
    if (length(first.line) > 0) {
      n.fields <- length(stringr::str_split(first.line, "\t")[[1]])
      if (n.fields < 12) {
        stop("User defined 'paf.file' has less then 12 expected tab-delimeted fields!!!")
      }
      close(con)
    } else {
      stop("User defined 'paf.file' seems to be empty!!!")
    }
  } else {
    stop("User defined 'paf.file' doesn't exist!!!")
  }
  
  ## Load PAF file ##
  if (file.exists(paf.file)) {
    ptm <- startTimedMessage(paste0("[readPaf] Loading PAF file: ", paf.file))
    ## Read PAF lines
    paf.lines <- readLines(paf.file)
    if (include.paf.tags) {
      fields <- stringr::str_split(paf.lines, "\t")
    } else {
      fields <- stringr::str_split(paf.lines, "\t", 13)
    }
    paf.fields <- S4Vectors::lapply(fields, "[", seq_len(12))
    field.names <- c("q.name", "q.len", "q.start", "q.end", "strand", "t.name", "t.len", "t.start", "t.end", "n.match", "aln.len", "mapq")
    
    for (i in seq_along(paf.fields)) {
      attr(paf.fields[[i]], "names") <- field.names
    }
    paf <- dplyr::bind_rows(paf.fields)
    cols.num <- c(2, 3, 4, 7:12)
    paf[cols.num] <- suppressWarnings(dplyr::bind_cols(S4Vectors::lapply(paf[cols.num], as.numeric)))
    
    if (include.paf.tags) {
      if (any(lengths(fields) > 12)) {
        paf.tags <- S4Vectors::lapply(fields, function(x) paste(x[13:length(x)]))
        paf <- dplyr::bind_cols(paf, processPafTags(paf.tags = paf.tags, restrict.paf.tags = restrict.paf.tags))
      }
    }
    stopTimedMessage(ptm)
    return(paf)
  } else {
    stop(paste0("PAF file ", paf.file, " doesn't exists !!!"))
    return(NULL)
  }
}

#' Process PAF specific alignment tags.
#'
#' @param paf.tags A \code{list} of PAF specific tag extracted from each alignment.
#' @inheritParams readPaf
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows
#' @return A \code{tibble} object
#' @author David Porubsky
#'
processPafTags <- function(paf.tags, restrict.paf.tags = c("NM", "cg")) {
  ## Make sure restrict.paf.tags takes only expected values
  allowed.tags <- c("tp", "cm", "s1", "s2", "NM", "MD", "AS", "SA", "ms", "nn", "ts", "cg", "cs", "dv", "de", "rl")
  restrict.paf.tags <- restrict.paf.tags[restrict.paf.tags %in% allowed.tags]
  if (length(restrict.paf.tags) == 0) {
    message(paste0("Submitted 'restrict.paf.tags' are not present in the allowed set of PAF tags: ", paste(allowed.tags, collapse = "; ")))
    restrict.paf.tags <- NULL
  }
  
  tags <- character()
  to.numeric <- integer()
  res <- list()
  n <- length(paf.tags)
  t.idx <- 0
  for (ali.idx in seq_along(paf.tags)) {
    split.tags <- stringr::str_split(paf.tags[[ali.idx]], ":")
    for (tag in split.tags) {
      if (!is.null(restrict.paf.tags) & tag[1] %in% restrict.paf.tags) {
        if (!(tag[1] %in% tags)) {
          t.idx <- t.idx + 1
          if (tag[2] %in% c("f", "H", "i")) {
            to.numeric <- c(to.numeric, t.idx)
          }
          res[[tag[1]]] <- rep(NA, n)
          tags <- c(tags, tag[1])
        }
        res[[tag[1]]][ali.idx] <- tag[3]
      }
    }
  }
  
  for (i in to.numeric) {
    res[[i]] <- as.numeric(res[[i]])
  }
  dplyr::bind_rows(res)
}

## function to read bed file:
read_bed <- function(path) {

bed <- read_delim(path, col_names = T)

return(bed)

}



## function to find windows to plot for summary
plot_coverage_windows <- function(pandepth_frame, All_blocks, start, end, scaffold, max_cov) 
{
  region_blocks <- All_blocks %>%
      filter(t.name == scaffold, ref_end >= start, ref_end <= end)
      
  if (end - start < 10000) {
  end <- end + 5000
  start <- start - 5000
  }

  center_plot <- (start + end)/2
  
  out_plt <- ggplot(pandepth_frame%>%
                      filter(Start >= start, 
                             End <= end, 
                             `#Chr` == scaffold),
                    aes(Start,
                        MeanDepth)) +
    geom_line(aes(color = asm), alpha = 0.5) +
    ylim(c(0, max_cov)) + 
    theme_bw() +
    geom_vline(data = region_blocks,
                aes(xintercept = ref_end, color = assembly),
                linetype = "dashed",
                alpha = 0.4,
                linewidth = 0.5) +
    scale_color_manual(values = custom_colors) +
    scale_x_continuous(labels = comma, breaks = center_plot, name = "Genomic Position") +
    coord_cartesian(xlim = c(start,end), expand = F)
  
  return(out_plt)
}

## Function to read pandepth infiles
read_pandepth <- function(pandepth_file, asm = "", yoff) {
  
  read_last <- function(file) {
    tail(readLines(file), n=1)
  }
  
  mead_depth <- sub(".*MeanDepth: *([0-9.]+).*", "\\1", read_last(pandepth_file))
  
  in_pandepth <- fread(pandepth_file, header = T)
  out_frame <- in_pandepth %>%
    mutate(asm = asm,
           NormDepth = MeanDepth / as.numeric(mead_depth),
           roundGC = round((as.numeric(`GC(%)`)/100), 2),
           yoff = yoff)
  
  return(out_frame)
}

contig_break_blocks <- function(paf_file,
                                asm = "",
                                y_number,
                                min_contig_size = 10000,
                                min_map_q = 50,
                                min_aln_len = 1000,
                                lookup = 10)
{
  
  paf.table <- readPaf(
    paf.file = paf_file, include.paf.tags = TRUE,
    restrict.paf.tags = "cg"
  )
  
  paf.filtered <- paf.table %>%
    group_by(q.name) %>%
    mutate(covered_nts = sum(aln.len)) %>%
    ## filter must be changed programmatically, as different assemblies will have different specs here
    filter(covered_nts >= min_contig_size, mapq >= min_map_q, aln.len >= min_aln_len) %>%
    group_by(q.name, t.name) %>%
    arrange(t.start, t.end) %>%
    summarise(ref_start = min(t.start),
              ref_end = max(t.end),
              aln.len = sum(aln.len)) %>%
    mutate(assembly = asm,
           yoff = y_number) %>%
    filter(ref_end - ref_start < 10 * aln.len) %>%
    ungroup() 
  
  min_after_leads <- function(col, n, default = Inf) {
    leads <- map(1:n, ~ lead(col, .x, default = default))
    reduce(leads, pmin, na.rm = TRUE)
  }
  
  max_before_lags <- function(col, n, default = -Inf) {
    lags <- map(1:n, ~ lag(col, .x, default = default))
    reduce(lags, pmax, na.rm = TRUE)
  }
  
  
  remove_shrapnel <- paf.filtered %>%
    group_by(t.name) %>% 
    arrange(ref_start) %>%
    mutate(
      max_before = max_before_lags(ref_end, lookup), 
      min_after = min_after_leads(ref_start, lookup),
      to_filter_tmp = max_before > ref_end,
      to_filter = to_filter_tmp | 
        (lag(!to_filter_tmp, 1, default = F) & 
           lead(!to_filter_tmp, 1, default = F) & 
           (ref_start < max_before & ref_end > min_after) & 
           (min_after - max_before < 5000)),
      surrounded = (lag(!to_filter_tmp, 1, default = F) & 
                      lead(!to_filter_tmp, 1, default = F))
    ) %>%
    filter(!to_filter)
  
  return(remove_shrapnel)
}




read_data <- function(input_path_paf, input_path_pandepth, bed_path) { 
      ## set up variables
      tmp_path <- paste(input_path_pandepth, "*_all.win.stat.gz", sep="/")
      in_files <- Sys.glob(tmp_path)
      pandepth_list <- list()
      
      print(in_files)
      ## read files
      for (i in seq_along(in_files)) {
        pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
        #lib_color = custom_colors[pattern][[1]]
        print(in_files[i])
        pandepth_list[[i]] <- read_pandepth(in_files[i], asm = pattern, yoff = i)
      }
      
      pandepth_list_long <- rbindlist(pandepth_list)
      
      ## set up variables
      tmp_path <- paste(input_path_paf, "*_to_ref.paf", sep="/")
      in_files <- Sys.glob(tmp_path)
      paf_list <- list()
      
      ## read files
      for (i in seq_along(in_files)) {
        pattern <- str_extract(in_files[i], paste(labels, collapse="|"))
        paf_list[[i]] <- contig_break_blocks(in_files[i], asm = pattern, y_number = i)
      }
      
      All_blocks <- rbindlist(paf_list)

      ## find interesting windows:
      ## Function to find biggest chromosome and select a good window
      print(All_blocks)
      
      biggest_scaffold <- All_blocks %>%
        group_by(t.name) %>%
        summarise(length = max(ref_end)) %>%
        ungroup() %>%
        filter(length == max(length))

      print(biggest_scaffold)
      ## This is the biggest scaffold
      biggest_scaffold_name <- unique(biggest_scaffold$t.name)

      ## if context window is large but assemblies are fragmented
      context_fraction <- All_blocks %>%
        filter(t.name == biggest_scaffold_name) %>%
        group_by(assembly) %>%
        summarise(mapped_contigs =n ()) %>%
        ungroup()

      ## Median of contigs mapped to largest reference contig:
      context_fraction_num <- mean(context_fraction$mapped_contigs)
      
      ## Context window:
      context_window <- min(biggest_scaffold$length/500, biggest_scaffold$length/(2*context_fraction_num))
      
      ## Avoid too small windows. Minimum size should be 50,000
      if (context_window < 50000) {
        context_window <- 50000
      }
      
      print("-------------------------")
      print("context_window")
      print(context_window)
      print("-------------------------")
      
      
      max_drops <- length(unique(All_blocks$assembly)) * 2.5
      
      ## Find a representative region (size ~ 1 Mb with all drop outs)
      windows_with_all_samples <- All_blocks %>%
        rowwise() %>%
        filter({
          valid_rows <- which(All_blocks$ref_end <= ref_end + context_window & 
                                All_blocks$ref_end >= ref_end - context_window & 
                                All_blocks$t.name == t.name)
          window_samples <- All_blocks$assembly[valid_rows]
          length(unique(window_samples)) + 0.5 * length(unique(window_samples)) >= length(unique(All_blocks$assembly)) & 
            length(window_samples) <= max_drops #& length(window_samples) >= length(unique(All_blocks$assembly))
        })%>%
        ungroup() %>%
        mutate(
          window_start = ref_end - (context_window + 10000),
          window_end = ref_end + (context_window + 10000),
          window_id = ceiling(window_start/context_window)
        )
    
      
      ## Windows to plot
      slices <- windows_with_all_samples %>%
        filter(t.name == biggest_scaffold_name) %>%
        group_by(window_id) %>%
        summarise(start = ifelse(min(window_start) <0, 0, min(window_start)),
                  end = max(window_end)) %>%
        select(-window_id)
    
      print(slices)
      
      ## Max Y value for coverage plots:
      max_cov <- pandepth_list_long %>%
        group_by(asm) %>%
        summarise(mean_depth = mean(MeanDepth))
      
      max_cov <- max(max_cov$mean_depth) * 1.75
         
      ## summarising factor to create at least 500 data points for lines
      ## window_size = 200, data_points   = 1000 -> 200*1000 = 20,000
      window_factor <- max(10^round(log10(context_window)) / 200 / 1000, 1)
      
      
      ## Bed file
      breakpoints_bed <- read_bed(bed_path)
      
      ## Output plots
      out_plots <- pmap(list(slices$start[1:min(50, length(slices$start))], slices$end[1:min(50, length(slices$start))]), 
                       ~{
        
        plot1 <- contig_breaks_zoom(All_blocks, .x, .y, biggest_scaffold_name)
        
        coverage_lines <- pandepth_list_long %>%
          filter(`#Chr` == biggest_scaffold_name, Start >= .x, End <= .y) %>%
          arrange(Start) %>% 
          mutate(window = ceiling(row_number() / window_factor)) %>%
          group_by(window, asm) %>%
          summarise(
            pos = min(End),
            mean_depth_yoff = mean(MeanDepth/(max_cov * 1.75) + yoff)) %>%
          ungroup()
        
        print(coverage_lines)
        
        plot2 <- plot1 + theme(axis.text.y = element_blank()) + 
          geom_line(data = coverage_lines %>%
                      mutate(ref_end = pos,
                             yoff = mean_depth_yoff), 
                    aes(x = ref_end, y = mean_depth_yoff, color = asm),
                    linewidth = 0.1)
        
        overview <- ggplot(All_blocks %>%
                             filter(t.name == biggest_scaffold_name), 
                           aes(y = yoff, 
                               yend = yoff,
                               x = ref_start, 
                               xend = ref_end)) +
          geom_segment(aes(color = assembly),
                       linewidth = 0.1, 
                       arrow = arrow(90, unit(1, "mm")))  +
          scale_color_manual(values = custom_colors) +
          scale_x_continuous(labels = comma) +
          theme_void() +
          theme(
            legend.position = "none",
            axis.line.x = element_line(colour = "black"),
            axis.ticks.x = element_line(color = "black",
                                        linewidth = 0.1),
            axis.ticks.length.x = unit(2, "mm"),
            axis.text.x = element_text(size = 10)) + 
          coord_cartesian(ylim = c(0.5,max(All_blocks$yoff) + 0.5), expand = F) +
          geom_rect(aes(xmin = .x - 10000,
                        xmax = .y + 10000,
                        ymin = min(yoff) - 0.3 ,
                        ymax = max(yoff) + 0.3),
                    fill = "transparent",
                    color = "black")
        
        bed_flt <- breakpoints_bed %>%
          filter(chr == biggest_scaffold_name, start >= .x, end <= .y ) %>%
          select(start, end)
        
        
        Cov_plots <- pmap(bed_flt, ~{
          plot_coverage_windows(pandepth_list_long, 
                                All_blocks, 
                                start = .x, 
                                end = .y, 
                                scaffold=biggest_scaffold_name, 
                                max_cov = max_cov)
          
        }) 
          
        final_Cov_plots <- wrap_plots(Cov_plots, guides = "collect", axis_titles = "collect", axes = "collect")
        
        plots <- overview + plot2 + final_Cov_plots
        
        layout <- plot_layout(design = "1\n2\n3", axis_titles = "collect", 
                              heights = c(0.5,0.75,1))
        plots + layout
        
      })

  return(out_plots)
}


## Write into function
contig_breaks_zoom <- function(blocks, start, end, contig) {
  zoom <- blocks %>%
    filter(t.name == contig)
  
  ggplot(zoom, aes(y = yoff, 
                   yend = yoff,
                   x = ref_start, 
                   xend = ref_end)) +
    geom_segment(aes(color = assembly),
                 linewidth = 1, 
                 arrow = arrow(90, unit(2, "mm"))) +
    scale_y_continuous(labels = unique(zoom$assembly),
                       breaks = unique(zoom$yoff)) + 
    scale_color_manual(values = custom_colors) +
    scale_x_continuous(labels = comma) +
    theme_void() +
    theme(
          legend.position = "none",
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(color = "black",
                                      linewidth = 0.1),
          axis.ticks.length.x = unit(2, "mm"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10)
          )+
    coord_cartesian(ylim = c(0.5,max(zoom$yoff) + 0.5),
                    xlim = c(start, end), expand = F) +
    xlab(paste0("reference contig: ", contig))
}

## output plots
out_summary <- read_data(input_path_paf, input_path_pandepth, input_bed_path)

#ggsave(output_summary_plot, plot = out_summary[[1]], device = "pdf", width = 12, height = 14, units = "in")
ggexport(out_summary, filename=output_summary_plot)


