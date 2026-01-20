## Script to calculate contig breaks from paf alignments against a reference
library(GenomeInfoDb)
library(tidyverse)
library(patchwork)
library(data.table)
library(GenomicRanges)
library(RColorBrewer)
library(S4Vectors)

## Safe (for colorblind people) Color Palette
safe <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

## Input
input_path_paf = snakemake@params[["paf_path"]]
input_names = snakemake@params[["sample_names"]]
in_colors = snakemake@params[["colors"]]

## output
output_bed = snakemake@output[["bed"]]
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


## Custom sort function in order to retain the order that we use in all other plots
custom_sort <- function(x, order_vec) {
  
  ## Not necessary in cases where there is only one string per sample, but included for 
  ## strings separated by ,
  n_parts <- lengths(regmatches(x, gregexpr(",", x))) + 1
  sort_key <- n_parts
  
  ## Put Ranks to letters (see what happens with PolA and PolB)
  first_letters <- sub(",.*", "", x)
  first_rank <- match(first_letters, order_vec)
  
  # Second letter for strings with ,
  second_letters <- ifelse(grepl(",", x), sub(".*,", "", x), NA_character_)
  second_rank <- match(second_letters, order_vec)
  
  ## lexicographic tiebreaker
  order(sort_key, first_rank, second_rank, x)
}


contig_break_bed <- function(paf_file,
                            asm = "",
                            min_contig_size = 10000,
                            min_map_q = 50,
                            min_aln_len = 1000,
                            lookup = 10
)
{
  
  paf.table <- readPaf(
    paf.file = paf_file, include.paf.tags = TRUE,
    restrict.paf.tags = "cg"
  )
  
  min_after_leads <- function(col, n, default = Inf) {
    leads <- map(1:n, ~ lead(col, .x, default = default))
    purrr::reduce(leads, pmin, na.rm = TRUE)
  }
  
  max_before_lags <- function(col, n, default = -Inf) {
    lags <- map(1:n, ~ lag(col, .x, default = default))
    purrr::reduce(lags, pmax, na.rm = TRUE)
  }
  
  
  contig_breaks <- paf.table %>%
    group_by(q.name) %>%
    mutate(covered_nts = sum(aln.len)) %>%
    filter(covered_nts >= min_contig_size, mapq >= min_map_q, aln.len >= min_aln_len) %>%
    ungroup() %>%
    group_by(q.name, t.name) %>%
    arrange(t.start, t.end) %>%
    summarise(start = min(t.start), 
              end = max(t.end),
              aln.len = sum(aln.len)) %>%
    filter(end - start < 10 * aln.len) %>%
    ungroup() %>%
    select(-q.name) %>%     
    group_by(t.name) %>%
    arrange(start) %>%
    mutate(
      max_before = max_before_lags(end, lookup), 
      min_after = min_after_leads(start, lookup),
      to_filter_tmp = max_before > end,
      to_filter = to_filter_tmp | 
        (lag(!to_filter_tmp, 1, default = F) & 
           lead(!to_filter_tmp, 1, default = F) & 
           (start < max_before & end > min_after) & 
           (min_after - max_before < 5000))
    ) %>%
    filter(!to_filter) %>%
    mutate(
      break_start = end - 2000,
      break_end = end + 2000
    ) %>%
    select(-to_filter_tmp, -to_filter, -max_before, -min_after) %>%
    filter(n() > 1) %>%
    ungroup()
  
  
  contig_bed <- contig_breaks %>%
    select(-start, -end) %>% #relocate(contig) -q.name
    mutate(polymerase = asm) %>%
    dplyr::rename(chr = t.name,
                  start = break_start,
                  end = break_end)
  
  return(contig_bed)
}


merge_bed <- function(path) {
  ## set up variables
  tmp_path <- paste(path, "*_to_ref.paf", sep="/") #if window size is a variable: paste(path, paste0("*.", window_size, ".win.stat.gz"), sep="/")
  in_files <- Sys.glob(tmp_path)
  paf_list <- list()
  
  
  ## read files (sort was added)
  for (i in seq_along(sort(in_files))) {
    pattern <- str_extract(in_files[i], paste(input_names, collapse="|"))
    print(pattern)
    paf_list[[i]] <- contig_break_bed(in_files[i], asm=pattern)
    names(paf_list)[i] <- pattern
  }
  
  ## Merge all bed files into one dataframe
  all_breaks <- rbindlist(paf_list)
  

  all_long <- all_breaks[start != 0][
    , .(end = min(end)), by = .(chr, start, polymerase)
  ][
    end - start <= 20000
  ]
  
  ## Create GRanges obejct from GenomicRanges package  
  all_bed <- with(all_long, GRanges(
    seqnames = chr,
    ranges = IRanges(start, end),
    polymerase = polymerase
  ))
  
  ## Merge as in bedtools merge
  merged_gr <- GenomicRanges::reduce(all_bed, with.revmap=TRUE)

  ## Add polymerase info back
  merged_gr$polymerase <- sapply(merged_gr$revmap, function(idx) {
    paste(sort(unique(mcols(all_bed[idx])$polymerase)), collapse=",")
  })
  
  ## Change into dataframe
  all_merged <- as.data.frame(merged_gr)[, c("seqnames", "start", "end","polymerase")]
  colnames(all_merged) <- c("chr", "start", "end", "polymerase")
  
  ###############################
  ##      output plot          ##
  ###############################
  
  ## Count dropout events:
  prepare_mat <- all_merged %>%
    group_by(polymerase) %>%
    summarise(dropouts = n())
  
  ## Get samples:
  sample_names <- prepare_mat %>%
    mutate(polymerase = str_split(polymerase, ',')) %>%
    unnest(polymerase)
  
  ## Create a matrix that will be filled
  samples <- unique(sample_names$polymerase)[custom_sort(unique(sample_names$polymerase), labels)]
  mat <- matrix(0, nrow = length(samples), ncol = length(samples))
  rownames(mat) <- colnames(mat) <- samples
  
  # Fill matrix from dropout counts in perpare_mat
  for(i in 1:nrow(prepare_mat)) {
    samps <- strsplit(prepare_mat$polymerase[i], ",")[[1]]
    principal_samp <- samps[1]
    for(s in samps) {
      if (length(samps) > 1 & s == principal_samp) {
        next
      }
      mat[principal_samp, s] <- mat[principal_samp, s] + prepare_mat$dropouts[i]
      current_index <- match(s,samps)
      if (length(samps) > 1 & current_index != length(samps)) {
        mat[s, samps[current_index + 1]] <- mat[s, samps[current_index + 1]] + prepare_mat$dropouts[i]
      }
    }
  }
  
  ## Full Matrix for barplots
  full_mat <- mat
  
  for(i in 1:nrow(full_mat)) {
    for(j in 1:ncol(full_mat)) {
      if(i != j) {
        full_mat[i,j] <- full_mat[j,i] <- max(full_mat[i,j], full_mat[j,i])
      }
    }
  }
  
  total_dropouts <- rowSums(full_mat)
  total_dropouts <- tibble(total_dropouts, polymerase = names(total_dropouts))
  
  ## Dataframe for tileplot
  tile_plot_df <- expand.grid(Row = samples, Col = samples) %>%
    mutate(value = as.vector(mat))
  
  
  tileplot <- (ggplot(total_dropouts, aes(polymerase, total_dropouts)) +
      geom_col(aes(fill = polymerase)) +
       theme_minimal() +
       theme(axis.text.x = element_blank(),
             axis.title = element_blank(),
             legend.position = "none") +
      scale_fill_manual(values = custom_colors)) /
    (ggplot(tile_plot_df, aes(x = Col, y = Row, fill = value)) +
       geom_tile(color = "white") +
       scale_fill_gradient2(low = "lightblue",
                            mid = "pink",
                            high = "red", 
                            transform = "log10", 
                            midpoint = 100) +
       geom_text(aes(label = ifelse(value!=0, value,"")), color = "black", size = 4) + 
       theme_minimal() +
       theme(axis.text.x = element_text(angle = 45, hjust = 1),
             axis.title = element_blank()) +
       labs(fill = "Contig Breaks")) +
    plot_layout(heights = c(0.4, 1))
  
  
  ## Output to return
  out_list <- list()
  out_list[[1]] <- all_merged
  out_list[[2]] <- tileplot
  
  return(out_list)
  
}

## Create output
## Bed file and Confusion Matrix of breaks
contig_break_out <- merge_bed(input_path_paf)

## Write Output
write_delim(contig_break_out[[1]], output_bed, delim = "\t", col_names = TRUE)
ggsave(output_plot, plot = contig_break_out[[2]], device = "pdf", width = 12, height = 14, units = "in")


