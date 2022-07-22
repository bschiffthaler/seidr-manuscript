#' ---
#' title: "Seidr assessment - Yeastract plots"
#' author: "Bastian Schiffthaler and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(scales)
  library(ggpubr)
  library(here)
  library(SeidRFile)
})

#' * Out dir
dir.create(here("analysis/seidr-plots"),recursive=TRUE,showWarnings=FALSE)

#' * Function
make_edge_tibble <- function(x, indices, algorithms, nodes, yeastract) {
  xi <- x$node_i %in% indices
  xj <- x$node_j %in% indices
  f <- xi | xj
  k <- xi[f] + xj[f]
  r <- x$rank[f]
  ranks <- do.call(rbind, r)
  colnames(ranks) <- algorithms
  tbl <- tibble(Source = nodes[x$node_i[f]],
                Target = nodes[x$node_j[f]],
                MatchType = k)
  cbind(tbl, ranks) %>% 
    left_join(yeastract, by = c("Source" = "Source", "Target" = "Target"))
}

extract_edges <- function(sf,yeastract){
  
  indices <- unique(match(unique(c(yeastract$Source, yeastract$Target)), nodes(sf)))  
  
  edges <- bind_rows(seidr_chunked_apply(
    sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
    FUN = make_edge_tibble, indices = indices, 
    algorithms = algorithms(sf), nodes = nodes(sf),
    yeastract=yeastract))
  
  edges$InYeastract <- !is.na(edges$Evidence)
  edges$Direction[is.na(edges$Direction)] <- "Not Connected"
  edges$Direction[edges$Direction == "N/A"] <- "Unknown"
  edges$Direction <- factor(edges$Direction, levels = c("Not Connected", "Positive", "Negative", "Unknown"))

  return(edges)
}

plot_summary <- function(edges,plt,
                         thm=theme(text = element_text(size = 20),
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.8),
                                   legend.position = 'none'),
                         lvls=c("Seidr", "CLR", "ELNET", "GENIE3", "LLR", "MI", "NARROMI", "PCOR",
                                "Pearson", "PLSNET", "Spearman", "TIGRESS", "TOMsimilarity", "TOM")){
  
  edges %>%
    pivot_longer(
      -c(Source, Target, MatchType, Sname, Tname,  N, Date, Evidence, Direction,
         Type, Experiment, InYeastract),
      names_to = "Algorithm") %>%
    filter(Algorithm != "ARACNE") %>%
    mutate(Algorithm = recode(Algorithm, irp = "Seidr",
                              SPEARMAN = "Spearman",
                              PEARSON = "Pearson",
                              TOMSimilarity = "TOMsimilarity")) %>%
    mutate(Algorithm = factor(Algorithm, levels = lvls)) %>%
    ggplot(aes(x = Algorithm, y = value, fill = Direction)) +
    geom_boxplot(notch = TRUE) +
    stat_compare_means(label = "p.signif",
                       method.args = list(alternative = "greater")) +
    theme_bw() +
    scale_fill_npg() +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
    thm +
    xlab("") +
    ylab("Edge Weight Rank")  +
    theme(legend.position = "right")
  
  rm(edges)
  gc()
  
  ggsave(plt, width = 16, height = 9)
}


#' * Variables
sf <- here("networks/scere/aggregate/irp.sf")
plt <- here("analysis/seidr-plots/sc_yeastract_edgerank.pdf")
rda <- here("networks/scere/edges_yeastract.rda")
yeastract <- here("gold_standard/yeast2019-with-regdata_all.tsv")

#' * Data
yeastract <- read_tsv(yeastract,
  col_names = c("Source", "Sname", "Target", "Tname", "N", 
                "Date", "Evidence", "Direction", "Type", 
                "Experiment"),show_col_types = FALSE)
yeastract <- bind_rows(yeastract, yeastract %>% rename(Source="Target",Target="Source"))

#' # Run
if (!file.exists(plt)){
  if (file.exists(rda)){
    load(rda)
  } else {
    edges <- extract_edges(SeidrFileFromPath(sf),yeastract)
    gc()
    save(edges,file=rda)
  }
  dev.null <- plot_summary(edges,plt)
  rm(edges)
  gc()
}
