#' ---
#' title: "Seidr assessment - BioGRID plots"
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
  library(here)  
  library(tidyverse)
  library(ggsci)
  library(scales)
  library(ggpubr)
  library(SeidRFile)
})

#' * Out dir
dir.create(here("analysis/seidr-plots"),recursive=TRUE,showWarnings=FALSE)

#' * Function
make_edge_tibble <- function(x, indices, algorithms, nodes, biogrid) {
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
    left_join(biogrid,by = c("Source" = "Source", "Target" = "Target"))
}

extract_edges <- function(sf,biogrid){
  biogrid <- read_tsv(
    biogrid,
    skip = 1, show_col_types = FALSE,
    col_names = c("Source", "Target", "Experiment", "Evidence")
  )
  biogrid <- bind_rows(biogrid, biogrid %>% rename(Source="Target",Target="Source"))
  
  
  indices <- unique(match(unique(c(biogrid$Source, biogrid$Target)), nodes(sf)))  
  edges <- bind_rows(seidr_chunked_apply(
    sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
    FUN = make_edge_tibble, indices = indices, 
    algorithms = algorithms(sf), nodes = nodes(sf),
    biogrid=biogrid))

  edges$InBiogrid <- !is.na(edges$Evidence)
  edges$Evidence[is.na(edges$Evidence)] <- "Not Connected"
  edges$Evidence <- recode(edges$Evidence, genetic = "Genetic", physical = "Physical")
  edges$Evidence <- factor(edges$Evidence, levels = c("Not Connected", "Genetic", "Physical"))

  return(edges)
}

plot_summary <- function(edges,plt,
             thm=theme(text = element_text(size = 20),
                       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.8),
                       legend.position = 'none'),
             lvls=c("Seidr", "CLR", "ELNET", "GENIE3", "LLR", "MI", "NARROMI", "PCOR",
                    "Pearson", "PLSNET", "Spearman", "TIGRESS", "TOMsimilarity")){
  edges %>%
    group_by(InBiogrid) %>%
    sample_n(sum(edges$InBiogrid)) %>%
    pivot_longer(
      -c(Source, Target, MatchType, Evidence, Experiment, InBiogrid),
      names_to = "Algorithm") %>%
    filter(Algorithm != "ARACNE") %>%
    mutate(Algorithm = recode(Algorithm, irp = "Seidr",
                              SPEARMAN = "Spearman",
                              PEARSON = "Pearson",
                              TOMSimilarity = "TOMsimilarity",
                              TOM="TOMsimilarity")) %>%
    mutate(Algorithm = factor(Algorithm, levels = lvls)) %>%
    ggplot(aes(x = Algorithm, y = value, fill = Evidence)) +
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

run_org <- function(plt,rda,sf,biogrid){
  if (!file.exists(plt)){
    if (file.exists(rda)){
      load(rda)
    } else {
      edges <- extract_edges(SeidrFileFromPath(sf),biogrid)
      gc()
      save(edges,file=rda)
    }
    dev.null <- plot_summary(edges,plt)
    rm(edges)
    gc()
  }
}

#' # Run
#' ## A. thaliana
dev.null <- run_org(here("analysis/seidr-plots/at_biogrid_edgerank.pdf"),
                    here("networks/athal/edges_biogrid.rda"),      
                    here("networks/athal/results/aggregate/irp.sf"),
                    here("gold_standard/athal_biogrid_all.tsv"))

#' ## D. melanogaster
dev.null <- run_org(here("analysis/seidr-plots/dm_biogrid_edgerank.pdf"),
                    here("networks/dmel/edges_biogrid.rda"),      
                    here("networks/dmel/results/aggregate/irp.sf"),
                    here("gold_standard/dmel_biogrid_all.tsv"))

#' ## S. cerevisiae
dev.null <- run_org(here("analysis/seidr-plots/sc_biogrid_edgerank.pdf"),
                    here("networks/scere/edges_biogrid.rda"),
                    here("networks/scere/aggregate/irp.sf"),
                    here("gold_standard/scere_biogrid_all.tsv"))

#' Session Info
#' ```{R session info, echo=FALSE}
#' sessionInfo()
#' ```
