#' ---
#' title: "Seidr assessment - KEGG plots"
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
  
  library(org.Dm.eg.db)
  library(org.At.tair.db)
  library(org.Sc.sgd.db)
  
  library(SeidRFile)
  library(KEGGREST)
  library(pheatmap)
})

#' * Out dir
dir.create(here("analysis/seidr-plots"),recursive=TRUE,showWarnings=FALSE)

#' * Functions
make_edge_tibble <- function(x, indices, algorithms, nodes, go, id) {
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
    left_join(go, by = c("Source" = id)) %>% 
    left_join(go, by = c("Target" = id))
}

extract_edges <- function(sf,db,key){
  go <- as_tibble(select(db, keys = nodes(sf),c(key, "PATH"), key))
  indices <- unique(match(go[[key]], nodes(sf)))
  
  return(bind_rows(seidr_chunked_apply(
    sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
    FUN = make_edge_tibble, 
    indices = indices, 
    algorithms = algorithms(sf),
    nodes = nodes(sf),
    go = go,
    id=key)))
}

plot_summary <- function(edges,filename,
                         lvls=c("Seidr", "CLR", "ELNET", "GENIE3", "LLR", 
                                "MI", "NARROMI", "PCOR","Pearson", "PLSNET", 
                                "Spearman", "TIGRESS", "TOMsimilarity", "TOM"),
                         thm=theme(text = element_text(size = 20),
                                    axis.text.x = element_text(angle = 45, 
                                                               hjust = 1, 
                                                               vjust=0.8),
                                    legend.position = 'none')
                         ){
  
  edges %>%
    # Check if both nodes are in the same PW
    mutate(InPath = PATH.x == PATH.y) %>% 
    # Drop nodes without any PW annotation
    filter(! is.na(InPath)) %>%
    # Consider edges with annotations in two separate pathways for both by
    # repeating the edge
    pivot_longer(cols = c(PATH.x, PATH.y), names_to = "PathSource",
                 values_to = "ReprPath") %>%
    # Pivot algorithms
    pivot_longer(-c(Source, Target, MatchType, InPath, PathSource, ReprPath),
                 names_to = "Algorithm", values_to = "Rank") %>%
    # Get rank means per pathway
    group_by(ReprPath, InPath, Algorithm) %>%
    filter(Algorithm != "ARACNE") %>%
    summarise(M = mean(Rank, na.rm = TRUE)) %>%
    mutate(Algorithm = recode(Algorithm, irp = "Seidr",
                              SPEARMAN = "Spearman",
                              PEARSON = "Pearson",
                              TOMSimilarity = "TOMsimilarity")) %>%
    mutate(Algorithm = factor(Algorithm, levels = lvls)) %>%
    ggplot(aes(y = M, fill = InPath, x = Algorithm)) +
    geom_boxplot(notch = TRUE) +
    stat_compare_means(label = "p.signif",
                       method.args = list(alternative = "greater")) +
    theme_bw() +
    scale_fill_npg() +
    scale_y_continuous(labels = comma) +
    theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
    thm +
    xlab("") +
    ylab("Mean Edge Weight Rank per Pathway")
  
  ggsave(filename, width = 16, height = 9)
}

run_org <- function(plt,rda,sf,org.db,key){
  if (!file.exists(plt)){
    if (file.exists(rda)){
      load(rda)
    } else {
      edges <- extract_edges(SeidrFileFromPath(sf),org.db,key)
      gc()
      save(edges,file=rda)
    }
    dev.null <- plot_summary(edges,plt)
    rm(edges)
    gc()
  }
}

#' # Data
#' ## S. cerevisiae
#' * Load or prep and plot the data
dev.null <- run_org(plt=here("analysis/seidr-plots/sc_kegg_edgerank.pdf"),
                    rda=here("networks/scere/edges.rda"),
                    sf=here("networks/scere/aggregate/irp.sf"),
                    org.db=org.Sc.sgd.db,
                    key="ENSEMBL")

#' ## A. thaliana
#' * Load the data
dev.null <- run_org(plt=here("analysis/seidr-plots/at_kegg_edgerank.pdf"),
                    rda=here("networks/athal/edges.rda"),
                    sf=here("networks/athal/results/aggregate/irp.sf"),
                    org.db=org.At.tair.db,
                    key="TAIR")

#' ## D. melanogaster
#' * Load the data
dev.null <- run_org(plt=here("analysis/seidr-plots/dm_kegg_edgerank.pdf"),
                    rda=here("networks/dmel/edges.rda"),
                    sf=here("networks/dmel/results/aggregate/irp.sf"),
                    org.db=org.Dm.eg.db,
                    key="FLYBASE")

#' Session Info
#' ```{R session info, echo=FALSE}
#' sessionInfo()
#' ```
