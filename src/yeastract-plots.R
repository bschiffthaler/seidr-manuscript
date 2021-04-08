library(tidyverse)
library(ggsci)
library(scales)
library(ggpubr)

library(org.Dm.eg.db)
library(org.At.tair.db)
library(org.Sc.sgd.db)

library(SeidRFile)
library(KEGGREST)
library(pheatmap)

th <- theme(text = element_text(size = 20), legend.position = 'none')
th45 <- theme(text = element_text(size = 20),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = 'none')
levels <- c("Seidr", "CLR", "ELNET", "GENIE3", "LLR", "MI", "NARROMI", "PCOR",
            "PEARSON", "PLSNET", "SPEARMAN", "TIGRESS", "TOMSimilarity", "TOM")
# go_terms <- tibble(Term = c("GO:0031497", "GO:0007049", "GO:0042254",
#                             "GO:0006113", "GO:0015979", "GO:0035282"),
#                    Name = c("Chromatin Assembly", "Cell Cycle",
#                             "Ribosome Biogenesis", "Fermentation",
#                             "Photosynthesis", "Segmentation"))

yeastract <- read_tsv(
  "/mnt/picea/projects/seidr/gold_standard/yeast2019-with-regdata_all.tsv",
  col_names = c("Source", "Sname", "Target", "Tname", "N", "Date", "Evidence", "Direction", "Type", "Experiment")
)

make_edge_tibble <- function(x, indices, algorithms, nodes) {
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
  cbind(tbl, ranks)
}

sc_sf <- SeidrFileFromPath("/mnt/picea/projects/seidr/networks/scere/aggregate/irp.sf")
sc_nodes <- nodes(sc_sf)
sc_go <- unique(c(yeastract$Source, yeastract$Target))
sc_indices <- unique(match(sc_go, sc_nodes))
chunks <- seidr_chunked_apply(
  sc_sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
  FUN = make_edge_tibble, indices = sc_indices, algorithms = algorithms(sc_sf),
  nodes = nodes(sc_sf)
)

chunks <-bind_rows(chunks)

yeastract2 <- yeastract
yeastract2$Source <- yeastract$Target
yeastract2$Target <- yeastract$Source
yeastract_symmetric <- rbind(yeastract, yeastract2)

sc_edges <- left_join(
  chunks, yeastract_symmetric, by = c("Source" = "Source", "Target" = "Target")
)

rm(chunks)
gc()

sc_edges$InYeastract <- !is.na(sc_edges$Evidence)
sc_edges$Direction[is.na(sc_edges$Direction)] <- "Not Connected"
sc_edges$Direction[sc_edges$Direction == "N/A"] <- "Unknown"
sc_edges$Direction <- factor(sc_edges$Direction, levels = c("Not Connected", "Positive", "Negative", "Unknown"))


sc_edges %>%
  pivot_longer(
    -c(Source, Target, MatchType, Sname, Tname,  N, Date, Evidence, Direction,
       Type, Experiment, InYeastract),
    names_to = "Algorithm"
  ) -> sc_edges_long
  
sc_edges_long %>%
  filter(Algorithm != "ARACNE") %>%
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(x = Algorithm, y = value, fill = Direction)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
  th +
  xlab("") +
  ylab("Edge Weight Rank")  +
  theme(legend.position = "right")
  
ggsave("~/seidr-plots/sc_yeastract_edgerank.pdf", width = 16, height = 9)
rm(ls())
gc()
