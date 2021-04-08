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

biogrid <- read_tsv(
  "/mnt/picea/projects/seidr/gold_standard/athal_biogrid_all.tsv",
  skip = 1,
  col_names = c("Source", "Target", "Experiment", "Evidence")
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

sc_sf <- SeidrFileFromPath("/mnt/picea/projects/seidr/networks/athal/results/aggregate/irp.sf")
sc_nodes <- nodes(sc_sf)
sc_go <- unique(c(biogrid$Source, biogrid$Target))
sc_indices <- unique(match(sc_go, sc_nodes))
chunks <- seidr_chunked_apply(
  sc_sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
  FUN = make_edge_tibble, indices = sc_indices, algorithms = algorithms(sc_sf),
  nodes = nodes(sc_sf)
)

chunks <-bind_rows(chunks)

biogrid2 <- biogrid
biogrid2$Source <- biogrid$Target
biogrid2$Target <- biogrid$Source
biogrid_symmetric <- rbind(biogrid, biogrid2)

sc_edges <- left_join(
  chunks, biogrid_symmetric, by = c("Source" = "Source", "Target" = "Target")
)

rm(chunks)
gc()

sc_edges$InBiogrid <- !is.na(sc_edges$Evidence)

sc_edges$Evidence[is.na(sc_edges$Evidence)] <- "Not Connected"
sc_edges$Evidence <- recode(sc_edges$Evidence, genetic = "Genetic", physical = "Physical")
sc_edges$Evidence <- factor(sc_edges$Evidence, levels = c("Not Connected", "Genetic", "Physical"))


sc_edges %>%
  group_by(InBiogrid) %>%
  sample_n(sum(sc_edges$InBiogrid)) %>%
  pivot_longer(
    -c(Source, Target, MatchType, Evidence, Experiment, InBiogrid),
    names_to = "Algorithm"
  ) -> sc_edges_long

sc_edges_long %>%
  filter(Algorithm != "ARACNE") %>%
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(x = Algorithm, y = value, fill = Evidence)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  th +
  xlab("") +
  ylab("Edge Weight Rank")  +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

ggsave("~/seidr-plots/at_biogrid_edgerank.pdf", width = 16, height = 9)
