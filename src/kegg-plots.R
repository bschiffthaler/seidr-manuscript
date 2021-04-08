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
sc_go <- as_tibble(select(org.Sc.sgd.db, keys = sc_nodes,
                          c("ENSEMBL", "PATH"), "ENSEMBL"))
sc_indices <- unique(match(sc_go$ENSEMBL, sc_nodes))
chunks <- seidr_chunked_apply(
  sc_sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
  FUN = make_edge_tibble, indices = sc_indices, algorithms = algorithms(sc_sf),
  nodes = nodes(sc_sf)
)

chunks <-bind_rows(chunks)
sc_edges <- left_join(
  left_join(chunks, sc_go, by = c("Source" = "ENSEMBL")),
  sc_go, by = c("Target" = "ENSEMBL")
)
rm(chunks)
gc()

sc_edges %>%
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
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(y = M, fill = InPath, x = Algorithm)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
  th +
  xlab("") +
  ylab("Mean Edge Weight Rank per Pathway")
ggsave("~/seidr-plots/sc_kegg_edgerank.pdf", width = 16, height = 9)

at_sf <- SeidrFileFromPath("/mnt/picea/projects/seidr/networks/athal/results/aggregate/irp.sf")
at_nodes <- nodes(at_sf)
at_go <- select(org.At.tair.db, keys = at_nodes, c("TAIR", "PATH"), "TAIR")
at_indices <- unique(match(at_go$TAIR, at_nodes))
chunks <- seidr_chunked_apply(
  at_sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
  FUN = make_edge_tibble, indices = at_indices, algorithms = algorithms(at_sf),
  nodes = nodes(at_sf)
)
chunks <- bind_rows(chunks)
at_edges <- left_join(
  left_join(chunks, at_go, by = c("Source" = "TAIR")),
  at_go, by = c("Target" = "TAIR")
)
rm(chunks)
gc()

at_edges %>%
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
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(y = M, fill = InPath, x = Algorithm)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
  th +
  xlab("") +
  ylab("Mean Edge Weight Rank per Pathway")
ggsave("~/seidr-plots/at_kegg_edgerank.pdf", width = 16, height = 9)


dm_sf <- SeidrFileFromPath("/mnt/picea/projects/seidr/networks/dmel/results/aggregate/irp.sf")
dm_nodes <- nodes(dm_sf)
dm_go <- select(org.Dm.eg.db, keys = dm_nodes, c("FLYBASE", "PATH"), "FLYBASE")
dm_indices <- unique(match(dm_go$FLYBASE, dm_nodes))
chunks <- seidr_chunked_apply(
  dm_sf, scores = FALSE, ranks = TRUE, edge_index = TRUE,
  FUN = make_edge_tibble, indices = dm_indices, algorithms = algorithms(dm_sf),
  nodes = nodes(dm_sf)
)
chunks <- bind_rows(chunks)
dm_edges <- left_join(
  left_join(chunks, dm_go, by = c("Source" = "FLYBASE")),
  dm_go, by = c("Target" = "FLYBASE")
)
rm(chunks)
gc()

dm_edges %>%
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
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(y = M, fill = InPath, x = Algorithm)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
  th +
  xlab("") +
  ylab("Mean Edge Weight Rank per Pathway")
ggsave("~/seidr-plots/dm_kegg_edgerank.pdf", width = 16, height = 9)

# 
# sc_edges$ReprTerm <- NA
# i <- which(is.na(match(sc_edges$GOALL.x, go_terms$Term)))
# sc_edges$ReprTerm[i] <- sc_edges$GOALL.y[i]
# i <- which(is.na(match(sc_edges$GOALL.y, go_terms$Term)))
# sc_edges$ReprTerm[i] <- sc_edges$GOALL.x[i]
# i <- which(match(sc_edges$GOALL.x, go_terms$Term) == match(sc_edges$GOALL.y, go_terms$Term))
# sc_edges$ReprTerm[i] <- sc_edges$GOALL.x[i]
# i <- which(match(sc_edges$GOALL.x, go_terms$Term) != match(sc_edges$GOALL.y, go_terms$Term))
# sc_edges$ReprTerm[i] <- sc_edges$GOALL.x[i]
# 
# 
# sc_edges %>%
#   head(10000) %>%
#   pivot_longer(-c(Source, Target, MatchType, GOALL.x, EVIDENCEALL.x,
#                   ONTOLOGYALL.x, GOALL.y, EVIDENCEALL.y, ONTOLOGYALL.y, ReprTerm)) %>%
#   filter(name != "ARACNE") %>%
#   group_by(ReprTerm, name, MatchType) %>%
#   summarise(Mean = median(value, na.rm = TRUE)) %>%
#   ungroup() %>%
#   group_by(ReprTerm, name) %>%
#   summarise(Diff = diff(Mean) * -1) %>%
#   ggplot(aes(x = name, y = Diff)) +
#   geom_point() +
#   facet_wrap(~ReprTerm)
# 
# 
# at_edges$ReprTerm <- NA
# i <- which(is.na(match(at_edges$GOALL.x, go_terms$Term)))
# at_edges$ReprTerm[i] <- at_edges$GOALL.y[i]
# i <- which(is.na(match(at_edges$GOALL.y, go_terms$Term)))
# at_edges$ReprTerm[i] <- at_edges$GOALL.x[i]
# i <- which(match(at_edges$GOALL.x, go_terms$Term) == match(at_edges$GOALL.y, go_terms$Term))
# at_edges$ReprTerm[i] <- at_edges$GOALL.x[i]
# i <- which(match(at_edges$GOALL.x, go_terms$Term) != match(at_edges$GOALL.y, go_terms$Term))
# at_edges$ReprTerm[i] <- at_edges$GOALL.x[i]
# 
# triu <- length(at_nodes) * (length(at_nodes) - 1) / 2
# 
# at_edges[sample(1:nrow(at_edges), 1000000), ] %>%
#   pivot_longer(-c(Source, Target, MatchType, GOALL.x, EVIDENCEALL.x,
#                   ONTOLOGYALL.x, GOALL.y, EVIDENCEALL.y, ONTOLOGYALL.y, ReprTerm)) %>%
#   filter(name != "ARACNE") %>%
#   group_by(ReprTerm, name, MatchType) %>%
#   summarise(Score = mean(value, na.rm = TRUE)) %>%
#   ungroup() %>%
#   group_by(ReprTerm, name) %>%
#   summarise(Diff = diff(Score) * -1) %>%
#   ggplot(aes(x = name, y = Diff)) +
#   geom_linerange(aes(ymin = 0, ymax = Diff)) +
#   geom_point() +
#   facet_wrap(~ReprTerm) +
#   theme_bw()
# 
# at_edges[sample(1:nrow(at_edges), 1000000), ] %>%
#   pivot_longer(-c(Source, Target, MatchType, GOALL.x, EVIDENCEALL.x,
#                   ONTOLOGYALL.x, GOALL.y, EVIDENCEALL.y, ONTOLOGYALL.y, ReprTerm)) %>%
#   filter(name != "ARACNE") %>%
#   ggplot(aes(x = name, y = value, fill = c("Out", "In")[MatchType])) +
#   geom_boxplot() +
#   facet_wrap(~ReprTerm) +
#   stat_compare_means(label = "p.signif") +
#   theme_bw()


as_tibble(sc_edges) %>%
  mutate(across(ARACNE:irp, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })) %>%
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
  group_by(ReprPath, Algorithm) %>%
  summarise(Diff = diff(M) * -1) -> tmp
tmp$Diff[!is.finite(tmp$Diff)] <- 0
tmp %>%
  mutate(Algorithm = recode(Algorithm, irp = "Seidr")) %>%
  mutate(Algorithm = factor(Algorithm, levels = levels)) %>%
  ggplot(aes(y = M, fill = InPath, x = Algorithm)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(label = "p.signif",
                     method.args = list(alternative = "greater")) +
  theme_bw() +
  scale_fill_npg() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_text(face = c('bold', rep('plain', 13)))) +
  th +
  xlab("") +
  ylab("Mean Edge Weight Rank per Pathway")
ggsave("~/seidr-plots/sc_kegg_edgerank.pdf", width = 16, height = 9)

# tmp %>%
#   ungroup() %>%
#   group_by(ReprPath, Algorithm) %>%
#   summarise(Diff = diff(M) * -1) %>%
#   pivot_wider(names_from = "Algorithm", values_from = "Diff") -> tmp2
# m <- as.matrix(tmp2[, -1])
# s <- sign(m)
# m <- log10(abs(m))
# m <- m * s
# m[!is.finite(m)] <- 0
# pheatmap(m, scale = "none")
# 
# colSums(m)
# tmp %>%
 

