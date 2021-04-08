library(tidyverse)
library(ggridges)
library(scales)
library(ggpubr)
library(GGally)
library(ggsci)
library(magrittr)
# 
drought_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
            "preprosessed/seidr/newaggregate/aggregated.centrality")
)
# drought_centrality <- read_tsv(
#   file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
#             "preprosessed/seidr/newbackbone/backbone-10-percentBS.centrality")
# )
drought_centrality$Condition <- "Drought"

unstressed_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/facility/seidr/results/",
            "newaggregated.centrality")
)
# unstressed_centrality <- read_tsv(
#   file.path("/mnt/picea/projects/spruce/facility/seidr/results/",
#             "backbone-10-percent.centrality")
# )
unstressed_centrality$Condition <- "Unstressed"

centrality <- rbind(drought_centrality, unstressed_centrality)

gsoi <- read_lines(
  file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
            "preprosessed/seidr/newbackbone/GOI.txt"),
  skip = 1
)

set.seed(314154354)
gsrand <- sample(drought_centrality$Node, length(gsoi))

centrality$Type <- "Other"
centrality$Type[which(centrality$Node %in% gsoi)] <- "GeneOfInterest"
centrality$Type[which(centrality$Node %in% gsrand)] <- "Random"

centrality <- pivot_longer(
  centrality, -c(Condition, Node, Type), names_to = "Centrality", values_to = "Value"
)

comps <- list(
  c("GeneOfInterest", "Other"),
  c("Random", "Other"),
  c("GeneOfInterest", "Random")
)

####
# Comparison box plot
# Computational issues in the calculation of Katz centrality, but results are
# consistent in all other measures, so not being re-done.
####

centrality %>%
  filter(! Centrality %in% c("Katz")) %>%
  group_by(Condition, Centrality) %>%
  mutate(z_value = (Value - min(Value)) / (max(Value) - min(Value))) %>%
  mutate(z_value = z_value + abs(min(z_value, na.rm = TRUE))) %>%
  ungroup() %>%
  ggplot(aes(y = z_value,
             x = Type,
             fill = Type,
             group = Type)) +
  geom_boxplot(notch = TRUE) +
  stat_compare_means(comparisons = comps) +
  facet_grid(Condition~Centrality, scales = "free") +
  theme_bw() +
  #scale_y_log10(labels = comma) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_npg() +
  xlab("Gene Set") +
  ylab("Centrality [standard-scores]") +
  theme(axis.text.x = element_blank())
ggsave("~/seidr-plots/gsoi_centrality.pdf", width = 16, height = 9)


centrality %>%
  filter(Centrality != "Katz") %>%
  group_by(Condition, Centrality) %>%
  mutate(z_value = (Value - min(Value)) / (max(Value) - min(Value))) %>%
  mutate(z_value = z_value + abs(min(z_value, na.rm = TRUE))) %>%
  summarise(means = mean(z_value)) %>%
  ungroup() %>%
  ggplot(aes(y = z_value,
             x = paste(Type, Centrality),
             fill = Type,
             pch = Condition)) +
  stat_summary(fun.data = mean_ci, geom = "point", aes(col = Centrality), size = 3) +
  stat_summary(fun.data = mean_ci, geom = "line", aes(group = paste(Condition, Centrality)), lty = 3) +
  #stat_compare_means(comparisons = comps) +
  #facet_grid(Condition~Centrality, scales = "free") +
  theme_bw() +
  #facet_wrap(~Type, nrow = 1) +
  #scale_y_log10(labels = comma) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_npg() +
  annotate("rect", xmin = 6.5, xmax = 12.5, ymin = 0, ymax = Inf, alpha = 0.1) +
  xlab("Gene Set") +
  ylab("Centrality [z-scores]") +
  theme(axis.text.x = element_blank())

# drought_centrality %>% 
#   select(PageRank, Strength, Eigenvector, Katz, Laplacian, Closeness, Betweenness) %>%
#   apply(2, rank) %>%
#   as_tibble() %>%
#   ggpairs() +
#   theme_bw()

centrality %>%
  filter(! Centrality %in% c("Katz")) %>%
  group_by(Condition, Centrality) %>%
  mutate(z_value = (Value - min(Value)) / (max(Value) - min(Value))) %>%
  ungroup() %>%
  group_by(Condition, Centrality, Type) %>%
  summarise(Mean = mean(z_value)) %>%
  mutate(Diff = Mean - Mean[Type == "Other"]) %>%
  mutate(DDiff = Diff - abs(Diff[Type == "Random"])) %>%
  filter(Type == "GeneOfInterest") %>%
  ungroup() %>%
  ggplot(aes(y = DDiff,
             x = Condition,
             group = Centrality)) +
  geom_point(size = 5, pch = 19, aes(col = Centrality)) +
  #geom_line(lty = 3) +
  theme_bw() +
  #scale_y_log10(labels = comma) +
 # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 18)) +
  scale_fill_npg() +
  xlab("Condition") +
  ylab(expression(Delta["GOI"] - paste("|",Delta["Random"],"|")))
ggsave("~/seidr-plots/gsoi_centrality_magnitude.pdf", width = 16, height = 9)



library(GSEABase)
library(fgsea)
library(GOfuncR)
library(tidyverse)
library(magrittr)


# Get centrality metrics
drought_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
            "preprosessed/seidr/newaggregate/aggregated.centrality")
)
unstressed_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/facility/seidr/results/",
            "newaggregated.centrality")
)
# Get basic GO annotations. Version 1.1, 11/18/2020
gene_term <- read_tsv(
  "/mnt/picea/storage/reference/Picea-abies/v1.1/gopher/gene_to_go.tsv",
  col_names = c("Gene", "Term")
)
terms <- str_split(gene_term$Term, "\\|")

# Also add parental annotations, names and domain GOfuncR 1.10.0, default graph
gene_term <- tibble(Gene = rep(gene_term$Gene, sapply(terms, length)),
                    Term = unlist(terms))  %>%
  mutate(Domain = get_names(Term)$root_node) %>%
  mutate(Name = get_names(Term)$go_name)
parents <- as_tibble(get_parent_nodes(gene_term$Term))
gene_term <- left_join(gene_term, parents, by = c("Term" = "child_go_id"))

# Pivot
gene_term %<>%
  pivot_longer(-c(Gene, Domain, Name, parent_name, distance),
               names_to = "Type", values_to = "Term")
gene_term %<>% filter(!is.na(Term))

# Now set up GSEA annotations. Split genes by pathways
gene_term_bp <- filter(gene_term, Domain == "biological_process")
annot <- split(gene_term_bp$Gene, gene_term_bp$Term)
annot <- lapply(annot, unique)

# Score type 'pos' since we want a one sided test with alt=greater
gseas_drought <- apply(dplyr::select(drought_centrality, -c(Node, Katz)), 2, function(s){
  r <- rank(s)
  names(r) <- drought_centrality$Node
  gsea <- fgsea(pathways = annot, stats = r, scoreType = "pos", nproc = 8,
                eps = 0)
  left_join(gsea, get_names(gsea$pathway), by = c("pathway" = "go_id"))
})
gseas_drought <- do.call(rbind, lapply(seq_along(gseas_drought), function(i){
  x <- gseas_drought[[i]]
  x$Centrality <- names(gseas_drought)[i]
  x
}))
# Unstressed
gseas_unstressed <- apply(dplyr::select(unstressed_centrality, -c(Node, Katz)), 2, function(s){
  r <- rank(s)
  names(r) <- unstressed_centrality$Node
  gsea <- fgsea(pathways = annot, stats = r, scoreType = "pos", nproc = 8,
                eps = 0)
  left_join(gsea, get_names(gsea$pathway), by = c("pathway" = "go_id"))
})
gseas_unstressed <- do.call(rbind, lapply(seq_along(gseas_unstressed), function(i){
  x <- gseas_unstressed[[i]]
  x$Centrality <- names(gseas_unstressed)[i]
  x
}))
gseas_drought$Condition = "Drought"
gseas_unstressed$Condition = "Unstressed"

arrange(gseas_unstressed, padj) %>%
  filter(padj < 0.01) -> x

writexl::write_xlsx(x = list(
  unstressed = arrange(gseas_unstressed, padj) %>%
    filter(padj < 0.01),
  stressed = arrange(gseas_drought, padj) %>%
    filter(padj < 0.01)
), col_names = TRUE, path = "~/seidr-plots/gseas.xlsx")

cat(apply(x[!duplicated(x$pathway), ][, c("pathway", "padj")], 1, paste, collapse = " "), sep = "\n")

###
# Now we can summarise significant terms by their GO Slim ancestor
# to get a high level overview of the processes.
# GO Slim: 11/19/2020
###
gseas_drought %>% 
  group_by(go_name, pathway) %>% 
  summarise(mpadj = median(padj)) %>% 
  filter(mpadj < 0.01) %>% arrange(mpadj) %>% 
  ungroup() %>% 
  dplyr::select(-go_name) -> x
slim <- goSlim(
  GOCollection(x$pathway),
  getOBOCollection(paste0("http://current.geneontology.org/ontology/subsets/",
                          "goslim_plant.obo")),
  "BP"
)
slim %>% arrange(desc(Percent))

gseas_unstressed %>% 
  group_by(go_name, pathway) %>% 
  summarise(mpadj = median(padj)) %>% 
  filter(mpadj < 0.01) %>% arrange(mpadj) %>% 
  ungroup() %>% 
  dplyr::select(-go_name) -> x
slim <- goSlim(
  GOCollection(x$pathway),
  getOBOCollection(paste0("http://current.geneontology.org/ontology/subsets/",
                          "goslim_plant.obo")),
  "BP"
)
slim %>% arrange(desc(Percent)) 


###
# Chek GSEA for the curated set of drought genes and a random sample of the
# same magnitude
###

# Check GOIs
annot2 = list(GeneOfInterest = gsoi, Random = gsrand)
# Drought
gseas_drought <- apply(dplyr::select(drought_centrality, -c(Node, Katz)), 2, function(s){
  r <- rank(s)
  names(r) <- drought_centrality$Node
  gsea <- fgsea(pathways = annot2, stats = r, scoreType = "pos", nproc = 8,
                eps = 0)
  as_tibble(gsea)
})
gseas_drought <- do.call(rbind, lapply(seq_along(gseas_drought), function(i){
  x <- gseas_drought[[i]]
  x$Centrality <- names(gseas_drought)[i]
  x
}))
# Unstressed
gseas_unstressed <- apply(dplyr::select(unstressed_centrality, -c(Node, Katz)), 2, function(s){
  r <- rank(s)
  names(r) <- unstressed_centrality$Node
  gsea <- fgsea(pathways = annot2, stats = r, scoreType = "pos", nproc = 8,
                eps = 0)
  as_tibble(gsea)
})
gseas_unstressed <- do.call(rbind, lapply(seq_along(gseas_unstressed), function(i){
  x <- gseas_unstressed[[i]]
  x$Centrality <- names(gseas_unstressed)[i]
  x
}))
gseas_drought$Condition = "Drought"
gseas_unstressed$Condition = "Unstressed"
rbind(gseas_drought, gseas_unstressed)  %>%
  ggplot(aes(x = Centrality, y = -log10(padj), col = pathway)) +
  geom_point(size = 5) +
  theme_bw() +
  geom_label(data = tibble(x = 1, y = 2.05, l = "P=0.01", Condition = "Drought"),
             aes(x = x, y = y, label = l),
             col = "black") +
  facet_wrap(~Condition) +
  geom_hline(yintercept = -log10(0.01), col = "red", lty = 2) +
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_npg(name = "Gene set")
ggsave("~/seidr-plots/gsea_goi.pdf", width = 16, height = 9)
###
# Effect size plots
###
# With crossbars
rbind(gseas_drought, gseas_unstressed) %>%
  group_by(Condition, Centrality) %>%
  summarise(y = diff(log10(padj))) %>%
  ggplot(aes(x = Centrality, ymin = 0, ymax = y, col = Condition)) +
  geom_linerange(width = 0.25, size = 1, position = position_dodge(width = 0.25)) +
  theme_bw() +
  theme(text = element_text(size = 18))
# Just lines
rbind(gseas_drought, gseas_unstressed) %>%
  group_by(Condition, Centrality) %>%
  summarise(y = diff(log10(padj))) %>%
  ungroup() %>%
  group_by(Centrality) %>%
  summarise(y = diff(rev(y))) %>%
  ggplot(aes(x = Centrality, ymin = 0, ymax = y)) +
  geom_linerange(width = 0.25, size = 1, position = position_dodge(width = 0.25)) +
  theme_bw() +
  theme(text = element_text(size = 18))


###
# Get centrality metrics
drought_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
            "preprosessed/seidr/newaggregate/aggregated.centrality")
)
annot <- read_tsv("ftp://plantgenie.org/Data/ConGenIE/Picea_abies/v1.0/Annotation/Pabies1.0-b2g_ID-Desc.tsv.gz")
annot %>%
  mutate(SeqName = str_replace(SeqName, "\\.\\d+$", "")) -> annot
drought_centrality %>%
  select(-c(Katz, Betweenness)) %>%
  mutate(
    PageRank = rank(PageRank),
    Strength = rank(Strength),
    Eigenvector = rank(Eigenvector),
    Laplacian = rank(Laplacian),
    Closeness = rank(Closeness)
  ) %>%
  pivot_longer(-Node) %>%
  group_by(Node) %>%
  summarise(MedCentrRank = median(value)) %>%
  arrange(desc(MedCentrRank)) %>%
  left_join(annot, by = c("Node" = "SeqName")) %>%
  mutate(MedCentrRank = rev(rank(MedCentrRank)))

unknowns <- c("MA_10432675g0010", "MA_69177g0010", "MA_10192193g0010", "MA_419923g0010")

drought_centrality <- read_tsv(
  file.path("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/",
            "preprosessed/seidr/newbackbone/backbone-10-percentBS.centrality")
)
annot <- read_tsv("ftp://plantgenie.org/Data/ConGenIE/Picea_abies/v1.0/Annotation/Pabies1.0-b2g_ID-Desc.tsv.gz")
annot %>%
  mutate(SeqName = str_replace(SeqName, "\\.\\d+$", "")) -> annot
drought_centrality %>%
  select(-c(Katz, Betweenness)) %>%
  mutate(
    PageRank = rank(PageRank),
    Strength = rank(Strength),
    Eigenvector = rank(Eigenvector),
    Laplacian = rank(Laplacian),
    Closeness = rank(Closeness)
  ) %>%
  pivot_longer(-Node) %>%
  group_by(Node) %>%
  summarise(MedCentrRank = median(value)) %>%
  arrange(desc(MedCentrRank)) %>%
  left_join(annot, by = c("Node" = "SeqName")) %>%
  mutate(MedCentrRank = rev(rank(MedCentrRank))) %>%
  filter(Node %in% unknowns)
