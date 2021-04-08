setwd("/mnt/picea/projects/spruce/vhurry/drought-stress-needles/preprosessed/seidr/")
load("../analysis/salmon/dds.rda")
library(DESeq2)
library(tidyverse)
vsd <- varianceStabilizingTransformation(dds)
vst <- t(assay(vsd))
colnames(vst) <- gsub("\\.\\d+$", "", colnames(vst))


clusters <- read_tsv(
  file.path(
    "newbackbone/backbone-10-infomap/",
    "backbone-10-percentBS_directed_for_infomap.resolved"
  ),
  col_names = c("Path", "C1", "C2", "C3", "C4", "Flow", "N", "Gene"),
  col_types = "cccccnic"
)


tmp <- prcomp(vst[, filter(clusters, C1 == "2")$Gene])
tmp <- cbind(tibble(Expression = tmp$x[, 1]), colData(dds))
# # Make a temporary matrix
# tmp <- as_tibble(scale(vst[, filter(clusters, C1 == "1")$Gene]))
# # Cbind metadataonto it
# tmp <- cbind(tmp, colData(dds))
# # Make into long format
# tmp <- pivot_longer(tmp, -c(SampleName, SciLifeID, Level, SampleID),
#                     names_to = "Gene", values_to = "Expression")
# Plot
ggplot(tmp, aes(x = Level, y = Expression, group = 1)) +
  stat_summary(fun.data = mean_cl_boot, geom = "line", lwd = 10, col = "#E64B35FF") +
  geom_point(size = 10, pch = 21, fill = "#4DBBD5FF") +
  theme_void() +
  theme(text = element_text(size = 36))


tmp <- scale(vst[, filter(clusters, C1 == "2")$Gene])
tmp <- as_tibble(cbind(tmp, colData(dds)))
tmp %>% 
  pivot_longer(-c(SampleName, SciLifeID, Level, SampleID)) %>%
  ggplot(aes(x = Level, y = value, group = name)) +
  stat_summary(fun.data = mean_cl_boot, geom = "line", alpha = 0.05)

source("https://raw.githubusercontent.com/UPSCb/UPSCb-common/master/src/R/gopher.R")

expr_genes <- scan("genes.tsv", "")

enr_res <- lapply(1:16, function(cl) {
  message(paste("Cl: ", cl))
  res <- gopher(
    clusters$Gene[clusters$C1 == cl], expr_genes, alpha = 1,
    host = "https://franklin.upsc.se", port = 5432, url = "pabies",
    endpoint = "enrichment"
  )
  res
})

save(enr_res, file = "~/seidr-plots/enr_go.RData")


library(rrvgo)

# Results from gofer2
load("~bastian/seidr-plots/enr_go.RData")

terms <- enr_res[[3]]$go$id
padj <- enr_res[[3]]$go$padj

# Setup matrices
simMatrix <- calculateSimMatrix(
  terms,
  orgdb="org.Pabies.eg.db",
  ont="BP",
  method="Wang"
)
scores <- setNames(-log10(padj), terms)
reducedTerms <- reduceSimMatrix(
  simMatrix,
  scores,
  threshold=0.8,
  orgdb="org.Pabies.eg.db"
)

# Various plots
heatmapPlot(
  simMatrix,
  reducedTerms,
  annotateParent=TRUE,
  annotationLabel="parentTerm",
  fontsize=6
)

# Need to increase this option to see all labels
options("ggrepel.max.overlaps" = Inf)
scatterPlot(simMatrix, reducedTerms, size = "score")

treemapPlot(reducedTerms)

wordcloudPlot(reducedTerms, min.freq=1, colors="black")


for( f in c(1:16) ) {
  message(f)
  df <- dplyr::filter(enr_res[[f]]$go, namespace == "BP" & padj <= 0.1)
  terms <- df$id
  padj <- df$padj
  # Setup matrices
  simMatrix <- calculateSimMatrix(
    terms,
    orgdb="org.Pabies.eg.db",
    ont="BP",
    method="Wang"
  )
  scores <- setNames(-log10(padj), terms)
  reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores,
    threshold=0.25,
    orgdb="org.Pabies.eg.db"
  )
  pdf(paste0("~/seidr-plots/clusters/treemap_f_0.25_", f, ".pdf"), width = 16, height = 16)
  # x <- tapply(reducedTerms$score, reducedTerms$parentTerm, sum)
  # x <- data.frame(word = names(x), freq = x)
  # wordcloud::wordcloud(x$word, x$freq, colors = ggsci::pal_npg()(10), rot.per = 0, scale = c(4, 1))
  treemapPlot(reducedTerms, fontsize.labels = 20)
  dev.off()
}


as_tibble(scale(vst[, c("MA_10430013g0020", "MA_76955g0010", "MA_135627g0010",
                  "MA_103341g0010", "MA_121534g0010", "MA_10426365g0010")])) %>%
  cbind(as_tibble(colData(dds))) %>%
  pivot_longer(-c(SampleName, SciLifeID, Level, SampleID)) %>%
  ggplot(aes(x = Level, y = value, group = name, col = name)) +
  stat_summary(fun.data = mean_cl_boot, geom = "line") + 
  theme_bw() +
  ylab("Mean expression [VST]") +
  theme(text = element_text(size = 18)) +
  scale_color_discrete(name = "")
ggsave("~/seidr-plots/tfs_2.pdf", width = 16, height = 9)

writexl::write_xlsx(lapply(enr_res, "[[", "go"), path = "~/seidr-plots/go_enr.xlsx")
writexl::write_xlsx(lapply(enr_res, "[[", "mapman"), path = "~/seidr-plots/mapman_enr.xlsx")
