library(readr)
library(DESeq2)
library(tximport)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
setwd('/mnt/picea/home/bastian/experimental_area/SRAsalmon2')
logs <- dir('.', pattern="salmon_quant.log", full.names=TRUE, recursive=TRUE)

get_mapping_rate <- function(x) {
  lines <- readLines(x)
  i <- which(str_detect(lines, 'Mapping rate'))
  if (length(i) < 1) { return(NA) }
  n <- as.numeric(str_replace(str_extract(lines[i], '\\d+\\.\\d+%'), '%', ''))
  return(n)
}

get_sample <- function(x) {
  return(basename(dirname(dirname(x))))
}

metadata <- tibble(Sample = unlist(sapply(logs, get_sample)),
                   Rate = unlist(sapply(logs, get_mapping_rate)))

ggplot(metadata, aes(x = Rate)) +
  geom_freqpoly(binwidth = 5, lwd = 1) + theme_bw()

metadata %>% filter(Rate > 75) -> f_metadata

in_samples <- file.path(f_metadata$Sample, 'quant.sf')

t2g <- read_tsv('tx2gene.tsv', col_names = c('Tx', 'Gene'))

txi <- tximport(in_samples, 'salmon', tx2gene = t2g)

meta <- data.frame(x = factor(1:length(in_samples)))

dds <- DESeqDataSetFromTximport(txi, meta, ~x)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

vss <- assay(vsd)

pca <- prcomp(t(vss))

plot(pca$x[, 1], pca$x[, 2])

vss <- apply(vss, 2, function(f) { return(f - median(f)) } )
vss <- t(vss)

write_tsv(as.data.frame(vss), 'expression_data.tsv')

save.image('scere.RData')
