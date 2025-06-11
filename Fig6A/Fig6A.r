library(tidyverse)
library(tximport)
library(DESeq2)
library(ComplexHeatmap)
library(khroma)
library(colorRamp2)
library(immunedeconv)
theme_set(theme_bw(base_size = 25))

dds <- readRDS("./files/HLA_neg_vs_pos.dds.rld.rds")
info <- readxl::read_excel("./files/1_頸癌bulk_metafile_v2.xlsx") %>%
  select(new_sample, histology, HPV, HPV4pattern, HLA_pos, met_pos, CD8_pos, infiltration, PDL1pos) %>%
  filter(histology != "MIX")
info[sapply(info, is.character)] <- lapply(info[sapply(info, is.character)], as.factor)
dds <- dds[,info$new_sample]
colData(dds) <- cbind(colData(dds), info)
gtf <- rtracklayer::import("../250327_Figure/files/gencode.v45.basic.annotation.success.sorted.gtf")
convert <- mcols(gtf[gtf$type == "gene"]) %>%
  as.data.frame() %>%
  .[, c("gene_id", "gene_name")] %>%
  tibble()

#heatmap
files <- list.files("../241227_CIBERSORT/res", full.names = TRUE)
data <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) %>%
  .$abundance
colnames(data) <- list.files("../241227_CIBERSORT/res") %>% str_remove_all("\\.genes.results")
old_ids <- rownames(data)
new_rownames <- convert$gene_name[match(old_ids, convert$gene_id)]
new_rownames[is.na(new_rownames)] <- old_ids[is.na(new_rownames)]
rownames(data) <- new_rownames

info <- readxl::read_excel("./files/1_頸癌bulk_metafile_v2.xlsx", sheet = "original") %>%
  select(newHLA, Histology, CD8_number, filt_absent, NewHLApos, met_cutoff, anyLOH, PDL1pos, HPV_positive) %>%
  dplyr::rename(new_sample = newHLA)

order_index <- order(info$Histology, info$CD8_number, info$NewHLApos)
info_ordered <- info[order_index, ]
info_ordered$NewHLApos <- str_replace_all(info_ordered$NewHLApos, c("neg" = "negative", "pos" = "positive"))
info_ordered$PDL1pos <- str_replace_all(info_ordered$PDL1pos, c("neg" = "negative", "pos" = "positive"))
info_ordered$filt_absent <- str_replace_all(info_ordered$filt_absent, c("0" = "absent", "1" = "infiltrated+excluded"))
info_ordered <- info_ordered %>%
  select(Histology, CD8_number, filt_absent, NewHLApos, met_cutoff, anyLOH, PDL1pos, HPV_positive)
info_ordered$Histology <- factor(info_ordered$Histology, levels = c("SCC", "AC", "GAS", "Small", "MIX"))

res.MCPcounter <- deconvolute_mcp_counter(gene_expression_matrix = data, feature_types = "HUGO_symbols")
mat <- t(t(res.MCPcounter) / colSums(res.MCPcounter))
mat <- cbind(mat, matrix(0, nrow = 10, ncol = nrow(info) - ncol(mat)))
colnames(mat) <- info$new_sample
mat_ordered <- mat[, order_index]

colnames(info_ordered)[2] <- "CD8_score"
colnames(info_ordered)[3] <- "CD8 infiltration"
colnames(info_ordered)[4] <- "HLA-I"
colnames(info_ordered)[5] <- "methylation"
colnames(info_ordered)[6] <- "LOH"
colnames(info_ordered)[7] <- "PDL1"
colnames(info_ordered)[8] <- "HPV"

pallet <- colour("light")(9)
colour_list <- list(Histology = c("SCC" = "red", "AC" = "green", "GAS" = "lightgreen", "Small" = "blue", "MIX" = "#FFC800"),
                    CD8_score = colorRamp2(c(0, 20), c("white", "#CC1669")),
                    `CD8 infiltration` = c("infiltrated+excluded" = "#CC1669", "absent" = "#325E9A"),
                    `HLA-I` = c("negative" = "#00aeda", "positive" = "#003b73"),
                    methylation = c("non_methyl" = "#F58F39", "methyl" = "#8220F5"),
                    LOH = c("LOH" = "steelblue", "ROH" = "tomato"),
                    PDL1 = c("negative" = "#38A89D", "positive" = "#E07A5F"),
                    HPV = c("positive" = "#F067A6", "negative" = "#406F79"))
annotation <- HeatmapAnnotation(df = info_ordered, col = colour_list)
hm <- Heatmap(matrix(nrow = 0, ncol = ncol(mat_ordered)),
              top_annotation = annotation,  column_split = info_ordered$Histology)
pdf(file = "heatmap_simple.pdf", width = 12, height = 4)
draw(hm, annotation_legend_side = "bottom")
dev.off()

