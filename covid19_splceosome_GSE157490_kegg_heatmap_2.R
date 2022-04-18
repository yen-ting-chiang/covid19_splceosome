setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/covid19_splceosome")
getwd()
library(dplyr)
library(readr)

filenames <- list.files(pattern="gsea_report", full.names=TRUE)
ldf <- lapply(filenames, read.csv, sep = "\t")
res <- lapply(ldf, summary)
names(ldf) <- substr(filenames, 22,100)
names(res) <- substr(filenames, 22,100)

infection_0h_kegg <- 
  bind_rows(ldf[["pos_0h.tsv"]], 
            ldf[["neg_0h.tsv"]])
infection_2h_kegg <- 
  bind_rows(ldf[["pos_2h.tsv"]], 
            ldf[["neg_2h.tsv"]])
infection_4h_kegg <- 
  bind_rows(ldf[["pos_4h.tsv"]], 
            ldf[["neg_4h.tsv"]])
infection_12h_kegg <- 
  bind_rows(ldf[["pos_12h.tsv"]], 
            ldf[["neg_12h.tsv"]])
infection_16h_kegg <- 
  bind_rows(ldf[["pos_16h.tsv"]], 
            ldf[["neg_16h.tsv"]])
infection_24h_kegg <- 
  bind_rows(ldf[["pos_24h.tsv"]], 
            ldf[["neg_24h.tsv"]])
infection_36h_kegg <- 
  bind_rows(ldf[["pos_36h.tsv"]], 
            ldf[["neg_36h.tsv"]])

DFlist <- list(infection_0h_kegg,
               infection_2h_kegg,
               infection_4h_kegg,
               infection_12h_kegg,
               infection_16h_kegg,
               infection_24h_kegg,
               infection_36h_kegg)

filtered_DFlist <- lapply(DFlist,
                          FUN = function(x)
                          {
                            x <- x %>% 
                              filter(NOM.p.val<=0.05) %>% 
                              filter(FDR.q.val<=0.25)
                          })

names(filtered_DFlist) <- c("infection_0h_kegg",
                            "infection_2h_kegg",
                            "infection_4h_kegg",
                            "infection_12h_kegg",
                            "infection_16h_kegg",
                            "infection_24h_kegg",
                            "infection_36h_kegg")

multi_full <- Reduce(
  function(...) {
    full_join(..., 
              by = c("NAME" = "NAME"),
              keep = FALSE)
  },
  filtered_DFlist
)


kegg_for_heatmap <- multi_full %>% 
  select(NAME, starts_with("NES"))

colnames(kegg_for_heatmap) = c("NAME","NES.0h",
                               "NES.2h",
                               "NES.4h",
                               "NES.12h",
                               "NES.16h",
                               "NES.24h",
                               "NES.36h")
kegg_for_heatmap <- kegg_for_heatmap %>% 
  arrange(desc(NES.4h)) %>% 
  filter(is.na(NES.4h) == FALSE)

write.csv(kegg_for_heatmap, 
          file = "kegg_for_heatmap_for_selection.csv")



#pheatmap----------------------------------------------------------
library(pheatmap)
library(RColorBrewer)
kegg_for_heatmap <- 
  read.csv(file = "kegg_for_heatmap_selected.csv", 
           header = T)

kegg_for_heatmap_tmp = 
  kegg_for_heatmap[,c(3:8)]
row.names(kegg_for_heatmap_tmp) = 
  kegg_for_heatmap[,10]

pheatmap(kegg_for_heatmap_tmp,
         color = colorRampPalette(rev(brewer.pal(n = 10, 
                                                 name = "RdYlBu")))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "grey60",
         na_col = "white",
         show_rownames = TRUE)

# install.packages("colorRampPalette")
# install.packages("brewer.pal")
# install.packages("pheatmap")
