setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/covid19_splceosome/tft")
getwd()
library(dplyr)
library(readr)

filenames <- list.files(pattern="gsea_report", full.names=TRUE)
ldf <- lapply(filenames, read.csv, sep = "\t")
res <- lapply(ldf, summary)
names(ldf) <- substr(filenames, 22,100)
names(res) <- substr(filenames, 22,100)

infection_0h_tft <- 
  bind_rows(ldf[["pos_0h.tsv"]], 
            ldf[["neg_0h.tsv"]])
infection_1h_tft <- 
  bind_rows(ldf[["pos_1h.tsv"]], 
            ldf[["neg_1h.tsv"]])
infection_2h_tft <- 
  bind_rows(ldf[["pos_2h.tsv"]], 
            ldf[["neg_2h.tsv"]])
infection_4h_tft <- 
  bind_rows(ldf[["pos_4h.tsv"]], 
            ldf[["neg_4h.tsv"]])
infection_12h_tft <- 
  bind_rows(ldf[["pos_12h.tsv"]], 
            ldf[["neg_12h.tsv"]])
infection_16h_tft <- 
  bind_rows(ldf[["pos_16h.tsv"]], 
            ldf[["neg_16h.tsv"]])
infection_24h_tft <- 
  bind_rows(ldf[["pos_24h.tsv"]], 
            ldf[["neg_24h.tsv"]])

DFlist <- list(infection_0h_tft,
               infection_2h_tft,
               infection_4h_tft,
               infection_12h_tft,
               infection_16h_tft,
               infection_24h_tft)

filtered_DFlist <- lapply(DFlist,
                          FUN = function(x)
                          {
                            x <- x %>% 
                              filter(NOM.p.val<=0.05) %>% 
                              filter(FDR.q.val<=0.25)
                          })

names(filtered_DFlist) <- c("infection_0h_tft",
                            "infection_2h_tft",
                            "infection_4h_tft",
                            "infection_12h_tft",
                            "infection_16h_tft",
                            "infection_24h_tft")

multi_full <- Reduce(
  function(...) {
    full_join(..., 
              by = c("NAME" = "NAME"),
              keep = FALSE)
  },
  filtered_DFlist
)


tft_for_heatmap <- multi_full %>% 
  select(NAME, starts_with("NES"))

colnames(tft_for_heatmap) = c("NAME","NES.0h",
                               "NES.2h",
                               "NES.4h",
                               "NES.12h",
                               "NES.16h",
                               "NES.24h")
write.csv(tft_for_heatmap, 
          file = "tft_for_heatmap.csv")



#pheatmap----------------------------------------------------------
library(pheatmap)
library(RColorBrewer)
tft_for_heatmap <- 
  read.csv(file = "tft_for_heatmap.csv", 
           header = T)
tft_for_heatmap <- tft_for_heatmap %>% 
  arrange(desc(NES.4h)) %>% 
  filter(is.na(NES.4h) == FALSE)

tft_for_heatmap_tmp = 
  tft_for_heatmap[,c(3:8)]
row.names(tft_for_heatmap_tmp) = 
  tft_for_heatmap[,2]

pheatmap(tft_for_heatmap_tmp,
         color = colorRampPalette(rev(brewer.pal(n = 10, 
                                                 name = "RdYlBu")))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "grey60",
         na_col = "white",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 1)

# install.packages("colorRampPalette")
# install.packages("brewer.pal")
# install.packages("pheatmap")
