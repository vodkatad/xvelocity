library(ggplot2)
library(ggsankey)
library(ComplexHeatmap)

load('~/fortsnePlotP43W3C.RData')
models <- readRDS('/scratch/trcanmed/connector/local/share/data/Models.RDs')

md <- merge(models, fortsnePlot, by="ShortID")

df <- md %>%
  make_long("col", "Cluster")

ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)


data <- models
data <- data[, c(1,2)]
colnames(data) <- c('model', 'cluster')

cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)

data <- merge(data, cris, by.x="model", by.y="genealogy", all.x=TRUE)

snv <- read.table('/scratch/trcanmed/connector/local/share/data/muts/dtb_snv_april2021_velocity_all_rePR.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
cn <- read.table('/scratch/trcanmed/connector/local/share/data/muts/dtb_cn_april2021_velocity.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
colnames(cn)[4] <- 'ERBB2'

cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv')
colnames(cetuxi) <- c('smodel', 'w3')
data <- merge(data, cetuxi, by.x="model", by.y="smodel", all.x=TRUE)
data$cetuxi <- ifelse(data$w3 < -50, 'OR', ifelse(data$w3 > 35, 'PD', 'SD'))

data <- merge(data, snv, by.x="model", by.y="CASE", all.x=TRUE)
data <- merge(data, cn, by.x="model", by.y="CASE", all.x=TRUE)
data$PIK3CA <- NULL
# https://journal.r-project.org/archive/2016/RJ-2016-032/RJ-2016-032.pdf

# https://rdrr.io/bioc/TRONCO/src/R/visualization.R


rownames(data) <- data$model
data$model <- NULL

data <- data[order(data$cluster),]

odata <- data[, c('KRAS', 'NRAS', 'BRAF', 'MET', 'EGFR', 'ERBB2')]
odata$KRAS <- ifelse(odata$KRAS == "wt", '', 'mut')
odata$BRAF <- ifelse(odata$BRAF == "wt", '', 'mut')
odata$NRAS <- ifelse(odata$NRAS == "wt", '', 'mut')
odata$EGFR <- ifelse(odata$EGFR == "wt", '', 'ampl')
odata$MET <- ifelse(odata$MET == "wt", '', 'ampl')
odata$ERBB2 <- ifelse(odata$ERBB2 == "wt", '', 'ampl')
annot <- data[, c('cluster', 'cris', 'cetuxi', 'w3'), drop=FALSE]
col = c("mut" = "forestgreen", 'ampl'="red")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  mut = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*0.33, 
              gp = gpar(fill = col["mut"], col = NA))
  },
  ampl =function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*0.9, 
              gp = gpar(fill = col["ampl"], col = NA))
  }
)

#odata <- ifelse(odata==1, TRUE, FALSE)
#odata[is.na(odata)] <- FALSE
#opd <- list(Mutated=t(odata))
opd <- t(odata)


colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F")#,"#FF7F00","#CAB2D6")
nc <- unique(c(annot$cluster))
nc <- nc[!is.na(nc)]
nc <- nc[order(nc)]
names(colors) <- nc
connector_col <- colors

colors <- c('darkorange2', 'brown3','darkblue', 'darkgreen', 'aquamarine2')
nc <- unique(c(annot$cris))
nc <- nc[!is.na(nc)]
nc <- nc[order(nc)]
names(colors) <- nc
cris_col <- colors

colors <- c('steelblue', 'red', 'orange')
nc <- unique(c(annot$cetuxi))
nc <- nc[!is.na(nc)]
nc <- nc[order(nc)]
names(colors) <- nc
cetuxi_col <- colors

get_recist <- function(x) {
  res <- ifelse(x < -50, 'steelblue', ifelse(x > 35, 'red', 'orange'))
  return(res)
}

colors <- list(CONNECTOR=connector_col, CRIS=cris_col)#, Cetuximab=cetuxi_col)
svg('CTGC_A3.svg', width=13.578, height=3.38)#, units="in")
oncoPrint(opd, alter_fun = alter_fun, col = col, 
          column_order = rownames(odata),
          row_order = colnames(odata),
          remove_empty_columns = FALSE, remove_empty_rows = FALSE, 
          pct_gp=gpar(fontsize=10),
          column_names_gp=gpar(fontsize=10),
          top_annotation = HeatmapAnnotation(CONNECTOR=annot$cluster, 
                                             CRIS=annot$cris, 
                                             Cetuximab = anno_barplot(annot$w3, gp = gpar(fill = get_recist(annot$w3))),
                                             col=colors,
                                             gap = unit(2, "mm"),
                                             height = unit(4, "cm")),
          right_annotation = NULL,
          heatmap_height = unit(8, "cm"))
graphics.off()


fisher_two_way <- function(data, col1, col2, w1, w2) {
  df <- data
  df$isone <- ifelse(df[,col1]==w1, 'yes', 'no')
  df$istwo <- ifelse(df[,col2]==w2, 'yes', 'no')
  twoway <- table(df$isone, df$istwo)
  print(twoway)
  fisher.test(twoway)
}

# OR in Aa:
fisher_two_way(data, 'cluster', 'cetuxi', 'Aa', 'OR')
fisher_two_way(data, 'cluster', 'cetuxi', 'Ab', 'PD')


# PD in all:
lapply(c('Aa','Ab','Ac', 'Ba', 'Bb', 'Bc', 'C'),  function(x) {fisher_two_way(data, 'cluster', 'cetuxi', x, 'PD')} )

# depletio on muts in responder clusters:
data$mut <- apply(data[,c('KRAS','NRAS','BRAF', 'MET','ERBB2')], 1, function(x) {any(x!='wt')})
data$mut <- ifelse(is.na(data$mut),  FALSE, data$mut)

fisher_two_way(data, 'cluster', 'mut', 'Aa', TRUE)
fisher_two_way(data, 'cluster', 'mut', 'Ab', TRUE)

data$isBC <- ifelse(data$cluster %in% c("Ba", "Bb", "Bc", "C"), 'yes', 'no')
fisher_two_way(data, 'isBC', 'mut', 'yes', TRUE)

# cris-C in Aa
fisher_two_way(data, 'cris', 'cluster', 'CRIS-C', 'Aa')
# cris-C in OR
fisher_two_way(data, 'cetuxi', 'cris', 'OR', 'CRIS-C')

# odd ratios CrisC- Aa  and OR

data$isAa <- ifelse(data$cluster == 'Aa', 'Aa', 'notAa')
data$isCrisC <- ifelse(data$cris == 'CRIS-C', 'CRIS-C', 'notCRIS-C')
data$isOR <- ifelse(data$cetuxi == 'OR', 'OR', 'ORnot') # not not OR to have the correct ordering of rows/cols in table
fisher.test(table(data$isAa, data$isCrisC))
fisher.test(table(data$isOR, data$isCrisC))

# calcolo i due odds ratio (non quelli MLE conditional riportati da R)
#https://www.biochemia-medica.com/en/journal/19/2/10.11613/BM.2009.011/fullArticle

conn <- as.data.frame.matrix(table(data$isAa, data$isCrisC))

reci <- as.data.frame.matrix(table(data$isOR, data$isCrisC))

(conn[1,1]*conn[2,2]) / (conn[1,2] * conn[2,1])
(reci[1,1]*reci[2,2]) / (reci[1,2] * reci[2,1])


