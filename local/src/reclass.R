library(ggplot2)
library(readxl)
bayes <- readRDS('/scratch/trcanmed/connector/local/share/data/second_round_bayes.rds')
table(bayes$col)
#cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
#data <- merge(data, cris, by.x="model", by.y="genealogy", all.x=TRUE)

cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/november2021/recist_3w.tsv', sep="\t", stringsAsFactors = F, header=T)
### madplayer recovering recistssss
cetuxinew <- read_excel('/scratch/trcanmed/pdx_methylation/local/share/data/local_dfs/DTB_Treatments_Update JUNE 2022_.xlsx', sheet=3)
# cetuxijanmods <- read.table('/scratch/trcanmed/pdx_methylation/local/share/data/local_dfs/modello_ctx_jan', header=FALSE)
# cetuxijanmods <- read_excel("/scratch/trcanmed/pdx_methylation/local/share/data/local_dfs/DTB_Treatments_Update 2022_Jan.xlsx")
cetuxinew <- as.data.frame(cetuxinew[-(324:nrow(cetuxinew)),])
write.table(cetuxinew, "/scratch/trcanmed/pdxopedia/local/share/data/treats/june2022/DTB_treats_cetuxi_june2022.tsv", sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

data <- merge(bayes, cetuxi, by.x="ShortID", by.y="smodel", all.x=TRUE)
na <- data[is.na(data$cetuxi_3w),]
# cetuxinew <- cetuxinew[,c(3,6,11,12,14,15)]
# names(cetuxinew) <- c("model","recist","cetuxi_3w","cetuxi_6w","cetuxi_3w.a","cetuxi_6w.a")
nacetuxi <- cetuxinew[cetuxinew$CASE %in% na$ShortID,]
### we got 12 new recist!!
for( m in 1:length(nacetuxi$CASE) ) {
  data[data$ShortID==nacetuxi$CASE[m],"recist"] <- nacetuxi[nacetuxi$CASE==nacetuxi$CASE[m],"Response Class 3wks"]
  data[data$ShortID==nacetuxi$CASE[m],"cetuxi_3w"] <- nacetuxi[nacetuxi$CASE==nacetuxi$CASE[m],"Column2"]
}
check <- data[is.na(data$recist),]
saveRDS(data, file="/scratch/trcanmed/connector/local/share/data/recistRecovered.rds")
# bayesrec <- readRDS('/scratch/trcanmed/connector/local/share/data/recistRecovered.rds')
###


df <- as.data.frame(table(data$col, data$recist))
colnames(df) <- c('connector', 'recist', 'n')
df$recist <- factor(df$recist, levels=c('OR', 'SD', 'PD'))
ggplot(data =  df, mapping = aes(x = recist, y = connector)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() + theme(legend.position = "none")

### OR/SD in Ac/Ad investigations: madplayer
ador <- data[(data$col=='Ad' | data$col=='Ac'),]
ador <- ador[order(ador$ShannonIndex),]
ador <- ador[!is.na(ador$recist),]
ador$recist <- factor(ador$recist, levels=c('OR', 'SD', 'PD'))

ggplot(data=ador, aes(x=recist, y=ShannonIndex)) +
  geom_boxplot() +
  geom_jitter(aes(color=col)) +
  theme_bw()

forplot <- data[!is.na(data$recist),]
forplot$recist <- factor(forplot$recist, levels=c('OR', 'SD', 'PD'))
ggplot(data=data, aes(x=recist, y=ShannonIndex)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=col)) +
  theme_bw()
### it seems that the shannon index in PDs tends to < 1 and many Ad/Ac that are misclassified
### with recist, has the same < 1 shannon
###

models <- readRDS('/scratch/trcanmed/connector/local/share/data/Models.RDs')

length(intersect(models$ShortID, bayes$ShortID))
length(setdiff(bayes$ShortID, models$ShortID))

common <- intersect(models$ShortID, bayes$ShortID)

md <- merge(models, bayes, by="ShortID")
df2 <- as.data.frame(table(md$Cluster, md$col))
colnames(df) <- c('', 'recist', 'n')
ggplot(data =  df, mapping = aes(x = recist, y = connector)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() + theme(legend.position = "none")

# nah this is less clusters
load('~/fortsnePlotP43W3C.RData')
bayes <- readRDS('/scratch/trcanmed/connector/local/share/data/second_round_bayes.rds')
length(intersect(fortsnePlot$ShortID, bayes$ShortID))
common <- intersect(fortsnePlot$ShortID, bayes$ShortID)

md <- merge(fortsnePlot, bayes, by="ShortID")
df2 <- as.data.frame(table(md$col.x, md$col.y))
colnames(df2) <- c('Old', 'New', 'n')
df2$New <- factor(df2$new, levels=unique(df2$Old))
ggplot(data =  df2, mapping = aes(x = Old, y = New)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() + theme(legend.position = "none")


#### irinotecan
bayes <- readRDS('/scratch/trcanmed/connector/local/share/data/irinoClassificationWithBayesAssociation.RDs')

table(bayes$col)
#cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
#data <- merge(data, cris, by.x="model", by.y="genealogy", all.x=TRUE)

irino <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/february2023/irino_feb23.txt', sep="\t", stringsAsFactors = F, header=T)
irino <- irino[, c('CASE','X3WKS', 'X6WKS', 'Response.Class.3wks')]
colnames(irino) <- c('CASE', 'p3w', 'p6w', 'recist')
data <- merge(bayes, irino, by.x="ShortID", by.y="CASE", all.x=TRUE)

df <- as.data.frame(table(data$col, data$recist))
colnames(df) <- c('connector', 'recist', 'n')
df$recist <- factor(df$recist, levels=c('PR', 'SD', 'PD'))
ggplot(data =  df, mapping = aes(x = recist, y = connector)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() + theme(legend.position = "none")

#### oncoprint
library(ComplexHeatmap)
bayes <- readRDS('/scratch/trcanmed/connector/local/share/data/second_round_bayes.rds')
class <- readRDS("/scratch/trcanmed/connector/local/share/data/OurClassification.RDs")

snv <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/dtb_mutations_update060623_expandedPRLM.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
snv <- snv[snv$X != "PR",]
cn <- read.table('/scratch/trcanmed/connector/local/share/data/muts/dtb_cn_april2021_velocity.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
colnames(cn)[4] <- 'ERBB2'


cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=T, stringsAsFactors = F)
#cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/november2021/recist_3w.tsv', header=T, sep="\t")
colnames(class) <- c('ShortID', 'w3', 'recist')
data <- merge(bayes, class, by="ShortID", all.x=TRUE)

data <- merge(data, snv, by.x="ShortID", by.y="CASE", all.x=TRUE)
data <- merge(data, cn, by.x="ShortID", by.y="CASE", all.x=TRUE)
data <- merge(data, cris, by.x="ShortID", by.y="genealogy", all.x=TRUE)
data$PIK3CA <- NULL

rownames(data) <- data$ShortID
data$model <- NULL

colnames(data)[2] <- 'cluster'
data <- data[order(data$cluster),]

odata <- data[, c('KRAS', 'NRAS', 'BRAF', 'MET', 'EGFR', 'ERBB2')]
odata$KRAS <- ifelse(odata$KRAS == "wt", '', 'mut')
odata$BRAF <- ifelse(odata$BRAF == "wt", '', 'mut')
odata$NRAS <- ifelse(odata$NRAS == "wt", '', 'mut')
odata$EGFR <- ifelse(odata$EGFR == "wt", '', 'ampl')
odata$MET <- ifelse(odata$MET == "wt", '', 'ampl')
odata$ERBB2 <- ifelse(odata$ERBB2 == "wt", '', 'ampl')
#annot <- data[, c('cluster', 'cris', 'cetuxi', 'w3'), drop=FALSE]
annot <- data[, c('cluster','cris', 'recist', 'w3'), drop=FALSE]
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


colors <- c("#a8cddf","#a6d271","#2579b5","#339b2a","#d0b8da","#ed9798","#e8151c","#ffbb64","#ffe2a3")
nc <- #unique(c(annot$cluster))
nc <- c('Aa', 'Ab', 'Ac', 'Ad', 'Ae', 'Ba', 'Bb', 'Bc', 'C')
# try to guess the missing one
#library(RColorBrewer)
#pp <-  brewer.pal(n=length(nc), 'Paired')

setdiff(pp, toupper(colors))
nc <- nc[!is.na(nc)]
nc <- nc[order(nc)]
names(colors) <- nc
connector_col <- colors

colors <- c('darkorange2', 'brown3','darkblue', 'darkgreen', 'aquamarine2', 'black')
nc <- c('CRIS-A', 'CRIS-B', 'CRIS-C', 'CRIS-D', 'CRIS-E', 'HET')
nc <- nc[!is.na(nc)]
nc <- nc[order(nc)]
names(colors) <- nc
cris_col <- colors
names(cris_col) <- 
# colors <- c('steelblue', 'red', 'orange')
# nc <- unique(c(annot$cetuxi))
# nc <- nc[!is.na(nc)]
# nc <- nc[order(nc)]
# names(colors) <- nc
# cetuxi_col <- colors

get_recist <- function(x) {
  res <- ifelse(x < -50, 'steelblue', ifelse(x > 35, 'red', 'orange'))
  return(res)
}

colors <- list(CONNECTOR=connector_col, CRIS=cris_col)#, Cetuximab=cetuxi_col)

svg('/scratch/trcanmed/connector/local/share/data/connector_eacr_paper.svg', width=13.578, height=3.38)#, units="in")
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
save.image('/scratch/trcanmed/connector/local/share/data/connector_eacr.RData')

