library(ggplot2)
load('~/fortsnePlotP43W3C.RData')
fortsnePlot$sha <- 1 + (  max(fortsnePlot$ShannonIndex)  - fortsnePlot$ShannonIndex)

cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE)
write.table(fortsnePlot, file= '/scratch/trcanmed/connector/local/share/data/tsne_cetuxi_3w.tsv', quote=F, sep="\t", row.names = F)

m_cris <- merge(fortsnePlot, cris, by.x="ShortID", all.x=TRUE, by.y="genealogy")


ggplot(data=m_cris, aes(x=tsne1, y=tsne2, size=sha, color=cris))+geom_point()+theme_bw()

# Check vs their data (muts mainly, can be done by hand on examples)
# can probalby be done with the merge with methy, that has lots of info.

# methy classes
methy <- read.table("/scratch/trcanmed/pdx_methylation/local/share/data/annotations/All_samples_info_final_070322_wclusters-cn_CRISvsd_enr.tsv", sep="\t", header=T)

m_methy <- merge(fortsnePlot, methy, by.x="ShortID", all.x=TRUE, by.y="model")
m_methy$cluster <- as.factor(m_methy$cluster)
ggplot(data=m_methy, aes(x=tsne1, y=tsne2, size=sha, color=cluster))+geom_point()+theme_bw()

# MUTS (only top hits from Marco's enrichments)


# DEG


# growth speed of PDO (? just a random idea)
ctg <- read.table('/scratch/trcanmed/biobanca/local/share/data/CTG_growth_X_CORR_VALUES.txt', sep="\t",header=T, comment.char="")
colnames(ctg)[1] <- 'smodel'
colnames(ctg)[3] <- 'CTG_75k'

m <- merge(fortsnePlot, ctg, by.x="ShortID", all.x=TRUE, by.y="smodel")
ggplot(data=m, aes(x=tsne1, y=tsne2, size=sha, color=CTG_75k))+geom_point()+theme_bw()
ggplot(data=m, aes(x=tsne1, y=tsne2, size=sha, color=log(CTG_75k)))+geom_point()+theme_bw()

# cetuximab in vivo
samples_S <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/samples_data', sep="\t", header=T)
samples_R <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_R/samples_data', sep="\t", header=T)

fortsnePlot$hasDE <- ifelse(fortsnePlot$ShortID %in% c(as.character(samples_S$sample), as.character(samples_R$sample)), 'yes', 'no')
ggplot(data=fortsnePlot, aes(x=tsne1, y=tsne2, size=sha, color=hasDE))+geom_point()+theme_bw()


##
tp <- read.table('/tmp/tp53', sep="\t", header=F)
tp$smodel <- substr(tp$V4, 0, 7)
tp$lfc <- tp$V5

m_tp53 <- merge(fortsnePlot, tp, by.x="ShortID", all.x=TRUE, by.y="smodel")


ggplot(data=m_tp53, aes(x=tsne1, y=tsne2, size=sha, color=lfc))+geom_point()+theme_bw()
ggplot(data=m_tp53, aes(x=col, y=lfc))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw()+ggtitle('TP53 lfc targeted sanger')


## enrich cris C in Aa vs C in OR
load('~/fortsnePlotP43W3C.RData')
fortsnePlot$sha <- 1 + (  max(fortsnePlot$ShannonIndex)  - fortsnePlot$ShannonIndex)

cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE)



cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/november2021/cetuxi_3w.tsv', sep="\t", header=F, stringsAsFactors = FALSE)
colnames(cetuxi) <- c('smodel', 'perc')
cetuxi$recist <- ifelse(cetuxi$perc < -50, 'OR', ifelse(cetuxi$perc > 35, 'PD', 'SD'))

m_cris_reci <- merge(cris, cetuxi, by.x="genealogy", by.y="smodel")

if (!length(unique(m_cris_reci$genealogy)) == nrow(m_cris_reci)) {
  stop('Something wrong!')
}


mm_cris_reci <- merge(m_cris_reci, fortsnePlot, by.x="genealogy", by.y="ShortID")

if (!length(unique(mm_cris_reci$genealogy)) == nrow(mm_cris_reci)) {
  stop('Something wrong!')
}

mm_cris_reci$isAa <- ifelse(mm_cris_reci$col == 'Aa', 'Aa', 'notAa')
mm_cris_reci$isCrisC <- ifelse(mm_cris_reci$cris == 'CRIS-C', 'CRIS-C', 'notCRIS-C')
mm_cris_reci$isOR <- ifelse(mm_cris_reci$recist == 'OR', 'OR', 'ORnot') # not not OR to have the correct ordering of rows/cols in table
fisher.test(table(mm_cris_reci$isAa, mm_cris_reci$isCrisC))
fisher.test(table(mm_cris_reci$isOR, mm_cris_reci$isCrisC))

# calcolo i due odds ratio (non quelli MLE conditional riportati da R)
#https://www.biochemia-medica.com/en/journal/19/2/10.11613/BM.2009.011/fullArticle

conn <- as.data.frame.matrix(table(mm_cris_reci$isAa, mm_cris_reci$isCrisC))

reci <- as.data.frame.matrix(table(mm_cris_reci$isOR, mm_cris_reci$isCrisC))

(conn[1,1]*conn[2,2]) / (conn[1,2] * conn[2,1])
(reci[1,1]*reci[2,2]) / (reci[1,2] * reci[2,1])

### selection high

data <- read.table('/home/mferri/allchrts_Pd_high_low_krt.tsv', sep="\t", header=TRUE)
high <- as.data.frame(apply(data, 1, function(x) {length(x[x=='high'])}))
low <- as.data.frame(apply(data, 1, function(x) {length(x[x=='low'])}))

selh <- head(high[order(-high[,1]),,drop=F], n=6)
sell <- head(low[order(-low[,1]),,drop=F], n=5)

data[rownames(data) %in% rownames(selh),]
data[rownames(data) %in% rownames(sell),]

high <- c('CRC0438LMX0A02001TUMR01000', 'CRC0480LMX0B02001TUMR02000', 'CRC1562LMX0A02001TUMR01000', 'CRC1586LMX0B02001TUMR01000')
low <- c('CRC1144LMX0A02002TUMR01000','CRC1342LMX0B02001TUMR01000', 'CRC0276LMX0A02001TUMR04000', 'CRC0420LMX0A02001TUMR06000')

smodels_h <- substr(high, 0, 7)
smodels_l <- substr(low, 0, 7)
mm_cris_reci[mm_cris_reci$genealogy %in% smodels_h,c('genealogy', 'cris', 'perc', 'col')]
mm_cris_reci[mm_cris_reci$genealogy %in% smodels_l,c('genealogy', 'cris', 'perc', 'col')]


# Marti's code
merged2$krtclass <- ifelse(merged2$genealogy %in% low, 'Low', ifelse(merged2$genealogy %in% high, 'High', 'Ignavi'))

ggplot(data = merged2, aes(x = col, y = KRT6A, color = krtclass)) + geom_jitter() + ylab('KRT6A') + scale_color_manual(values = c("Low"="steelblue", "Ignavi"="gray", "High"="red"))

ggplot(data = merged2, aes(x = col, y = KRT80, color = krtclass)) + geom_jitter() + ylab('KRT80') + scale_color_manual(values = c("Low"="steelblue", "Ignavi"="gray", "High"="red"))


###
final2 <- read.table('/home/mferri/allchrts_Pd_high_low_krt.tsv', sep ="\t", header=T)
class <- read.table("/home/mferri/Pd_high_low_krt.tsv", sep = "\t", header=T)


high <- as.data.frame(apply(final2, 1, function(x) {length(x[x=='high'])}))
low <- as.data.frame(apply(final2, 1, function(x) {length(x[x=='low'])}))

selh <- head(high[order(-high[,1]),,drop=F], n=6)
#sell <- head(low[order(-low[,1]),,drop=F], n=8)

final2[rownames(final2) %in% rownames(selh),]
#final2[rownames(final2) %in% rownames(sell),]

high <- c('CRC0438LMX0A02001TUMR01000', 'CRC0480LMX0B02001TUMR02000', 'CRC1562LMX0A02001TUMR01000', 'CRC1586LMX0B02001TUMR01000')

colnames(low) <- 'nlow'
mm <- merge(low, class, by="row.names")

low <- c('CRC1144LMX0A02002TUMR01000','CRC1342LMX0B02001TUMR01000', 'CRC0276LMX0A02001TUMR04000', 'CRC0420LMX0A02001TUMR06000')

class[rownames(class) %in% high,]
class[rownames(class) %in% low,]

cris <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/cris_vsd_ok_prediction_result.tsv', sep="\t", header=T)

cris$smodel <- substr(cris$sample.names, 0,7)
class$smodel <- substr(rownames(class), 0,7)

m1 <- merge(cris, class, by.x="sample.names", by.y="row.names")
m2 <- merge(m1, fortsnePlot, by.x="smodel.x", by.y="ShortID")

recist_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3.tsv"
recist <- read.table(recist_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
recist$classification <- NA

for (i in seq(length(recist$case))) {
  if (recist[i, "perc"] < -50.0){
    recist[i, "classification"] <- "OR"
  } else if (recist[i, "perc"] > -50.0 & recist[i, "perc"] < 35.0) {
    recist[i, "classification"] <- "SD"
  } else {
    recist[i, "classification"] <- "PD"
  }
}
names(recist)[names(recist) == "case"] <- "ShortID"
colnames(recist)[3] <- 'recist'
m3 <- merge(m2, recist, by.x="smodel.x", by.y="ShortID")
res <- m3[m3$recist=="PD",]


res$isCRISB <- ifelse(res$predict.label2 == "CRIS-B", 'yes', 'no')
res$isKRThigh <- ifelse(res$KRT_classification == "KRT_high", 'yes', 'no')
table(res$isCRISB, res$isKRThigh)


d <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/cris_deg_krt_connector/krt_cris_smodel_collapse.tsv', sep="\t", header=T, stringsAsFactors = F)


data <- read.table('/home/mviviani/misc_tmp/CONNECTOR/Jun2022/June2022_cluster-prepared.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
data <- data[, c(1,2,3,4)]

snv <- read.table('/scratch/trcanmed/connector/local/share/data/muts/dtb_snv_april2021_velocity_all_rePR.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
cn <- read.table('/scratch/trcanmed/connector/local/share/data/muts/dtb_cn_april2021_velocity.tsv', sep="\t", stringsAsFactors = FALSE, header=TRUE)
colnames(cn)[4] <- 'ERBB2'

cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv')
colnames(cetuxi) <- c('smodel', 'w3')
data <- merge(data, cetuxi, by.x="model", by.y="smodel", all.x=TRUE)

data <- merge(data, snv, by.x="model", by.y="CASE", all.x=TRUE)
data <- merge(data, cn, by.x="model", by.y="CASE", all.x=TRUE)
data$PIK3CA <- NULL

m <- merge(data, d, by.x="row.names", by.y="model")
