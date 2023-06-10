
library(ComplexHeatmap)

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


#colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6")
#colors <- c("#33a02c","#fb9a99","#a6cee3","#1f78b4","#e31a1c","#ff7f00","#b2df8a","#f48114","#fdbf6f")
colors <- c("#a8cddf","#a6d271","#2579b5","#339b2a","#d0b8da","#ed9798","#e8151c","#ffbb64","#ffe2a3")

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

svg('/scratch/trcanmed/connector/local/share/data/CTGC_A5_paper.svg', width=13.578, height=3.38)#, units="in")

colors <- list(CONNECTOR=connector_col, CRIS=cris_col)#, Cetuximab=cetuxi_col)
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


data$mut <- apply(data[,c('KRAS','NRAS','BRAF', 'MET','ERBB2')], 1, function(x) {any(x!='wt')})
data$mut <- ifelse(is.na(data$mut),  FALSE, data$mut)

fisher.test(table(data$cluster, data$mut))
chisq.test(table(data$cluster, data$mut))

fisher_two_way <- function(data, col1, col2, w1, w2) {
  df <- data
  df$isone <- ifelse(df[,col1]==w1, 'yes', 'no')
  df$istwo <- ifelse(df[,col2]==w2, 'yes', 'no')
  twoway <- table(df$isone, df$istwo)
  print(twoway)
  fisher.test(twoway)
}

fisher_two_way(data, 'cluster', 'cetuxi', 'Aa', 'OR')

lapply(c('Aa','Ab','Ac', 'Ba', 'Bb', 'Bc', 'C'),  function(x) {fisher_two_way(data, 'cluster', 'cetuxi', x, 'PD')} )

fisher_two_way(data, 'cluster', 'mut', 'Aa', TRUE)
fisher_two_way(data, 'cluster', 'mut', 'Ab', TRUE)

fisher_two_way(data, 'cluster', 'cris', 'Aa', 'CRIS-C')
fisher_two_way(data, 'cetuxi', 'cris', 'OR', 'CRIS-C')


data$isBC <- ifelse(data$cluster %in% c("Ba", "Bb", "Bc", "C"), 'yes', 'no')
fisher_two_way(data, 'isBC', 'mut', 'yes', TRUE)

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

## begin volcano
#egrassi@godot:/scratch/trcanmed/DE_RNASeq/dataset/connector_for_all$ for f in  *_*/*deseq2*tsv; do bawk -v f=$f '{print f, $1,$2, $3,$6}' $f; done | grep -v baseMean > volcano

vd <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/connector_Aa_vsall/volcano', sep="\t", header=F, stringsAsFactors = F)
colnames(vd) <- c('v','gene', 'base', 'lfc', 'padj')


vd$vv <- sapply(strsplit(vd$v, '/'), function(x) {x[1]})
vd$vs <- sapply(strsplit(vd$vv, '_'), function(x) {x[2]})
vd$gene <- gsub('H_','', vd$gene, fixed=T)

library(reshape)

padj <- vd[vd$vs == "Bb", ]
data <- vd[, c('gene', 'lfc', 'vs')]
odata <- padj[, c('gene', 'padj', 'vs')]

m <- merge(data, odata, by="gene")
go <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/connector_Aa_vsall/GO0031424_keratinization.tsv', sep="\t", header=T)


final2 <- read.table('/home/mferri/allchrts_Pd_high_low_krt.tsv', sep ="\t", header=T)
krt_h <- colnames(final2)
## upto here volcanos preparation 



m <- m[order(m$padj),]

p <- ggplot(data=m, aes(x=lfc, y=-log10(padj), color=vs.x))+geom_point(size=0.5, alpha=0.7)+theme_bw()+
  scale_color_manual(values=colors)
mm <- m[1:50,]
#mmm <- mm[mm$vs.x=="Bb",]
p + geom_text_repel(data=mm, aes(label=gene))

m2 <- m[m$padj < 0.1,]
p <- ggplot(data=m2, aes(x=lfc, y=-log10(padj), color=vs.x))+geom_point(size=0.5, alpha=0.7)+theme_bw()+
  scale_color_manual(values=colors)
mm <- m2[1:50,]
#mmm <- mm[mm$vs.x=="Bb",]
p + geom_text_repel(data=mm, aes(label=gene))

m3 <- m[m$`vs.x` %in% c('Ac','Ba','Bb'),]
m3$GO <- ifelse(m3$gene %in% go$SYMBOL, 'Keratinization', 'Other')
p <- ggplot(data=m3, aes(x=lfc, y=-log10(padj), color=vs.x, size=GO))+geom_point(alpha=0.7)+theme_bw()+
  scale_size_manual(values=c(1.5, 0.3))+
  scale_color_manual(values=colors)

mmm <- m3[ m3$gene %in% krt_h,]
p + geom_text_repel(data=mmm, aes(label=gene), size=4)

###################3 (non la facciamo perchè mettiamo il plot con tutti Ac etc)  fare questa per ogni confronto Aa vs almeno 5 campioni, occhio che il codice qui ha
# solo i padj di Bb quindi va cambiato
m3 <- m[m$`vs.x` %in% c('Bb'),]
m3$GO <- ifelse(m3$gene %in% go$SYMBOL, 'Keratinization', 'Other')
lfct <- 0.5849625
padjt <- 0.01
m3$col <- ifelse(abs(m3$lfc) > lfct & m3$padj < padjt, "both", ifelse(abs(m3$lfc) > lfct, "LFC",
                 ifelse(m3$padj < padjt, "padj", "NS")))
m3$col <- factor(m3$col, levels=c("LFC", "padj", "both", "NS"))
m3$lfc <- -m3$lfc
p <- ggplot(data=m3, aes(x=lfc, y=-log10(padj), color=col, size=GO))+geom_point(alpha=0.7)+theme_bw()+
  scale_size_manual(values=c(1.5, 0.3))+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE)

mmm <- m3[ m3$gene %in% krt_h,]
p + geom_text_repel(data=mmm, aes(label=gene), size=4, color="black")+ xlab('Log2 fold change')+ylab('-log10(padj)')+ggtitle('CTGC-Bb vs Aa')+ theme(plot.title = element_text(hjust = 0.5))

ggsave('/scratch/trcanmed/connector/Bb_Aa.svg', height=17, width=17, units="cm")

## Bb only label all GO
p <- ggplot(data=m3, aes(x=lfc, y=-log10(padj), color=col))+geom_point(alpha=0.7)+theme_bw()+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE)

mmm <- m3[ m3$gene %in% go$SYMBOL,]
p + geom_text_repel(data=mmm, aes(label=gene), size=4, color="black")


m$krt <- grepl('KRT', m$gene, fixed=F)


ggplot(data=m, aes(x=lfc, y=-log10(padj), color=vs.x, size=krt))+geom_point(alpha=0.7)+theme_bw()+scale_color_manual(values=colors$CONNECTOR)+scale_size_manual(values=c(0.3, 1.5))
final2 <- read.table('/home/mferri/allchrts_Pd_high_low_krt.tsv', sep ="\t", header=T)

krt_h <- colnames(final2)





mmm <- m[ m$gene %in% krt_h,]
p + geom_text_repel(data=mmm, aes(label=gene))

#### the WINNER CHOSEN ONE
tmpd <- read.table('/home/mviviani/misc_tmp/CONNECTOR/Jun2022/June2022_cluster-prepared.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)

nctgc <- unique(tmpd$cluster)
names(connector_col) <- nctgc[order(nctgc)]

vdkrt <- vd[vd$gene %in%  go$SYMBOL,]
vdkrt$lfc <- - vdkrt$lfc
#vdkrt <- vdkrt[vdkrt$vs %in% c('Ac', 'Ba', 'Bb', 'Bc', 'C'),]
vdkrt <- vdkrt[vdkrt$vs %in% c('Ac', 'Ba', 'Bb'),]

vdkrtBb <- vdkrt[vdkrt$vs == "Bb",]
vdkrtBb <- vdkrtBb[order(-vdkrtBb$lfc),]
krt_hh <- head(vdkrtBb$gene, n=5)

mmm <- vdkrt[vdkrt$gene %in% krt_hh,]

ggplot(data=vdkrt, aes(x=lfc, y=-log10(padj), color=vs))+geom_point(alpha=0.7, size=0.5)+
scale_color_manual(values=connector_col)+
geom_text_repel(data=mmm, aes(label=gene, color=vs), size=1.5, show.legend = FALSE)+xlab('Log2 FoldChange') + ylab('-Log10(padj)')+
ggtitle('Keratinization genes in CTGC DEG')+
guides(color=guide_legend(title="CTGC-* vs Aa"))+gtheme


ggsave('/scratch/trcanmed/connector/threeVolcano_leg.pdf', height=6, width=6, units="cm")

ggplot(data=vdkrt, aes(x=lfc, y=-log10(padj), color=vs))+geom_point(alpha=0.7, size=0.5)+
  scale_color_manual(values=connector_col)+
  geom_text_repel(data=mmm, aes(label=gene, color=vs), size=1.5, show.legend = FALSE)+xlab('Log2 FoldChange') + ylab('-Log10(padj)')+
  ggtitle('Keratinization genes in CTGC DEG')+
  gtheme+
  theme(legend.position="none")

ggsave('/scratch/trcanmed/connector/threeVolcano.pdf', height=6, width=6, units="cm")


### third panel
load('/scratch/trcanmed/DE_RNASeq/dataset/krt_hl_outcon/pippo.Rdata')
gene <- 'H_HOPX'
data <-plotCounts(dds, gene, intgroup='class', returnData=T)
data <- data[order(data$count),]
data <- data[data$class %in% c('H','L'),]
data$class <- ifelse(data$class == "H",'keratin-high', 'keratin-low')
data$class <- factor(data$class, levels=c('keratin-low', 'keratin-high'))
e2 <- ggplot(data, aes_string(x = what, y = "count"))
gene <- gsub('H_','', gene)
e3 <- e2 + geom_jitter(aes_string(color = what),   position = position_jitter(0.2),size = 0.5) +
  stat_summary( aes_string(color = what), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.1, color="darkgreen")+
  scale_y_continuous(trans='log10')+labs(color = "Group", x="Group", y="Log10(nreads)")+ggtitle("HOPX expression")+
  gtheme

e3
ggsave('/scratch/trcanmed/connector/hopx_legend.pdf', height=6, width=6, units="cm")

e3 <- e2 + geom_jitter(aes_string(color = what),   position = position_jitter(0.2),size = 0.5) +
  stat_summary( aes_string(color = what), fun.data="mean_sdl",  fun.args = list(mult=1),  geom = "pointrange",  size = 0.1, color="darkgreen")+
  scale_y_continuous(trans='log10')+labs(color = "Group", x="Group", y="Log10(nreads)")+ggtitle("HOPX expression")+
  gtheme+theme(legend.position="none")

e3
ggsave('/scratch/trcanmed/connector/hopx.pdf', height=6, width=6, units="cm")

### B PD?
hl <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/krt_hl_outcon/samples_data', sep="\t", header=T)
# not managing replicates
m <- merge(data, hl, by.x="row.names", by.y="model")
