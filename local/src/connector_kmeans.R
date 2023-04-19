library(ggplot2)
textSize <- 2
current_theme <-
  theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    axis.title = element_text(size = rel(textSize)),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0,
                               size = rel(textSize)),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(size = rel(textSize), face = "italic"),
    legend.text = element_text(size = rel(textSize)),
    plot.title = element_text(
      face = "bold",
      size = rel(textSize),
      hjust = 0.5
    )
  )



load('/scratch/trcanmed/stash/PDXCLusterGenomicPhysician.RData')

cet3 <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t")
colnames(cet3) <- c('smodel', 'perc_3w')

cet6 <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi6w_CRC0078_PR.tsv', sep="\t")
colnames(cet6) <- c('smodel', 'perc_6w')

m1 <- merge(GenCet, cet3, all.x=TRUE, by.x="ShortID", by.y="smodel")
m2 <- merge(m1, cet6, all.x=TRUE, by.x="ShortID", by.y="smodel")

# chiedere ad Eugy di 11 senza dati risposta cetuxi:
# m2[is.na(m2$perc_3w), ]
# per lo piu` modelli 'alti'

# we remove models without a cluster CRC0456, CRC2407 (?)
m3 <- m2[!is.na(m2$Cluster),]
# and those without a cetuxi response
pd <- m3[!is.na(m3$perc_3w),]

pd$confmatrix <- as.factor(paste0(pd$Cluster, "+", pd$PhyCluster))

# 1 is PD
# 2 is OR
# 3 is SD
pd$error <- as.factor(ifelse(pd$confmatrix %in% c('2+PD','2+SD', '3+PD'), 'phy_worse', ifelse(pd$confmatrix %in% c('1+OR','1+SD', '3+OR'), 'phy_better', 'right')))
pd$Cluster <- as.factor(pd$Cluster)
levels(pd$Cluster) <- c('1_PD','2_OR','3_SD')
ggplot(data=pd, aes(x=reorder(ShortID,-perc_3w), y=perc_3w, fill=Cluster, color=error))+
geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+labs(fill='Confusion matrix')+
scale_color_manual(values=c('blue','red','white'))+scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
geom_hline(yintercept=-50)+geom_hline(yintercept=35)

# TODO 6w
