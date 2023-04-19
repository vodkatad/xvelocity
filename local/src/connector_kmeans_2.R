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
pd$error_code <- as.factor(ifelse(pd$error == "phy_worse", "Worse", ifelse(pd$error == "phy_better", 'Better', '')))
pd$just <- ifelse(pd$perc_3w > 0, -0.5, 1.1)
pd$pospoint <- ifelse(pd$perc_3w > 0, pd$perc_3w+10, pd$perc_3w-10)

levels(pd$Cluster) <- c('1_growing','2_shrinking','3_stable')

ggplot(data=pd, aes(x=reorder(ShortID,-perc_3w), y=perc_3w, fill=Cluster, color=error))+
geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+labs(fill='Connector cluster')+
scale_color_manual(values=c('blue','red','white'))+scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
geom_hline(yintercept=-50)+geom_hline(yintercept=35)

ggplot(data=pd, aes(x=reorder(ShortID,-perc_3w), y=perc_3w, fill=Cluster))+
  geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+labs(fill='Confusion matrix')+
  scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
  geom_hline(yintercept=-50)+geom_hline(yintercept=35)+
  geom_text(aes(label = error_code, color=error_code, vjust=just), 
            size = 5,  face = "bold")+
  scale_color_manual(values=c('white','darkblue','darkred'))

ggplot(data=pd, aes(x=reorder(ShortID,-perc_3w), y=perc_3w, fill=Cluster))+
  geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+
  scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
  geom_hline(yintercept=-50)+geom_hline(yintercept=35)+
  geom_point(aes(y=pospoint, shape = error_code, color=error_code), 
            size = 2)+
  scale_color_manual(values=c('white','darkblue','darkred'), guide="none")+scale_shape_manual(values=c(46,25,24))+
  labs(fill='Connector cluster', shape="Recist Classification", color=NULL)

# 6w
pd2 <- pd[!is.na(pd$perc_6w),]
pd2$pospoint <- ifelse(pd2$perc_6w > 0, pd2$perc_6w+10, pd2$perc_6w-10)


ggplot(data=pd2, aes(x=reorder(ShortID,-perc_6w), y=perc_6w, fill=Cluster, color=error))+
  geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+labs(fill='Confusion matrix')+
  scale_color_manual(values=c('blue','red','white'))+scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
  geom_hline(yintercept=-50)+geom_hline(yintercept=35)



ggplot(data=pd2, aes(x=reorder(ShortID,-perc_6w), y=perc_6w, fill=Cluster))+
  geom_col()+current_theme+xlab('models')+ylab('3 weeks Cetuximab response')+
  scale_fill_manual(values=c('#ff6666','#009933','#d6d6c2'))+
  geom_hline(yintercept=-50)+geom_hline(yintercept=35)+
  geom_point(aes(y=pospoint, shape = error_code, color=error_code), 
             size = 2)+
  scale_color_manual(values=c('white','darkblue','darkred'), guide="none")+scale_shape_manual(values=c(46,25,24))+
  labs(fill='Connector cluster', shape="Recist Classification", color=NULL)

## exp g
data <- read.table(gzfile('~/velocity/cetuxi6w_long_header.tsv.gz'), sep="\t", header=TRUE, stringsAsFactors = FALSE)
data$smodel <- substr(data$longen, 0, 7)

expg <- as.data.frame(sapply(unique(data$smodel), function(x) { length(unique(data[data$smodel ==x,'exp_group'])) }))
colnames(expg) <- 'n_exp_g'

mexpg <- merge(pd, expg, by.x="ShortID", by.y="row.names")
nrow(pd) == nrow(mexpg)

mexpg$error <- as.factor(ifelse(mexpg$confmatrix %in% c('2+PD','2+SD', '3+PD'), 'phy_worse', ifelse(mexpg$confmatrix %in% c('1+OR','1+SD', '3+OR'), 'phy_better', 'right')))

ggplot(data=mexpg, aes(x=confmatrix, y=n_exp_g, color=error))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+current_theme+theme(axis.text.x = element_text(size = rel(textSize), angle = 90, vjust = 0.5, hjust=1))+ scale_color_manual(values=c('blue','red','black'))

data$start_date_d <- as.Date(data$start_date)
data$end_date_d <- as.Date(data$end_date)
data$measure_date_d <- as.Date(data$measure_date)


sd_mice <- function(expg, data) {
  # we get the sd between same measure date
  # TODO limit to dates between start/end dat
  d <- data[data$exp_group == expg,]
  d <- d[d$measure_date_d+1 >= d$start_date_d & d$measure_date_d <= d$end_date_d,]
  if (nrow(d) > 1) {
    sds <- sapply(unique(d$measure_date), function(x) { sd(d[d$measure_date==x,'volume']) } )
    res <- mean(sds, na.rm=TRUE)
  #names(res) <- unique(d$smodel) 
  } else {
    res <- NA
  }
  return(res)
}

sd_expg <- sapply(unique(data$exp_group), sd_mice, data)
sd_smodel <- sapply(unique(data$exp_group), function(x) { unique(data[data$exp_group==x, 'smodel'])} )

sd_df <- data.frame(expg=names(sd_expg), sd=sd_expg, smodel=sd_smodel)


mean_sd_smodel <- as.data.frame(sapply(unique(sd_df$smodel), function(x) { mean(sd_df[sd_df$smodel==x,'sd'], na.rm=TRUE)}))
colnames(mean_sd_smodel) <- 'mean_sd_expg'
rownames(mean_sd_smodel) <- unique(sd_df$smodel)


msd <- merge(pd, mean_sd_smodel, by.x="ShortID", by.y="row.names")
nrow(pd) == nrow(msd)

msd$error <- as.factor(ifelse(msd$confmatrix %in% c('2+PD','2+SD', '3+PD'), 'phy_worse', ifelse(msd$confmatrix %in% c('1+OR','1+SD', '3+OR'), 'phy_better', 'right')))

ggplot(data=msd, aes(x=confmatrix, y=mean_sd_expg, color=error))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+current_theme+theme(axis.text.x = element_text(size = rel(textSize), angle = 90, vjust = 0.5, hjust=1))+ scale_color_manual(values=c('blue','red','black'))


### n topi

n_mice <- function(expg, data) {
  # we get the sd between same measure date
  # TODO limit to dates between start/end dat
  d <- data[data$exp_group == expg,]
  d <- d[d$measure_date_d+1 >= d$start_date_d & d$measure_date_d <= d$end_date_d,]
  return(length(unique(d$longen)))
}

n_expg <- sapply(unique(data$exp_group), n_mice, data)
n_smodel <- sapply(unique(data$exp_group), function(x) { unique(data[data$exp_group==x, 'smodel'])} )

n_df <- data.frame(expg=names(n_expg), n=n_expg, smodel=sd_smodel)


mean_n_smodel <- as.data.frame(sapply(unique(n_df$smodel), function(x) { mean(n_df[n_df$smodel==x,'n'], na.rm=TRUE)}))
colnames(mean_n_smodel) <- 'mean_n_expg'
rownames(mean_n_smodel) <- unique(n_df$smodel)


msd <- merge(pd, mean_n_smodel, by.x="ShortID", by.y="row.names")
nrow(pd) == nrow(msd)

msd$error <- as.factor(ifelse(msd$confmatrix %in% c('2+PD','2+SD', '3+PD'), 'phy_worse', ifelse(msd$confmatrix %in% c('1+OR','1+SD', '3+OR'), 'phy_better', 'right')))

ggplot(data=msd, aes(x=confmatrix, y=mean_n_expg, color=error))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+current_theme+theme(axis.text.x = element_text(size = rel(textSize), angle = 90, vjust = 0.5, hjust=1))+ scale_color_manual(values=c('blue','red','black'))



### n misure

n_mice <- function(expg, data) {
  # we get the sd between same measure date
  # TODO limit to dates between start/end dat
  d <- data[data$exp_group == expg,]
  d <- d[d$measure_date_d+1 >= d$start_date_d & d$measure_date_d <= d$end_date_d,]
  return(length(unique(d$measure_date_d)))
}

n_expg <- sapply(unique(data$exp_group), n_mice, data)
n_smodel <- sapply(unique(data$exp_group), function(x) { unique(data[data$exp_group==x, 'smodel'])} )

n_df <- data.frame(expg=names(n_expg), n=n_expg, smodel=sd_smodel)


mean_n_smodel <- as.data.frame(sapply(unique(n_df$smodel), function(x) { mean(n_df[n_df$smodel==x,'n'], na.rm=TRUE)}))
colnames(mean_n_smodel) <- 'mean_n_expg'
rownames(mean_n_smodel) <- unique(n_df$smodel)


msd <- merge(pd, mean_n_smodel, by.x="ShortID", by.y="row.names")
nrow(pd) == nrow(msd)

msd$error <- as.factor(ifelse(msd$confmatrix %in% c('2+PD','2+SD', '3+PD'), 'phy_worse', ifelse(msd$confmatrix %in% c('1+OR','1+SD', '3+OR'), 'phy_better', 'right')))

ggplot(data=msd, aes(x=confmatrix, y=mean_n_expg, color=error))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+current_theme+theme(axis.text.x = element_text(size = rel(textSize), angle = 90, vjust = 0.5, hjust=1))+ scale_color_manual(values=c('blue','red','black'))

