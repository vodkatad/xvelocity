bayes <- readRDS('/scratch/trcanmed/connector/local/share/data/second_round_bayes.rds')
table(bayes$col)
#cris <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
#data <- merge(data, cris, by.x="model", by.y="genealogy", all.x=TRUE)

cetuxi <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/november2021/recist_3w.tsv', sep="\t", stringsAsFactors = F, header=T)
data <- merge(bayes, cetuxi, by.x="ShortID", by.y="smodel", all.x=TRUE)

df <- as.data.frame(table(data$col, data$recist))
colnames(df) <- c('connector', 'recist', 'n')
df$recist <- factor(df$recist, levels=c('OR', 'SD', 'PD'))
ggplot(data =  df, mapping = aes(x = recist, y = connector)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() + theme(legend.position = "none")


##

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
