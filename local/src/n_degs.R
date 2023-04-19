library(ggplot2)
n_deg <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/connector_for_all/n_deg.txt' , sep = "\t", header=F, stringsAsFactors = F)
n_deg_krt <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/connector_for_all/n_deg_krt.txt' , sep = "\t", header=F, stringsAsFactors = F)

colnames(n_deg) <- c('direction', 'comparison', 'n')
colnames(n_deg_krt) <- c('kdirection', 'kcomparison', 'kn')

n_deg$v <- sapply(strsplit(n_deg$comparison, '/'), function(x) {x[1]})
n_deg$vs <- sapply(strsplit(n_deg$v, '_'), function(x) {x[2]})


n_deg_krt$kv <- sapply(strsplit(n_deg_krt$kcomparison, '/'), function(x) {x[1]})
n_deg_krt$kvs <- sapply(strsplit(n_deg_krt$kv, '_'), function(x) {x[2]})


deg <- cbind(n_deg, n_deg_krt)
deg$nKRT <- deg$n - deg$kn

pdata <- data.frame(vs=c(deg$vs, deg$kvs), direction=c(deg$direction, deg$kdirection), 
                    n_deg=c(deg$nKRT, deg$kn), class=c(rep('deg', nrow(deg)), rep('keratinization', nrow(deg))))

pdata$ln_deg <- log(pdata$n_deg+1)

ggplot(data=pdata[pdata$direction=="UP",], aes(x=vs, y=ln_deg,fill=class))+geom_col()+theme_bw()+
  geom_col(data=pdata[pdata$direction=="DOWN",], aes(x=vs, y=-ln_deg,fill=class))+scale_fill_manual(values=c('goldenrod', ' forestgreen'))
  

ggplot(data=pdata[pdata$direction=="UP",], aes(x=vs, y=n_deg,fill=class))+geom_col()+theme_bw()+
  geom_col(data=pdata[pdata$direction=="DOWN",], aes(x=vs, y=-n_deg,fill=class))+scale_fill_manual(values=c('goldenrod', ' forestgreen'))


load('/scratch/trcanmed/DE_RNASeq/dataset/connector_for_all/Aa_Bb/GO.Rdata')

