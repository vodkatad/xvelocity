#!/usr/bin/env Rscript

set.seed(42)
library(getopt)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'cet_in', 'c', 1, 'character',
  'irino_in', 'i', 1, 'character',
  'one_out', 'o', 1, 'character',
  'trad_out', 't', 1, 'character',
  'dic_out', 'd', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

cet <- read.table(gzfile(opt$cet_in), sep='\t', header=TRUE)
# cet$model <- substr(cet$longen,1,10)
cet$drug <- "cetuxi"
### model
# cetuniq <- cet[!duplicated(cet$model),]
### longen
cetuniq <- cet[!duplicated(cet$longen),]

iri <- read.table(gzfile(opt$irino_in), sep='\t', header=TRUE)
iri$drug <- "irino"
# iri$model <- substr(iri$longen,1,10)
### model
# iriuniq <- iri[!duplicated(iri$model),]
### longen
iriuniq <- iri[!duplicated(iri$longen),]

one <- rbind(cet,iri)
write.table(one, gzfile(opt$one_out), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

one_uniq <- rbind(cetuniq, iriuniq)

### model
# md <- one_uniq[,c("model","drug")]
# md <- md[order(md$model),]
### longen
md <- one_uniq[,c("longen","drug")]
md <- md[order(md$longen),]
md$trad <- paste0("CRC",seq(0000, by=1, length.out = nrow(md)))
write.table(md, opt$dic_out, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

### model
# final <- merge(one, md, by='model', all.x=TRUE)
### longen
final <- merge(one, md, by='longen', all.x=TRUE)

tothem <- final[,c("trad","start_date","end_date","measure_date","volume")]
names(tothem)[1] <- "model"
write.table(tothem, gzfile(opt$trad_out), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
