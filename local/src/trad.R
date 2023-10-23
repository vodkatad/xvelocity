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
cet$model <- substr(cet$longen,1,7)
cet$type <- substr(cet$longen,1,10)
cet$drug <- "cetuxi"
### --- ###
iri <- read.table(gzfile(opt$irino_in), sep='\t', header=TRUE)
iri$model <- substr(iri$longen,1,7)
iri$type <- substr(iri$longen,1,10)
iri$drug <- "irino"

one <- rbind(cet,iri)
write.table(one, gzfile(opt$one_out), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

new_cet <- NULL
### 1169 & 3116 have both LMX and PRX, so we need to manage this to not loose one during the unique opn the model.
mods <- unique(cet$model)
for( i in 1:length(mods) ) {
  tmp <- cet[cet$model==mods[i],]
  if( any(grepl("PRX", tmp$type)) ) {
    tmp <- tmp[!duplicated(tmp$type),]
  } else {
    tmp <- tmp[!duplicated(tmp$model),]
  }
  new_cet <- rbind(new_cet,tmp)
}

new_iri <- NULL
mods <- unique(iri$model)
for( i in 1:length(mods) ) {
  tmp <- iri[iri$model==mods[i],]
  if( any(grepl("PRX", tmp$type)) ) {
    tmp <- tmp[!duplicated(tmp$type),]
  } else {
    tmp <- tmp[!duplicated(tmp$model),]
  }
  new_iri <- rbind(new_iri,tmp)
}

new <- rbind(new_cet, new_iri)

md <- new[,c("model","drug","type")]
md <- md[order(md$model),]
md$trad <- paste0("CRC",seq(1000, by=10, length.out = nrow(md)))
out  <- md[,c("model","drug","trad")]
write.table(out, opt$dic_out, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

### 196, 515 & 3210 have both cetuxi and chemio treat, so we need to manage this as well,
### so I merge each df alone and bind them together after
cetcet <- merge(md[md$drug=='cetuxi',], cet, by='type')
cetcet <- cetcet[,c(5,6,7,8,9,10,11,12,4)]
names(cetcet)[c(7,8)] <- c("model", "drug")
iriri <- merge(md[md$drug=='irino',], iri, by='type')
iriri <- iriri[,c(5,6,7,8,9,10,11,12,4)]
names(iriri)[c(7,8)] <- c("model", "drug")

final <- rbind(cetcet, iriri)
tothem <- final[,c("trad","start_date","end_date","measure_date","volume")]
names(tothem)[1] <- "model"
write.table(tothem, gzfile(opt$trad_out), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
