library(GEOquery, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(survival)
setwd("~")
load("GSE5851.RData")
gset <- gse[[1]]


#entrez <- c(3853)
#names(entrez) <- c("KRT6A")
entrez <- 84525
names(entrez) <- c("HOPX")

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  exprlist <- lapply(entrez, function(x) {log(exprs(gset[which(pData(featureData(gset))[,"ENTREZ_GENE_ID"]==x)]))/log(2)})
} else {
  exprlist <- lapply(entrez, function(x) {exprs(gset[which(pData(featureData(gset))[,"ENTREZ_GENE_ID"]==x)])})
}

geoSamples <- pData(phenoData(gset))[,c("title","geo_accession")]
#metagenedf <- cbind(expr=exprlist$KRT6A[1,], geoSamples) # TODO check correspondence of clinical data with excel!
metagenedf <- cbind(expr=exprlist$HOPX[1,], geoSamples) # TODO check correspondence of clinical data with excel!

responsesOrig <- read.table("Khambata-Ford.txt", sep="\t", header=T)
x <- gsub("FL_","_", responsesOrig[,"Affymetrix.ID"] )
responsesOrig$Affymetrix.ID <- x
responses <- data.frame(response=responsesOrig$Best.Clinical.Response.Assessment, ids=x)
responsesPruned <- responses[responses$ids != "",]

pfs <- data.frame(response=responsesOrig$Progression.free.survival....of.days, ids=responsesOrig$Affymetrix.ID)
pfsPruned <- pfs[pfs$ids != "",]
pfs_data <- merge(metagenedf, pfs, by.x="title", by.y="ids")
pfs_data$event <- rep(1, nrow(pfs_data))
x <- pfs_data
median <- median(x$expr)
x$group <- as.factor(ifelse(x$expr > median, "High HOPX", "Low HOPX"))
surv <- Surv(time=x$response, event=x$event)
sd <-survdiff(surv~x$group)
pvalue <- 1 - pchisq(sd$chisq , length(sd$n) - 1)
plot(survfit(surv ~ x$group), col=c("red","blue"), lwd=2, lty=1:1, xlab = "time",
     ylab = "S(t)", main=paste("HOPX median"), sub=paste("pvalue:",pvalue))
legend("bottomright", levels(x$group), col=c("red","blue"), lty=1:1, lwd=2)



pfs <- data.frame(response=responsesOrig$Best.Clinical.Response.Assessment, ids=responsesOrig$Affymetrix.ID)
pfsPruned <- pfs[pfs$ids != "",]

data <- merge(metagenedf, pfs, by.x="title", by.y="ids")
data$response <- factor(data$response, levels=c("CR","PR","SD","PD","UTD"))
ggplot(data, aes(x = response, y=expr, fill=response)) + geom_boxplot() +geom_jitter()

pd <- data[data$response=="PD",]
x <- x[x$title %in% pd$title,]
median <- median(x$expr)
x$group <- as.factor(ifelse(x$expr > median, "High HOPX", "Low HOPX"))
surv <- Surv(time=x$response, event=x$event)
sd <-survdiff(surv~x$group)
pvalue <- 1 - pchisq(sd$chisq , length(sd$n) - 1)
plot(survfit(surv ~ x$group), col=c("red","blue"), lwd=2, lty=1:1, xlab = "time",
     ylab = "S(t)", main=paste("HOPX median"), sub=paste("pvalue:",pvalue))
legend("bottomright", levels(x$group), col=c("red","blue"), lty=1:1, lwd=2)


# last resort HOPX+K6

entrez <- c(84525, 3853)
names(entrez) <- c("HOPX", "KRT6A")

metagene <- sapply(seq(1,80), function(x) { median(log(exprs(gset[which(pData(featureData(gset))[,"ENTREZ_GENE_ID"] %in% entrez)])[,x])/log(2))})
metagenedf <- data.frame(expr=metagene)
geoSamples <- pData(phenoData(gset))[,c("title","geo_accession")]
metagenedf <- cbind(metagenedf, geoSamples, stringsAsFactors=F) # TODO check correspondence of clinical data with excel!


responsesOrig <- read.table("Khambata-Ford.txt", sep="\t", header=T)
x <- gsub("FL_","_", responsesOrig[,"Affymetrix.ID"] )
responsesOrig$Affymetrix.ID <- x
responses <- data.frame(response=responsesOrig$Best.Clinical.Response.Assessment, ids=x)
responsesPruned <- responses[responses$ids != "",]

pfs <- data.frame(response=responsesOrig$Progression.free.survival....of.days, ids=responsesOrig$Affymetrix.ID)
pfsPruned <- pfs[pfs$ids != "",]

pfs_data <- merge(metagenedf, pfsPruned, by.x="title", by.y="ids")
pfs_data$event <- rep(1, nrow(pfs_data))
x <- pfs_data
median <- median(x$expr)
x$group <- as.factor(ifelse(x$expr > median, "High", "Low"))
surv <- Surv(time=x$response, event=x$event)
sd <-survdiff(surv~x$group)
pvalue <- 1 - pchisq(sd$chisq , length(sd$n) - 1)
plot(survfit(surv ~ x$group), col=c("red","blue"), lwd=2, lty=1:1, xlab = "time",
     ylab = "S(t)", main=paste("median"), sub=paste("pvalue:",pvalue))
legend("bottomright", levels(x$group), col=c("red","blue"), lty=1:1, lwd=2)

