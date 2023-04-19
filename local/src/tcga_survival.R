# downloaded z score RNAseq vs all samples and clinical data on 27th/07/22, 
#https://www.cbioportal.org/study/clinicalData?id=coadread_tcga
#https://www.cbioportal.org/results/download?cancer_study_list=coadread_tcga&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=rna_seq_v2_mrna_median_all_sample_Zscores&case_set_id=coadread_tcga_all&gene_list=HOPX%250AKRT6A%250AKRT6B%250AKRT17%250AKRT16&geneset_list=%20&tab_index=tab_visualize&Action=Submit&comparison_subtab=clinical
# firehose legacy
library(survival)
library(survminer)
library(ggplot2)

setwd('/scratch/trcanmed/connector/local/share/data/tcga')
expr <- read.table('expr.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
clin <- read.table('coadread_tcga_clinical_data.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)


dim(expr)
length(substr(expr$SAMPLE_ID, 1, 12))
length(unique(substr(expr$SAMPLE_ID, 1, 12)))

# we remove the patients with duplicated samples, they are few (4)
# alternative: average expr of the duplicated samples
dup <- as.data.frame(table(substr(expr$SAMPLE_ID, 1, 12)))

expr$patient_id <- substr(expr$SAMPLE_ID, 1, 12)
expr_nodup <- expr[expr$patient_id %in% dup[dup$Freq == 1,'Var1'],]

expr_nodup <- expr_nodup[apply(expr_nodup[, c(3,4,5,6,7)], 1, function(x) {all(x!="NP")}),]

for (i in c(3,4,5,6,7)) {
  expr_nodup[,i] <- as.numeric(expr_nodup[,i])
}

clin <- clin[, c("Patient.ID", 'Disease.Free..Months.','Disease.Free.Status', 'Overall.Survival..Months.', 'Overall.Survival.Status')]


mdata <- merge(clin, expr_nodup, by.x="Patient.ID", by.y="patient_id")

mdata2 <- mdata[!is.na(mdata$Overall.Survival.Status),]
table(sapply(strsplit(mdata2$Overall.Survival.Status, ":"), `[[`, 2))

mdata2 <- mdata[!is.na(mdata$Disease.Free.Status),]
table(sapply(strsplit(mdata2$Disease.Free.Status, ":"), `[[`, 2))

median_all <- apply(mdata[, c('HOPX', 'KRT6A' ,'KRT6B')], 1, max)

data <- data.frame(row.names=mdata$Patient.ID, 
                   expr=mdata$HOPX, 
                   #expr=median_all,
                   event=as.numeric(sapply(strsplit(mdata$Overall.Survival.Status, ":"), `[[`, 1)),
                   time=mdata$Overall.Survival..Months.
                   #event=as.numeric(sapply(strsplit(mdata$Disease.Free.Status, ":"), `[[`, 1)),  # works for HOPX 0.035
                   #time=mdata$Disease.Free..Months.
                  )

table(is.na(data$event), is.na(data$time))

x <- data[!is.na(data$event) & !is.na(data$time),]

#median <- quantile(x$expr, 0.75)
#median <- median(x$expr)
cutp <- surv_cutpoint(
  x,
  time = "time",
  event = "event",
  variables = "expr",
  minprop = 0.1,
  progressbar = TRUE
)
median <- as.numeric(cutp$cutpoint[1])

plot(cutp, "expr", palette = "npg")

x$group <- as.factor(ifelse(x$expr > median, "High", "Low"))
surv <- Surv(time=x$time, event=x$event)
sd <-survdiff(surv~x$group)
pvalue <- 1 - pchisq(sd$chisq , length(sd$n) - 1)
plot(survfit(surv ~ x$group), col=c("red","blue"), lwd=2, lty=1:1, xlab = "time",
    ylab = "S(t)", main=paste("HOPX median"), sub=paste("pvalue:",pvalue))
legend("bottomright", levels(x$group), col=c("red","blue"), lty=1:1, lwd=2)

fit <- survfit(surv ~ x$group, data = x)

ggsurvplot(
  fit,
  data = x,
  size = 1,                 # change line size
  palette =
    c("#E69F00", "#56B4E9"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

## corrected pvalue

res.cat <- surv_categorize(cutp)
head(res.cat)

# 4. Fit survival curves and visualize

fit2 <- survfit(surv ~ expr, data = res.cat)

ggsurvplot(
  fit2,
  data = res.cat,
  size = 1,                 # change line size
  palette =
    c("#E69F00", "#56B4E9"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

library(maxstat)

ms <- maxstat.test(surv ~ expr, data = x, smethod="LogRank",  pmethod="condMC", minprop=0.1, maxprop=0.9)

median == ms$estimate

ms$p.value


ggsurvplot(
  fit2,
  data = res.cat,
  size = 1,                 # change line size
  palette =
    c("#E69F00", "#56B4E9"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = paste('pval=', as.character(ms$p.value)),              # Add p-value
  risk.table = FALSE,
  legend.labs =
    c("High HOPX", "Low HOPX"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = gtheme # Change ggplot2 theme,
)

ggsave('/scratch/trcanmed/connector/hopx_surv_leg.pdf', height=6, width=6, units="cm")

ggsurvplot(
  fit2,
  data = res.cat,
  size = 1,                 # change line size
  palette =
    c("#E69F00", "#56B4E9"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = paste('pval=', as.character(ms$p.value)),              # Add p-value
  risk.table = FALSE,
  legend="none", # Add risk table
  legend.labs =
    c("High HOPX", "Low HOPX"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = gtheme # Change ggplot2 theme
)

ggsave('/scratch/trcanmed/connector/hopx_surv.pdf', height=6, width=6, units="cm")

