d <- read.table("~/work/stash/velocity/CRC_phys_standard_cetux_long_expg.tsv", sep="\t", header=F)
colnames(d) <- c('gen','expg','date','vol')

# Filtriamo via i gruppi sperimentali con numero di topi/misure estremi
nums <- as.data.frame(table(d$expg))

q <- quantile(nums$Freq, probs=c(0.01, 0.99))
nums$Var1 <- as.character(nums$Var1)
remove <- c(nums[nums$Freq <= q[1],'Var1'], nums[nums$Freq >= q[2],'Var1'])
dd <- d[! d$expg %in% remove,]

topos <- sapply(unique(dd$expg), function(x) {y <- dd[dd$expg==x,]; length(unique(y$gen))})
#table(topos)
#topos[topos == 3]
ddd <- dd[dd$expg !="ORGANOIDI_CRC0177LMX0D.2019-05-08.Cetuximab Standard.0",]

ddd$realdate <- as.Date(as.character(ddd$date), format="%Y-%m-%d")

# Blocco fit lineari

fit_all_mice <-  function(expg_d) {
  expg_d <- expg_d[order(expg_d$realdate),]
  vels <- sapply(unique(expg_d$gen), function(x) fit_one_mouse(expg_d[expg_d$gen==x,]) )
  res <- t(vels)
  colnames(res) <- c('slope','pval','R2', 'n_measures','logLik')
  res <- as.data.frame(res)
  res$expg <- unique(expg_d$expg)
  res
  }

lmp <- function (sm) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- sm$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

fit_one_mouse <- function(mouse_d) {
  res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'vol'] <- mouse_d[1,'vol']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'vol'] <- mouse_d[i, 'vol']
    res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
  }
  model <- lm(data=res, formula="vol~day")
  sm <- summary(model)
  r2 <- sm$r.squared
  pval <- lmp(sm) 
  ll <- logLik(model)
  #ggplot(res,aes(day,vol))+geom_point()+geom_smooth(method='lm')
  #ggsave(paste0(unique(mouse_d$gen),".png"))
  # do we get only some R2/pvals?
  vel <- model$coefficients[2]
  attributes(vel) <- NULL
  return(c(vel, pval, r2, nrow(mouse_d),ll))
}


all_vels <- lapply(unique(ddd$expg), function(x) fit_all_mice(ddd[ddd$expg==x,]) )
dfvel <- do.call(rbind, all_vels)
dfvel <- dfvel[!is.na(dfvel$pval),]

dfvel$padj <- p.adjust(dfvel$pval)
dfvel_sign <- dfvel[dfvel$padj < 0.05,]
dfvel_fil <- dfvel[dfvel$pval<0.05 & dfvel$R2 > 0.7,]

summarize_expg <-  function(dfvel_g) {
  avg <- mean(dfvel_g$slope)
  sd <- sd(dfvel_g$slope)
  n <- nrow(dfvel_g)
  return(c(avg, sd, n))
}

get_mice <- function(dfvel_g) {
  unique(substr(rownames(dfvel_g),0,14))
}

# gp!
all_avg <- sapply(unique(dfvel_fil$expg), function(x) summarize_expg(dfvel_fil[dfvel_fil$expg==x,]) )

avg <- t(all_avg)
colnames(avg) <- c('mean','sd','nmice')

m <- sapply(unique(dfvel_fil$expg), function(x) get_mice(dfvel_fil[dfvel_fil$expg==x,]) )
avg <- as.data.frame(avg)
avg$mouse <- m
write.table(avg, file="average_slope.tsv", sep="\t", quote=F)
write.table(dfvel, file="all_models.tsv", sep="\t", quote=F)

dd <- as.data.frame(table(substr(avg$mouse, 0,10)))
re <- dd[dd$Freq>1,'Var1']
avg$m <- substr(avg$mouse, 0, 10)
stable <- avg[avg$m %in% as.character(re),]


ggplot(data=stable, aes(x=m, y=mean,))+geom_point()+coord_flip()

avgnn <- avg[!is.na(avg$sd),]
avgnn$lower <- avgnn$mean - avgnn$sd
avgnn$upper <- avgnn$mean + avgnn$sd


ggplot(avgnn, aes(x=reorder(mouse,mean), y=mean)) +
geom_point(position=position_dodge(1), stat="identity") +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(1))+theme_bw()+ggtitle('avg linear fit slope vol ~ day')+theme(axis.text.x=element_blank())


### (vmax - vmin)/ (tmax-tmin)
# 
# notfit_all_mice <-  function(expg_d) {
#   expg_d <- expg_d[order(expg_d$realdate),]
#   vels <- sapply(unique(expg_d$gen), function(x) notfit_one_mouse(expg_d[expg_d$gen==x,]) )
#   return(vels)
# }
# 
# 
# notfit_one_mouse <- function(mouse_d) {
#   res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
#   res[1,'vol'] <- mouse_d[1,'vol']
#   for (i in seq(2, nrow(mouse_d))) {
#     res[i,'vol'] <- mouse_d[i, 'vol']
#     res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
#   }
#   vmax <- max(res$vol)
#   vmin <- min(res$vol)
#   imax <- match(vmax, res$vol)
#   imin <- match(vmin, res$vol)
#   tmax <- res[imax,'day']
#   tmin <- res[imin,'day']
#   return((vmax-vmin)/(tmax-tmin))
# }
# 
# all_vels_nofit <- lapply(unique(ddd$expg), function(x) notfit_all_mice(ddd[ddd$expg==x,]) )
# dfvel_nofit <- as.data.frame(unlist(all_vels_nofit))
# m <- merge(dfvel_nofit, dfvel, by="row.names")
# cor.test(m$slope, m$`unlist(all_vels_nofit)`)
# 
# colnames(m)[2] <- 'vmaxvmin'
# colnames(m)[3] <- 'slopefit'
# ggplot(m,aes(vmaxvmin,slopefit))+geom_point()+geom_smooth(method='lm')

### exp fit with log transf of y
#https://stackoverflow.com/questions/31851936/exponential-curve-fitting-in-r

efit_all_mice <-  function(expg_d) {
  expg_d <- expg_d[order(expg_d$realdate),]
  vels <- sapply(unique(expg_d$gen), function(x) efit_one_mouse(expg_d[expg_d$gen==x,]) )
  res <- t(vels)
  colnames(res) <- c('slope','pval','R2', 'n_measures', 'logLik')
  res <- as.data.frame(res)
  res$expg <- unique(expg_d$expg)
  res
}


efit_one_mouse <- function(mouse_d) {
  res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'vol'] <- mouse_d[1,'vol']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'vol'] <- mouse_d[i, 'vol']
    res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
  }
  res$vol <- log(res$vol)
  model <- lm(data=res, formula="vol~day")
  sm <- summary(model)
  r2 <- sm$r.squared
  pval <- lmp(sm) 
  ll <- logLik(model)
  #ggplot(res,aes(day,vol))+geom_point()+geom_smooth(method='lm')
  #ggsave(paste0(unique(mouse_d$gen),".png"))
  # do we get only some R2/pvals?
  vel <- model$coefficients[2]
  attributes(vel) <- NULL
  return(c(vel, pval, r2, nrow(mouse_d),ll))
}

all_expvels <- lapply(unique(ddd$expg), function(x) efit_all_mice(ddd[ddd$expg==x,]) )
dfvel_exp <- do.call(rbind, all_expvels)

cor.test(m2$slope.y, m2$slope.y)
dfvel_exp$kind <- 'exp'
dfvel$kind <- 'linear'
m2 <- merge(dfvel_exp, dfvel, by="row.names")
cor.test(m2$slope.y, m2$slope.y)
ggplot(m2,aes(slope.x,slope.y))+geom_point()+geom_smooth(method='lm')
m3 <- m2[,c('slope.x','pval.x','R2.x','kind.x', 'logLik.x')]
m4 <- m2[,c('slope.y','pval.y','R2.y','kind.y', 'logLik.y')]
colnames(m3) <- c('slope','pval','R2','kind', 'logLik')
colnames(m4) <- c('slope','pval','R2','kind', 'logLik')
m5 <- rbind(m3,m4)
ggplot(m5, aes(x=R2, color=kind))+geom_density()
ggplot(m5, aes(x=-log10(pval), color=kind))+geom_density()
ggplot(m5, aes(x=logLik, color=kind))+geom_density()


dfvel_exp_fil <- dfvel_exp[dfvel_exp$pval<0.05 & dfvel_exp$R2 > 0.7,]
dfvel_exp_fil <- dfvel_exp_fil[!is.na(dfvel_exp_fil$pval),]
# gp!
all_avg <- sapply(unique(dfvel_exp_fil$expg), function(x) summarize_expg(dfvel_exp_fil[dfvel_exp_fil$expg==x,]) )

avg <- t(all_avg)
colnames(avg) <- c('mean','sd','nmice')

m <- sapply(unique(dfvel_exp_fil$expg), function(x) get_mice(dfvel_exp_fil[dfvel_exp_fil$expg==x,]) )
avg <- as.data.frame(avg)
avg$mouse <- m
#write.table(avg, file="average_slope_exp.tsv", sep="\t", quote=F)
#write.table(dfvel, file="all_models_exp.tsv", sep="\t", quote=F)

dd <- as.data.frame(table(substr(avg$mouse, 0,10)))
re <- dd[dd$Freq>1,'Var1']
avg$m <- substr(avg$mouse, 0, 10)
stable <- avg[avg$m %in% as.character(re),]


ggplot(data=stable, aes(x=m, y=mean,))+geom_point()+coord_flip()

avgnn <- avg[!is.na(avg$sd),]
avgnn$lower <- avgnn$mean - avgnn$sd
avgnn$upper <- avgnn$mean + avgnn$sd
avgnn$mouse <- as.character(avgnn$mouse)
ggplot(avgnn, aes(x=reorder(mouse,mean), y=mean)) +
        geom_point(position=position_dodge(1), stat="identity") +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(1))+theme_bw()+ggtitle('avg fit slope log(vol) ~ day')+theme(axis.text.x=element_blank())


## correlate with cetuxi and chemio

cetuxi <- read.table('/home/data/Dropbox/work/strata/pdxopedia/localsharedata/treats/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w.tsv', sep="\t", quote="", header=FALSE)
head(cetuxi)
colnames(cetuxi) <- c("case", "perc")
avgnn$case <- substr(avgnn$mouse, 0, 7)
head(avgnn)
irino <- read.table('/home/data/Dropbox/work/strata/pdxopedia/localsharedata/treats/CHEMIO_WATERFALL_PLOT_Eugy_2020maggio_w3.txt', sep="\t", quote="", header=TRUE)
mirino <- merge(avgnn, irino, by="case")
mcetuxi <- merge(avgnn, cetuxi, by="case")
dim(mirino)
dim(mcetuxi)
dim(cetuxi)
dim(irino)
#ggplot(mirino,aes(x,y))+geom_point()+geom_smooth(method='lm')
head(mirin)
head(mirino)
ggplot(mirino,aes(mean,perc))+geom_point()+geom_smooth(method='lm')
summary(lm(data=mirino, formula=as.formula(perc~lower)))
summary(lm(data=mirino, formula=as.formula(perc~mean)))
cor.test(mirino$mean, mirino$perc)
ggplot(mcetuxi,aes(mean,perc))+geom_point()+geom_smooth(method='lm')
summary(mcetuxi$mean)
summary(mcetuxi$perc)
cor.test(mcetuxi$mean, mcetuxi$perc)
summary(lm(data=mcetuxi, formula=as.formula(perc~mean)))
dim(cetuxi)
dim(mcetuxi)
