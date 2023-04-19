# file esteso da formato wide a long da process_las.pl - il formato wide esce da query del LAS
# che pero` in produzione non funziona (da chiedere ad Ale, stanno cmq indagando).
d <- read.table("~/work/stash/velocity/CRC_phys_standard_cetux_long_expg.tsv", sep="\t", header=F)
colnames(d) <- c('gen','expg','date','vol')

# Filtriamo via i gruppi sperimentali con numero di topi/misure estremi
nums <- as.data.frame(table(d$expg))

png('histogram_nmeasuresmice_expg.png')
hist(nums$Freq)
graphics.off()

q <- quantile(nums$Freq, probs=c(0.01, 0.99))
nums$Var1 <- as.character(nums$Var1)
remove <- c(nums[nums$Freq <= q[1],'Var1'], nums[nums$Freq >= q[2],'Var1'])
dd <- d[! d$expg %in% remove,]

# filtriamo via un singolo esperimento con tre soli topi trattati con la fisiologica
topos <- sapply(unique(dd$expg), function(x) {y <- dd[dd$expg==x,]; length(unique(y$gen))})
#table(topos)
remove <- names(topos[topos == 3])
cat(remove)
ddd <- dd[dd$expg != remove,]
# assegno una data di tipo Date per calcolare comodamente i delta in giorni
ddd$realdate <- as.Date(as.character(ddd$date), format="%Y-%m-%d")

## Blocco fit lineari ###############################


# fit_all_mice <-  function(expg_d) {
#   expg_d <- expg_d[order(expg_d$realdate),]
#   vels <- sapply(unique(expg_d$gen), function(x) fit_one_mouse(expg_d[expg_d$gen==x,]) )
#   res <- t(vels)
#   colnames(res) <- c('slope','pval','R2', 'n_measures','logLik')
#   res <- as.data.frame(res)
#   res$expg <- unique(expg_d$expg)
#   res
# }
# 
# lmp <- function (sm) {
#   #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
#   #f <- summary(modelobject)$fstatistic
#   f <- sm$fstatistic
#   p <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(p) <- NULL
#   return(p)
# }
# 
# fit_one_mouse <- function(mouse_d) {
#   res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
#   res[1,'vol'] <- mouse_d[1,'vol']
#   for (i in seq(2, nrow(mouse_d))) {
#     res[i,'vol'] <- mouse_d[i, 'vol']
#     res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
#   }
#   model <- lm(data=res, formula="vol~day")
#   sm <- summary(model)
#   r2 <- sm$r.squared
#   pval <- lmp(sm) 
#   ll <- logLik(model)
#   #ggplot(res,aes(day,vol))+geom_point()+geom_smooth(method='lm')
#   #ggsave(paste0(unique(mouse_d$gen),".png"))
#   # do we get only some R2/pvals?
#   vel <- model$coefficients[2]
#   attributes(vel) <- NULL
#   return(c(vel, pval, r2, nrow(mouse_d),ll))
# }
# 
# 
#
# all_vels <- lapply(unique(ddd$expg), function(x) fit_all_mice(ddd[ddd$expg==x,]) )
# dfvel <- do.call(rbind, all_vels)
# dfvel <- dfvel[!is.na(dfvel$pval),]
# 
# dfvel$padj <- p.adjust(dfvel$pval)
# dfvel_sign <- dfvel[dfvel$padj < 0.05,]
# dfvel_fil <- dfvel[dfvel$pval<0.05 & dfvel$R2 > 0.7,]

summarize_expg <-  function(dfvel_g) {
  avg <- mean(dfvel_g$slope)
  sd <- sd(dfvel_g$slope)
  n <- nrow(dfvel_g)
  return(c(avg, sd, n))
}

get_mice <- function(dfvel_g) {
  unique(substr(rownames(dfvel_g),0,14))
}
# 
# 
# all_avg <- sapply(unique(dfvel_fil$expg), function(x) summarize_expg(dfvel_fil[dfvel_fil$expg==x,]) )
# 
# avg <- t(all_avg)
# colnames(avg) <- c('mean','sd','nmice')
# 
# m <- sapply(unique(dfvel_fil$expg), function(x) get_mice(dfvel_fil[dfvel_fil$expg==x,]) )
# avg <- as.data.frame(avg)
# avg$mouse <- m
# write.table(avg, file="average_slope.tsv", sep="\t", quote=F)
# write.table(dfvel, file="all_models.tsv", sep="\t", quote=F)
# 
# dd <- as.data.frame(table(substr(avg$mouse, 0,10)))
# re <- dd[dd$Freq>1,'Var1']
# avg$m <- substr(avg$mouse, 0, 10)
# stable <- avg[avg$m %in% as.character(re),]
# 
# 
# ggplot(data=stable, aes(x=m, y=mean,))+geom_point()+coord_flip()
# 
# avgnn <- avg[!is.na(avg$sd),]
# avgnn$lower <- avgnn$mean - avgnn$sd
# avgnn$upper <- avgnn$mean + avgnn$sd
# 
# 
# ggplot(avgnn, aes(x=reorder(mouse,mean), y=mean)) +
#   geom_point(position=position_dodge(1), stat="identity") +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(1))+theme_bw()+ggtitle('avg linear fit slope vol ~ day')+theme(axis.text.x=element_blank())

## blocco (vmax - vmin)/ (tmax-tmin) ###########################################

notfit_all_mice <-  function(expg_d) {
  expg_d <- expg_d[order(expg_d$realdate),]
  vels <- sapply(unique(expg_d$gen), function(x) notfit_one_mouse(expg_d[expg_d$gen==x,]) )
  res <- data.frame(slope=vels, expg=unique(expg_d$expg))
  return(res)
}


notfit_one_mouse <- function(mouse_d) {
  res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'vol'] <- mouse_d[1,'vol']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'vol'] <- mouse_d[i, 'vol']
    res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
  }
  vmax <- max(res$vol)
  vmin <- min(res$vol)
  imax <- match(vmax, res$vol) # in caso di parimerito tra volumi max prendiamo il primo in ordine di tempo
  imin <- match(vmin, res$vol) # idem per il min
  tmax <- res[imax,'day']
  tmin <- res[imin,'day']
  return((vmax-vmin)/(tmax-tmin))
}

# otteniamo i risultati di tutti i singoli topi dei vari gruppi sperimentali
all_vels_nofit <- lapply(unique(ddd$expg), function(x) notfit_all_mice(ddd[ddd$expg==x,]) )
dfvel_nofit <- do.call(rbind, all_vels_nofit)
# rimuovo cmq il dato del topo con sole due misure ma non faccio altri filtri
tt <- as.data.frame(table(ddd$gen))
remove <- as.character(tt[tt$Freq==2,'Var1'])

dfvel_fil <- dfvel_nofit[rownames(dfvel_nofit) != remove,]

# e ne calcoliamo la media/sd della velocita` risultante
all_avg <- sapply(unique(dfvel_fil$expg), function(x) summarize_expg(dfvel_fil[dfvel_fil$expg==x,]) )

avg <- t(all_avg)
colnames(avg) <- c('mean','sd','nmice')

m <- sapply(unique(dfvel_fil$expg), function(x) get_mice(dfvel_fil[dfvel_fil$expg==x,]) )
avg <- as.data.frame(avg)
avg$mouse <- as.character(m)
write.table(avg, file="average_vmaxvmin.tsv", sep="\t", quote=F)
#write.table(dfvel, file="all_models_vmaxvmin.tsv", sep="\t", quote=F)

# vediamo i singoli modelli come vanno ignorando XA/XB, se son simili o meno quando ci son piu` experimental group per loro
dd <- as.data.frame(table(substr(avg$mouse, 0,10)))
re <- dd[dd$Freq>1,'Var1']
avg$m <- substr(avg$mouse, 0, 10)
stable <- avg[avg$m %in% as.character(re),]

ggplot(data=stable, aes(x=m, y=mean))+geom_point()+coord_flip()
ggsave('differences_sameshortgenealogy_differentaarm.png')


avg$lower <- avg$mean - avg$sd
avg$upper <- avg$mean + avg$sd


ggplot(avg, aes(x=reorder(mouse,mean), y=mean)) +
  geom_point(position=position_dodge(1), stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(1))+theme_bw()+ggtitle('vmax-vmin estimate mean/sd')+theme(axis.text.x=element_blank())
ggsave('overall_picture.png')


# m <- merge(dfvel_nofit, dfvel, by="row.names")
# cor.test(m$slope, m$`unlist(all_vels_nofit)`)
# 
# colnames(m)[2] <- 'vmaxvmin'
# colnames(m)[3] <- 'slopefit'
# ggplot(m,aes(vmaxvmin,slopefit))+geom_point()+geom_smooth(method='lm')

### blocco fit  exp fit with log transf of y e confronto con i lineari ###################################################3
#https://stackoverflow.com/questions/31851936/exponential-curve-fitting-in-r

# efit_all_mice <-  function(expg_d) {
#   expg_d <- expg_d[order(expg_d$realdate),]
#   vels <- sapply(unique(expg_d$gen), function(x) efit_one_mouse(expg_d[expg_d$gen==x,]) )
#   res <- t(vels)
#   colnames(res) <- c('slope','pval','R2', 'n_measures', 'logLik')
#   res <- as.data.frame(res)
#   res$expg <- unique(expg_d$expg)
#   res
# }
# 
# 
# efit_one_mouse <- function(mouse_d) {
#   res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
#   res[1,'vol'] <- mouse_d[1,'vol']
#   for (i in seq(2, nrow(mouse_d))) {
#     res[i,'vol'] <- mouse_d[i, 'vol']
#     res[i,'day'] <- mouse_d[i, 'realdate'] - mouse_d[1,'realdate']
#   }
#   res$vol <- log(res$vol)
#   model <- lm(data=res, formula="vol~day")
#   sm <- summary(model)
#   r2 <- sm$r.squared
#   pval <- lmp(sm) 
#   ll <- logLik(model)
#   vel <- model$coefficients[2]
#   attributes(vel) <- NULL
#   return(c(vel, pval, r2, nrow(mouse_d),ll))
# }
# 
# all_expvels <- lapply(unique(ddd$expg), function(x) efit_all_mice(ddd[ddd$expg==x,]) )
# dfvel_exp <- do.call(rbind, all_expvels)
# 
# all_expvels <- lapply(unique(ddd$expg), function(x) efit_all_mice(ddd[ddd$expg==x,]) )
# dfvel_exp <- do.call(rbind, all_expvels)
# cor.test(m2$slope.y, m2$slope.y)
# dfvel_exp$kind <- 'exp'
# dfvel$kind <- 'linear'
# m2 <- merge(dfvel_exp, dfvel, by="row.names")
# cor.test(m2$slope.y, m2$slope.y)
# ggplot(m2,aes(slope.x,slope.y))+geom_point()+geom_smooth(method='lm')
# m3 <- m2[,c('slope.x','pval.x','R2.x','kind.x', 'logLik.x')]
# m4 <- m2[,c('slope.y','pval.y','R2.y','kind.y', 'logLik.y')]
# colnames(m3) <- c('slope','pval','R2','kind', 'logLik')
# colnames(m4) <- c('slope','pval','R2','kind', 'logLik')
# m5 <- rbind(m3,m4)
# ggplot(m5, aes(x=R2, color=kind))+geom_density()
# ggplot(m5, aes(x=-log10(pval), color=kind))+geom_density()
# ggplot(m5, aes(x=logLik, color=kind))+geom_density()
# 
# 
# dfvel_exp_fil <- dfvel_exp[dfvel_exp$pval<0.05 & dfvel_exp$R2 > 0.7,]
# dfvel_exp_fil <- dfvel_exp_fil[!is.na(dfvel_exp_fil$pval),]
# # gp!
# all_avg <- sapply(unique(dfvel_exp_fil$expg), function(x) summarize_expg(dfvel_exp_fil[dfvel_exp_fil$expg==x,]) )
# 
# avg <- t(all_avg)
# colnames(avg) <- c('mean','sd','nmice')
# 
# m <- sapply(unique(dfvel_exp_fil$expg), function(x) get_mice(dfvel_exp_fil[dfvel_exp_fil$expg==x,]) )
# avg <- as.data.frame(avg)
# avg$mouse <- m
# #write.table(avg, file="average_slope_exp.tsv", sep="\t", quote=F)
# #write.table(dfvel, file="all_models_exp.tsv", sep="\t", quote=F)
# 
# dd <- as.data.frame(table(substr(avg$mouse, 0,10)))
# re <- dd[dd$Freq>1,'Var1']
# avg$m <- substr(avg$mouse, 0, 10)
# stable <- avg[avg$m %in% as.character(re),]
# 
# 
# ggplot(data=stable, aes(x=m, y=mean,))+geom_point()+coord_flip()
# 
# avgnn <- avg[!is.na(avg$sd),]
# avgnn$lower <- avgnn$mean - avgnn$sd
# avgnn$upper <- avgnn$mean + avgnn$sd
# avgnn$mouse <- as.character(avgnn$mouse)
# ggplot(avgnn, aes(x=reorder(mouse,mean), y=mean)) +
#   geom_point(position=position_dodge(1), stat="identity") +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(1))+theme_bw()+ggtitle('avg fit slope log(vol) ~ day')+theme(axis.text.x=element_blank())
# 
