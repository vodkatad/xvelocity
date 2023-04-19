setwd('~/work/stash/velocity/new')
cetuxi <- read.table(gzfile('cetuxi6w_long_header.tsv.gz'), sep="\t", header=T)
pla <- read.table(gzfile('placebo_long_header.tsv.gz'), sep="\t", header=T)
head(cetuxi)
head(pla)
m <- merge(cetuxi, pla, by="exp_group")
head(m)
dim(m[m$start_date.x!=m$start_date.y,])
dim(m[m$start_date.x != m$start_date.y,])
dim(m[m$start_date.x == m$start_date.y,])
dim(m)
head(m[m$start_date.x == m$start_date.y,])
m[m$start_date.x == m$start_date.y,]
m$start_date.x <- as.Date(as.character(m$start_date.x), format="%Y-%m-%d")
m$start_date.y <- as.Date(as.character(m$start_date.y), format="%Y-%m-%d")
m$end_date.x <- as.Date(as.character(m$end_date.x), format="%Y-%m-%d")
m$end_date.y <- as.Date(as.character(m$end_date.y), format="%Y-%m-%d")

pla$start_date <- as.Date(as.character(pla$start_date), format="%Y-%m-%d")
cetuxi$start_date <- as.Date(as.character(cetuxi$start_date), format="%Y-%m-%d")
pla$end_date <- as.Date(as.character(pla$end_date), format="%Y-%m-%d")
cetuxi$end_date <- as.Date(as.character(cetuxi$end_date), format="%Y-%m-%d")
pla$measure_date <- as.Date(as.character(pla$measure_date), format="%Y-%m-%d")
cetuxi$measure_date <- as.Date(as.character(cetuxi$measure_date), format="%Y-%m-%d")

### placebo not dropping


fit_all_mice <-  function(expg_d) {
  expg_d <- expg_d[order(expg_d$measure_date),]
  vels <- sapply(unique(expg_d$longen), function(x) fit_one_mouse(expg_d[expg_d$longen==x,]) )
  res <- t(vels)
  colnames(res) <- c('slope','pval','R2', 'n_measures','logLik')
  res <- as.data.frame(res)
  res$exp_group <- unique(expg_d$exp_group)
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
  res = data.frame(volume=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'volume'] <- mouse_d[1,'volume']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'volume'] <- mouse_d[i, 'volume']
    res[i,'day'] <- mouse_d[i, 'measure_date'] - mouse_d[1,'measure_date']
  }
  model <- lm(data=res, formula="volume~day")
  sm <- summary(model)
  r2 <- sm$r.squared
  pval <- lmp(sm) 
  ll <- logLik(model)
  ggplot(res,aes(day,volume))+geom_point()+geom_smooth(method='lm')
  ggsave(paste0(unique(mouse_d$longen),".png"))
  # do we get only some R2/pvals?
  vel <- model$coefficients[2]
  attributes(vel) <- NULL
  return(c(vel, pval, r2, nrow(mouse_d),ll))
}


####

all_vels <- lapply(unique(pla$exp_group), function(x) fit_all_mice(pla[pla$exp_group==x,]) )
dfvel <- do.call(rbind, all_vels)
dfvel <- dfvel[!is.na(dfvel$pval),]
