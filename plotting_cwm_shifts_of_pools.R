## plotting function for CWM shifts over succession of different demographic pools

pdf(file='cwm_bam.pdf')
par(mfrow=c(2,2), mar=c(4,4,1,1))
for(trt in 1:length(res)){
  col <- rgb(0.6,0,0,1)
  lcol <- rgb(0.8,0,0,1)
  p <- 'baM'
  t <- names(res)[trt]
  m <- 'cwm'
  vals <- res[[t]][[m]][,p]
  sites <- rownames(res[[t]][[m]])
  x <- vals[match(as.character(census$SITE.CENSUS), sites)]  
  censusyr <- census$ages + census$CENSUS
  x <- cbind(x, censusyr, census$SITE)
  x <- x[!is.na(rowSums(x)),]
  rownames(x) <- NULL
  colnames(x) <- c('y','age','site')
  x <- as.data.frame(x)
  x$age2 <- x$age^2
  lmod <- lmer(y ~ age + (1 + age|site), data=x)
  qmod <- lmer(y ~ age + age2 + (1 + age|site), data=x)
  mod <- list(lmod, qmod)[[ifelse(AIC(qmod)+2 < AIC(lmod), 2, 1)]]
  plot_t(ylab=paste(t, m, p), pool=p, trait=t, metric=m, census=census, col=col)
  if(AIC(qmod)+2 < AIC(lmod)){
    for(i in 1:nrow(coef(qmod)$site)){
      nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
      y <- coef(qmod)$site[i,1] + (nd$age * coef(qmod)$site[i,2]) + (nd$age2 * coef(qmod)$site[i,3])
    }
    nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
    y <- fixef(mod)[1] + (nd$age * fixef(mod)[2]) + (nd$age2 * fixef(mod)[3])
    points(nd$age, y, type='l', col=lcol, lwd=4)
  } else {
    for(i in 1:nrow(coef(lmod)$site)){
      nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
      y <- coef(lmod)$site[i,1] + (nd$age * coef(lmod)$site[i,2])
    }
    nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
    y <- fixef(mod)[1] + (nd$age * fixef(mod)[2])
    points(nd$age, y, type='l', col=lcol, lwd=4)
  }
  r2 <- round(r.squaredGLMM(mod),2)
  mtext(bquote('R'['m']^2~'='~.(r2[1])~'; '~'R'['m']^2~'='~.(r2[2])), cex=.8, line=-2)
}
dev.off()












col <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))

pdf(file='cwm_ba0_color.pdf')
par(mfrow=c(2,2), mar=c(4,4,1,1))
for(trt in 1:length(res)){
  for(pool in 1){  
    #  for(pool in 1:4){
    p <- c('ba0','baG','baR','baM')[pool]
    #col <- c(rgb(0,0,0,1),rgb(0,0.6,0,1),rgb(0,0,0.6,1),rgb(0.6,0,0,1))[pool]
    #lcol <- c(rgb(0,0,0,.7),rgb(0,0.8,0,1),rgb(0,0,0.8,1),rgb(0.8,0,0,1))[pool]  
    t <- names(res)[trt]
    m <- 'cwm'
    vals <- res[[t]][[m]][,p]
    sites <- rownames(res[[t]][[m]])
    x <- vals[match(as.character(census$SITE.CENSUS), sites)]  
    censusyr <- census$ages + census$CENSUS
    x <- cbind(x, censusyr, census$SITE)
    x <- x[!is.na(rowSums(x)),]
    rownames(x) <- NULL
    colnames(x) <- c('y','age','site')
    x <- as.data.frame(x)
    x$age2 <- x$age^2
    lmod <- lmer(y ~ age + (1 + age|site), data=x)
    qmod <- lmer(y ~ age + age2 + (1 + age|site), data=x)
    mod <- list(lmod, qmod)[[ifelse(AIC(qmod)+2 < AIC(lmod), 2, 1)]]
    if(pool==1) {
      #plot_t(ylab=paste(t, m, p), pool=p, trait=t, metric=m, census=census)
      plot_t(ylab=paste(t, m, p), pool=p, trait=t, metric=m, census=census, col=col)
    }
    if(AIC(qmod)+2 < AIC(lmod)){
      for(i in 1:nrow(coef(qmod)$site)){
        nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
        y <- coef(qmod)$site[i,1] + (nd$age * coef(qmod)$site[i,2]) + (nd$age2 * coef(qmod)$site[i,3])
        #      points(nd$age, y, type='l', col=lcol, lwd=4)
      }
      nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
      y <- fixef(mod)[1] + (nd$age * fixef(mod)[2]) + (nd$age2 * fixef(mod)[3])
      points(nd$age, y, type='l', col=lcol, lwd=4)
    } else {
      for(i in 1:nrow(coef(lmod)$site)){
        nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
        y <- coef(lmod)$site[i,1] + (nd$age * coef(lmod)$site[i,2])
        #      points(nd$age, y, type='l', col=lcol, lwd=4)
      }
      nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
      y <- fixef(mod)[1] + (nd$age * fixef(mod)[2])
      points(nd$age, y, type='l', col=lcol, lwd=4)
    }
  }}
dev.off()


