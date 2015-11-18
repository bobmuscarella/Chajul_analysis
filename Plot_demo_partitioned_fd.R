##################################################
### PLOT DEMOGRAPHICALLY-PARTIIONED CWM AND FD ###
##################################################
library(lme4)
library(MuMIn)

# feve, fdis, fric, rao, cwm, ba
pdf('feve.pdf')

par(mfrow=c(3,3), mar=c(4,4,2,2), oma=c(4,4,0,0))
for(setpool in 1:4){
  pool <- setpool
  metric <- 'feve'

for(trt in 1:length(res)){
  p <- c('ba0','baG','baR','baM')[pool]
  col <- c(rgb(0,0,0,1), rgb(0,0.6,0,1), rgb(0,0,0.6,1), rgb(0.6,0,0,1))[pool]
  t <- names(res)[trt]
  m <- metric
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
    nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
    y <- fixef(mod)[1] + (nd$age * fixef(mod)[2]) + (nd$age2 * fixef(mod)[3])
    est <- quantile(fixef(sim(mod, n.sims=500))[,2], probs=c(0.025, 0.975))
    lty <- ifelse(sign(est[1]) == sign(est[2]), 1, 2)
    lcol <- ifelse(lty==1, 2, 1)
    lwd <- ifelse(lty==1, 3, 1) 
    points(nd$age, y, type='l', col=lcol, lwd=lwd, lty=lty)
  } else {  
    nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
    y <- fixef(mod)[1] + (nd$age * fixef(mod)[2])
    est <- quantile(fixef(sim(mod, n.sims=500))[,2], probs=c(0.025, 0.975))
    lty <- ifelse(sign(est[1]) == sign(est[2]), 1, 2)
    lcol <- ifelse(lty==1, 2, 1)
    lwd <- ifelse(lty==1, 3, 1) 
    points(nd$age, y, type='l', col=lcol, lwd=lwd, lty=lty)  
  }  
  r2 <- round(r.squaredGLMM(mod),2)
  mtext(paste('R2m =',r2[1],',','R2c =',r2[2]), cex=.8)
}
plot(0,0,cex=NA,axes=F,ylab='',xlab='')
plot(0,0,cex=NA,axes=F,ylab='',xlab='')
plot(0,0,cex=NA,axes=F,ylab='',xlab='')
}
dev.off()






##########################################
### PLOT DEMOGRAPHICALLY-PARTIIONED BA ###
##########################################
library(lme4)
library(MuMIn)

plot_t2 <- function(metric='cwm', pool='ba0', trait='SLA', col=NULL, ylab='trait', census, add=F){
  #ylim <- range(res[[trait]][[metric]], na.rm=T)
  vals <- res[[trait]][[metric]][,pool]
  sites <- rownames(res[[trait]][[metric]])
  x <- vals[match(as.character(census$SITE.CENSUS), sites)]  
  if(is.null(col)){
    col <- rainbow(length(unique(census$SITE)))
  } else {
    col <- rep(col, length(unique(census$SITE))) 
  }
  censusyr <- census$ages + census$CENSUS
  if(add==F){
    #plot(censusyr, x, xlab='Years', ylab=ylab, col=0, ylim=ylim)
    plot(censusyr, x, xlab='Years', ylab=ylab, col=0)    
  }
  for (i in 1:length(unique(census$SITE))){
    site <- unique(census$SITE)[i]
    points(censusyr[census$SITE==site], x[census$SITE==site], type='l', col=col[i], lwd=2)
  }
}



pdf('BA.pdf')

par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(4,4,0,0))
for(setpool in 1:4){
  pool <- setpool
  metric <- 'ba'
  
  #  for(trt in 1:length(res)){
  for(trt in 1){
    p <- c('ba0','baG','baR','baM')[pool]
    col <- c(rgb(0,0,0,1), rgb(0,0.6,0,1), rgb(0,0,0.6,1), rgb(0.6,0,0,1))[pool]
    t <- names(res)[trt]
    m <- metric
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
    plot_t2(ylab=paste(t, m, p), pool=p, trait=t, metric=m, census=census, col=col)
    if(AIC(qmod)+2 < AIC(lmod)){
      nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
      y <- fixef(mod)[1] + (nd$age * fixef(mod)[2]) + (nd$age2 * fixef(mod)[3])
      est <- quantile(fixef(sim(mod, n.sims=500))[,2], probs=c(0.025, 0.975))
      lty <- ifelse(sign(est[1]) == sign(est[2]), 1, 2)
      lcol <- ifelse(lty==1, 2, 1)
      lwd <- ifelse(lty==1, 3, 1) 
      points(nd$age, y, type='l', col=lcol, lwd=lwd, lty=lty)
    } else {  
      nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
      y <- fixef(mod)[1] + (nd$age * fixef(mod)[2])
      est <- quantile(fixef(sim(mod, n.sims=500))[,2], probs=c(0.025, 0.975))
      lty <- ifelse(sign(est[1]) == sign(est[2]), 1, 2)
      lcol <- ifelse(lty==1, 2, 1)
      lwd <- ifelse(lty==1, 3, 1) 
      points(nd$age, y, type='l', col=lcol, lwd=lwd, lty=lty)  
    }  
    r2 <- round(r.squaredGLMM(mod),2)
    mtext(paste('R2m =',r2[1],',','R2c =',r2[2]), cex=.8)
  }
}
dev.off()





