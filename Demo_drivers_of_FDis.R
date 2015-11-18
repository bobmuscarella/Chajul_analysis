###########################################################
#### CALCULATE DEMOGRAPHICALLY-PARTITIONED FD, CWM, BA ####
###########################################################

library(RColorBrewer)
library(FD)

###################
#### LOAD DATA ####
###################
setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/Neoselvas/Chajul/DATA")
load("Chajul_data_processed_wtraits_4.27.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
tdata <- read.csv('Mean_traits_4.27.15.csv')
head(data)
head(census)

### CHECK FOR OUTLIERS...
outliers <- vector()
for(i in -30:30){
  outliers[i] <- length(data$uID[!is.na(data$growth) & data$growth > i])
}
plot(outliers, type='l', ylim=c(0,30))
# data <- data[!is.na(data$growth) & data$growth < 8,]

### FUNCTION TO CONVERT DBH TO BA
dbh2ba <- function(dbh){
  return(((pi / 40000) * (dbh^2)))
}


###################
#### PREP DATA ####
###################
### IDENTIFY LAST CENSUS AT EACH SITE (NO GROWTH/MORTALITY DATA)
lastcensus <- tapply(census$year, census$SITE2, max)
yr <- ifelse(substring(lastcensus,3,3)==0,substring(lastcensus, 4, 4), substring(lastcensus, 3, 4))
lastcensus <- paste(names(lastcensus), yr, sep='.')

### IDENTIFY FIRST CENSUS AT EACH SITE (NO RECRUITMENT DATA)
firstcensus <- tapply(census$year, census$SITE2, min)
yr <- ifelse(substring(firstcensus,3,3)==0,substring(firstcensus, 4, 4), substring(firstcensus, 3, 4))
firstcensus <- paste(names(firstcensus), yr, sep='.')

### EXCLUDE SITES WITH FEWER THAN 2 CENSUSES
number_censuses <- tapply(data$SITE.CENSUS, data$SITE, function(x) length(unique(x)))
focal_sites <- names(number_censuses)[number_censuses>1]
data <- data[data$SITE %in% focal_sites,]

### REMOVE INDIVIDUALS WITH NO DBH RECORDED (???)
dbhlist <- split(data$DBH, data$uID)
tmp <- unlist(lapply(dbhlist, function(x) sum(is.na(x))!=length(x)))
data <- data[data$uID %in% names(tmp)[tmp==T],]
data <- droplevels(data)
dbhlist <- split(data$DBH, data$uID)

### GET THE RECRUITMENT YEAR FOR EACH INDIVIDUAL
yearlist <- split(data$year, data$uID)
year <- unlist(lapply(dbhlist, function(x) min(which(!is.na(x)))-1))
year.recruit <- rep(NA, length(yearlist))
for(i in 1:length(yearlist)){
	if(year[i] == 0){ year.recruit[i] <- NA	} else { year.recruit[i] <- yearlist[[i]][year[i]]	}
}
names(year.recruit) <- unique(data$uID)
data$year.recruit <- year.recruit[match(data$uID, names(year.recruit))]
data$year.recruit <- data$year.recruit+1
data$recruit <- data$year == data$year.recruit
data$recruit[is.na(data$recruit)] <- FALSE

### REMOVE DATA FROM FIRST AND LAST CENSUSES (OR LAST CENSUS ONLY)
data <- data[! data$SITE.CENSUS %in% c(firstcensus, lastcensus), ]
data <- droplevels(data)


###########################
#### MAKE CALCULATIONS ####
###########################
res <- list()

tdata$log.SLA <- log(tdata[,'SLA'])
tdata$log.sFtP <- log(tdata[,'sFtP'])
tdata$log.SV <- log(tdata[,'SV'])
tdata$log.P <- log(tdata[,'P'])

focaltraits <- which(names(tdata) %in% c('log.SLA','log.sFtP','WD','log.SV','log.P','QY'))

for (t in 1:length(focaltraits)){
  trait <- names(tdata)[focaltraits[t]]
  print(paste('*** Now workingo on trait ',trait,'***'))
	n <- length(unique(data$SITE.CENSUS))
  results.FDis <- as.data.frame(matrix(ncol=10, nrow=n))
  ba <- cwm <- results.FRic <- results.RaoQ <- results.FEve <- results.FDis

  for (i in 1:length(unique(data$SITE.CENSUS))){	
  	site <- unique(data$SITE.CENSUS)[i]
  	d <- data[data$SITE.CENSUS %in% site, ]
### Get initial basal area
  	ba0 <- tapply(d$ba[d$recruit == F], d$SPECIES[d$recruit == F], sum, na.rm=T)
  	ba0[is.na(ba0)] <- 0
### Get increase of basal area due to growth
  	d$baG <- dbh2ba(d$DBH + d$growth) - dbh2ba(d$DBH)
  	baG <- tapply(d$baG, d$SPECIES, sum, na.rm=T)
  	baG[is.na(baG)] <- 0
### Get increase of basal area due to recruitment
  	baR <- tapply(d$ba[d$recruit == T], d$SPECIES[d$recruit == T], sum, na.rm=T)
  	baR[is.na(baR)] <- 0
### Get decrease of basal area due to mortality
  	baM <- tapply(d$ba[d$survive %in% 0], d$SPECIES[d$survive %in% 0], sum, na.rm=T)
  	baM[is.na(baM)] <- 0
### Get other BA combinations
    ba0G <- ba0 + baG
  	ba0R <- ba0 + baR
    ba0GR <- ba0 + baR + baG
    ba0M <- ba0 - baM	
    baGR <- baG + baR
    ba1 <- ba0 + baR + baG - baM
    trt <- tdata[,trait]
    names(trt) <- as.character(tdata$species)
    trt <- trt[match(names(ba0), names(trt))]
    names(trt) <- names(ba0)
    x <- as.data.frame(cbind(trt, ba0, ba0G, ba0R, ba0GR, ba0M, baG, baR, baGR, baM, ba1))
### EXCLUDE SPECIES NOT OBSERVED IN SITE AND THOSE WITHOUT TRAIT DATA
  	x <- x[rowSums(x[,-1]) != 0,]
  	x <- x[!is.na(x$trt),]
### CALCULATE THE METRICS
    if(nrow(x)>0){
  		ba[i,] <- colSums(x[,-1])
      cwm[i,] <- colSums(sweep(x[,-1], MARGIN=1, x[,1], '*')) / colSums(x[,-1])
  		tmp <- x[,1]
      names(tmp) <- rownames(x)
      d <- dist(tmp)
  		omits <- which(colSums(x) == 0)
  	  x[,omits] <- 1
### NEED TO SET NEGATIVE GROWTH TO 0 TO AVOID ERRORS...
      x[x<0] <- 0
  	  tmpres <- dbFD(d, t(x[,-1]), messages=F)
      tmpres$RaoQ[omits] <- tmpres$FDis[omits] <- tmpres$FEve[omits] <- tmpres$FRic[omits] <- NA
  		results.FDis[i,] <- tmpres$FDis
  		results.FEve[i,] <- tmpres$FEve
  		results.FRic[i,] <- tmpres$FRic
      results.RaoQ[i,] <- tmpres$RaoQ
  	}
    print(paste(i, 'out of', length(unique(data$SITE.CENSUS))))
	}
  rownames(cwm) <- rownames(ba) <- rownames(results.FRic) <- rownames(results.RaoQ) <- rownames(results.FEve) <- rownames(results.FDis) <- unique(data$SITE.CENSUS)
  res[[t]] <- list(ba, cwm, results.FDis, results.FEve, results.FRic, results.RaoQ)
  names(res[[t]]) <- c('ba','cwm','fdis','feve','fric','rao')
  names(res)[[t]] <- trait
  for (j in 1:length(res[[t]])){
    names(res[[t]][[j]]) <- c('ba0','ba0G','ba0R','ba0GR','ba0M','baG','baR','baGR','baM','ba1')  
  }
}

# save(res, file="demo_partioning_of_FD.RDA")
load("demo_partioning_of_FD.RDA")
str(res)




########################################################
#### SOME PLOTTING OF DEMOG-PARTITIONED FD, CWM, BA ####
########################################################

plot_t <- function(metric='cwm', pool='ba0', trait='SLA', col=NULL, ylab='trait', census, add=F){
  ylim <- range(res[[trait]][[metric]], na.rm=T)
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
    plot(censusyr, x, xlab='Years', ylab=ylab, col=0, ylim=ylim)
#    plot(censusyr, x, xlab='Years', ylab=ylab, col=0)    
  }
  for (i in 1:length(unique(census$SITE))){
    site <- unique(census$SITE)[i]
    points(censusyr[census$SITE==site], x[census$SITE==site], type='l', col=col[i], lwd=2)
  }
}







par(mfrow=c(3,3), mar=c(4,4,2,2), oma=c(4,4,0,0))

library(lme4)
library(MuMIn)

for(trt in 1:length(res)){
p <- 'ba0'
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

plot_t(ylab=paste(t, m, p), pool=p, trait=t, metric=m, census=census, col=rgb(0,0,0,1))

if(AIC(qmod)+2 < AIC(lmod)){
  for(i in 1:nrow(coef(qmod)$site)){
    nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
    y <- coef(qmod)$site[i,1] + (nd$age * coef(qmod)$site[i,2]) + (nd$age2 * coef(qmod)$site[i,3])
    points(nd$age, y, type='l', col=rgb(1,0,0,0.5))
  }
  nd <- data.frame(age=(-1:35), age2=((-1:35)^2), site=rep(1,length((-1:35))))
  y <- fixef(mod)[1] + (nd$age * fixef(mod)[2]) + (nd$age2 * fixef(mod)[3])
  points(nd$age, y, type='l', col=4, lwd=3)
} else {
  for(i in 1:nrow(coef(lmod)$site)){
    nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
    y <- coef(lmod)$site[i,1] + (nd$age * coef(lmod)$site[i,2])
    points(nd$age, y, type='l', col=rgb(0,0,1,0.5))
  }
  nd <- data.frame(age=(-1:35), site=rep(1,length((-1:35))))
  y <- fixef(mod)[1] + (nd$age * fixef(mod)[2])
  points(nd$age, y, type='l', col=2, lwd=3)
}

r2 <- round(r.squaredGLMM(mod),2)
mtext(paste('R2m =',r2[1],',','R2c =',r2[2]), cex=.8)
}











######################
###   Plotting CWM   #####
######################
col <- rainbow(length(unique(census$SITE)))

par(mfrow=c(2,2))

	trait <- names(tdata)[-1][t]												# Not using log traits
	# trait <- names(data)[c(23,28:32, 35, 36:44)][t]				# Using log traits

	census$cwm0 <- results[,1][match(as.character(census$SITE.CENSUS), rownames(results))]
	census$cwmG <- results[,2][match(as.character(census$SITE.CENSUS), rownames(results))]
	census$cwmR <- results[,3][match(as.character(census$SITE.CENSUS), rownames(results))]
	census$cwmM <- results[,4][match(as.character(census$SITE.CENSUS), rownames(results))]

# PLOT CWM0
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwm0, col=0, xlab='Years since abandonment', ylab='Trait value', 
			main=paste('Mean',trait,'at Start',sep=' '), ylim=c(0.2,0.8))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$cwm0[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT CWMG
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwmG, col=0, xlab='Years since abandonment', ylab='Trait value', 
			main=paste('Mean',trait,'of Growth',sep=' '), ylim=c(0.2,0.8))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$cwmG[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT CWMR
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwmR, col=0, xlab='Years since abandonment', ylab='Trait value', main=paste('Mean',trait,'of Recruitment',sep=' '), ylim=c(0.2,0.8))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$cwmR[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT CWMM
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwmM, col=0, xlab='Years since abandonment', ylab='Trait value',
			main=paste('Mean',trait,'of Mortality',sep=' '), ylim=c(0.2,0.8))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$cwmM[census$SITE==site], type='l', col=col[i], lwd=2)
	}





#####################
###   Plotting FDis   #####
#####################
col <- rainbow(length(unique(census$SITE)))

par(mfrow=c(2,2))

	trait <- names(tdata)[-1][t]												# Not using log traits
	# trait <- names(data)[c(23,28:32, 35, 36:44)][t]				# Using log traits

	res <- list(results.FRic, results.FDis, results.FEve)[[1]]

	census$FD0 <- res[,1][match(as.character(census$SITE.CENSUS), rownames(res))]
	census$FDG <- res[,2][match(as.character(census$SITE.CENSUS), rownames(res))]
	census$FDR <- res[,3][match(as.character(census$SITE.CENSUS), rownames(res))]
	census$FDM <- res[,4][match(as.character(census$SITE.CENSUS), rownames(res))]

# PLOT FDis0
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$FD0, col=0, xlab='Years since abandonment', ylab='FD', 
			main=paste('FD of Start',sep=' '))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$FD0[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT FDisG
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$FDG, col=0, xlab='Years since abandonment', ylab='FD', 
			main=paste('FD of Growth',sep=' '))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$FDG[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT FDisR
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$FDR, col=0, xlab='Years since abandonment', ylab='FD', 
			main=paste('FD of Recruitment',sep=' '))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$FDR[census$SITE==site], type='l', col=col[i], lwd=2)
	}

# PLOT FDisM
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$FDM, col=0, xlab='Years since abandonment', ylab='FD',
			main=paste('FD of Mortality',sep=' '))
	for (i in 1:length(unique(census$SITE))){
		site <- unique(census$SITE)[i]
		points(censusyr[census$SITE==site], census$FDM[census$SITE==site], type='l', col=col[i], lwd=2)
	}




