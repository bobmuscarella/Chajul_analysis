library(RColorBrewer)
library(FD)

###################
#### LOAD DATA ####
###################
setwd("/Users/Bob/Projects/Neoselvas/Chajul/DATA")
load("Chajul_data_processed_wtraits_4.27.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
head(data)
head(census)

#########################
#### LOAD FUNCTIONS ####
#########################
DemographicCWM <- function(ba0, trait, ba.g, ba.r, ba.m){
	ba1 <- ba0 + ba.g + ba.r - ba.m
	x <- as.data.frame(cbind(trait, ba0, ba.g, ba.r, ba.m, ba1))
	x <- x[!is.nan(x$trait),]
	cwm0 <- sum(x$ba0 * x$trait) / sum(x$ba0)
	cwm.g <- sum(x$ba.g * x$trait) / sum(x$ba.g)
	cwm.g <- ifelse(is.nan(cwm.g), 0, cwm.g)
	cwm.r <- sum(x$ba.r * x$trait) / sum(x$ba.r)
	cwm.r <- ifelse(is.nan(cwm.r), 0, cwm.r)
	cwm.m <- sum(x$ba.m * x$trait) / sum(x$ba.m)
	cwm.m <- ifelse(is.nan(cwm.m), 0, cwm.m)
	cwm1 <- sum(x$ba1 * x$trait) / sum(x$ba1)
	cwm.flux <- (cwm0*sum(x$ba0)) + (cwm.g * sum(x$ba.g)) + (cwm.r * sum(x$ba.r)) + (cwm.m * sum(x$ba.m))
	t0.prop <- (cwm0 * sum(x$ba0)) / cwm.flux
	g.prop <- (cwm.g * sum(x$ba.g)) / cwm.flux
	r.prop <- (cwm.r * sum(x$ba.r)) / cwm.flux
	m.prop <- (cwm.m * sum(x$ba.m)) / cwm.flux
	out <- round(c(cwm0, cwm1, cwm.g, cwm.r, cwm.m, t0.prop, g.prop, r.prop, m.prop), 3)
	names(out) <- c('cwm0', 'cwm1', 'cwm.g', 'cwm.r', 'cwm.m', 't0', 'Growth', 'Recruitment', 'Mortality')
	return(out)
}


PlotDemographicCWM <- function(ba0, trait, ba.g, ba.r, ba.m){
	cwm0 <- sum(ba0/sum(ba0) * trait)
	cwm.g <- sum(ba.g * trait) / sum(ba.g)
	cwm.g <- ifelse(is.nan(cwm.g), 0, cwm.g)
	cwm.r <- sum(ba.r * trait) / sum(ba.r)
	cwm.r <- ifelse(is.nan(cwm.r), 0, cwm.r)
	cwm.m <- sum(ba.m * trait) / sum(ba.m)
	cwm.m <- ifelse(is.nan(cwm.m), 0, cwm.m)
	ba1 <- ba0 + ba.g + ba.r - ba.m
	cwm1 <- sum(ba1/sum(ba1) * trait)
	ylim <- range(cwm0, cwm1, cwm.g, cwm.r, cwm.m)
	xlim <- c(min(ba0, ba.g, ba.r, ba.m), max((ba0+sum(ba0)), (ba1+sum(ba1)), (ba.g+sum(ba.g)), (ba.r+sum(ba.r)), ba.m+sum(ba.m)))
	plot(c(sum(ba0), sum(ba1)), c(cwm0,cwm1), pch=21, ylim=ylim, xlim=xlim, cex=3, bg=c(0,'grey'), xlab='Basal Area', ylab='CWM trait value')
	arrows(sum(ba0), cwm0, sum(ba1), cwm1, len=.1, lwd=3)
	abline(v=sum(ba0), h=cwm0, lty=3)
	points(c(sum(ba.m), sum(ba.g), sum(ba.r)), c(cwm.m, cwm.g, cwm.r), pch=21, bg=2:4, cex=3)
	cwm.flux <- (cwm0*sum(ba0)) + (cwm.g * sum(ba.g)) + (cwm.r * sum(ba.r)) + (cwm.m * sum(ba.m))
	t0.prop <- (cwm0 * sum(ba0)) / cwm.flux
	g.prop <- (cwm.g * sum(ba.g)) / cwm.flux
	r.prop <- (cwm.r * sum(ba.r)) / cwm.flux
	m.prop <- (cwm.m * sum(ba.m)) / cwm.flux
	mtext(paste("t0", round(t0.prop, 3), "| Growth: ", round(g.prop, 3), "| Recruitment: ", round(r.prop, 3), "| Mortality: ", round(m.prop, 3)))
	abline(h=c(cwm.m, cwm.g, cwm.r), v=c(sum(ba.m),sum(ba.g), sum(ba.r)), col=2:4, lwd=.5)
}

dbh2ba <- function(dbh){
	return((pi / 40000) * (dbh^2))
}


### WORKING ON HOW TO PARTITION FUNCTIONAL DIVERSITY....

DemographicFD <- function(ba0, trait, ba.g, ba.r, ba.m){

	ba1 <- ba0 + ba.g + ba.r - ba.m
	ba0g <- ba0 + ba.g
	ba0r <- ba0 + ba.r
	ba0m <- ba0 - ba.m


# Temp fixes to avoid negative biomass...
	ba0 <- ifelse(ba0<0, 0, ba0)
	ba.g <- ifelse(ba.g<0, 0, ba.g)
	ba.r <- ifelse(ba.r<0, 0, ba.r)
	ba.m <- ifelse(ba.m<0, 0, ba.m)
	ba0g <- ifelse(ba0g<0, 0, ba0g)
	ba0r <- ifelse(ba0r<0, 0, ba0r)
	ba0m <- ifelse(ba0m<0, 0, ba0m)



#	a <- t(cbind(ba0, ba.g, ba.r, ba.m))
	a <- t(cbind(ba0, ba0g, ba0r, ba0m, ba.g, ba.r, ba.m, ba1))

	colnames(a) <- 1:ncol(a)

	if(class(trait) == 'numeric') { names(trait) <- 1:length(trait) }
	if(class(trait) == 'matrix') { rownames(trait) <- 1:nrow(trait) }

	d <- dist(trait)
	fd <- dbFD(d, a, stand.FRic=T, calc.CWM=F)

# fd

# fd$FRic
# fd$FRic['ba.g'] * sum(ba.g)/sum(ba0)

	# t0.prop <- (fd[[1]][[1]] * sum(ba0)) / sum(ba1)
	# g.prop <- (fd[[1]][[2]] * sum(ba.g)) / sum(ba1)
	# r.prop <- (fd[[1]][[3]] * sum(ba.r)) / sum(ba1)
	# m.prop <- (fd[[1]][[4]] * sum(ba.m)) / sum(ba1)

	FRic <- fd[[3]]
	FEve <- fd[[5]]
	FDiv <- fd[[6]]
	FDis <- fd[[7]]

	x <- as.data.frame(rbind(FRic, FEve, FDiv, FDis))
	x$Metric <- rownames(x)
	rownames(x) <- NULL

#	out <- c(round(cwm0,3), round(cwm1,3), round(t0.prop, 3), round(g.prop, 3), round(r.prop, 3), round(m.prop, 3))
#	names(out) <- c('cwm0', 'cwm1', 't0', 'Growth', 'Recruitment', 'Mortality')
#	return(out)
	return(x)
}








# COMPUTE DEMOGRAPHIC DRIVERS OF CWM SHIFTS

out <- matrix(nrow=length(unique(data$SITE.CENSUS)), ncol=5)
FDout <- data.frame()

tdata <- tdata[order(tdata$species),]
trait <- 'WD'

for (i in 1:length(unique(data$SITE.CENSUS))){

	site <- unique(data$SITE.CENSUS)[i]
	d <- data[data$SITE.CENSUS == site, ]

	# Get initial basal area
	ba0 <- tapply(d$ba[d$recruit == F], d$SPECIES[d$recruit == F], sum, na.rm=T)
	ba0[is.na(ba0)] <- 0
	
	# Get increase of basal area due to growth
	d$ba.g <- dbh2ba(d$DBH + d$growth) - dbh2ba(d$DBH)
	ba.g <- tapply(d$ba.g, d$SPECIES, sum, na.rm=T)
	ba.g[is.na(ba.g)] <- 0

	# Get increase of basal area due to recruitment
	ba.r <- tapply(d$ba[d$recruit == T], d$SPECIES[d$recruit == T], sum, na.rm=T)
	ba.r[is.na(ba.r)] <- 0
	
	# Get decrease of basal area due to mortality
	ba.m <- tapply(d$ba[d$survive %in% 0], d$SPECIES[d$survive %in% 0], sum, na.rm=T)
	ba.m[is.na(ba.m)] <- 0

	trt <- tdata[,trait]
	names(trt) <- as.character(tdata$species)
	trt <- trt[match(names(ba0), names(trt))]
	names(trt) <- names(ba0)

	x <- as.data.frame(cbind(trt, ba0, ba.g, ba.r, ba.m))
	x <- x[!is.na(x$trt),]

	x <- x[rowSums(x[, colnames(x) != 'trt']) > 0,]
	x$site <- NULL

	if(sum(colSums(x) > 0) == 5){	# error trap for zero species communities...

	x <- DemographicFD(x$ba0, x$trait, x$ba.g, x$ba.r, x$ba.m)
	x$site <- site
	FDout <- rbind(FDout, x)
	out[i,] <- c(site, round(c(sum(x$ba0), sum(x$ba.g), sum(x$ba.r), sum(x$ba.m)), 3))
	print(i)

	}

}


out <- as.data.frame(out)
names(out) <- c('SITE.CENSUS', names(DemographicCWM(x$ba0, x$trait, x$ba.g, x$ba.r, x$ba.m)), 'BA.0','BA.g',"BA.r","BA.m")







#############################
#### SIMULATION SANDBOX ####
#############################
nsp <- 30
ba0=rep(1, nsp)
# trait=abs(rnorm(nsp))  # for single trait analysis
trait=cbind(abs(rnorm(nsp)), abs(rnorm(nsp)), abs(rnorm(nsp)), abs(rnorm(nsp))) # for multivariate analyses
ba.g=abs(rnorm(nsp))
ba.r=abs(rnorm(nsp))
ba.m=runif(nsp, 0, 9)

ba.g[sample(1:nsp, runif(nsp, 1, nsp/2))] <- 0
ba.r[sample(1:nsp, runif(nsp, 1, nsp/2))] <- 0
ba.m[sample(1:nsp, runif(nsp, 1, nsp))] <- 0





PlotDemographicCWM(ba0, trait, ba.g, ba.r, ba.m)


DemographicCWM(ba0, trait, ba.g, ba.r, ba.m)
tba <- sum(ba0, ba.g, ba.r, ba.m)
plot(DemographicCWM(ba0, trait, ba.g, ba.r, ba.m)[3:6], c(sum(ba0), sum(ba.g), sum(ba.r), sum(ba.m))/tba, pch=21, bg=c(1,3,4,2), cex=3)
abline(0,1)


########################################
#### TESTING SOME DIVERSITY PARTITIONING IDEAS ####
########################################

nsp <- 4
ba0=rep(1, nsp)
trait=t(t(c(.1,.2,.3,.4)))  # for single trait analysis
# trait=abs(rnorm(nsp))  # for single trait analysis
# trait=cbind(abs(rnorm(nsp)), abs(rnorm(nsp)), abs(rnorm(nsp)), abs(rnorm(nsp))) # for multivariate analyses
ba.g=abs(rnorm(nsp))
ba.r=abs(rnorm(nsp))
ba.m=runif(nsp, 0, .1)

ba0 <- t(as.matrix(ba0))
colnames(ba0) <- letters[1:4]
rownames(trait) <- letters[1:4]
res.ba0 <- dbFD(trait, ba0)
res.ba0

bag <- t(as.matrix(ba.g))
colnames(bag) <- letters[1:4]
res.bag <- dbFD(trait, bag)
res.bag

bam <- t(as.matrix(ba.m))
colnames(bam) <- letters[1:4]
res.bam <- dbFD(trait, bam)
res.bam

ba1 <- ba0+bag-bam
res.ba1 <- dbFD(trait, ba1)
res.ba1


res <- cbind(unlist(res.ba0), unlist(res.bag), unlist(res.bam), unlist(res.ba1))
tmp <- t(rbind(ba0, bag, bam, ba1))
tmp <- rbind(tmp, colSums(tmp))
rownames(tmp)[5] <- 'totalba'
res <- rbind(tmp, res)
colnames(res) <- c('ba0','bag','bam','ba1')
res

rba0 <- res['totalba',1]/res['totalba',4]
rbag <- res['totalba',2]/res['totalba',4]
rbam <- res['totalba',3]/res['totalba',4]

(res['FEve',1] * rba0) + (res['FEve',2] * rbag) - (res['FEve',3] * rbam)
res['FEve',4]


(res['FDiv',1] * rba0) + (res['FDiv',2] * rbag) - (res['FDiv',3] * rbam)
res['FDiv',4]

(res['FDis',1] * rba0) + (res['FDis',2] * rbag) - (res['FDis',3] * rbam)
res['FDis',4]

all.equal((res['CWM.x',1] * rba0) + (res['CWM.x',2] * rbag) - (res['CWM.x',3] * rbam), res['CWM.x',4])


########################################
#### GENERATE PLOTTING ARGUMENTS ####
########################################
col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
col <- c(col, col)

###########################################
#### CWM TRAITS WITH AGE (ABUNDANCE) ####
###########################################
par(mfrow=c(4,4), mar=c(2,2,2,2))
for (t in 1:16){
	trait <- names(tdata)[-1][t]												# Not using log traits
	# trait <- names(data)[c(23,28:32, 35, 36:44)][t]				# Using log traits
	d <- data[data$status=='alive',]
	cwm <- tapply(d[,trait], d$SITE.CENSUS, mean, na.rm=T)
	census$cwm <- cwm[match(as.character(census$SITE.CENSUS), names(cwm))]
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwm, col=0, xlab='Years since abandonment', ylab='Stems per plot', main=trait)
		for (i in 1:length(unique(census$SITE))){
			site <- unique(census$SITE)[i]
			points(censusyr[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=2)
		}
}

##########################################
#### CWM TRAITS WITH AGE (BASAL AREA) ####
##########################################
quartz()
par(mfrow=c(4,4), mar=c(2,2,2,2))
for (t in 1:16){
	trait <- names(tdata)[-1][t]												# Not using log traits
	# trait <- names(data)[c(23,28:32, 35, 36:44)][t]				# Using log traits
	tmpcwm <- data[,trait] * data$relBA
	cwm <- tapply(tmpcwm, data$SITE.CENSUS, sum, na.rm=T)
	census$cwm <- cwm[match(as.character(census$SITE.CENSUS), names(cwm))]
	censusyr <- census$ages+census$CENSUS
	plot(censusyr, census$cwm, col=0, xlab='Years since abandonment', ylab='Trait value', main=trait)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=2)
	}
}

######################################################
#### TOTAL BASAL AREA INCREASES WITH SUCCESSION ####
######################################################

censusyr <- census$ages+census$CENSUS
plot(censusyr, census$totalBA, col=0, xlab='Years since abandonment', ylab='Total BA')
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], census$totalBA[census$SITE==site], type='l', col=col[i], lwd=2)
}

######################################################
#### FORMAT DATA FOR DEMOGRAPHIC CWM ANALYSIS ####
######################################################

# IDENTIFY LAST CENSUS AT EACH SITE (NO GROWTH/MORTALITY DATA)
lastcensus <- tapply(census$year, census$SITE2, max)
yr <- ifelse(substring(lastcensus,3,3)==0,substring(lastcensus, 4, 4), substring(lastcensus, 3, 4))
lastcensus <- paste(names(lastcensus), yr, sep='.')

# IDENTIFY FIRST CENSUS AT EACH SITE (NO RECRUITMENT DATA)
firstcensus <- tapply(census$year, census$SITE2, min)
yr <- ifelse(substring(firstcensus,3,3)==0,substring(firstcensus, 4, 4), substring(firstcensus, 3, 4))
firstcensus <- paste(names(firstcensus), yr, sep='.')

# EXCLUDE SITES WITH FEWER THAN 3 CENSUSES
number_censuses <- tapply(data$SITE.CENSUS, data$SITE, function(x) length(unique(x)))
focal_sites <- names(number_censuses)[number_censuses>2]
tmpdata <- data[data$SITE %in% focal_sites,]

# REMOVE INDIVIDUALS WITH NO DBH RECORDED (???)
dbhlist <- split(tmpdata$DBH, tmpdata$uID)
tmp <- unlist(lapply(dbhlist, function(x) sum(is.na(x))!=length(x)))
tmpdata <- tmpdata[tmpdata$uID %in% names(tmp)[tmp==T],]
tmpdata <- droplevels(tmpdata)
dbhlist <- split(tmpdata$DBH, tmpdata$uID)

# GET THE RECRUITMENT YEAR FOR EACH INDIVIDUAL
yearlist <- split(tmpdata$year, tmpdata$uID)
year <- unlist(lapply(dbhlist, function(x) min(which(!is.na(x)))-1))
year.recruit <- rep(NA, length(yearlist))
for(i in 1:length(yearlist)){
	if(year[i] == 0){ year.recruit[i] <- NA	} else { year.recruit[i] <- yearlist[[i]][year[i]]	}
}
names(year.recruit) <- unique(tmpdata$uID)
tmpdata$year.recruit <- year.recruit[match(tmpdata$uID, names(year.recruit))]
tmpdata$year.recruit <- tmpdata$year.recruit+1
tmpdata$recruit <- tmpdata$year == tmpdata$year.recruit
tmpdata$recruit[is.na(tmpdata$recruit)] <- FALSE

# REMOVE DATA FROM FIRST AND LAST CENSUSES
tmp <- tmpdata[! tmpdata$SITE.CENSUS %in% c(firstcensus, lastcensus), ]
tmp <- droplevels(tmp)

# GET SPECIES TRAIT VALUES
trait <- tapply(tmp[,'SLA'], tmp$SPECIES, mean, na.rm=T)

# COMPUTE DEMOGRAPHIC DRIVERS OF CWM SHIFTS

out <- matrix(nrow=length(unique(tmp$SITE.CENSUS)), ncol=14)
for (i in 1:length(unique(tmp$SITE.CENSUS))){

	site <- unique(tmp$SITE.CENSUS)[i]
	d <- tmp[tmp$SITE.CENSUS == site, ]

	# Get initial basal area
	ba0 <- tapply(d$ba[d$recruit == F], d$SPECIES[d$recruit == F], sum, na.rm=T)
	ba0[is.na(ba0)] <- 0
	
	# Get increase of basal area due to growth
	d$ba.g <- dbh2ba(d$DBH + d$growth) - dbh2ba(d$DBH)
	ba.g <- tapply(d$ba.g, d$SPECIES, sum, na.rm=T)
	ba.g[is.na(ba.g)] <- 0

	# Get increase of basal area due to recruitment
	ba.r <- tapply(d$ba[d$recruit == T], d$SPECIES[d$recruit == T], sum, na.rm=T)
	ba.r[is.na(ba.r)] <- 0
	
	# Get decrease of basal area due to mortality
	ba.m <- tapply(d$ba[d$survive %in% 0], d$SPECIES[d$survive %in% 0], sum, na.rm=T)
	ba.m[is.na(ba.m)] <- 0

	x <- as.data.frame(cbind(trait, ba0, ba.g, ba.r, ba.m))
	x <- x[!is.na(x$trait),]

	out[i,] <- c(site, DemographicCWM(x$ba0, x$trait, x$ba.g, x$ba.r, x$ba.m), round(c(sum(x$ba0), sum(x$ba.g), sum(x$ba.r), sum(x$ba.m)), 3))


}
out <- as.data.frame(out)
names(out) <- c('SITE.CENSUS', names(DemographicCWM(x$ba0, x$trait, x$ba.g, x$ba.r, x$ba.m)), 'BA.0','BA.g',"BA.r","BA.m")
for(i in 2:ncol(out)) out[,i] <- as.numeric(as.character(out[,i]))
out$GR <- out$Growth + out$Recruitment

tmp <- cbind(census, out[match(census$SITE2.CENSUS, out$SITE.CENSUS),])

plot(out$BA.0 / (out$BA.0 + out$BA.g + out$BA.r + out$BA.m), out$t0, ylim=c(0,1), pch=16)
points(out$BA.g / (out$BA.0 + out$BA.g + out$BA.r + out$BA.m), out$Growth, col=3, pch=16)
points(out$BA.r / (out$BA.0 + out$BA.g + out$BA.r + out$BA.m), out$Recruitment, col=4, pch=16)
points(out$BA.m / (out$BA.0 + out$BA.g + out$BA.r + out$BA.m), out$Mortality, col=2, pch=16)

summary(lm(out$BA.0 / (out$BA.0 + out$BA.g + out$BA.r + out$BA.m) ~ out$t0))

(out$BA.0 + out$BA.g + out$BA.r + out$BA.m)


rowSums(out[,8:11])



plot(tmp$CENSUS + tmp$AGE00, )


censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$cwm0, col=0, xlab='Years since abandonment', ylab='CWM value', ylim=c(10,20))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$cwm0[census$SITE==site], type='l', col=1)
	points(censusyr[census$SITE==site], tmp$cwm.g[census$SITE==site], type='l', col=3)
	points(censusyr[census$SITE==site], tmp$cwm.r[census$SITE==site], type='l', col=2)
	points(censusyr[census$SITE==site], tmp$cwm.m[census$SITE==site], type='l', col=4)
	points(censusyr[census$SITE==site], tmp$cwm1[census$SITE==site], type='l', col=5)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$cwm0[census$SITE==site], pch=21, bg=1)
			points(censusyr[census$SITE==site], tmp$cwm.g[census$SITE==site], pch=21, bg=3)
			points(censusyr[census$SITE==site], tmp$cwm.r[census$SITE==site], pch=21, bg=4)
			points(censusyr[census$SITE==site], tmp$cwm.m[census$SITE==site], pch=21, bg=2)
			points(censusyr[census$SITE==site], tmp$cwm1[census$SITE==site], pch=21, bg=5)
		}
}



###########################################
#### DO SOME PLOTTING OF THE RESULTS ####
###########################################

########################################################
#### FIT AGE VS. % CONTRIBUTION CURVES ON ONE PLOT ####
########################################################

censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$t0, col=0, xlab='Years since abandonment', ylab='% Contribution to CWM change', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$t0[census$SITE==site], type='l', col=1)
	# points(censusyr[census$SITE==site], tmp$GR[census$SITE==site], type='l', col=3)
	points(censusyr[census$SITE==site], tmp$Growth[census$SITE==site], type='l', col=3)
	points(censusyr[census$SITE==site], tmp$Mortality[census$SITE==site], type='l', col=2)
	points(censusyr[census$SITE==site], tmp$Recruitment[census$SITE==site], type='l', col=4)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$t0[census$SITE==site], pch=21, bg=1)
			# points(censusyr[census$SITE==site], tmp$GR[census$SITE==site], pch=21, bg=3)
			points(censusyr[census$SITE==site], tmp$Growth[census$SITE==site], pch=21, bg=3)
			points(censusyr[census$SITE==site], tmp$Mortality[census$SITE==site], pch=21, bg=2)
			points(censusyr[census$SITE==site], tmp$Recruitment[census$SITE==site], pch=21, bg=4)
		}
}

y <- tmp$t0
x <- censusyr 
dat <- data.frame(y=y[!is.na(y)], x=x[!is.na(y)])
mod <- nls(y ~ Vm * x / (K+x), data=dat, start=list(K=0.5, Vm=1))
summary(mod)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], predict(mod, new), type='l', lwd=3, col=1)

y <- tmp$Growth
x <- censusyr
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=3)
summary(mod)

y <- tmp$Recruitment
x <- censusyr
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y[y>0]) ~ x[y>0])
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=4)
summary(mod)

y <- tmp$Mortality
x <- censusyr
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=2)
summary(mod)


##########################################################
#### FIT AGE VS. % CONTRIBUTION CURVES ON MULTI-PLOT ####
##########################################################
par(mfrow=c(2,2), mar=c(4,4,2,2))

censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$t0, col=0, xlab='Years since abandonment', ylab='Standing Biomass % Contribution', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$t0[census$SITE==site], type='l', col=col[i], lwd=3)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$t0[census$SITE==site], pch=16, col=col[i])
		}
}
y <- tmp$t0
x <- censusyr 
dat <- data.frame(y=y[!is.na(y)], x=x[!is.na(y)])
mod <- nls(y ~ Vm * x / (K+x), data=dat, start=list(K=0.5, Vm=1))
summary(mod)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], predict(mod, new), type='l', lwd=2)


censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$Growth, col=0, xlab='Years since abandonment', ylab='Growth % Contribution', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$Growth[census$SITE==site], type='l', col=col[i], lwd=3)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$Growth[census$SITE==site], pch=16, col=col[i])
		}
}
y <- tmp$Growth
x <- censusyr
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=2, col=1)
summary(mod)


plot(censusyr, tmp$Recruitment, col=0, xlab='Years since abandonment', ylab='Recruitment % Contribution', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$Recruitment[census$SITE==site], type='l', col=col[i], lwd=3)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$Recruitment[census$SITE==site], pch=16, col=col[i])
		}
}
y <- tmp$Recruitment
x <- censusyr
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y[y>0]) ~ x[y>0])
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=2, col=1)
summary(mod)


plot(censusyr, tmp$Mortality, col=0, xlab='Years since abandonment', ylab='Mortality % Contribution', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$Mortality[census$SITE==site], type='l', col=col[i], lwd=3)
		if(sum(census$SITE==site)==3){
			points(censusyr[census$SITE==site], tmp$Mortality[census$SITE==site], pch=16, col=col[i])
		}
}
y <- tmp$Mortality
x <- censusyr
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,50,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=2)
summary(mod)






#############################################################
#### FIT TOTAL BA VS. % CONTRIBUTION CURVES ON ONE PLOT ####
#############################################################

plot(tmp$totalBA, tmp$t0, col=0, xlab='Total Basal Area (m2)', ylab='% Contribution to CWM change', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
		points( tmp$totalBA[census$SITE==site], tmp$t0[census$SITE==site], pch=21, bg=1)
		# points( tmp$totalBA[census$SITE==site], tmp$GR[census$SITE==site], pch=21, bg=3)
		points( tmp$totalBA[census$SITE==site], tmp$Growth[census$SITE==site], pch=21, bg=3)
		points( tmp$totalBA[census$SITE==site], tmp$Mortality[census$SITE==site], pch=21, bg=2)
		points( tmp$totalBA[census$SITE==site], tmp$Recruitment[census$SITE==site], pch=21, bg=4)
}

y <- tmp$t0
x <- tmp$totalBA
dat <- data.frame(y=y[!is.na(y)], x=x[!is.na(y)])
mod <- nls(y ~ Vm * x / (K+x), data=dat, start=list(K=0.5, Vm=1))
summary(mod)
new <- data.frame(x=seq(-5,5,.01))
points(new[,1], predict(mod, new), type='l', lwd=3, col=1)

y <- tmp$Growth
x <- tmp$totalBA
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,5,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=3)
summary(mod)

y <- tmp$Recruitment
x <- tmp$totalBA
x <- x[y>0]; y <- y[y>0]
mod <- lm(log(y[y>0]) ~ x[y>0])
new <- data.frame(x=seq(-5,5,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=4)
summary(mod)

y <- tmp$Mortality
x <- tmp$totalBA
mod <- lm(log(y) ~ x)
new <- data.frame(x=seq(-5,5,.01))
points(new[,1], exp(predict(mod, new)), type='l', lwd=3, col=2)
summary(mod)







tmp2 <- tmp[tmp$Growth>0,]
censusyr <- tmp2$ages+tmp2$CENSUS

yg <- tmp2$Growth / (tmp2$Growth + tmp2$Recruitment + tmp2$Mortality)
yr <- tmp2$Recruitment / (tmp2$Growth + tmp2$Recruitment + tmp2$Mortality)
ym <- tmp2$Mortality / (tmp2$Growth + tmp2$Recruitment + tmp2$Mortality)

plot(censusyr, y, col=0, ylab='', ylim=c(0,1), cex.axis=1.5)
for (i in 1:length(unique(tmp2$SITE))){
	site <- unique(tmp2$SITE)[i]
	points(censusyr[tmp2$SITE==site], yg[tmp2$SITE==site], col=col[i], type='l', lwd=2)
#	points(censusyr[tmp2$SITE==site], yr[tmp2$SITE==site], col=4, type='l')
#	points(censusyr[tmp2$SITE==site], ym[tmp2$SITE==site], col=2, type='l')
#	points(censusyr[tmp2$SITE==site], y[tmp2$SITE==site], col=col[i], type='l')
	}
abline(h=.5, lty=3)






censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$t0, col=0, xlab='Years since abandonment', ylab='% Contribution to CWM change', ylim=c(0,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	points(censusyr[census$SITE==site], tmp$t0[census$SITE==site], pch=21, bg=1)
	points(censusyr[census$SITE==site], tmp$Growth[census$SITE==site], pch=21, bg=3)
	points(censusyr[census$SITE==site], tmp$Mortality[census$SITE==site], pch=21, bg=2)
	points(censusyr[census$SITE==site], tmp$Recruitment[census$SITE==site], pch=21, bg=4)
}


censusyr <- tmp$ages+tmp$CENSUS
plot(censusyr, tmp$t0, col=0, xlab='Years since abandonment', ylab='% Contribution to CWM change', ylim=c(-1,1))
for (i in 1:length(unique(tmp$SITE))){
	site <- unique(census$SITE)[i]
	g <- tmp$Growth[census$SITE==site]
	r <- tmp$Recruitment[census$SITE==site]
	points(censusyr[census$SITE==site], ((g / (r + g))-0.5)*2, type='l', col=col[i], lwd=2)
}




