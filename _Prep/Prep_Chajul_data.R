library(reshape2)
library(RColorBrewer)

setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
data <- read.csv("Chajul_data_4.27.15.csv")
data$HT12 <- as.numeric(as.character(data$HT12))
plots <- read.csv("plot_ages.csv")

#head(data)
#head(plots)

data$uID <- paste(as.character(data$SITE2), as.character(data$TAG), as.character(data$NUMBER), sep='.')
data <- data[!data$TAG=="",]  # DROP RECORDS WITH NO TAG
outs <- names(table(data$uID))[table(data$uID)==2]
data <- data[!data$uID %in% outs,]
data$AGE00 <- plots$AGE00[match(data$SITE, plots$SITE)]

### Change some error in species names
data$SPECIES[data$SPECIES=="trema laxiflora"] <- "Trema laxiflora"
data$SPECIES[data$SPECIES=='Bursera simarouba'] <- 'Bursera simaruba'
data$SPECIES[data$SPECIES=='Guettarda anomala'] <- 'Guatteria anomala'
data$SPECIES[data$SPECIES=='Phylostillum subcesile'] <- 'Phylostilon subsessile'
data$SPECIES[data$SPECIES=='Senna papilosa'] <- 'Senna papillosa'
data <- droplevels(data)

### Melt the DBH data...
d <- data[,!(colnames(data) %in% paste('HT',0:14,sep=''))] 

d <- melt(d, id.vars = c("uID","OLDSITE","SITE2","NUMBER",
                         "YEAR_RECRUIT","TAG","SITE","QUAD",
                         "SUBQUAD","TAG","SPECIES","NOTES",'AGE00'),
          value.name='DBH')
d <- d[order(d$uID),]
d$census <- 0:14

### Melt and add the HT data...
d2 <- data[,!(colnames(data) %in% paste('DIA',0:14,sep=''))] 
d2 <- melt(d2, id.vars = c("uID","OLDSITE","SITE2","NUMBER",
                         "YEAR_RECRUIT","TAG","SITE","QUAD",
                         "SUBQUAD","TAG","SPECIES","NOTES",'AGE00'),value.name='HT')
d2 <- d2[order(d2$uID),]
d$HT <- d2$HT

### FIX HT ERRORS (IDENTIFIED WITH ELBOW GREASE!)
d$HT[d$uID %in% 'RAF1.784.13815' & d$census %in% '6'] <- 640
d$HT[d$uID %in% 'RAF0.714.15758' & d$census %in% '11'] <- 1180
d$HT[d$uID %in% 'PEPED2.47.12134' & d$census %in% '13'] <- 1600
d$HT[d$uID %in% 'PED2.843.11828' & d$census %in% '13'] <- 1100
d$HT[d$uID %in% 'HEC17.91.6451' & d$census %in% '10'] <- 1380
d$HT[d$uID %in% 'HEC17.91.6451' & d$census %in% '11'] <- 1380
d$HT[d$uID %in% 'HEC17.91.6451' & d$census %in% '12'] <- 1450
d$HT[d$uID %in% 'HEC17.253.6725' & d$census %in% '12'] <- 2200
d$HT[d$uID %in% 'HEC17.253.6725' & d$census %in% '13'] <- 2500
d$HT[d$uID %in% 'HEC1.1864.8636' & d$census %in% '13'] <- 1500
d$HT[d$uID %in% 'GUM3.613.5406' & d$census %in% '13'] <- 1600
d$HT[d$uID %in% 'GUM25.268.3616' & d$census %in% '13'] <- 1400
d$HT[d$uID %in% 'FER4.1570.4732' & d$census %in% '10'] <- 1750
d$HT[d$uID %in% 'FER4.1570.4732' & d$census %in% '13'] <- 289
d$HT[d$uID %in% 'PEGUM3.1.6265' & d$census %in% '10'] <- 900
d$HT[d$uID %in% 'SAN8.346.17949' & data$census %in% '13'] <- 1115
d$HT[d$uID %in% 'FER4.231.4429' & d$census %in% '10'] <- 1750
d$HT[d$uID %in% 'FER4.1570.4732' & d$census %in% '10'] <- 230
d$HT[d$uID %in% 'FER1.1097.3062' & d$census %in% '4'] <- 1030
d$HT[d$uID %in% 'HEC1.1505.7989' & d$census %in% '8'] <- 1250
d$HT[d$uID %in% 'HEC17.479.7007' & d$census %in% '10'] <- 1640
d$HT[d$uID %in% 'HEC17.807.7047' & d$census %in% '3'] <- 270
d$HT[d$uID %in% 'PED2.1139.12053' & d$census %in% '13'] <- 330
d$HT[d$uID %in% 'RAF1.108.13136' & d$census %in% '9'] <- 1100


### THESE ARE QUESTIONABLE BUT LOOKS LIKE BREAKAGE:
#d[d$uID %in% 'SAN8.346.17949' & d$census %in% '11',]
#d[d$uID %in% 'HEC1.1155.8818' & d$census %in% '13',]
#d[d$uID %in% 'PED2.881' & d$census %in% 'PED2.8',]
#d[d$uID %in% 'PEFER4.32.4770' & d$census %in% 'PEFER4.11',]

### ADD IN CENSUS INFO
census <- read.csv('site.censuses.csv')
censuses <- as.character(census$SITE.CENSUS[census$TAKEN=="TRUE"])
census$year <- substring(as.character(as.Date(as.character(census$DATE), format='%m/%d/%y')),1,4)
head(census)

### PLOT CENSUS VERSUS AGE
#census$DATE <- as.Date(census$DATE, format='%m/%d/%y')
#plot(census$DATE[census$PE==F], census$AGE00[census$PE==F]+census$CENSUS[census$PE==F], col=0, axes=F, ylab="Years Post Abandonment", xlab="Calendar Year", ylim=c(0,30))	
#for(i in 1:length(unique(census$SITE2))){
#	site <- census[census$SITE2 == unique(census$SITE2)[i],]
#	site <- site[site$PE==FALSE,]
#	col <- c(brewer.pal(12,"Set3"),brewer.pal(8,"Set2"), brewer.pal(9,"Set1"))
#	points(site$DATE, jitter(site$AGE00+site$CENSUS), type='l', col=col[i], lwd=2, pch=21)
#	points(site$DATE, jitter(site$AGE00+site$CENSUS), type='p', bg=col[i], pch=21)
#}
#axis(1, labels=1999:2014, at=as.Date(paste(1,1,2000:2015,sep='/'), format='%m/%d/%Y'), las=2)
#axis(2)

d$SITE.CENSUS <- paste(d$SITE2, d$census, sep='.')
census$DATE <- as.Date(census$DATE, format='%m/%d/%y')
d$DATE <- census$DATE[match(d$SITE.CENSUS, as.character(census$SITE.CENSUS))]

d$growth <- NA
d$ht.growth <- NA
d$YEAR_RECRUIT2 <- NA
d$survive <- NA
d$year <- 2000:2014
d$status <- 	ifelse(!d$SITE.CENSUS %in% censuses, 'unknown', ifelse(is.na(d$DBH), 'dead', 'alive'))

### Correct for problematic measurements with straightforward error (determined by hand)...
fix <- read.csv("fix_diams.csv")
d$DBH[paste(d$uID, d$variable) %in% paste(fix$uID, fix$variable)] <- fix$NEW
me <- read.csv("remove_for_measurement_error.csv")
d <- d[!d$uID %in% me$uID,]

d <- d[order(d$census, d$uID),]
g <- vector()
htg <- vector()
s <- vector()
int <- vector()
for(census in 0:14){
	c1 <- d[d$census == census,]
	c2 <- d[d$census == (census+1),]
	g <- c(g, c2$DBH - c1$DBH)
	htg <- c(htg, c2$HT - c1$HT)
  if(census == 14){
    g <- c(g, rep(NA, times=nrow(c1)))
    htg <- c(htg, rep(NA, times=nrow(c1)))
  }
	s <- c(s, ifelse(c1$status=='alive' & c2$status=='alive', 1, 
					ifelse(c1$status=='alive' & c2$status=='dead', 0, NA)))
	if(census == 14){
    s <- c(s, rep(NA, times=nrow(c1)))
	}
	int <- c(int, c2$DATE - c1$DATE)
	if(census == 14){
    int <- c(int, rep(NA, times=nrow(c1)))
	}
}

d$growth <- g
d$ht.growth <- htg
d$survive <- s
d$int <- int
d <- d[order(d$uID),]

tmp <- d[d$SITE.CENSUS %in% censuses,]

names(tmp)
tmp2 <- tmp[,!colnames(tmp) %in% c("OLDSITE","NUMBER","YEAR_RECRUIT","TAG","NOTES","variable") ]
head(tmp2)

gdata <- tmp2[!is.na(tmp2$growth),]
gdata <- gdata[order(gdata$uID),]

# # par(mfrow=c(2,2))
# # hist(gdata$growth, xlim=c(-5,5), breaks=1000)
# # abline(v=0,col=2,lwd=3)
# # hist(gdata$growth, breaks=1000)
# # abline(v=0,col=2,lwd=1)

# checks <- gdata$uID[gdata$growth < -5]
# checks <- c(gdata$uID[gdata$growth > 5], checks)
# checks <- unique(checks)
# write.csv(d[d$uID %in% checks,], file='check.growth_NEW.csv')
# write.csv(d[d$uID %in% checks,], file='check.growth.csv')


# PERHAPS SOME MORE GROWTH OUTLIER CLEANING...
# tapply(gdata$growth, gdata$SPECIES, sd, na.rm=T)* 10

# # head(gdata)
# # plot(jitter(gdata$AGE00+gdata$census), gdata$growth, ylim=c(-2,5), pch=16, col=grey(.7,.1), xlab='Years since abandonment', ylab='Diameter growth (cm / yr)')
# # points(1:30, tapply(gdata$growth, (gdata$AGE00+gdata$census), mean), pch=21, bg=5, cex=1.5)
# # abline(h=0,col=2,lwd=1)

# # plot(1:30, tapply(gdata$growth, (gdata$AGE00+gdata$census), mean), pch=21, bg=5, cex=1.5, ylim=c(-2,5), xlab='Years since abandonment', ylab='Diameter growth (cm / yr)')
# # abline(h=0,col=2,lwd=1)
# # for(i in 1:30)	{
	# # x <- subset(gdata, (gdata$AGE00+gdata$census)==i)
	# # y <- tapply(x$growth, x$SPECIES, mean, na.rm=T)
	# # y <- y[!is.na(y)]
	# # points(jitter(rep(i, length(y)),.5), y, pch=16, col=grey(.5,.5), cex=.75) 
# # }
# # points(1:30, tapply(gdata$growth, (gdata$AGE00+gdata$census), mean), pch=21, bg=5, cex=1.5)


# REMOVE SPECIES WITH NO GROWTH
x <- is.na(tapply(!is.na(gdata$growth), gdata$SPECIES, sum)==0)
gdata <- gdata[!gdata$SPECIES %in% names(x)[x==T],]
gdata <- droplevels(gdata)

gout <- vector()
for(i in 1:30)	{
	x <- subset(gdata, (gdata$AGE00+gdata$census)==i)
	y <- tapply(x$growth, x$SPECIES, mean, na.rm=T)
	gout <- cbind(gout, y)
}

abund <- log(tapply(gdata$growth, gdata$SPECIES, length)/50)+.2

# # plot(1:30,gout[1,],ylim=c(-2,5), type='l', col=grey(0,.2), lwd=abund[1], xlab='Years since abandonment', ylab='Diameter growth (cm / yr)')
# # for(i in 2:nrow(gout)){
# # points(1:30,gout[i,], type='l',col=grey(0,.2), lwd=abund[i])
# # }
# # abline(h=0,lty=2,col=2,lwd=2)

head(tmp2)
tmp2$ba <- 0.00007854 * (tmp2$DBH^2)
ba <- tapply(tmp2$ba, tmp2$SITE.CENSUS, sum, na.rm=T)
table(tmp2$AGE00+tmp2$census)


census <- read.csv('site.censuses.csv')
census <- census[census$TAKEN==TRUE,]
census$year <- substring(as.character(as.Date(as.character(census$DATE), format='%m/%d/%y')),1,4)
census$DATE <- as.Date(census$DATE, format='%m/%d/%y')
census$ba <- ba[match(as.character(census$SITE.CENSUS), names(ba))]
head(census)

ages <- tapply(tmp2$AGE00, tmp2$SITE2, mean)

census$ages <- ages[match(as.character(census$SITE), names(ages))]

# # col <- c(brewer.pal(9,"Set1")[c(1:5,7:9)],brewer.pal(8,"Set2"), brewer.pal(12,"Set3"))
# # plot(census$ages+census$CENSUS, census$ba, type='l', col=0, xlab='Years since abandonment', ylab='Total basal area (m^2)')
# # for(i in 1:length(unique(census$SITE))){
	# # x <- census[census$SITE == unique(census$SITE)[i],]
	# # print(x$ages+x$CENSUS)
	# # points(x$ages+x$CENSUS, x$ba, type='l', col=col[i], lwd=3)
# # }


stemdens <- tapply(tmp2$SITE[tmp2$status=='alive'], tmp2$SITE.CENSUS[tmp2$status=='alive'], length)
census$stemdens <- stemdens[match(census$SITE.CENSUS, names(stemdens))]

# # plot(census$ages+census$CENSUS, census$stemdens, col=0, xlab='Years since abandonment', ylab='Stems per plot')
# # col <- c(brewer.pal(9,"Set1")[c(1:5,7:9)],brewer.pal(8,"Set2"), brewer.pal(12,"Set3"))
# # for (i in 1:length(unique(census$SITE))){
	# # site <- unique(census$SITE)[i]
	# # points((census$ages+census$CENSUS)[census$SITE==site], census$stemdens[census$SITE==site], type='l', col=col[i], lwd=3)	
# # }

data <- tmp2

setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
#save(data, file="Chajul_data_processed_notraits_11.20.15.RDA")
#save(census, file="Chajul_census_processed_11.20.15.RDA")
#save(data, file="Chajul_data_processed_notraits_4.27.15.RDA")
#save(census, file="Chajul_census_processed_8.25.15.RDA")
load("Chajul_data_processed_notraits_11.20.15.RDA")
head(data)

## READ IN TRAIT DATA
tdata <- read.csv("Mean_traits_4.27.15.csv")
wood <- read.csv("wood_traits_4.27.15.csv")
tdata$WDMC <- wood$WDMC[match(tdata$species, wood$species)]
data <- cbind(data, tdata[match(data$SPECIES, tdata$species),])

#par(mfrow=c(4,4), mar=c(2,2,2,2))
#for(t in 2:16) {hist(tdata[,t], main=names(tdata)[t])}
#quartz()
#par(mfrow=c(4,4), mar=c(2,2,2,2))
#for(t in 2:16) {hist(log(tdata[,t]), main=names(tdata)[t])}

data$log.LA <- log(data$LA)
data$log.SLA <- log(data$SLA)
data$log.LDMCsimple <- log(data$LDMCsimple)
data$log.LT <- log(data$LT)
data$log.F2P <- log(data$F2P)
data$log.sFtP <- log(data$sFtP)
data$log.Ft <- log(data$Ft)
data$log.P <- log(data$P)
data$log.N <- log(data$N)

#save(data, file="Chajul_data_processed_wtraits_11.20.15.RDA")
#save(data, file="Chajul_data_processed_wtraits_4.27.15.RDA")

### ADD CHAVE WD AND GENUS-LEVEL AVERAGE WD DATA... 
setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
load("Chajul_data_processed_wtraits_11.20.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")
tdata <- read.csv('Mean_traits_4.27.15.csv')

nowd <- sort(as.character(unique(data$SPECIES[is.na(data$WD)])))
chave <- read.csv('ChaveWoodDensity.csv')
chave$sp <- paste(chave$GENUS, chave$SPECIES)
chavetmp <- chave[chave$sp %in% nowd,c('sp','WSG')]
newwd <- chavetmp$WSG[match(data$SPECIES, chavetmp$sp)]
data$all.WD <- ifelse(!is.na(data$WD), data$WD, newwd)

wd <- tapply(data$all.WD, data$SPECIES, mean)
genus.list <- unlist(lapply(strsplit(names(wd), " "), function(x) x[[1]]))
genus.wd <- tapply(wd, genus.list, mean, na.rm=T)
data$genus <- unlist(lapply(strsplit(as.character(data$SPECIES), " "), function(x) x[[1]]))
newwd <- genus.wd[match(data$genus, names(genus.wd))]
data$all.WD <- ifelse(!is.na(data$all.WD), data$all.WD, newwd)

# PROBABLY SOME MORE HEIGHT ERRORS TO CLEAN BUT GOOD FOR NOW...
range(data$ht.growth,na.rm=T)
hist(data$ht.growth, breaks=1000)
checks <- data$uID[!is.na(data$ht.growth) & data$ht.growth<(-500)]
checks <- c(checks, data$uID[!is.na(data$ht.growth) & data$ht.growth>(500)])
unique(checks)

### CALCULATE AGB BASED ON CHAVE 2014 AND FULL WD DATA
data$agb <- 0.0673 * ((data$all.WD * data$DBH^2 * data$HT/100)^0.976) 
  
#save(data, file="Chajul_data_processed_wtraits_11.20.15.RDA")
load("Chajul_data_processed_wtraits_11.20.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")


###########################
### END DATA PROCESSING ###
###########################

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #

##################################
### START EXPLORATORY PLOTTING ###
##################################

par(mfrow=c(4,4), mar=c(2,2,2,2))
col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
col <- c(col, col)
for (t in 1:16){
trait <- names(data)[c(23,28:32, 35, 36:44)][t]			# Use log traits
# trait <- names(tdata)[-1][t]			# Don't use log traits
tmpcwm <- tapply(data[,trait][data$status=='alive'], data$SITE.CENSUS[data$status=='alive'], mean, na.rm=T)
census$cwm <- tmpcwm[match(as.character(census$SITE.CENSUS), names(tmpcwm))]
plot(census$ages+census$CENSUS, census$cwm, col=0, xlab='Years since abandonment', ylab='Stems per plot', main=trait)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points((census$ages+census$CENSUS)[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=1)
	}
}

###  PLOT SLA, WD
par(mfrow=c(2,2), mar=c(2,2,2,2))
col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
col <- c(col, col)

data$LMA <- 1/data[,'SLA'] * 1000
data$log.LMA <- log(data$LMA)

for (t in 1:4){
trait <- c('WD', 'LMA', 'Chlorophyll', 'log.F2P')[t]
tmpcwm <- tapply(data[,trait][data$status=='alive'], data$SITE.CENSUS[data$status=='alive'], mean, na.rm=T)
census$cwm <- tmpcwm[match(as.character(census$SITE.CENSUS), names(tmpcwm))]
y <- census$ages+census$CENSUS
# y <- census$ba
plot(y, census$cwm, col=0, xlab='Years since abandonment', ylab='Stems per plot', main=trait)
#plot(, census$cwm, col=0, xlab='Years since abandonment', ylab='Stems per plot', main=trait)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points(y[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=3)
#	points(y[census$SITE==site], census$cwm[census$SITE==site], col=1, pch=16)
	}
}



head(data)
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
head(data)


#### TOTAL BASAL AREA INCREASES WITH SUCCESSION ####
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
plot(census$ages+census$CENSUS, census$totalBA, col=0, xlab='Years since abandonment', ylab='Total BA')
col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
col <- c(col, col)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points((census$ages+census$CENSUS)[census$SITE==site], census$totalBA[census$SITE==site], type='l', col=col[i], lwd=3)	
}
######################################################


quartz()
par(mfrow=c(4,4), mar=c(2,2,2,2))
for (t in 1:16){
	trait <- names(tdata)[-1][t]
	tmpcwm <- data[,trait] * data$relBA
	cwm <- tapply(tmpcwm[data$status=='alive'], data$SITE.CENSUS[data$status=='alive'], sum, na.rm=T)
	census$cwm <- cwm[match(as.character(census$SITE.CENSUS), names(cwm))]
	plot(census$ages+census$CENSUS, census$cwm, col=0, xlab='Years since abandonment', ylab='Trait value', main=trait)
	col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
	col <- c(col, col)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points((census$ages+census$CENSUS)[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=3)
	}
}



#### USE LOG TRAITS WHEN NON NORMAL ####
quartz()
par(mfrow=c(4,4), mar=c(2,2,2,2))
for (t in 1:16){
	trait <- names(data)[c(23,28:32, 35, 36:44)][t]
	tmpcwm <- data[,trait] * data$relBA
	cwm <- tapply(tmpcwm[data$status=='alive'], data$SITE.CENSUS[data$status=='alive'], sum, na.rm=T)
	census$cwm <- cwm[match(as.character(census$SITE.CENSUS), names(cwm))]
	plot(census$ages+census$CENSUS, census$cwm, col=0, xlab='Years since abandonment', ylab='Trait value', main=trait)
	col <- c(brewer.pal(12,"Set3")[c(1,3:8,10:12)],brewer.pal(8,"Set2"), brewer.pal(9,"Set1")[c(1:5,7:8)])
	col <- c(col, col)
for (i in 1:length(unique(census$SITE))){
	site <- unique(census$SITE)[i]
	points((census$ages+census$CENSUS)[census$SITE==site], census$cwm[census$SITE==site], type='l', col=col[i], lwd=3)
	}
}
###########################################





head(data)
data$SITE.CENSUS.QUAD <- paste(data$SITE.CENSUS, data$QUAD, sep=".")



library(FD)

# rownames(tdata) <- tdata$species

trt <- cbind(log(tdata[,c(2,3,5,14,15)]), tdata[,c(4,7,8,9,10,11,12,13,16)]) 
rownames(trt) <- tdata$species

fdiv <- list()
# fdiv.quad <- list()

for(t in 1:ncol(trt)){
	trait <- colnames(trt)[t]
	message(paste('Now working on', trait))
	x <- tdata[,trait]
	names(x) <- rownames(tdata)
	x <- x[!is.na(x)]

	mat <- table(data$SITE.CENSUS[data$status=='alive'], data$species[data$status=='alive'])
	mat <- as.data.frame.matrix(mat)
	mat <- as.matrix(mat)
	mat <- mat[,colnames(mat) %in% names(x)]
	mat <- mat[,colSums(mat) > 0]
	mat <- mat[rowSums(mat) > 0,]

	# qmat <- table(data$SITE.CENSUS.QUAD[data$status=='alive'], data$species[data$status=='alive'])
	# qmat <- as.data.frame.matrix(qmat)
	# qmat <- as.matrix(qmat)
	# qmat <- qmat[,colnames(qmat) %in% names(x)]
	# qmat <- qmat[,colSums(qmat) > 0]
	# qmat <- qmat[rowSums(qmat) > 0,]

	x <- x[names(x) %in% colnames(mat)]
	x <- x[match(colnames(mat), names(x))]

	# qx <- x[names(x) %in% colnames(qmat)]
	# qx <- x[match(colnames(qmat), names(x))]

	fdiv[[t]] <- dbFD(x, mat, corr="cailliez")
	# fdiv[[t]] <- dbFD(x, mat)
	# fdiv.quad[[t]] <- dbFD(qx, qmat)
	names(fdiv)[t] <- trait
}


####################################
### IF DOING THE QUAD ANALYSIS...
# expand census to have 125 rows per site.census (one for each quad.site.census)
qcensus <- census[rep(seq_len(nrow(census)), each=125), ]
qcensus$quad <- 1:125
qcensus$SITE.CENSUS.QUAD <- paste(qcensus$SITE.CENSUS, qcensus$quad, sep='.')
qcensus$SITE.QUAD <- paste(qcensus$SITE, qcensus$quad, sep='.')

		m=2
		vals <- list(res$FRic, res$FEve, res$FDis, res$RaoQ, res$CWM)[[m]]
		if(m != 5) qcensus$y <- vals[match(as.character(qcensus$SITE.CENSUS.QUAD), names(vals))]
		if(m == 5) qcensus$y <- vals[match(as.character(qcensus$SITE.CENSUS.QUAD), rownames(vals)),]
		qcensus <- qcensus[!is.na(qcensus$y),]
		plot(qcensus$ages + qcensus$CENSUS, qcensus$y, col=0, xlab='Years since abandonment', ylab="metric")
		rcols <- rainbow(length(unique(qcensus$SITE)), alpha=0.5)
			for (i in 1:length(unique(qcensus$SITE.QUAD))){
				quad <- unique(qcensus$SITE.QUAD)[i]
				points((qcensus$ages+qcensus$CENSUS)[qcensus$SITE.QUAD==quad], qcensus$y[qcensus$SITE.QUAD==quad], type='l', col=rcols[as.factor(qcensus$SITE)[i]], lwd=1)
			}
####################################


setwd("~/Desktop")
pdf(file="Functional_Diversity.pdf", onefile=T)
for (t in 1:ncol(trt)){
#	quartz()
	par(mfrow=c(2,3), mar=c(2,4,1,1), oma=c(20,1, 2, 2))
	for (m in 1:5){
		metric <- c('FRic','FEve','FDis','RaoQ','CWM')[m]
		vals <- list(fdiv[[t]]$FRic, fdiv[[t]]$FEve, fdiv[[t]]$FDis, fdiv[[t]]$RaoQ, fdiv[[t]]$CWM)[[m]]
		if(m != 5) y <- vals[match(as.character(census$SITE.CENSUS), names(vals))]
		if(m == 5) y <- vals[match(as.character(census$SITE.CENSUS), rownames(vals)),]
		plot(census$ba, y, col=0, xlab='Years since abandonment', ylab=metric)
#		plot(census$ages+census$CENSUS, y, col=0, xlab='Years since abandonment', ylab=metric)
			for (i in 1:length(unique(census$SITE))){
				site <- unique(census$SITE)[i]
				points(census$ba[census$SITE==site], y[census$SITE==site], type='l', col=col[i], lwd=1)
#				points((census$ages+census$CENSUS)[census$SITE==site], y[census$SITE==site], type='l', col=col[i], lwd=1)
			}
}
mtext(colnames(trt)[-1][t], line=16)
}
dev.off()














head(data)



census$LMA.WD.corr <- NA
for (i in 1:length(unique(data$SITE.CENSUS))){
	focsite <- unique(data$SITE.CENSUS)[[i]]
	tmpdat <- data[data$SITE.CENSUS %in% focsite,]
	tmpdat <- tmpdat[tmpdat$status=='alive',]
	splist <- as.character(unique(tmpdat$SPECIES))
	tmptdat <- tdata[tdata$species %in% splist,]
	if(nrow(tmptdat) > 0) census$LMA.WD.corr[i] <- cor(1/tmptdat$SLA, tmptdat$WD, use='c')
}

y <- census$LMA.WD.corr






head(tdata)

trt <- cbind(log(tdata[,c(2,3,5,14,15)]), tdata[,c(4,7,8,9,10,11,12,13,16)]) 
rownames(trt) <- tdata$species
trtcorr <- list()
for (i in 1:length(unique(data$SITE.CENSUS))){
	focsite <- unique(data$SITE.CENSUS)[[i]]
	tmpdat <- data[data$SITE.CENSUS %in% focsite,]
	tmpdat <- tmpdat[tmpdat$status=='alive',]
	tmptdat <- trt[rownames(trt) %in% as.character(unique(tmpdat$SPECIES)),]
	if(nrow(tmptdat) > 2) trtcorr[[i]] <- cor(tmptdat, use='pairwise.complete.obs')
	if(nrow(tmptdat) < 3) trtcorr[[i]] <- cor(trt, use='pairwise.complete.obs') * NA
}



colnames(trtcorr[[1]])
y <- unlist(lapply(trtcorr, function(x) x["WD", "SV"]))
ba <- census$ba
age <- census$ages+census$CENSUS
ba <- ba[census$SITE.CENSUS %in% data$SITE.CENSUS]
age <- age[census$SITE.CENSUS %in% data$SITE.CENSUS]
par(mfrow=c(2,2))
main <- round(summary(lm(y ~ age))$r.sq, 2)
plot(age, y, col=0, xlab='Years since abandonment', main=main)
	for (i in 1:length(unique(census$SITE))){
				site <- unique(census$SITE)[i]
				points(age[census$SITE==site], y[census$SITE==site], type='l', col=col[i], lwd=2)
}
abline(lm(y ~ age), lwd=3)
abline(h=0, lty=3)
main <- round(summary(lm(y ~ ba))$r.sq, 2)
plot(ba, y, col=0, xlab='stand basal area', main=main)
	for (i in 1:length(unique(census$SITE))){
				site <- unique(census$SITE)[i]
				points(ba[census$SITE==site], y[census$SITE==site], type='l', col=col[i], lwd=2)
}
abline(lm(y ~ ba), lwd=3)
abline(h=0, lty=3)


# What functional strategies are arguably correlated but their decoupling could reflect niche partitioning?
# Trait-gradient analysis to a successional gradient
# How does the strength of trait correlations of co-occurring species change during succession?
#		- a higher decree of niche partitioning in odler forests may be reflected by lower correlation across species trait values




head(census)








