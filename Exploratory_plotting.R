### EXPLORATORY PLOTS FOR CHAJUL DATA

library(plyr)
Rescale01 <- function(x){
   (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

sdata <- data

###########################
# PLOT AVERAGE SURVIVAL BY CALENDER YEAR
head(sdata)
sdata$age <- sdata$AGE00 + sdata$census

df <- ddply(sdata, c("SITE2","year"), function(df)mean(df$survive))
df <- df[df$SITE2 %in% names(table(df$SITE2)[table(df$SITE2)>1]),]
df <- droplevels(df)
#col <- c(brewer.pal(9,"Set1")[c(1:5,7:9)],brewer.pal(8,"Set2"))
col <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))
plot(df$year, df$V1, type='l', col=0, xlab='Year', ylab='mean survival', ylim=c(0.4,1))
for(i in 1:length(unique(df$SITE2))){
	tmp <- df[df$SITE2 == unique(df$SITE2)[i],]
	lty <- ifelse(substring(as.character(unique(tmp$SITE2)),1,2)=='PE', 2,1)
	points(tmp$year, tmp$V1, type='l', col=col[i], lwd=3, lty=lty)
}
points(1998:2013,Rescale01(map), pch=21, col=rgb(0,.5,1,.4), lwd=10, type='l')
points(map)


###########################
### PLOT AVERAGE SURVIVAL BY FOREST AGE
df <- ddply(sdata, c("SITE2","age"), function(df)mean(df$survive, na.rm=T))
df <- df[df$SITE2 %in% names(table(df$SITE2)[table(df$SITE2)>1]),]
df <- droplevels(df)
col <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))
plot(df$age, df$V1, type='l', col=0, xlab='Years Post Abandonment', ylab='mean survival', ylim=c(0.4,1))
for(i in 1:length(unique(df$SITE2))){
	tmp <- df[df$SITE2 == unique(df$SITE2)[i],]
	lty <- ifelse(substring(as.character(unique(tmp$SITE2)),1,2)=='PE', 1,1)
	points(tmp$age, tmp$V1, type='l', col=col[i], lwd=3, lty=lty)
}
points(tapply(sdata$survive, sdata$age, mean, na.rm=T), type='l', lwd=20, col=grey(.5,.75))
# add points to show relative sample size
# cex <- Rescale01 (tapply(gdata$growth, gdata$age, length))*2+.5
# points(tapply(sdata$survive, sdata$age, mean), cex=cex, pch=21, bg=1)






######################################################
######################################################
### LOOK AT SPECIES AVERAGE GROWTH AND SURVIVAL, RELATE TO TRAITS

head(data)
sdata <- data
sdata$age <- sdata$AGE00 + sdata$census

head(sdata)
sdf <- ddply(sdata, c("age","SPECIES"), function(df)mean(df$survive, na.rm=T))
gdf <- ddply(sdata, c("age","SPECIES"), function(df)mean(df$growth, na.rm=T))
sdf <- droplevels(sdf)
gdf <- droplevels(gdf)

sdf <- sdf[!is.nan(sdf$V1),]
gdf <- gdf[!is.nan(gdf$V1),]
# sdf$wd <- tdata$WD[match(sdf$SPECIES, rownames(tdata))]
# gdf$wd <- tdata$WD[match(gdf$SPECIES, rownames(tdata))]


pdf(file='test.pdf', onefile=T)
for (i in 2:ncol(tdata)){
	sdf$wd <- tdata[,i][match(sdf$SPECIES, rownames(tdata))]
	gdf$wd <- tdata[,i][match(gdf$SPECIES, rownames(tdata))]

# sdf$wd <- tdata$WD[match(sdf$SPECIES, rownames(tdata))]
# gdf$wd <- tdata$WD[match(gdf$SPECIES, rownames(tdata))]

sn <- table(sdata$SPECIES[!is.na(sdata$survive)])
sn <- sn[sn>0]
gn <- table(sdata$SPECIES[!is.na(sdata$growth)])
gn <- gn[gn>0]

ys <- tapply(sdf$V1, sdf$SPECIES, mean, na.rm=T)
yg <- tapply(gdf$V1, gdf$SPECIES, mean, na.rm=T)
wds <- tapply(sdf$wd, sdf$SPECIES, mean)
wdg <- tapply(gdf$wd, gdf$SPECIES, mean)
ns <- sn[match(names(y), names(sn))]
ng <- gn[match(names(y), names(gn))]
sdat <- data.frame(ys, wds, ns)
gdat <- data.frame(yg, wdg, ng)
sdat <- sdat[!is.na(rowSums(sdat)),]
gdat <- gdat[!is.na(rowSums(gdat)),]

trait <- colnames(tdata)[i]

par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(sdat$wd, sdat$y, cex=sqrt(sdat$ns)/20, pch=21, bg=grey(.5,.5), ylab='species average survival', xlab=trait)
#abline(lm(ys ~ wds, weights=ns, data=sdat))

plot(gdat$wd, gdat$y, cex=sqrt(gdat$ng)/20, pch=21, bg=grey(.5,.5), ylab='species average growth', xlab= trait)
#abline(lm(yg ~ wdg, weights=ng, data=gdat))

col <- rainbow(101)[round((Rescale01(sdat$wd)*100)+1, 0)]
plot(gdat$yg, sdat$ys, cex=sqrt(sdat$ns)/20, pch=21, bg=col, ylab='average survival', xlab='average growth')
summary(lm(sdat$ys ~ gdat$yg, weights=sdat$ns))
#abline(lm(sdat$ys ~ gdat$yg, weights=sdat$ns))
}
dev.off()

cor(tdata[,-1], use='complete.obs')



res <- data.frame()
for (i in 1:length(unique(sdata$SITE.CENSUS))){
	site <- unique(sdata$SITE.CENSUS)[i]
	tmp <- sdata[sdata$SITE.CENSUS == site,]
	tmp <- tmp[tmp$status=='alive',]
	tmp2 <- tmp[!duplicated(tmp$species),]
	if(nrow(tmp) > 0){


cor(tmp[,20:35], use='c')[lower.tri(cor(tmp[,20:35], use='c'))]




names(tmp[,20:35])




sla.wd.ab <- cor(tmp$log.SLA, tmp$WD, use='c')
sla.wd.pa <- cor(tmp2$log.SLA, tmp2$WD, use='c')

chl.sla.ab <- cor(tmp$Chlorophyll, tmp$log.SLA, use='c')
chl.sla.pa <- cor(tmp2$Chlorophyll, tmp2$log.SLA, use='c')

f2p.sla.ab <- cor(tmp$log.F2P, tmp$log.SLA, use='c')
f2p.sla.pa <- cor(tmp2$log.F2P, tmp2$log.SLA, use='c')

ln.sla.ab <- cor(tmp$log.N, tmp$log.SLA, use='c')
ln.sla.pa <- cor(tmp2$log.N, tmp2$log.SLA, use='c')

thk.sla.ab <- cor(tmp$log.LT, tmp$log.SLA, use='c')
thk.sla.pa <- cor(tmp2$log.LT, tmp2$log.SLA, use='c')

age <- unique(tmp$age)
site <- strsplit(unique(tmp$SITE.CENSUS), '\\.')[[1]][1]
stems <- nrow(tmp)
rich <- length(unique(tmp$SPECIES))

tmp <- data.frame(sla.wd.ab, sla.wd.pa, chl.sla.ab, chl.sla.pa, f2p.sla.ab, f2p.sla.pa, ln.sla.ab, ln.sla.pa, thk.sla.ab, thk.sla.pa, age, site, stems, rich)
res <- rbind(res, tmp)
}}


par(mfrow=c(2,2))

plot(res[,c(11,5)]); abline(h=0, lty=2)

plot(res[,c(3,2)])
abline(h=0, lty=2)


pdf(file='trait_corrs.pdf', onefile=T)
par(mfrow=c(2,2), mar=c(4,4,2,2))
for(pair in 1:10){
plot(res[,c(11,pair)], col=NA)
col <- rainbow(length(unique(res$site)))
for (i in 1:length(unique(res$site))){
	site <- unique(res$site)[i]
	tmp <- res[res$site==site,]
	points(tmp[,11], tmp[,pair], type='l', col=col[i])
	if(i==length(unique(res$site))) abline(h=0, lty=2)
}}
dev.off()





plot(res[,c(3,2)], col=NA)
col <- rainbow(length(unique(res[,4])))
for (i in 1:length(unique(res[,4]))){
	site <- unique(res[,4])[i]
	tmp <- res[res[,4]==site,]
	points(tmp[,3], tmp[,2], type='l', col=col[i])
}
abline(h=0, lty=2)


summary(lm(res[,2] ~ res[,3]))
abline(lm(res[,2] ~ res[,3]))



i=1
mods <- data.frame()
for (i in 1:length(unique(sdata$SITE.CENSUS))){
	site <- unique(sdata$SITE.CENSUS)[i]
	tmp <- sdata[sdata$SITE.CENSUS == site,]
	tmp <- tmp[tmp$status=='alive',]
	tmp <- tmp[!is.na(tmp$growth),]
	if(nrow(tmp) > 1){
		dcwm <- abs(tmp$WD - sum(tmp$relBA * tmp$WD, na.rm=T))
		mods <- rbind(mods, coef(lm(tmp$growth ~ tmp$WD + dcwm)))
		}
}

names(mods) <- c('int','wd','dcwm')
head(mods)

plot(mods$dcwm, mods$wd); abline(h=0)


	plot(dcwm, tmp$growth)
	summary(lm(tmp$growth ~ dcwm))
	abline(lm(tmp$growth ~ dcwm))

	plot(tmp$WD, tmp$growth)
	abline(lm(tmp$growth ~ tmp$WD + dcwm))

# How does performance vary as a function of trait values and the difference between the species trait and the traits of the local community?

#
plot(jitter(sdata$growth[-nrow(sdata)]), jitter(sdata$survive[-1]))


# Do neighborhood interactions become increasingly important determinants of composition as succession proceeds?
#  How does the strength of trait-performance relationships vary during succession?
	# In early succession, resources are abundant and traits that provide fast growth should do so unencoumbered by local interactions and limited resources
	# Later in succession, light becomes increasingly limited and the link between traits and growth will be weaker.

# Do priority effects explain variation in individual site trajectories?
# Do landscape conditions explain variation in trajectories?
# How does functional hypervolume vary during succession?

# Is functional space filled or expanded during succession?

# How can we elaborate those patterns using individual survival/growth rates?


S <- as.matrix(tapply(sdata$survive, sdata$SPECIES, mean, na.rm=T))
head(tdata)
tdata$survive <- S[match(rownames(tdata), rownames(S))]


head(tdata)
plot(tdata$Chlorophyll,(tdata$survive))



library(hypervolume)


head(tdata)

focsp <- as.character(sdata$SPECIES[sdata$SITE.CENSUS == unique(sdata$SITE.CENSUS)[1]])
tmp <- match(focsp, rownames(tdata))
tmp <- tmp[!is.na(tmp)]
dat <- scale(tdata[tmp, c('WD','SLA','P','F2P')])
hv1 <- hypervolume(dat, bandwidth=.5)

focsp <- as.character(unique(sdata$SPECIES[sdata$SITE.CENSUS == unique(sdata$SITE.CENSUS)[1]]))
dat <- scale(tdata[rownames(tdata) %in% focsp, c('WD','SLA','P','F2P')])
hv1 <- hypervolume(dat, bandwidth=1)
plot(hv1)



head(sdata)
focsp <- as.character(unique(sdata$SPECIES[sdata$age == 1]))
dat <- scale(tdata[rownames(tdata) %in% focsp, c('WD','SLA','P','F2P')])
hv1 <- hypervolume(dat, bandwidth=1)
plot(hv1)


head(sdata)
focsp <- as.character(unique(sdata$SPECIES[sdata$age == 10]))
dat <- scale(tdata[rownames(tdata) %in% focsp, c('WD','SLA','P','F2P')])
hv2 <- hypervolume(dat, bandwidth=1)
plot(hv2)



plot(hv2@RandomUniformPointsThresholded[,1:2])
points(hv2@Data[,1:2], col=2, pch=16)


######################################################
######################################################






###########################
### PLOT AVERAGE GROWTH BY FOREST AGE
head(gdata)
gdata$age <- gdata$AGE00 + gdata$census
df <- ddply(gdata, c("SITE2","age"), function(df)mean(df$growth))
df <- df[df$SITE2 %in% names(table(df$SITE2)[table(df$SITE2)>1]),]
df <- droplevels(df)
#col <- c(brewer.pal(9,"Set1")[c(1:5,7:9)],brewer.pal(8,"Set2"))
col <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Set3"))
plot(df$age, df$V1, type='l', col=0, xlab='Years Post Abandonment', ylab='mean growth')
for(i in 1:length(unique(df$SITE2))){
	tmp <- df[df$SITE2 == unique(df$SITE2)[i],]
	lty <- ifelse(substring(as.character(unique(tmp$SITE2)),1,2)=='PE', 1,1)
	points(tmp$age, tmp$V1, type='l', col=col[i], lwd=3, lty=lty)
}
points(tapply(gdata$growth, gdata$age, mean), type='l', lwd=20, col=grey(.5,.75))
# add points to show relative sample size
# cex <- Rescale01 (tapply(gdata$growth, gdata$age, length))*2+.5
# points(tapply(gdata$growth, gdata$age, mean), cex=cex, pch=21, bg=1)




###########################
### HISTOGRAMS OF STEM DIAMETER
hist(sdata$DBH, prob=T)
sdata <- droplevels(sdata)
ages <- tapply(sdata$AGE00, sdata$SITE2, mean)

par(mfrow=c(4,4))
for(i in 1:length(unique(sdata$SITE2))){
		tmp <- sdata[sdata$SITE2 == unique(sdata$SITE2)[i],]
		if(i == 17) {quartz(); par(mfrow=c(4,4))}
		hist(tmp$DBH, add=F, col=col[i], prob=T, main=paste(unique(sdata$SITE2)[i], ages[i]), xlim=c(0,40))
}





###########################
### PLOT REL ABUND OF INDIVIDUAL SPECIES WITH FOREST AGE
par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(2,2,1,1))
for(s in 1:9){
focsp <- names(rev(sort(table(sdata$SPECIES))))[s]
tmp <- sdata[sdata$SPECIES == focsp,]
tmp <- table(tmp$SITE2, tmp$age)
b <- sdata[sdata$SPECIES != focsp,]
b <- table(b$SITE2, b$age)
a <- b*0
a[rownames(bb) %in% rownames(tmp),colnames(bb) %in% colnames(tmp)] <- tmp
a <- round(a/(a+b), 3)
plot(c(1,30),c(0,1), col=NA, xlab='', ylab='', cex.lab=1.25, ylim=c(0,max(a, na.rm=T)+.1))
mtext(focsp, 3, -1.25, cex=.75, font=2)
if(s %in% c(1,4,7)) mtext('Relative abundance', 2, 2.25)
if(s %in% c(7,8,9)) mtext('Plot age', 1, 2.5)
for(i in 1:nrow(a)){
	aa <- a[i,]
	points(aa, type='l', col=col[i], lwd=2, lty=lty)
	}
}


###########################
### PLOT SPECIES RICHNESS BY FOREST AGE
library(vegan)
tmp <- table(gdata$age, gdata$SITE2 , gdata$SPECIES)
dim(tmp)
indv <- list()
for (i in 1:dim(tmp)[1]){ indv[[i]] <- tmp[i,,]}
indv <- lapply(indv, rowSums)
sr <- list()
for (i in 1:dim(tmp)[1]){ sr[[i]] <- tmp[i,,]>0}
sr <- lapply(sr, rowSums)

r <- list()
for (i in 1:length(sr)){ r[[i]] <- sr[[i]]/indv[[i]]}

sr <- matrix(nrow=22, ncol=29, dimnames=list(rownames(tmp[1,,]),1:29))
for (i in 1:length(r)){ sr[,i] <- r[[i]] }
sr[is.nan(sr)] <- NA
meansr <- apply(sr, 2, mean, na.rm=T)

plot(0,0, ylim=c(0,max(sr, na.rm=T)), xlim=c(0,30), type='l', col=0, xlab='Years Post Abandonment', ylab='Species Per Stem')
for(i in 1:nrow(sr)){
	lty <- ifelse(substring(as.character(rownames(sr)[i]),1,2)=='PE', 1,1)
	points(1:29, sr[i,], type='l', col=col[i], lwd=3, lty=lty)
}













