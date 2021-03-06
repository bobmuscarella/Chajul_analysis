#####################################################################
### HIERARCICHAL MODEL OF TRAIT-MEDIATED GROWTH DURING SUCCESSION ###
#####################################################################
library(jagsUI)
library(lme4)
library(MuMIn)
library(arm)

###################
#### LOAD DATA ####
###################
setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
load("Chajul_data_processed_wtraits_12.3.15.RDA")
load("Chajul_census_processed_12.3.15.RDA")
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
tdata <- read.csv('Mean_traits_4.27.15.csv')
data$log.SV <- log(data$SV)

### FUNCTION TO CONVERT DBH TO BA
dbh2ba <- function(dbh){
  return(((pi / 40000) * (dbh^2)))
}

### FUNCTION TO Z-TRANSFORM DATA
z.score <- function (data, center=T) {
  xm <- ifelse(center==T, mean (data, na.rm=TRUE), 0)
  xsd <- sd(data, na.rm=TRUE)
  xtrans <- (data - xm) / (2 * xsd) 
  return(xtrans)
}

###################
#### PREP DATA ####
###################
age <- census$ages + census$CENSUS
data$age <- as.vector(age[match(data$SITE.CENSUS, census$SITE.CENSUS)])
plotBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$plotBA <- plotBA[match(data$SITE.CENSUS, names(plotBA))]
plotAGB <- tapply(data$agb, data$SITE.CENSUS, sum, na.rm=T)
data$plotAGB <- plotAGB[match(data$SITE.CENSUS, names(plotAGB))]

data <- data[!is.na(data$DBH),]
data <- data[!is.na(data$growth),]
data <- droplevels(data)

### CHECK FOR GROWTH OUTLIERS...
#outliers <- vector()
#for(i in -30:30){
#  outliers[i] <- length(data$uID[!is.na(data$annual.growth) & data$annual.growth > i])
#}

### CONDIT METHOD TO EXCLUDE OUTLIERS
Growth.Include <- (data$annual.growth < 7.5 & data$annual.growth > (data$DBH * -0.25))

### BCI METHOD TO EXCLUDE STEMS...
# Growth.Include <- (((data$growth*10) > 4 * (-(0.0062 * data$DBH + 0.904))) & ((data$growth*10) < 75))

### REMOVE STEMS BASED ON SD CUTOFF
# glim <- tapply(data$growth, data$SPECIES, sd, na.rm=T) * 5
# glim <- glim[!is.na(glim)]
# data$glim <- glim[match(data$species,names(glim))]
# sum(abs(data$growth) > data$glim, na.rm=T)
# data <- data[abs(data$growth) < data$glim,]


data <- data[Growth.Include==T,]
data <- droplevels(data)


####################################
#### SUBSET DATA FOR TESTING... ####
# data <- data[sample(1:nrow(data), 10000),]
####################################

### SELECT A TRAIT
focal.trait <- 'log.SV'
data <- data[!is.na(data[,focal.trait]),]
data <- droplevels(data)
data <- data[order(data$SPECIES, data$uID),]

### SCALE (+ center) DATA, OVERALL DATA
data$growth.z <- as.vector(scale(data$log.growth, center=F))
data$age.z <- z.score(data$age)
data$plotAGB.z <- z.score(data$plotAGB)
data$plotBA.z <- z.score(data$plotBA)
data$log.dbh.z <- z.score(log(data$DBH))

### SCALE AND CENTER DATA, WITHIN SPECIES
# data$log.dbh.z <- unlist(tapply(log(data$DBH), data$species, z.score))


### Center trait values (and convert SLA to LMA)
if(focal.trait=='log.sFtP'){  
trait.z <- as.vector(z.score(tapply(1/data[,focal.trait], data$SPECIES, mean))) 
} else {
  trait.z <- as.vector(z.score(tapply(data[,focal.trait], data$SPECIES, mean))) 
}

### MAKE INPUT DATA
d <- list(
    ntree = nrow(data),
    nindiv = length(unique(data$uID)),
    nspecies = length(unique(data$SPECIES)),
    nplot = length(unique(data$SITE2)),
    trait = trait.z,
    species = as.numeric(data$SPECIES),
    indiv = as.numeric(as.factor(data$uID)),
    growth = data$growth.z,
    dbh = data$log.dbh.z,
    age = data$age.z,
    plotba = data$plotBA.z,
    plotagb = data$plotAGB.z,
    plot = as.numeric(data$SITE),
    census = as.numeric(as.factor(data$SITE.CENSUS))
)


################################################
##################### LMER #####################
################################################
growth <- data$growth.z
indiv <- d$indiv
trait <- trait.z[d$species]
species <- d$species
dbh <- d$dbh
plot <- d$plot
plotba <- d$plotba

m1 <- lmer(growth ~ plotba + trait + plotba * trait + dbh 
           + (1|plot) + (1|species) + (1|indiv))

mod <- m1

est <- apply(fixef(sim(mod, n.sims=500)), 2, quantile, probs=c(0.025, 0.5, 0.975))

#est <- est[,-1]
pch <- ifelse(sign(est[1,])==sign(est[3,]), 16, 1)
#par(mar=c(4,10,4,6))
#pdf(file='lma.coeff.pdf')
par(mar=c(20,15,4,5))
plot(est[2,], ncol(est):1, xlim=c(ifelse(min(est)>0,-.1,min(est)), ifelse(max(est)<0,.1,max(est))), pch=pch, axes=F, ylab='', xlab='Standard Effect')
segments(est[1,], ncol(est):1, est[3,], lwd=3)
abline(v=0, lty=2)
axis(2, labels=colnames(est), at=ncol(est):1, las=2)
axis(1)
box()
#dev.off()


# PLOT INTERACTION PREDICTION
newx <- seq(min(plotba), max(plotba), length.out=50)
newdata.hightrait <- data.frame(plotba=newx, 
                                trait=rep(max(trait), length(newx)), 
                                dbh=rep(0,length(newx)),                                 
                                plot=factor(rep(1,length(newx))),
                                cen=factor(rep(1,length(newx))),
                                species=factor(rep(1,length(newx))),
                                indiv=factor(rep(1,length(newx))))
newdata.lowtrait <- data.frame(plotba=newx, 
                               trait=rep(min(trait), length(newx)), 
                               dbh=rep(0,length(newx)), 
                               plot=factor(rep(0,length(newx))),
                               cen=factor(rep(1,length(newx))),
                               species=factor(rep(1,length(newx))),
                               indiv=factor(rep(1,length(newx))))
y.pred.hightrait <- predict(mod, newdata.hightrait, allow.new.levels=T)
y.pred.lowtrait <- predict(mod, newdata.lowtrait, allow.new.levels=T)

#pdf(file='lma_interaction.pdf')
par(mar=c(6,6,4,4))
ylim <- range(y.pred.lowtrait, y.pred.hightrait)
plot(newx, y.pred.hightrait, type='l', lwd=4, 
     ylab='', xlab='Plot BA', axes=F, cex.lab=1.5, ylim=ylim)
labs <- seq(0, 5, by=0.2)
axis(1, labels=labs, at=(labs - mean(data$plotBA)) / sd(data$plotBA), cex.axis=1.25)
labs2 <- seq(-2,2, by=0.2)
rms <- sqrt(sum(data$annual.growth^2)/(length(data$annual.growth)-1)) # TO BACK TRANSFORM GROWTH PREDICTIONS TO CM / YEAR
axis(2, labels=labs2, at=(labs2 / rms), las=2, cex.axis=1.25)
points(newx, y.pred.lowtrait, type='l', col=2, lwd=4)
legend('topright',legend=c('High trait', 'Low trait'), lty=1, col=1:2, bty='n', lwd=4, cex=1.25)
abline(h=0,lty=3)
mtext(bquote('Predicted Growth Rate ( mm/yr'^1~')'), 2, 3.5, cex=1.5)
# dev.off()




#/////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\
######################################################
##################### JAGS MODEL #####################
######################################################

sink("Growth_Age_x_Trait.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- exp(z[i])    
    
    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * (plotba[i])
    + beta.3[species[i]] * (dbh[i])
    + indiv.effect[indiv[i]]
    + plot.effect[plot[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])  # Trait-mediated avg growth
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])  # Trait x age response
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])  # DBH Effect
    }
    
    for( i.a in 1:nindiv ) {
     indiv.effect[i.a] ~ dnorm(0, tau[5])
    }
    
    for( p.a in 1:nplot ) {
      plot.effect[p.a] ~ dnorm(0, tau[6])
    }
    
    for( b in 1:2 ) {
      beta.t[b] ~ dnorm(0, 1E-4)
    }

    for( m in 1:3 ) {
      mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:6 ) {
    tau[t] ~ dgamma(1E3, 1E3)
    }
    
    sigma <- 1 / sqrt(tau)
    }",
fill=TRUE)
sink()



sink("OneLevel_Growth_Age_x_Trait.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    growth[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- exp(z[i])    
    
    z[i] <- mu.beta[1]
    + mu.beta[2] * (plotba[i])
    + mu.beta[3] * (trait[species[i]])
    + mu.beta[4] * (trait[species[i]]) * (plotba[i])
    + beta.1[species[i]] * (dbh[i])
    + indiv.effect[indiv[i]]
    + plot.effect[plot[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[5], tau[2])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[3])
    }
    
    for( p.a in 1:nplot ) {
    plot.effect[p.a] ~ dnorm(0, tau[4])
    }
    
    for( m in 1:5 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:4 ) {
    tau[t] ~ dgamma(1E3, 1E3)
    }
    
    sigma <- 1 / sqrt(tau)
    }",
fill=TRUE)
sink()



################################################
### Set initial values, monitors, iterations and run model ###
################################################

inits <- function (){
  list(
    beta.t = rnorm(2),
    mu.beta = rnorm(3),
    tau = rgamma(6, 1E3, 1E3)
  )
}



inits <- function (){
  list(
    mu.beta = rnorm(5),
    tau = rgamma(4, 1E3, 1E3)
  )
}



# Set monitors
#params <- c("beta.1","beta.2","beta.3","beta.t","mu.beta","sigma")
#params <- c("beta.t","mu.beta","sigma")
params <- c("beta.1","mu.beta","sigma")

# Run model
adapt <- 500
iter <- 500
burn <- 250
thin <- 2
chains <- 2

mod2 <- jagsUI::jags(d, inits, params, 
                    "OneLevel_Growth_Age_x_Trait.bug",
#                    "Growth_Age_x_Trait.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=F)

mod

plot(mod)

mod <- update(mod, n.iter=1000)


x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma',rownames(x)),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3], lwd=2)
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()




m <- mod2
x <- cbind(unlist(m$q50),unlist(m$q2.5),unlist(m$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma',rownames(x)),]
x <- x[grep('mu',rownames(x)),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3], lwd=2)
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()


















