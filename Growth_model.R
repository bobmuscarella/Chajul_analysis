#####################################################################
### HIERARCICHAL MODEL OF TRAIT-MEDIATED GROWTH DURING SUCCESSION ###
#####################################################################
library(jagsUI)
library(FD)

###################
#### LOAD DATA ####
###################
setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
load("Chajul_data_processed_wtraits_11.20.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")
#load("Chajul_data_processed_wtraits_4.27.15.RDA")
#load("Chajul_census_processed_8.25.15.RDA")
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
tdata <- read.csv('Mean_traits_4.27.15.csv')

head(data)

age <- census$ages + census$CENSUS
data$age <- as.vector(age[match(data$SITE.CENSUS, census$SITE.CENSUS)])
quadBA <- tapply(data$ba, paste(data$SITE.CENSUS, data$QUAD), sum, na.rm=T)
data$quadBA <- quadBA[match(paste(data$SITE.CENSUS, data$QUAD), names(quadBA))]


#plot(log(data$quadBA[data$quadBA!=0]+1), jitter(data$age[data$quadBA!=0]))

cor(log(data$quadBA+1), data$growth, use='p')
cor(data$age, data$growth, use='p')



### FUNCTION TO CONVERT DBH TO BA
dbh2ba <- function(dbh){
  return(((pi / 40000) * (dbh^2)))
}

data <- data[!is.na(data$DBH),]
data <- data[!is.na(data$growth),]
data <- droplevels(data)


###################
#### PREP DATA ####
###################

### CHECK FOR GROWTH OUTLIERS...
outliers <- vector()
for(i in -30:30){
  outliers[i] <- length(data$uID[!is.na(data$growth) & data$growth > i])
}
#data <- data[!is.na(data$growth) & data$growth < 8,]

### BCI METHOD
Growth.Include <- (((data$growth*10) > 4 * (-(0.0062 * data$DBH + 0.904))) & ((data$growth*10) < 75))
data <- data[Growth.Include==T,]
data <- droplevels(data)

### REMOVE STEMS BASED ON SD CUTOFF
#glim <- tapply(data$growth, data$SPECIES, sd, na.rm=T) * 5
#glim <- glim[!is.na(glim)]
#data$glim <- glim[match(data$species,names(glim))]
#sum(abs(data$growth) > data$glim, na.rm=T)
#data <- data[abs(data$growth) < data$glim,]

####################################
#### SUBSET DATA FOR TESTING... ####
#data <- data[sample(1:nrow(data), 2000),]

####################################

### GET PLOT AGES
age <- census$ages + census$CENSUS

### SELECT A TRAIT
focal.trait <- 'log.SLA'
data <- data[!is.na(data[,focal.trait]),]
data <- droplevels(data)

data <- data[order(data$uID, data$SPECIES, data$SITE.CENSUS),]

### MAKE INPUT DATA
d <- list(
    ntree = nrow(data),
    nindiv = length(unique(data$uID)),
    nspecies = length(unique(data$SPECIES)),
    nplot = length(unique(data$SITE2)),
    trait = as.vector(scale(tapply(data[,focal.trait], data$SPECIES, mean))),
    species = as.numeric(data$SPECIES),
    indiv = as.numeric(as.factor(data$uID)),
    growth = as.vector(scale(data$growth, center=F)),
    dbh = as.vector(scale(log(data$DBH))),
#    age = as.vector(scale(age[match(data$SITE.CENSUS, census$SITE.CENSUS)])),
    quadba = as.vector(scale(log(data$quadBA))),
    plot = as.numeric(data$SITE2)
  )


##############################
#### Write the model file ####
##############################
sink("Growth_Age_x_Trait.bug")

cat(" model {
    
    for( i in 1:ntree ) {

      growth[i] ~ dnorm(mu[i], tau[1])

      mu[i] <- exp(z[i])    

      z[i] <- beta.1[species[i]]
              + beta.2[species[i]] * (quadba[i])
              + beta.3[species[i]] * (dbh[i])
              + indiv.effect[indiv[i]]
              + plot.effect[plot[i]]
    }
    
    for( j in 1:nspecies ) {
      beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])  # Avg Growth
      beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])  # Trait Effect
      beta.3[j] ~ dnorm(mu.beta[3], tau[4])                         # DBH Effect
    }

    for( i.a in 1:nindiv ) {
      indiv.effect[i.a] ~ dnorm(0, tau[5])
    }
    
    for( p.a in 1:nplot ) {
      plot.effect[p.a] ~ dnorm(0, tau[6])
    }

    beta.t[1] ~ dnorm(0, 1E-4)
    beta.t[2] ~ dnorm(0, 1E-4)
    
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



################################################
### Set initial values, monitors, iterations and run model ###
################################################

#### TRY WITH A MIXED MODEL FIRST...
library(lme4)
library(MuMIn)
library(arm)

growth <- d$growth
trait <- d$trait[d$species]
species <- d$species
age <- d$age
dbh <- d$dbh
plot <- d$plot
cen <- d$census
quadba <- d$quadba


m1 <- lm(growth ~ quadba + trait + quadba*trait + dbh)
m2 <- lmer(growth ~ quadba + trait + quadba*trait + dbh + (1|plot))
m3 <- lmer(growth ~ age + trait + age*trait + dbh + (1|cen))

summary(m1)
summary(m2)
summary(m3)

mod <- m2

summary(mod)

est <- apply(fixef(sim(m2, n.sims=1000)), 2, quantile, probs=c(0.025, 0.5, 0.975))

r.squaredGLMM(mod)

pch <- ifelse(sign(est[1,])==sign(est[3,]), 16, 1)
par(mar=c(4,10,4,6))
plot(est[2,], ncol(est):1, xlim=c(min(est), max(est)), pch=pch, axes=F, cex=2, ylab='', xlab='Standard Effect')
segments(est[1,], ncol(est):1, est[3,], lwd=3, col=2)
abline(v=0, lty=2)
axis(2, labels=colnames(est), at=ncol(est):1, las=2)
axis(1)
box()



############

# Set initial values
inits <- function (){
  list(
    beta.t = rnorm(2),
    mu.beta = rnorm(3),
    tau = rgamma(6, 1E3, 1E3)
  )
}

# Set monitors
#params <- c("beta.1","beta.2","beta.3","beta.t","mu.beta","sigma")
params <- c("beta.t","mu.beta","sigma")

# Run model
adapt <- 500
iter <- 1000
burn <- 500
thin <- 2
chains <- 3

mod <- jagsUI::jags(d, inits, params, 
                    "Growth_Age_x_Trait.bug", 
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
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()


















