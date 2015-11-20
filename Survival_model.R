#######################################################################
### HIERARCICHAL MODEL OF TRAIT-MEDIATED SURVIVAL DURING SUCCESSION ###
#######################################################################
library(jagsUI)
library(FD)

###################
#### LOAD DATA ####
###################
setwd("/Users/Bob/Projects/Postdoc/Demo Drivers of FD/DATA")
load("Chajul_data_processed_wtraits_4.27.15.RDA")
load("Chajul_census_processed_8.25.15.RDA")
totalBA <- tapply(data$ba, data$SITE.CENSUS, sum, na.rm=T)
data$totalBA <- totalBA[match(data$SITE.CENSUS, names(totalBA))]
data$relBA <- data$ba/data$totalBA
census$totalBA <- totalBA[match(census$SITE.CENSUS, names(totalBA))]
tdata <- read.csv('Mean_traits_4.27.15.csv')

# FIX AN INCORRECT DATE
data$int[data$SITE.CENSUS %in% c('PEHEC17.12')] <- 22+365


###################
#### PREP DATA ####
###################
data <- data[!is.na(data$survive),]

##############################
### SUBSET DATA FOR TESTING...
focsp <- sort(unique(data$species))[as.logical(rowSums((table(data$species, data$survive)==0) > 0))]
data <- data[data$species %in% focsp,]

data <- data[sample(1:nrow(data), 500),]
##############################

### GET PLOT AGES
age <- census$ages + census$CENSUS

### SELECT A TRAIT
focal.trait <- 'WD'
data <- data[!is.na(data[,focal.trait]),]
data <- droplevels(data)

data <- data[order(as.numeric(as.factor(data$uID)), data$species, data$SITE),]

### MAKE INPUT DATA
d <- list(
  ntree = nrow(data),
  nindiv = length(unique(data$uID)),
  nspecies = length(unique(data$SPECIES)),
  nplot = length(unique(data$SITE2)),
  trait = as.vector(scale(tapply(data[,focal.trait], data$SPECIES, mean))),
  species = as.numeric(data$SPECIES),
  indiv = as.numeric(as.factor(data$uID)),
  days = as.vector(data$int),
  survive = data$survive,
  dbh = as.vector(scale(log(data$DBH))),
  age = as.vector(scale(age[match(data$SITE.CENSUS, census$SITE.CENSUS)])),
  plot = as.numeric(data$SITE2)
)



##############################
#### Write the model file ####
##############################
sink("Survive_Age_x_Trait.bug")

cat(" model {
    
    for( i in 1:ntree ) {
    
    survive[i] ~ dinterval(t[i], days[i])

    t[i] ~ dweib(r, mu[i])

    mu[i] <- exp(z[i])

    z[i] <- beta.1[species[i]]
    + beta.2[species[i]] * (age[i])
    + beta.3[species[i]] * (dbh[i])
    + indiv.effect[indiv[i]]
    + plot.effect[plot[i]]

    }

    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[1])
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[2])
    beta.3[j] ~ dnorm(mu.beta[3], tau[3])
    }
    
    for( i.a in 1:nindiv ) {
    indiv.effect[i.a] ~ dnorm(0, tau[4])
    }
    
    for( p.a in 1:nplot ) {
    plot.effect[p.a] ~ dnorm(0, tau[5])
    }

    beta.t[1] ~ dnorm(0, 1E-4)
    beta.t[2] ~ dnorm(0, 1E-4)
    
    for( m in 1:3 ) {
    mu.beta[m] ~ dnorm(0, 1E-4)
    }
    
    for( t in 1:5 ) {
    tau[t] ~ dgamma(1E3, 1E3)
    }

    r ~ dgamma(1, 1E-3)
    sigma <- 1 / sqrt(tau)
    
    }"
    , fill=TRUE)
sink()



################################################
### Set initial values, monitors, iterations and run model ###
################################################

inits <- function (eps=0.1){  
  list(beta.t = rnorm(2),
       mu.beta = rnorm(3),
       tau = rgamma(5, 1E3, 1E3),
       r = 2,
       t = with(d, days + ifelse(survive, eps, -eps)))
}

# Set monitors
#params <- c("beta.1","beta.2","beta.3","beta.t","mu.beta","sigma")
params <- c("beta.t","mu.beta","sigma")


# Run model
adapt <- 250#2500
iter <- 500#5000
burn <- 250#2500
thin <- 1#5
chains <- 3

mod <- jagsUI::jags(d, inits, params, 
                    "Survive_Age_x_Trait.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=F)

table(data$species, data$survive)


mod
mod.sim <- update(mod.sim, n.iter=250)
plot(mod)

#save(mod, file='wd_Survive_Age_x_Trait_11.9.15.RDA')


x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma',rownames(x)),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()




x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma',rownames(x)),]
x <- x[(nrow(x)-4):(nrow(x)),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()




par(mar=c(15,15,1,15))
x <- x[c('beta.t1','beta.t2','mu.beta2','mu.beta3'),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5, xlim=c(0,5))
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
abline(h=0,lty=2)
axis(1, labels=c("Trait x mean effect",
                 "Trait x age effect", 
                 "Age effect","DBH effect"), at=1:nrow(x), las=2)
axis(2); box()





x <- x[paste('beta.1',1:d$nspecies,sep=''),]

plot(d$trait, x[,2])









