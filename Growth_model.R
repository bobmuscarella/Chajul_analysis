
#####################################################################
### HIERARCICHAL MODEL OF TRAIT-MEDIATED GROWTH DURING SUCCESSION ###
#####################################################################
library(jagsUI)
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

### FUNCTION TO CONVERT DBH TO BA
dbh2ba <- function(dbh){
  return(((pi / 40000) * (dbh^2)))
}

###################
#### PREP DATA ####
###################

### CHECK FOR GROWTH OUTLIERS...
outliers <- vector()
for(i in -30:30){
  outliers[i] <- length(data$uID[!is.na(data$growth) & data$growth > i])
}
data <- data[!is.na(data$growth) & data$growth < 8,]


### SUBSET DATA FOR TESTING...
# data <- data[sample(1:nrow(data), 1000),]

### GET PLOT AGES
age <- census$ages + census$CENSUS

### SELECT A TRAIT
focal.trait <- 'WD'
data <- data[!is.na(data[,focal.trait]),]
data <- droplevels(data)

### MAKE INPUT DATA
d <- list(
  ntree = nrow(data),
  nindiv = length(unique(data$uID)),
  nspecies = length(unique(data$SPECIES)),
  nplot = length(unique(data$SITE2)),
  trait = as.vector(scale(tapply(data[,focal.trait], data$SPECIES, mean))),
  species = as.numeric(data$SPECIES),
  indiv = as.numeric(as.factor(data$uID)),
  growth = data$growth,
  dbh = as.vector(scale(log(data$DBH))),
  age = as.vector(scale(age[match(data$SITE.CENSUS, census$SITE.CENSUS)])),
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
    + beta.2[species[i]] * (age[i])
    + beta.3[species[i]] * (dbh[i])
    + indiv.effect[indiv[i]]
    + plot.effect[plot[i]]
    }
    
    for( j in 1:nspecies ) {
    beta.1[j] ~ dnorm(mu.beta[1] + beta.t[1] * trait[j], tau[2])
    beta.2[j] ~ dnorm(mu.beta[2] + beta.t[2] * trait[j], tau[3])
    beta.3[j] ~ dnorm(mu.beta[3], tau[4])
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
    tau[t] ~ dgamma(1E-3, 1E-3)
    }
    
    sigma <- 1 / sqrt(tau)
    
    }"
    , fill=TRUE)
sink()



################################################
### Set initial values, monitors, iterations and run model ###
################################################

# Set initial values
inits <- function (){
  list(
    beta.t = rnorm(2),
    mu.beta = rnorm(3),
    tau = rgamma(6, 1E-3, 1E-3) + 1E-5)
}

# Set monitors
params <- c("beta.t","mu.beta","sigma")

# Run model
adapt <- 100
iter <- 250
burn <- 150
thin <- 2
chains <- 3

mod <- jagsUI::jags(d, inits, params, 
                    "Growth_Age_x_Trait.bug", 
                    n.chains=chains, n.adapt=adapt, n.iter=iter, 
                    n.burnin=burn, n.thin=thin, parallel=F)

plot(mod)

x <- cbind(unlist(mod$q50),unlist(mod$q2.5),unlist(mod$q97.5))
x <- x[-nrow(x),]
x <- x[-grep('sigma',rownames(x)),]
pch <- ifelse(sign(x[,2])==sign(x[,3]), 16, 1)
plot(x[,1],ylim=c(min(x), max(x)), pch=pch, axes=F, ylab='Std. Effect', xlab='', cex=1.5)
segments(1:nrow(x), x[,2], 1:nrow(x), x[,3])
abline(h=0,lty=2)
axis(1, labels=rownames(x), at=1:nrow(x), las=2)
axis(2); box()


















