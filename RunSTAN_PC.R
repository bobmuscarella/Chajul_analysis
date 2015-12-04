setwd("K:/Bob/DemoDriversofFD/GIT/")

n = nrow(data)
nindiv = length(unique(data$uID))
nspecies = length(unique(data$SPECIES))
nplot = length(unique(data$SITE2))
trait = trait.z
species = as.numeric(data$SPECIES)
indiv = as.numeric(as.factor(data$uID))
growth = as.vector(data$growth.z)
dbh = data$log.dbh.z
age = data$age.z
plotba = data$plotBA.z
plotagb = data$plotAGB.z
plot = as.numeric(data$SITE)
census = as.numeric(as.factor(data$SITE.CENSUS))



dataSTAN <- list(growth = growth, 
                 species = species, 
                 nspecies=max(as.numeric(as.factor(as.character(species)))), 
                 plot = plot, 
                 nplot = max(as.numeric(as.factor(as.character(plot)))), 
                 indiv = indiv, 
                 nindiv = max(as.numeric(as.factor(as.character(indiv)))),
                 n = nrow(data), 
                 trait = trait, 
                 dbh = dbh,
                 plotagb = plotagb)

library(rstan)

fit <- stan(file="Bob.stan", 
            data=dataSTAN, 
            iter=1000, 
            chains=4, cores=4,
            thin=2
            )

# Check convergence looking at the last column
print(fit,c("tau","beta.mu","beta.t"))

#posterior distribution
tau1<-extract(fit,"tau[1]")
tau2<-extract(fit,"tau[2]")
tau3<-extract(fit,"tau[3]")
tau4<-extract(fit,"tau[4]")
tau5<-extract(fit,"tau[5]")
beta.mu1<-extract(fit,"beta.mu[1]")
beta.mu2<-extract(fit,"beta.mu[2]")
beta.mu3<-extract(fit,"beta.mu[3]")
beta.t1<-extract(fit,"beta.t[1]")
beta.t2<-extract(fit,"beta.t[2]")

res_q<-apply(cbind(beta.mu1[[1]],beta.mu2[[1]],beta.mu3[[1]],beta.t1[[1]],beta.t2[[1]],tau1[[1]],tau2[[1]],tau3[[1]],tau4[[1]],tau5[[1]]),2,quantile,c(0.025,0.05,0.5,0.95,0.975))
colnames(res_q)<-c("beta.mu3","beta.mu2","beta.mu3","beta.t1","beta.t2","tau1","tau2","tau3","tau3","tau5")

res_q

#extract the predictions
ypred<-extract(fit,"ynew")
ypred_q<-apply(ypred[[1]],2,quantile,c(0.025,0.05,0.5,0.95,0.975))

