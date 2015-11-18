######################################
######################################
######################################
dbh2ba <- function(dbh){
  return((pi / 40000) * (dbh^2))
}



# START SIMULATION BY SAYING NUMBER OF SPECIES IN POOL
#nsp <- 50
nsp <- 10

# GENERATE A VECTOR OF TRAIT VALUES
trait <- abs(rnorm(nsp, 0, 5))
trait2 <- abs(rnorm(nsp, 0, 1))
#trait <- c(1:nsp)

# TWO SPECIES INITIAL BASAL AREA
ba0 <- runif(nsp, 2, 5) * (1/trait)
#ba0 <- c(20,20)

# INITIAL CWM
cwm0 <- sum(ba0/sum(ba0) * trait)

# AMOUNT OF BA GAINED BY GROWTH
ba.growth <- abs(rnorm(nsp)) * trait
#ba.growth <- c(1,1)

# AVERAGE TRAIT VALUE OF GROWTH BA
cwm.g <- sum(ba.growth * trait) / sum(ba.growth)
#cwm.g <- sum(ba.growth * trait)

# AMOUNT OF BA LOST TO MORTALITY
mort <- rlnorm(nsp, .1)
#ba.mort <- c(1,1)
ba.mort <- ifelse(mort < ba0, mort, 0)

# AVERAGE TRAIT VALUE OF DEAD BA
cwm.m <- sum(ba.mort * trait) / sum(ba.mort)
#cwm.m <- sum(ba.mort * trait)

# NEW BASAL AREA (OLD BA + GROWTH - MORTALITY)
ba1 <- ba0 + ba.growth - ba.mort

# NEW CWM
cwm1 <- sum(ba1/sum(ba1) * trait)

# WHAT IS THE PROPORTION CONTRIBUTION OF GROWTH AND MORTALITY TO THE OBSERVED CHANGE IN CWM?

# CHANGE IN CWM
# delta.cwm <- cwm1 - cwm0
# delta.ba <- sum(ba1 - ba0)

# cwm.flux <- (cwm.g * sum(ba.growth)) + (cwm.m * sum(ba.mort))
# cwm.flux <- (cwm0*sum(ba0)) + (cwm.g * sum(ba.growth)) + (cwm.m * sum(ba.mort))
# ba.flux <- sum(ba.growth) + sum(ba.mort)

# t0.prop <- (cwm0 * sum(ba0)) / cwm.flux
# g.prop <- (cwm.g * sum(ba.growth)) / cwm.flux
# m.prop <- (cwm.m * sum(ba.mort)) / cwm.flux


# plot(c(sum(ba0), sum(ba1)), c(cwm0,cwm1), ylim=range(c(cwm.g, cwm.m, cwm0, cwm1)), pch=21, xlim=c(0, max(sum(ba0), sum(ba1))), cex=3, bg=c(0,'grey'))
# arrows(sum(ba0), cwm0, sum(ba1), cwm1, len=.1, lwd=3)
# abline(h=cwm0, lty=3)
# points(sum(ba.growth), cwm.g, pch=21, bg=3, cex=3)
# points(sum(ba.mort), cwm.m, pch=21, bg=2, cex=3)
# mtext(paste("t0", round(t0.prop, 3), "| Growth: ", round(g.prop, 3), "| Mortality: ", round(m.prop, 3)))
#mtext(paste("Growth: ", round(g.prop, 3), "| Mortality: ", round(m.prop, 3)))

# plot(c(sum(ba0), sum(ba1)), c(cwm0,cwm1), ylim=range(trait), pch=21, xlim=c(min(ba0, ba.growth, ba.mort),max(sum(ba0+ba.growth), sum(ba0+ba.mort))), cex=c(3,2), bg=c(0,'grey'))
# points(ba.growth+sum(ba.growth), trait, pch=22, bg=3, cex=2)
# points(ba.mort+sum(ba.mort), trait, pch=ifelse(ba.mort>0, 22, NA), bg=2)
# points(ba0+sum(ba0), trait, pch=22, bg=0, cex=2)
# points(ba1+sum(ba1), trait, pch=22, bg='grey', cex=1)
# arrows(sum(ba0), cwm0, sum(ba1), cwm1, len=.1, lwd=3)
# abline(h=cwm0, lty=3)
# points(sum(ba.growth), cwm.g, pch=21, bg=3, cex=2)
# points(sum(ba.mort), cwm.m, pch=21, bg=2, cex=1)
# mtext(paste("t0", round(t0.prop, 3), "| Growth: ", round(g.prop, 3), "| Mortality: ", round(m.prop, 3)))

# cwm1 == ((cwm0*sum(ba0)) + (cwm.g * sum(ba.growth)) - (cwm.m * sum(ba.mort))) / (sum(ba1))

# How to attribute proportion of CWM change due to growth vs. mortality?
# How does a little bit of growth but imposing a big trait value compare to a lot of growth and trait value near the current CWM?
# What does that mean biologically?
# How can we tease apart absolute magnitude of growth / mortality versus relative contribution to CWM shifts?


#########################################
#### Partition change in functional composition to components of growth, recruitment, mortality, colonization and extinction.
#########################################
# GROWTH: Increase in basal area for all individuals between t0 and t1 that were present in t0.
# RECRUITMENT: Increase in basal area from new individuals of species present in t0 and t1.
# MORTALITY: Change in basal area for all individuals that died between t0 - t1 for species present in both times.
# COLONIZATION: Addition of basal area from new recruits of species absent in t0 but present in t1.
# EXTINCTION: Reduction in basal area for all species present in t0 but absent in t1.

ba <- as.data.frame(cbind(ba0, ba.growth, ba.mort, ba1))
rownames(ba) <- letters[1:nsp]
names(trait) <- letters[1:nsp]
names(trait2) <- letters[1:nsp]
trt <- cbind(trait, trait2)
trt <- apply(trt, 2, scale)
rownames(trt) <- letters[1:nsp]

ba$ba0g <- ba$ba0 + ba$ba.growth
ba$ba0m <- ba$ba0 - ba$ba.mort
ba$ba0gm <- ba$ba0 + ba$ba.growth - ba$ba.mort

res <- dbFD(trt, t(ba))



par(mfrow=c(3,3), mar=c(2,2,2,0))
plot(trt[,1], trt[,2], cex=ba[,'ba0'], col=NA, main="Time 1")
cwm0 <- colSums(sweep(trt, 1, ba0, "*"))/sum(ba0)
#for(i in 1:nsp){
  lty <- ifelse(ba[,'ba0']>0, 1, 0)
  segments(cwm0[1], cwm0[2], trt[i,1], trt[i,2], col=1, lty=lty)
}
points(cwm0[1], cwm0[2], pch='+', cex=1.5)
points(trt[,1], trt[,2], pch=21, cex=ba[,'ba0'], bg=rgb(0,0,0,.4))
draw.circle(cwm0[1], cwm0[2], radius=res$FDiv[1], lwd=2, border=5, lty=2)
draw.circle(cwm0[1], cwm0[2], radius=res$FDis[1], lwd=2, border=6, lty=2)
Plot_ConvexHull(trt[,1], trt[,2], lty=3, lwd=2)

plot(trt[,1], trt[,2], cex=ba[,'ba.growth'], pch=21, bg=rgb(0,.6,0,.9), main="Growth")
cwm.g <- colSums(sweep(trt, 1, ba.growth, "*"))/sum(ba.growth)
#for(i in 1:nsp){
#  lty <- ifelse(ba[,'ba.growth']>0, 1, 0)
#  segments(cwm.g[1], cwm.g[2], trt[i,1], trt[i,2], col=3, lty=lty[i])
#}
points(cwm.g[1], cwm.g[2], pch='+', cex=1.5, col=rgb(0,.3,0,.9))
draw.circle(cwm.g[1], cwm.g[2], radius=res$FDiv[2], lwd=2, border=5, lty=2)
draw.circle(cwm.g[1], cwm.g[2], radius=res$FDis[2], lwd=2, border=6, lty=2)
Plot_ConvexHull(trt[,1][ba[,'ba.growth']>0], trt[,2][ba[,'ba.growth']>0], lty=3, lwd=2)

plot(trt[,1], trt[,2], cex=ba[,'ba.mort'], pch=21, bg=rgb(1,0,0,.9), main="Mortality")
cwm.m <- colSums(sweep(trt, 1, ba.mort, "*"))/sum(ba.mort)
#for(i in 1:nsp){
#  lty <- ifelse(ba[,'ba.mort']>0, 1, 0)
#  segments(cwm.m[1], cwm.m[2], trt[i,1], trt[i,2], col=2, lty=lty[i])
#}
points(cwm.m[1], cwm.m[2], pch='+', cex=1.5, col=rgb(1,0,0,.9))
draw.circle(cwm.m[1], cwm.m[2], radius=res$FDiv[3], lwd=2, border=5, lty=2)
draw.circle(cwm.m[1], cwm.m[2], radius=res$FDis[3], lwd=2, border=6, lty=2)
Plot_ConvexHull(trt[,1][ba[,'ba.mort']>0], trt[,2][ba[,'ba.mort']>0], lty=3, lwd=2)

plot(trt[,1], trt[,2], cex=ba[,'ba1'], pch=21, bg=rgb(0,0,1,.8), main="Time 2")
cwm1 <- colSums(sweep(trt, 1, ba1, "*"))/sum(ba1)
for(i in 1:nsp){
  lty <- ifelse(ba[,'ba1']>0, 1, 0)
  segments(cwm1[1], cwm1[2], trt[i,1], trt[i,2], col=rgb(0,0,1,.8), lty=lty[i])
}
points(cwm1[1], cwm1[2], pch='+', cex=1.5, col=rgb(0,0,1,.8))
draw.circle(cwm1[1], cwm1[2], radius=res$FDiv[4], lwd=2, border=5, lty=2)
draw.circle(cwm1[1], cwm1[2], radius=res$FDis[4], lwd=2, border=6, lty=2)
Plot_ConvexHull(trt[,1][ba[,'ba1']>0], trt[,2][ba[,'ba1']>0], lty=3, lwd=2)


plot(c(cwm0[1], cwm1[1]), c(cwm0[2], cwm1[2]), pch=21, col=c(1,4), xlim=range(trt[,1]), ylim=range(trt[,2]), cex=2, main="Shift in CWM")
points(cwm.m[1], cwm.m[2], pch='+', cex=1.5, col=rgb(1,0,0,.9))
points(cwm.g[1], cwm.g[2], pch='+', cex=1.5, col=rgb(0,1,0,.9))
arrows(cwm0[1],  cwm0[2], cwm1[1], cwm1[2], col=1, len=0.05, lwd=1)

#bg <- ifelse(ba[,'ba1'] - ba[,'ba0'] < 0, rgb(1,0,0,.8), rgb(0,1,0,.8))
#plot(trt[,1], trt[,2], cex=abs(ba[,'ba1'] - ba[,'ba0']), pch=21, bg=bg, main="Winner/loosers")

barplot(res$FRic[c(1,5,6,4)], main='FRic')
barplot(res$FEve[c(1,5,6,4)], main='FEve', col=0)
barplot(res$FDiv[c(1,5,6,4)], main='FDiv', col=5)
barplot(res$FDis[c(1,5,6,4)], main='FDis', col=6)



Plot_ConvexHull<-function(xcoord, ycoord, col=1, lwd=1, lty=1){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = col, lwd=lwd, lty=lty)
}  




#==================================

par(mfrow=c(3,3), mar=c(2,2,2,0))
plot(trt[,1], trt[,2], cex=ba[,'ba0'], col=NA, main="Time 1")
cwm0 <- colSums(sweep(trt, 1, ba0, "*"))/sum(ba0)
points(cwm0[1], cwm0[2], pch='+', cex=1.5)
points(trt[,1], trt[,2], pch=21, cex=ba[,'ba0'], bg=rgb(0,0,0,.4))
Plot_ConvexHull(trt[,1], trt[,2], lty=3, lwd=2)
draw.circle(cwm[1], cwm[2], radius=res$FDiv[1], lwd=2, border=5, lty=2)
draw.circle(cwm[1], cwm[2], radius=res$FDis[1], lwd=2, border=6, lty=2)

plot(trt[,1], trt[,2], cex=ba$ba0g, pch=21, bg=rgb(0,.6,0,.9), main="T1 + Growth")
cwm <- colSums(sweep(trt, 1, ba$ba0g, "*"))/sum(ba$ba0g)
points(cwm[1], cwm[2], pch='+', cex=1.5)
Plot_ConvexHull(trt[,1][ba$ba0g>0], trt[,2][ba$ba0g>0], lty=3, lwd=2)
draw.circle(cwm[1], cwm[2], radius=res$FDiv[5], lwd=2, border=5, lty=2)
draw.circle(cwm[1], cwm[2], radius=res$FDis[5], lwd=2, border=6, lty=2)

plot(trt[,1], trt[,2], cex=ba$ba0m, pch=21, bg=rgb(1,0,0,.7), main="T1 - Mort")
cwm <- colSums(sweep(trt, 1, ba$ba0m, "*"))/sum(ba$ba0m)
points(cwm[1], cwm[2], pch='+', cex=1.5)
Plot_ConvexHull(trt[,1][ba$ba0m>0], trt[,2][ba$ba0m>0], lty=3, lwd=2)
draw.circle(cwm[1], cwm[2], radius=res$FDiv[6], lwd=2, border=5, lty=2)
draw.circle(cwm[1], cwm[2], radius=res$FDis[6], lwd=2, border=6, lty=2)

plot(trt[,1], trt[,2], cex=ba$ba1, pch=21, bg=rgb(0,0,1,.7), main="Time 2")
cwm <- colSums(sweep(trt, 1, ba$ba1, "*"))/sum(ba$ba1)
points(cwm[1], cwm[2], pch='+', cex=1.5)
Plot_ConvexHull(trt[,1][ba$ba0m>0], trt[,2][ba$ba0m>0], lty=3, lwd=2)
draw.circle(cwm[1], cwm[2], radius=res$FDiv[4], lwd=2, border=5, lty=2)
draw.circle(cwm[1], cwm[2], radius=res$FDis[4], lwd=2, border=6, lty=2)

plot(1,1,cex=NA,axes=F)

barplot(res$FRic[c(1,5,6,4)], main='FRic')
barplot(res$FDiv[c(1,5,6,4)], main='FDiv', col=5)
barplot(res$FDis[c(1,5,6,4)], main='FDis', col=6)
barplot(res$FEve[c(1,5,6,4)], main='FEve', col=0)




























demoCWM(ba0,trait,ba.g,ba.m)
plot_demoCWM(ba0, trait, ba.g, ba.m)

demoCWM <- function(ba0, trait, ba.g, ba.m){
	cwm0 <- sum(ba0 * trait) / sum(ba0)
	cwm.g <- sum(ba.g * trait) / sum(ba.g)
	cwm.m <- sum(ba.m * trait) / sum(ba.m)
	cwm.m <- ifelse(is.nan(cwm.m), 0, cwm.m)
	ba1 <- ba0 + ba.g - ba.m
	cwm1 <- sum(ba1 * trait) / sum(ba1)
	cwm.flux <- (cwm0*sum(ba0)) + (cwm.g * sum(ba.g)) + (cwm.m * sum(ba.m))
	t0.prop <- (cwm0 * sum(ba0)) / cwm.flux
	g.prop <- (cwm.g * sum(ba.g)) / cwm.flux
	m.prop <- (cwm.m * sum(ba.m)) / cwm.flux
	out <- c(round(cwm0,3), round(cwm1,3), round(t0.prop, 3), round(g.prop, 3), round(m.prop, 3))
	names(out) <- c('cwm0', 'cwm1', 't0', 'Growth', 'Mortality')
	return(out)
}


plot_demoCWM <- function(ba0, trait, ba.g, ba.m){
	cwm0 <- sum(ba0/sum(ba0) * trait)
	ba1 <- ba0 + ba.g - ba.m
	cwm1 <- sum(ba1/sum(ba1) * trait)
	cwm.m <- sum(ba.m * trait) / sum(ba.m)
	cwm.m <- ifelse(is.nan(cwm.m), 0, cwm.m)
	cwm.g <- sum(ba.g * trait) / sum(ba.g)
	ylim <- range(cwm0, cwm1, cwm.g, cwm.m)
	# ylim <- c(min(trait)-sd(trait), max(trait)+sd(trait))
	xlim <- c(min(ba0, ba.g, ba.m), max((ba0+sum(ba0)), (ba1+sum(ba1)), (ba.g+sum(ba.g)), ba.m+sum(ba.m)))
	plot(c(sum(ba0), sum(ba1)), c(cwm0,cwm1), pch=21, ylim=ylim, xlim=xlim, cex=3, bg=c(0,'grey'))
	# points(ba0+sum(ba0), trait, pch=23, bg=0, cex=1)
	# points(ba1+sum(ba1), trait, pch=22, bg='grey', cex=1)
	# points(ba.g, trait, pch=23, bg=3, cex=1)
	# points(ba.m, trait, pch=ifelse(ba.m>0, 23, NA), bg=2)
	arrows(sum(ba0), cwm0, sum(ba1), cwm1, len=.1, lwd=3)
	# arrows(ba0+sum(ba0), trait, ba1+sum(ba1), trait, len=.1, lwd=1)
	abline(v=sum(ba0), h=cwm0, lty=3)
	points(sum(ba.g), cwm.g, pch=21, bg=3, cex=3)
	points(sum(ba.m), cwm.m, pch=21, bg=2, cex=3)
	cwm.flux <- (cwm0*sum(ba0)) + (cwm.g * sum(ba.g)) + (cwm.m * sum(ba.m))
	t0.prop <- (cwm0 * sum(ba0)) / cwm.flux
	g.prop <- (cwm.g * sum(ba.g)) / cwm.flux
	m.prop <- (cwm.m * sum(ba.m)) / cwm.flux
	mtext(paste("t0", round(t0.prop, 3), "| Growth: ", round(g.prop, 3), "| Mortality: ", round(m.prop, 3)))
	abline(h=cwm.m, v=sum(ba.m), col=2, lwd=.5)
	abline(h=cwm.g, v=sum(ba.g), col=3, lwd=.5)
}


nsp <- 20
ba0 <- rlnorm(nsp)
trait <- abs(rnorm(nsp))
ba.g <- abs(rnorm(nsp))
ba.m <- vector()
for(i in 1:nsp){ ba.m[i] <- runif(1, 0, (ba0[i] + ba.g[i])) * 0.5}

demoCWM(ba0, trait, ba.g, ba.m)

plot_demoCWM(ba0, trait, ba.g, ba.m)

plot(trait, ba0, ylim=c(0, max(ba0)), pch=16)
abline(v=trait, col=rgb(0,0,0,.2))
points(trait, ba.g, col=3, pch=16)
points(trait, ba.m, col=2, pch=16)
arrows(trait, ba.g, trait, ba.m, col=ifelse(ba.g>ba.m, 3, 2), len=.1)

plot(trait, ba.g-ba.m, pch=16, col=ifelse(ba.g>ba.m, 3, 2), cex=2); abline(h=0, lty=3)


############################
###  DYNAMIC SIMULATION  ###
############################
nsp <- 100
nruns <- 10
out <- matrix(nrow=nruns, ncol=5)
ba0 <- rlnorm(nsp)
trait <- abs(rnorm(nsp))
for(i in 1:nruns){
	ba.g <- abs(rnorm(nsp))
	mort <- runif(nsp,0,3)
	ba.m <- ifelse(mort > (ba.g + ba1), 0, mort)
	if( i == 1) {
		out[i,] <- demoCWM(ba0, trait, ba.g, ba.m)
		ba1 <- ba0 + ba.g - ba.m
		}
	if( i > 1){
		out[i,] <- demoCWM(ba0=ba1, trait, ba.g, ba.m)
 		ba1 <- ba1 + ba.g - ba.m
 		}
}
colnames(out) <- names(demoCWM(ba0, trait, ba.g, ba.m))
out

# plot(out[,1], type='l', ylim=c(0,max(out)))
# points(out[,3], type='l', col=4)
plot(out[,3], type='l', ylim=c(0,1))
points(out[,4], type='l', col=3)
points(out[,5], type='l', col=2)
legend('topleft', colnames(out[,c(3:5)]), pch=21, pt.bg=c(1,3,2))




ba0 <- rlnorm(100)
g <- abs(rnorm(100))
m <- runif(100,0,3)
trait <- abs(rnorm(100))
demographic_CWM(ba0, g, m, trait)


ba0 <- c(1,1)
g <- c(2,1)
m <- c(1,2)
trait <- c(1,5)
demographic_CWM(ba0, g, m, trait)









# HOW MUCH DID CWM TRAIT VALUES CHANGE BETWEEN T0 AND T1?
	# - how do CWM trait values change 
	
# WHAT IS THE RELATIVE CONTRIBUTION TO THE CHANGE IN CWM FROM GROWTH VERSUS MORTALITY?
# HOW DO THESE VALUES DIFFER IF WE ADD COLONIZATION AND EXTINCTION?
# So, if growth is most important contributor to change in CWM, who is growing?  We also need to relate average demographic rates to species traits at each time step (i.e., high SLA species are growing fast early on and then have high mortality later in succession)

#We see a large amount of variability between years and between sites in studies of forest succession.  Do we see regularity in the demographic drivers of functional shifts?

#Linking theory of forest succession with functional traits requires a stronger link between the demographic drivers of shifts in functional composition of communities.  These have largely been infered.  Here, we introduce a novel method for partitioning the relative contribution of observed shifts in CWM trait values to different demographic process.  This method will help better connect trait-based ecology with existing theoretical understanding of forest dynamics.




demographic_CWM <- function(ba0, g, r, c, m, e, trait){
	rba0 <- ba0/sum(ba0)
	cwm0 <- sum(rba0 * trait)

	rba.growth <- g/sum(g)
	cwm.g <- sum(rba.growth * trait)
	rba.recruit <- r/sum(r)
	cwm.r <- sum(rba.recruit * trait)
	rba.colon <- c/sum(c)
	cwm.c <- sum(rba.colon * trait)
	rba.mort <- m/sum(m)
	cwm.m <- sum(rba.mort * trait)
	rba.ext <- e/sum(e)
	cwm.e <- sum(rba.ext * trait)

	ba1 <- ba0 + g + r + c - m - e
	rba1 <- ba1/sum(ba1)
	cwm1 <- sum(rba1 * trait)

	cwm.flux <- (cwm0*sum(ba0)) + (cwm.g * sum(g)) + (cwm.r * sum(r)) + (cwm.c * sum(c)) + (cwm.m * sum(m)) + (cwm.e * sum(e))

	t0.prop <- (cwm0 * sum(ba0)) / cwm.flux
	g.prop <- (cwm.g * sum(g)) / cwm.flux
	r.prop <- (cwm.r * sum(r)) / cwm.flux
	c.prop <- (cwm.c * sum(c)) / cwm.flux
	m.prop <- (cwm.m * sum(m)) / cwm.flux
	e.prop <- (cwm.e * sum(e)) / cwm.flux

	out <- c(round(cwm0,3), 
				round(cwm1,3), 
				round(t0.prop, 3), 
				round(g.prop, 3), 
				round(r.prop, 3), 
				round(c.prop, 3), 
				round(m.prop, 3),
				round(e.prop, 3))
	names(out) <- c('cwm0', 'cwm1', 't0', 'Growth','Recruitment','Colonization','Mortality','Extinction')
	return(out)
}


ba0 <- rlnorm(100)
g <- abs(rnorm(100))
r <- abs(rnorm(100, .1))
c <- abs(rnorm(100, .3))
m <- runif(100,0,3)
e <- sample(c(0,1), 100, c(.9,.1), replace=T)
trait <- abs(rnorm(100))
demographic_CWM(ba0, g, r, c, m, e, trait)










#====================
# OTHER STUFF...


nsp <- 10
ba0 <- rlnorm(nsp)
trait <- abs(rnorm(nsp))
ba.g <- abs(rnorm(nsp))
ba.m <- vector()
for(i in 1:nsp){ ba.m[i] <- runif(1, 0, (ba0[i] + ba.g[i]))}


ba0 <- c(200,500)
trait <- c(0,1)
ba.g <- c(100,500)
ba.m <- c(1,1)

sum(ba.m*trait) / sum(ba.m)
sum(ba0*trait) / sum(ba0)

x0 <- sum(ba0)
y0 <- sum(ba0 * trait) / sum(ba0)

xg <- sum(ba0) + sum(ba.g)
yg <- sum(ba.g * trait) / sum(ba.g)

xm <- sum(ba0) - sum(ba.m)
ym <- sum(ba.m * trait) / sum(ba.m)
#ym <- (y0-ym) + y0

xr <- sum(ba0 + ba.g - ba.m)
yr <- sum(yg + (y0-ym) + y0 - y0)

#plot(c(x0, xm, xg, (xg - sum(ba.m))), c(y0, ym, yg, (y0 + yg - ym)), pch=21, bg=1:4, cex=c(3,2,2,2))
plot(c(x0, xm, xg, xr), c(y0, ym, yg, yr), pch=21, bg=1:4, cex=c(4,2,2,3))

arrows(x0, y0, xg, yg, col=3, len=.1, lwd=2)
arrows(x0, y0, xm, ym, col=2, len=.1, lwd=2)
abline(v=x0, h=y0, lty=3)
#arrows(x0, y0, (xg - sum(ba.m)), (y0 + yg - ym), len=.1, col=4)
arrows(x0, y0, xr, yr, len=.1, col=4, lwd=3)
arrows(xg, yg, xr, yr, col=2, len=.1, lwd=2)

xgcontrib <- abs(x0 - xg) / abs(x0 - (x0 + xg - xm))
xmcontrib <- abs(x0 - xm) / abs(x0 - (x0 + xg - xm))
ygcontrib <- abs(y0 - yg) / abs(y0 - (y0 + yg - ym))
ymcontrib <- abs(y0 - ym) / abs(y0 - (y0 + yg - ym))

gcon <- (xgcontrib + ygcontrib) / (xgcontrib + ygcontrib + xmcontrib + ymcontrib)
mcon <- (xmcontrib + ymcontrib) /(xgcontrib + ygcontrib + xmcontrib + ymcontrib)
mtext(paste("Growth:", round(gcon,3), "| Mortality:", round(mcon,3)))




#========================================
#========================================
#========================================


CWMdemog <- function(ba0, trait, ba.g, ba.m){

	x0 <- sum(ba0)
	y0 <- sum(ba0 * trait) / sum(ba0)
	
	xg <- sum(ba0) + sum(ba.g)
	yg <- sum(ba.g * trait) / sum(ba.g)
	
	xm <- sum(ba0) - sum(ba.m)
	ym <- sum(ba.m * trait) / sum(ba.m)

	xr <- sum(ba0 + ba.g - ba.m)
#	yr <- sum(yg + (y0-ym) + y0 - y0)
	ba1 <- ba0 + ba.g - ba.m
	yr <- sum(ba1/sum(ba1) * trait)

	xgcontrib <- abs(x0 - xg) / abs(x0 - (x0 + xg - xm))
	xmcontrib <- abs(x0 - xm) / abs(x0 - (x0 + xg - xm))
	ygcontrib <- abs(y0 - yg) / abs(y0 - (y0 + yg - ym))
	ymcontrib <- abs(y0 - ym) / abs(y0 - (y0 + yg - ym))
	gcon <- (xgcontrib + ygcontrib) / (xgcontrib + ygcontrib + xmcontrib + ymcontrib)
	mcon <- (xmcontrib + ymcontrib) /(xgcontrib + ygcontrib + xmcontrib + ymcontrib)

	return(list(cwm0=y0, cwm1=yr, ba0=x0, ba1=xr, growth.contribution=gcon, mortality.contribution=mcon))

# THIS PROVIDES CONTRIBUTION IN UNITS OF TRAITS
	# abs(y0-yr) * gcon
	# abs(y0-yr) * mcon

}

nsp <- 100
ba0 <- rlnorm(nsp)
trait <- abs(rnorm(nsp))
ba.g <- abs(rnorm(nsp))
ba.m <- vector()
for(i in 1:nsp){ ba.m[i] <- runif(1, 0, (ba0[i] + ba.g[i]) * 0.0001) }

res <- CWMdemog(ba0, trait, ba.g, ba.m)
res

	

	
	
	#plot(c(x0, xm, xg, (xg - sum(ba.m))), c(y0, ym, yg, (y0 + yg - ym)), pch=21, bg=1:4, cex=c(3,2,2,2))
plot(c(x0, xm, xg, xr), c(y0, ym, yg, yr), pch=21, bg=1:4, cex=c(4,2,2,3))

arrows(x0, y0, xg, yg, col=3, len=.1, lwd=2)
arrows(x0, y0, xm, ym, col=2, len=.1, lwd=2)
abline(v=x0, h=y0, lty=3)
#arrows(x0, y0, (xg - sum(ba.m)), (y0 + yg - ym), len=.1, col=4)
arrows(x0, y0, xr, yr, len=.1, col=4, lwd=3)
arrows(xg, yg, xr, yr, col=2, len=.1, lwd=2)

xgcontrib <- abs(x0 - xg) / abs(x0 - (x0 + xg - xm))
xmcontrib <- abs(x0 - xm) / abs(x0 - (x0 + xg - xm))
ygcontrib <- abs(y0 - yg) / abs(y0 - (y0 + yg - ym))
ymcontrib <- abs(y0 - ym) / abs(y0 - (y0 + yg - ym))

gcon <- (xgcontrib + ygcontrib) / (xgcontrib + ygcontrib + xmcontrib + ymcontrib)
mcon <- (xmcontrib + ymcontrib) /(xgcontrib + ygcontrib + xmcontrib + ymcontrib)
mtext(paste("Growth:", round(gcon,3), "| Mortality:", round(mcon,3)))



