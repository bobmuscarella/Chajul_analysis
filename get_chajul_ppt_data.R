setwd("/Users/Bob/Projects/Neoselvas/Data/chajul/chajul ppt data")
library(raster)

# These give precipitation in mm/10 day period
# # # month <- unlist(lapply(strsplit(list.files(), "_"), function(x) x[2]))
# # # year <- unlist(lapply(strsplit(list.files(), "_"), function(x) x[3]))
# # # days <- unlist(lapply(lapply(strsplit(list.files(), "_"), function(x) x[1]), function(x) substring(x, 5,100)))

# # # coords <- cbind(-90.99825926, 16.09176667)

# # # ppt <- vector()
# # # for (i in 1:length(list.files())){
	# # # ppt <- c(ppt, extract(raster(list.files()[[i]]), coords))
# # # }

# # # rows <- 1:length(days)
# # # x <- as.data.frame(cbind(days, month, year, ppt), row.names=rows)
# # # x$ppt <- as.numeric(as.character(x$ppt))

# # # m <- as.data.frame(cbind(1:12, c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')))
# # # x$m <- m$V1[match(x$month, m$V2)]

# # # d <- as.data.frame(cbind(c(1,2,3,3,3,3), c('1-10','11-20','21-28','21-29','21-30','21-31')))
# # # x$d <- d$V1[match(x$day, d$V2)]

# # # x <- x[order(x$year, x$m, x$d),]
# # # rownames(x) <- NULL

# # # x$ym <- paste(x$year, x$m, 1, sep='-')
# # # x$date <- as.Date(paste(x$year, x$m, x$d, sep='-'))

# # # head(x)
# # # ppt <- x
# # # save(ppt, file="Chajul_ppt_data.RDA")

setwd("/Users/Bob/Projects/Neoselvas/Chajul/DATA/chajul ppt data")
load("Chajul_ppt_data.RDA")
head(ppt)

monthlysums <- tapply(ppt$ppt, ppt$ym, sum)
monthlysums <- monthlysums[order(as.Date(names(monthlysums)))]
yearlysums <- tapply(ppt$ppt, ppt$year, sum)

map <- tapply(ppt$ppt, ppt$year, sum)

labs <- paste("1/1/", 1998:2013, sep='')
plot(monthlysums, type='l', ylim=c(0, max(yearlysums)), axes=F, xlab="", ylab='Precipitation (mm)', main="Monthly and annual precipitation at Chajul")
axis(1, labels=labs, at=seq(1, length(monthlysums),12), las=2)
axis(2)
points(cbind(seq(13,193,12), map), type='p', bg=1, pch=21)
polygon(x=c(13, seq(13,193,12), 193), y=c(0,map,0), col=rgb(0,0,1,.75))
polygon(x=c(1,1:length(monthlysums),length(monthlysums)), y=c(0, monthlysums, 0), col=3)

range(map)

# abline(v=seq(1,300,12), col=2)
segments(x0=seq(13,300,12), y0=0, x1=seq(13,300,12), y1=map, lwd=1)

map[3:16]











