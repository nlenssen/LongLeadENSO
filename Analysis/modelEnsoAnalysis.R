modelNameVec <- c('CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5',
	'CESM2', 
	'GISS-E21G', 
	'CESM1-NMME', 
	'CCSM4-NMME',
	'CM2.1-NMME',
	'CM2.5-NMME',
	'Observations')

nModels <- length(modelNameVec)

# Set the experiment parameters to play nice with exisiting code
initModelNameVec <- modelNameVec
initExperimentName <- 'tos_zos_NEW'
trialFields <- c('tos', 'zos')

# the bounds for the Nino/Nina Hovemollers
hovLonBounds <- c(120, 270)
hovLatBounds <- c(-3,3)

timeBounds <- c(-12, 24)
timeBoundVec <- min(timeBounds):max(timeBounds)

plotdir <- '~/Desktop/scratchFigures/ensoPaper/Hovemollers'
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)


###############################################################################
# Loop over all models to get the rate of 2x LN events
###############################################################################

nYearsVec     <- rep(NA, nModels)
djfEventCount <- matrix(NA, 3, nModels)
twoYearProb   <- matrix(NA, 2, nModels)
threeYearProb <- matrix(NA, 2, nModels)

for(mod in 1:nModels){
	# set the two model names
	modelName <- modelNameVec[mod]
	initModelName <- initModelNameVec[mod]

	# run the analyses
	source('Namelists/baseNamelist.Rnl')

	load(sprintf('%s/step01_KeyInfo.Rda',wdir))

	djfEvents <- eventVec[which(timeMap[,2]==1)]
	nYears    <- length(djfEvents)

	# save number of years
	nYearsVec[mod] <- nYears

	djfEventCount[,mod] <- table(djfEvents)/length(djfEvents)

	for(j in c(1,3)){
		# get the out vec ind
		ind <- ifelse(j==1,1,2)
		
		# get all of the years where two year events happen
		eventInds <- which(djfEvents==j-2)
		count <- 0
		for(k in eventInds){
			if(k+1 > length(djfEvents)) next()
			if(djfEvents[k+1] == (j-2)) count <- count + 1
		}

		twoYearProb[ind,mod] <- count/nYears

		# get all of the years where three year events happen
		count <- 0
		for(k in eventInds){
			if(k+1 > length(djfEvents)) next()
			if(k+2 > length(djfEvents)) next()
			if(djfEvents[k+1] == (j-2) & djfEvents[k+2] == (j-2)) count <- count + 1
		}

		threeYearProb[ind,mod] <- count/nYears
	}
}


vizProb <- t(djfEventCount)
rownames(vizProb) <- modelNameVec
colnames(vizProb) <- c('LN', 'NN', 'EN')
round(vizProb,3)

viz2 <- t(twoYearProb)
rownames(viz2) <- modelNameVec
colnames(viz2) <- c('2x LN', '2x EN')
round(viz2,3)

viz3 <- t(threeYearProb)
rownames(viz3) <- modelNameVec
colnames(viz3) <- c('3x LN', '3x EN')
round(viz3,3)

###############################################################################
# Loop over all models to do the hovemollers
###############################################################################
zr <- pm(2.5)
pal <- rev(brewer.pal(11,'RdBu'))



for(mod in 1:nModels){

# set the two model names
modelName <- modelNameVec[mod]
initModelName <- initModelNameVec[mod]

# run the analyses
source('Namelists/baseNamelist.Rnl')



# load the step 1 info
load(sprintf('%s/step01_KeyInfo.Rda',wdir))


# get the lon lat inds (inclusive for now)
lonInds <- which(lon >= min(hovLonBounds) & lon <= max(hovLonBounds))
latInds <- which(lat >= min(hovLatBounds) & lat <= max(hovLatBounds))

# load in the subset data
handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,'tos',modelName))

sstField <- ncvar_get(handle,'tos_anom', start=c(min(lonInds),min(latInds),1),
										 count=c(length(lonInds),length(latInds),-1))

nc_close(handle)


# get all the dfj EN events
enInds <- which(timeMap[,2] == 1 & eventVec == 1)

# get the evolution of all of these from 12 months prior to 24 months afterwards
compositeArray <- array(NA, dim=c(nrow(sstField), ncol(sstField), length(timeBoundVec), length(enInds)))

for(i in 1:length(enInds)){
	tempInds <- (enInds[i] + min(timeBounds)):(enInds[i] + max(timeBounds))

	if(min(tempInds) < 1 | max(tempInds) > nrow(timeMap)) next()

	compositeArray[,,,i] <- sstField[,,tempInds]
}

enHovemoller <- apply(compositeArray, c(1,3), mean, na.rm=T)




# get all the dfj EN events
lnInds <- which(timeMap[,2] == 1 & eventVec == -1)

# get the evolution of all of these from 12 months prior to 24 months afterwards
compositeArray <- array(NA, dim=c(nrow(sstField), ncol(sstField), length(timeBoundVec), length(enInds)))

for(i in 1:length(lnInds)){
	tempInds <- (lnInds[i] + min(timeBounds)):(lnInds[i] + max(timeBounds))

	if(min(tempInds) < 1 | max(tempInds) > nrow(timeMap)) next()

	compositeArray[,,,i] <- sstField[,,tempInds]
}

lnHovemoller <- apply(compositeArray, c(1,3), mean, na.rm=T)






pdf(sprintf('%s/hov_%s.pdf', plotdir, modelName),10,6)
par(mfrow = c(1,2))
image.plot(lon[lonInds], timeBoundVec, enHovemoller, zlim=zr, col=pal,
	yaxp = c(-12,24, 6),
	xlab='Longitude', ylab='Months',
	main=sprintf('El Niño Evolution (%s piControl)', modelName))

abline(h=c(-12,0,12,24), lty=c(3,1,3,3))

image.plot(lon[lonInds], timeBoundVec, lnHovemoller, zlim=zr, col=pal,
	yaxp = c(-12,24, 6),
	xlab='Longitude', ylab='Months',
	main=sprintf('La Niña Evolution (%s piControl)', modelName))

abline(h=c(-12,0,12,24), lty=c(3,1,3,3))
dev.off()

}





