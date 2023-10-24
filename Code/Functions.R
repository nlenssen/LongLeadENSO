###############################################################################
# Analogue Forecasting Functions
###############################################################################

# Error function for the simple MSE error of Ding et al. (2018)
simpleMseError <- function(ind, trialMat, targetMat){
	
	sum(weightingMat*(targetMat/targetSigma[ind] - trialMat/trialSigma[ind])^2,na.rm=T)

}


# The big function that creates an analogue forecast off of a single initial
# state as passed in by the initalState list of the form
# initialState <- list(lon=lon[lonInds],
# 					 lat=lat[latInds],
# 					 fields=initialStateFields, # an array of dim [lon,lat,field]
# 					 fieldNames=c('tos', 'zos'),
# 					 timeMap=library$timeMap,
# 					 mask=maximalMask)

createAnalogueForecastGeneral <- function(forecastInd, initialState, libraryInds,
											errorFunction = simpleMseError,
											fieldWeights = NULL,
											forecastFieldNames = NULL,
											libIndsExclude = NULL,
											obsHindcast = FALSE,
											verificationModelName = NULL)
											{
	##########
	# Set up some obs
	##########
	if(is.null(fieldWeights)) fieldWeights <- rep(1, length(initialState$fieldNames))

	if(is.null(forecastFieldNames)) forecastFieldNames <- initialState$fieldNames

	if(is.null(verificationModelName)) verificationModelName <- modelName

	nTestFields <- length(initialState$fieldNames)
	nForecastFields <- length(forecastFieldNames)


	################
	# Tests on input
	################

	if(length(initialState$lon) != dim(initialState$fields)[1] |
	   length(initialState$lat) != dim(initialState$fields)[2]){
		stop('Mismatch of lon/lat size in initial state')
	}

	if(nTestFields != dim(initialState$fields)[3]){
		stop('Mismatch of number of fields in initial state')
	}


	##########
	# Determine the inds to search the library over
	##########

	# get the max ind to pull from to avoid NA values by forecasting off the tail of the data
	clippedLibraryInds <- libraryInds[libraryInds < nrow(timeMap) - forecastLength]

	# find the closest anagloues over the indo pacific box used in Ding 2018 method
	seasonOfInterest <- initialState$timeMap[forecastInd,2]
	loopInds <- intersect(which(timeMap[,2] == seasonOfInterest), clippedLibraryInds)

	# remove possible library states to search (mostly not in use as the libaryInds is passed in)
	removeInds <- which(loopInds %in% libIndsExclude)
	if(length(removeInds) > 0){
		loopInds <- loopInds[-removeInds]
	}

	##########
	# Get cosine weighting over the reduced domain
	##########

	# need the cos weighting mat over the domain
	cosMatFull <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180)) *
									initialState$mask
	cosMat <- cosMatFull[lonInds,latInds]
	weightingMat <- cosMat/sum(cosMat,na.rm=T)


	##########
	# Loop over the library to determine the errorVec corresponding to the loopInds variable
	##########

	errorVec <- rep(NA, length(loopInds))

	for(i in 1:length(loopInds)){
		
		partialError <- rep(NA, nTestFields)
		
		for(j in 1:nTestFields){

			varName <- initialState$fieldNames[j]

			handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',
				ddir,modelName,varName,modelName))

			trialMat <- ncvar_get(handle,sprintf('%s_anom',varName),
								start=c(min(lonInds),min(latInds),loopInds[i]),
								count=c(length(lonInds),length(latInds),1))*maximalMask[lonInds,latInds]
			nc_close(handle)

			# possibility to add a ... or other args here to pass more things to the error function
			partialError[j] <- do.call(errorFunction, list(ind=j,
														   trialMat=trialMat,
														   targetMat=initialState$fields[,,j]))
		}

		errorVec[i] <- sum(partialError*fieldWeights)
	}


	##########
	# Pull the forecast fields
	##########


	# select the ensemble indices that have the least error and place in order
	rankedError <-rank(errorVec)
	ensembleInds <- rep(NA, ensembleSize)
	for(i in 1:ensembleSize) ensembleInds[i] <- which(rankedError==i)

	# pull the sst fields of the ensemble and get the ensembleMean
	# 5-dimension array of[lon,lat,forecastTime,ensemble,field]
	forecastArray <- array(NA,
		dim=c(length(lonInds),length(latInds),2*forecastLength+1,ensembleSize,nForecastFields))

	for(i in 1:ensembleSize){
		for(j in 1:nForecastFields){

			startInd <- max(1, loopInds[ensembleInds[i]] - (forecastLength+ 1))

			handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,forecastFieldNames[j],modelName))
			forecastArray[,,,i,j] <- ncvar_get(handle,sprintf('%s_anom',forecastFieldNames[j]),
								start=c(min(lonInds),min(latInds),startInd),
								count=c(length(lonInds),length(latInds),forecastLength*2+1))
			nc_close(handle)

		}
	}

	# get the ensemble mean forecast as well
	ensembleMeanForecast <- array(NA,
		dim=c(length(lonInds),length(latInds),2*forecastLength+1,nForecastFields))

	for(j in 1:nForecastFields){
		ensembleMeanForecast[,,,j] <- arrMult(apply(forecastArray[,,,,j],c(1,2,3),mean),maximalMask[lonInds,latInds])
	}



	##########
	# Pull the verification fields (if running against a model)
	##########

	if(!obsHindcast){

		# get some important metadata
		handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,verificationModelName,forecastFieldNames[1],modelName))
		nt <- length(ncvar_get(handle,'time'))
		nc_close(handle)


		# get the corresponding verification fields as well
		verificationLength <- min(length(forecastInd:nt),forecastLength)

		targetFields <- array(NA,
			dim=c(length(lonInds),length(latInds),2*forecastLength+1,nForecastFields))

		# loop through all verificiation fields
		for(j in 1:nForecastFields){

				startInd <- max(1, loopInds[ensembleInds[i]] - (forecastLength+ 1))


				handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,forecastFieldNames[j],modelName))
				targetFields[,,,j] <- ncvar_get(handle,sprintf('%s_anom',forecastFieldNames[j]),
									start=c(min(lonInds),min(latInds),startInd),
									count=c(length(lonInds),length(latInds),2*verificationLength+1))
				nc_close(handle)

		}

	}


	##########
	# Pull the verification fields (if running observational hindcasts, need to debug potentially)
	##########

	if(obsHindcast){
		# get the corresponding verification fields as well
		verificationMaxInd <- min(forecastInd+forecastLength, nrow(obsDataList$timeMap))

		if(length(forecastInd:obsMaxInd)==1){
			targetFields         <- obsDataList$tos[lonInds,latInds,forecastInd:obsMaxInd] *
										 maximalMask[lonInds,latInds]
		} else{
			targetFields         <- arrMult(obsDataList$tos[lonInds,latInds,forecastInd:obsMaxInd],
										maximalMask[lonInds,latInds])
		}

		ensembleMeanForecast <- arrMult(apply(forecastArray,c(1,2,3),mean),maximalMask[lonInds,latInds])
	}


	##########
	# function output
	##########

	# package and send out
	outList <- list(targetFields=targetFields,ensembleMeanForecast=ensembleMeanForecast,forecastArray=forecastArray,
						errorVec=errorVec[ensembleInds],ensembleInds=ensembleInds)


	return(outList)


}

# The big function that creates an analogue forecast off of a single initial
# state as passed in by the initalState list of the form
# initialState <- list(lon=lon[lonInds],
# 					 lat=lat[latInds],
# 					 fields=initialStateFields, # an array of dim [lon,lat,field]
# 					 fieldNames=c('tos', 'zos'),
# 					 timeMap=library$timeMap,
# 					 mask=maximalMask)


# forecastInd = testingInds[i]
# initialState = initialState
# libraryInds = libraryInds
# errorFunction = simpleMseError
# forecastFieldNames = forecastFieldNames
# libIndsExclude = testingInds[i]
# analogueSize = analogueSize
# fieldWeights = NULL
createAnalogueForecastMinimal <- function(forecastInd, initialState, libraryInds,
											errorFunction = simpleMseError,
											fieldWeights = NULL,
											forecastFieldNames = NULL,
											libIndsExclude = NULL,
											analogueSize = 50)
											{
	##########
	# Set up some obs
	##########
	if(is.null(fieldWeights)) fieldWeights <- rep(1, length(initialState$fieldNames))

	if(is.null(forecastFieldNames)) forecastFieldNames <- initialState$fieldNames

	nTestFields <- length(initialState$fieldNames)
	nForecastFields <- length(forecastFieldNames)


	################
	# Tests on input
	################

	if(length(initialState$lon) != dim(initialState$fields)[1] |
	   length(initialState$lat) != dim(initialState$fields)[2]){
		stop('Mismatch of lon/lat size in initial state')
	}

	if(nTestFields != dim(initialState$fields)[3]){
		stop('Mismatch of number of fields in initial state')
	}


	##########
	# Determine the inds to search the library over
	##########

	# get the max ind to pull from to avoid NA values by forecasting off the tail of the data
	clippedLibraryInds <- libraryInds[libraryInds < nrow(timeMap) - forecastLength]

	# find the closest anagloues over the indo pacific box used in Ding 2018 method
	seasonOfInterest <- initialState$timeMap[forecastInd,2]
	
	loopInds <- intersect(which(timeMap[,2] == seasonOfInterest), clippedLibraryInds)

	# remove possible library states to search (mostly not in use as the libaryInds is passed in)
	removeInds <- which(loopInds %in% libIndsExclude)
	if(length(removeInds) > 0){
		loopInds <- loopInds[-removeInds]
	}

	##########
	# Get cosine weighting over the reduced domain
	##########

	# need the cos weighting mat over the domain
	cosMatFull <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180)) *
									initialState$mask
	cosMat <- cosMatFull[lonInds,latInds]
	weightingMat <- cosMat/sum(cosMat,na.rm=T)


	##########
	# Loop over the library to determine the errorVec corresponding to the loopInds variable
	##########

	errorVec <- rep(NA, length(loopInds))

	for(i in 1:length(loopInds)){
		
		partialError <- rep(NA, nTestFields)
		
		for(j in 1:nTestFields){

			varName <- initialState$fieldNames[j]

			handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',
				ddir,modelName,varName,modelName))

			trialMat <- ncvar_get(handle,sprintf('%s_anom',varName),
								start=c(min(lonInds),min(latInds),loopInds[i]),
								count=c(length(lonInds),length(latInds),1))*maximalMask[lonInds,latInds]
			nc_close(handle)

			# possibility to add a ... or other args here to pass more things to the error function
			partialError[j] <- do.call(errorFunction, list(ind=j,
														   trialMat=trialMat,
														   targetMat=initialState$fields[,,j]))


		}

		# rm(trialMat)
		# gc()

		errorVec[i] <- sum(partialError*fieldWeights)
	}



	# select the ensemble indices that have the least error and place in order
	rankedError <-rank(errorVec)

	# get the position of the ensemble in the full library
	errorVecOut <- rep(NA, analogueSize)
	ensembleInds <- rep(NA, analogueSize)

	# get the top analogueSize matches and save their error value as well as
	# the index of their position in the library
	for(i in 1:analogueSize){
		tempInd <- which(rankedError==i)

		# get the error
		errorVecOut[i] <- errorVec[tempInd]
		
		# get the position in the library
		ensembleInds[i] <- loopInds[tempInd]

	}


	##########
	# function output
	##########

	# package and send out
	outList <- list(errorVec=errorVecOut, ensembleInds=ensembleInds)



	return(outList)


}


###############################################################################
# ENSO Statistics Functions
###############################################################################
probExceedence <- function(obj, threshold, greater=TRUE){
	if(greater){
		sum(obj > threshold, na.rm=T)/sum(!is.na(obj))
	} else{
		sum(obj < threshold, na.rm=T)/sum(!is.na(obj))
	}
}

###############################################################################
# Initialized Forecast Processing Functions
###############################################################################

###############################################################################
# Forecast Synthesis/Verification Functions
###############################################################################

# calculate the ROC of a forecast/obs pair
rocCalculationSingle <- function(forecastMat, obsVec, probStep=0.05){
	uniqueProbs <- seq(0,1,by=probStep)

	# get some length values
	nForecasts <- nrow(forecastMat)

	# build objects to save things
	rocArray     <- array(NA, dim=c(length(uniqueProbs),3,2))
	numWarnMat <- array(NA, dim=c(length(uniqueProbs),3))

	# these don't depend on probability
	numEvents    <- rep(NA, 3)
	numNonEvents <- rep(NA, 3)

	for(j in 1:3){
		# make sure we only have good inds
		goodInds <- which(!is.na(forecastMat[,j]) & !is.na(obsVec))
		

		# pull the forecasts for a given event
		tempFcasts <- forecastMat[goodInds,j]
		tempObs    <- obsVec[goodInds]

		# calculate the actual event counts (e e' in Mason Graham 2002)
		numEvents[j]    <- sum(tempObs == (j-2))
		numNonEvents[j] <- length(tempObs) - numEvents[j]

		for(i in 1:length(uniqueProbs)){
			warnInds <- which(tempFcasts >= uniqueProbs[i])
			numWarnMat[i,j] <- length(warnInds)

			if(length(warnInds)==0) next()

			hr <- sum(tempObs[warnInds] == (j-2))/sum(tempObs == (j-2))
			far <- sum(tempObs[warnInds] != (j-2))/sum(tempObs != (j-2))
			rocArray[i,j,] <- c(hr,far)
		}
		
	}


	# calculate the roc score (area under the curve with hr on y, far on x)
	rocScoreVec <- rep(NA, 3)

	for(j in 1:3){
		x <- c(rocArray[,j,2],0)
		y <- c(rocArray[,j,1],0)

		x <- x[which(!is.na(x))]
		y <- y[which(!is.na(y))]

		deltaX <- -diff(x)
		deltaY <- -diff(y)
		rocScoreVec[j] <- sum(deltaX * y[-1] + 1/2 * deltaX * deltaY) - 0.5
	}

	# some signifgance work following equations (7) and (12) of Mason and Graham (2002)
	pvalVec <- rep(NA, 3)
	sdVec   <- rep(NA, 3)

	for(j in 1:3){
		e      <- numEvents[j]
		ePrime <- numNonEvents[j]

		A      <- rocScoreVec[j] + 0.5

		U <- e * ePrime * (1-A)

		# null distribution
		uNullMean <- e*ePrime/2
		uNullSd <- sqrt(e*ePrime*(e+ePrime+1)/12)

		sdVec[j]   <- uNullSd
		pvalVec[j] <- pnorm(U,uNullMean,uNullSd)
	}

	return(list(roc=rocArray,rocScore=rocScoreVec,pval=pvalVec,sd=sdVec, n=nForecasts))
}

# calculation of the rpss
rpssCalculationSingle <- function(forecastMat, obsVec){
	nForecasts <- nrow(forecastMat)

	# store the RPS values at each point
	rpsForecast <- rep(NA, nForecasts)
	rpsBaseline <- rep(NA, nForecasts)

	# climatology forecast
	climatology <- c(1/4,1/2,1/4)
	cumClimatology <- cumsum(climatology)

	# loop through  everything
	for(tInd in 1:nForecasts){
		tempForecast <- forecastMat[tInd,]
		
		tempObs <- rep(0, 3)	
		tempObs[obsVec[tInd] + 2] <- 1

		cumForecast   <- cumsum(tempForecast)
		cumOutcome    <- cumsum(tempObs)

		rpsForecast[tInd] <- sum((cumForecast - cumOutcome)^2)

		rpsBaseline[tInd] <- sum((cumClimatology - cumOutcome)^2)
	}


	# calculate the RPSS
	rpss <- 1 - mean(rpsForecast)/mean(rpsBaseline)

	return(list(rpss=rpss, n=nForecasts))
}


# calculate the RPSS and decomposition for a forecastMat and obsVec
rpssDecomp <- function(forecastMat, obsVec, probStep=0.05){
	nForecasts <- nrow(forecastMat)


	# climatology forecast
	climForecast <- matrix(c(1/4,1/2,1/4),nrow=nrow(forecastMat),ncol=3,byrow=T)


	####################################################################
	# Decompose the Brier score into three parts following [Murphy 1973]
	# using the notation of [Stephenson et al. 2008]
	#
	# Note: this agrees with the RPS calculation in the RPSS function
	####################################################################
	# Break forecast into two brier calculations [Make outcome 0/1!!!!!]
	lowForecast <- cbind(forecastMat[,1],forecastMat[,2] + forecastMat[,3])
	lowClim     <- cbind(climForecast[,1],climForecast[,2] + climForecast[,3])
	lowOutcome  <- obsVec == 0 | obsVec == 1

	highForecast <- cbind(forecastMat[,1] + forecastMat[,2], forecastMat[,3])
	highClim     <- cbind(climForecast[,1] + climForecast[,2], climForecast[,3])
	highOutcome  <- obsVec == 1

	lowDecomp     <- brierDecomp(lowForecast,lowOutcome,probStep)
	highDecomp    <- brierDecomp(highForecast,highOutcome,probStep)	

	lowClimDecomp  <- brierDecomp(lowClim,lowOutcome,probStep)
	highClimDecomp <- brierDecomp(highClim,highOutcome,probStep)

	rps0 <- lowDecomp + highDecomp
	rpsClim <- lowClimDecomp + highClimDecomp

	rpsSeries <- rps0[1]
	relSeries <- rps0[2]
	resSeries <- rps0[3]
	uncSeries <- rps0[4]

	rpssSeries   <- 1 - rps0[1]/rpsClim[1]
	relssSeries <- 1 - rps0[2]/rpsClim[2]

	decompList <- list(rps = rpsSeries, rel = relSeries, res = resSeries, unc = uncSeries)
	skillList    <- list(rpss=rpssSeries, relss=relssSeries)

	outList <- list(decomp=decompList,skill=skillList)

	return(outList)
}

# calculate the BSS, reliability, resolution, and uncertainty for a 
# single binary forecast
# Brier decomposition according to Murphy 1973 following the notation of
# Stephenson et al. 2008.
brierDecomp <- function(forecast,outcome, probStep){
	n <- nrow(forecast)

	# calculate dBar (or oBar in 2008 notation)
	dBar <- 1/sum(outcome)

	# make bins
	bins <- seq(0,1,by=probStep)
	centers <- bins[-1]-(probStep)/2

	# determine the bin counts
	nk <- rep(NA, length(centers))
	oBarK <- rep(NA, length(centers))


	# Do first bin outside loop to include zero
	tempInds <- which(forecast[,2]>=bins[1] & forecast[,2]<=bins[2])
	nk[1] <- length(tempInds)
	oBarK[1] <- mean(outcome[tempInds])

	# loop through the rest of the bins
	for(k in 2:length(centers)){
		tempInds <- which(forecast[,2]>bins[k] & forecast[,2]<=bins[k+1])
		nk[k] <- length(tempInds)
		oBarK[k] <- mean(outcome[tempInds])
	}

	# fix no value bins
	oBarK[is.na(oBarK)] <- 0

	# calculate the traditional three components
	oBar <- 1/n*sum(nk*oBarK)

	rel <- 1/n * sum(nk*(centers - oBarK)^2)
	res <- 1/n * sum(nk*(oBarK-oBar)^2)
	unc <- oBar * (1-oBar)

	brier <- rel - res + unc

	outVec <- c(brier,rel,res,unc)
	names(outVec) <- c('Brier', 'REL', 'RES', 'UNC')


	return(outVec)
}


# calculate the necessary parts for a reliability diagram
reliabilityCalc <- function(fcast,obs,binSize=0.05){
	
	bins <- seq(0-(binSize)/2,1+(binSize)/2,by=binSize)
	centers <- bins[-1]-(binSize)/2

	nk       <- matrix(NA, length(centers), 3)
	fcastBar <- matrix(NA, length(centers), 3)
	obsBar   <- matrix(NA, length(centers), 3)

	for(terc in 1:3){
		# figure out which obs were of this tercile
		obsInds <- which(obs==(terc-2))

		# loop through the rest of the bins
		for(i in 1:length(centers)){
			# figure out which forecasts gave this probabiliy
			tempInds <- which(fcast[,terc]>bins[i] & fcast[,terc]<=bins[i+1])
			
			# determine how many forecasts gave this probability (for sharpness)
			nk[i,terc] <- length(tempInds)

			# calculate the average forecast given within this probability
			fcastBar[i,terc] <- mean(fcast[tempInds,terc])

			# calculate the occurence of this terc WHEN the outcome actually occurs
			obsBar[i,terc] <- length(intersect(tempInds,obsInds))/length(tempInds)
		}

	}
	# return as one matrix with 10 columns:
	# 1 centers
	# 2:4 fcastBar
	# 5:7 obsBar
	# 8:10 nk

	return(cbind(centers,fcastBar,obsBar,nk))
}



# function to plot a reliability diagram
reliabilityDiagram <- function(relabilityMat,
		yr=NULL,textSize=1.5,scaleFactor=0,...){

	# upack the reliability object
	nk <- relabilityMat[,8:10]
	fcastBar <- relabilityMat[,2:4]
	obsBar <- relabilityMat[,5:7]

	# plotting params
	pal <- c('blue','green3','red')
	transPal <- adjustcolor(pal,alpha=0.3)
	transPal2 <- adjustcolor(pal,alpha=0.5)

	# rescale the bars to have a better axis
	barArrange <- nk/10^scaleFactor

	if(is.null(yr)){
		yr <- c(0,max(barArrange))*1.25
	}

	# plot the barplot of counts
	par(mar=c(5,5,5,6) + 0.1)
	barplot(t(barArrange),beside=T,axes=FALSE,ylim=yr,col=transPal)
	axis(4,cex.axis=textSize,cex.lab=textSize)
	
	par(new=TRUE) 

	# plot the reliability curves
	plot(x=NULL,xlim=c(0,1),ylim=c(0,1),
		xlab='Forecast Probability',ylab='Empirical Probability',
		cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, ...)
	

	mtext(sprintf('# of Forecasts Issued (x 10^%d)',scaleFactor), side=4, line=3,cex=1.25)

	grid()
	for(i in 1:3){
		xVec <- fcastBar[,i]
		yVec <- obsBar[,i]

		goodInds <- which(!is.na(xVec))

		points(xVec[goodInds], yVec[goodInds],type='b',col=pal[i],lwd=1.5)

		abline(lm(  yVec[goodInds] ~ xVec[goodInds],weights = nk[goodInds,i]),lty=2,lwd=2,col=transPal2[i])

	}
	abline(0,1)


	legend('topleft',
		c('La Niña', 'Neutral', 'El Niño'),
		lwd=2, col=pal,bg='white', cex = textSize, bty='n')

}


# a function that converts the indices in a testingInds vec to the indices in
# forecast/obs array points in a simple way
init2forecastInds <- function(indVec, testingInds){
	outInds <- rep(NA, length(indVec))

	for(i in 1:length(indVec)){
		outInds[i] <- which(testingInds == indVec[i])
	}

	return(outInds)
}


criticalCor <- function(n, alpha = .05, twoSided=FALSE) {
  if(twoSided) alpha <- alpha/2

  df <- n - 2
  critical.t <- qt( alpha, df, lower.tail = F )
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(critical.r)
}

# A helper function for plotting
fillInVec <- function(y){
	outVec <- rep(NA, length(xVec))
	for(b in 1:nBins){
		# get the leads in the give bin
		leadInds <- (leadBins[b]):(leadBins[b+1])
		
		outVec[leadInds] <- y[b]
	}
	return(outVec)			
}