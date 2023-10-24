if(restart) load(sprintf('%s/step01_KeyInfo.Rda',wdir))

###############################################################################
# Find `analogueSize` top analogues
###############################################################################

# the inds of the init/verif model to use (also probably all of them eventually)
testingInds <- 1:(nrow(initTimeMap) - 100)

if(initModelName == 'Observations') testingInds <- 1:696
nInit <- length(testingInds)

# the inds of the library model to use (probably use all of them)
libraryInds <- 1:nrow(timeMap)

##########
# helper function to pass to foreach (i is index of the testingInds and is of length nInit)
##########

createForecast <- function(i){
	require(ncdf4)

	##########
	# Pull needed data for the state at forecast initialization
	# NOTE: This could be from a different model than the library of analogues
	##########
	
	initialStateFields <- array(NA, dim=c(length(lonInds),length(latInds),length(trialFields)))

	for(j in 1:length(trialFields)){
		handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,initModelName,trialFields[j],initModelName))
		initialStateFields[,,j] <-  ncvar_get(handle,sprintf('%s_anom',trialFields[j]),
						start=c(min(lonInds),min(latInds),testingInds[i]),
						count=c(length(lonInds),length(latInds),1))*maximalMask[lonInds,latInds]
		nc_close(handle)
	}


	# build the initial state list to feed into the forecast function
	initialState <- list(lon=lon[lonInds],
						 lat=lat[latInds],
						 fields=initialStateFields,
						 fieldNames=trialFields,
						 timeMap=initTimeMap,
						 mask=maximalMask)


	# make sure to exclude the initialization if testing perfect model
	libIndsExclude <- NULL
	if(perfectModel){
		libIndsExclude <- testingInds[i]
	}

	##########
	# Run the function to find the best analogues and return output
	##########
	outList <- createAnalogueForecastMinimal(forecastInd = testingInds[i],
											 initialState = initialState,
											 libraryInds = libraryInds,
											 errorFunction = simpleMseError,
											 forecastFieldNames = forecastFieldNames,
											 libIndsExclude = libIndsExclude,
											 analogueSize = analogueSize)

	return(outList)
}

##########
# Run the forecast code (parallelized, chunked to make memory leak less of a problem)
##########
chunkSize <- 5e2
startVec <- seq(1,nInit,by=chunkSize)

outList <- vector(mode='list', length=nInit)

pb   <- txtProgressBar(1, length(startVec), style=3)
for(j in 1:length(startVec)){
	setTxtProgressBar(pb, j)

	loopInds <- startVec[j]:min(nInit, startVec[j]+(chunkSize-1))

	cl <- makeCluster(nCores_step02)
	registerDoParallel(cl)
	tempList <- foreach(i=loopInds) %dopar% createForecast(i)
	stopCluster(cl)
	gc()	

	for(k in 1:length(loopInds)){
		outList[[loopInds[k]]] <- tempList[k]
	}
	rm(tempList)
	gc()	
}




##########
# process the loop output and save to other objects (think about the best way
# to do this in the future, this is probably fine though)
##########

# stuff to save as we loop over all initializations
errorMat     <- matrix(NA,analogueSize, nInit)
ensembleTimePointsMat <- matrix(NA,analogueSize, nInit)

for(i in 1:nInit){
	ensembleTimePointsMat[,i] <- outList[[i]][[1]]$ensembleInds
	errorMat[,i]              <- outList[[i]][[1]]$errorVec
}


##########
# process the loop output and save to other objects (think about the best way
# to do this in the future, this is probably fine though)
##########
save(ensembleTimePointsMat, errorMat, testingInds, nInit, libraryInds, file=sprintf('%s/step02_KeyResults.Rda',wdir))
