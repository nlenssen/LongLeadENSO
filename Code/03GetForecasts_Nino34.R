if(restart){
	load(sprintf('%s/step01_KeyInfo.Rda',wdir))
	load(sprintf('%s/step02_KeyResults.Rda',wdir))
}

###############################################################################
# Pull the relevant info from the time steps according to the top K indices as 
# set by `ensembleSize`

# NOTE: This is where you now pull fields, etc to try to separate finding the
# closest matching states from pulling the variabile(s) of interest.
###############################################################################


# get the raw nino34 index [ensemble member, initialization, lead]
nino34ObsArray <- array(NA, dim=c(nInit, forecastLength+1))
nino34ForecastArray <- array(NA, dim=c(nInit, forecastLength+1, ensembleSize))

# get the event forecast [ensemble member, initialization, lead]
eventObsArray      <- array(NA, dim=c(nInit, forecastLength+1))
eventForecastArray <- array(NA, dim=c(nInit, forecastLength+1, ensembleSize))


# loop through the initializations and get the forecast/obs pair objects
for(i in 1:nInit){

	# get the time points in terms of the full piControl
	testInd <- testingInds[i]
	forecastInds <- ensembleTimePointsMat[1:ensembleSize, i]

	# get the TRUE nino34 progression for each of these
	nino34ObsArray[i,] <- trueNino34[testInd:(testInd+forecastLength)]
	eventObsArray[i,]  <- trueEventVec[testInd:(testInd+forecastLength)]
	

	for(j in 1:ensembleSize){
		# get the location of the forecast 
		tempInd <- forecastInds[j]

		nino34ForecastArray[i,,j] <- nino34[tempInd:(tempInd+forecastLength)]
		eventForecastArray[i,,j]  <- eventVec[tempInd:(tempInd+forecastLength)]
	}
}

# make probablistic forecasts [LN, N, EN]
nino34ProbForecastArray <- array(NA, dim=c(nInit, forecastLength+1, 3))

for(i in 1:nInit){
	forecastMatTemp <- eventForecastArray[i,,]

	for(j in 1:(forecastLength + 1)){
		fcastVec <- rep(NA, 3)

		for(k in -1:1){
			fcastVec[k+2] <- sum(forecastMatTemp[j,]==k)/ensembleSize
		}
		
		nino34ProbForecastArray[i,j,] <- fcastVec
	}
}

###############################################################################
# Create some important metadata objects to make it easier to verify the
# forecasts
###############################################################################
forecastTimeArray <- array(NA, dim=c(dim(nino34ObsArray), 3))

for(i in 1:nInit){
	initInd <- testingInds[i]
	forecastTimeArray[i,,] <- initTimeMap[initInd:(initInd+forecastLength),]
}


save(forecastTimeArray, nino34ForecastArray, nino34ProbForecastArray, eventObsArray,nino34ObsArray,
	file = sprintf('%s/step03_KeyResults.Rda', wdir))