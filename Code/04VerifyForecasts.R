
if(restart){
	load(sprintf('%s/step03_KeyResults.Rda',wdir))
}

###############################################################################
# Verify the forecast aggregating over all forecast conditions
# Probabilistic Verification (ROC and RPSS)
###############################################################################

# First, group forecasts with the same target season and lead

# initialize a nested list (thanks stack overflow!)
rocList    <- lapply(initHistList<-vector(mode = 'list',12),function(x) x<-vector(mode='list',forecastLength+1))
rpssList   <- lapply(initHistList<-vector(mode = 'list',12),function(x) x<-vector(mode='list',forecastLength+1))
rpssDecompList <- lapply(initHistList<-vector(mode = 'list',12),function(x) x<-vector(mode='list',forecastLength+1))


rocScoreArray <- array(NA, dim=c(12, forecastLength + 1, 3))
pValArray     <- array(NA, dim=c(12, forecastLength + 1, 3))

rpssArray     <- array(NA, dim=c(12, forecastLength + 1))
reliabilityArray    <- array(NA, dim=c(12, forecastLength + 1))
resolutionArray    <- array(NA, dim=c(12, forecastLength + 1))

reliabilityDiagArray <- array(NA, dim=c(21, 10, 12, forecastLength + 1))

# loop over target season and lead
for(i in 1:12){
	for(j in 1:(forecastLength + 1)){
		# This now should be getting all array inds where month = i 
		# and lead = j
		subInds <- which(forecastTimeArray[,j,2] == i)

		# Get the forecast mat and the observed vec
		forecastMat <- nino34ProbForecastArray[subInds,j,]
		obsEventVec <- eventObsArray[subInds,j]

		# ROC Score
		rocTemp  <- rocCalculationSingle(forecastMat, obsEventVec, probStep = 0.05)
		rocScoreArray[i,j,] <- rocTemp$rocScore
		pValArray[i,j,]     <- rocTemp$pval

		# RPSS
		rpssTemp <- rpssCalculationSingle(forecastMat, obsEventVec)
		rpssArray[i,j] <- rpssTemp$rpss

		# RPSS decomposition
		decompTemp <- rpssDecomp(forecastMat, obsEventVec)
		reliabilityArray[i,j] <- decompTemp$decomp$rel
		resolutionArray[i,j]  <- decompTemp$decomp$res

		# relaibility diagram stuff!
		reliabilityDiagArray[,,i,j] <- reliabilityCalc(forecastMat, obsEventVec)
	}
}

###############################################################################
# Deterministic Verification (AC and MSESS)
###############################################################################
climatologyForecast <- rep(NA,12)

for(i in 1:12){
	inds <- which(forecastTimeArray[,1,2] == i)
	climatologyForecast[i] <- mean(trueNino34[inds],na.rm=T)
}

acMat <- matrix(NA, forecastLength+1, 12)
acSigMat <- matrix(NA, forecastLength+1, 12)
sigMat <- matrix(NA, forecastLength+1, 12)
msessMat <- matrix(NA, forecastLength+1, 12)
ampBiasMat <- matrix(NA, forecastLength+1, 12)

for(i in 1:(forecastLength+1)){
	for(j in 1:12){
		subInds <- which(forecastTimeArray[,i,2] == j)

		# pull the nino34 series from the matrices
		obsSeries      <- nino34ObsArray[subInds,i]
		forecastSeries <- apply(nino34ForecastArray[subInds,i,],1,mean,na.rm=T)

		goodInds <- which(!is.na(obsSeries) & !is.na(forecastSeries))
		obsSeries      <- obsSeries[goodInds]
		forecastSeries <- forecastSeries[goodInds]

		# calculate the error (epsilon following DelSole+Tippett 2014)
		errorForecast    <- forecastSeries - obsSeries
		errorClimatology <- climatologyForecast[j] - obsSeries
		
		# calculate the MSEs
		forecastMSE <- mean((errorForecast)^2)
		climMSE     <- mean((errorClimatology)^2)

		# significance test
		sigMat[i,j] <- wilcox.test(errorForecast^2, errorClimatology^2,alternative='less')$p.value < 0.05
		# calculate the skill score
		msessMat[i,j] <- 1 - forecastMSE/climMSE


		# calculate anomaly correlation as well
		acMat[i,j] <- cor(forecastSeries,obsSeries)
		criticalValue <- criticalCor(length(forecastSeries))
		acSigMat[i,j] <- acMat[i,j] > criticalValue

		# calculate the amplitude bias
		ampBiasMat[i,j] <- (cor(forecastSeries,obsSeries) - sd(forecastSeries)/sd(obsSeries))^2
	}
}

msessList <- list(acSigMat=acSigMat, acMat=acMat,
	msessMat=msessMat, sigMat=sigMat, ampBiasMat=ampBiasMat)

save(rocScoreArray, pValArray,msessList,rpssArray,reliabilityArray,resolutionArray,
	file = sprintf('%s/step04_AllForecastVerification.Rda', wdir))


