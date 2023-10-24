runAnalysis <- FALSE

source('Code/Functions.R')
# stuff
modelNameVec <- c('CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5', 
	'CESM2', 
	# 'GISS-E3', # too little data
	'GISS-E21G', 
	# 'GISS-E21H', # too little data
	'CESM1-NMME', 
	'CCSM4-NMME',
 	'CM2.1-NMME', # rerunning overnight as it is huge (4000 years)
	'CM2.5-NMME')
nModels <- length(modelNameVec)

baseplotdir <- '~/Desktop/scratchFigures/ensoPaper'


# Other relevant vars from the namelist
forecastLength <- 5*12-1


# get the initialization stuff from /Scripts/runMultipleModels_Cross.R


initExperimentName <- 'tos_zos_NEW'
trialFields <- c('tos', 'zos')
restart <- TRUE

# the different initializations to run over
initModelNameVec <- c('CESM1-1-CAM5-CMIP5', 'GISS-E21G', 'Observations')

###############################################################################
# Loop over the different initializations
###############################################################################
if(runAnalysis){

for(it in 1:length(initModelNameVec)){

initModelName <- initModelNameVec[it]



plotdir <- sprintf('%s/crossModelFigures_Binned_TargetBinned/%s',baseplotdir, initModelName)
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

ofdir <- sprintf('%s/plotData',plotdir)
if(!dir.exists(ofdir)) dir.create(ofdir, recursive = TRUE)


# lead bins of interest
leadBins <- c(1,2,7,13,19,25,31,37)
nBins <- length(leadBins) - 1

# quantile bounds
quantileBounds <- c(0, 0.75, 0.95, 1)
# quantileBounds <- c(0, 0.05, 0.25, 0.75, 0.95, 1)

nquants <- length(quantileBounds) - 1


###############################################################################
# Loop over each model GIVEN the initialization model
###############################################################################

masterRpssArray <- array(NA, dim=c(nquants, nBins, 3, nModels))

for(mod in 1:nModels){


	modelName <- modelNameVec[mod]

	# source the namelist to properly link/set all objects
	source('Namelists/baseNamelist.Rnl')


	# load the key data from the experiment
	load(sprintf('%s/step01_KeyInfo.Rda',wdir))
	load(sprintf('%s/step02_KeyResults.Rda',wdir))
	load(sprintf('%s/step03_KeyResults.Rda',wdir))
	load(sprintf('%s/step04_AllForecastVerification.Rda', wdir))


	###############################################################################
	# Look at the  predictiability based on the Nino34 magnitude at time of
	# initialization 
	###############################################################################

	# do some prep work to the event and nino vectors (true forecasts)
	indsDJF <- which(initTimeMap[,2]==1)
	nino34DJF   <- trueNino34[indsDJF]
	eventVecDJF <- trueEventVec[indsDJF]


	# calculate the quantiles over the testing period
	quantMat <- matrix(NA, 12, length(quantileBounds))
	for(i in 1:12){
		tempInds <- testingInds[which(forecastTimeArray[,1,2] == i)]
		quantMat[i,] <- quantile(trueNino34[tempInds], quantileBounds, na.rm=T)	
	}


	rpssScoreArrayLN2  <- array(NA, dim=c(nquants, nBins,3))
	reliabilityArrayLN2    <- array(NA, dim=c(nquants, nBins,3))
	resolutionArrayLN2     <- array(NA, dim=c(nquants, nBins,3))

	# loop over all quantiles and all lead bins
	for(i in 1:nquants){
		for(targ in -1:1){
			for(b in 1:nBins){

				# get the leads in the give bin
				leadInds <- (leadBins[b]):(leadBins[b+1]-1)

				# initialize the forecast/obs objects to collect over the bins
				forecastMat <- matrix(NA, nrow=0, ncol=3)
				obsVec      <- c()

				# loop over all leads in the block to collect forecast/obs pairs
				for(j in leadInds){
					
					djfInds <- which(forecastTimeArray[,j,2] == 1)
					
					# Find all of the inds for the given quant/lead 
					# AND given post state at that lead
					subInds <- c()
					for(k in 1:length(djfInds)){
						tempNino34 <- trueNino34[which(initTimeMap[,3] == forecastTimeArray[djfInds[k],1,3])]

						targetEvent <- trueEventVec[which(initTimeMap[,3] == forecastTimeArray[djfInds[k],j,3])]
						
						if(is.na(tempNino34)) next()
						if(is.na(targetEvent)) next()

						quants <- quantMat[forecastTimeArray[djfInds[k],1,2],]

						# check to see if it falls within the specific quant
						if(tempNino34 > quants[i] & tempNino34 <= quants[i+1]){
							# check to see the state at initialization
							if(targetEvent == targ){
								subInds <- c(subInds, djfInds[k])
							}
						}
					}

					if(is.null(subInds)) next()
					# Append the forecast mat and the observed vec

					forecastMat <- rbind(forecastMat, nino34ProbForecastArray[subInds,j,])
					obsVec      <- c(obsVec, eventObsArray[subInds,j])
				}

				if(nrow(forecastMat) == 0) next()
				
				tempRpss   <- rpssCalculationSingle(forecastMat, obsVec)
				tempDecomp <- rpssDecomp(forecastMat, obsVec)

				rpssScoreArrayLN2[i,b,targ+2]   <- tempRpss$rpss
				reliabilityArrayLN2[i,b,targ+2] <- tempDecomp$decomp$rel
				resolutionArrayLN2[i,b,targ+2]  <- tempDecomp$decomp$res
			}
		}		
	}

	# keep track of all 
	masterRpssArray[,,,mod] <- rpssScoreArrayLN2
	###############################################################################
	# RPSS Plotting (3 quants)
	###############################################################################
	
	pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

	# pal <- c(rev(brewer.pal(9,'YlGnBu')[c(6,8)]), 'grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

	# target stratification colors for second plot attempt
	pal2 <- c('blue', 'green', 'red')

	xVec <- 0:(max(leadBins)-1)
	maxInd <- max(xVec)

	elNameVec <- c('La Niña', 'Neutral', 'El Niño')
	initStratNameVec <- c('No EN', 'Weak EN', 'Strong EN')

	# line plot of rpss (binned)
	pdf(sprintf('%s/crossModel_%s_Rpss_DoubleDecomp01.pdf',plotdir,modelName),16,5)
	
	par(mfrow=c(1,3))
	
	for(k in 1:3){
		plot(NULL, xlab='Lead (months)', ylab='RPSS',
				xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-1, 1.0),
				main=sprintf('DJF %s Forecast (%s, %s Init.)',elNameVec[k], modelNameVec[mod], initModelName))

		abline(h=0)
		abline(v=c(12,24),lty=3)

		for(i in 1:nquants){
			# Plot DJF ROC Scores as a function of lead
			points(xVec, fillInVec(rpssScoreArrayLN2[i,,k]), col=pal[i], type='l', lwd=2)

		}
		legend('topright', rev(initStratNameVec), col=rev(pal), lwd=2, bty='n')

		# plot(NULL, xlab='Lead (months)', ylab='RPS Decomposition',
		# 		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-1,0.3),
		# 		main='RPS Decomposition')
		# abline(h=0)

		# for(i in 1:nquants){
		# 	# Plot DJF ROC Scores as a function of lead
		# 	points(xVec, fillInVec(-reliabilityArrayLN2[i,,k]), col=pal[i], type='l', lwd=2, lty=3)
		# 	points(xVec, fillInVec(resolutionARrayLN2[i,,k]), col=pal[i], type='l', lwd=2, lty=1)
		# }
		# legend('bottomleft', c('Resolution', 'Reliability'), col='black', lwd=2, lty=c(1,3), bty='n')
	}

	dev.off()


	# Plot all three
	pdf(sprintf('%s/crossModel_%s_Rpss_DoubleDecomp02.pdf',plotdir,modelName),16,5)
	
	par(mfrow=c(1,3))
	for(i in rev(1:nquants)){
	
		plot(NULL, xlab='Lead (months)', ylab='RPSS',
				xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-1, 1.0),
				main=sprintf('%s Init (%s, %s Init.)',initStratNameVec[i], modelNameVec[mod], initModelName))

		abline(h=0)
		abline(v=c(12,24),lty=3)

		for(k in 1:3){
			# Plot DJF ROC Scores as a function of lead
			points(xVec, fillInVec(rpssScoreArrayLN2[i,,k]), col=pal2[k], type='l', lwd=2)

		}
		legend('topright', paste(elNameVec, 'Target'), col=pal2, lwd=2, bty='n')
	}

	dev.off()

	# Plot just EN/LN (as neutral has generally low RPSS due to poor reliability in particular)
	subELInds <- c(1,3)
	pdf(sprintf('%s/crossModel_%s_Rpss_DoubleDecomp03.pdf',plotdir,modelName),16,5)
	
	par(mfrow=c(1,3))
	for(i in rev(1:nquants)){
	
		plot(NULL, xlab='Lead (months)', ylab='RPSS',
				xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
				main=sprintf('%s Init (%s, %s Init.)',initStratNameVec[i], modelNameVec[mod], initModelName))

		abline(h=0)
		abline(v=c(12,24),lty=3)

		for(k in subELInds){
			# Plot DJF ROC Scores as a function of lead
			points(xVec, fillInVec(rpssScoreArrayLN2[i,,k]), col=pal2[k], type='l', lwd=2)

		}
		legend('topright', paste(elNameVec[subELInds], 'Target'), col=pal2[subELInds], lwd=2, bty='n')
	}

	dev.off()

	save(xVec, rpssScoreArrayLN2, file=sprintf('%s/RPSS_%s.Rda',ofdir, modelName))


}

save(masterRpssArray, file=sprintf('%s/RPSS_AllModels.Rda',ofdir))

}
}

###############################################################################
# Final Figure 2 Work
###############################################################################

# lead bins of interest
leadBins <- c(1,2,7,13,19,25,31,37)
nBins <- length(leadBins) - 1

# quantile bounds
quantileBounds <- c(0, 0.75, 0.95, 1)
# quantileBounds <- c(0, 0.05, 0.25, 0.75, 0.95, 1)

nquants <- length(quantileBounds) - 1


plotdirFinal <- '~/Desktop/scratchFigures/ensoPaper/finalFigures'

# dim=c(nquants, nBins,ElLa,modelName,initModelName)
plotArray <- array(NA, dim=c(nquants,nBins,3,2,2))

# Stuff to make paths work
modelNameVec <- c('GISS-E21G','CESM1-1-CAM5-CMIP5')
initModelNameVec <- c('GISS-E21G','CESM1-1-CAM5-CMIP5')

for(i in 1:2){
	for(j in 1:2){
		plotdir <- sprintf('%s/crossModelFigures_Binned_TargetBinned/%s',baseplotdir, initModelNameVec[j])
		ofdir <- sprintf('%s/plotData',plotdir)
		load(sprintf('%s/RPSS_%s.Rda',ofdir, modelNameVec[i]))


		plotArray[,,,i,j] <- rpssScoreArrayLN2
	}
}

pal2 <- c('blue', 'green', 'red')

elNameVec <- c('La Niña', 'Neutral', 'El Niño')
initStratNameVec <- c('No EN', 'Weak EN', 'Strong EN')

xVec <- 0:(max(leadBins)-1)
maxInd <- max(xVec)

# only look at LN/EN Targets (dim 3 of array)
subELInds <- c(1,3)


mod <- 1
it  <- 1
inds <- 3:1
for(i in 1:length(inds)){
	pdf(sprintf('%s/Figure_02_%s.pdf',plotdirFinal,letters[i]),7,7)

	plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main=sprintf('%s at Initialization (%s Target)', initStratNameVec[inds[i]], initModelNameVec[it]))

	abline(h=0)
	abline(v=c(12,24),lty=3)

	for(k in subELInds){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, fillInVec(plotArray[inds[i],,k,mod,it]), col=pal2[k], type='l', lwd=2)

		points(xVec, fillInVec(plotArray[inds[i],,k,mod,it+1]), col=pal2[k], type='l', lwd=2,lty=2)

	}

	legend('topright', c(paste(elNameVec[subELInds], 'Target (Perfect)'),
						 paste(elNameVec[subELInds], 'Target (Cross-Model)')),
				col=pal2[subELInds], lwd=2, lty=c(1,1,2,2), bty='n')

	legend(-3.34, 1.08, sprintf('(%s)',letters[i]) , bty='n')

	dev.off()
}



mod <- 2
it  <- 2
inds <- 3:1
for(i in 1:length(inds)){
	pdf(sprintf('%s/Figure_02_%s.pdf',plotdirFinal,letters[i+3]),7,7)

	plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main=sprintf('%s at Initialization (%s Target)', initStratNameVec[inds[i]], initModelNameVec[it]))

	abline(h=0)
	abline(v=c(12,24),lty=3)

	for(k in subELInds){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, fillInVec(plotArray[inds[i],,k,mod,it]), col=pal2[k], type='l', lwd=2)

		points(xVec, fillInVec(plotArray[inds[i],,k,mod,it-1]), col=pal2[k], type='l', lwd=2,lty=2)

	}

	legend('topright', c(paste(elNameVec[subELInds], 'Target (Perfect)'),
						 paste(elNameVec[subELInds], 'Target (Cross-Model)')),
				col=pal2[subELInds], lwd=2, lty=c(1,1,2,2), bty='n')

	legend(-3.34, 1.08, sprintf('(%s)',letters[i+3]) , bty='n')

	dev.off()
}



