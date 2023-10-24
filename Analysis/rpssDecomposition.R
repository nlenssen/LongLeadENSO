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
	# 'CM2.1-NMME', # rerunning overnight as it is huge (4000 years)
	'CM2.5-NMME')
nModels <- length(modelNameVec)

baseplotdir <- '~/Desktop/scratchFigures/ensoPaper'


# Other relevant vars from the namelist
forecastLength <- 5*12-1


# get the initialization stuff from /Scripts/runMultipleModels_Cross.R
# initModelName <- 'CESM1-1-CAM5-CMIP5'
initModelName <- 'Observations'
initExperimentName <- 'tos_zos_NEW'
trialFields <- c('tos', 'zos')
restart <- TRUE


plotdir <- sprintf('%s/crossModelFigures/%s',baseplotdir, initModelName)
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)


###############################################################################
# Loop over each model GIVEN the initialization model
###############################################################################
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
	# Look at the LN predictiability based on the EN magnitude at time of
	# initialization (3 quantiles to start)
	###############################################################################
	quantileBounds <- c(0, 0.75, 0.95, 1)
	quantileBounds <- c(0, 0.05, 0.25, 0.75, 0.95, 1)

	nquants <- length(quantileBounds) - 1

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


	rocListLN2  <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',forecastLength+1))
	rpssListLN2 <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',forecastLength+1))
	rpssDecompListLN2 <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',forecastLength+1))

	for(i in 1:nquants){
		for(j in 1:(forecastLength + 1)){
			
			djfInds <- which(forecastTimeArray[,j,2] == 1)
			
			subInds <- c()
			for(k in 1:length(djfInds)){
				tempNino34 <- trueNino34[which(initTimeMap[,3] == forecastTimeArray[djfInds[k],1,3])]

				if(is.na(tempNino34)) next()

				quants <- quantMat[forecastTimeArray[djfInds[k],1,2],]

				if(tempNino34 > quants[i] & tempNino34 <= quants[i+1]){
					subInds <- c(subInds, djfInds[k])
				}
			}

	 

			# Get the forecast mat and the observed vec
			forecastMat <- nino34ProbForecastArray[subInds,j,]
			obsVec      <- eventObsArray[subInds,j]

			rocListLN2[[i]][[j]]  <- rocCalculationSingle(forecastMat, obsVec)
			rpssListLN2[[i]][[j]] <- rpssCalculationSingle(forecastMat, obsVec)
			rpssDecompListLN2[[i]][[j]] <- rpssDecomp(forecastMat, obsVec)
		}

	}


	# extract things from the ROC list for easy plotting
	rocScoreArrayLN2  <- array(NA, dim=c(nquants, forecastLength + 1, 3))
	pValArrayLN2      <- array(NA, dim=c(nquants, forecastLength + 1, 3))

	rpssScoreArrayLN2  <- array(NA, dim=c(nquants, forecastLength + 1))

	rpssScoreCheckArrayLN2 <- array(NA, dim=c(nquants, forecastLength + 1))
	reliabilityArrayLN2    <- array(NA, dim=c(nquants, forecastLength + 1))
	resolutionARrayLN2     <- array(NA, dim=c(nquants, forecastLength + 1))

	for(i in 1:nquants){
		for(j in 1:(forecastLength + 1)){
			rocScoreArrayLN2[i,j,] <- rocListLN2[[i]][[j]]$rocScore
			pValArrayLN2[i,j,]     <- rocListLN2[[i]][[j]]$pval

			# get the rpss score
			rpssScoreArrayLN2[i,j]     <- rpssListLN2[[i]][[j]]$rpss

			# unpack the decomposition
			rpssScoreCheckArrayLN2[i,j] <- rpssDecompListLN2[[i]][[j]]$skill$rpss
			reliabilityArrayLN2[i,j]    <- rpssDecompListLN2[[i]][[j]]$decomp$rel
			resolutionARrayLN2[i,j]     <- rpssDecompListLN2[[i]][[j]]$decomp$res

		}
	}


	###############################################################################
	# RPSS Plotting (5 quants)
	###############################################################################
	pal <- c(rev(brewer.pal(9,'YlGnBu')[c(6,8)]), 'grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

	maxInd <- 36

	# only plotting la nina here
	xVec <- 0:maxInd
	subInds <- which(0:forecastLength < maxInd+1)

	# bar plot of rpss
	inds <- c(11,13,23,25)
	plotMat <- rpssScoreArrayLN2[,inds]
	colnames(plotMat) <- c('10', '12', '22','24')

	pdf(sprintf('%s/crossModel_%s_Rpss_ENStrat_Bar.pdf',plotdir, modelName),11,5)

	barplot(plotMat,col = pal, beside = TRUE, ylim=c(0, 1.0),
		xlab='lead (months)', ylab='RPSS',
		main='DJF Target Forecast')
	legend('topleft', rev(c('Strong LN', 'Weak LN', 'Neutral', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=15, bty='n')
	dev.off()


	# line plot of rpss
	pdf(sprintf('%s/crossModel_%s_Rpss_ENStrat_Combo.pdf',plotdir,modelName),11,10)
	par(mfrow=c(2,2))
	plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main='DJF Target Forecast')

	abline(h=0)
	abline(v=c(inds-1),lty=2)

	points(xVec, rpssArray[1,subInds], col='grey40', type='l', lwd=1.5, lty=3)

	for(i in 1:nquants){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, rpssScoreArrayLN2[i,subInds], col=pal[i], type='l', lwd=2)

	}
	legend('topright', rev(c('Strong LN', 'Weak LN', 'Neutral', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=2, bty='n')

	plot(xVec, rpssScoreArrayLN2[5,subInds]- rpssScoreArrayLN2[3,subInds],
		xlab='Lead (months)', ylab='Strong EN - No EN',type='l', lwd=2)
	abline(h=0)
	abline(v=c(inds-1),lty=2)

	plot(NULL, xlab='Lead (months)', ylab='RPS Decomposition',
			xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.55,0.3),
			main='DJF Target Forecast RPS Decomposition')
	abline(h=0)
	abline(v=c(inds-1),lty=2)

	for(i in 1:nquants){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, -reliabilityArrayLN2[i,subInds], col=pal[i], type='l', lwd=2, lty=3)
		points(xVec, resolutionARrayLN2[i,subInds], col=pal[i], type='l', lwd=2, lty=1)
	}
	legend('bottomleft', c('Resolution', 'Reliability'), col='black', lwd=2, lty=c(1,3), bty='n')

	barplot(plotMat,col = pal, beside = TRUE, ylim=c(-0.2, 1.0),
		xlab='lead (months)', ylab='RPSS',
		main='DJF Target Forecast')
	legend('topleft', rev(c('Strong LN', 'Weak LN', 'Neutral', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=15, bty='n')

	dev.off()
}