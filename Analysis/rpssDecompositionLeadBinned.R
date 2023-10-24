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

plotdirFinal <- '~/Desktop/scratchFigures/ensoPaper/finalFigures'

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
for(it in 1:length(initModelNameVec)){

initModelName <- initModelNameVec[it]



plotdir <- sprintf('%s/crossModelFigures_Binned/%s',baseplotdir, initModelName)
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)


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
if(runAnalysis){
masterRpssArray <- array(NA, dim=c(nquants, nBins, nModels))

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

	# prep objects to save into
	rocScoreArrayLN2  <- array(NA, dim=c(nquants, nBins, 3))
	pValArrayLN2      <- array(NA, dim=c(nquants, nBins, 3))

	rpssScoreArrayLN2  <- array(NA, dim=c(nquants, nBins))

	rpssScoreCheckArrayLN2 <- array(NA, dim=c(nquants, nBins))
	reliabilityArrayLN2    <- array(NA, dim=c(nquants, nBins))
	resolutionARrayLN2     <- array(NA, dim=c(nquants, nBins))

	# loop over all quantiles and all lead bins
	for(i in 1:nquants){
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
				subInds <- c()
				for(k in 1:length(djfInds)){
					tempNino34 <- trueNino34[which(initTimeMap[,3] == forecastTimeArray[djfInds[k],1,3])]

					if(is.na(tempNino34)) next()

					quants <- quantMat[forecastTimeArray[djfInds[k],1,2],]

					if(tempNino34 > quants[i] & tempNino34 <= quants[i+1]){
						subInds <- c(subInds, djfInds[k])
					}
				}

				# Append the forecast mat and the observed vec
				forecastMat <- rbind(forecastMat, nino34ProbForecastArray[subInds,j,])
				obsVec      <- c(obsVec, eventObsArray[subInds,j])
			}
			
			# now calculate the skill scores for the given chunk of leads
			rocListTemp  <- rocCalculationSingle(forecastMat, obsVec)
			rpssListTemp <- rpssCalculationSingle(forecastMat, obsVec)
			rpssDecompListTemp <- rpssDecomp(forecastMat, obsVec)

			# save the ROC stuff
			rocScoreArrayLN2[i,b,] <- rocListTemp$rocScore
			pValArrayLN2[i,b,]     <- rocListTemp$pval

			# save the rpss score
			rpssScoreArrayLN2[i,b]     <- rpssListTemp$rpss

			# unpack and save the decomposition
			rpssScoreCheckArrayLN2[i,b] <- rpssDecompListTemp$skill$rpss
			reliabilityArrayLN2[i,b]    <- rpssDecompListTemp$decomp$rel
			resolutionARrayLN2[i,b]     <- rpssDecompListTemp$decomp$res


		}		
	}

	# save things to an full output array
	masterRpssArray[,,mod] <- rpssScoreArrayLN2

	###############################################################################
	# RPSS Plotting (5 quants)
	###############################################################################
	pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])
	# pal <- c(rev(brewer.pal(9,'YlGnBu')[c(6,8)]), 'grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

	xVec <- 0:(max(leadBins)-1)
	maxInd <- max(xVec)


	# line plot of rpss (binned)
	pdf(sprintf('%s/crossModel_%s_Rpss_ENStrat_Combo.pdf',plotdir,modelName),11,5)
	par(mfrow=c(1,2))
	plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main=sprintf('DJF Target Forecast (%s: Obs. Init.)',modelNameVec[mod]))

	abline(h=0)

	for(i in 1:nquants){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, fillInVec(rpssScoreArrayLN2[i,]), col=pal[i], type='l', lwd=2)

	}
	legend('topright', rev(c('No EN', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=2, bty='n')

	
	plot(NULL, xlab='Lead (months)', ylab='RPS Decomposition',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.3,0.3),
			main='RPS Decomposition')
	abline(h=0)

	for(i in 1:nquants){
		# Plot DJF ROC Scores as a function of lead
		points(xVec, fillInVec(-reliabilityArrayLN2[i,]), col=pal[i], type='l', lwd=2, lty=3)
		points(xVec, fillInVec(resolutionARrayLN2[i,]), col=pal[i], type='l', lwd=2, lty=1)
	}
	legend('bottomleft', c('Resolution', 'Reliability'), col='black', lwd=2, lty=c(1,3), bty='n')

	dev.off()
}

save(masterRpssArray, file=sprintf('Data/binnedAnalysis_Init_%s.Rda',initModelName))
}
load(sprintf('Data/binnedAnalysis_Init_%s.Rda',initModelName))
###############################################################################
# plot the difference between the quants over lead time for all 5 models at all leads
###############################################################################
pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])
# pal <- c(rev(brewer.pal(9,'YlGnBu')[c(6,8)]), 'grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

xVec <- 0:(max(leadBins)-1)
maxInd <- max(xVec)

pdf(sprintf('%s/crossModel_Differences_Final_01.pdf',plotdir),7,7)
# strong EN - no EN
plot(NULL, xlab='Lead (months)', ylab='Strong EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.4,1),
		main=sprintf('Extra Strong EN Skill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[3,,i]) - fillInVec(masterRpssArray[1,,i]), col=i, type='l',lwd=2)
}
abline(h=0)
abline(v=c(12,24),lty=3)

legend('topright', modelNameVec, col=1:nModels, lwd=2, bty='n', cex=0.75)
dev.off()

# weak EN - no EN
pdf(sprintf('%s/crossModel_Differences_Final_02.pdf',plotdir),7,7)
plot(NULL, xlab='Lead (months)', ylab='Weak EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.4,1),
		main=sprintf('Extra Weak EN SKill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[2,,i]) - fillInVec(masterRpssArray[1,,i]), col=i, type='l',lwd=2)
}
abline(h=0)
abline(v=c(12,24),lty=3)
legend('topright', modelNameVec, col=1:nModels, lwd=2, bty='n', cex=0.75)
dev.off()

}

###############################################################################
# Figure 03 Paper Plots
###############################################################################

# GISS First
initModelName <- 'GISS-E21G'
load(sprintf('Data/binnedAnalysis_Init_%s.Rda',initModelName))


pdf(sprintf('%s/Figure_03_b.pdf',plotdirFinal),7,7)
# strong EN - no EN
plot(NULL, xlab='Lead (months)', ylab='Strong EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=pm(0.6),
		main=sprintf('Extra Strong EN Skill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[3,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)

legend(-3.34, 0.68, '(b)', bty='n')
dev.off()

# weak EN - no EN
pdf(sprintf('%s/Figure_03_c.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='Weak EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=pm(0.6),
		main=sprintf('Extra Weak EN SKill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[2,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)


legend(-3.34, 0.68, '(c)', bty='n')

dev.off()


# CESM next

initModelName <- 'CESM1-1-CAM5-CMIP5'
load(sprintf('Data/binnedAnalysis_Init_%s.Rda',initModelName))


pdf(sprintf('%s/Figure_03_e.pdf',plotdirFinal),7,7)
# strong EN - no EN
plot(NULL, xlab='Lead (months)', ylab='Strong EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=pm(0.6),
		main=sprintf('Extra Strong EN Skill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[3,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)

legend(-3.34, 0.68, '(e)', bty='n')
dev.off()

# weak EN - no EN
pdf(sprintf('%s/Figure_03_f.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='Weak EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=pm(0.6),
		main=sprintf('Extra Weak EN SKill (%s Init.)',initModelName))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[2,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)


legend(-3.34, 0.68, '(f)', bty='n')

dev.off()


