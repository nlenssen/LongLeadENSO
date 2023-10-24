# keep false if not using new data as this will make it run in ~2 sec instead of 2 hours
runAnalysis <- FALSE

library(vioplot)

modelNameVec <- c(
	'CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5', # already ran
	'CESM2', 
	# 'GISS-E3', # too little data
	'GISS-E21G', 
	# 'GISS-E21H', # too little data
	'CESM1-NMME', 
	'CCSM4-NMME',
	'CM2.1-NMME', # rerunning overnight as it is huge (4000 years)
	'CM2.5-NMME'
	)

nModels <- length(modelNameVec)

# model names same for perfect
initModelNameVec <- modelNameVec
initExperimentName <- 'tos_zos_NEW'

# set other things
trialFields <- c('tos', 'zos')

# thing pasted from the namelist just to make the code work
forecastLength <- 5*12-1

plotdir <- '~/Desktop/scratchFigures/ensoPaper/finalFigures'
if(!dir.exists(plotdir)) dir.create(plotdir)

###############################################################################
# Load observational skill from step04
###############################################################################
obsRpssArray <- array(NA, dim=c(12, forecastLength + 1, nModels))

# Eventual loop over models
for(mod in 1:nModels){

	# set the two model names
	modelName <- modelNameVec[mod]
	initModelName <- 'Observations'

	# run the analyses
	source('Namelists/baseNamelist.Rnl')


	# load in needed data
	load(sprintf('%s/step04_AllForecastVerification.Rda',wdir))


	obsRpssArray[,,mod] <- rpssArray
}


###############################################################################
# Load observational skill from binned analysis
###############################################################################
load('Data/binnedAnalysis_Init_Observations.Rda')
obsBinnedRpssArray <- masterRpssArray

###############################################################################
# Stuff from runMultipleModel_Perfect
###############################################################################

# The parameter that is currently controlling how long of a chunk to verify
# from the full piControl run.
#
# NOTE: hard coded from obs right now
initSampleLength <- 109*12


# number of times to sample from the piControl
nSamples <- 200


###############################################################################
# Stuff from rpssDecompositionLeadBinned
###############################################################################

# lead bins of interest
leadBins <- c(1,2,7,13,19,25,31,37)
nBins <- length(leadBins) - 1

# quantile bounds
quantileBounds <- c(0, 0.75, 0.95, 1)
# quantileBounds <- c(0, 0.05, 0.25, 0.75, 0.95, 1)

nquants <- length(quantileBounds) - 1


###############################################################################
# Loop over models and set up the workspace
###############################################################################
if(runAnalysis){
rpssArraySampling        <- array(NA, dim=c(12, forecastLength + 1,nSamples, nModels))
reliabilityArraySampling <- array(NA, dim=c(12, forecastLength + 1,nSamples, nModels))
resolutionArraySampling  <- array(NA, dim=c(12, forecastLength + 1,nSamples, nModels))

# also track the number of each event in the sample to show how that relates
# to better skill
nEventsMat <- array(NA, dim=c(3, nSamples, nModels))

# save all the binned RPSS
masterBinnedRpssArray <- array(NA, dim=c(nquants, nBins, nSamples, nModels))


# Eventual loop over models
for(mod in 1:nModels){

	# set the two model names
	modelName <- modelNameVec[mod]
	initModelName <- 'CESM1-1-CAM5-CMIP5'

    print(sprintf('%02d/%02d Model: %s, Init: %s',mod, nModels, modelName, initModelName))

	# run the analyses
	source('Namelists/baseNamelist.Rnl')

	# load in needed data
	load(sprintf('%s/step01_KeyInfo.Rda',wdir))
	load(sprintf('%s/step02_KeyResults.Rda',wdir))
	load(sprintf('%s/step03_KeyResults.Rda',wdir))

	###############################################################################
	# Loop over the number of samples and perform various verifications
	###############################################################################

	# create a sample padding enough to make sure every forecast can be complete
	sampleStartInds <- sample(1:(nInit - forecastLength - 100), nSamples, replace=T)



	for(s in 1:nSamples){
		##########
		# 01 Verify the forecast aggregating over all forecast conditions
		#			(as in 04VerifyForecasts.R)
		##########
		tempSampleInds <- sampleStartInds[s]:(sampleStartInds[s] + initSampleLength - 1)

		# get the event counts over the record
		subInds0 <- which(initTimeMap[,2] == 1)
		subInds  <- intersect(subInds0, tempSampleInds)
		nEventsMat[,s,mod] <- table(trueEventVec[subInds])


		# loop over target season and lead
		for(i in 1:12){
			for(j in 1:(forecastLength + 1)){
				# This now should be getting all array inds where month = i 
				# and lead = j
				subInds0 <- which(forecastTimeArray[,j,2] == i)
				subInds  <- intersect(subInds0, tempSampleInds)

				# Get the forecast mat and the observed vec
				forecastMat <- nino34ProbForecastArray[subInds,j,]
				obsEventVec <- eventObsArray[subInds,j]
				
				rpssTempList   <- rpssCalculationSingle(forecastMat, obsEventVec)
				decompTempList <- rpssDecomp(forecastMat, obsEventVec)

				rpssArraySampling[i,j,s,mod] <- rpssTempList$rpss
				reliabilityArraySampling[i,j,s,mod] <- decompTempList$decomp$rel
				resolutionArraySampling[i,j,s,mod]  <- decompTempList$decomp$res
			}
		}



		##########
		# 02) Verify the forecast for the quantiles and time bins
		# 			(as in rpssDecompositionLeadBinned.R)
		##########

		# calculate the quantiles over the testing period
		quantMat <- matrix(NA, 12, length(quantileBounds))
		for(i in 1:12){
			tempInds <- testingInds[intersect(which(forecastTimeArray[,1,2] == i),tempSampleInds)]
			quantMat[i,] <- quantile(trueNino34[tempInds], quantileBounds, na.rm=T)	
		}

		# prep objects to save into
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
					
					djfInds <- intersect(which(forecastTimeArray[,j,2] == 1),tempSampleInds)
					
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
				rpssListTemp <- rpssCalculationSingle(forecastMat, obsVec)
				rpssDecompListTemp <- rpssDecomp(forecastMat, obsVec)

				# save the rpss score
				rpssScoreArrayLN2[i,b]     <- rpssListTemp$rpss

				# unpack and save the decomposition
				rpssScoreCheckArrayLN2[i,b] <- rpssDecompListTemp$skill$rpss
				reliabilityArrayLN2[i,b]    <- rpssDecompListTemp$decomp$rel
				resolutionARrayLN2[i,b]     <- rpssDecompListTemp$decomp$res


			}		
		}

		# save things to an full output array
		masterBinnedRpssArray[,,s,mod] <- rpssScoreArrayLN2

	}
}

save(rpssArraySampling, reliabilityArraySampling, resolutionArraySampling,
	masterBinnedRpssArray,
	file='Data/samplingUncAnalysisResults.R')
}
###############################################################################
# Some plotting
###############################################################################



load('Data/samplingUncAnalysisResults.R')
############
# Figure 4a
############

# subset the lead time
maxInd <- 36
xVec <- 0:maxInd
subInds <- which(0:forecastLength < maxInd+1)

# set plotting parameters
pal <- 1:nModels

pdf(sprintf('%s/Figure_04_a.pdf',plotdir), 7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Obs. Hindcast 1901-2009)')

abline(h=c(0))
abline(v=c(12,24), lty=3)

for(i in 1:nModels){
	points(xVec, obsRpssArray[1,subInds,i],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(18,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')
legend(-3.34, 1.08, '(a)', bty='n')

dev.off()


############
# Figure 4b
############

# DJF focus for now
i <- 1

# subset the lead time
maxInd <- 36
xVec <- 0:maxInd
subInds <- which(0:forecastLength < maxInd+1)

pal <- c('black', 'blue')

# remove the perfect model skill from this plot
modelLoopVec <- 1:nModels[-c(5,8)]


# get the full sampling CI
samplingCI_Cross <- apply(rpssArraySampling[i,subInds,,modelLoopVec],1, quantile, probs=c(0.025, 0.5, 0.975))

# get the spread of the obs
spreadSkill_Obs <- apply(obsRpssArray[i,subInds,],1,quantile, probs=c(0, 0.5, 1))


pdf(sprintf('%s/Figure_04_b.pdf',plotdir), 7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main=sprintf('Observational Hindcast Sampling Uncertainty'))

abline(h=c(0))
abline(v=c(12,24), lty=3)

# plot the sampling CI estimates from CESM1.1 Cross
points(xVec, samplingCI_Cross[2,], type='l', lwd=2, col=pal[1])
polygon(c(xVec,rev(xVec)), c(samplingCI_Cross[1,], rev(samplingCI_Cross[3,])), col=adjustcolor(pal[1],alpha=0.25), border=NA)

# plot the
points(xVec, spreadSkill_Obs[2,], type='l', lwd=2, col=pal[2])
polygon(c(xVec,rev(xVec)), c(spreadSkill_Obs[1,], rev(spreadSkill_Obs[3,])), col=adjustcolor(pal[2],alpha=0.25), border=NA)

legend(9,0.7, c('Cross-Model Sampling Unc. Skill (CESM1.1)', 'Observational Hindcasts (1901-2009)'), col=pal, lwd=c(2), bty='n')
legend(-3.34, 1.08, '(b)', bty='n')
dev.off()


############
# Figure 4c attempt02 (VIOLIN)
############
shiftFactor <- 1.5
squishFactor <- 1

pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])
pal2 <- adjustcolor(pal, alpha=0.2)

xVec <- c(0, ((leadBins[-1] - diff(leadBins)/2)[-1])-1) - shiftFactor


pdf(sprintf('%s/Figure_04_c.pdf',plotdir), 7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(-1,maxInd), ylim=c(-0.2, 1.0),
		main=sprintf('Observational Hindcast State Dependent Skill'))

# Add key lines here
abline(h=c(0,1))
abline(v=c(12,24), lty=3)

# plot the sampling dist
for(i in 1:3){
	for(b in 1:nBins){
		vioplot(masterBinnedRpssArray[i,b,,modelLoopVec], at=xVec[b] + (i-1)*squishFactor, add=T,
			wex=2, col=pal2[i], colMed2=pal2[i], colMed='black', pchMed=21)
	}
}

# plot the values calculated in obs
for(i in 1:3){
	for(b in 1:nBins){
		boxplot(obsBinnedRpssArray[i,b,modelLoopVec], at=xVec[b] + (i-1)*squishFactor, add=T,
			col=pal[i], outline=FALSE, width=0.5)
	}
}


legend(26,0.97, rev(c('No EN', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=10, bty='n')
legend(-4, 1.08, '(c)', bty='n')

dev.off()


############
# Figure 4c attempt01 (NEED TO MAKE VIOLIN)
############

# # subset the lead time
# xVec <- 0:(max(leadBins)-1)

# # use the same pal as previous figures

# # remove the perfect model skill from this plot
# modelLoopVec <- 1:nModels[-c(5,8)]

# # get the full sampling CI
# binnedCI_Cross <- apply(masterBinnedRpssArray[,,,modelLoopVec],c(1,2), quantile, probs=c(0.025, 0.5, 0.975))

# # get the spread of the obs
# binnedSpread_Obs <- apply(obsBinnedRpssArray[,,modelLoopVec],c(1,2),quantile, probs=c(0.025, 0.5, 0.975))


# pdf(sprintf('%s/Figure_04_c_Attempt01.pdf',plotdir), 7,7)
# plot(NULL, xlab='Lead (months)', ylab='RPSS',
# 		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
# 		main=sprintf('Initial State Dependent Skill (Obs. 1901-2009)'))

# abline(h=0)

# # plot the sampling CI estimates from CESM1.1 Cross
# for(i in 1:3){
# 	points(xVec, fillInVec(binnedCI_Cross[2,i,]), type='l', lwd=2, col=pal[i])
# 	polygon(c(xVec,rev(xVec)), c(fillInVec(binnedCI_Cross[1,i,]), rev(fillInVec(binnedCI_Cross[3,i,]))), col=adjustcolor(pal[i],alpha=0.25), border=NA)
# }

# # plot the obs estimates as well (dashed?)
# for(i in 1:3){
# 	points(xVec, fillInVec(binnedSpread_Obs[2,i,]), type='l', lty=2, lwd=2, col=pal[i])
# 	polygon(c(xVec,rev(xVec)), c(fillInVec(binnedSpread_Obs[1,i,]), rev(fillInVec(binnedSpread_Obs[3,i,]))), lty=2, col=adjustcolor(pal[i],alpha=0.25), border=NA)
# }


# legend(6,0.8, c('Cross-Model Sampling Unc. Skill (CESM1.1)', 'Observational Hindcasts (1901-2009)'), col=pal, lwd=c(2), bty='n')

# dev.off()

