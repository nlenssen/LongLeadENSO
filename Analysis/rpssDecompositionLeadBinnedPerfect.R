runAnalysis <- FALSE

source('Code/Functions.R')

# stuff
modelNameVec <- c('CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5', 
	'CESM2', 
	'GISS-E21G', 
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

# lead bins of interest
leadBins <- c(1,2,7,13,19,25,31,37)
nBins <- length(leadBins) - 1

# quantile bounds
quantileBounds <- c(0, 0.75, 0.95, 1)
# quantileBounds <- c(0, 0.05, 0.25, 0.75, 0.95, 1)

nquants <- length(quantileBounds) - 1

plotdir <- sprintf('%s/crossModelFigures_Binned_Perfect',baseplotdir)
if(!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

# Stuff for plotting
xVec <- 0:(max(leadBins)-1)
maxInd <- max(xVec)


###############################################################################
# Loop over the different initializations
###############################################################################


if(runAnalysis){

masterRpssArray <- array(NA, dim=c(nquants, nBins, nModels))

for(mod in 1:nModels){









###############################################################################
# Loop over each model GIVEN the initialization model
###############################################################################


	initModelName <- modelNameVec[mod]
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


	rocListLN2  <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',nBins))
	rpssListLN2 <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',nBins))
	rpssDecompListLN2 <- lapply(initHistList<-vector(mode = 'list',nquants),function(x) x<-vector(mode='list',nBins))

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
			rocListLN2[[i]][[b]]  <- rocCalculationSingle(forecastMat, obsVec)
			rpssListLN2[[i]][[b]] <- rpssCalculationSingle(forecastMat, obsVec)
			rpssDecompListLN2[[i]][[b]] <- rpssDecomp(forecastMat, obsVec)
		}		
	}



	# extract things from the ROC list for easy plotting
	rocScoreArrayLN2  <- array(NA, dim=c(nquants, nBins, 3))
	pValArrayLN2      <- array(NA, dim=c(nquants, nBins, 3))

	rpssScoreArrayLN2  <- array(NA, dim=c(nquants, nBins))

	rpssScoreCheckArrayLN2 <- array(NA, dim=c(nquants, nBins))
	reliabilityArrayLN2    <- array(NA, dim=c(nquants, nBins))
	resolutionARrayLN2     <- array(NA, dim=c(nquants, nBins))

	for(i in 1:nquants){
		for(j in 1:nBins){
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

	masterRpssArray[,,mod] <- rpssScoreArrayLN2
	###############################################################################
	# RPSS Plotting (5 quants)
	###############################################################################
	pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])
	# pal <- c(rev(brewer.pal(9,'YlGnBu')[c(6,8)]), 'grey40',brewer.pal(9,'YlOrRd')[c(5,8)])

	# line plot of rpss (binned)
	pdf(sprintf('%s/crossModel_%s_Rpss_ENStrat_Combo.pdf',plotdir,modelName),11,5)
	par(mfrow=c(1,2))
	plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main=sprintf('DJF Target Forecast (%s Perfect)',modelNameVec[mod]))

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


save(masterRpssArray, file='Data/binnedAnalysis_Init_Perfect.Rda')



}

load('Data/binnedAnalysis_Init_Perfect.Rda')

###############################################################################
# plot the difference between the quants over lead time for all 5 models at all leads
###############################################################################

pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])



pdf(sprintf('%s/crossModel_Differences_Final.pdf',plotdir),14,6)
par(mfrow=c(1,2))
# strong EN - no EN
plot(NULL, xlab='Lead (months)', ylab='Strong EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.4,1),
		main=sprintf('Extra Strong EN Skill (Perfect)'))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[3,,i]) - fillInVec(masterRpssArray[1,,i]), col=i, type='l',lwd=2)
}
abline(h=0)
abline(v=c(12,24),lty=3)

legend('topright', modelNameVec, col=1:nModels, lwd=2, bty='n', cex=0.75)

# weak EN - no EN
plot(NULL, xlab='Lead (months)', ylab='Weak EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.4,1),
		main=sprintf('Extra Weak EN SKill (Perfect)'))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[2,,i]) - fillInVec(masterRpssArray[1,,i]), col=i, type='l',lwd=2)
}
abline(h=0)
abline(v=c(12,24),lty=3)
legend('topright', modelNameVec, col=1:nModels, lwd=2, bty='n', cex=0.75)
dev.off()




###############################################################################
# Figure 01b Paper Plots
###############################################################################
plotdirFinal <- '~/Desktop/scratchFigures/ensoPaper/finalFigures'

# get the model inds First
inds <- c(which(modelNameVec=='CESM1-1-CAM5-CMIP5'), which(modelNameVec=='GISS-E21G'))

pal <- c('grey40',brewer.pal(9,'YlOrRd')[c(5,8)])


pdf(sprintf('%s/Figure_01_b.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
			xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
			main=sprintf('DJF Target Forecast (Perfect Model)'))

abline(h=0)
abline(v=c(12,24),lty=3)

for(i in 1:nquants){
	# Plot DJF ROC Scores as a function of lead
	points(xVec, fillInVec(masterRpssArray[i,,inds[1]]), col=pal[i], type='l', lwd=2)
	points(xVec, fillInVec(masterRpssArray[i,,inds[2]]), col=pal[i], type='l', lwd=2,lty=2)

}

legend('topright', rev(c('No EN', 'Weak EN', 'Strong EN')), col=rev(pal), lwd=2, bty='n')

legend(0.5, -0.05, c('CESM1-1-CAM5-CMIP5   ', 'GISS-E21G'), lwd=2, lty=c(1,2),
	horiz=T, bty='n', seg.len=3)

legend(-3.34, 1.08, '(b)', bty='n')

dev.off()

###############################################################################
# Figure 01c Paper Plot
###############################################################################

pdf(sprintf('%s/Figure_01_c.pdf',plotdirFinal),7,7)

plot(NULL, xlab='Lead (months)', ylab='Strong EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2,1),
		main=sprintf('Extra Strong EN Skill (Perfect Model)'))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[3,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)


legend(-3.34, 1.08, '(c)', bty='n')

dev.off()


###############################################################################
# Figure 01d Paper Plot
###############################################################################

# weak EN - no EN
pdf(sprintf('%s/Figure_01_d.pdf',plotdirFinal),7,7)

plot(NULL, xlab='Lead (months)', ylab='Weak EN RPSS - No EN RPSS',
		xaxp = c(0, maxInd, maxInd/12 * 2), xlim=c(0,maxInd), ylim=c(-0.2,1),
		main=sprintf('Extra Weak EN SKill (Perfect Model)'))
for(i in 1:nModels){
	points(xVec, fillInVec(masterRpssArray[2,,i]) - fillInVec(masterRpssArray[1,,i]),
		col=i, type='o',lwd=2, pch=rep(c(rep(NA,3),i,rep(NA,2)),maxInd/6))
}
abline(h=0)
abline(v=c(12,24),lty=3)

legend(-3.34, 1.08, '(d)', bty='n')

dev.off()


