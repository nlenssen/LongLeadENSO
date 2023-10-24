modelNameVec <- c('CM4',
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
	'CM2.5-NMME')
nModels <- length(modelNameVec)

plotdir <- '~/Desktop/scratchFigures/ensoPaper/Fig01'
if(!dir.exists(plotdir)) dir.create(plotdir)

plotdirFinal <- '~/Desktop/scratchFigures/ensoPaper/finalFigures'


# Stuff for file paths from namelist
ddir   <- '/Users/nale4362/Documents/AnalogueData'
trialFields <- c('tos', 'zos')

# Other relevant vars from the namelist
forecastLength <- 5*12-1

###############################################################################
# Load Perfect Model Verification
###############################################################################

rpssListPerfect <- list()
relListPerfect <- list()
resListPerfect <- list()

rocListPerfect  <- list()

for(i in 1:nModels){
	modelName <- modelNameVec[i]

	# load in the results

	experimentName <- 'tos_zos'
	wdir <- sprintf('%s/%s/experiments/%s_NEW',ddir, modelName, experimentName)

	load(sprintf('%s/step04_AllForecastVerification.Rda', wdir))

	# save in the output list
	rpssListPerfect[[i]]  <- rpssArray
	
	relListPerfect[[i]] <- reliabilityArray
	resListPerfect[[i]] <- resolutionArray

	rocListPerfect[[i]]  <- rocScoreArray
}

###############################################################################
# Load Cross Model Verification #1
###############################################################################

rpssListCross <- list()
relListCross <- list()
resListCross <- list()

rocListCross  <- list()

for(i in 1:nModels){
	modelName <- modelNameVec[i]

	initModelName <- 'CESM1-1-CAM5-CMIP5'

	experimentName <- sprintf('CrossPerfectModel_%s_%s', initModelName, paste(trialFields, collapse='_'))

	# Make working directory and sub-figure directory if they don't exist
	wdir <- sprintf('%s/%s/experiments/%s_NEW',ddir, modelName, experimentName)


	# load in the result
	load(sprintf('%s/step04_AllForecastVerification.Rda', wdir))

	# save in the output list
	rpssListCross[[i]] <- rpssArray

	relListCross[[i]] <- reliabilityArray
	resListCross[[i]] <- resolutionArray

	rocListCross[[i]]  <- rocScoreArray
}


###############################################################################
# Load Cross Model Verification #2
###############################################################################

rpssListCross2 <- list()
relListCross2 <- list()
resListCross2 <- list()

rocListCross2  <- list()

for(i in 1:nModels){
	modelName <- modelNameVec[i]

	initModelName <- 'GISS-E21G'

	experimentName <- sprintf('CrossPerfectModel_%s_%s', initModelName, paste(trialFields, collapse='_'))

	# Make working directory and sub-figure directory if they don't exist
	wdir <- sprintf('%s/%s/experiments/%s_NEW',ddir, modelName, experimentName)


	# load in the result
	load(sprintf('%s/step04_AllForecastVerification.Rda', wdir))

	# save in the output list
	rpssListCross2[[i]] <- rpssArray

	relListCross2[[i]] <- reliabilityArray
	resListCross2[[i]] <- resolutionArray

	rocListCross2[[i]]  <- rocScoreArray
}



###############################################################################
# Figure 1: Perfect and Cross DJF Nino3.4 Skill (RPSS)
###############################################################################

# subset the lead time
maxInd <- 36
xVec <- 0:maxInd
subInds <- which(0:forecastLength < maxInd+1)

# set plotting parameters
pal <- 1:nModels


pdf(sprintf('%s/01RpssFigs.pdf',plotdir),21,12)
par(mfrow=c(2,3))

##########
# plotting the perfect model skill
##########
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Perfect Model)')

abline(h=0)
for(i in 1:nModels){
	points(xVec, rpssListPerfect[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(24,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')

##########
# cross model #1
##########
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Cross-Model, CESM1.1-CAM5-CMIP5)')

abline(h=0)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(24,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')

##########
# cross model #2
##########
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Cross-Model, GISS-E21G)')

abline(h=0)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross2[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(24,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')


##########
# plotting the perfect model decomp
##########
plot(NULL, xlab='Lead (months)', ylab='Resolution/Reliability',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.1, 0.5),
		main='RPSS Decomposition (Perfect Model)')

abline(h=0)
for(i in 1:nModels){
	points(xVec, resListPerfect[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
	points(xVec, - relListPerfect[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6),lty=2)
}

legend(24,0.525,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')

##########
# cross model #1 decomp
##########
plot(NULL, xlab='Lead (months)', ylab='Resolution/Reliability',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.1, 0.5),
		main='RPSS Decomposition (Cross-Model, CESM1.1-CAM5-CMIP5)')

abline(h=0)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, resListCross[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
	points(xVec, - relListCross[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6),lty=2)
}

legend(24,0.525,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')

##########
# cross model #2 decomp
##########
plot(NULL, xlab='Lead (months)', ylab='Resolution/Reliability',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.1, 0.5),
		main='RPSS Decomposition (Cross-Model, CESM1.1-CAM5-CMIP5)')

abline(h=0)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, resListCross2[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
	points(xVec, - relListCross2[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6),lty=2)
}

legend(24,0.525,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')


dev.off()




##########
# New Figure 1a
##########
maxInd <- 36
xVec <- 0:maxInd
subInds <- which(0:forecastLength < maxInd+1)

# set plotting parameters
pal <- 1:nModels

pdf(sprintf('%s/Figure_01_a.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Perfect Model)')

abline(h=0)
abline(v=c(12,24),lty=3)

for(i in 1:nModels){
	points(xVec, rpssListPerfect[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(14,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')
legend(-3.34, 1.08, '(a)', bty='n')
dev.off()



##########
# New Figure 3a
##########

pdf(sprintf('%s/Figure_03_a.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Cross-Model, GISS-E21G)')

abline(h=0)
abline(v=c(12,24),lty=3)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross2[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(18,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')
legend(-3.34, 1.08, '(a)', bty='n')
dev.off()



##########
# New Figure 3d
##########
pdf(sprintf('%s/Figure_03_d.pdf',plotdirFinal),7,7)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Cross-Model, CESM1.1-CAM5-CMIP5)')

abline(h=0)
abline(v=c(12,24),lty=3)

for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross[[i]][1,subInds],type='o',lwd=2,col=pal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(18,1.05,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')
legend(-3.34, 1.08, '(d)', bty='n')
dev.off()



###############################################################################
# Figure 1: Scatter of Perfect model v Cross model skill at various leads
###############################################################################

# leadVec <- (0:24)+1
# plotLeads <- seq(3,24,by=3) + 1
# rSqVec <- rep(NA, length(leadVec))

# par(mfrow=c(2,4))

# for(j in 1:length(leadVec)){

# 	plotMat <- matrix(NA, nModels,2)
# 	for(i in 1:nModels){
# 		plotMat[i,1] <- rpssListPerfect[[i]][1,leadVec[j]]
# 		plotMat[i,2] <- rpssListCross[[i]][1,leadVec[j]]
		
# 	}	
# 	rSqVec[j] <- cor(plotMat)[1,2]^2
		
# 	if(j %in% plotLeads){
# 		plot(NULL, xlab='Perfect Model RPSS', ylab='Cross-model RPSS',
# 		main=sprintf('Lead %02d',leadVec[j]-1),
# 		xlim=c(-0.2,1), ylim=c(-0.2,1))

# 		points(plotMat,col=1:nModels,pch=1:nModels)

# 		abline(h=0, lty=3)
# 		abline(v=0, lty=3)
# 		abline(0,1,lty=1)
# 		legend(-0.2,1,modelNameVec,col=pal, pch=1:nModels, lwd=2,bty='n')
# 		legend('bottomright',sprintf('R-Squared: %.02f', rSqVec[j]),bty='n')
# 	}
# }

# check 
# plot(leadVec-1, rSqVec, type='b')
