modelNameVec <- c('CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5',
	'CESM2', 
	'GISS-E21G', 
	'CESM1-NMME', 
	'CCSM4-NMME',
	'CM2.1-NMME', 
	'CM2.5-NMME')
nModels <- length(modelNameVec)




if(!exists('baseplotdir')){
	baseplotdir <- '~/Dropbox/CUBoulder/EnsoPaper'
	if(!dir.exists(baseplotdir)) dir.create(baseplotdir, recursive = TRUE)
}

if(!exists('plotdirFinal')){
	plotdirFinal <- '~/Dropbox/CUBoulder/EnsoPaper/finalFigures'
	if(!dir.exists(plotdirFinal)) dir.create(plotdirFinal, recursive = TRUE)
}

# Stuff for file paths from namelist
ddir   <- '/Users/nale4362/Data/AnalogueData'
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


##########
# New Figure 1a
##########
maxInd <- 36
xVec <- 0:maxInd
subInds <- which(0:forecastLength < maxInd+1)

textSize <- 1.5

# set plotting parameters


pdf(sprintf('%s/Figure_01_a.pdf',plotdirFinal),7,7)
# b,l,t,r Default: c(5, 4, 4, 2) + 0.1
par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF Probabilistic ENSO Skill (Perfect Model)',
		cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

abline(h=0)
abline(v=c(12,24),lty=3)

for(i in 1:nModels){
	points(xVec, rpssListPerfect[[i]][1,subInds],type='o',lwd=3, cex=1.5,
		col=qualPal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(17,1.06,modelNameVec,col=qualPal, pch=1:nModels, lwd=3,bty='n', cex = 1.2)
legend(-3.34, 1.08, '(a)', bty='n', cex = textSize)
dev.off()



##########
# New Figure 3a
##########

pdf(sprintf('%s/Figure_03_a.pdf',plotdirFinal),7,7)
par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF ENSO Skill (Cross-Model, GISS-E21G)',
		cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

abline(h=0)
abline(v=c(12,24),lty=3)
for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross2[[i]][1,subInds],type='o',lwd=3, cex=1.5,
		col=qualPal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(17,1.06,modelNameVec,col=qualPal, pch=1:nModels, lwd=2,bty='n', cex = 1.2)
legend(-3.34, 1.08, '(a)', bty='n', cex = textSize)
dev.off()



##########
# New Figure 3d
##########
pdf(sprintf('%s/Figure_03_d.pdf',plotdirFinal),7,7)
par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(NULL, xlab='Lead (months)', ylab='RPSS',
		xaxp = c(0, maxInd, 6), xlim=c(0,maxInd), ylim=c(-0.2, 1.0),
		main='DJF ENSO Skill (Cross-Model, CESM1.1)',
		cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

abline(h=0)
abline(v=c(12,24),lty=3)

for(i in 1:nModels){
	# if(modelNameVec[[i]] == 'CESM1-NMME' | modelNameVec[[i]] == 'CESM1-1-CAM5-CMIP5') next()
	points(xVec, rpssListCross[[i]][1,subInds],type='o',lwd=3, cex=1.5,
		col=qualPal[i],pch=rep(c(i,rep(NA,5)),maxInd/6))
}

legend(17,1.06,modelNameVec,col=qualPal, pch=1:nModels, lwd=2,bty='n', cex = 1.2)
legend(-3.34, 1.08, '(d)', bty='n', cex = textSize)
dev.off()

