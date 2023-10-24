library(colorout)

modelNameVec <- c('CM4',
	'ESM4',
	'CanESM5',
	'MIROC6',
	'CESM1-1-CAM5-CMIP5',
	'CESM2', 
	# 'GISS-E3', 
	'GISS-E21G', 
	# 'GISS-E21H', 
	'CESM1-NMME', 
	'CCSM4-NMME',
	'CM2.1-NMME', 
	'CM2.5-NMME'
	)


nModels <- length(modelNameVec)


# initModelNameVec <- 'CM2.5-NMME'
initModelNameVec <- 'Observations'

initExperimentName <- 'tos_zos_NEW'

# set other things
trialFields <- c('tos', 'zos')

# set the steps to run
runSteps <- c(1:4)

# Use true if not running from step 1 through steps of interest
restart <- TRUE

# loop over models to create the step 1 thing
for(mod in 2:nModels){
	# set the two model names
	modelName <- modelNameVec[mod]
    initModelName <- initModelNameVec

    print(sprintf('%02d/%02d Model: %s, Init: %s',mod, nModels, modelName, initModelName))

    # run the analyses
    source('Namelists/baseNamelist.Rnl')
    source('Scripts/runSingleModel.R')

    gc()
}