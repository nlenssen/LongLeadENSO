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
	'CM2.5-NMME')

nModels <- length(modelNameVec)

# same for perfect
initModelNameVec <- modelNameVec
initExperimentName <- 'tos_zos_NEW'

# set other things
trialFields <- c('tos', 'zos')


# loop over models to create the step 1 thing
for(mod in 1:nModels){
	modelName <- modelNameVec[mod]
    
    # using this b/c we know it works
    initModelName <- 'CESM1-1-CAM5-CMIP5'

    source('Namelists/baseNamelist.Rnl')

	source('Code/01ProcessData.R')
}