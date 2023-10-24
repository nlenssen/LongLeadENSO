if(!exists('runSteps')){
	runSteps <- 1:2
}


###############################################################################
# Run the process
###############################################################################

timeVec <- rep(NA, 5)

# Step 1: Process the model and obs data (< 10 Seconds)
# 	- calculate anomalies
# 	- get time info into timeMap format
# 	- calculate nino3.4 indices
if(sum(runSteps==1) > 0){
	timeVec[1] <- system.time(source(sprintf('%s/01ProcessData.R',cdir)))[3]
}

# Step 2: Get the top K analogues given: (minutes to hours, parallelized)
# 	- library (set by `initModelName` with size given by `libraryInds`)
# 	- initial states
# 	- fields of interest
# 	- weighting/penalty function
if(sum(runSteps==2) > 0){
	timeVec[2] <- system.time(source(sprintf('%s/02BestAnalogues.R',cdir)))[3]
}


# Step 3: Get the forecast evolution of interest (< 10 Seconds)
# 	- full field vs. Nino3.4 (Currently only getting Nino3.4)
# 	- other variables?
# 	- both prior and post the initial condition?
if(sum(runSteps==3) > 0){
	timeVec[3] <- system.time(source(sprintf('%s/03GetForecasts_Nino34.R',cdir)))[3]
}

# Step 4: Verify the forecast in some way (< 10 Seconds)
# 	- ROC for event detection
# 	- correlation/MSESS for mean evolution
# 	- Other probablistic forecast metric? Something that accounts for sharpness maybe?
if(sum(runSteps==4) > 0){
	timeVec[4] <- system.time(source(sprintf('%s/04VerifyForecasts.R',cdir)))[3]
}

# Step 5: Plotting/Other
if(sum(runSteps==5) > 0){
	timeVec[5] <- system.time(source(sprintf('%s/05VisulizeResults.R',cdir)))[3]
}


# Step 6: Bonus analysis 01
if(sum(runSteps==6) > 0){
	timeVec[6] <- system.time(source(sprintf('%s/06StratificationSkillAnalysis.R',cdir)))[3]
}
