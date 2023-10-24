###############################################################################
# Load in and clean data
###############################################################################

# get the grid info
handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,trialFields[1],modelName))

lon <- ncvar_get(handle,'lon')

lat <- ncvar_get(handle,'lat')
time <- ncvar_get(handle,'time')

sstField <- ncvar_get(handle,sprintf('%s_anom',trialFields[1]),start=c(1,1,1),count=c(-1,-1,1))
nlon <- length(lon)
nlat <- length(lat)
nt <- length(time)
nYear <- nt/12

nc_close(handle)

# figure out lon/lat Inds
lonInds <- which(lon > min(targetLon) & lon < max(targetLon))
latInds <- which(lat > min(targetLat) & lat < max(targetLat))

# make a maximal mask
maximalMask <- ifelse(!is.na(sstField),1,NA)

# make the observation/target sigma calculations for the prediction
cosMatFull <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180)) *
								maximalMask

cosMat <- cosMatFull[lonInds,latInds]
weightingMat <- cosMat/sum(cosMat,na.rm=T)


# make the library/trial sigma calculations for the prediction
trialSigma <- rep(NA, length(trialFields))

for(i in 1:length(trialFields)){
	handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,trialFields[i],modelName))
	trialArrayTemp <- arrMult(ncvar_get(handle,sprintf('%s_anom',trialFields[i]),
					start=c(min(lonInds),min(latInds),1),
					count=c(length(lonInds),length(latInds),-1)),maximalMask[lonInds,latInds])
	nc_close(handle)

	trialSigma[i] <- weighted.sd(c(trialArrayTemp),rep(c(weightingMat),dim(trialArrayTemp)[3]))

	rm(trialArrayTemp)
	gc()

}


###############################################################################
# Calculate ENSO climatology stuff (Define events)
###############################################################################
nino34LonInds <- which(lon > min(nino34Lon) & lon < max(nino34Lon))
nino34LatInds <- which(lat > min(nino34Lat) & lat < max(nino34Lat))

handle <- nc_open(sprintf('%s/%s/%s_anomaly_%s.nc',ddir,modelName,trialFields[1],modelName))

sstField <- ncvar_get(handle,'tos_anom',
	start=c(min(nino34LonInds),min(nino34LatInds),1),
	count=c(length(nino34LonInds),length(nino34LatInds),-1))

nc_close(handle)

nino34Monthly <- rep(NA, nt)

for(i in 1:nt){
	nino34Monthly[i] <- calculateIndex(sstField[,,i],lon[nino34LonInds], lat[nino34LatInds],
		lonRange=nino34Lon, latRange= nino34Lat)
}

nino34 <- rollmean(nino34Monthly,3,na.pad=T)

# Could do this differently!
eventCutoffs <- matrix(NA, 12, 2)

for(i in 1:12){
	monthInds <- which(timeMap[,2] == i)

	eventCutoffs[i,] <- quantile(nino34[monthInds], c(0.25, 0.75), na.rm=T) 	
}

# Define events based off of this

eventVec <- rep(0, nt)

for(i in 1:12){
	monthInds <- which(timeMap[,2] == i)

	eventVec[monthInds[which((nino34[monthInds] < eventCutoffs[i,1]))]] <- -1
	eventVec[monthInds[which((nino34[monthInds] > eventCutoffs[i,2]))]] <-  1

}

# Save just the stuff needed for initialization
save(timeMap, nino34, eventVec, trialSigma, file=sprintf('%s/step01_InitVerifInfo.Rda',wdir))

# save to make sure we have it (remove the restart flag in case it gets overwritten)
# Rethink this, messes with some important things on reload...
save(lon,lat,nlon,nlat,nt,nYear,lonInds,latInds,
	maximalMask,cosMatFull,cosMat,weightingMat,trialSigma,
	nino34LonInds, nino34LatInds, nino34Monthly, nino34, eventCutoffs, eventVec,
	file=sprintf('%s/step01_KeyInfo.Rda',wdir))





