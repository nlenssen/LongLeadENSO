ddir <- '/Users/nale4362/Documents/AnalogueData/Observations/Raw/CERA20C'
ofdir <- '/Users/nale4362/Documents/AnalogueData/Observations'
# variable names
var_fname  <- c('sosstsst','sossheig')
var_ofname <- c('tos', 'zos')
dlNameVec <- c('SST Anomaly', 'SSH Anomaly')
unitVec   <- c('DegC', 'cm')

# useful testing ind (January 1998 El Nino)
# i <- 1165

# Stuff needed for regridding
regrid <- FALSE
targetGrid <- '/Users/nale4362/Documents/AnalogueData/mygrid_2x2'

# hard code in time stuff
tYear <- 1901:2009
timeMap <- cbind(rep(tYear,each=12),1:12,NA)
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12

# first loop over all files and make them 2x2
if(regrid){
	files <- system(sprintf('ls %s/NOAA_1x1/*',ddir), intern=T)

	for(i in 1:length(files)){
		ofname <- sprintf('%s/2x2/%s_2x2.nc',ddir,substr(files[i],74,94))
		cdoRemapcon(targetGrid, files[i], ofname)
	}

}
# pull one file to get grid info
handle <- nc_open(sprintf('%s/2x2/sosstsst_c20c_1m_1901_2x2.nc',ddir))
lon  <- ncvar_get(handle, 'lon')
lat  <- ncvar_get(handle, 'lat')

nlon <- length(lon)
nlat <- length(lat)

nc_close(handle)

# loop over variables and process
for(v in 1:length(var_fname)){

	# First, pull the raw data
	files <- system(sprintf('ls %s/2x2/%s*.nc',ddir, var_fname[v]), intern=T)
	nFiles <- length(files)

	dataArray <- array(NA, dim=c(nlon,nlat,nFiles*12))
	nt <- dim(dataArray)[3]

	for(i in 1:nFiles){

		handle <- nc_open(files[i])
		tempData  <- ncvar_get(handle, var_fname[v])
		nc_close(handle)

		tempInds <- ((i-1)*12+1):(i*12)
		dataArray[,,tempInds] <- tempData
	}

	# Now, calculate anomalies using "sliding climatology"
	dataAnomalies <- array(NA, dim=dim(dataArray))

	# loop over tYear
	for(i in 1:length(tYear)){

		# get the inds of this year
		yearInds <- ((i-1)*12+1):(i*12)

		# get the moving climatology period
		if(tYear[i] < 1931){
			climYears <- 1901:1930
		} else{
			climYears <- (tYear[i]-30):(tYear[i]-1)
		}

		# calculate the anomaly
		for(j in 1:12){
			tempInds <- which(timeMap[,2]==j & timeMap[,1] %in% climYears)
			climatology <- apply(dataArray[,,tempInds],c(1,2),mean)

			dataAnomalies[,,yearInds[j]] <- dataArray[,,yearInds[j]] - climatology
		}

	}


	# write the netcdf file in the original format
	londim   <- ncdim_def("lon","degrees_east",as.double(lon)) 
	latdim   <- ncdim_def("lat","degrees_north",as.double(lat)) 
	timeDim   <- ncdim_def("time",'months from 1901',as.double(1:nrow(timeMap)),unlim=TRUE)

	# define anomaly variable
	fillvalue <- 1e32
	dlname <- dlNameVec[v]
	anom_def <- ncvar_def(sprintf('%s_anom',var_ofname[v]),unitVec[v],list(londim,latdim,timeDim),
							fillvalue,dlname,prec="single")

	# create netCDF file and put arrays
	ncfname <- sprintf('%s/%s_anomaly_Observations.nc',ofdir, var_ofname[v])
	if(file.exists(ncfname)) file.remove(ncfname)
	ncout <- nc_create(ncfname,list(anom_def),force_v4=T)

	# put variables
	ncvar_put(ncout,anom_def,dataAnomalies)

	# put additional attributes into dimension and data variables
	ncatt_put(ncout,"lon","axis","X")
	ncatt_put(ncout,"lat","axis","Y")
	ncatt_put(ncout,"time","axis","T")

	# add global attributes
	ncatt_put(ncout,0,"title",sprintf('CERA20C Sea Surface %s (Rolling Clim)',dlname))
	ncatt_put(ncout,0,"institution",'GISS/IRI')
	history <- paste("N. Lenssen", date(), sep=", ")
	ncatt_put(ncout,0,"history",history)

	# close the file, writing data to disk
	nc_close(ncout)

}



