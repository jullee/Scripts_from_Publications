########################## GOALS AND NOTES ###########################

### Project: Freeman et al. GEB 2018

### Notes: Functions used to calcualte change in area (within a study region) resulting from species' elevational shifts 
## An example is included but to reproduce the full dataset requires downloading files from source and processing with "dem_processing_script.R"

### Date: July 11, 2017

### Version of R:  R 3.3.3  

########################### END SECTION ##############################

############################ LIBRARIES ###############################

library(maptools)
library(raster)

########################### END SECTION ##############################

###################### FUNCTION DEFINITIONS ##########################


### 1) Function to convert DEM into binary raster showing only cells of specified elevational range


elev_subset<-function(dem, min, max){
	dem[dem<min|dem>max]<-NA
	dem[!(is.na(dem))]<-1
	return(dem)
}

### 2) Function to calculate differences in area

# name = species' name
# mountain_dem = projected and cropped dem of the survey mountain
# hist_warm = species' historical min
# hist_cold = species' historical max
# modern_warm = species' modern min
# modern_cold = species' modern max

area_change<-function(mountain_dem,proj,hist_warm,hist_cold,modern_cold,modern_warm,tmpdir){
	
	rasterOptions(tmpdir=tmpdir)
	
	# project raster
	mountain_dem_ea<-projectRaster(mountain_dem,crs=proj,filename="./proj_temp/projectedRaster")
	
	# cell area
	cellarea<-prod(res(mountain_dem_ea))
	
	# historic area
	hist_range<-elev_subset(mountain_dem_ea,min=hist_warm,max=hist_cold)
	hist_cells<-cellStats(hist_range,'sum',na.rm=TRUE)
	hist_area<-hist_cells*cellarea
	rm(hist_range)
	do.call(file.remove,list(normalizePath(list.files(path=tmpdir,full.names =TRUE))))
	
	# modern area
	modern_range<-elev_subset(mountain_dem_ea,min=modern_warm,max=modern_cold)
	modern_cells<-cellStats(modern_range,'sum',na.rm=TRUE)
	modern_area<-modern_cells*cellarea
	rm(modern_range)
	do.call(file.remove,list(normalizePath(list.files(path=tmpdir,full.names =TRUE))))
	
	do.call(file.remove,list(normalizePath(list.files(path=paste(getwd(),"/proj_temp",sep=""),full.names =TRUE))))
	
	# change in area
	delta_cells<-modern_cells-hist_cells
	delta_area<-modern_area-hist_area
	
	return(cbind(hist_cells,modern_cells,delta_cells,hist_area,modern_area,delta_area))
		
}

########################### END SECTION ##############################

########################## DATA AND ANALYSIS #########################

## Example using the script ##

setwd() # set to location of processed rasters

rast<-raster("./Example/Freeman2014_Karkar.tif")
prj<-"+proj=aea +lat_1=4 +lat_2=-8 +lat_0=0 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
area_change(mountain_dem=rast,proj=prj,hist_warm=800,hist_cold=1800,modern_cold=1600,modern_warm=800,tmpdir=tempdir())


########################### END SECTION ##############################

########################## FINAL COMMENTS ############################

########################### END SCRIPT ###############################
