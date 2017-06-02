########################## GOALS AND NOTES ###########################

### Project: lyrata niche models

### Goal: Build niche model using Maxent. 

### Notes  

# This script details the basic maxent analysis with setting of features and regularization parameters
# This step follows the tuning step used to find the optimal features and regularization parameter 
# Path names and file names need to be adjusted when working with the figshare data

### Date: July 20, 2015

### Version of R:  R 3.2.0  

########################### END SECTION ##############################

############################ LIBRARIES ###############################

library(dismo)
library(maptools)
library(rgdal)
library(raster)

########################### END SECTION ##############################

###################### FUNCTION DEFINITIONS ##########################

### 1. Function to mask the first raster of a stack by a specified background extent and return the stack
# where:
# background = either a shapefile or raster of the desired background extent
# layers = the intial raster stack; the first layer of which will get masked

maskLayersToBackground<-function(background, layers){
	masked<-mask(layers[[1]], background)
	newStack<-stack(masked,layers[[2:nlayers(layers)]])
	names(newStack)[1]<-names(layers)[1]
	return(newStack)
}

### 2. Function to build niche models with k-fold cross validation

## Need copy of maxent jar file in R dismo library folder (see dismo documentation)
## Pre-processing of input data to same extent, projection etc. should have already been completed

## Input Variables:
# loc = the locality data (2 columns of x and y) 
# k = desired number of folds
# calibrationRasters = a raster stack to be used during model calibration (e.g. after processing with maskLayersToBackground) 
# out_dir = path to directory where results are to be stored (default is current directory)
# out_name = prefix for each output folder and file (e.g. 'range_MCP_ModelOutput'; default is "modelOutput")
# basicArgs = c("-J","-P","removeDuplicates=TRUE","writebackgroundpredictions","maximumiterations=5000") #a vector defining basic arguments to pass to maxent
# featureOmission = "" #a vector defining feature types to be omitted (e.g. c("noproduct","nohinge")
# reg = 1.0 #a number (double) defining the regularization parameter 
# predictRaster = can either be the same as 'calibrationRasters' (default) or a new stack of rasters (e.g. representing different time or place) to transfer the maxent model to
# factors = a vector of the names of variables which should be treated as categorical (default is empty vector)

## Outputs:
# This function outputs a folder for each of the k iterations of the modelling process, each of which contains all the standard output of maxent (including the maxent.html, plots and the species.lambdas files) as well as: 
# a) the coordinates for both the background points and locality data used in the modelling; 
# b) the maxent object from R (can be recalled using dget); 
# c) the raw and logistic raster surfaces
# d) the evaluation object (can be recalled using dget)

# The function also writes the following to the output directory: 
# a) a csv file of the coordinates of the background points used in model evaluation; 
# b) a csv file with the fold information for all localities (these two files are useful if you want to evaluate another set of models based on the same input localities using the same evaluation dataset); 
# c) the model evaluation object with various metrics 
# d) a table with the model evaluation results for each iteration of modelling; 
# e) a list of the input arguments


maxent_kfoldCrossEval<-function (loc,k,calibrationRasters,out_dir =".",out_name="ModelOutput",basicArgs = c("-J","-P","removeDuplicates=TRUE","writebackgroundpredictions","maximumiterations=5000"), featureOmission = "", reg = 1.0,predictRasters=calibrationRasters,factors=vector()){
	
	# Step 1: Process the locality data into folds, get sample of background points
		fold<-kfold(loc,k)
		backgr<-randomPoints(calibrationRasters, length(loc[fold==1,][,1]), p=loc) # exclude presences from background sample
	
	# Step 2: Initial a list to store evaluation results and set the arguments for input into Maxent	
		replicateEval1<-list() # initialize a list to store AUC and evaluation data
		max_args<-c(basicArgs,featureOmission, paste("betamultiplier=",reg,sep=""))
		
	# Step 3: Build model and make a prediction surface for every combination of k-1 folds	
		for (i in 1:k){
			dir<-paste(out_dir,"/",out_name,i,sep="") 
			train<-loc[fold!=i,]
			test<-loc[fold==i,]
				# run Maxent
				if(length(factors)==0){
					me<-maxent(calibrationRasters,train, args=max_args, path=dir) 
				} else {
					me<-maxent(calibrationRasters,train, factors=factors, args=max_args, path=dir) 
				}
				save(me,file=paste(dir,"/maxentModelObject_",i,sep="")) # save the maxent object
				
				# evaluate the models and put AUC results in list 
				replicateEval1[[i]]<-evaluate(p=test, a=backgr, me, x=calibrationRasters) 
				dput(replicateEval1[[i]],file=paste(dir,"/modelEvalObject_",i,sep="")) # save evaluation object so can use other metrics besides AUC
				
				# write testing and background data used in evaluation
				write.csv(backgr,file=paste(dir,"/backgroundPts_Model",i,sep=""),row.names=FALSE)
				write.csv(test,file=paste(dir,"/testLocs_Model",i,sep=""),row.names=FALSE)
				
				# make prediction rasters (both logistic and raw)
				predict(me, predictRasters, progress="", filename=paste(dir,"/",out_name,i,"_prediction_logistic",sep=""),format="GTiff")
				predict(me, predictRasters,args="outputformat=raw", progress="", filename=paste(dir,"/",out_name,i,"_prediction_raw",sep=""),format="GTiff")
				}
				
	# Step 4: Write AUC results, background data points used across all models, locality fold information and input parameters to file
	
		ModelAUC1<-sapply(replicateEval1,function(x){slot(x,'auc')})
		table1<-as.data.frame(ModelAUC1)
		table1$Model<-row.names(table1)
		table1<-table1[,c(2,1)]
		write.csv(table1,file=paste(out_dir,"/",out_name,"_ModelEvaluationResults.csv",sep=""),row.names=FALSE)
        write.csv(backgr,file=paste(out_dir,"/",out_name,"_backgr_evalPts.csv",sep=""),row.names=FALSE)   
        locFolds<-cbind(loc,fold)
        write.csv(locFolds,file=paste(out_dir,"/",out_name,"_localityFolds_evalPts.csv",sep=""),row.names=FALSE)
        
        # also write maxent input parameters
        write("Maxent_Arguments:",file=paste(out_dir,"/",out_name,"_maxentParameters.txt",sep=""))
        write(max_args,file=paste(out_dir,"/",out_name,"_maxentParameters.txt",sep=""),append=TRUE)
}


########################### END SECTION ##############################

########################## DATA AND ANALYSIS #########################

# modify this path as appropriate
setwd("/Users/Julie/Documents/Research/Post-Doc_Switzerland/lyrata_NicheModel/Data_and_Analysis") 
path<-getwd()

### Locality Data

myPts<-read.csv("./LocalityData/allSites_aea_onshore_noveltyThinned_MahalanobisDistances.csv",header=TRUE)
loc<-myPts[,c("Long","Lat")]

### Environmental Layers
# final environmental layers (e.g. after correlation analysis etc.)

setwd("./EnvironmentalLayers")
files<-list.files(pattern="aea.tif")
elayers<-stack(files)
setwd(path)


### Background Extent

# Ecoregions
ECO_range<-raster("./BackgroundExtent/ecoregions_background_aea.tif")


### Set optimal features and regularization parameter for each background
# This is based on tuning results

ECO_features<-c("nothreshold","nohinge")
ECO_reg=0.5 

### Setting up parameters for function

range<-get(ECO_range)
f<-get(ECO_features)
r<-get(ECO_reg)

calibrationRasters<-maskLayersToBackground(background=range,layers=elayers)
od<-paste(path,"/Maxent_Models/All_Pops_Model/ModelResults/ECO_background",sep="")
dir.create(od)
on<-"ModelOutput_ECO_"
	
### Run models

maxent_kfoldCrossEval(loc=loc,k=10,calibrationRasters=calibrationRasters,out_dir=od,out_name=on,featureOmission=f,reg=r,predictRasters=elayers)

########################### END SECTION ##############################

########################## FINAL COMMENTS ############################

# NA

########################### END SCRIPT ###############################