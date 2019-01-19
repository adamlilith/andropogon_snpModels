### ANDROPOGON PHENO-GENETIC-BIOGEOGRAPHICAL MODELING
### Adam B Smith | Missouri Botanical Garden | 2018
### source('C:/Ecology/Drive/Research/Andropogon/Andropogon - Pheno-geno Modeling/Code/00 Andropogon Pheno-Genetic-Biogeographic Modeling.r')

	print(Sys.time())
	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

################
### CONTENTS ###
################

### libraries ###
### constants ###
### functions ###

### copy and collate pheno/geno data ###
### create study region extents and masks ###
### collate environmental data ###
### copy and modify species distribution data ###
### match species distribution data to environmental data ###
### construct PCA for environmental data ###
### match pheno/genotypic data to environmental data ###
### add populations to spatial polygon of species occurrences ###

### model species distribution using bayesLopod ###
### model future distribution of species from bayesLopod output ###
### map bayesLopod present and future prediction ###

### model SNPs ###
### assess SNP model performance and sensitivity to predictors ###
### EDA of SNP model performance and predictor importance ###
### identify alleles most likely affected by climate ###
### make maps of alleles most likely responding to climate ###
### make map of allelic diversity ###


# ### copy and process raw response data ###
# ### create data frame with rows representing INDIVIDUALS ###
# ### create data frame with rows representing POPULATIONS ###
# ### create study region ###
# ### collate environmental data ###
# ### extract environmental data ###
# ### merge shapefile with environmental data with shapefile with specimen data ###
# ### merge occurrence data and population-level phenotype and genotype data ###
# ### model species distribution using bayesLopod ###
# ### GJAM @ population level ###
# ### plot GJAM @ population level ###

#################
### libraries ###
#################

	library(geosphere)
	library(rgeos)
	library(raster)
	library(dismo)
	# library(bayesLopod)
	library(gjam)
	library(MuMIn)

	library(omnibus)
	library(enmSdm)
	library(legendary)
	library(statisfactory)
	
	library(scales)
	library(viridisLite)

	library(openxlsx)

	options(stringsAsFactors=FALSE)
	rasterOptions(format='GTiff', overwrite=TRUE)

	setwd('C:/ecology/Drive/Research/Andropogon/Andropogon - Pheno-geno Modeling')

#################
### constants ###
#################

	phenoVars <- c('height_cm', 'leafWidth_cm', 'biomass_g', 'spad')

	bioclims <- c(1, 12)
	# predictors <- c('bio02', 'bio05', 'bio07', 'bio12', 'solarRadiation')
	predictors <- c('bio01', 'bio12')
			
	ll <- c('longitude', 'latitude')

	rcps <- c('4pt5', '8pt5')
	
	buffSize <- 200 # size of buffer of focal region for plotting in km

#################
### functions ###
#################

	# stack environmental rasters
	stackCurrentEnv <- function(predictors) {
	
		# predictors	character vector of predictors matching raster file names

		env <- raster(paste0('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000/', predictors[1], '.tif'))
		
		if (length(predictors) > 1) {
			env <- stack(env, paste0('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000/', predictors[2:length(predictors)], '.tif'))
		}
		
		env
	
	}

	# stack environmental rasters
	stackFutureEnv <- function(predictors, rcp) {
	
		# predictors	character vector of predictors matching raster file names
		# rcp			RCP in format as per '4pt5'
	
		gcms <- read.csv('./Data/WORLDCLIM Ver 1pt4 Rel 3/GCMs.csv')
		gcms <- gcms$gcm
		allEnv <- list()

		for (gcm in gcms) {
		
			env <- raster(paste0('./Data/WORLDCLIM Ver 1pt4 Rel 3/RCP ', rcp, ' ', gcm, ' 2061 to 2080/', predictors[1], '.tif'))
			
			if (length(predictors) > 1) {
				env <- stack(env, paste0('./Data/WORLDCLIM Ver 1pt4 Rel 3/RCP ', rcp, ' ', gcm, ' 2061 to 2080/', predictors[2:length(predictors)], '.tif'))
			}
			
			allEnv[[length(allEnv) + 1]] <- env
			name <- gsub(gcm, pattern='-', replacement='')
			name <- gsub(name, pattern='!', replacement='')
			name <- tolower(name)
			names(allEnv)[length(allEnv)] <- name
			
		}
		
		allEnv
	
	}

	# get names of coefficients (excluding intercept, including 2-way interactions)
	getCoeffNames <- function(predictors) {

		# predictors	vector of predictor names
	
		coeffNames <- predictors
		for (countPred1 in 1:(length(predictors) - 1)) {
			for (countPred2 in 2:length(predictors)) {
				
				pred1 <- predictors[countPred1]
				pred2 <- predictors[countPred2]
				ia <- paste0(pred1, 'X', pred2)
				coeffNames <- c(coeffNames, ia)
			
			}
		}
		
		coeffNames
		
	}

	# nice name and unit of each response variable
	respNice <- function(x) {

		if (x == 'biomass_g') {
			nice <- 'Biomass'
			unit <- ' (g)'
		} else if (x == 'leafWidth_cm') {
			nice <- 'Leaf Width'
			unit <- ' (cm)'
		} else if (x == 'height_cm') {
			nice <- 'Height'
			unit <- ' (cm)'
		} else if (x == 'spad') {
			nice <- 'SPAD'
			unit <- ''
		} else if (x == 'structureGroup1') {
			nice <- 'STRUCTURE Population 1'
			unit <- ''
		} else if (x == 'structureGroup2') {
			nice <- 'STRUCTURE Population 2'
			unit <- ''
		} else if (x == 'structureGroup3') {
			nice <- 'STRUCTURE Population 3'
			unit <- ''
		} else if (x == 'structureGroup4') {
			nice <- 'STRUCTURE Population 4'
			unit <- ''
		}
		
		list(nice=nice, unit=unit)
		
	}

	# nice name and unit of each response variable
	predNice <- function(x) {

		if (x == 'bio01') {
			nice <- 'BIO01'
			name <- 'Mean Annual Temperature'
			unit <- ' (deg C)'
		} else if (x == 'bio02') {
			nice <- 'BIO02'
			name <- 'Diurnal Temperature Range'
			unit <- ' (deg C)'
		} else if (x == 'bio05') {
			nice <- 'BIO05'
			name <- 'Maximum Temperature'
			unit <- ' (deg C)'
		} else if (x == 'bio07') {
			nice <- 'BIO07'
			name <- 'Temperature Range'
			unit <- ' (deg C)'
		} else if (x == 'bio12') {
			nice <- 'BIO12'
			name <- 'Mean Annual Precipitation'
			unit <- ' (mm)'
		} else if (x == 'solarRadiation') {
			nice <- 'Solar Radiation'
			name <- 'Solar Radiation'
			unit <- ' (MJ)'
		} else if (grepl(x, pattern='X')) {
			xIsWhere <- regexpr(x, pattern='X')[1]
			x1 <- substr(x, 1, xIsWhere - 1)
			x2 <- substr(x, xIsWhere + 1, nchar(x))
			
			out1 <- predNice(x1)
			out2 <- predNice(x2)
			
			nice <- paste(out1$nice, '*', out2$nice)
			name <- paste(out1$name, '*\n', out2$name)
			unit <- ''
		}
		
		list(nice=nice, name=name, unit=unit)
		
	}

	# aggregate individual-level site data to population level
	aggPhenoGenoByPop <- function(phenoGeno, sp=TRUE) {
		
		# phenoGeno		data frame of individual-level geno/pheno data
		# sp			TRUE: output is SPatialPointsDataFrame
		
		phenoGenoPops <- aggregate(phenoGeno[ , c('population', ll)], by=list(phenoGeno$population), mean)
		phenoGenoPops$population <- NULL
		names(phenoGenoPops)[1] <- 'population'
		nas <- naRows(phenoGenoPops[ , ll])
		if (length(nas) > 0) phenoGenoPops <- phenoGenoPops[-nas, ]
		
		if (sp) phenoGenoPops <- SpatialPointsDataFrame(phenoGenoPops[ , ll], data=phenoGenoPops, proj4string=getCRS('wgs84', TRUE))
		phenoGenoPops
		
	}
		
	# get polygon to represent extent of plotting region
	getPlotFocus <- function(buffSize) {
	
		# buffSize	size of buffer around populations in km
		
		load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
		
		# aggregate geno/pheno data to population
		phenoGenoPops <- aggPhenoGenoByPop(phenoGeno, sp=TRUE)
		
		# get focal region extent
		phenoGenoPopsEa <- sp::spTransform(phenoGenoPops, getCRS('albersNA', TRUE))
		focus <- gBuffer(phenoGenoPopsEa, width=buffSize * 1000)
		focus <- sp::spTransform(focus, getCRS('wgs84', TRUE))
		focus <- extent(focus)
		focus <- as(focus, 'SpatialPolygons')
		projection(focus) <- getCRS('wgs84')
	
		focus
	
	}

# say('########################################')
# say('### copy and collate pheno/geno data ###')
# say('########################################')

	# file.copy('./Data/Phenotypic & Genotypic Data/RAW/gbs2014 combined genotype phenotype data 01OCT2018.xlsx', './Data/Phenotypic & Genotypic Data/00 gbs2014 combined genotype phenotype data 01OCT2018.xlsx', overwrite=TRUE)

	# phenoGeno <- read.xlsx('./Data/Phenotypic & Genotypic Data/00 gbs2014 combined genotype phenotype data 01OCT2018.xlsx', sheet='Sheet1')

	# # rename some columns
	# phenoGenoNames <- names(phenoGeno)

	# thisName <- which(phenoGenoNames %in% 'pop')
	# names(phenoGeno)[thisName] <- 'population'
	
	# thisName <- which(phenoGenoNames %in% 'height.(cm).2014')
	# names(phenoGeno)[thisName] <- 'height_cm'
	
	# thisName <- which(phenoGenoNames %in% 'spad_avg_2014')
	# names(phenoGeno)[thisName] <- 'spad'
	
	# thisName <- which(phenoGenoNames %in% 'leaf_width_avg_2014')
	# names(phenoGeno)[thisName] <- 'leafWidth_cm'
	
	# thisName <- which(phenoGenoNames %in% 'veg.wt_2014')
	# names(phenoGeno)[thisName] <- 'vegMass_g'
	
	# thisName <- which(phenoGenoNames %in% 'boot.wt_2014')
	# names(phenoGeno)[thisName] <- 'bootMass_g'
	
	# thisName <- which(phenoGenoNames %in% 'seed.wt_2014')
	# names(phenoGeno)[thisName] <- 'seedMass_g'
	
	# thisName <- which(phenoGenoNames %in% 'total.weight_2014')
	# names(phenoGeno)[thisName] <- 'biomass_g'
	
	# phenoGenoNames <- names(phenoGeno)
	# for (i in seq_along(phenoGenoNames)) {
		
		# name <- phenoGenoNames[i]
		# nameStart <- substr(name, 1, 1)
		# asNum <- as.numeric(nameStart)
		
		# if (!is.na(asNum)) phenoGenoNames[i] <- paste0('snp', name)
	# }
	
	# names(phenoGeno) <- phenoGenoNames
			
	# save(phenoGeno, file='./Data/Phenotypic & Genotypic Data/01 Phenotypic & Genotypic Data - Renamed Columns.RData')

# say('#############################################')
# say('### create study region extents and masks ###')
# say('#############################################')

	# outDir <- './Study Region'
	# dirCreate(outDir)

	# ### study region polygons
	# usa1 <- raster::getData('GADM', country='USA', level=1, path=outDir)
	# can1 <- raster::getData('GADM', country='CAN', level=1, path=outDir)
	# mex1 <- raster::getData('GADM', country='MEX', level=1, path=outDir)
	
	# nam1 <- rbind(can1, usa1, mex1)
	# nam1 <- nam1[-which(nam1$NAME_1 %in% c('Alaska', 'Hawaii')), ]

	# usa0 <- raster::getData('GADM', country='USA', level=0, path=outDir)
	# can0 <- raster::getData('GADM', country='CAN', level=0, path=outDir)
	# mex0 <- raster::getData('GADM', country='MEX', level=0, path=outDir)
	
	# nam0 <- rbind(can0, usa0, mex0)
	# nam0 <- crop(nam0, nam1)

	# save(nam1, file=paste0(outDir, '/North America Level 1 Sans Alaska.RData'))
	# save(nam0, file=paste0(outDir, '/North America Level 0 Sans Alaska.RData'))
	
	# ### study region mask raster
	# elev <- raster('D:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/Elevation - 10 arcmin/elevation.tif')
	# mask <- crop(elev, nam0)
	# mask <- rasterize(nam0, elev)
	# mask <- crop(mask, nam0)
	# mask <- 1 + 0 * mask

	# names(mask) <- 'mask'
	# writeRaster(mask, './Study Region/mask')
	
# say('##################################')
# say('### collate environmental data ###')
# say('##################################')

	# mask <- raster('./Study Region/mask.tif')

	# ### elevation
	# #############
		
		# outDir <- './Data/Topography - WORLDCLIM 1pt4'
		# dirCreate(outDir)
		
		# elev <- raster('D:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/Elevation - 10 arcmin/elevation.tif')
		# elev <- crop(elev, mask)
		# elev <- elev * mask
		# names(elev) <- 'elevation'
		# writeRaster(elev, paste0(outDir, '/elevation'))
		
	# ### current climate (1970-2010)
	# ###############################

		# dirCreate('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000')
		
		# ### BIOCLIMS
		# ############
		
		# for (bio in bioclims) {
			
			# rast <- raster(paste0('D:/ecology/Climate/WORLDCLIM Ver 2 Rel June 1 2016/10 arcmin 1970 to 2000/wc2.0_bio_10m_', prefix(bio, 2), '.tif'))
			# rast <- crop(rast, mask)
			# rast <- rast * mask
			# rast <- setMinMax(rast)
			# names(rast) <- paste0('bio', prefix(bio, 2))
			# writeRaster(rast, paste0('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000/bio', prefix(bio, 2)))
			
		# }
		
	# ### future climate
	# ##################
	
	# gcmInfo <- read.csv('D:/Ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/GCM Names, Scenarios Modeled, and File Codes.csv')
	# gcmsInEns <- read.csv('D:/Ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/RCP4pt5 - 10 arcmin/2041 to 2060 - 10 arcmin - RCP4pt5 - !ENSEMBLE/!GCMs Used in Ensemble.csv')
	# gcms <- c(gcmsInEns[ , 2])
	# gcms <- c('!ENSEMBLE', gcms)
	
	# write.csv(data.frame(gcm=gcms), './Data/WORLDCLIM Ver 1pt4 Rel 3/GCMs.csv', row.names=FALSE)
	
	# for (rcp in rcps) {

		# for (gcm in gcms) {
	
			# outDir <- paste0('./Data/WORLDCLIM Ver 1pt4 Rel 3/RCP ', rcp, ' ', gcm, ' 2061 to 2080')
			# dirCreate(outDir)
		
			# rcpShort <- gsub(rcp, pattern='pt', replacement='')
			# gcmCode <- tolower(gcmInfo$code[gcmInfo$gcm == gcm])
			
			# say(rcp, ' ', gcm)
	
			# for (bio in bioclims) {
		
				# rast <- raster(paste0('D:/Ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/RCP', rcp, ' - 10 arcmin/2061 to 2080 - 10 arcmin - RCP', rcp, ' - ', gcm, '/World/', gcmCode, rcpShort, 'bi70', prefix(bio, 2), '.tif'))
				# rast <- crop(rast, mask)
				# rast <- rast * mask
				# if (bio < 12) rast <- rast / 10
				# rast <- setMinMax(rast)
				# names(rast) <- paste0('bio', prefix(bio, 2))
				# writeRaster(rast, paste0(outDir, '/bio', prefix(bio, 2)))
				
			# }
			
		# } # next GCM
	
	# } # next RCP
	

# say('#################################################')
# say('### copy and modify species distribution data ###')
# say('#################################################')

	# outDir <- './Data/Andropogon Occurrences'
	# dirCreate(outDir)
	
	# file.copy('C:/Ecology/Drive/Research/Andropogon/Analysis - Phenotype Modeling/Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata', paste0(outDir, '/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
	
	# load(paste0(outDir, '/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
	
	# # clean names and remove old columns
	# ag <- gadm
	
	# ag@data$areaKm2 <- NULL
	# ag@data$elevation <- NULL
	# ag@data$solarRadiation <- NULL
	# ag@data$wc01 <- NULL
	# ag@data$wc02 <- NULL
	# ag@data$wc05 <- NULL
	# ag@data$wc07 <- NULL
	# ag@data$wc12 <- NULL
	# ag@data$longitude <- NULL
	# ag@data$site <- NULL
	
	# unwanted <- which(grepl(names(ag), pattern='height'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='height'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='leafwidth'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='biomass'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='spad'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='julian'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='boot'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='structure'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='freqlocus'))
	# ag@data[ , unwanted] <- NULL

	# unwanted <- which(grepl(tolower(names(ag)), pattern='geneticsamplesize'))
	# ag@data[ , unwanted] <- NULL

	# agNames <- names(ag)
	# this <- which(agNames %in% 'NAME_0')
	# names(ag)[this] <- 'country'
	
	# this <- which(agNames %in% 'NAME_1')
	# names(ag)[this] <- 'stateProv'
	
	# this <- which(agNames %in% 'NAME_2')
	# names(ag)[this] <- 'county'
	
	# this <- which(agNames %in% 'numCrd1')
	# names(ag)[this] <- 'agRecordsQual1'
	
	# this <- which(agNames %in% 'numCrd2')
	# names(ag)[this] <- 'agRecordsQual2'
	
	# this <- which(agNames %in% 'numCrd3')
	# names(ag)[this] <- 'agRecordsQual3'
	
	# this <- which(agNames %in% 'anyAg1to3')
	# names(ag)[this] <- 'agRecordsQual1to3'
	
	# this <- which(agNames %in% 'poaRec')
	# names(ag)[this] <- 'poaceaeRecords'
	
	# this <- which(agNames %in% 'agDensity')
	# names(ag)[this] <- 'agRecordsQual1to3_numPerKm2'
	
	# this <- which(agNames %in% 'poaDensity')
	# names(ag)[this] <- 'poaceaeRecords_numPerKm2'
	
	# save(ag, file=paste0(outDir, '/!14 Andropogon and Poaceae Records on GADM.RData'))

# say('#############################################################')
# say('### match species distribution data to environmental data ###')
# say('#############################################################')

	# outDir <- './Data/Andropogon Occurrences'

	# load(paste0(outDir, '/!14 Andropogon and Poaceae Records on GADM.RData'))
	
	# envRasts <- raster::stack(
		# c(
			# listFiles('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000', pattern='.tif'),
			# './Data/Topography - WORLDCLIM 1pt4/elevation.tif'
		# )
	# )
	
	# env <- raster::extract(envRasts, ag, fun=mean)
	# colnames(env) <- names(envRasts)
	# ag <- cbind(ag, env)
	
	# save(ag, file=paste0(outDir, '/!15 Andropogon and Poaceae Records on GADM with BIOCLIM Variables.RData'))
	
# say('############################################')
# say('### construct PCA for environmental data ###')
# say('############################################')

	# outDir <- './Data/Andropogon Occurrences'
	# load(paste0(outDir, '/!15 Andropogon and Poaceae Records on GADM with BIOCLIM Variables.RData'))
	
	# agOnly <- ag[!is.na(ag$agRecordsQual1to3), ]
	# agOnly <- agOnly[agOnly$agRecordsQual1to3 > 0, ]
	# bios <- agOnly@data[ , predictors]
	# nas <- naRows(bios)
	# if (length(nas) > 0) bios <- bios[-nas, ]

	# # PCA
	# pca <- princomp(bios, cor=TRUE)
	# pcPrediction <- predict(pca, ag@data)
	# colnames(pcPrediction) <- paste0('pc', seq_along(predictors))
	
	# # save
	# ag <- cbind(ag, pcPrediction)
	# save(ag, file=paste0(outDir, '/!15 Andropogon and Poaceae Records on GADM with BIOCLIM Variables and PCA.RData'))
	
	# dirCreate('./PCA')
	# sink('./PCA/PCA.txt', split=TRUE)
		# say(paste0('PCA on ', paste(predictors, collapse=' '), ' in counties with Andropogon gerardii occurrences'), post=2)
		# print(summary(pca))
	# sink()
	
	# save(pca, file='./PCA/PCA on BIOCLIM Variables across Counties with Andropogon gerardii Occurrences.RData')
	
# say('########################################################')
# say('### match pheno/genotypic data to environmental data ###')
# say('########################################################')

	# load('./Data/Phenotypic & Genotypic Data/01 Phenotypic & Genotypic Data - Renamed Columns.RData')
	# load('./PCA/PCA on BIOCLIM Variables across Counties with Andropogon gerardii Occurrences.RData')
	
	# envRasts <- raster::stack(
		# c(
			# listFiles('./Data/WORLDCLIM Ver 2 Rel June 1 2016/1970-2000', pattern='.tif'),
			# './Data/Topography - WORLDCLIM 1pt4/elevation.tif'
		# )
	# )
	
	# env <- raster::extract(envRasts, phenoGeno[ , ll], fun=mean)
	# colnames(env) <- names(envRasts)
	
	# pcPrediction <- predict(pca, env)
	# colnames(pcPrediction) <- paste0('pc', seq_along(predictors))
	
	# phenoGeno <- insertCol(x=env, into=phenoGeno, at='snp120025_ref_T_alt_C')
	# phenoGeno <- insertCol(x=pcPrediction, into=phenoGeno, at='snp120025_ref_T_alt_C')
	
	# save(phenoGeno, file='./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	
# say('#################################################################')
# say('### add populations to spatial polygon of species occurrences ###')
# say('#################################################################')

	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	# load('./Data/Andropogon Occurrences/!15 Andropogon and Poaceae Records on GADM with BIOCLIM Variables and PCA.RData')

	# # aggregate geno/pheno data to population
	# phenoGenoPops <- aggPhenoGenoByPop(phenoGeno, sp=TRUE)
	
	# overs <- over(ag, phenoGenoPops)
	# ag@data <- insertCol(overs, into=ag@data, at='bio01')
	
	# save(ag, file='./Data/Andropogon Occurrences/!16 Andropogon and Poaceae Records with Sampled Populations.RData')
	
# say('###################################################')
# say('### model species distribution using bayesLopod ###')
# say('###################################################')

	# outDir <- './Models - Species/bayesLopod'
	# dirCreate(outDir)
	
	# load('./Data/Andropogon Occurrences/!15 Andropogon and Poaceae Records on GADM with BIOCLIM Variables and PCA.RData')
	
	# ag <- ag[-which(is.na(ag$poaceaeRec)), ]
	
	# gadmLopod <- shapeLopodData(ag, fieldN='poaceaeRec', fieldY='agRecordsQual1to3', Adjacency=TRUE, keepFields=TRUE)
	
	# modelVarP <- modelLopod(gadmLopod, varP=TRUE, CAR=TRUE)
	# modelConstantP <- modelLopod(gadmLopod, varP=FALSE, CAR=TRUE)
	
	# save(lopodVarP, file=paste0(outDir, '/Model - bayesLopod with Variable Detection.RData'))
	# save(lopodConstantP, file=paste0(outDir, '/Model - bayesLopod with Constant Detection.RData'))
	
	# lopodPredictVarP <- lopodShape(modelVarP, param='psi_i', extrapolate=TRUE, metric='mean')
	# lopodPredictConstantP <- lopodShape(modelConstantP, param='psi_i', extrapolate=TRUE, metric='mean')
	
	# save(lopodPredictVarP, file=paste0(outDir, '/Prediction - bayesLopod with Variable Detection.RData'))
	# save(lopodPredictConstantP, file=paste0(outDir, '/Prediction - bayesLopod with Constant Detection.RData'))

# say('###################################################################')
# say('### model future distribution of species from bayesLopod output ###')
# say('###################################################################')

	# # generalization
	# rcp <- '8pt5'
	# lopodBuffSize <- 400 # size of buffer around populations to crop LOPOD shape (in km)

	# # load('./Models - Species/bayesLopod/Prediction - bayesLopod with Constant Detection.Rdata')
	# load('./Models - Species/bayesLopod/Prediction - bayesLopod with Variable Detection.Rdata')
	
	# lopodPredictVarP$freqLocusS1_2050023allele110 <- lopodPredictVarP$freqLocusS1_27189028allele120 <- lopodPredictVarP$freqLocusS1_33964042allele100 <- lopodPredictVarP$freqLocusS1_73120031allele110 <- lopodPredictVarP$geneticSampleSize <- lopodPredictVarP$sampEffort <- NULL
	
	# # crop to Midwest
	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	# focus <- getPlotFocus(buffSize=buffSize)
	# lopodPredictVarP <- crop(lopodPredictVarP, focus)	
	
	# # setup for modeling
	# resp <- logitAdj(lopodPredictVarP@data[ , 'psi_i'], epsilon=0)
	# wc01 <- scale(lopodPredictVarP@data[ , 'wc01'])
	# wc12 <- scale(lopodPredictVarP@data[ , 'wc12'])
	# weight <- log10(lopodPredictVarP@data$areaKm2)
	# weight <- weight - min(weight) + 0.001 * diff(range(weight))
	
	# data <- data.frame(resp=resp, wc01=wc01, wc12=wc12, weight=weight)
	
	# model <- trainGlm(data=data, resp='resp', preds=c('wc01', 'wc12'), w='weight', family='gaussian')
	
	# center <- c(attr(wc01, 'scaled:center'), attr(wc12, 'scaled:center'))
	# scale <- c(attr(wc01, 'scaled:scale'), attr(wc12, 'scaled:scale'))

	# ### predict to future
	# env <- stackFutureEnv(predictors=c('bio01', 'bio12'), rcp=rcp)
	# thisEnv <- raster::extract(env[[1]], lopodPredictVarP, fun=mean)
	# thisEnv <- as.data.frame(thisEnv)
	# names(thisEnv) <- c('wc01', 'wc12')

	# thisEnv <- scale(thisEnv, center=center, scale=scale)
	# thisEnv <- as.data.frame(thisEnv)
	
	# predictions <- predict(model, thisEnv)
	# predictions <- probitAdj(predictions, epsilon=0)
	
	# psi_ensembleGcm_rcp85_2070s <- data.frame(psi_ensembleGcm_rcp85_2070s=predictions)
	# lopodPredictVarP@data <- cbind(lopodPredictVarP@data, psi_ensembleGcm_rcp85_2070s)
	
	# save(lopodPredictVarP, file='./Models - Species/bayesLopod/Prediction - bayesLopod with Variable Detection with Future Predictions.RData')

# say('####################################################')
# say('### map bayesLopod present and future prediction ###')
# say('####################################################')
	
	# # generalization
	# thold <- 0.90 # threshold psi above which to highlight a county
	
	# # load geno/pheno and GIS data
	# load('./Study Region/North America Level 1 Sans Alaska.RData')
	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	# load('./Models - Species/bayesLopod/Prediction - bayesLopod with Variable Detection with Future Predictions.RData')

	# # aggregate geno/pheno data to population
	# phenoGenoPops <- aggPhenoGenoByPop(phenoGeno, sp=TRUE)
	# focus <- getPlotFocus(buffSize=buffSize)
	
	# nam1 <- crop(nam1, focus)
	# lopodPredictVarP <- crop(lopodPredictVarP, focus)
	
	# colors
	# cols <- colorRampPalette(c('white', 'green4'))
	# cols <- cols(101)
	# for (i in seq_along(cols)) cols[i] <- alpha(cols[i], 0.85)
		
	# dirCreate('./Figures & Tables/Species Model - Maps')

	# plot
	# png(paste0('./Figures & Tables/Species Model - Maps/bayesLopod x Linear Model - Present and Future Predictions.png'), width=2 * 1200, height=1000, res=300)

		# par(mfrow=c(1, 2), oma=0.1 * c(1, 4, 1, 26), mai=0 * c(1, 1, 1, 1), mar=rep(0, 4))

		# ### present
		# ###########
		
		# preds <- lopodPredictVarP@data$psi_i
		# preds <- round(100 * preds) + 1
		# theseCols <- cols[preds]
		
		# minPred <- min(preds, na.rm=TRUE)
		# maxPred <- max(preds, na.rm=TRUE)

		# coreQuant <- quantile(lopodPredictVarP@data$psi_i, thold, na.rm=TRUE)
		# whichCore <- which(lopodPredictVarP@data$psi_i >= coreQuant)
		# rangeCore <- lopodPredictVarP[whichCore, ]
		# rangeCore <- gUnaryUnion(rangeCore)
		
		# plot(lopodPredictVarP, col=theseCols, breaks=minPred:maxPred, border=NA, axes=FALSE, box=FALSE)
		# plot(rangeCore, col='darkgreen', border='black', add=TRUE)
		# plot(nam1, add=TRUE)
		
		# cents <- gCentroid(lopodPredictVarP, byid=TRUE)
		# cents <- cents[lopodPredictVarP@data$anyAg1to3 > 0, ]
		# points(cents, pch=16, cex=0.2)
		
		# # title
		# usr <- par('usr')
		# x <- usr[1]
		# y <- usr[4] - 0.08 * (usr[4] - usr[3])
		# name <- paste0('a) Current climate')
		# text(x, y, labels=name, xpd=NA, cex=0.7, pos=4, font=1)

		# ### future
		# ###########
		
		# preds <- lopodPredictVarP@data$psi_ensembleGcm_rcp85_2070s
		# preds <- round(100 * preds) + 1
		# theseCols <- cols[preds]
		
		# minPred <- min(preds, na.rm=TRUE)
		# maxPred <- max(preds, na.rm=TRUE)
		
		# whichCore <- which(lopodPredictVarP@data$psi_ensembleGcm_rcp85_2070s >= coreQuant)
		# rangeCore <- lopodPredictVarP[whichCore, ]
		# rangeCore <- gUnaryUnion(rangeCore)
		
		# plot(lopodPredictVarP, col=theseCols, breaks=minPred:maxPred, border=NA, axes=FALSE, box=FALSE)
		# plot(rangeCore, col='darkgreen', border='black', add=TRUE)
		# plot(nam1, add=TRUE)
		
		# cents <- gCentroid(lopodPredictVarP, byid=TRUE)
		# cents <- cents[lopodPredictVarP@data$anyAg1to3 > 0, ]
		# points(cents, pch=16, cex=0.2)
		
		# # title
		# usr <- par('usr')
		# x <- usr[1]
		# y <- usr[4] - 0.08 * (usr[4] - usr[3])
		# name <- paste0('b) Future climate in 2070s under RCP8.5')
		# text(x, y, labels=name, xpd=NA, cex=0.7, pos=4, font=1)
		
		# # legend
		# labs <- c('0', '0.25', '0.50', '0.75', sprintf('%.2f', round(coreQuant, 2)))

		# legendGrad(
			# x='right',
			# inset = -0.10,
			# vert = TRUE,
			# width = 0.12,
			# height = 1,
			# labels = labs,
			# labAdj = 0.3,
			# cex = 0.62,
			# col = c('white', 'green4'),
			# border = 'black',
			# title = 'Occupancy',
			# titleAdj = c(0.5, 0.865),
			# adjX = c(0, 0.3),
			# adjY = c(0.115, 0.72),
			# boxBg = par('bg'),
			# boxBorder = NULL,
			# swatches = list(list(swatchAdjY=c(0.77, 0.815), col='darkgreen', border='black', labels='Range\ncore'))
		# )	

		# title(sub=date(), cex.sub=0.35, line=-1, outer=TRUE)
		
	# dev.off()
	
# say('##################')
# say('### model SNPs ###')
# say('##################')	
	
	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	
	# # output
	# dirCreate('./Models - SNPs')
	
	# # SNP names
	# snps <- names(phenoGeno)
	# snps <- snps[grepl(pattern='snp', x=snps)]
	
	# # make formulae
	# formulae <- paste('resp ~', paste(predictors, collapse= ' + '))
	# formulae <- as.formula(formulae)
	# formulae <- makeFormulae(formulae, quad=FALSE, ia=TRUE, returnFx=as.character)

	# # for each SNP
	# for (countSnp in seq_along(snps)) {
	# # for (countSnp in 1:1000) {
	
		# snp <- snps[countSnp]
	
		# say('Modeling ', snp, ' (model ', countSnp, ' of ', length(snps), ') on ', date())
	
		# ### collate data
		
		# # wanting one column for number of reference allele copies and one for alternative allele copies
		
		# snpData <- phenoGeno[ , c('population', ll, predictors, 'pc1', 'pc2', snp)]

		# start <- regexpr(snp, pattern='ref')[1]
		# refBase <- substr(snp, start + 4, start + 4)
		# altBase <- substr(snp, nchar(snp), nchar(snp))
	
		# snpCode <- snpData[ , snp]
		# refAllele <- numGeno <- rep(NA, length(snpCode))
		
		# ### tally number of alleles in each population
		# for (i in seq_along(refAllele)) {
		
			# if (is.na(snpCode[i])) {
				# refAllele[i] <- NA
				# numGeno[i] <- 0
			# } else if (snpCode[i] == -1) {
				# refAllele[i] <- NA
				# numGeno[i] <- 0
			# } else if (snpCode[i] == 0) {
				# refAllele[i] <- 2
				# numGeno[i] <- 1
			# } else if (snpCode[i] == 1) {
				# refAllele[i] <- 1
				# numGeno[i] <- 1
			# } else if (snpCode[i] == 2) {
				# refAllele[i] <- 0
				# numGeno[i] <- 1
			# }
			
		# }
		
		# snpData$refAllele <- refAllele
		# snpData$altAllele <- 2 - refAllele
		# snpData$numGeno <- numGeno
		
		# snpCounts <- stats::aggregate(x=snpData[ , c('numGeno', 'refAllele', 'altAllele')], by=list(snpData$population), FUN=sum, na.rm=TRUE)
		# snpMean <- aggregate(snpData, by=list(snpData$population), mean, na.rm=TRUE)
		
		# snpMean$population <- NULL
		# names(snpMean)[1] <- 'population'

		# ### merge data
		# snpByPop <- snpMean
		# snpByPop[ , snp] <- NULL
		# snpByPop$numGeno <- snpCounts$numGeno
		# snpByPop$refAllele <- snpCounts$refAllele
		# snpByPop$altAllele <- snpCounts$altAllele
		
		# if (any(snpByPop$numGeno == 0)) snpByPop <- snpByPop[-which(snpByPop$numGeno == 0), ]
	
		# # scale to 0 mean and unit sd
		# scales <- scale(snpByPop[ , predictors])
		# snpByPopScaled <- snpByPop
		# snpByPopScaled[ , predictors] <- scale(snpByPopScaled[ , predictors])
	
		# ### full model
		# resp <- as.matrix(snpByPop[ , c('refAllele', 'altAllele')])

		# allSitesModels <- data.frame()
		# for (i in seq_along(formulae)) {
		
			# thisForm <- as.formula(formulae[i])
			# thisModel <- glm(thisForm, data=snpByPopScaled, weights=snpByPop$numGeno, family='binomial')
			# AICc <- AICc(thisModel)
			
			# allSitesModels <- rbind(
				# allSitesModels,
				# data.frame(
					# formula = formulae[i],
					# AICc = AICc
				# )
			# )
			
		# }
		
		# allSitesModels <- allSitesModels[order(allSitesModels$AICc), ]
		
		# allSitesModelFormula <- as.formula(allSitesModels$formula[1])
		# allSitesModel <- glm(allSitesModelFormula, data=snpByPopScaled, weights=snpByPop$numGeno, family='binomial', na.action='na.fail')
		
		# ### LOO cross-validation
		# looModels <- list()
		# for (countPop in seq_along(snpByPop$population)) {
		
			# pop <- snpByPopScaled$population[countPop]
		
			# looSnpByPopScaled <- snpByPopScaled[(snpByPopScaled$population != pop), ]
			# resp <- as.matrix(looSnpByPopScaled[ , c('refAllele', 'altAllele')])
			# looModel <- glm(allSitesModelFormula, data=looSnpByPopScaled, weights=looSnpByPopScaled$numGeno, family='binomial', na.action='na.fail')
			
			# looModels[[countPop]] <- list()
			# names(looModels)[countPop] <- pop
			# looModels[[countPop]]$looPop <- pop
			# looModels[[countPop]]$looModel <- looModel
			
		# }
		
		# ### remember
		# snpModel <- list()
		# snpModel$snp <- snp
		# snpModel$refBase <- refBase
		# snpModel$altBase <- altBase
		# snpModel$snpByPop <- snpByPop
		# snpModel$predictorMeans <- attributes(scales)$`scaled:center`
		# snpModel$predictorSds <- attributes(scales)$`scaled:scale`
		# snpModel$allSitesModelFormula <- allSitesModelFormula
		# snpModel$allSitesModels <- allSitesModels
		# snpModel$allSitesModel <- allSitesModel
		# snpModel$looModels <- looModels
		
		# save(snpModel, file=paste0('./Models - SNPs/Model for SNP ', snp, '.RData'))

	# }
	
# say('##################################################################')
# say('### assess SNP model performance and sensitivity to predictors ###')
# say('##################################################################')	

	# # generalization
	# iters <- 1000 # times to iterate permutation test

	# dirCreate('./Figures & Tables/SNP Models')

	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	
	# # SNP names
	# snps <- names(phenoGeno)
	# snps <- snps[grepl(pattern='snp', x=snps)]

	# # output
	# perform <- data.frame()

	# # by SNP
	# for (countSnp in seq_along(snps)) {
	# # for (countSnp in 1:30) {
	
		# snp <- snps[countSnp]
		
		# say(countSnp, ' ', snp)
		
		# load(paste0('./Models - SNPs/Model for SNP ', snp, '.RData'))
		
		# snpByPop <- snpModel$snpByPop

		# predictorMeans <- snpModel$predictorMeans
		# predictorSds <- snpModel$predictorSds
		
		# ### pseudo R2
		# #############
		
		# dev <- snpModel$allSitesModel$deviance
		# nullDev <- snpModel$allSitesModel$null.deviance
		# pseudoR2 <- (nullDev - dev) / nullDev

		# ### variation in alleles across populations
		# ###########################################
		
		# snpByPop <- snpModel$snpByPop
		# alleleFreq <- snpByPop$refAllele / rowSums(snpByPop[ , c('refAllele', 'altAllele')])
		# sdAlleleFreq <- sd(alleleFreq)
	
		# ### variable importance
		# #######################
		
		# ### predict to all-sites
		# allSitesModel <- snpModel$allSitesModel
		# predictionAllSites <- fitted(allSitesModel, type='response')
	
		# ### create predictor matrix
		# invioableRandPreds <- snpByPop[ , predictors]
		
		# # add column for intercept
		# intercept <- data.frame(intercept = rep(1, nrow(snpByPop)))
		# invioableRandPreds <- insertCol(intercept, into=invioableRandPreds, at=1)

		# # rescale
		# for (pred in predictors) {
			# invioableRandPreds[ , pred] <- (invioableRandPreds[ , pred] - predictorMeans[[pred]]) / predictorSds[[pred]]
		# }
		
		# # add 2-way interaction terms
		# for (countPred1 in 1:(length(predictors) - 1)) {
			# for (countPred2 in 2:length(predictors)) {
				
				# pred1 <- predictors[countPred1]
				# pred2 <- predictors[countPred2]
				
				# ia <- invioableRandPreds[ , pred1] * invioableRandPreds[ , pred2]
				# ia <- data.frame(ia=ia)
				# names(ia) <- paste0(pred1, 'X', pred2)
				# invioableRandPreds <- cbind(invioableRandPreds, ia)
			
			# }
		# }
		
		# # get columns for predictors that are in model
		# coeffsInModel <- sprintf('%.48f', coef(allSitesModel))
		# coeffsInModel <- as.numeric(coeffsInModel)
		# coeffInModelNames <- names(coefficients(allSitesModel))
		# names(coeffsInModel) <- coeffInModelNames
		# coeffInModelNames <- coeffInModelNames[-which(coeffInModelNames == '(Intercept)')]
		# coeffInModelNames <- gsub(coeffInModelNames, pattern=':', replacement='X')
		
		# randPreds <- invioableRandPreds[ , 'intercept', drop=FALSE]
		
		# for (thisCoeff in coeffInModelNames) {
		
			# addCol <- if (thisCoeff %in% names(invioableRandPreds)) {
				# invioableRandPreds[ , thisCoeff, drop = FALSE]
			# }
			
			# randPreds <- cbind(randPreds, addCol)
			
		# }

		# # to store importances
		# coeffNames <- getCoeffNames(predictors)
		# importances <- importancesLower <- importancesUpper <- rep(NA, length(coeffNames))
		# names(importances) <- names(importancesLower) <- names(importancesUpper) <- coeffNames

		# # test effect of each predictor
		# for (thisCoeff in coeffNames) {
		
			# thisRandPreds <- randPreds

			# # coefficient isn't in model
			# if (!(thisCoeff %in% coeffInModelNames)) {
			
				# corImp <- 0
				
			# # model with more than just intercept
			# } else {
			
				# corImp <- rep(NA, iters)
			
				# for (iter in 1:iters) {
				
					# thisRandPreds[ , thisCoeff] <- sample(thisRandPreds[ , thisCoeff], nrow(thisRandPreds))

					# predictionAllSitesRand <- c(as.matrix(thisRandPreds) %*% coeffsInModel)
					# predictionAllSitesRand <- allSitesModel$family$linkinv(predictionAllSitesRand)

					# thisCorImp <- cor(predictionAllSites, predictionAllSitesRand)
					# corImp[iter] <- if (is.na(thisCorImp)) {
						# 0
					# } else {
						# 1 - thisCorImp
					# }
				
				# }
				
			# }
			
			# importances[[thisCoeff]] <- mean(corImp)
			# importancesLower[[thisCoeff]] <- quantile(corImp, 0.025)
			# importancesUpper[[thisCoeff]] <- quantile(corImp, 0.975)
			
		# } # next predictor
	
		# names(importancesLower) <- paste0(names(importancesLower), '_0pt025')
		# names(importancesUpper) <- paste0(names(importancesUpper), '_0pt975')

		# importances <- c(importances, importancesLower, importancesUpper)
		# importances <- importances[order(names(importances))]
		# importances <- t(as.data.frame(importances))
		# colnames(importances) <- paste0('import_', colnames(importances))
	
		# ### LOO cross-validation
		# numPops <- length(snpModel$looModels)
		# looDiffs <- rep(NA, numPops)
		# names(looDiffs) <- names(snpModel$looModels)
		
		# for (countPop in 1:numPops) {
		
			# pop <- names(snpModel$looModels)[countPop]
			
			# ### predict to LOO population
			# looPop <- snpByPop[snpByPop$population == pop, ]
			
			# # rescale
			# for (pred in predictors) {
				# looPop[ , pred] <- (looPop[ , pred] - predictorMeans[[pred]]) / predictorSds[[pred]]
			# }

			# looModel <- snpModel$looModels[[pop]]$looModel
			# predictionLoo <- predict(looModel, newdata=looPop, type='response')

			# obsFreq <- looPop$refAllele / (looPop$refAllele + looPop$altAllele)
			# looDiffs[countPop] <- predictionLoo - obsFreq
			
		# }
		
		# meanLooDiff <- mean(looDiffs)
		# looDiff_0pt025 <- quantile(looDiffs, 0.025)
		# looDiff_0pt975 <- quantile(looDiffs, 0.975)
		# meanAbsLooDiff <- mean(abs(looDiffs))
		# meanAbsLooDiff_0pt025 <- quantile(abs(looDiffs), 0.025)
		# meanAbsLooDiff_0pt975 <- quantile(abs(looDiffs), 0.975)
			
		# ### observed frequencies
		# meanRefAlleleFreq <- sum(snpByPop$refAllele) / (sum(snpByPop$refAllele) + sum(snpByPop$altAllele))
			
		# ### remember
		# perform <- rbind(
			# perform,
			# cbind(
				# data.frame(
					# snp=snp,
					# refBase=snpModel$refBase,
					# altBase=snpModel$altBase,
					# meanRefAlleleFreq=meanRefAlleleFreq,
					# sdAlleleFreq=sdAlleleFreq,
					# numPops=nrow(snpModel$snpByPop),
					# numIndivid=sum(snpByPop$numGeno),
					# pseudoR2=pseudoR2,
					# looDiff_0pt025=looDiff_0pt025,
					# meanLooDiff=meanLooDiff,
					# looDiff_0pt975=looDiff_0pt975,
					# meanAbsLooDiff_0pt025=meanAbsLooDiff_0pt025,
					# meanAbsLooDiff=meanAbsLooDiff,
					# meanAbsLooDiff_0pt975=meanAbsLooDiff_0pt975
				# ),
				# importances
			# )
		# )
		
		# # print(perform)
	# }

	# rownames(perform) <- 1:nrow(perform)
	# write.csv(perform, './Figures & Tables/SNP Models - Summaries/SNP Model Performance and Predictor Importance.csv', row.names=FALSE)
	
# say('#############################################################')
# say('### EDA of SNP model performance and predictor importance ###')
# say('#############################################################')

	# perform <- read.csv('./Figures & Tables/SNP Models - Summaries/SNP Model Performance and Predictor Importance.csv')

	# png('./Figures & Tables/SNP Models - Summaries/SNP Models - EDA.png', width=2000, height=1200, res=200)
		
		# par(mfrow=c(2, 4))
		
		# x <- perform$sdAlleleFreq
		# x <- na.omit(x)
		# x <- x[!is.infinite(x)]
		# hist(x, breaks=seq(0, 1, by=0.025), xlab='SD', main='Standard deviation of\nallelic frequencies\nacross populations')
		# abline(v=rep(mean(x, na.rm=TRUE), 2), lwd=2, col='firebrick1')
		
		# x <- perform$pseudoR2
		# x <- na.omit(x)
		# x <- x[!is.infinite(x)]
		# hist(x, breaks=seq(0, 1, by=0.025), xlab='Pseudo R2', main='Pseudo R2')
		# abline(v=rep(mean(x, na.rm=TRUE), 2), lwd=2, col='firebrick1')
		
		# x <- perform$meanLooDiff
		# x <- na.omit(x)
		# x <- x[!is.infinite(x)]
		# hist(x, breaks=40, xlab='Mean LOO Difference', main='LOO DIfference')
		# abline(v=rep(mean(x, na.rm=TRUE), 2), lwd=2, col='firebrick1')
		
		# x <- perform$meanAbsLooDiff
		# x <- na.omit(x)
		# x <- x[!is.infinite(x)]
		# hist(x, breaks=40, xlab='Mean Absolute LOO Difference', main='Absolute LOO DIfference')
		# abline(v=rep(mean(x, na.rm=TRUE), 2), lwd=2, col='firebrick1')
		
		# coeffNames <- getCoeffNames(predictors)
		# imports <- perform[ , paste0('import_', coeffNames)]
	
		# step <- 0.05
		# xmax <- max(as.matrix(imports))
		# xmax <- step * ceiling(xmax / step)
		
		# breaks <- seq(0, xmax, by=step)

		# for (thisCoeff in coeffNames) {
		
			# x <- perform[ , paste0('import_', thisCoeff)]
			# x <- na.omit(x)
			# nice <- predNice(thisCoeff)$nice
			# name <- predNice(thisCoeff)$name
			# hist(x, breaks=breaks, xlab=nice, main=paste0('Importance of\n', name), xlim=c(0, xmax))
			# abline(v=rep(mean(x), 2), lwd=2, col='firebrick1')
			
		# }

	# dev.off()

# say('########################################################')
# say('### identify alleles most likely affected by climate ###')
# say('########################################################')
	
	# # generalization
	# thold <- 0.50 # get alleles that fall within best "thold" percentile
	
	# # load model performance data
	# perform <- read.csv('./Figures & Tables/SNP Models - Summaries/SNP Model Performance and Predictor Importance.csv')
	
	# # get alleles most strongly associated with climate
	# bestAlleles <- perform[
		# perform$sdAlleleFreq >= 0.1 &
		# perform$pseudoR2 >= 0.5 &
		# perform$meanLooDiff >= quantile(perform$meanLooDiff, 0.5 - 0.5 * thold) &
		# perform$meanLooDiff <= quantile(perform$meanLooDiff, 0.5 + 0.5 * thold) &
		# perform$meanAbsLooDiff <= quantile(perform$meanAbsLooDiff, thold) & 
		# (perform$import_bio01 >= 1 - thold |
		# perform$import_bio12 >= 1 - thold |
		# perform$import_bio01Xbio12 >= 1 - thold),
	# ]

	# sink('./Figures & Tables/SNP Models - Summaries/SNPs Most Sensitive to Climate.txt', split=TRUE)
		
		# # summarize
		# say('Using threshold of ', thold, ' to define "best" models.')
		
		# say('There are ', nrow(bestAlleles), ' alleles most likely responding to climate.')

		# n <- sum(bestAlleles$import_bio01 >= 1 - thold & bestAlleles$import_bio12 < 1 - thold & bestAlleles$import_bio01Xbio12 < 0.75)
		# say(n , ' alleles sensitive to just BIO01.')
		
		# n <- sum(bestAlleles$import_bio01 < 1 - thold & bestAlleles$import_bio12 >= 1 - thold & bestAlleles$import_bio01Xbio12 < 0.75)
		# say(n , ' alleles sensitive to just BIO12.')
		
		# n <- sum(bestAlleles$import_bio01 >= 1 - thold & bestAlleles$import_bio12 >= 1 - thold & bestAlleles$import_bio01Xbio12 < 0.75)
		# say(n , ' alleles sensitive to BIO01 and BIO12 (but not interaction).')
		
		# n <- sum(bestAlleles$import_bio01 < 1 - thold & bestAlleles$import_bio12 < 1 - thold & bestAlleles$import_bio01Xbio12 >= 0.75)
		# say(n , ' alleles sensitive to interaction between BIO01 and BIO12.')
		
	# sink()
		
	# write.csv(bestAlleles, './Figures & Tables/SNP Models - Summaries/SNPs Most Sensitive to Climate.csv')

# say('##############################################################')
# say('### make maps of alleles most likely responding to climate ###')
# say('##############################################################')

	# # generalization
	# rcp <- '8pt5' # RCP
	
	# bestAlleles <- read.csv('./Figures & Tables/SNP Models - Summaries/SNPs Most Sensitive to Climate.csv')

	# # load geno/pheno and GIS data
	# load('./Study Region/North America Level 1 Sans Alaska.RData')
	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	# sqEnv <- stackCurrentEnv(predictors)
	# futEnv <- stackFutureEnv(predictors, rcp=rcp)
	# load('./Models - Species/bayesLopod/Prediction - bayesLopod with Variable Detection with Future Predictions.Rdata')

	# # aggregate geno/pheno data to population
	# phenoGenoPops <- aggPhenoGenoByPop
	# focus <- getPlotFocus(buffSize=buffSize)
	
	# sqEnv <- crop(sqEnv, focus)
	# futEnv <- lapply(futEnv, crop, y=focus)
	# nam1 <- crop(nam1, focus)
	# lopodPredictVarP <- crop(lopodPredictVarP, focus)
	
	# # colors
	# cols <- cividis(101)
	
	# cols <- colorRampPalette(c('darkgreen', 'green', 'purple', 'red', 'darkred'))
	# cols <- colorRampPalette(c('chartreuse', 'purple', 'firebrick3'))
	# cols <- colorRampPalette(c('darkred', 'goldenrod', 'navyblue'))
	# cols <- cols(101)
		
	# maxUncert <- sd(c(rep(100, 4), rep(0, 5)))
	# maxUncert <- 5 * round(maxUncert / 5)
		
	# colsUncert <- colorRampPalette(c('gray80', 'red'))
	# colsUncert <- colsUncert(maxUncert)

	# dirCreate('./Figures & Tables/SNP Models - Maps')

	# ### by SNP
	# for (countSnp in 1:nrow(bestAlleles)) {
	# # for (countSnp in 1) {
	
		# snp <- bestAlleles$snp[countSnp]
		# say(snp)
	
		# # load model
		# load(paste0('./Models - SNPs/Model for SNP ', snp, '.RData'))
	
		# # rescale
		# predictorMeans <- snpModel$predictorMeans
		# predictorSds <- snpModel$predictorSds
		
		# thisSqEnv <- sqEnv
		# layerNames <- names(thisSqEnv)
		# for (pred in predictors) {
			# thisSqEnv[[pred]] <- (thisSqEnv[[pred]] - predictorMeans[[pred]]) / predictorSds[[pred]]
		# }
		# names(thisSqEnv) <- layerNames
	
		# thisFutEnv <- futEnv
		# for (gcm in names(thisFutEnv)) {
			# for (pred in predictors) {
				# thisFutEnv[[gcm]][[pred]] <- (thisFutEnv[[gcm]][[pred]] - predictorMeans[[pred]]) / predictorSds[[pred]]
			# }
			# names(thisFutEnv[[gcm]]) <- layerNames
		# }
	
		# ### current predict
		# sqPrediction <- predict(thisSqEnv, snpModel$allSitesModel, type='response')
		# sqPrediction <- round(100 * sqPrediction)
		# minPredSq <- minValue(sqPrediction)
		# maxPredSq <- maxValue(sqPrediction)
		# theseColsSq <- cols[(1 + minPredSq):(1 + maxPredSq)]
	
		# ### future predict
		# futPredictionList <- lapply(thisFutEnv, predict, model=snpModel$allSitesModel, type='response')
		# futPrediction <- futPredictionList[[1]]
		# for (i in 2:length(futPredictionList)) futPrediction <- stack(futPrediction, futPredictionList[[i]])
		# names(futPrediction) <- names(futPredictionList)
		
		# futPredictionMean <- mean(futPrediction)
		# futPredictionMean <- round(100 * futPredictionMean)
		
		# minPredFut <- min(minValue(futPredictionMean))
		# maxPredFut <- max(maxValue(futPredictionMean))
		# theseColsFut <- cols[(1 + minPredFut):(1 + maxPredFut)]
	
		# # uncertainty around future
		# futUncertain <- calc(futPrediction, sd)
		# futUncertain <- round(100 * futUncertain)
		# maxPredUncert <- cellStats(futUncertain, 'max')
		# theseColsUncert <- colsUncert[1:(1 + maxPredUncert)]
	
		# # get allelic frequencies
		# snpByPop <- snpModel$snpByPop
		# freqRef <- snpByPop$refAllele / rowSums(snpByPop[ , c('refAllele', 'altAllele')])
		# freqAlt <- 1 - freqRef
	
		# # plot
		# png(paste0('./Figures & Tables/SNP Models - Maps/', snp, '.png'), width=3 * 1200, height=1000, res=300)

			# par(mfrow=c(1, 3), oma=0.1 * c(4, 8, 1, 1), mai=0.1 * c(1, 1, 1, 1))

			# ### present
			# ###########
			
			# plot(sqPrediction, col=theseColsSq, breaks=minPredSq:maxPredSq, legend=FALSE, axes=FALSE, box=FALSE)
			# plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_i), border=NA, add=TRUE)
			# plot(nam1, add=TRUE)
			
			# # add populations
			# for (countPop in 1:nrow(snpByPop)) {

				# freqRef <- snpByPop$refAllele[countPop] / (snpByPop$refAllele[countPop] + snpByPop$altAllele[countPop])
				# popCol <- cols[round(100 * freqRef)]
				# points(x=snpByPop[countPop, ll[1]], y=snpByPop[countPop, ll[2]], pch=21, bg=popCol, cex=1.2)
				
			# }
			
			# # title
			# usr <- par('usr')
			# x <- usr[1]
			# y <- usr[4] - 0.02 * (usr[4] - usr[3])
			# underscore <- regexpr(snp, pattern='_')
			# snpName <- substr(snp, 5, underscore[1] - 1)
			# snpName <- paste0('a) SNP ', snpName, ': Current Climate')
			# text(x, y, labels=snpName, xpd=NA, cex=1.1, pos=4, font=2)
			
			# # legend
			# refBase <- snpModel$refBase
			# altBase <- snpModel$altBase
		
			# labs <- c(
				# paste0('100% ', refBase),
				# paste0('50% ', refBase, ' / 50% ', altBase),
				# paste0('100% ', altBase)
			# )
		
			# legendGrad(
				# x='bottom',
				# y = NULL,
				# inset = -0.1,
				# vert = FALSE,
				# width = 0.999,
				# height = 0.15,
				# labels = labs,
				# labAdj = -0.33,
				# cex = 0.9,
				# col = rev(cols),
				# border = 'black',
				# title = '',
				# titleAdj = c(0.5, 0.9),
				# adjX = c(0, 1),
				# adjY = c(0.8, 1),
				# boxBg = par('bg'),
				# boxBorder = NULL,
				# swatches = NULL
			# )	

			# # stats
			# allele <- perform[perform$snp == snp, ]

			# lab <- paste0('Pseudo-R2: ', sprintf('%.2f', allele$pseudoR2), '\n')
			# lab <- paste0(lab, 'MAD (LOO): ', sprintf('%.2f', allele$meanAbsLooDiff), '\n')
			
			# imports <- names(allele)[which(grepl(names(allele), pattern='import_'))]
			# importsMean <- imports[-which(grepl(imports, pattern='_0pt025') | grepl(imports, pattern='_0pt975'))]
			# importsLower <- imports[which(grepl(imports, pattern='_0pt025'))]
			# importsUpper <- imports[which(grepl(imports, pattern='_0pt975'))]
			
			# for (i in seq_along(importsMean)) {
			
				# thisImport <- importsMean[i]
				# thisImportSimple <- gsub(thisImport, pattern='import_', replacement='')
				# nice <- predNice(thisImportSimple)$nice
				
				# meanVal <- sprintf('%.2f', allele[ , thisImport])
				# lowVal <- sprintf('%.2f', allele[ , paste0(thisImport, '_0pt025')])
				# highVal <- sprintf('%.2f', allele[ , paste0(thisImport, '_0pt975')])
				
				# thisLab <- paste0('Importance of ', nice, ': ', meanVal, ' (', lowVal, ' - ', highVal, ')\n')
				
				# lab <- paste0(lab, thisLab)
				
			# }
			
			# text(usr[1], usr[3] + 0.125 * (usr[4] - usr[3]), labels=lab, xpd=NA, cex=0.67, pos=4, font=2)
			
			# ### future
			# ###########
			
			# plot(futPredictionMean, col=theseColsFut, breaks=minPredFut:maxPredFut, legend=FALSE, axes=FALSE, box=FALSE)
			# plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_ensembleGcm_rcp85_2070s), border=NA, add=TRUE)
			# plot(nam1, add=TRUE)
			
			# # add populations
			# for (countPop in 1:nrow(snpByPop)) {

				# freqRef <- snpByPop$refAllele[countPop] / (snpByPop$refAllele[countPop] + snpByPop$altAllele[countPop])
				# popCol <- cols[round(100 * freqRef)]
				# points(x=snpByPop[countPop, ll[1]], y=snpByPop[countPop, ll[2]], pch=21, bg=popCol, cex=1.2)
				
			# }
			
			# # title
			# usr <- par('usr')
			# x <- usr[1]
			# y <- usr[4] - 0.02 * (usr[4] - usr[3])
			# underscore <- regexpr(snp, pattern='_')
			# snpName <- substr(snp, 5, underscore[1] - 1)
			# snpName <- paste0('b) Mean prediction for SNP ', snpName, ' under RCP8.5 for 2070s')
			# text(x, y, labels=snpName, xpd=NA, cex=1.1, pos=4, font=2)
			
			# # legend
			# refBase <- snpModel$refBase
			# altBase <- snpModel$altBase
		
			# labs <- c(
				# paste0('100% ', refBase),
				# paste0('50% ', refBase, ' / 50% ', altBase),
				# paste0('100% ', altBase)
			# )
		
			# legendGrad(
				# x='bottom',
				# y = NULL,
				# inset = -0.1,
				# vert = FALSE,
				# width = 0.999,
				# height = 0.15,
				# labels = labs,
				# labAdj = -0.33,
				# cex = 0.9,
				# col = rev(cols),
				# border = 'black',
				# title = '',
				# titleAdj = c(0.5, 0.9),
				# adjX = c(0, 1),
				# adjY = c(0.8, 1),
				# boxBg = par('bg'),
				# boxBorder = NULL,
				# swatches = NULL
			# )	

			# ### uncertainty
			# ###############
			
			# plot(futUncertain, col=theseColsUncert, breaks=0:maxPredUncert, legend=FALSE, axes=FALSE, box=FALSE)
			# plot(nam1, add=TRUE)
			
			# # add populations
			# points(x=snpByPop[ , ll[1]], y=snpByPop[ , ll[2]], pch=1, cex=1.2)
			
			# # title
			# usr <- par('usr')
			# x <- usr[1]
			# y <- usr[4] - 0.02 * (usr[4] - usr[3])
			# underscore <- regexpr(snp, pattern='_')
			# snpName <- paste('c) Uncertainty in allelic frequencies under RCP8.5 for 2070s')
			# text(x, y, labels=snpName, xpd=NA, cex=1.1, pos=4, font=2)
			
			# # legend
			# labs <- paste0(c(0, round(maxUncert / 2, 1), maxUncert), '%')
		
			# legendGrad(
				# x='bottom',
				# y = NULL,
				# inset = -0.1,
				# vert = FALSE,
				# width = 0.999,
				# height = 0.15,
				# labels = labs,
				# labAdj = -0.33,
				# cex = 0.9,
				# col = colsUncert,
				# border = 'black',
				# title = '',
				# titleAdj = c(0.5, 0.9),
				# adjX = c(0, 1),
				# adjY = c(0.8, 1),
				# boxBg = par('bg'),
				# boxBorder = NULL,
				# swatches = NULL
			# )	

			# title(sub=date(), cex.sub=0.5, line=-0.5, outer=TRUE)
			
		# dev.off()

	# } # next allelle

# say('#####################################')
# say('### make map of allelic diversity ###')
# say('#####################################')

	# # generalization
	# rcp <- '8pt5' # RCP
	# thold <- 0.9 # topmost quantile of psi above which is designated as range "core"
	
	# bestAlleles <- read.csv('./Figures & Tables/SNP Models - Summaries/SNPs Most Sensitive to Climate.csv')

	# # load geno/pheno and GIS data
	# load('./Study Region/North America Level 1 Sans Alaska.RData')
	# load('./Data/Phenotypic & Genotypic Data/02 Phenotypic & Genotypic Data - Environmental Data & PCA.RData')
	# sqEnv <- stackCurrentEnv(predictors)
	# futEnv <- stackFutureEnv(predictors, rcp=rcp)
	# load('./Models - Species/bayesLopod/Prediction - bayesLopod with Variable Detection with Future Predictions.Rdata')

	# # aggregate geno/pheno data to population
	# phenoGenoPops <- aggPhenoGenoByPop(phenoGeno, sp=TRUE)

	# # crop to focal region extent
	# focus <- getPlotFocus(buffSize=buffSize)
	# sqEnv <- crop(sqEnv, focus)
	# futEnv <- lapply(futEnv, crop, y=focus)
	# nam1 <- crop(nam1, focus)
	# lopodPredictVarP <- crop(lopodPredictVarP, focus)
	
	# ### by SNP
	# for (countSnp in 1:nrow(bestAlleles)) {
	# # for (countSnp in 1:5) {
	
		# snp <- bestAlleles$snp[countSnp]
		# say(snp)
	
		# # load model
		# load(paste0('./Models - SNPs/Model for SNP ', snp, '.RData'))
	
		# # rescale
		# predictorMeans <- snpModel$predictorMeans
		# predictorSds <- snpModel$predictorSds
		
		# thisSqEnv <- sqEnv
		# layerNames <- names(thisSqEnv)
		# for (pred in predictors) {
			# thisSqEnv[[pred]] <- (thisSqEnv[[pred]] - predictorMeans[[pred]]) / predictorSds[[pred]]
		# }
		# names(thisSqEnv) <- layerNames
	
		# thisFutEnv <- futEnv
		# for (gcm in names(thisFutEnv)) {
			# for (pred in predictors) {
				# thisFutEnv[[gcm]][[pred]] <- (thisFutEnv[[gcm]][[pred]] - predictorMeans[[pred]]) / predictorSds[[pred]]
			# }
			# names(thisFutEnv[[gcm]]) <- layerNames
		# }
	
		# ### current predict
		# thisSqPred <- predict(thisSqEnv, snpModel$allSitesModel, type='response')
		# thisSqDiversity <- thisSqPred * (1 - thisSqPred)
		# names(thisSqDiversity) <- snp
	
		# ### future predict
		# futPredictionList <- lapply(thisFutEnv, predict, model=snpModel$allSitesModel, type='response')
		# thisFutPred <- futPredictionList[[1]]
		# for (i in 2:length(futPredictionList)) thisFutPred <- stack(thisFutPred, futPredictionList[[i]])
		# names(thisFutPred) <- names(thisFutPred)
		# thisFutPred <- mean(thisFutPred)
		# thisFutDiversity <- thisFutPred * (1 - thisFutPred)
		# names(thisFutDiversity) <- snp
	
		# # get allelic frequencies
		# snpByPop <- snpModel$snpByPop
		# freqRef <- snpByPop$refAllele / rowSums(snpByPop[ , c('refAllele', 'altAllele')])
		# freqAlt <- 1 - freqRef
		# thisObsDiversity <- freqRef * freqAlt
		
		# if (countSnp == 1) {
			# sqDiversity <- thisSqDiversity
			# futDiversity <- thisFutDiversity
			# obsDiversity <- matrix(thisObsDiversity, ncol=1)
			# colnames(obsDiversity) <- snp
		# } else {
			# sqDiversity <- stack(sqDiversity, thisSqDiversity)
			# futDiversity <- stack(futDiversity, thisFutDiversity)
			# obsDiversity <- cbind(obsDiversity, matrix(thisObsDiversity, ncol=1))
			# colnames(obsDiversity)[countSnp] <- snp
		# }
	
	# } # next SNP
	
	# ### calculate diversity
	# #######################
	
	# sqDiversity <- sum(sqDiversity / 0.5^2) / nlayers(sqDiversity)
	# futDiversity <- sum(futDiversity / 0.5^2) / nlayers(futDiversity)
	
	# obsDiversity <- rowSums(obsDiversity / 0.5^2) / ncol(obsDiversity)
	
	# ### plot
	# ########
	
	# # colors
	# cols <- cividis(101)
		
	# dirCreate('./Figures & Tables/SNP Models - Genetic Diversity')

	# coreQuant <- quantile(lopodPredictVarP@data$psi_i, thold, na.rm=TRUE)
	
	# # plot
	# png(paste0('./Figures & Tables/SNP Models - Genetic Diversity/Genetic Diversity.png'), width=2 * 1200, height=1000, res=300)

		# par(mfrow=c(1, 2), oma=0.1 * c(1, 4, 1, 26), mai=0 * c(1, 1, 1, 1), mar=rep(0, 4))

		# ### present
		# ###########

		# # diversity
		# plot(sqDiversity, col=cols, breaks=seq(0, 1, by=0.01), legend=FALSE, axes=FALSE, box=FALSE)
		# plot(nam1, add=TRUE)

		# # range core
		# coreQuant <- quantile(lopodPredictVarP$psi_i, thold, na.rm=TRUE)
		# whichCore <- which(lopodPredictVarP@data$psi_i >= coreQuant)
		# rangeCore <- lopodPredictVarP[whichCore, ]
		# rangeCore <- gUnaryUnion(rangeCore)
		# plot(rangeCore, col=NA, border='white', lwd=1, add=TRUE)

		# # fade by probability of occurrence
		# plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_i), border=NA, add=TRUE)
		
		# # add populations
		# popCol <- cols[round(100 * obsDiversity)]
		# points(phenoGenoPops, pch=21, bg=popCol, cex=1.2)
		
		# # title
		# usr <- par('usr')
		# x <- usr[1]
		# y <- usr[4] - 0.1 * (usr[4] - usr[3])
		# text(x, y, labels='a) Current climate', xpd=NA, cex=0.7, pos=4, font=2)
		
		# ### future
		# ###########

		# # diversity
		# plot(futDiversity, col=cols, breaks=seq(0, 1, by=0.01), legend=FALSE, axes=FALSE, box=FALSE)
		# plot(nam1, add=TRUE)

		# # range core
		# whichCore <- which(lopodPredictVarP@data$psi_ensembleGcm_rcp85_2070s >= coreQuant)
		# rangeCore <- lopodPredictVarP[whichCore, ]
		# rangeCore <- gUnaryUnion(rangeCore)
		# plot(rangeCore, col=NA, border='white', lwd=1.5, add=TRUE)

		# # fade by probability of occurrence
		# plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_ensembleGcm_rcp85_2070s), border=NA, add=TRUE)
		
		# # add populations
		# popCol <- cols[round(100 * obsDiversity)]
		# points(phenoGenoPops, pch=21, bg=popCol, cex=1.2)

		# # title
		# usr <- par('usr')
		# x <- usr[1]
		# y <- usr[4] - 0.1 * (usr[4] - usr[3])
		# text(x, y, labels='b) 2070s RCP8.5', xpd=NA, cex=0.7, pos=4, font=2)
		
		# # legend
		# legendGrad(
			# x='right',
			# inset = -0.17,
			# vert = TRUE,
			# width = 0.12,
			# height = 1,
			# labels = c(0, 0.25, 0.5, 0.75, 1),
			# labAdj = 0.35,
			# cex = 0.65,
			# col = cols,
			# border = 'black',
			# title = 'Diversity',
			# titleAdj = c(0.5, 0.855),
			# adjX = c(0, 0.35),
			# # adjY = c(0.23, 0.775),
			# adjY = c(0.13, 0.775),
			# boxBg = par('bg'),
			# boxBorder = NULL
			# # swatches = list(list(swatchAdjY=c(0.13, 0.175), col=NA, border='orange', lwd=1.5, labels='Range\ncore'))
		# )	

		# title(sub=date(), cex.sub=0.35, line=-1, outer=TRUE)
		
	# dev.off()
	


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
# # say('##########################################')
# # say('### define LOO cross-validation blocks ###')
# # say('##########################################')

	# # say('I am wanting to define leave-one-out neighborhoods based on spatial autocorrelation of the predictors. In particular, I will use these neighborhoods to create LOO cross-validation blocks and to weight performance scores obtained from LOO cross-validation. The blocks will be defined based on a spatial variogram calculated across the region containing the training sites. I will use a region defined by a buffer applied to counties with populations with a width equal to the mean pairwise distance between sampled counties.', breaks=80)
	
	# # load('./Data/Andropogon Occurrences/!16 Andropogon and Poaceae Records with Sampled Populations.RData')
	
	# # # create buffer around sampled sites
	# # sampledAg <- ag[!is.na(ag$population), ]
	# # cents <- gCentroid(sampledAg, byid=TRUE)
	# # nnDists <- distm(cents)
	# # diag(nnDists) <- NA
	# # meanDist <- mean(nnDists, na.rm=TRUE)
	# # sampledAg <- sp::spTransform(sampledAg, getCRS('albersNA', TRUE))
	# # focalArea <- gBuffer(sampledAg, width=meanDist)
	# # focalArea <- sp::spTransform(focalArea, getCRS('wgs84', TRUE))
	# # agIn <- ag[focalArea, ]

	# # # get data for variogram
	# # cents <- gCentroid(agIn, byid=TRUE)
	# # cents <- SpatialPointsDataFrame(coordinates(cents), data=agIn@data[ , c('pc1', 'pc2')], proj4=getCRS('wgs84', TRUE))
	# # cents <- sp::spTransform(cents, getCRS('albersNA', TRUE))
	
	# # # fit variograms
	# # vario1 <- automap::autofitVariogram(pc1 ~ 1, cents)
	# # vario2 <- automap::autofitVariogram(pc2 ~ 1, cents)
	
	# # agIn <- sp::spTransform(agIn, getCRS('albersNA', TRUE))
	# # ak <- automap::autoKrige(pc1 ~ 1, cents, agIn)
	
	# # # varioMod1 <- afvmod(pc1 ~ 1, cents, miscFitOptions=list(equal.width.bins = TRUE))
	# # # varioMod2 <- afvmod(pc1 ~ 1, cents, miscFitOptions=list(equal.np.bins = TRUE))
	





















	
	
	
	
# say('##########################################')
# say('### copy and process raw response data ###')
# say('##########################################')

	# ### create "stage 00" copies of all raw data files
	# #################################################
		
		# # NEUTRAL genetic data (inc STRUCTURE -- OLD VERSION!!! using all alleles) by INDIVIDUAL
		# data <- read.xlsx('./Pheno-Geno-Biogeographic Data/RESTARTS/gbs2014_filterCS STRUCTURE analysis groupings for ADAM.xlsx')
		# names(data)[names(data) == 'pop'] <- 'site'
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv', row.names=FALSE)
		
		# # NEUTRAL genetic data (inc STRUCTURE -- REVISED!!! using non-candidate alleles) by INDIVIDUAL
		# data <- read.xlsx('C:/ecology/Dropbox/Andropogon/Modeling genotype-phenotype/gbs2014_REV-feb27-2018-filterCS_neutral allele frequency (mkuch2).xlsx', sheet='neutral genetic grouping K-4')
		# names(data)[names(data) == 'pop'] <- 'site'
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data REVISED [ORIG gbs2014_REV-feb27-2018-filterCS_neutral allele frequency (mkuch2)].csv', row.names=FALSE)
		
		# # PHENOTYPE by INDIVIDUAL (has *some* overlap in records with NEUTRAL genetic data)
		# data <- read.csv('./Phenotypic Analysis/Data/01 GBS_multi-pop_2014_Data 2016-07-11 By Individual.csv')
		# names(data)[names(data) == 'population'] <- 'site'
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Individual - Phenotype [ORIG GBS_multi-pop_2014_Data 2016-07-11].csv', row.names=FALSE)
	
		# # PHENOTYPE by POPULATION
		# data <- read.csv('./Phenotypic Analysis/Data/01 GBS_multi-pop_2014_Data 2015-07-08 By Population.csv')
		# data <- data[1:35, ]
		# names(data)[names(data) == 'population'] <- 'site'
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Population - Phenotype [ORIG GBS_multi-pop_2014_Data 2015-07-08].csv', row.names=FALSE)
	
		# # CANDIDATE ALLELIC FREQUENCIES by POPULATION
		# for (locus in loci) {
			# data <- read.xlsx('./Pheno-Geno-Biogeographic Data/RESTARTS/Candidate Allele Frequencies.xlsx', sheet=locus)
			# write.csv(data, paste0('C:/ecology/Drive/Research/Andropogon/Pheno-Geno-Biogeographic Data/00 By Population - Candidate Allele Frequencies for Locus ', locus, ' [ORIG Candidate Allele Frequencies].csv'), row.names=FALSE)
		# }
		
		# # PHENOTYPE and GENOTYPE by INDIVIDUAL
		# data <- read.xlsx('C:/ecology/Dropbox/Andropogon/2017-12-26 Figs for adam/gbs2014 phenotype and genotype samples.xlsx', sheet='Sheet1')
		# names(data)[names(data) == 'site'] <- 'siteNum'
		# names(data)[names(data) == 'pop'] <- 'site'
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Individual - IDs for ALL Individuals [ORIG gbs2014 phenotype and genotype samples].csv', row.names=FALSE)
		
		# data <- read.xlsx('C:/ecology/Dropbox/Andropogon/2017-12-26 Figs for adam/gbs2014 phenotype and genotype samples.xlsx', sheet='notes')
		# write.csv(data, './Pheno-Geno-Biogeographic Data/00 By Individual - IDs ALL Individuals Metadata [ORIG gbs2014 phenotype and genotype samples].csv', row.names=FALSE)
		
# say('############################################################')
# say('### create data frame with rows representing INDIVIDUALS ###')	
# say('############################################################')
	
	# # want one data frame for numGeno with fields re STRUCTURE-assigned population group, phenotypic measurements, vegetation class of habitat

	# ind <- read.csv('Pheno-Geno-Biogeographic Data/00 By Individual - IDs for ALL Individuals [ORIG gbs2014 phenotype and genotype samples].csv')
	# ind$gbs2014_phenotype_sort <- ind$gbs2014_filterCS_pop_sort <- ind$kuch <- ind$gbs_SEQUENCE <- NULL
	# names(ind)[names(ind) == 'gbs_name_2'] <- 'individual'
	
	# ind$site[ind$site == 'NE-4' | ind$site == 'NE-5'] <- 'NE-4/NE-5'
	# ind$site[ind$site == 'NE-3' | ind$site == 'NE-7'] <- 'NE-3/NE-7'

	# ### match numGeno to STRUCTURE groups
	# #########################################
		
		# # population-level genetic data
		# fromData <- read.csv('Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv')
		# fromData$individual <- paste0(substr(fromData$individual, 1, 2), '-', substr(fromData$individual, 3, 4), '-', substr(fromData$individual, 5, 6))
		
		# # are all numGeno in genetic data frame also in master numGeno data frame?
		# stopifnot(any(fromData$individual %in% ind$individual))
		
		# for (countGenGroup in 1:4) {
			
			# ind$DUMMY <- NA
			# names(ind)[ncol(ind)] <- paste0('structureGroup', countGenGroup)
			
			# for (countIndivid in 1:nrow(fromData)) {
			
				# individual <- fromData$individual[countIndivid]
			
				# ind[ind$individual == fromData$individual[countIndivid], paste0('structureGroup', countGenGroup)] <-
					# fromData[countIndivid, paste0('GeneticGroup', countGenGroup)]
			
			# }
		# }

	# ### match numGeno to PHENOTYPE (only numGeno that were also GENOTYPED)
	# ##############################################################################
		
		# # individual-level neutral genetic and phenotypic data
		# fromData <- read.csv('./Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv', as.is=TRUE)
		
		# fromData$individual <- paste0(substr(fromData$individual, 1, 2), '-', substr(fromData$individual, 3, 4), '-', substr(fromData$individual, 5, 6))
		
		# for (col in c('height', 'spad', 'leaf', 'veg_wgt', 'boot_wgt', 'seed_wgt', 'biomass')) {
		
			# newCol <- if (col == 'height') {
				# 'heightCm'
			# } else if (col == 'spad') {
				# 'spad'
			# } else if (col == 'leaf') {
				# 'leafWidthCm'
			# } else if (col == 'veg_wgt') {
				# 'biomassVegetativeG'
			# } else if (col == 'boot_wgt') {
				# 'biomassBootG'
			# } else if (col == 'seed_wgt') {
				# 'biomassSeedG'
			# } else if (col == 'biomass') {
				# 'biomassG'
			# }

			# ind$DUMMY <- NA
			# names(ind)[ncol(ind)] <- newCol
			
			# for (i in 1:nrow(fromData)) {
			
				# individual <- fromData$individual[i]
				# ind[ind$individual == individual, newCol] <- fromData[fromData$individual == individual, col]
			
			# }
			
		# }
	
	# ### match numGeno that were not genotyped to PHENOTYPE
	# ##########################################################
	
		# say('NOTE: As of yet there is no way to match numGeno that were not genotyped with their phenotype because there is no common ID between data sets!', level=2)
		
	# ### match numGeno to Kuchler vegetation groups
	# ##################################################
	
		# # individual-level phenotypic data
		# fromData <- read.csv('Pheno-Geno-Biogeographic Data/00 By Individual - Phenotype [ORIG GBS_multi-pop_2014_Data 2016-07-11].csv')
		
		# ind$kuchlerType <- ind$kuchlerForm <- NA # Form encompasses Type
		
		# for (site in unique(ind$site)) {
		
			# matchRow <- which(fromData$site == site)
			
			# stopifnot(same(fromData$Kuchler_Vegetation_type[matchRow]))
			# ind$kuchlerType[ind$site == site] <- fromData$Kuchler_Vegetation_type[matchRow[1]]
			# ind$kuchlerForm[ind$site == site] <- fromData$Kuchler_description[matchRow[1]]
			
		# }
		
	
	# write.csv(ind, 'Pheno-Geno-Biogeographic Data/01 By Individual - Pheno-Geno-Biogeographic Data.csv', row.names=FALSE)
		
# say('############################################################')
# say('### create data frame with rows representing POPULATIONS ###')
# say('############################################################')
	
	# # want one data frame for populations with fields re averaged STRUCTURE-assigned population group, candidate allelic frequencies, phenotypic measurements, vegetation class of habitat
	# pop <- read.csv('Pheno-Geno-Biogeographic Data/00 By Population - Phenotype [ORIG GBS_multi-pop_2014_Data 2015-07-08].csv')
	# names(pop)[which(names(pop) == 'population')] <- 'site'
	# pop$Annual_Precipitation_cm <- NULL
	
	# for (col in c('height_cm', 'jul_height_TIP', 'spad_avg', 'jul_spad_TIP', 'leaf_width_avg', 'jul_leaf_width_TIP', 'biomass', 'jul_biomass_TIP', 'boot_discrete', 'jul_boot_date', 'jul_boot_TIP', 'jul_transplant_date', 'height_cmLower', 'height_cmUpper', 'spad_avgLower', 'spad_avgUpper', 'leaf_width_avgLower', 'leaf_width_avgUpper', 'biomassLower', 'biomassUpper')) {
	
		# newCol <- if (col == 'height_cm') {
			# 'heightCm'
		# } else if (col == 'jul_height_TIP') {
			# 'heightJulianTIP'
		# } else if (col == 'spad_avg') {
			# 'spad'
		# } else if (col == 'jul_spad_TIP') {
			# 'spadJulianTIP'
		# } else if (col == 'leaf_width_avg') {
			# 'leafWidthCm'
		# } else if (col == 'jul_leaf_width_TIP') {
			# 'leafWidthJulianTIP'
		# } else if (col == 'biomass') {
			# 'biomassG'
		# } else if (col == 'jul_biomass_TIP') {
			# 'biomassJulianTIP'
		# } else if (col == 'boot_discrete') {
			# 'booted'
		# } else if (col == 'jul_boot_date') {
			# 'bootDateJulian'
		# } else if (col == 'jul_boot_TIP') {
			# 'bootJulianTIP'
		# } else if (col == 'jul_transplant_date') {
			# 'transplantJulian'
		# } else if (col == 'height_cmLower') {
			# 'heightCmLower'
		# } else if (col == 'height_cmUpper') {
			# 'heightCmUpper'
		# } else if (col == 'spad_avgLower') {
			# 'spadLower'
		# } else if (col == 'spad_avgUpper') {
			# 'spadUpper'
		# } else if (col == 'leaf_width_avgLower') {
			# 'leafWidthCmLower'
		# } else if (col == 'leaf_width_avgUpper') {
			# 'leafWidthCmUpper'
		# } else if (col == 'biomassLower') {
			# 'biomassGLower'
		# } else if (col == 'biomassUpper') {
			# 'biomassGUpper'
		# }
		
		# names(pop)[names(pop) == col] <- newCol
		
	# }

	# ### match populations with *MEAN* STRUCTURE-assigned population probabilities
	# #############################################################################
		
		# fromData <- read.csv('./Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv')
		# fromData$group <- fromData$individual <- fromData$individual.1 <- NULL
		# fromData <- aggregate(fromData, by=list(fromData$site), mean, na.rm=TRUE)
		# fromData$site <- NULL
		# names(fromData)[1] <- 'site'

		# fromData$site[fromData$site == 'IL-02'] <- 'IL-2'
		# fromData$site[fromData$site == 'NE-4' | fromData$site == 'NE-5'] <- 'NE-4/NE-5'
		# fromData$site[fromData$site == 'NE-3' | fromData$site == 'NE-7'] <- 'NE-3/NE-7'
		
		# pop$structureGroup4 <- pop$structureGroup3 <- pop$structureGroup2 <- pop$structureGroup1 <- NA
		
		# for (countPop in 1:nrow(fromData)) {
		
			# fromSite <- fromData$site[countPop]
			# pop$structureGroup1[pop$site == fromSite] <- fromData$GeneticGroup1[countPop]
			# pop$structureGroup2[pop$site == fromSite] <- fromData$GeneticGroup2[countPop]
			# pop$structureGroup3[pop$site == fromSite] <- fromData$GeneticGroup3[countPop]
			# pop$structureGroup4[pop$site == fromSite] <- fromData$GeneticGroup4[countPop]
			
		# }

	# ### match populations with candidate allelic frequency data
	# ###########################################################
	
		# for (locus in loci) {

			# # construct coherent data frame for this locus
			# origFreqs <- read.csv(paste0('C:/ecology/Drive/Research/Andropogon/Pheno-Geno-Biogeographic Data/00 By Population - Candidate Allele Frequencies for Locus ', locus, ' [ORIG Candidate Allele Frequencies].csv'))
			
			# fromData <- data.frame(
				# site=c(origFreqs$Stratum1, origFreqs$Stratum2),
				# allele=c(origFreqs$Allele1, origFreqs$Allele2),
				# freq=c(origFreqs$Frequency1, origFreqs$Frequency2)
			# )
			
			# if (anyNA(fromData$site)) fromData <- fromData[!is.na(fromData$site), ]

			# stopifnot(length(unique(na.omit(fromData$allele))) == 2)
			# alleles <- sort(unique(na.omit(fromData$allele))) 
			# allele <- alleles[1]
			# altAllele <- alleles[2]
			
			# # if allele was completely missing from a population then it does not appear as a line in the data
			# # add new line with 0 frequency
			# for (site in unique(fromData$site)) {
				
				# rowsWithSiteData <- which(fromData$site == site)
				
				# if (length(rowsWithSiteData) == 1) {
				
					# missingAllele <- if (fromData$allele[rowsWithSiteData] == allele) {
						# altAllele
					# } else {
						# allele
					# }
				
					# fromData <- rbind(
						# fromData,
						# data.frame(
							# site=site,
							# allele=missingAllele,
							# freq=0
						# )
					# )
				# }
				
			# }
			
			# fromData <- fromData[fromData$allele == allele, ]
			
			# # assign allelic frequencies to populations
			# pop$DUMMY <- NA
			# fieldName <- paste0('freqLocus', locus, 'allele', allele)
			# names(pop)[ncol(pop)] <- fieldName
			
			# for (site in unique(pop$site)) {

				
				# allelicFreq <- if (site == 'NE-4/NE-5') {
					# fromData$freq[fromData$site == 'NE-5'] # only NE-5 sampled
				# } else if (site == 'NE-3/NE-7') { # neither sampled
					# NA
				# } else {
					# allelicFreq <- fromData$freq[fromData$site == site]
				# }
				
				# if (length(allelicFreq) == 0) allelicFreq <- NA
				
				# pop[pop$site == site, fieldName] <- allelicFreq
				
			# }
			
		# }

	# ### match populations with genetic sample size
	# ##############################################
		
		# fromData <- read.csv('./Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv', as.is=TRUE)
		
		# fromData$site[fromData$site == 'IL-02'] <- 'IL-2'
		# fromData$site[fromData$site == 'NE-4' | fromData$site == 'NE-5'] <- 'NE-4/NE-5'
		# fromData$site[fromData$site == 'NE-3' | fromData$site == 'NE-7'] <- 'NE-3/NE-7'

		# pop$geneticSampleSize <- NA
		
		# for (site in unique(fromData$site)) {
		
			# nPops <- length(which(fromData$site == site))
			# stopifnot(nPops > 0)
			# pop$geneticSampleSize[pop$site == site] <- nPops
			
		# }

	# write.csv(pop, './Pheno-Geno-Biogeographic Data/01 By Population - Pheno-Geno-Biogeographic Data.csv', row.names=FALSE)

# say('###########################')
# say('### create study region ###')
# say('###########################')		

	# say('Doing this in R maxes out the memory on a 24-GB machine, so I will do it in ArcGIS. The steps in R are shown below, nonetheless.')

	# gadm <- shapefile('C:/ecology/Political Geography/GADM/ver2/gadm2_northAmericaAndCentralAmerica_sansAlaska')
	
	# gadmAg <- readRDS('Species Records V1/!13b GADM Ver 2 - Dissolved to Level 2 - Multipart - WGS84 - Training Region - Tallied Poaceae in Each County.rds')
	# gadmAg <- gadmAg[gadmAg$anyAg1to3 > 0, ]
	
	# gadmAg <- gUnaryUnion(gadmAg)
	
	# # create 2000-km buffer around any county with AG
	# gadmAg <- sp::spTransform(gadmAg, getCRS('albersNA'))
	# agBuffer <- gBuffer(gadmAg, width=1500000)
	
	# # crop
	# studyRegion <- gIntersection(agBuffer, gadm, byid = TRUE, drop_lower_td = TRUE)

# say('##################################')
# say('### extract environmental data ###')
# say('##################################')		

	# say('Using county-level mean of each predictor following ....[what was that paper?].')
	
	# gadm <- readRDS('./Extents/GADM Ver 2 - Dissolved to Level 2 - Multipart - WGS84 - Clipped to Projection Region.rds')
	
	# ### environmental data
	
	# envStack <- stack(c('./WORLDCLIM Ver 2 Rel June 1 2016/elevation.tif', listFiles('./WORLDCLIM Ver 2 Rel June 1 2016/1970-2000')))
	# env <- raster::extract(envStack, gadm, df=TRUE, weights=TRUE)
	
	# for (thisPred in c('elevation', 'solarRadiation', 'wc01', 'wc02', 'wc05', 'wc07', 'wc12')) {
		
		# say(thisPred)
		
		# weighted <- env
		# weighted$weightedVal <- weighted[ , thisPred] * weighted$weight
		# weighted <- aggregate(weighted, by=list(weighted$ID), sum, na.rm=TRUE)
		# weighted$ID <- NULL
		# names(weighted)[1] <- 'ID'
		
		# gadm$DUMMY <- NA
		# names(gadm)[ncol(gadm)] <- thisPred
		# gadm@data[ , ncol(gadm)] <- weighted$weightedVal
		
	# }
	
	# ### county area
	# gadmProj <- sp::spTransform(gadm, getCRS('climateNA'))
	# countyArea <- gArea(gadmProj, byid=TRUE)
	# countyArea <- countyArea / (1000^2)
	# gadm$area_km2 <- countyArea
	
	# save(gadm, file='./Extents/GADM Ver 2 - Dissolved to Level 2 - Multipart - WGS84 - Clipped to Projection Region with WORLDCLIM Ver 2 Rel June 1 2016.RData')

# say('#################################################################################')
# say('### merge shapefile with environmental data with shapefile with specimen data ###')
# say('#################################################################################')
	
	# load('./Extents/GADM Ver 2 - Dissolved to Level 2 - Multipart - WGS84 - Clipped to Projection Region with WORLDCLIM Ver 2 Rel June 1 2016.RData')
	
	# gadmAg <- shapefile('./Species Records V2/!13c_GADM_Ver_2_-_Dissolved_to_Level_2_-_Multipart_-_WGS84_-_Training_Region_-_RE-Tallied_Poaceae_in_Each_County')
	
	# gadmNames <- paste0(gadm$NAME_0, gadm$NAME_1, gadm$NAME_2)
	# gadmAgNames <- paste0(gadmAg$NAME_0, gadmAg$NAME_1, gadmAg$NAME_2)
	
	# for (field in c('numCrd1', 'numCrd2', 'numCrd3', 'anyAg1to3', 'poaRec')) {
	
		# gadm$DUMMY <- as.integer(rep(NA, nrow(gadm)))
		# names(gadm)[ncol(gadm)] <- field
		
		# for (i in 1:nrow(gadmAg)) {
			# index <- which(gadmNames == gadmAgNames[i])
			# gadm@data[index, field] <- as.integer(gadmAg@data[i , field])
		# }
	
	# }
	
	# gadm$agDensity <- gadm$anyAg1to3 / gadm$areaKm2
	# gadm$poaDensity <- gadm$poaRec / gadm$areaKm2
	
	# dirCreate('./Species Records V3')
	
	# write.csv('This data file is the merging of 13c in the "Species Records V2" folder AND the RData file "GADM Ver 2 - Dissolved to Level 2 - Multipart - WGS84 - Clipped to Projection Region with WORLDCLIM Ver 2 Rel June 1 2016.RData" in the "Extents" folder. Note that counties with no Poa records were not tallied (they may indeed have Poacea records.)', './Species Records V3/!README.txt', row.names=FALSE)
	
	# save(gadm, file='./Species Records V3/!13c GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Records.RData')
	
# say('##############################################################################')
# say('### merge occurrence data and population-level phenotype and genotype data ###')
# say('##############################################################################')

	# load('./Species Records V3/!13c GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Records.RData')
	# pop <- read.csv('./Pheno-Geno-Biogeographic Data/01 By Population - Pheno-Geno-Biogeographic Data.csv')

	# popSp <- SpatialPointsDataFrame(cbind(pop$longitude, pop$latitude), data=pop, proj4string=CRS(projection(gadm)))
	# countiesWithPop <- sp::over(popSp, gadm)

	# gadmNames <- paste0(gadm$NAME_0, gadm$NAME_1, gadm$NAME_2)
	# popNames <- paste0(countiesWithPop$NAME_0, countiesWithPop$NAME_1, countiesWithPop$NAME_2)

	# gadmIndex <- match(popNames, gadmNames)
	
	# # genetic/phenotypic data
	# for (field in c('site', 'longitude', 'heightCm', 'heightCmLower', 'heightCmUpper', 'nHeight', 'heightJulianTIP', 'spad', 'spadLower', 'spadUpper', 'nSpad', 'spadJulianTIP', 'leafWidthCm', 'leafWidthCmLower', 'leafWidthCmUpper', 'nLeafWidth', 'leafWidthJulianTIP', 'biomassG', 'biomassGLower', 'biomassGUpper', 'nBiomass', 'biomassJulianTIP', 'booted', 'bootDateJulian', 'bootJulianTIP', 'transplantJulian', 'structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4', 'freqLocusS1_2050023allele110', 'freqLocusS1_27189028allele120', 'freqLocusS1_33964042allele100', 'freqLocusS1_73120031allele110', 'geneticSampleSize')) {
	
		# gadm@data$DUMMY <- NA
		# names(gadm@data)[ncol(gadm@data)] <- field
		# gadm@data[gadmIndex, field] <- pop[ , field]
		
	# }
	
	# # remove counties that are Great Lakes
	# gadm <- gadm[-which(gadm$NAME_2 == 'Lake Superior'), ]
	# gadm <- gadm[-which(gadm$NAME_2 == 'Lake Michigan'), ]
	# gadm <- gadm[-which(gadm$NAME_2 == 'Lake Hurron'), ]
	# gadm <- gadm[-which(gadm$NAME_2 == 'Lake Erie'), ]
	# gadm <- gadm[-which(gadm$NAME_2 == 'Lake Ontario'), ]

	# gadm <- gadm[-which(gadm$NAME_2 == 'Kenora'), ]
		
	# save(gadm, file='./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')
	
# say('###############################')
# say('### GJAM @ population level ###')
# say('###############################')

	# dirCreate('./Models - PBGMs/GJAM @ Population Level')

	# load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')

	# ### responses
		
		# basePhenoResponses <- c('biomassG', 'leafWidthCm', 'heightCm', 'spad')
		# phenoResponses <- c('biomassGLower', 'biomassG', 'biomassGUpper', 'leafWidthCmLower', 'leafWidthCm', 'leafWidthCmUpper', 'heightCmLower', 'heightCm', 'heightCmUpper', 'spadLower', 'spad', 'spadUpper')
		# structureGroupResponses <- c('structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')
		# # responses <- c(phenoResponses, structureGroupResponses)
		# responses <- c(phenoResponses, structureGroupResponses)
		
		# typeNames <- c(rep('CA', length(phenoResponses)), rep('FC', length(structureGroupResponses)))
		
		# y <- gadm@data[!is.na(gadm@data$site), responses]

	# ### effort
			
		# effort <- NA * y
		# for (resp in responses) {
			
			# effortCol <- if (resp %in% c('biomassG', 'biomassGLower', 'biomassGUpper')) {
				# 'nBiomass'
			# } else if (resp %in% c('leafWidthCm', 'leafWidthCmLower', 'leafWidthCmUpper')) {
				# 'nLeafWidth'
			# } else if (resp %in% c('heightCm', 'heightCmLower', 'heightCmUpper')) {
				# 'nHeight'
			# } else if (resp %in% c('spad', 'spadLower', 'spadUpper')) {
				# 'nSpad'
			# } else if (resp %in% c('freqLocusS1_2050023allele110', 'freqLocusS1_27189028allele120', 'freqLocusS1_33964042allele100', 'freqLocusS1_73120031allele110', 'structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')) {
				# 'geneticSampleSize'
			# }
			
			# effort[!is.na(y[ , resp]), resp] <- gadm@data[!is.na(gadm@data[ , effortCol]), effortCol]
			# if (anyNA(effort[ , resp])) effort[is.na(effort[ , resp]), resp] <- 0
			
		# }
		
		# effort <- list(columns=names(effort), values=as.matrix(effort))

	# ### predictors
			
		# x <- gadm@data[!is.na(gadm@data$site), predictors]
		# pca <- princomp(x, cor=TRUE)
		# # biplot(pca)
		# x <- as.data.frame(pca$scores)
		# names(x) <- paste0('pc', prefix(seq_along(predictors), 2))
		# save(pca, file='./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.RData')

	# ### model: PHENOTYPE and STRUCTURE
			
		# form <- as.formula(~ pc01 + pc02 + pc01:pc02 + I(pc01^2) + I(pc02^2))
		
		# tuning <- list(ng = 5000, burnin = 1000, typeNames = typeNames, effort = effort)

		# model <- gjam(form, xdata = x, ydata = y, modelList = tuning)
		# summary(model)
		
		# save(model, file='./Models - PBGMs/GJAM @ Population Level/GJAM on Phenotype & STRUCTURE Groups.RData')
		
# say('####################################')
# say('### plot GJAM @ population level ###')
# say('####################################')

	# ### generalization
	# ##################

		# # states
		# states <- readRDS('./Extents/GADM Ver 2 - Dissolved to Level 1 - Multipart - WGS84 - Clipped to Projection Region.rds')
		
		# # shapefile with predictors/response data
		# load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')
		
		# # GJAM
		# load('./Models - PBGMs/GJAM @ Population Level/GJAM on Phenotype & STRUCTURE Groups.RData')
		
		# # LOPOD
		# load('./Models - PBGMs/bayesLopod with Variable Detection/LOPOD Prediction.RData')
		# load('./Models - PBGMs/bayesLopod with Constant Detection/LOPOD Prediction.RData')

		# # responses
		# basePhenoResponses <- c('biomassG', 'leafWidthCm', 'heightCm', 'spad')
		# phenoResponses <- c('biomassGLower', 'biomassG', 'biomassGUpper', 'leafWidthCmLower', 'leafWidthCm', 'leafWidthCmUpper', 'heightCmLower', 'heightCm', 'heightCmUpper', 'spadLower', 'spad', 'spadUpper')
		# structureGroupResponses <- c('structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')
		# responses <- c(phenoResponses, structureGroupResponses)

	# # say('### observed vs predicted')
	# # say('#########################')

		# # png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted (GJAM).png', height=1400, width=2000)
		
			# # layout(matrix(c(1:4, rep(5, 4)), ncol=4, byrow=TRUE))
			# # par(cex=1.1)
		
			# # ### pheno: observed vs expected
			
			# # par(pty='s', cex=1.4)
			# # for (resp in basePhenoResponses) {
				
				# # niceName <- respNice(resp)

				# # y <- model$inputs$y
				
				# # # plot limits
				# # resps <- c(y[ , resp], y[ , paste0(resp, 'Lower')], y[ , paste0(resp, 'Upper')], model$prediction$ypredMu[ , resp], model$prediction$ypredMu[ , paste0(resp, 'Lower')], model$prediction$ypredMu[ , paste0(resp, 'Upper')])
				# # nice <- pretty(resps, 3)
				# # lims <- range(nice)
				
				# # # colors and R2
				# # cols <- c('#1b7837', 'white', '#762a83')
				# # r2 <- rep(NA, 3)
				# # names(cols) <- names(r2) <- c('Upper', 'Mean', 'Lower')
				
				# # # plot
				# # for (type in c('Mean', 'Lower', 'Upper')) {
					
					# # if (type == 'Mean') {
						# # typeName <- ''
						# # col <- cols[typeName]
					# # } else {
						# # typeName <- type
						# # col <- alpha(cols[typeName], 0.6)
					# # }
					
					# # obs <- y[ , paste0(resp, typeName)]
					# # pred <- model$prediction$ypredMu[ , paste0(resp, typeName)]
					
					# # r2[type] <- cor(obs, pred)
				
					# # if (type == 'Mean') {
						# # plot(obs, pred, xlim=lims, ylim=lims, xlab=paste('Observed', niceName$unit), ylab=paste('Predicted', niceName$unit), main=niceName$nice, pch=21, cex=1.9, bg=col)
						# # abline(0, 1)
						# # points(obs, pred, pch=21, cex=1.9)
					# # } else {
						# # points(obs, pred, bg=col, pch=21, cex=1.9)
					# # }
					
					
				# # }
				
				# # legend('bottomright', inset=0.02, legend=c('97.5% Quantile', 'Population Mean', '2.5% Quantile'), pt.bg=cols, pch=21, pt.cex=1.2)
				
				# # # R2
				# # box <- par('usr')
				# # xPos <- box[1] + 0.02 * (box[2] - box[1])
				# # yPos <- box[4] - 0.02 * (box[4] - box[3]) - 0.05 * (1:4) * (box[4] - box[3])
				# # text(xPos, yPos, labels=paste(c('Pseudo-R2', 'Upper:', 'Mean:', 'Lower:'), c('', sprintf('%.2f', r2))), pos=4)
				
			# # }
		
			# # ### structure groups: left bar is observed, right is predicted
		
			# # par(pty='m')
			# # structObs <- t(y[ , structureGroupResponses])
			# # structPred <- t(model$prediction$ypredMu[ , structureGroupResponses])

			# # sites <- gadm@data$site[!is.na(gadm@data$site)]
			# # obsPred <- matrix(NA, ncol=1, nrow=length(structureGroupResponses))
			# # for (i in seq_along(sites)) obsPred <- cbind(obsPred, cbind(structObs[ , i], structPred[ , i]))
			# # obsPred <- obsPred[ , -1]
			
			# # for (i in 1:ncol(obsPred)) if (!anyNA(obsPred[ , i])) obsPred[ , i] <- obsPred[ , i] / sum(obsPred[ , i])

			# # sites <- gadm@data$site[!is.na(gadm@data$site)]
			# # cols <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a')
			# # barplot(obsPred, ylim=c(0, 1), space=c(0.7, 0), yaxt='n', col=cols, main='STRUCTURE Groups')
			# # axis(2, pos=-0.5)
			# # text(-5, 0.5, xpd=NA, labels='Probability', srt=90, cex=1.1)
			# # text(-0.5 + 2.7 * seq_along(sites), -0.02, labels=sites, xpd=NA, srt=90, adj=c(1, 0))
			# # text(mean(-0.5 + 2.7 * seq_along(sites)), 1.03, labels='For each site the left bar is observed and the right is predicted. Missing observed bars were not genotyped.', xpd=NA, cex=1)
			# # legend('bottom', inset=c(0.5, -0.19), xpd=NA, ncol=4, fill=cols, legend=paste('Group', 1:4), cex=0.85)

			# # title(main=paste0('GJAM using 1st Two Axes of PCA on BIOCLIM Variables\n', date()), outer=TRUE, line=-3.5, cex.main=1.6)
			
		# # dev.off()

	# say('### map of predicted phenotypes')
	# say('###############################')

		# # # generalization
		# # # lopod <- lopodPredictVarP
		# # lopod <- lopodPredictConstantP
	
		# # # define map extent
		# # gadm <- sp::spTransform(gadm, getCRS('albersNA', TRUE))
		# # mapExtent <- gadm[!is.na(gadm$heightCm), ]
		# # mapExtent <- gBuffer(mapExtent, width=200000)
		# # mapExtent <- as(extent(mapExtent), 'SpatialPolygons')
		# # projection(mapExtent) <- getCRS('albersNA')
		
		# # lopod <- sp::spTransform(lopod, getCRS('albersNA'))
		
		# # load('./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.RData')
		# # pcaPred <- as.data.frame(predict(pca, gadm@data))
		# # names(pcaPred) <- paste0('pc', prefix(seq_along(predictors), 2))
		# # gadm@data <- cbind(gadm@data, pcaPred)
	
		# # # predict GJAM to GADM (present)
		# # pred <- gjamPredict(model, newdata=list(xdata=gadm@data[ , paste0('pc', prefix(seq_along(predictors), 2))]))
		# # gadmPred <- gadm
		# # for (resp in basePhenoResponses) gadmPred@data[ , resp] <- pred$sdList$yMu[ , resp]
		# # for (resp in basePhenoResponses) gadmPred@data[ , paste0(resp, 'Lower')] <- pred$sdList$yMu[ , paste0(resp, 'Lower')]
		# # for (resp in basePhenoResponses) gadmPred@data[ , paste0(resp, 'Upper')] <- pred$sdList$yMu[ , paste0(resp, 'Upper')]
	
		# # states <- sp::spTransform(states, getCRS('albersNA', TRUE))
	
		# # countiesInMap <- crop(gadmPred, mapExtent)
		# # statesInMap <- crop(states, mapExtent)
		# # lopodInMap <- crop(lopod, mapExtent)
		
		# # colors
		# cols <- c('#ff7f00', '#6a3d9a', 'darkgreen')
		# pal <- colorRampPalette(cols)
		# pal <- pal(101)
		
		# png(paste0('./Models - PBGMs/GJAM @ Population Level/Phenotype Maps for 1970-2000.png'), height=2 * 900, width=2600, res=450)

			# # two columns, one row: current mean prediction and current predicted range
			# par(mfrow=c(2, 2), cex=1.1, mar=0.2 * c(1, 2, 2, 6))

			# count <- 1
			
			# # by RESPONSE
			# for (resp in basePhenoResponses) {
			# # for (resp in basePhenoResponses[1]) {
				
				# say(resp)
				# niceName <- respNice(resp)

				# ### plot mean response
				# ######################
				
				# # range of observed values
				# obs <- na.omit(gadm@data[ , resp])
				# minOfObs <- min(obs, na.rm=TRUE)
				# maxOfObs <- max(obs, na.rm=TRUE)
				# maxOfObs <- 1.1 * maxOfObs
				# minOfObs <- max(0, 0.9 * minOfObs)
				# prettyObs <- pretty(c(minOfObs, maxOfObs), 3)
				
				# # predicted values
				# vals <- countiesInMap@data[ , resp]
				# valsScaled <- vals - minOfObs
				# valsScaled <- valsScaled / (maxOfObs - minOfObs)
				# if (any(valsScaled < 0)) valsScaled[valsScaled < 0] <- 0
				# if (any(valsScaled > 1)) valsScaled[valsScaled > 1] <- 1
				# col <- pal[round(100 * valsScaled) + 1]
				
				# # fade by probability of presence
				# prPres <- lopodInMap$psi_i
				
				# plot(countiesInMap, col=alpha(col, prPres), border=NA)
				# plot(statesInMap, col=NA, border='black', add=TRUE, lwd=0.3)

				# pos <- par('usr')
				# text(x=pos[1] + 0.01 * (pos[2] - pos[1]), y=pos[4] + 0.0 * (pos[4] - pos[3]), labels=paste0(letters[count], ') ', niceName$nice), xpd=NA, cex=0.5, xpd=NA, pos=4)

				# legendGrad('right', inset=-0.07, width=0.1, height=1.1, gradAdjX=c(0.1, 0.5), gradAdjY=c(0.08, 0.83), labels=prettyObs, col=cols, title=paste0(niceName$nice, '\n', niceName$unit), boxBorder=NULL, xpd=NA, cex=0.4, labAdj=0.4, lwd=0.5)
				
				# # observed
				# gadmObs <- gadm[!is.na(gadm@data[ , resp]), ]
				# vals <- gadmObs@data[ , resp]
				# valsScaled <- vals - minOfObs
				# valsScaled <- valsScaled / (maxOfObs - minOfObs)
				# if (any(valsScaled < 0)) valsScaled[valsScaled < 0] <- 0
				# if (any(valsScaled > 1)) valsScaled[valsScaled > 1] <- 1
				# col <- pal[round(100 * valsScaled) + 1]
				
				# prPres <- rep(NA, nrow(gadmObs))
				# for (i in 1:nrow(gadmObs)) {
					# prPres[i] <- lopod$psi_i[which(lopod$NAME_1 == gadmObs$NAME_1[i] & lopod$NAME_2 == gadmObs$NAME_2[i])]
				# }
				
				# plot(gadmObs, col='white', border=NA, add=TRUE)
				# plot(gadmObs, col=alpha(col, prPres), border='black', lwd=1, add=TRUE)
			
				# count <- count + 1
				
			# }
			
			# title(sub=date(), outer=TRUE, line=-0.85, cex.sub=0.2)
			
		# dev.off()
				
	# say('### map of STRUCTURE groups')
	# say('###########################')

		# # generalization
		# lopod <- lopodPredictVarP
		# # lopod <- lopodPredictConstantP
	
		# # define map extent
		# gadm <- sp::spTransform(gadm, getCRS('albersNA', TRUE))
		# mapExtent <- gadm[!is.na(gadm$heightCm), ]
		# mapExtent <- gBuffer(mapExtent, width=200000)
		# mapExtent <- as(extent(mapExtent), 'SpatialPolygons')
		# projection(mapExtent) <- getCRS('albersNA')
		
		# lopod <- sp::spTransform(lopod, getCRS('albersNA'))
		
		# load('./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.RData')
		# pcaPred <- as.data.frame(predict(pca, gadm@data))
		# names(pcaPred) <- paste0('pc', prefix(seq_along(predictors), 2))
		# gadm@data <- cbind(gadm@data, pcaPred)
	
		# # predict GJAM to GADM (present)
		# pred <- gjamPredict(model, newdata=list(xdata=gadm@data[ , paste0('pc', prefix(seq_along(predictors), 2))]))
		# gadmPred <- gadm
		# for (resp in 1:4) gadmPred@data[ , paste0('structureGroup', 1:4)] <- pred$sdList$yMu[ , paste0('structureGroup', 1:4)]
	
		# states <- sp::spTransform(states, getCRS('albersNA', TRUE))
	
		# countiesInMap <- crop(gadmPred, mapExtent)
		# statesInMap <- crop(states, mapExtent)
		# lopodInMap <- crop(lopod, mapExtent)
		
		# cols <- c('yellow', 'darkblue', 'black', 'red')
		
		# png(paste0('./Models - PBGMs/GJAM @ Population Level/STRUCTURE Maps for 1970-2000 (GJAM).png'), height=900, width=1400, res=450)

			# # two columns, one row: current mean prediction and current predicted range
			# par(cex=1.1, mar=0.2 * c(1, 2, 2, 10))

			# ### predicted
			
			# # colors
			# vals <- as.data.frame(countiesInMap[ , paste0('structureGroup', 1:4)])
			# valsToCols <- apply(vals, 1, colorFrom4Vector, cols=cols)
			# predCols <- matrix(NA, ncol=3, nrow=nrow(vals))
			# xy <- matrix(NA, ncol=2, nrow=nrow(vals))
			
			# for (i in 1:nrow(vals)) {
				# predCols[i, ] <- valsToCols[[i]]$col
				# xy[i, ] <- valsToCols[[i]]$xy
			# }
			
			# colnames(predCols) <- c('h', 's', 'v')
			# colnames(xy) <- c('x', 'y')
			
			# # fade by probability of presence
			# prPres <- lopodInMap$psi_i
			# if (anyNA(prPres)) prPres[is.na(prPres)] <- 0
			
			# plot(lopodInMap, col=hsv(h=predCols[ , 'h'], s=predCols[ , 's'], v=predCols[ , 'v'], alpha=prPres), border=NA)
			# plot(statesInMap, col=NA, border='black', add=TRUE, lwd=0.5)

			# pos <- par('usr')
			# text(x=pos[1] + 0.01 * (pos[2] - pos[1]), y=pos[4] + 0.0 * (pos[4] - pos[3]), labels=paste0('a) STRUCTURE groups'), xpd=NA, cex=0.3, xpd=NA, pos=4)

			# legendQuad('topright', inset=c(-0.13, 0.038), width=0.17, height=0.25, swatchAdjX=c(0.13, 0.87), swatchAdjY=c(0.18, 0.82), labels=LETTERS[1:4], title='', boxBorder=NULL, xpd=NA, cex=0.3, labelAdj=1.1, cols=cols, aspect=FALSE, lwd=0.5)
			
			# ### observed
			# obs <- gadm[ , paste0('structureGroup', 1:4)]
			# nas <- naRows(obs)
			# obs <- obs[-nas, ]
			# where <- coordinates(gCentroid(obs, byid=TRUE))
			
			# for (i in 1:nrow(where)) pies(unlist(as.data.frame(obs[i, ])), add=TRUE, xPos=where[i, 1], yPos=where[i, 2], radius=40000, col=cols, aspect=FALSE, lwd=0.5)
			
			# title(sub=date(), outer=TRUE, line=-0.85, cex.sub=0.1)
			
		# dev.off()
				





				
	
	# say('### model diagnostics: phenotype: observed vs expected in environmental space')
	
	# png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted in Environmental Space (GJAM).png', height=1800, width=1600, res=300)

		# par(mfrow=c(4, 3), mar=0.3 * c(7.5, 6, 3, 5.5) + 0.1, cex.main=0.8, cex.lab=0.7, cex.axis=0.5, mgp=c(3, 0.2, 0))
	
		# inc <- 20 # number of increments into which to divide each environmental axis
		
		# pcs <- c('pc01', 'pc02')
	
		# pc01 <- seq(1.05 * min(x$pc01), 1.05 * max(x$pc01), length.out=inc)
		# pc02 <- seq(1.05 * min(x$pc02), 1.05 * max(x$pc02), length.out=inc)
		
		# xSpace <- expand.grid(pc01, pc02)
		# names(xSpace) <- as.matrix(pcs)
		# newdata <- list(xdata=xSpace)
		
		# preds <- gjamPredict(model, newdata=newdata); say('')

		# for (resp in phenoResponses) {
		
			# niceName <- respNice(resp)
			# niceName <- paste0(niceName$nice, niceName$unit)
		
			# plot(0, type='n', xlim=range(xSpace$pc01), ylim=range(xSpace$pc02), pty='s', main=niceName, xaxt='n', yaxt='n', bty='n', pty='s', xlab='', ylab='')
			# axis(1, tck=-0.02)
			# axis(2, tck=-0.02)
			# title(xlab=toupper(pcs[1]), line=1)
			# title(ylab=toupper(pcs[2]), line=1)

			# # background
			# predObs <- y[ , resp]
			# predBg <- preds$sdList$yMu[ , resp]
			# cols <- predBg - min(predBg, predObs)
			# cols <- cols / max(cols)
			
			# for (i in 1:(nrow(xSpace) - inc)) {
			
				# xs <- c(xSpace$pc01[i], xSpace$pc01[i+1], xSpace$pc01[i+1], xSpace$pc01[i])
				# ys <- c(xSpace$pc02[i], xSpace$pc02[i], xSpace$pc02[i+inc], xSpace$pc02[i+inc])
			
				# polygon(
					# x=xs,
					# y=ys,
					# col=alpha('darkgreen', cols[i]),
					# border=NA
				# )
			
			# }

			# # legend
			# legendGrad('bottomright', inset=c(-0.08, 0), width=0.1, height=1, labels=sprintf('%.1f', round(c(min(predBg), mean(predBg), max(predBg)), 1)), col=c('white', 'darkgreen'), xpd=NA, boxBorder=NA, gradAdjY=c(0.0345, 0.9645), labAdj=0.35, cex=0.6, boxBg=NA)
			
			# # loadings
			# mag <- 3.5
			# loads <- as.matrix(pca$loadings)
			# for (pred in predictors) {
				# xs <- mag * loads[pred, 'Comp.1']
				# ys <- mag * loads[pred, 'Comp.2']
				# arrows(x0=0, y0=0, x1=xs, y1=ys, angle=25, length=0.05, lwd=0.5)
				# text(xs, ys, labels=predNice(pred)$nice, pos=4, cex=0.6, lwd=0.5, xpd=NA)
			# }
			
			# # observed
			# cols <- predObs - min(predBg, predObs)
			# cols <- cols / max(cols)
			
			# points(x[ , pcs], pch=16, col='white', cex=1, xpd=NA)
			# points(x[ , pcs], pch=21, bg=alpha('darkgreen', cols), cex=1, xpd=NA)
			
		# }
		
		# title(sub=paste0('GJAM Using 1st Two PC Axes vs Population-Level Phenotypes | ', date()), outer=TRUE, line=-0.9, cex.sub=0.6)
		
	# dev.off()
	
# say('### model diagnostics: STRUCTURE in environmental space')
# # source('C:/ecology/Drive/Research/Andropogon/Scripts - PGBMs/TEMP.r')

	# png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted STRUCTURE Population in Environmental Space (GJAM).png', height=1800, width=1600, res=300)

		# par(mfrow=c(1, 1), mar=0.3 * c(7.5, 6, 3, 5.5) + 0.1, cex.main=0.8, cex.lab=0.7, cex.axis=0.5, mgp=c(3, 0.2, 0))
	
		# # inc <- 100 # number of increments into which to divide each environmental axis
		# inc <- 20 # number of increments into which to divide each environmental axis
		
		# pcs <- c('pc01', 'pc02')
	
		# pc01 <- seq(1.05 * min(x$pc01), 1.05 * max(x$pc01), length.out=inc)
		# pc02 <- seq(1.05 * min(x$pc02), 1.05 * max(x$pc02), length.out=inc)
		
		# xSpace <- expand.grid(pc01, pc02)
		# names(xSpace) <- as.matrix(pcs)
		# newdata <- list(xdata=xSpace)
		
		# # predSpace <- gjamPredict(model, newdata=newdata); say('')

		# # plot
		# plot(0, type='n', xlim=range(xSpace$pc01), ylim=range(xSpace$pc02), pty='s', main='STRUCTURE Populations', xaxt='n', yaxt='n', bty='n', pty='s', xlab='', ylab='')
		# axis(1, tck=-0.02)
		# axis(2, tck=-0.02)
		# title(xlab=toupper(pcs[1]), line=1)
		# title(ylab=toupper(pcs[2]), line=1)

		# # background 
		# predBg <- predObs <- list()
		# for (i in 1:4) {
			# predObs[[i]] <- y[ , paste0('structureGroup', i)]
			# predBg[[i]] <- predSpace$sdList$yMu[ , paste0('structureGroup', i)]
		# }
		
		# for (i in 1:(nrow(xSpace) - inc)) {
		
			# xs <- c(xSpace$pc01[i], xSpace$pc01[i+1], xSpace$pc01[i+1], xSpace$pc01[i])
			# ys <- c(xSpace$pc02[i], xSpace$pc02[i], xSpace$pc02[i+inc], xSpace$pc02[i+inc])

			# # for (group in 1:4) {
			# for (group in 1:2) {
			
				# # rgb(1, 1, 1) is all black
				
				# col <- rgb(predBg[[1]][i], predBg[[2]][i], predBg[[3]][i], 1 - predBg[[4]][i])
				# col <- if (group == 1) {
					# rgb(predBg[[1]][i], 0, 0, 1 - predBg[[1]][i])
				# } else if (group == 2) {
					# rgb(0, predBg[[2]][i], 0, 1 - predBg[[2]][i])
				# } else if (group == 3) {
					# rgb(0, 0, predBg[[3]][i], 1 - predBg[[3]][i])
				# } else if (group == 4) {
					# rgb(1, 1, 1, 1 - predBg[[4]][i])
				# }
				
				# polygon(
					# x=xs,
					# y=ys,
					# col=col,
					# border=NA
				# )
				
			# }
		
		# }

		# # legend
		# col <- c('white', 'blue', 'green', 'red')
		# legendGrad('bottomright', inset=c(-0.08, 0), width=0.1, height=1, labels=4:1, col=col, xpd=NA, boxBorder=NA, gradAdjY=c(0.0345, 0.9645), labAdj=0.35, cex=0.6, boxBg=NA)
		
		# # loadings
		# mag <- 3.5
		# loads <- as.matrix(pca$loadings)
		# for (pred in predictors) {
			# xs <- mag * loads[pred, 'Comp.1']
			# ys <- mag * loads[pred, 'Comp.2']
			# arrows(x0=0, y0=0, x1=xs, y1=ys, angle=25, length=0.05, lwd=0.5)
			# text(xs, ys, labels=predNice(pred)$nice, pos=4, cex=0.6, lwd=0.5, xpd=NA)
		# }
		
		# # observed
		# isNa <- is.na(predObs[[1]])
		# predObsNoNa <- lapply(predObs, na.omit)
		# col <- rgb(predObsNoNa[[1]], predObsNoNa[[2]], predObsNoNa[[3]], 1 - predObsNoNa[[4]])
		
		# points(x[isNa, pcs], pch=0, cex=3, xpd=NA)
		
		# sitesNoNa <- sites[!isNa]
		
		# for (i in 1:sum(!isNa)) {

			# cols <- c('red', 'green', 'blue', 'white')
		
			# pies(x=c(predObsNoNa[[1]][i], predObsNoNa[[2]][i], predObsNoNa[[1]][3], predObsNoNa[[4]][i]), add=TRUE, xPos=x[!isNa, pcs[1]][i], yPos=x[!isNa, pcs[2]][i], radius=0.15, col=cols)
			
			# text(x[!isNa, pcs[1]][i], x[!isNa, pcs[2]][i], labels=sitesNoNa[i], adj=c(-0.5, -0.5), xpd=NA, cex=0.8)
		
		# }

		# title(sub=paste0('GJAM Using 1st Two PC Axes vs Population-Level STRUCTURE Groups | ', date()), outer=TRUE, line=-0.9, cex.sub=0.6)
		
	# dev.off()
	
	
# say('######################################################################')
# say('### GAM @ population level vs STRUCTURE groups & candidate alleles ###')
# say('######################################################################')

	# dir <- 'GAM @ Population Level vs STRUCTURE & Candidate Alleles'
	# dirCreate('./Models - PBGMs/', dir)

	# load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')

	# ### responses and effort
		
		# responsesStruct <- c('structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')
		# responsesAlleles <- c('freqLocusS1_2050023allele110', 'freqLocusS1_27189028allele120', 'freqLocusS1_33964042allele100', 'freqLocusS1_73120031allele110')
		
		# data <- gadm@data[!is.na(gadm@data$site), c(responsesStruct, responsesAlleles, predictors, 'geneticSampleSize')]
		# data <- data[-naRows(data), ]
		# data$effort <- data[ , 'geneticSampleSize'] / max(data[ , 'geneticSampleSize'])
		
	# ### predictors
			
		# pca <- princomp(data[ , predictors], cor=TRUE)
		# # biplot(pca)
		# x <- as.data.frame(pca$scores)
		# names(x) <- paste0('pc', prefix(seq_along(predictors), 2))
		# data <- cbind(data, x)
		# save(pca, file=paste0('./Models - PBGMs/', dir, '/PCA on Environmental Predictors at Sampled Populations.RData'))
		
		# data$pc01pc02 <- data$pc01 * data$pc02

	# ### model: STRUCTURE
		
		# for (i in seq_along(responsesStruct)) {
		
			# resp <- paste0('structureGroup', i)
			# model <- trainGam(data=data, resp=resp, preds=c('pc01', 'pc02', 'pc01pc02'), family=betar, interaction=NULL, w=data$effort, verbose=FALSE)
			# save(model, file=paste0('./Models - PBGMs/', dir, '/GAM vs STRUCTURE Population ', i, '.RData'))
			
		# }
		
	# ### model: alleles
		
		# for (i in seq_along(responsesAlleles)) {
		
			# resp <- responsesAlleles[i]
			# model <- trainGam(data=data, resp=resp, preds=c('pc01', 'pc02', 'pc01pc02'), family=betar, interaction=NULL, w=data$effort, verbose=FALSE)
			# save(model, file=paste0('./Models - PBGMs/', dir, '/GAM vs Allele ', resp, '.RData'))
			
		# }
		
		
# say('##########################################################################')
# say('### map GAM @ population level vs STRUCTURE groups & candidate alleles ###')
# say('##########################################################################')

	# # ### generalization
	# # ##################

		# dir <- 'GAM @ Population Level vs STRUCTURE & Candidate Alleles'

		# # states
		# states <- readRDS('./Extents/GADM Ver 2 - Dissolved to Level 1 - Multipart - WGS84 - Clipped to Projection Region.rds')
		
		# # shapefile with predictors/response data
		# load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')
		
		# # GJAM
		# load('./Models - PBGMs/GJAM @ Population Level/GJAM on Phenotype & STRUCTURE Groups.RData')
		
		# # LOPOD
		# load('./Models - PBGMs/bayesLopod with Variable Detection/LOPOD Prediction.RData')
		# load('./Models - PBGMs/bayesLopod with Constant Detection/LOPOD Prediction.RData')

		# # responses
		# responsesStruct <- c('structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')
		# responsesAlleles <- c('freqLocusS1_2050023allele110', 'freqLocusS1_27189028allele120', 'freqLocusS1_33964042allele100', 'freqLocusS1_73120031allele110')

		# # generalization
		# lopod <- lopodPredictVarP
		# # lopod <- lopodPredictConstantP
	
		# # define map extent
		# gadm <- sp::spTransform(gadm, getCRS('albersNA', TRUE))
		# mapExtent <- gadm[!is.na(gadm$heightCm), ]
		# mapExtent <- gBuffer(mapExtent, width=200000)
		# mapExtent <- as(extent(mapExtent), 'SpatialPolygons')
		# projection(mapExtent) <- getCRS('albersNA')
		
		# lopod <- sp::spTransform(lopod, getCRS('albersNA'))
		
		# load('./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.RData')
		# pcaPred <- as.data.frame(predict(pca, gadm@data))
		# names(pcaPred) <- paste0('pc', prefix(seq_along(predictors), 2))
		# gadm@data <- cbind(gadm@data, pcaPred)
		# gadm@data$pc01pc02 <- gadm@data$pc01 * gadm@data$pc02
	
	# say('### map of STRUCTURE groups')
	# say('###########################')

		# responses <- responsesStruct
	
		# # predict GAM to GADM (present)
		# pred <- matrix(NA, ncol=length(responses), nrow=nrow(gadm))
		# colnames <- responses
		
		# for (i in seq_along(responses)) {
		
			# load(paste0('./Models - PBGMs/', dir, '/GAM vs STRUCTURE Population ', i, '.RData'))
			# pred[ , i] <- predict(model, newdata=gadm@data, type='response')
			
		# }
		
		# # standardize
		# sums <- rowSums(pred)
		# for (i in seq_along(responses)) pred[ , i] <- pred[ , i] / sums
		
		# gadmPred <- gadm
		# gadmPred@data[ , responses] <- pred
	
		# states <- sp::spTransform(states, getCRS('albersNA', TRUE))
	
		# countiesInMap <- crop(gadmPred, mapExtent)
		# statesInMap <- crop(states, mapExtent)
		# lopodInMap <- crop(lopod, mapExtent)
		
		# cols <- c('yellow', 'darkblue', 'black', 'red')
		
		# png(paste0('./Models - PBGMs/', dir, '/STRUCTURE Maps for 1970-2000 (Individual GAMs).png'), height=900, width=1400, res=450)

			# par(cex=1.1, mar=0.2 * c(1, 2, 2, 10))

			# ### predicted
			
			# # colors
			# vals <- as.data.frame(countiesInMap[ , responses])
			# valsToCols <- apply(vals, 1, colorFrom4Vector, cols=cols)
			# predCols <- matrix(NA, ncol=3, nrow=nrow(vals))
			# xy <- matrix(NA, ncol=2, nrow=nrow(vals))
			
			# for (i in 1:nrow(vals)) {
				# predCols[i, ] <- valsToCols[[i]]$col
				# xy[i, ] <- valsToCols[[i]]$xy
			# }
			
			# colnames(predCols) <- c('h', 's', 'v')
			# colnames(xy) <- c('x', 'y')
			
			# # fade by probability of presence
			# prPres <- lopodInMap$psi_i
			# if (anyNA(prPres)) prPres[is.na(prPres)] <- 0
			
			# plot(lopodInMap, col=hsv(h=predCols[ , 'h'], s=predCols[ , 's'], v=predCols[ , 'v'], alpha=prPres), border=NA)
			# plot(statesInMap, col=NA, border='black', add=TRUE, lwd=0.5)

			# pos <- par('usr')
			# text(x=pos[1] + 0.01 * (pos[2] - pos[1]), y=pos[4] + 0.0 * (pos[4] - pos[3]), labels=paste0('a) STRUCTURE groups'), xpd=NA, cex=0.3, xpd=NA, pos=4)

			# legendQuad('topright', inset=c(-0.13, 0.038), width=0.17, height=0.25, swatchAdjX=c(0.13, 0.87), swatchAdjY=c(0.18, 0.82), labels=LETTERS[1:4], title='', boxBorder=NULL, xpd=NA, cex=0.3, labelAdj=1.1, cols=cols, aspect=FALSE, lwd=0.5)
			
			# ### observed
			# obs <- gadm[ , responses]
			# nas <- naRows(obs)
			# obs <- obs[-nas, ]
			# where <- coordinates(gCentroid(obs, byid=TRUE))
			
			# for (i in 1:nrow(where)) pies(unlist(as.data.frame(obs[i, ])), add=TRUE, xPos=where[i, 1], yPos=where[i, 2], radius=40000, col=cols, aspect=FALSE, lwd=0.5)
			
			# title(sub=date(), outer=TRUE, line=-0.85, cex.sub=0.1)
			
		# dev.off()

	# say('### maps of candidate allele frequencies')
	# say('########################################')

		# responses <- responsesAlleles
	
		# # predict GAM to GADM (present)
		# pred <- matrix(NA, ncol=length(responses), nrow=nrow(gadm))
		# colnames <- responses
		
		# for (i in seq_along(responses)) {
		
			# resp <- responses[i]
			# load(paste0('./Models - PBGMs/', dir, '/GAM vs Allele ', resp, '.RData'))
			# pred[ , i] <- predict(model, newdata=gadm@data, type='response')
			
		# }
		
		# gadmPred <- gadm
		# gadmPred@data[ , responses] <- pred
	
		# states <- sp::spTransform(states, getCRS('albersNA', TRUE))
	
		# countiesInMap <- crop(gadmPred, mapExtent)
		# statesInMap <- crop(states, mapExtent)
		# lopodInMap <- crop(lopod, mapExtent)
		
		# cols <- c('blue', 'red')
		# pal <- colorRampPalette(cols)
		# pal <- pal(101)

		# png(paste0('./Models - PBGMs/', dir, '/Maps of Candidate Allelic Frequency for 1970-2000 (Individual GAMs).png'), height=1800, width=2800, res=450)

			# par(mfrow=c(2, 2), cex=1.1, mar=0.2 * c(1, 2, 2, 10))

			# for (i in seq_along(responses)) {
				
				# resp <- responses[i]
				# say(resp)
				
				# ### predicted
				
				# # colors
				# vals <- as.data.frame(countiesInMap)[ , resp]
				# vals <- 1 + round(100 * vals)

				# # fade by probability of presence
				# prPres <- lopodInMap$psi_i
				# if (anyNA(prPres)) prPres[is.na(prPres)] <- 0
				# mapCols <- alpha(pal[vals], prPres)
				
				# plot(countiesInMap, col=mapCols, border=NA)
				# plot(statesInMap, col=NA, border='black', add=TRUE, lwd=0.5)

				# locus <- substr(resp, 10, 19)
				# if (locus == 'S1_2050023') {
					# allele <- '110'
					# altAllele <- '130'
				# } else if (locus == 'S1_27189028') {
					# allele <- '120'
					# altAllele <- '130'
				# } else if (locus == 'S1_33964042') {
					# allele <- '100'
					# altAllele <- '130'
				# } else if (locus == 'S1_73120031') {
					# allele <- '110'
					# altAllele <- '130'
				# }	
				
				# pos <- par('usr')
				# text(x=pos[1] + 0.01 * (pos[2] - pos[1]), y=pos[4] + 0.0 * (pos[4] - pos[3]), labels=paste0({letters}[i], ') Locus ', locus), xpd=NA, cex=0.45, xpd=NA, pos=4)

				# legendGrad('right', inset=c(-0.06, 0), width=0.13, height=0.92, labels=paste0('allele ', c(altAllele, allele), '\n(100%)'), labAdj=0.2, col=cols, gradAdjX=c(0.2, 0.4), gradAdjY=c(0, 1), cex=0.35, xpd=NA, boxBorder=NA, boxBg=NA, lwd=0.6)
				
				# ### observed
				# obs <- gadm[ , resp]
				# obs <- obs[-nas, ]
				# where <- coordinates(gCentroid(obs, byid=TRUE))
				# props <- cbind(1 - as.data.frame(obs), as.data.frame(obs))
				
				# for (i in 1:nrow(where)) pies(unlist(props[i, ]), add=TRUE, xPos=where[i, 1], yPos=where[i, 2], radius=40000, col=cols, aspect=FALSE, lwd=0.5)

				# title(sub=date(), outer=TRUE, line=-0.85, cex.sub=0.1)

			# } # next allele
			
		# dev.off()

# say('############################################################################################')
# say('### observed vs predicted GAM @ population level vs STRUCTURE groups & candidate alleles ###')
# say('#############################################################################################')

		# png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted.png', height=1400, width=2000)
		
			# layout(matrix(c(1:4, rep(5, 4)), ncol=4, byrow=TRUE))
			# par(cex=1.1)
		
			# ### pheno: observed vs expected
			
			# par(pty='s', cex=1.4)
			# for (resp in basePhenoResponses) {
				
				# niceName <- respNice(resp)

				# y <- model$inputs$y
				
				# # plot limits
				# resps <- c(y[ , resp], y[ , paste0(resp, 'Lower')], y[ , paste0(resp, 'Upper')], model$prediction$ypredMu[ , resp], model$prediction$ypredMu[ , paste0(resp, 'Lower')], model$prediction$ypredMu[ , paste0(resp, 'Upper')])
				# nice <- pretty(resps, 3)
				# lims <- range(nice)
				
				# # colors and R2
				# cols <- c('#1b7837', 'white', '#762a83')
				# r2 <- rep(NA, 3)
				# names(cols) <- names(r2) <- c('Upper', 'Mean', 'Lower')
				
				# # plot
				# for (type in c('Mean', 'Lower', 'Upper')) {
					
					# if (type == 'Mean') {
						# typeName <- ''
						# col <- cols[typeName]
					# } else {
						# typeName <- type
						# col <- alpha(cols[typeName], 0.6)
					# }
					
					# obs <- y[ , paste0(resp, typeName)]
					# pred <- model$prediction$ypredMu[ , paste0(resp, typeName)]
					
					# r2[type] <- cor(obs, pred)
				
					# if (type == 'Mean') {
						# plot(obs, pred, xlim=lims, ylim=lims, xlab=paste('Observed', niceName$unit), ylab=paste('Predicted', niceName$unit), main=niceName$nice, pch=21, cex=1.9, bg=col)
						# abline(0, 1)
						# points(obs, pred, pch=21, cex=1.9)
					# } else {
						# points(obs, pred, bg=col, pch=21, cex=1.9)
					# }
					
					
				# }
				
				# legend('bottomright', inset=0.02, legend=c('97.5% Quantile', 'Population Mean', '2.5% Quantile'), pt.bg=cols, pch=21, pt.cex=1.2)
				
				# # R2
				# box <- par('usr')
				# xPos <- box[1] + 0.02 * (box[2] - box[1])
				# yPos <- box[4] - 0.02 * (box[4] - box[3]) - 0.05 * (1:4) * (box[4] - box[3])
				# text(xPos, yPos, labels=paste(c('Pseudo-R2', 'Upper:', 'Mean:', 'Lower:'), c('', sprintf('%.2f', r2))), pos=4)
				
			# }
		
			# ### structure groups: left bar is observed, right is predicted
		
			# par(pty='m')
			# structObs <- t(y[ , structureGroupResponses])
			# structPred <- t(model$prediction$ypredMu[ , structureGroupResponses])

			# sites <- gadm@data$site[!is.na(gadm@data$site)]
			# obsPred <- matrix(NA, ncol=1, nrow=length(structureGroupResponses))
			# for (i in seq_along(sites)) obsPred <- cbind(obsPred, cbind(structObs[ , i], structPred[ , i]))
			# obsPred <- obsPred[ , -1]
			
			# for (i in 1:ncol(obsPred)) if (!anyNA(obsPred[ , i])) obsPred[ , i] <- obsPred[ , i] / sum(obsPred[ , i])

			# sites <- gadm@data$site[!is.na(gadm@data$site)]
			# cols <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a')
			# barplot(obsPred, ylim=c(0, 1), space=c(0.7, 0), yaxt='n', col=cols, main='STRUCTURE Groups')
			# axis(2, pos=-0.5)
			# text(-5, 0.5, xpd=NA, labels='Probability', srt=90, cex=1.1)
			# text(-0.5 + 2.7 * seq_along(sites), -0.02, labels=sites, xpd=NA, srt=90, adj=c(1, 0))
			# text(mean(-0.5 + 2.7 * seq_along(sites)), 1.03, labels='For each site the left bar is observed and the right is predicted. Missing observed bars were not genotyped.', xpd=NA, cex=1)
			# legend('bottom', inset=c(0.5, -0.19), xpd=NA, ncol=4, fill=cols, legend=paste('Group', 1:4), cex=0.85)

			# title(main=paste0('GJAM using 1st Two Axes of PCA on BIOCLIM Variables\n', date()), outer=TRUE, line=-3.5, cex.main=1.6)
			
		# dev.off()

		
say('DONE!!!!!!!!!!!', level=1)
