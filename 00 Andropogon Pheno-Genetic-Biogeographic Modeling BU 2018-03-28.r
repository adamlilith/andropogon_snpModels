### ANDROPOGON PHENO-GENETIC-BIOGEOGRAPHICAL MODELING
### Adam B Smith | Missouri Botanical Garden | 2018
### source('C:/ecology/Drive/Research/Andropogon/Scripts - PGBMs/00 Andropogon Pheno-Genetic-Biogeographic Modeling.r')

print(Sys.time())
rm(list=ls())
memory.limit(memory.limit() * 2^30)
gc()

library(omnibus)
library(legendary)
library(enmSdm)
library(bayesLopod)
library(scales)
library(rgeos)
library(raster)
library(dismo)
library(gjam)
library(openxlsx)

options(stringsAsFactors=FALSE)
rasterOptions(format='GTiff', overwrite=TRUE)

setwd('C:/ecology/Drive/Research/Andropogon')

################
### CONTENTS ###
################

### constants ###
### functions ###
### copy and process raw response data ###
### create data frame with rows representing INDIVIDUALS ###
### create data frame with rows representing POPULATIONS ###
### create study region ###
### collate environmental data ###
### extract environmental data ###
### merge shapefile with environmental data with shapefile with specimen data ###
### merge occurrence data and population-level phenotype and genotype data ###
### model species distribution using bayesLopod ###A
### GJAM @ population level ###
### plot GJAM @ population level ###


#################
### constants ###
#################

predictors <- c('wc02', 'wc05', 'wc07', 'wc12', 'solarRadiation')
		
# names of candidate loci
loci <- c('S1_2050023', 'S1_27189028', 'S1_33964042', 'S1_73120031')

#################
### functions ###
#################

	# nice name and unit of each response variable
	respNice <- function(x) {

		if (x == 'biomassG') {
			nice <- 'Biomass'
			unit <- ' (g)'
		} else if (x == 'biomassGLower') {
			nice <- 'Biomass (Lowest 2.5%)'
			unit <- ' (g)'
		} else if (x == 'biomassGUpper') {
			nice <- 'Biomass (Upper 97.5%)'
			unit <- ' (g)'
		} else if (x == 'leafWidthCm') {
			nice <- 'Leaf Width'
			unit <- ' (cm)'
		} else if (x == 'leafWidthCmLower') {
			nice <- 'Leaf Width (Lowest 2.5%)'
			unit <- ' (cm)'
		} else if (x == 'leafWidthCmUpper') {
			nice <- 'Leaf Width (Upper 97.5%)'
			unit <- ' (cm)'
		} else if (x == 'heightCm') {
			nice <- 'Height'
			unit <- ' (cm)'
		} else if (x == 'heightCmLower') {
			nice <- 'Height (Lowest 2.5%)'
			unit <- ' (cm)'
		} else if (x == 'heightCmUpper') {
			nice <- 'Height (Upper 97.5%)'
			unit <- ' (cm)'
		} else if (x == 'spad') {
			nice <- 'SPAD'
			unit <- ''
		} else if (x == 'spadLower') {
			nice <- 'SPAD (Lowest 2.5%)'
			unit <- ''
		} else if (x == 'spadUpper') {
			nice <- 'SPAD (Upper 97.5%)'
			unit <- ''
		} else if (x == 'freqLocusS1_2050023allele110') {
			nice <- 'Locus S1_2050023 (Allele 110)'
			unit <- ''
		} else if (x == 'freqLocusS1_27189028allele120') {
			nice <- 'Locus S1_27189028 (Allele 120)'
			unit <- ''
		} else if (x == 'freqLocusS1_33964042allele100') {
			nice <- 'Locus S1_33964042 (Allele 100)'
			unit <- ''
		} else if (x == 'freqLocusS1_73120031allele110') {
			nice <- 'S1_73120031 (Allele 110)'
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

		if (x == 'wc01') {
			nice <- 'WC01'
			name <- 'Mean Annual Temperature'
			unit <- ' (deg C)'
		} else if (x == 'wc02') {
			nice <- 'WC02'
			name <- 'Diurnal Temperature Range'
			unit <- ' (deg C)'
		} else if (x == 'wc05') {
			nice <- 'WC05'
			name <- 'Maximum Temperature'
			unit <- ' (deg C)'
		} else if (x == 'wc07') {
			nice <- 'WC07'
			name <- 'Temperature Range'
			unit <- ' (deg C)'
		} else if (x == 'wc12') {
			nice <- 'WC12'
			name <- 'Mean Annual Precipitation'
			unit <- ' (mm)'
		} else if (x == 'solarRadiation') {
			nice <- 'Solar Radiation'
			name <- 'Solar Radiation'
			unit <- ' (MJ)'
		}
		
		list(nice=nice, name=name, unit=unit)
		
	}

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
	
	# # want one data frame for individuals with fields re STRUCTURE-assigned population group, phenotypic measurements, vegetation class of habitat

	# ind <- read.csv('Pheno-Geno-Biogeographic Data/00 By Individual - IDs for ALL Individuals [ORIG gbs2014 phenotype and genotype samples].csv')
	# ind$gbs2014_phenotype_sort <- ind$gbs2014_filterCS_pop_sort <- ind$kuch <- ind$gbs_SEQUENCE <- NULL
	# names(ind)[names(ind) == 'gbs_name_2'] <- 'individual'
	
	# ind$site[ind$site == 'NE-4' | ind$site == 'NE-5'] <- 'NE-4/NE-5'
	# ind$site[ind$site == 'NE-3' | ind$site == 'NE-7'] <- 'NE-3/NE-7'

	# ### match individuals to STRUCTURE groups
	# #########################################
		
		# # population-level genetic data
		# fromData <- read.csv('Pheno-Geno-Biogeographic Data/00 By Individual - Neutral Genetic Data [ORIG gbs2014_filterCS STRUCTURE analysis groupings for ADAM].csv')
		# fromData$individual <- paste0(substr(fromData$individual, 1, 2), '-', substr(fromData$individual, 3, 4), '-', substr(fromData$individual, 5, 6))
		
		# # are all individuals in genetic data frame also in master individuals data frame?
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

	# ### match individuals to PHENOTYPE (only individuals that were also GENOTYPED)
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
	
	# ### match individuals that were not genotyped to PHENOTYPE
	# ##########################################################
	
		# say('NOTE: As of yet there is no way to match individuals that were not genotyped with their phenotype because there is no common ID between data sets!', level=2)
		
	# ### match individuals to Kuchler vegetation groups
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
# say('### collate environmental data ###')
# say('##################################')

	# dirCreate('./WORLDCLIM Ver 2 Rel June 1 2016')
	
	# mask <- raster('./Masks/mask_northAmerica_10arcmin.tif')
	
	# # gadm <- readRDS('./Extents/GADM Ver 2 - Dissolved to Level 0 - Multipart - WGS84 - Clipped to Projection Region.rds')
	
	# elev <- raster('F:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/10 arcmin/Elevation - 10 arcmin/elevation.tif')
	# elev <- crop(elev, mask)
	# elev <- elev * mask
	# names(elev) <- 'elevation'
	
	# writeRaster(elev, './WORLDCLIM Ver 2 Rel June 1 2016/elevation')
	
	# wcs <- c(1, 2, 5, 7, 12)

	# ### CURRENT (1970-2010)
	# #######################

		# dirCreate('WORLDCLIM Ver 2 Rel June 1 2016/1970-2000')
		
		# ### BIOCLIMS
		# ############
		
		# for (wc in wcs) {
			
			# rast <- raster(paste0('F:/ecology/Climate/WORLDCLIM Ver 2 Rel June 1 2016/10 arcmin 1970 to 2000/wc2.0_bio_10m_', prefix(wc, 2), '.tif'))
			# rast <- crop(rast, mask)
			# rast <- rast * mask
			# rast <- setMinMax(rast)
			# names(rast) <- paste0('wc', prefix(wc, 2))
			# writeRaster(rast, paste0('./WORLDCLIM Ver 2 Rel June 1 2016/1970-2000/wc', prefix(wc, 2)))
			
		# }
		
		# ### SOLAR RADIATION
		# ###################
		
		# srad <- raster::stack(listFiles('F:/ecology/Climate/WORLDCLIM Ver 2 Rel June 1 2016/10 arcmin 1970 to 2000', pattern='srad'))
		# rast <- sum(srad)
		# rast <- crop(rast, mask)
		# rast <- rast * mask
		# rast <- setMinMax(rast)
		# names(rast) <- 'solarRadiation'
		# writeRaster(rast, paste0('WORLDCLIM Ver 2 Rel June 1 2016/1970-2000/solarRadiation'))
		
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
	
# say('###################################################')
# say('### model species distribution using bayesLopod ###')
# say('###################################################')

	# load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')
	
	# gadm <- gadm[-which(is.na(gadm$poaRec)), ]
	
	# gadmLopod <- shapeLopodData(gadm, fieldN='poaRec', fieldY='anyAg1to3', Adjacency=TRUE, keepFields=TRUE)
	
	# modelVarP <- modelLopod(gadmLopod, varP=TRUE, CAR=TRUE)
	# modelConstantP <- modelLopod(gadmLopod, varP=FALSE, CAR=TRUE)
	
	# dirCreate('./Models - PBGMs/bayesLopod with Variable Detection')
	# dirCreate('./Models - PBGMs/bayesLopod with Constant Detection')
	
	# save(lopodVarP, file='./Models - PBGMs/bayesLopod with Variable Detection/LOPOD Model.Rdata')
	# save(lopodConstantP, file='./Models - PBGMs/bayesLopod with Constant Detection/LOPOD Model.Rdata')
	
	# lopodPredictVarP <- lopodShape(modelVarP, param='psi_i', extrapolate=TRUE, metric='mean')
	# lopodPredictConstantP <- lopodShape(modelConstantP, param='psi_i', extrapolate=TRUE, metric='mean')
	
	# save(lopodPredictVarP, file='./Models - PBGMs/bayesLopod with Variable Detection/LOPOD Prediction.Rdata')
	# save(lopodPredictConstantP, file='./Models - PBGMs/bayesLopod with Constant Detection/LOPOD Prediction.Rdata')
	
say('###############################')
say('### GJAM @ population level ###')
say('###############################')

	dirCreate('./Models - PBGMs/GJAM @ Population Level')

	load('./Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.RData')

	### responses
		
		basePhenoResponses <- c('biomassG', 'leafWidthCm', 'heightCm', 'spad')
		phenoResponses <- c('biomassGLower', 'biomassG', 'biomassGUpper', 'leafWidthCmLower', 'leafWidthCm', 'leafWidthCmUpper', 'heightCmLower', 'heightCm', 'heightCmUpper', 'spadLower', 'spad', 'spadUpper')
		structureGroupResponses <- c('structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')
		# responses <- c(phenoResponses, structureGroupResponses)
		responses <- c(phenoResponses, structureGroupResponses)
		
		typeNames <- c(rep('CA', length(phenoResponses)), rep('FC', length(structureGroupResponses)))
		
		y <- gadm@data[!is.na(gadm@data$site), responses]

	### effort
			
		effort <- NA * y
		for (resp in responses) {
			
			effortCol <- if (resp %in% c('biomassG', 'biomassGLower', 'biomassGUpper')) {
				'nBiomass'
			} else if (resp %in% c('leafWidthCm', 'leafWidthCmLower', 'leafWidthCmUpper')) {
				'nLeafWidth'
			} else if (resp %in% c('heightCm', 'heightCmLower', 'heightCmUpper')) {
				'nHeight'
			} else if (resp %in% c('spad', 'spadLower', 'spadUpper')) {
				'nSpad'
			} else if (resp %in% c('freqLocusS1_2050023allele110', 'freqLocusS1_27189028allele120', 'freqLocusS1_33964042allele100', 'freqLocusS1_73120031allele110', 'structureGroup1', 'structureGroup2', 'structureGroup3', 'structureGroup4')) {
				'geneticSampleSize'
			}
			
			effort[!is.na(y[ , resp]), resp] <- gadm@data[!is.na(gadm@data[ , effortCol]), effortCol]
			if (anyNA(effort[ , resp])) effort[is.na(effort[ , resp]), resp] <- 0
			
		}
		
		effort <- list(columns=names(effort), values=as.matrix(effort))

	### predictors
			
		x <- gadm@data[!is.na(gadm@data$site), predictors]
		pca <- princomp(x, cor=TRUE)
		# biplot(pca)
		x <- as.data.frame(pca$scores)
		names(x) <- paste0('pc', prefix(seq_along(predictors), 2))
		save(pca, file='./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.Rdata')

	### model: PHENOTYPE and STRUCTURE
			
		form <- as.formula(~ pc01 + pc02 + pc01:pc02 + I(pc01^2) + I(pc02^2))
		
		tuning <- list(ng = 5000, burnin = 1000, typeNames = typeNames, effort = effort)

		model <- gjam(form, xdata = x, ydata = y, modelList = tuning)
		summary(model)
		
		save(model, file='./Models - PBGMs/GJAM @ Population Level/GJAM on Phenotype & STRUCTURE Groups.RData')
		
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

		# # png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted.png', height=1400, width=2000)
		
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
		
		# # load('./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.Rdata')
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
		
		# load('./Models - PBGMs/PCA on Environmental Predictors at Sampled Populations.Rdata')
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
		
		# png(paste0('./Models - PBGMs/GJAM @ Population Level/STRUCTURE Maps for 1970-2000.png'), height=900, width=1400, res=450)

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
	
	# png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted in Environmental Space.png', height=1800, width=1600, res=300)

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
	
# print(NON)
# say('### model diagnostics: STRUCTURE in environmental space')
# # source('C:/ecology/Drive/Research/Andropogon/Scripts - PGBMs/TEMP.r')

	# png('./Models - PBGMs/GJAM @ Population Level/Observed vs Predicted STRUCTURE Population in Environmental Space.png', height=1800, width=1600, res=300)

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
	
	
	
say('DONE!!!!!!!!!!!', level=1)
