library(diagram)
library(lubridate)
library(exactextractr)
library(maptools)
library(phylotate)
library(seraphim)

source("MCC_tree_extraction.r")

e_Cambodia_1 = extent(101, 109, 9, 16)
e_Cambodia_2 = extent(101.5, 108, 10, 15)

# 8. Testing the impact of environmental factors on lineage dispersal velocity

analysis = "N_genes"; localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900

showingPlots = FALSE; nberOfCores = 10; OS = "Unix"; randomisations = FALSE
fourCells = FALSE; nberOfRandomisations = 0; randomProcedure = 3
c = 0; envVariables = list(); resistances = list(); avgResistances = list()
envVariableNames = c("Land_cover_forests_RABV_CA_04","Land_cover_savannas_RABV_CA_04","Land_cover_grasslands_RABV_CA_04","Land_cover_croplands_RABV_CA_04",
"Land_cover_water_RABV_CA_04","Human_pop_density_RABV_CA_04","Elevation_RABV_CA_04","Annual_mean_temperature_RABV_CA_04","Annual_precipitation_RABV_CA_04")
for (k in c(10,100,1000))
	{
		for (i in 1:length(envVariableNames))
			{
				c = c+1
				rast = raster(paste("Environmental_rasters/",envVariableNames[i],".asc",sep=""))
				rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
				names(rast) = paste(envVariableNames[i],"_k",k,sep="")
				envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableNames[i],"_k",k,sep="")
				resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
			}
		for (i in 1:length(envVariableNames))
			{
				c = c+1
				rast = raster(paste("Environmental_rasters/",envVariableNames[i],".asc",sep=""))
				rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
				names(rast) = paste(envVariableNames[i],"_k",k,sep="")
				envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableNames[i],"_k",k,sep="")
				resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
			}
	}
# pathModel = 2; outputName = paste0("All_seraphim_results/",analysis,"_LC_extractions")
# spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  # nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)
pathModel = 3; outputName = paste0(analysis,"_CS_extractions")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)	
# pathModel = 2; outputName = paste0("All_seraphim_results/",analysis,"_LC_simulations")
# spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  # nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	
# pathModel = 3; outputName = paste0(analysis,"_CS_simulations")
# spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  # nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	

