library(sp)
library(fields)
library(raster)
library(geometry)
library(RColorBrewer)

# Defining the parameters

samplingYears = 2012:2017
samplesPerYear = 35
nberOfIterations = 100000
nberOfRepetitions = 1000
cambodia = getData("GADM", country="KHM", level=0)
A = 0.9158

# Loading the data

colourScale = colorRampPalette(brewer.pal(11,"Spectral"))(11)[c(1,3,5,7,9,11)]
cities = c(); spatialPoints = c(); years = c(); cols = c(); coords = c()
for (i in 1:length(samplingYears))
	{
		csv = read.csv(paste0(samplingYears[i],"_UTM_sampling_frame.csv"), header=F)[2:4]
		colnames(csv) = c("city","x_coords","y_coords")
		cities = c(cities, as.character(csv$city))
		spatialPoints = rbind(spatialPoints, cbind(csv$x_coords,csv$y_coords))
		years = c(years, rep(samplingYears[i], dim(csv)[1]))
		cols = c(cols, rep(colourScale[i], dim(csv)[1]))
	}
spatialPoints_1 = data.frame(lon=as.numeric(spatialPoints[,1]), lat=as.numeric(spatialPoints[,2]))
UTM_48N = CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")
WGS84 = CRS("+init=epsg:4326")
coordinates(spatialPoints_1) = c("lon", "lat")
proj4string(spatialPoints_1) = UTM_48N
spatialPoints_2 = spTransform(spatialPoints_1, WGS84)
spatialPoints_1 = spatialPoints_1@coords
spatialPoints_2 = spatialPoints_2@coords
# plot(cambodia); points(spatialPoints_2, pch=3, col=cols, cex=0.8)
for (i in 1:dim(spatialPoints)[1])
	{
		coords = c(coords, paste0(spatialPoints_2[i,1],"_",spatialPoints_2[i,2]))
	}
allSamples = cbind(cities, spatialPoints_1, spatialPoints_2, years, coords, cols)
colnames(allSamples) = c("city","UTM_x","UTM_y","lon","lat","year","coords","col")

# Script of the subsampling selection (to become a function?)

	# 1. Computation of the pairwise great-circle distances between all sampling points

allPairwiseDistances = matrix(0, nrow=dim(allSamples)[1], ncol=dim(allSamples)[1])
for (i in 2:dim(allSamples)[1])
	{
		for (j in 1:(i-1))
			{
	 			x1 = cbind(as.numeric(allSamples[i,"lon"]), as.numeric(allSamples[i,"lat"]))
	 			x2 = cbind(as.numeric(allSamples[j,"lon"]), as.numeric(allSamples[j,"lat"]))
	 			d = rdist.earth(x1, x2, miles=F, R=NULL)
				allPairwiseDistances[i,j] = d; allPairwiseDistances[j,i] = d
			}
	}

	# 2. Starting the different repetitions (with different initial subsamplings)

for (r in 1:nberOfRepetitions)
	{
		selectedSamples0 = which(as.numeric(allSamples[,"year"])=="2017") # to select all the samples from 2017
		points = cbind(as.numeric(allSamples[,"lon"]), as.numeric(allSamples[,"lat"]))
		hull = chull(points); selectedSamples0 = c(selectedSamples0, hull) # to select all the most distant samples
		# hull = c(hull,hull[1]); plot(points); lines(points[hull,])
		for (i in 1:(length(samplingYears)-1))
			{
				indices = which(as.numeric(allSamples[,"year"])==samplingYears[i])
				remainingNberOfSamplesToSelect = samplesPerYear-sum(allSamples[selectedSamples0,"year"]==samplingYears[i])
				for (j in 1:remainingNberOfSamplesToSelect)
					{
						newCoords = FALSE
						while (newCoords == FALSE)
							{
								newSample = sample(indices, 1, replace=F)
								if (!allSamples[newSample,"coords"]%in%allSamples[selectedSamples0,"coords"])
									{
										newCoords = TRUE; indices = indices[which(indices!=newSample)]
										selectedSamples0 = c(selectedSamples0, newSample)
									}
							}
					}
			}
		selectedSamples = selectedSamples0
		# print(length(unique(allSamples[selectedSamples,"coords"])))
		selectedPoints = cbind(as.numeric(allSamples[selectedSamples,"lon"]),as.numeric(allSamples[selectedSamples,"lat"]))		
		triangulation = delaunayn(selectedPoints)
		adjacentsDistances1 = sum(allPairwiseDistances[triangulation[,1], triangulation[,2]])
		adjacentsDistances2 = sum(allPairwiseDistances[triangulation[,2], triangulation[,3]])
		adjacentsDistances3 = sum(allPairwiseDistances[triangulation[,1], triangulation[,3]])
		D0 = sum(adjacentsDistances1, adjacentsDistances2, adjacentsDistances3)
		Di_1 = D0; Dmax = D0
		if (r == 1)
			{
				# dir.create(file.path(map_directory), showWarnings=F)
				# files = list.files(map_directory)
				# if (length(files) > 0)
					# {
						# for (i in 1:length(files))
							# {
								# file.remove(paste0(map_directory,"/",files[i]))
							# }
					# }
				DMAX = Dmax; finalSamples = selectedSamples
			}
		# pdf(paste0(map_directory,"/Initial_subsampling_",r,"_D0=",round(D0,0),".pdf"), width=6, height=6); par(oma=c(0,0,0,0), mar=c(0,3,0,3))
		# plot(cambodia, lwd=0.1, border="gray10", col="gray90")
		# indices = which(!c(1:dim(allSamples)[1])%in%selectedSamples0)
		# points1 = cbind(as.numeric(allSamples[indices,"lon"]),as.numeric(allSamples[indices,"lat"]))
		# points(points1, pch=3, col="green3", cex=0.4, lwd=0.5)
		# indices = selectedSamples0
		# points2 = cbind(as.numeric(allSamples[indices,"lon"]),as.numeric(allSamples[indices,"lat"]))
		# points(points2, pch=3, col="red", cex=0.5, lwd=1)
		# title(main=paste0("Initial random subsampling ",r,"   -   D0 = ",round(D0,0)), cex.main=0.6, col.main="gray30", line=-3)
		# dev.off()

		# Starting the different iterations for a particular repetition r:
		for (i in 1:nberOfIterations)
			{
				new_selectedSamples = selectedSamples
				for (j in 1:(length(samplingYears)-1))
					{
						indices = which(allSamples[,"year"] == samplingYears[j])
						indices1 = indices[which(indices%in%selectedSamples)]
						indices2 = indices[which(!indices%in%selectedSamples)]
						indices1_without_hull = indices1[which(!indices1%in%hull)]
						newCoords = FALSE
						while (newCoords == FALSE)
							{
								index1 = sample(indices1_without_hull, 1, replace=F)
								index2 = sample(indices2, 1, replace=F)
								if (!allSamples[index2,"coords"]%in%allSamples[new_selectedSamples,"coords"])
									{
										newCoords = TRUE
										new_selectedSamples[new_selectedSamples[]==index1] = index2
									}
							}
					}
				# print(length(unique(allSamples[new_selectedSamples,"coords"])))
				selectedPoints = cbind(as.numeric(allSamples[new_selectedSamples,"lon"]),as.numeric(allSamples[new_selectedSamples,"lat"]))		
				triangulation = delaunayn(selectedPoints)
				adjacentsDistances1 = sum(allPairwiseDistances[triangulation[,1], triangulation[,2]])
				adjacentsDistances2 = sum(allPairwiseDistances[triangulation[,2], triangulation[,3]])
				adjacentsDistances3 = sum(allPairwiseDistances[triangulation[,1], triangulation[,3]])
				Di = sum(adjacentsDistances1, adjacentsDistances2, adjacentsDistances3)
				if (Dmax < Di)
					{
						Dmax = Di
						print(paste0("          new d_max = ",round(Dmax,0)," km  -  iteration ",i," of repetition ",r))
					}
				if (Di_1 <= Di)
					{
						selectedSamples = new_selectedSamples; Di_1 = Di
					}	else	{
						# SA = nberOfIterations^A
						SA = i^A
						p = exp(((Di/Dmax)-(Di_1/Dmax))*SA)
						# print(paste0("p = ",round(p,4)))
						R = runif(1,0,1)
						if (R < p)
							{
								selectedSamples = new_selectedSamples; Di_1 = Di
							}
					}
			}
		print(paste0("     final d_max = ",round(Dmax,0)," for repetition ",r))
		if (DMAX < Dmax)
			{
				DMAX = Dmax; finalSamples = selectedSamples
				print(paste0("New overall D_max = ",round(Dmax,0)," for repetition ",r))
			}
	}
write.csv(allSamples[finalSamples,c("city","UTM_x","UTM_y","lon","lat","year")], "Final_samples_selection_NEW.csv", quote=F, row.names=F)
pdf(paste0("Final_samples_selection_NEW.pdf"), width=6, height=6); par(oma=c(0,0,0,0), mar=c(0,3,0,3))
plot(cambodia, lwd=0.1, border="gray10", col="gray90")
indices = which(!c(1:dim(allSamples)[1])%in%finalSamples)
points1 = cbind(as.numeric(allSamples[indices,"lon"]),as.numeric(allSamples[indices,"lat"]))
points(points1, pch=3, col="green3", cex=0.4, lwd=0.5)
indices = finalSamples
points2 = cbind(as.numeric(allSamples[indices,"lon"]),as.numeric(allSamples[indices,"lat"]))
points(points2, pch=3, col="red", cex=0.5, lwd=1)
title(main=paste0("Final random subsampling   -   Dmax = ",round(DMAX,0)), cex.main=0.6, col.main="gray30", line=-3)
dev.off()
