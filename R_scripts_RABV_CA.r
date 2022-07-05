library(adephylo)
library(diagram)
library(geometry)
library(lubridate)
library(exactextractr)
library(maptools)
library(phylotate)
library(seraphim)

# 1. Generating the figure based on overall maximum likelihood RABV tree
# 2. Selection of samples to sequence through a Markov chain procedure
# 3. Generating the figure based on the discrete phylogeographic analysis
# 4. Testing the correlation between the patristic and geographic distances
# 5. Preparing the environmental rasters for the landscape phylogeographic analyses
# 6. Plotting the environmental rasters for the landscape phylogeographic analyses
# 7. Extracting the spatio-temporal information embedded in posterior and MCC trees
# 8. Plotting the dispersal history of RABV lineages in Cambodia (for both analyses)
# 9. Estimating dispersal statistics based on continuous phylogeographic analyes
# 10. Continuous phylogeographic reconstruction for RABV in Tanzania (Brunker et al. 2018)
# 11. Comparing the weighted lineage dispersal velocity estimates with other datasets
# 12. Generating a null dispersal model for the landscape phylogeographic analyses
# 13. Testing the impact of environmental factors on lineage dispersal locations
# 14. Testing the impact of environmental factors on lineage dispersal velocity

wd = getwd(); source("MCC_tree_extraction.r")
e_Cambodia_1 = extent(101, 109, 9, 16)
e_Cambodia_2 = extent(101.5, 108, 10, 15)
purple = "#8151A1"; orange = "#FAA51A"

# 1. Generating the figure based on overall maximum likelihood RABV tree

tre = read.nexus("Global_RABV_ML.tree"); tab = read.csv("Global_RABV_ML.csv", head=T, sep=";")
locations = c("Cambodia","Southeast Asia","Asia","Caribbean","Central America","Europe",
			  "Middle East","North Africa","North America","South America","Sub-Saharan Africa")
colours = c(purple,orange,colorRampPalette(brewer.pal(12,"Paired"))(12)[c(1:8,12)])
cols = rep(NA, length(tre$tip.label))
for (i in 1:length(tre$tip.label))
	{
		index1 = which(tab[,"full_name"]==gsub("'","",tre$tip.label[i]))
		if (length(index1) != 1)
			{
				print(c(i))
			}	else	{
				index2 = which(locations==tab[index1,"continent"])
				cols[i] = colours[index2]
				if (tab[index1,"country"] == "Cambodia")
					{
						cols[i] = colours[1]
					}
			}
	}
pdf(paste0("Global_RABV_M_NEW.pdf"), width=7, height=6.5); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1, col="gray30")
plot(tre, type="fan", show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray50", edge.color="gray50")
for (i in 1:dim(tre$edge)[1])
	{
		if (!tre$edge[i,2]%in%tre$edge[,1])
			{
				nodelabels(node=tre$edge[i,2], pch=16, cex=0.65, col=cols[tre$edge[i,2]])
				nodelabels(node=tre$edge[i,2], pch=1, cex=0.65, col="gray50", lwd=0.5)
			}
	}
add.scale.bar(x=0.00, y=-0.015, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
legend(-0.1, 0.0, locations, col=colours, text.col="gray30", pch=16, pt.cex=1.3, box.lty=0, cex=0.75, y.intersp=1.1)
dev.off()

# 2. Selection of samples to sequence through a Markov chain procedure

setwd(paste0(wd,"/Selection_of_samples"))

	# 2.1. Defining the parameters

samplingYears = 2012:2017
samplesPerYear = 35
nberOfIterations = 100000
nberOfRepetitions = 1000
cambodia = getData("GADM", country="KHM", level=0)
A = 0.9158

	# 2.2. Loading the data

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
for (i in 1:dim(spatialPoints)[1])
	{
		coords = c(coords, paste0(spatialPoints_2[i,1],"_",spatialPoints_2[i,2]))
	}
allSamples = cbind(cities, spatialPoints_1, spatialPoints_2, years, coords, cols)
colnames(allSamples) = c("city","UTM_x","UTM_y","lon","lat","year","coords","col")

	# 2.3. Computation of the pairwise great-circle distances between all sampling points

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

	# 2.4. Starting the different repetitions (with different initial subsamplings)

for (r in 1:nberOfRepetitions)
	{
		selectedSamples0 = which(as.numeric(allSamples[,"year"])=="2017") # to select all the samples from 2017
		points = cbind(as.numeric(allSamples[,"lon"]), as.numeric(allSamples[,"lat"]))
		hull = chull(points); selectedSamples0 = c(selectedSamples0, hull) # to select all the most distant samples
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
		selectedPoints = cbind(as.numeric(allSamples[selectedSamples,"lon"]),as.numeric(allSamples[selectedSamples,"lat"]))		
		triangulation = delaunayn(selectedPoints)
		adjacentsDistances1 = sum(allPairwiseDistances[triangulation[,1], triangulation[,2]])
		adjacentsDistances2 = sum(allPairwiseDistances[triangulation[,2], triangulation[,3]])
		adjacentsDistances3 = sum(allPairwiseDistances[triangulation[,1], triangulation[,3]])
		D0 = sum(adjacentsDistances1, adjacentsDistances2, adjacentsDistances3)
		Di_1 = D0; Dmax = D0
		if (r == 1)
			{
				DMAX = Dmax; finalSamples = selectedSamples
			}	# (starting the different iterations for a particular repetition r:)
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
						SA = i^A
						p = exp(((Di/Dmax)-(Di_1/Dmax))*SA)
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
selectedSamples = read.csv("Final_samples_selection.csv", head=T)
points1 = cbind(as.numeric(allSamples[,"lon"]),as.numeric(allSamples[,"lat"]))
points(points1, pch=3, col="gray50", cex=0.4, lwd=0.5)
points2 = cbind(as.numeric(selectedSamples[,"lon"]),as.numeric(selectedSamples[,"lat"]))
points(points2, pch=3, col=purple, cex=0.5, lwd=1)
title(main=paste0("Final random subsampling   -   Dmax = 41846453"), cex.main=0.6, col.main="gray30", line=-3)
dev.off()

# 3. Generating the figure based on the discrete phylogeographic analysis

setwd(wd); mcc_tre = readAnnotatedNexus("BEAST_DTA_SEA.tree"); rootHeight = max(nodeHeights(mcc_tre))

	# 3.1. Comparing patristic distances between sequences sampled within the same country or not

times = as.matrix(distTips(mcc_tre, method="patristic"))
times1 = times[which(!grepl("Cambodia",row.names(times))),which(!grepl("Cambodia",colnames(times)))]
times2 = times[which(grepl("Cambodia",row.names(times))),which(grepl("Cambodia",colnames(times)))]
embeddedCambodianSequences = c("EU086170_Canis-lupus-familiaris_1997_Cambodia_NA_NA",
							   "KM366268_Dog_2005-09-07_Cambodia-PhnomPenh-MeanChey-PrekPra-OuAndaung1_11.492_104.952",
							   "KM366301_Dog_2009-08-19_Cambodia-KgChhnang-RoleaPhaAir-RoleaPhaAir-PreyKhmer_12.165_104.665")
times3 = times2[which(!gsub("'","",row.names(times2))%in%embeddedCambodianSequences),
			    which(!gsub("'","",colnames(times2))%in%embeddedCambodianSequences)]
times1a = c(); times1b = c(); countries = c(); times2 = times2[lower.tri(times2)]; times3 = times3[lower.tri(times3)]
for (i in 2:dim(times1)[1])
	{
		name1 = row.names(times1)[i]; country1 = unlist(strsplit(name1,"_"))[4]
		country1 = unlist(strsplit(country1,"-"))[1]; countries = c(countries, country1)
		for (j in 1:(i-1))
			{
				name2 = colnames(times1)[j]; country2 = unlist(strsplit(name2,"_"))[4]
				country2 = unlist(strsplit(country2,"-"))[1]; countries = c(countries, country2)
				if (country1 == country2) times1a = c(times1a, times1[i,j])
				if (country1 != country2) times1b = c(times1b, times1[i,j])
			}
	} # countries = unique(countries)
pdf(paste0("BEAST_DTA_SEA_NEW1.pdf"), width=8, height=2.2); # dev.new(width=7, height=2.2)
par(mar=c(1.5,1.7,0.5,0.5), oma=c(0,0,0,0), mgp=c(0,0.3,0), lwd=0.2, bty="o", col="gray30")
plot(density(times2), col=NA, xlim=c(0,700), ylim=c(0,0.025), axes=F, ann=F, main=NA)
polygon(density(times3), col=paste0(purple,"30"), border=NA)
polygon(density(times1a), col=paste0(orange,"30"), border=NA)
polygon(density(times1b), col=paste0("#4D4D4D30"), border=NA)
lines(density(times3), col=purple, lwd=1); lines(density(times1a), col=orange, lwd=1); lines(density(times1b), col="gray30", lwd=1)
axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.025, col.axis="gray30", mgp=c(0,0.04,0))
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.025, col.axis="gray30", mgp=c(1,0.30,0))
dev.off()

	# 3.2. Plotting the MCC tree of the discrete phylogeographic analysis
	
tipLabels = mcc_tre$tip.label; dates = rep(NA, length(mcc_tre$tip.label))
for (i in 1:length(tipLabels))
	{
		text = unlist(strsplit(tipLabels[i],"_"))[3]
		if (length(unlist(strsplit(text,"-"))) == 1) dates[i] = as.numeric(text)
		if (length(unlist(strsplit(text,"-"))) == 2) dates[i] = decimal_date(ymd(paste0(text,"-15")))
		if (length(unlist(strsplit(text,"-"))) == 3) dates[i] = decimal_date(ymd(text))
	}
mostRecentSamplingDatum = max(dates); tMRCA = mostRecentSamplingDatum-rootHeight
purple = "#8151A1"; orange = "#FAA51A"; cols = rep(orange,dim(mcc_tre$edge)[1])
for (i in 1:dim(mcc_tre$edge)[1])
	{
		if (mcc_tre$edge[i,1]%in%mcc_tre$edge[,2])
			{
				index = which(mcc_tre$edge[,2]==mcc_tre$edge[i,1])
				if ((mcc_tre$annotations[[index]]$location=="Cambodia")&(mcc_tre$annotations[[i]]$location=="Cambodia")) cols[i] = purple
			}
	}
pdf(paste0("BEAST_DTA_SEA_NEW2.pdf"), width=10, height=3.5); # dev.new(width=10, height=3.5)
par(mar=c(1,1,1,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o")
plot(mcc_tre, direction="downwards", show.tip.label=F, show.node.label=F, edge.width=0.45, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
for (i in 1:dim(mcc_tre$edge)[1])
	{
		if ((!mcc_tre$edge[i,2]%in%mcc_tre$edge[,1])&(tree$annotations[[i]]$location != "Cambodia"))
			{
				nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.5, col="white")
				nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.5, col=paste0(orange,"90"))
				nodelabels(node=mcc_tre$edge[i,2], pch=1, cex=0.5, col="gray30", lwd=0.2)
			}
		if ((!mcc_tre$edge[i,2]%in%mcc_tre$edge[,1])&(tree$annotations[[i]]$location == "Cambodia"))
			{
				nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.5, col="white")
				nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.5, col=paste0(purple,"90"))
				nodelabels(node=mcc_tre$edge[i,2], pch=1, cex=0.5, col="gray30", lwd=0.2)
			}
	}
selectedDates = c(1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050); selectedLabels = selectedDates
axis(lwd=0.2, at=mostRecentSamplingDatum-selectedDates, labels=selectedLabels, cex.axis=0.6, mgp=c(0,-0.8,-1), lwd.tick=0.2,
	 col.lab="gray30", col.axis="gray30", col.ticks="gray30", col="gray30", tck=-0.012, side=2)
dev.off()

	# 3.3. Age of the MRCA of the main Cambodian clade

age = mostRecentSamplingDatum-48.892 # obtained through FigTree
hpd95 = mostRecentSamplingDatum-c(57.360,40.906) # obtained through FigTree
cat(round(age,1),", 95% HPD = [",round(hpd95[1],1),"-",round(hpd95[2],1),"]",sep="")

# 4. Testing the correlation between the patristic and geographic distances
	
 	# H0: sequences closer to the Thai border do not tend to be more related to Thai sequences
 	# H1: sequences closer to the Thai border tend to be more related to Thai sequences
	# Spearman correlation computed between patristic distances and the geographic distance
	# (between each tip and the closest point of the Thai border)
	# --> H0: for 1000 posterior trees for which the geographic coordinates assigned to the tips are swapped
	# --> H1: for the same 1000 posterior trees but without swapping the sampling coordinates among tips

trees = read.nexus("BEAST_DTA_SEA.trees")
trees = trees[(ceiling(0.1*length(trees))+1):length(trees)]
interval = floor(length(trees)/1000)
trees = trees[seq(interval,1000*interval,interval)]
embeddedCambodianSequences = c() # if also including the Cambodian sequences embedded in the other clades
embeddedCambodianSequences = c("EU086170_Canis-lupus-familiaris_1997_Cambodia_NA_NA",
							   "KM366268_Dog_2005-09-07_Cambodia-PhnomPenh-MeanChey-PrekPra-OuAndaung1_11.492_104.952",
							   "KM366301_Dog_2009-08-19_Cambodia-KgChhnang-RoleaPhaAir-RoleaPhaAir-PreyKhmer_12.165_104.665")
fourTargetCambodianSequences = c("KM366298_Dog_2002-03-16_Cambodia-BanteayMeanchey-SereySophorn-ToeukThla-KampongSvay_13.592_102.976",
							   "KM366291_Dog_2003-12-22_Cambodia-BanteayMeanchey-ThmorPuok-ThmorPuok-ThmorPuok_13.939_103.058",
							   "KM366271_Dog_2001-09-25_Cambodia-SiemReap-AngloungVeng-TrapaingPrey-PrambeyRuoy_14.226_103.942",
							   "KM366209_Dog_2002-05-16_Cambodia-UdoorMeanchey-AnloungVeng-AnloungVeng-TuolKandal_14.229_104.104")
								# --> note: those 4 sequences are not in basal position in the MCC tree
fourTargetCambodianSequences = c() # if including all Cambodian sequences
rS_inferred = rep(NA, length(trees)); rS_swapped = rep(NA, length(trees))
cambodia = getData("GADM", country="KHM", level=0)
thailand = getData("GADM", country="THA", level=0)
coords = c(); geoDistsObs = c(); IDs = c()
for (j in 1:length(trees[[1]]$tip.label))
	{
		if ((grepl("Cambodia",trees[[1]]$tip.label[j]))&(!trees[[1]]$tip.label[j]%in%embeddedCambodianSequences)&(!trees[[1]]$tip.label[j]%in%fourTargetCambodianSequences))
			{
				minDist = 9999
				x = as.numeric(unlist(strsplit(trees[[1]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[1]]$tip.label[j],"_")))])
				y = as.numeric(unlist(strsplit(trees[[1]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[1]]$tip.label[j],"_")))-1])
				if ((!is.na(x))&(!is.na(y)))
					{
						x1 = cbind(x,y); coords = rbind(coords, x1)
						for (k in 1:length(thailand@polygons[[1]]@Polygons))
							{
								x2 = thailand@polygons[[1]]@Polygons[[k]]@coords
								d = min(rdist.earth(x1, x2, miles=F, R=NULL))
								if (d < minDist) minDist = d
							}
						geoDistsObs = c(geoDistsObs, minDist)
						IDs = c(IDs, trees[[1]]$tip.label[j])
					}
			}
	}
row.names(coords) = IDs; names(geoDistsObs) = IDs
for (i in 1:length(trees))
	{
		print(i)
		geoDists = c(); patDists = c(); IDs = c(); times = as.matrix(distTips(trees[[i]], method="patristic"))
		for (j in 1:length(trees[[i]]$tip.label))
			{
				if ((grepl("Cambodia",trees[[i]]$tip.label[j]))&(!trees[[i]]$tip.label[j]%in%embeddedCambodianSequences)&(!trees[[i]]$tip.label[j]%in%fourTargetCambodianSequences))
					{
						x = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))])
						y = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))-1])
						if ((!is.na(x))&(!is.na(y)))
							{
								index = which(names(geoDistsObs)==trees[[i]]$tip.label[j])
								geoDists = c(geoDists, geoDistsObs[index])
								index = which(row.names(times)==trees[[i]]$tip.label[j])
								indices = which(grepl("Thailand",colnames(times)))
								patDists = c(patDists, mean(times[index,indices]))
								IDs = c(IDs, trees[[1]]$tip.label[j])
							}
					}
			}
		rS_inferred[i] = cor(geoDists, patDists, method="spearman")
		if (i == 1)
			{
				pdf(paste0("BEAST_DTA_SEA_NEW3.pdf"), width=6, height=3); # dev.new(width=6, height=3)
				par(mfrow=c(1,2), oma=c(0,1,0,1), mar=c(0,0,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30")
				cambodia = gSimplify(getData("GADM", country="KHM", level=0), 0.01)
				plot(cambodia, col="gray90", border="gray70", lwd=0.75)
				cols = paste0(colorRampPalette(brewer.pal(11,"BrBG"))(121)[20:121],"90")
				cols = cols[(((patDists-min(patDists))/(max(patDists)-min(patDists)))*100)+1]
				dS = 1/patDists; cexs = (((dS-min(dS))/(max(dS)-min(dS)))*3)+0.3; positions = coords[IDs,]
				for (j in 1:dim(positions)[1])
					{
						points(positions[j,1], positions[j,2], pch=16, cex=cexs[j], col=cols[j])
						points(positions[j,1], positions[j,2], pch=1, cex=cexs[j], col="gray30", lwd=0.2)
					}
				# indices = which(cexs>2.3); text(positions[indices,], labels=IDs[indices], cex=cexs[indices]/10, col="gray30")
				rast = raster(as.matrix(c(min(patDists),max(patDists))))
				plot(rast, legend.only=T, add=T, col=colorRampPalette(brewer.pal(11,"BrBG"))(121)[20:121], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.84,0.85,0.25,0.40), adj=3,
					 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.7, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.3,0)), alpha=1, side=3)
				collectionDates = rep(NA, length(patDists)); c = 0
				for (j in 1:length(trees[[i]]$tip.label))
					{
						if ((grepl("Cambodia",trees[[i]]$tip.label[j]))&(!trees[[1]]$tip.label[j]%in%embeddedCambodianSequences)&(!trees[[i]]$tip.label[j]%in%fourTargetCambodianSequences))
							{
								x = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))])
								y = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))-1])
								if ((!is.na(x))&(!is.na(y)))
									{
										text = unlist(strsplit(trees[[i]]$tip.label[j],"_"))[3]; c = c+1
										if (length(unlist(strsplit(text,"-"))) == 1) collectionDates[c] = as.numeric(text)
										if (length(unlist(strsplit(text,"-"))) == 2) collectionDates[c] = decimal_date(ymd(paste0(text,"-15")))
										if (length(unlist(strsplit(text,"-"))) == 3) collectionDates[c] = decimal_date(ymd(text))
									}
							}
					}
				plot(cambodia, col="gray90", border="gray70", lwd=0.75)
				cols = paste0(colorRampPalette(brewer.pal(11,"PuOr"))(101)[1:101],"90")
				cols = cols[(((collectionDates-min(collectionDates))/(max(collectionDates)-min(collectionDates)))*100)+1]
				dS = 1/patDists; cexs = (((dS-min(dS))/(max(dS)-min(dS)))*3)+0.3; positions = coords[IDs,]
				for (j in 1:dim(positions)[1])
					{
						points(positions[j,1], positions[j,2], pch=16, cex=cexs[j], col=cols[j])
						points(positions[j,1], positions[j,2], pch=1, cex=cexs[j], col="gray30", lwd=0.2)
					}
				rast = raster(as.matrix(c(min(collectionDates),max(collectionDates))))
				plot(rast, legend.only=T, add=T, col=colorRampPalette(brewer.pal(11,"PuOr"))(121)[20:121], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.84,0.85,0.25,0.40), adj=3,
					 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.7, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.3,0)), alpha=1, side=3)
				dev.off() # plot(collectionDates, patDists) # --> almost a perfect correlation... (*)
			}
		geoDists = c(); patDists = c(); geoDistsRan = geoDistsObs
		names(geoDistsRan) = sample(names(geoDistsRan), length(geoDistsRan), replace=F)
		for (j in 1:length(trees[[i]]$tip.label))
			{
				if ((grepl("Cambodia",trees[[i]]$tip.label[j]))&(!trees[[i]]$tip.label[j]%in%embeddedCambodianSequences)&(!trees[[i]]$tip.label[j]%in%fourTargetCambodianSequences))
					{
						x = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))])
						y = as.numeric(unlist(strsplit(trees[[i]]$tip.label[j],"_"))[length(unlist(strsplit(trees[[i]]$tip.label[j],"_")))-1])
						if ((!is.na(x))&(!is.na(y)))
							{
								index = which(names(geoDistsRan)==trees[[i]]$tip.label[j])
								geoDists = c(geoDists, geoDistsRan[index])
								index = which(row.names(times)==trees[[i]]$tip.label[j])
								indices = which(grepl("Thailand",colnames(times)))
								patDists = c(patDists, mean(times[index,indices]))
							}
					}
			}
		rS_swapped[i] = cor(geoDists, patDists, method="spearman")
	}
p = sum(rS_inferred[!is.na(rS_inferred)]>rS_swapped[!is.na(rS_swapped)])/length(rS_inferred[!is.na(rS_inferred)]); BF = (p/(1-p))/(0.5/(1-0.5))
	# BF = 141.9 when including all Cambodian sequences, and BF = ~15.7 when discarding embedded and targetted Cambodian sequences
	# Interpretation: those results mostly reflect an impact of the heterogeneous sampling dates as most sequences close to the Thai border are among the oldest (*)

# 5. Preparing the environmental rasters for the landscape phylogeographic analyses

source("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/Seraphim/Divers_R/Divers_R/decreaseResolution.r")
source("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/Seraphim/Divers_R/Divers_R/landCoverRasters.r")
elevation = crop(raster("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/GIS_files/Original_rasters/Elevation.tif"), e_Cambodia_1)
writeRaster(elevation, "Environmental_rasters/Elevation_RABV_CA_008.asc", overwrite=T)
writeRaster(decreaseResolution(elevation, 5), "Environmental_rasters/Elevation_RABV_CA_04.asc", overwrite=T)
temperature = crop(raster("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/GIS_files/WorldClim_2.0/WC2.0_bio_2.5m_01.tif"), e_Cambodia_1)
writeRaster(temperature, "Environmental_rasters/Annual_mean_temperature_RABV_CA_04.asc", overwrite=T)
precipitation = crop(raster("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/GIS_files/WorldClim_2.0/WC2.0_bio_2.5m_12.tif"), e_Cambodia_1)
writeRaster(precipitation, "Environmental_rasters/Annual_precipitation_RABV_CA_04.asc", overwrite=T)
population = crop(raster("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/GIS_files/WorldPop_files/Population_count_2020_1km.tif"), e_Cambodia_1)
writeRaster(decreaseResolution(population, 5), "Environmental_rasters/Human_pop_density_RABV_CA_04.asc", overwrite=T)
land_cover = crop(raster("/Users/simondellicour/Dropbox/Simon_repos/Projects/Viruses/GIS_files/Original_rasters/Cover_land.tif"), e_Cambodia_1)
landCoverRasters(land_cover, 5, "RABV_CA_04") # selected variables: croplands, forests, grasslands, savannas, water

lakes = crop(shapefile("All_natural_Earth_files/Natural_Earth_lakes.shp"), e_Cambodia_1)
files = list.files("Environmental_rasters"); files = files[which((grepl(".asc",files))&(!grepl("water",files)))]
for (i in 1:length(files))
	{
		rast = raster(paste0("Environmental_rasters/",files[i]))
		pol = lakes@polygons[[1]]@Polygons[[1]] # because only one polygon
		p = Polygon(pol@coords); ps = Polygons(list(p),1)
		sps = SpatialPolygons(list(ps)); pol = sf::st_as_sfc(sps)
		extraction = exact_extract(rast, pol, include_cell=T)[[1]]
		extraction = extraction[which(extraction[,"coverage_fraction"]>0.5),]
		rast[extraction[,"cell"]] = NA
		writeRaster(rast, paste0("Environmental_rasters/",files[i]), overwrite=T)
	}

# 6. Plotting the environmental rasters for the landscape phylogeographic analyses

borders = crop(shapefile("All_natural_Earth_files/International_borders.shp"), e_Cambodia_1)
envVariableFiles = c("Land_cover_forests","Land_cover_savannas","Land_cover_grasslands","Land_cover_croplands","Land_cover_water",
					 "Human_pop_density","Elevation","Annual_mean_temperature","Annual_precipitation")
envVariableNames1 = c("Forests","Savannas","Grasslands","Croplands","Water","Human pop.","Elevation","Annual","Annual")
envVariableNames2 = c("","","","","","density (log10)","","mean temp.","precipitation")
rS = list(); cols = list(); croppingRaster = FALSE; colour1 = "white"
for (i in 1:length(envVariableFiles)) rS[[i]] = raster(paste0("Environmental_rasters/",envVariableFiles[i],"_RABV_CA_04.asc"))
rS[[6]][!is.na(rS[[6]][])] = log10(rS[[6]][!is.na(rS[[6]][])]+1); rS[[7]][rS[[7]][]<0] = 0
rS[[9]][] = rS[[9]][]/100 # legend: temperature in °C and precipitation in meters
colour2 = "olivedrab3"; r = rS[[1]]
cols[[1]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour2 = "burlywood3"; colour2 = "brown"; r = rS[[2]]
cols[[2]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour2 = "chartreuse4"; r = rS[[3]]
cols[[3]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour2 = "navajowhite4"; r = rS[[4]] 
cols[[4]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour2 = colorRampPalette(brewer.pal(9,"Blues"))(9)[5]; r = rS[[5]]
cols[[5]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
cols[[6]] = c("#FFFFFF",colorRampPalette(brewer.pal(9,"BuPu"))(100))
cols[[7]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
cols[[8]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(140)[1:100]
cols[[9]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
xmin = 102; xmax = 108; ymin = 10; ymax = 15
labsX = c(expression(102*degree*E), expression(108*degree*E))
labsY = c(expression(10*degree*N), expression(15*degree*N))

pdf("Environmental_factors_NEW.pdf", width=7.3, height=6.7)
par(mfrow=c(3,3), oma=c(2,2.5,1,0.0), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2, col="gray30"); colNA="#D0D0D0"
for (i in 1:length(rS))
	{
		plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA=colNA)
		axis(1,c(xmin,xmax),labels=labsX,pos=ymin(rS[[i]]),col="gray30",cex.axis=0.6,col.axis="gray30",lwd=0,lwd.tick=0.2,col.tick="gray30",tck=-0.02,mgp=c(0,0.10,0))
		axis(2,c(ymin,ymax),labels=labsY, pos=xmin(rS[[i]]),col="gray30",cex.axis=0.6,col.axis="gray30",lwd=0,lwd.tick=0.2,col.tick="gray30",tck=-0.02,mgp=c(0,0.27,0))
		plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.040,0.055,0.12,0.35), adj=3,
		axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		plot(borders, col="gray30", lwd=0.5, add=T)
		if (nchar(envVariableNames1[i] > 0)) mtext(envVariableNames1[i], side=3, adj=0, line=-2.4, at=101.3, cex=0.65, font=1, col="gray30")
		if (nchar(envVariableNames2[i] > 0)) mtext(envVariableNames2[i], side=3, adj=0, line=-3.3, at=101.3, cex=0.65, font=1, col="gray30")
		rect(xmin(rS[[i]]), ymin(rS[[i]]), xmax(rS[[i]]), ymax(rS[[i]]), lwd=0.2, border="gray30")
	}
dev.off()

# 7. Extracting the spatio-temporal information embedded in posterior and MCC trees

localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900; mostRecentSamplingDatum = decimal_date(ymd("2017-12-20"))
allTrees = scan(paste0("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.trees"), what="", sep="\n", quiet=T)
burnIn = 101; randomSampling = FALSE; nberOfTreesToSample = nberOfExtractionFiles; coordinateAttributeName = "gps"; nberOfCores = 5
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
mcc_tre = readAnnotatedNexus(paste0("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.tree"))
mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, paste0("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.csv"), row.names=F, quote=F)

localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900; mostRecentSamplingDatum = decimal_date(ymd("2017-12-20"))
allTrees = scan(paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.trees"), what="", sep="\n", quiet=T)
burnIn = 101; randomSampling = FALSE; nberOfTreesToSample = nberOfExtractionFiles; coordinateAttributeName = "gps"; nberOfCores = 5
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
mcc_tre = readAnnotatedNexus(paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.tree"))
mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.csv"), row.names=F, quote=F)

pdf(paste0("BEAST_DTA_SEA_NEW4.pdf"), width=3.5, height=3.5); # dev.new(width=3.5, height=3.5)
par(oma=c(0,0,0,0.0), mar=c(0,3,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30")
cambodia = gSimplify(getData("GADM", country="KHM", level=0), 0.01); purple = "#8151A1"
topography = mask(crop(raster("All_Natural_Earth_files/Gray_background.tif"), cambodia), cambodia); r = topography
lakes = crop(shapefile("All_natural_Earth_files/Natural_Earth_lakes.shp"), e_Cambodia_1)
pol = lakes@polygons[[1]]@Polygons[[1]]; p = Polygon(pol@coords); ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps)); pol = sf::st_as_sfc(sps); extraction = exact_extract(r, pol, include_cell=T)[[1]]
extraction = extraction[which(extraction[,"coverage_fraction"]>0.5),]; r[extraction[,"cell"]] = NA; topography = r
cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
mcc_tab = read.csv(paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.csv"), head=T)
plot(topography, col=cols, ann=F, frame=F, box=F, axes=F, legend=F); plot(cambodia, col=NA, border="gray70", lwd=0.75, add=T)
points(mcc_tab[which(!mcc_tab[,"node2"]%in%mcc_tab[,"node1"]),c("endLon","endLat")], cex=0.4, pch=1, col="gray30", lwd=0.2)
points(mcc_tab[which(!mcc_tab[,"node2"]%in%mcc_tab[,"node1"]),c("endLon","endLat")], cex=0.4, pch=16, col=paste0(purple,"90"))
dev.off()

pdf(paste0("BEAST_DTA_SEA_NEW5.pdf"), width=3.5, height=3.5); # dev.new(width=3.5, height=3.5)
par(oma=c(0,0,0,0.0), mar=c(0,3,0,0), mgp=c(0,0,0), lwd=0.2, col="gray30")
cambodia = gSimplify(getData("GADM", country="KHM", level=0), 0.01)
human_pop = mask(crop(raster(paste0("Environmental_rasters/Human_pop_density_RABV_CA_04.asc")), cambodia), cambodia)
human_pop[!is.na(human_pop[])] = log10(human_pop[!is.na(human_pop[])]+1); cols = colorRampPalette(brewer.pal(9,"Greys"))(121)[1:101]
plot(human_pop, col=cols, ann=F, frame=F, box=F, axes=F, legend=F); plot(cambodia, col=NA, border="gray70", lwd=0.75, add=T)
plot(human_pop, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.74,0.75,0.25,0.40), adj=3,
	 axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.7, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.3,0)), alpha=1, side=3)
dev.off()

# 8. Plotting the dispersal history of RABV lineages in Cambodia (for both analyses)

localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900
mcc = read.csv(paste0("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.csv"), head=T)
bsg = read.table(paste0("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.txt"), head=T)

localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900
mcc = read.csv(paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.csv"), head=T)
bsg = read.table(paste0("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.txt"), head=T)

mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]; mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)

	# 8.1. Loading the different GIS files used for the graphic

background = crop(raster("Environmental_rasters/Elevation_RABV_CA_008.asc"), e_Cambodia_2)
lakes = crop(shapefile("All_natural_Earth_files/Natural_Earth_lakes.shp"), e_Cambodia_2)
borders = crop(shapefile("All_Natural_Earth_files/International_borders.shp"), e_Cambodia_2)
coastLines = crop(shapefile("All_Natural_Earth_files/Coastline_borders.shp"), e_Cambodia_2)

	# 8.2. Estimating the 80% HPD regions for successive time slices

prob = 0.80; startDatum = min(mcc[,"startYear"]); precision = 5
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

	# 8.3. Defining the colour scale (for tree nodes and HPD regions)

cols = gsub("FF","",viridis::viridis(101)[1:101]); minYearColours = min(mcc[,"startYear"])
endYearsM = ((mcc[,"endYear"]-minYearColours)/(max(mcc[,"endYear"])-minYearColours)*100)+1; endYearsM[endYearsM<1] = 1
cols_mcc = cols[endYearsM]; col_start = cols[1]; legend = raster(as.matrix(cbind(0,0)))
legend[1] = minYearColours; legend[2] = max(mcc[,"endYear"]); cols_pol = list()
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		pol_index = round((((date-minYearColours)/(max(mcc[,"endYear"])-minYearColours))*100)+1)
		if (pol_index < 1) pol_index = 1
		cols_pol[[i]] = paste0(cols[pol_index],"20")
	}

	# 8.4. Generating and saving the spread graphic in a PDF format

dev.new(width=6.5, height=4.7); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(0,1,0,0), mgp=c(1,0.2,0), lwd=0.3)
plot(background, main="", cex.main=0.8, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="gray90", colNA="white")
plot(borders, lwd=0.5, add=T, col="gray50"); plot(lakes, add=T, border=NA, col="white")
for (i in 1:length(polygons))
	{
		plot(polygons[[i]], axes=F, col=cols_pol[[i]], add=T, border=NA)
	}
for (i in dim(mcc)[1]:1)
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.15, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=cols_mcc[i], cex=0.5)
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.5, lwd=0.15)
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=col_start, cex=0.5)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=0.5, lwd=0.15)
			}
	}
rect(e_Cambodia_2@xmin, e_Cambodia_2@ymin, e_Cambodia_2@xmax, e_Cambodia_2@ymax, xpd=T, lwd=0.5, border="gray30")
plot(legend, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.880,0.890,0.042,0.956), adj=3,
	 axis.args=list(cex.axis=0.55, lwd=0, lwd.tick=0.5, tck=-0.6, col="gray30", col.lab="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)

	# 8.5. Generating and saving the skygrid reconstruction graphs

bsg = bsg[which(bsg[,"Time"]>=minYear),]; bsg[,"Median"] = log10(bsg[,"Median"])
bsg[,"Mean"] = log10(bsg[,"Mean"]); bsg[,"Lower"] = log10(bsg[,"Lower"]); bsg[,"Upper"] = log10(bsg[,"Upper"])

dev.new(width=6.5, height=3.5); par(oma=c(1,1,0,0), mar=c(0,1,0,0), mgp=c(1,0.2,0), lwd=0.3)
minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
xMin = min(bsg[,"Time"]); xMax = max(bsg[,"Time"]); xMin = minYear; timeSlice = bsg[1,"Time"]-bsg[2,"Time"]
yMin = min(log10(bsg[,"Lower"])); yMax = max(log10(bsg[,"Upper"]))
yMin = 0.5; yMax = 3.0 # for the skygrid based on the full genomes
yMin = -0.3; yMax = 3.2 # for the skygrid based on the N genes
colours = cols[(((bsg[,c("Time")]-minYear)/(maxYear-minYear))*100)+1]
plot(bsg[,"Time"], bsg[,"Median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30", axes=F, xlab=NA, ylab=NA, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
xx_l = c(bsg[,c("Time")],rev(bsg[,c("Time")])); yy_l = c(bsg[,"Lower"],rev(bsg[,"Upper"]))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=rgb(187/255,187/255,187/255,0.25),border=0)
for (i in 1:length(bsg[,"Time"]))
	{
		x1 = bsg[i,"Time"]-(timeSlice/2); x2 = bsg[i,"Time"]+(timeSlice/2)
		y1 = bsg[i,"Lower"]-0.2; y2 = bsg[i,"Upper"]+0.2
		polygon(c(x1,x2,x2,x1), c(y1,y1,y2,y2), col=colours[i], border=NA)
	}
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l,yy_l,col=NA,border="gray30")
lines(bsg[,"Time"], bsg[,"Median"], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col="gray30")
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col="gray30", col.axis="gray30", mgp=c(0,0.04,0), at=seq(1950,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.017, col="gray30", col.axis="gray30", mgp=c(1,0.30,0), at=seq(-1,4,1))
mtext("log10(Ne)", side=2, col="gray30", cex=0.80, line=1.2, las=3)

# 9. Estimating dispersal statistics based on continuous phylogeographic analyes

localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900
timeSlices = 50; onlyTipBranches = FALSE; showingPlots = FALSE; outputName = "All_dispersal_statistics/Genomes"; nberOfCores = 5; slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	# median value of mean branch dispersal velocity = 16.1 km/year, 95% HPD = [13.2-22.4]
	# median value of weighted branch dispersal velocity = 6.7 km/year, 95% HPD = [6.1-7.4]
	# median value of original diffusion coefficient = 114.9 km2/year, 95% HPD = [96.0-144.9]
	# median value of weighted diffusion coefficient = 89.9 km2/year, 95% HPD = [78.0, 106.8]

localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900
timeSlices = 50; onlyTipBranches = FALSE; showingPlots = FALSE; outputName = "All_dispersal_statistics/N_genes"; nberOfCores = 5; slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	# median value of mean branch dispersal velocity = 30.2 km/year, 95% HPD = [23.5-61.8]
	# median value of weighted branch dispersal velocity = 12.6 km/year, 95% HPD = [11.0-14.3]
	# median value of original diffusion coefficient = 603.5 km2/year, 95% HPD = [400.7-1877.3]
	# median value of weighted diffusion coefficient = 367.5 km2/year, 95% HPD = [319.8-418.5]

# 10. Continuous phylogeographic reconstruction for RABV in Tanzania (Brunker et al. 2018)

localTreesDirectory = "Tree_extraction_files/Tanzania"; nberOfExtractionFiles = 900; mostRecentSamplingDatum = 2013.67945205479
allTrees = scan(paste0("BEAST_RRW_analysis/Tanzania/RABV_Brunker_et_al_Tanzania_j001.trees"), what="", sep="\n", quiet=T)
burnIn = 101; randomSampling = FALSE; nberOfTreesToSample = nberOfExtractionFiles; coordinateAttributeName = "gps"; nberOfCores = 5
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
mcc_tre = readAnnotatedNexus(paste0("BEAST_RRW_analysis/Tanzania/RABV_Brunker_et_al_Tanzania_j001.tree"))
mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, paste0("BEAST_RRW_analysis/Tanzania/RABV_Brunker_et_al_Tanzania_j001.csv"), row.names=F, quote=F)

localTreesDirectory = "Tree_extraction_files/Tanzania"; nberOfExtractionFiles = 900
mcc = read.csv(paste0("BEAST_RRW_analysis/Tanzania/RABV_Brunker_et_al_Tanzania_j001.csv"), head=T)
mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]; mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
tanzania = gSimplify(getData("GADM", country="TZA", level=0), 0.01)
background = mask(raster("Environmental_rasters/Elevation_RABV_TA_008.asc"), tanzania)
lakes = crop(shapefile("All_natural_Earth_files/Natural_Earth_lakes.shp"), tanzania)
prob = 0.80; startDatum = min(mcc[,"startYear"]); precision = 5
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
cols = gsub("FF","",viridis::viridis(101)[1:101]); minYearColours = min(mcc[,"startYear"])
endYearsM = ((mcc[,"endYear"]-minYearColours)/(max(mcc[,"endYear"])-minYearColours)*100)+1; endYearsM[endYearsM<1] = 1
cols_mcc = cols[endYearsM]; col_start = cols[1]; legend = raster(as.matrix(cbind(0,0)))
legend[1] = minYearColours; legend[2] = max(mcc[,"endYear"]); cols_pol = list()
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		pol_index = round((((date-minYearColours)/(max(mcc[,"endYear"])-minYearColours))*100)+1)
		if (pol_index < 1) pol_index = 1
		cols_pol[[i]] = paste0(cols[pol_index],"20")
	}

pdf("BEAST_RRW_analysis/Tanzania/RABV_Brunker_et_al_Tanzania_j001.pdf", width=6, height=5.2)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(0.5,0.5,0.5,0), mgp=c(1,0.2,0), lwd=0.3)
plot(background, main="", cex.main=0.8, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="gray90", colNA="white")
for (i in 1:length(polygons))
	{
		plot(polygons[[i]], axes=F, col=cols_pol[[i]], add=T, border=NA)
	}
for (i in dim(mcc)[1]:1)
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.1, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.02, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=cols_mcc[i], cex=0.1)
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.1, lwd=0.02)
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=col_start, cex=0.1)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=0.1, lwd=0.02)
			}
	}
plot(legend, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.900,0.910,0.042,0.956), adj=3,
	 axis.args=list(cex.axis=0.55, lwd=0, lwd.tick=0.5, tck=-0.6, col="gray30", col.lab="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
dev.off()

# 11. Comparing the weighted lineage dispersal velocity estimates with other datasets

nberOfExtractionFiles = 900; maximumDistance = 1000; registerDoMC(cores=1)
datasets = c("Genomes","N_genes","Previous/RABV_AF_900t","Tanzania","Previous/RABV_IR_900t","Previous/RABV_YU_900t")
colours1 = rep(NA, length(datasets)); colours2 = rep(NA, length(datasets))
colours1[1] = rgb(150,150,150,255,maxColorValue=255); colours2[1] = rgb(150,150,150,100,maxColorValue=255) # light grey
colours1[1] = rgb(76,76,76,255,maxColorValue=255); colours2[1] = rgb(60,60,60,100,maxColorValue=255) # dark grey
colours1[2] = rgb(129,81,161,255,maxColorValue=255); colours2[2] = rgb(129,81,161,100,maxColorValue=255) # purple
colours1[3] = rgb(250,165,26,255,maxColorValue=255); colours2[3] = rgb(250,165,26,100,maxColorValue=255) # orange
colours1[4] = rgb(139,101,8,255,maxColorValue=255); colours2[4] = rgb(139,101,8,100,maxColorValue=255) # brown
colours1[5] = rgb(222,67,39,255,maxColorValue=255); colours2[5] = rgb(222,67,39,100,maxColorValue=255) # red
colours1[6] = rgb(70,118,187,255,maxColorValue=255); colours2[6] = rgb(70,118,187,100,maxColorValue=255) # blue
for (i in 1:length(datasets))
	{
		localTreesDirectory = paste0("Tree_extraction_files/",datasets[i]); buffer = list()
		for (t in 1:nberOfExtractionFiles)
			{
				data = read.csv(paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep=""), h=T)
				dists = rdist.earth(data[,c("startLon","startLat")], data[,c("endLon","endLat")], miles=F)
				distances = diag(dists)
				if ((i == 1)&(t == 1))
					{
						minDistance = min(distances); maxDistance = max(distances)
					}	else	{
						if (minDistance > min(distances)) minDistance = min(distances)
						if (maxDistance < max(distances)) maxDistance = max(distances)
					}
			}
	}
minDistance = 0; maxDistance = 1000; interval = 5; wldvs_list1 = list(); wldvs_list2 = list()
for (h in 1:length(datasets))
	{
		localTreesDirectory = paste0("Tree_extraction_files/",datasets[h]); wldvs_temp = list()
		for (i in 1:3)
			{
				buffer = foreach(t = 1:nberOfExtractionFiles) %dopar% {				
				# for (t in 1:nberOfExtractionFiles) {
						data = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",t,".csv"), h=T)
						dists = rdist.earth(data[,c("startLon","startLat")], data[,c("endLon","endLat")], miles=F)
						distances = diag(dists); durations = data[,"endYear"]-data[,"startYear"]
						wldv = matrix(nrow=maxDistance/interval, ncol=2)
						wldv[,1] = seq(interval, maxDistance, interval)
						for (j in 1:dim(wldv)[1])
							{
								if (i == 1) # WLDV list for the maximum geographic distances
									{
										indices = which(distances<=wldv[j,1])
									}
								if (i == 2) # WLDV list for the sliding window of 20 km
									{
										indices = which((distances>=(wldv[j,1]-10))&(distances>=(wldv[j,1]+10)))
									}
								if (i == 3) # WLDV list for the sliding window of 50 km
									{
										indices = which((distances>=(wldv[j,1]-25))&(distances>=(wldv[j,1]+25)))
									}
								wldv[j,2] = sum(distances[indices])/sum(durations[indices])
							}
						colnames(wldv) = c("maxDist", "wldv")
						wldv
					}
				wldvs = matrix(nrow=maxDistance/interval, ncol=nberOfExtractionFiles)
				row.names(wldvs) = seq(interval, maxDistance, interval)
				for (t in 1:length(buffer)) wldvs[,t] = buffer[[t]][,2]
				wldvs_temp[[i]] = wldvs
			}
		wldvs_list1[[h]] = wldvs_temp
	}
for (h in 1:length(wldvs_list1))
	{
		wldvs_temp = list()
		for (i in 1:length(wldvs_list1[[h]]))
			{
				median = matrix(nrow=dim(wldvs_list1[[h]][[i]])[1], ncol=1)
				lower = matrix(nrow=dim(wldvs_list1[[h]][[i]])[1], ncol=1)
				upper = matrix(nrow=dim(wldvs_list1[[h]][[i]])[1], ncol=1)
				for (j in 1:dim(wldvs_list1[[h]][[i]])[1])
					{
						quantiles = quantile(wldvs_list1[[h]][[i]][j,], probs=c(0.025,0.975), na.rm=T)
						quantiles = quantile(wldvs_list1[[h]][[i]][j,], probs=c(0.1,0.9), na.rm=T)
						median[j,1] = median(wldvs_list1[[h]][[i]][j,], na.rm=T)
						lower[j,1] = quantiles[1]; upper[j,1] = quantiles[2]
					}
				tab = matrix(nrow=dim(wldvs_list1[[h]][[i]])[1], ncol=4)
				colnames(tab) = c("maxDistance","median","lower_95HPD","higher_95HPD")
				tab[,1] = as.numeric(row.names(wldvs_list1[[h]][[i]]))
				tab[,2] = median[,1]; tab[,3] = lower[,1]; tab[,4] = upper[,1]
				wldvs_temp[[i]] = tab
			}
		wldvs_list2[[h]] = wldvs_temp
	}
pdf(paste0("WLDV_distances_NEW.pdf"), width=7.3, height=2.5); H = 1
par(oma=c(0,0,0,0), mar=c(1.5,3,0,1), mgp=c(1,0.2,0), lwd=0.2, col="gray30")
layout(matrix(c(1,1,1,1,1,1,1,1,1,2,3,4,2,3,4), ncol=5)); cutOffs = c(50, 100, 200)
plottingDashedLinesForHPDintervals = FALSE; indices = c(6,5,4,3,2,1); c = 0
for (i in indices)
	{
		c = c+1; tab = wldvs_list2[[i]][[H]]
		if (c == 1)
			{
				if (H == 1) plot(tab[,c("maxDistance","median")], xlim=c(30,750), ylim=c(1.3,33), col=NA, axes=F, ann=F)
				if (H == 2) plot(tab[,c("maxDistance","median")], xlim=c(30,550), ylim=c(1.3,500), col=NA, axes=F, ann=F)
				if (H == 3) plot(tab[,c("maxDistance","median")], xlim=c(30,550), ylim=c(1.3,500), col=NA, axes=F, ann=F)
			}
		xx_l = c(tab[,c("maxDistance")],rev(tab[,c("maxDistance")]))
		yy_l = c(tab[,"lower_95HPD"],rev(tab[,"higher_95HPD"]))
		getOption("scipen"); opt = options("scipen"=20)
		polygon(xx_l, yy_l, col=colours2[i], border=0)
	}
for (i in indices)
	{
		if (plottingDashedLinesForHPDintervals == TRUE)
			{
				tab = wldvs_list2[[i]][[H]]
				lines(tab[,c("maxDistance","lower_95HPD")], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col=colours1[i], lty=2)
				lines(tab[,c("maxDistance","higher_95HPD")], lwd=0.7, type="l", cex.axis=0.8, cex.lab=0.8, col=colours1[i], lty=2)
			}
	}
for (i in indices)
	{
		tab = wldvs_list2[[i]][[H]]
		lines(tab[,c("maxDistance","median")], lwd=1.5, type="l", cex.axis=0.8, cex.lab=0.8, col=colours1[i])
	}
for (i in 1:length(cutOffs)) abline(v=cutOffs[i], col="gray60", lwd=0.75, lty=2)
axis(side=1, lwd.tick=0.2, cex.axis=0.70, lwd=0.2, tck=-0.017, col="gray30", col.axis="gray30", mgp=c(0,0.20,0), at=seq(0,800,100))
axis(side=2, lwd.tick=0.2, cex.axis=0.70, lwd=0.2, tck=-0.020, col="gray30", col.axis="gray30", mgp=c(1,0.40,0), at=seq(-5,35,5))
title(ylab="weighted lineage dispersal velocity (km/year)", cex.lab=0.9, mgp=c(1.7,0,0), col.lab="gray30")
title(xlab="geographic distance (km)", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
for (h in 1:length(cutOffs))
	{
		c = 0
		for (i in indices)
			{
				c = c+1
				vS = wldvs_list1[[i]][[H]][which(as.numeric(row.names(wldvs_list1[[i]][[H]]))==cutOffs[h]),]
				if (c == 1) plot(density(vS), xlim=c(0,25), ylim=c(0,1.2), col=NA, axes=F, ann=F, mar=c(3,2,0,1))
				polygon(density(vS), border=NA, col=colours2[i])
			}
		for (i in indices)
			{	
				vS = wldvs_list1[[i]][[H]][which(as.numeric(row.names(wldvs_list1[[i]][[H]]))==cutOffs[h]),]
				lines(density(vS), lwd=1.0, col=colours1[i])
			}
		axis(side=2, lwd.tick=0.2, cex.axis=0.70, lwd=0.2, tck=-0.06, col="gray30", col.axis="gray30", mgp=c(1,0.40,0), at=seq(-0.5,1.5,0.5))
		title(ylab="density", cex.lab=0.9, mgp=c(1.7,0,0), col.lab="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.70, lwd=0.2, tck=-0.06, col="gray30", col.axis="gray30", mgp=c(0,0.20,0), at=seq(-5,30,5))
		if (h == 3) title(xlab="weighted lineage dispersal velocity (km/year)", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
		mtext(paste0("< ",cutOffs[h]," km"), side=3, line=-2.0, at=23.5, cex=0.6)
	}
dev.off()

# 12. Generating a null dispersal model for the landscape phylogeographic analyses

analysis = "Genomes"; localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900
log = read.table("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.log", header=T)[102:1001,]
trees = readAnnotatedNexus("BEAST_RRW_analysis/Genomes/Compiled_Genomes_aligned_gamma.trees")[102:1001]

analysis = "N_genes"; localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900
log = read.table("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.log", header=T)[102:1001,]
trees = readAnnotatedNexus("BEAST_RRW_analysis/N_genes/Compiled_N_genes_aligned_gamma.trees")[102:1001]

if (!file.exists(paste0("Environmental_rasters/Minimum_convex_hull_",analysis,"_04.asc")))
	{
		points1 = c(); points2 = c()
		for (i in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
				points1 = rbind(points1, tab[,c("startLon","startLat")])
				points2 = rbind(points2, tab[,c("endLon","endLat")])
			}
		colnames(points1) = c("lon","lat"); colnames(points2) = c("lon","lat")
		points = rbind(points1, points2); jitter = 0
		if (jitter > 0)
			{
				points1 = points; points2 = points; points3 = points; points4 = points
				points1[,1] = points1[,1]-jitter; points2[,1] = points2[,1]+jitter
				points3[,2] = points3[,2]-jitter; points4[,2] = points4[,2]+jitter
				points = rbind(points, points1, points2, points3, points4)
			}
		hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
		ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		rast1 = raster("Environmental_rasters/Elevation_RABV_CA_04.asc")
		rast1[!is.na(rast1[])] = 0; rast2 = rasterize(sps, rast1, getCover=T, fun=avg)
		rast2[rast2[]==0] = NA; rast2[is.na(rast1[])] = NA; backgroundRaster = rast2
		writeRaster(backgroundRaster, paste0("Environmental_rasters/Minimum_convex_hull_",analysis,"_04.asc"), overwrite=T)
	}	else	{
		backgroundRaster = raster(paste0("Environmental_rasters/Minimum_convex_hull_",analysis,"_04.asc"))
	}
showingPlots = TRUE; nodesOnly = FALSE; newPlot = TRUE
showingPlots = FALSE; nodesOnly = FALSE; newPlot = FALSE
model = "gamma"; n1 = 100; n2 = 10; mostRecentSamplingDatum = 0
envVariables = list(backgroundRaster)
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"))
		if (mostRecentSamplingDatum < max(tab[,"endYear"])) mostRecentSamplingDatum = max(tab[,"endYear"])
	}
for (i in 1:nberOfExtractionFiles)
	{
		if (!file.exists(paste(localTreesDirectory,"/TreeSimulations_",i,".csv",sep=""))) {
		tree = trees[[i]]				
		tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		ancestID = which(!tab[,"node1"]%in%tab[,"node2"])[1]
		ancestPosition = c(tab[ancestID,"startLon"], tab[ancestID,"startLat"])
		rates = c(); geoDists = matrix(nrow=dim(tab)[1], ncol=1)
		for (j in 1:length(tree$annotations))
			{
				rates = c(rates, tree$annotations[[j]]$gps.rate)
			}
		for (j in 1:dim(tab)[1])
			{
				x1 = cbind(tab[j,"startLon"], tab[j,"startLat"])
				x2 = cbind(tab[j,"endLon"], tab[j,"endLat"])
				geoDists[j,1] = rdist.earth(x1, x2, miles=F, R=NULL)
			}
		col11 = log[i,"treeLengthPrecision1"]
		col12 = log[i,"treeLengthPrecision3"]
		col22 = log[i,"treeLengthPrecision2"]
		my_prec = c(col11, col12, col12, col22)
		halfDF = log[i,"gps.halfDF"]
		if (model == "cauchy") reciprocalRates = TRUE
		if (model == "gamma") reciprocalRates = TRUE
		if (model == "logN") reciprocalRates = FALSE
		tab = tab[order(tab[,"startYear"]),]
		cor = cor((tab[,"endLon"]-tab[,"startLon"]),(tab[,"endLat"]-tab[,"startLat"]))
		my_var = solve(matrix(my_prec,nrow=2))		
		sigma1 = sqrt(my_var[1,1]); sigma2 = sqrt(my_var[2,2])
		sigmas = c(sigma1, sigma2)
		cor = my_var[1,2]/(sqrt(my_var[1,1])*sqrt(my_var[2,2]))
		localTreesDirectory = localTreesDirectory
		simulation = simulatorRRW1(tree, rates, sigmas, cor, envVariables, mostRecentSamplingDatum,
							       ancestPosition, reciprocalRates, n1, n2, showingPlots, newPlot)
		if (showingPlots == TRUE)
			{
				dev.new(); plot(backgroundRaster, col="gray90", box=F, ann=F, axes=F, legend=F)
				for (j in 1:dim(tab)[1])
					{
						segments(tab[j,"startLon"], tab[j,"startLat"], tab[j,"endLon"], tab[j,"endLat"], col="red", lwd=0.2)
					}
				dev.new(); plot(backgroundRaster, col="gray90", box=F, ann=F, axes=F, legend=F)
				for (j in 1:dim(simulation)[1])
					{
						segments(simulation[j,"startLon"], simulation[j,"startLat"], simulation[j,"endLon"], simulation[j,"endLat"], col="red", lwd=0.2)
					}
			}
		file = as.character(paste(localTreesDirectory,"/TreeSimulations_",i,".csv",sep=""))
		write.csv(simulation, file, row.names=F, quote=F)
	}}

i = 100; sim = read.csv(paste(localTreesDirectory,"/TreeSimulations_",i,".csv",sep=""), head=T)
dev.new(width=5, height=4.7); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(0,0,0,0), mgp=c(0,0,0), lwd=0.3)
convex_hull = raster(paste0("Environmental_rasters/Minimum_convex_hull_",analysis,"_04.asc"))
plot(convex_hull, main="", cex.main=0.8, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="gray90", colNA="white")
for (i in dim(sim)[1]:1)
	{
		curvedarrow(cbind(sim[i,"startLon"],sim[i,"startLat"]), cbind(sim[i,"endLon"],sim[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.15, dr=NA, endhead=F)
	}
for (i in dim(sim)[1]:1)
	{
		points(sim[i,"endLon"], sim[i,"endLat"], pch=16, col="red", cex=0.3)
		if (i == 1) points(sim[i,"startLon"], sim[i,"startLat"], pch=16, col="red", cex=0.3)
	}

# 13. Testing the impact of environmental factors on lineage dispersal locations

analysis = "Genomes"; localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900
analysis = "N_genes"; localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900

envVariables = list(); envVariableNames = c("Land_cover_forests_RABV_CA_04","Land_cover_savannas_RABV_CA_04","Land_cover_grasslands_RABV_CA_04","Land_cover_croplands_RABV_CA_04",
					   "Land_cover_water_RABV_CA_04","Human_pop_density_RABV_CA_04","Elevation_RABV_CA_04","Annual_mean_temperature_RABV_CA_04","Annual_precipitation_RABV_CA_04")for (i in 1:length(envVariableNames)) envVariables[[i]] = raster(paste0("Environmental_rasters/",envVariableNames[i],".asc"))
for (i in 1:nberOfExtractionFiles)
	{
		obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		sim = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",i,".csv"), header=T)
		envValues_obs = matrix(nrow=dim(obs)[1], ncol=length(envVariables))
		envValues_sim = matrix(nrow=dim(sim)[1], ncol=length(envVariables))
		colnames(envValues_obs) = envVariableNames; colnames(envValues_sim) = envVariableNames
		for (j in 1:length(envVariables))
			{
				if (dim(envVariables[[j]])[3] > 1)
					{
						time_intervals = matrix(nrow=length(names(envVariables[[j]])), ncol=2)
						for (k in 1:length(names(envVariables[[j]])))
							{
								time_intervals[k,1] = as.numeric(unlist(strsplit(names(envVariables[[j]])[k],"_"))[2])
								time_intervals[k,2] = as.numeric(unlist(strsplit(names(envVariables[[j]])[k],"_"))[3])
							}
					}
				envValues_obs[,j] = raster::extract(envVariables[[j]], SpatialPoints(obs[,c("endLon","endLat")]))
				envValues_sim[,j] = raster::extract(envVariables[[j]], SpatialPoints(sim[,c("endLon","endLat")]))
			}
		write.csv(envValues_obs, paste0(localTreesDirectory,"/EnvValues_obs_",i,".csv"), row.names=F, quote=F)
		write.csv(envValues_sim, paste0(localTreesDirectory,"/EnvValues_sim_",i,".csv"), row.names=F, quote=F)
	}
BFs = matrix(nrow=length(envVariableNames), ncol=2)
row.names(BFs) = envVariableNames; colnames(BFs) = c("lower","higher")
meanEnvValues_obs_list = list(); meanEnvValues_sim_list = list()
for (i in 1:length(envVariableNames))
	{
		lowerEnvValues = 0; meanEnvValues_obs_list = list(); meanEnvValues_sim_list = list()
		meanEnvValues_obs = rep(NA, nberOfExtractionFiles); meanEnvValues_sim = rep(NA, nberOfExtractionFiles)
		for (j in 1:nberOfExtractionFiles)
			{
				meanEnvValues_obs[j] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_obs_",j,".csv"))[,gsub("-",".",envVariableNames[i])], na.rm=T)
				meanEnvValues_sim[j] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_sim_",j,".csv"))[,gsub("-",".",envVariableNames[i])], na.rm=T)
				if (meanEnvValues_obs[j] < meanEnvValues_sim[j]) lowerEnvValues = lowerEnvValues+1				
			}
		p = lowerEnvValues/nberOfExtractionFiles; BFs[i,"lower"] = round((p/(1-p))/(0.5/(1-0.5)),1)
		p = (1-(lowerEnvValues/nberOfExtractionFiles)); BFs[i,"higher"] = round((p/(1-p))/(0.5/(1-0.5)),1)
		meanEnvValues_obs_list[[i]] = meanEnvValues_obs; meanEnvValues_sim_list[[i]] = meanEnvValues_sim
	}
write.csv(BFs, paste0("All_seraphim_results/",analysis,"_support_dispersal_location.csv"), quote=F)

# 14. Testing the impact of environmental factors on lineage dispersal velocity

nberOfExtractionFiles = 900; source("spreadFactors_mod.r")
analysis = "Genomes"; localTreesDirectory = "Tree_extraction_files/Genomes"
showingPlots = FALSE; nberOfCores = 1; OS = "Unix"; randomisations = FALSE
fourCells = FALSE; nberOfRandomisations = 1; randomProcedure = 2
c = 0; envVariables = list(); resistances = list(); avgResistances = list()
envVariableNames = c("Land_cover_forests_RABV_CA_04")
for (k in c(10))
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
	}
pathModel = 2; outputName = paste0("TEST_randomProcedure2")
spreadFactors_mod(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
				  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)

analysis = "Genomes"; localTreesDirectory = "Tree_extraction_files/Genomes"; nberOfExtractionFiles = 900
analysis = "N_genes"; localTreesDirectory = "Tree_extraction_files/N_genes"; nberOfExtractionFiles = 900

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
pathModel = 2; outputName = paste0("All_seraphim_results/",analysis,"_LC_extractions")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)
pathModel = 3; outputName = paste0("All_seraphim_results/",analysis,"_CS_extractions")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F,randomisations=F)	
pathModel = 2; outputName = paste0("All_seraphim_results/",analysis,"_LC_simulations")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	
pathModel = 3; outputName = paste0("All_seraphim_results/",analysis,"_CS_simulations")
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	

extractions = list(); simulations = list()
extractions[[1]] = read.table(paste0("All_seraphim_results/",analysis,"_CS_extractions_LR_results.txt"), header=T)
simulations[[1]] = read.table(paste0("All_seraphim_results/",analysis,"_CS_simulations_LR_results.txt"), header=T)
allResults = matrix(nrow=length(envVariableNames)*1*2*3, ncol=6); kS = c(10,100,1000); CR = c("C","R"); L = 0
colnames(allResults) = c("Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
pathModels = c("Circuitscape path model")
for (i in 1:length(pathModels))
	{
		for (j in 1:length(envVariableNames))
			{
				for (k in 1:length(CR))
					{
						for (l in 1:length(kS))
							{
								L = L+1; c = 0; rasterName = gsub(".asc","",envVariableNames[j])
								allResults[L,2] = kS[l]; allResults[L,1] = paste0(rasterName," (",CR[k],")")
								index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(rasterName,colnames(extractions[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(rasterName,colnames(extractions[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(rasterName,colnames(simulations[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
								R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]
								allResults[L,3] = paste0(round(median(R2),3)," [",round(quantile(R2,0.025),3)," - ",round(quantile(R2,0.975),3),"]")
								allResults[L,4] = paste0(round(median(Qe),3)," [",round(quantile(Qe,0.025),3)," - ",round(quantile(Qe,0.975),3),"]")
								allResults[L,5] = round(sum(Qe>0)/nberOfExtractionFiles, 2)
								if ((sum(Qe>0)/length(Qe)) > 0.9)
									{
										for (m in 1:length(Qe))
											{
												if (Qs[m] < Qe[m]) c = c+1
											}
										p = c/length(Qe); BF = p/(1-p)
										allResults[L,6] = round(BF,1)
									}	else	{
										allResults[L,6] = "-"
									}
							}
					}
			}
	}
write.csv(allResults, paste0("All_seraphim_results/",analysis,"_support_dispersal_velocity.csv"), row.names=F, quote=F)

