# Uncovering the endemic circulation of rabies in Cambodia

Full text is available [here]().

## Abstract

In epidemiology, endemicity characterises sustained pathogen circulation in a geographical area, which involves a circulation that is not being maintained by external introductions. Because it could potentially shape the design of public health interventions, there is an interest in fully uncovering the endemic pattern of a disease. Here, we use a phylogeographic approach to investigate the endemic signature of rabies virus (RABV) circulation in Cambodia. Cambodia is located in one of the most affected regions in the world, but RABV circulation between and within Southeast Asian countries remains understudied in the area. Our analyses are based on a new comprehensive data set of 199 RABV genomes collected between 2014 and 2017 as well as previously published Southeast Asian RABV sequences. We show that most Cambodian sequences belong to a distinct clade that has been circulating almost exclusively in Cambodia. Our results thus point toward rabies circulation in Cambodia that does not rely on external introductions, which could have concrete implications in terms of mitigation strategy in the region. More globally, our study illustrates how phylogeographic investigations can be performed to assess and characterise viral endemicity in a context of relatively limited data.

## Analyses

This repository contains input data and codes necessary to reproduce the results presented in the original paper. The R code `R_scripts_RABV_CA.r` notably contains the procedure for sample selection, the post-processing steps of the continuous phylogeographic analyses, and all the steps of the landscape phylogeographic analysis.    

### Sample selection

We could sequence only 208 RABV positive dog brain samples out of the 627 that were available at the Institut Pasteur du Cambodge at the time of the study. In order to maximise the spatiotemporal coverage of the final dataset for the phylogeographic analyses, we subsampled the 208 sequences using a Markov process that maximises the sum D of great-circle distances along the edges of a Delaunay triangulation network connecting all selected samples. The available and selected samples are listed in the folder `Selection_of_samples`.

### Maximum likelihood tree of the Cambodian sequences and a representative sample of the canine RABV genetic diversity 

We identified the canine RABV clades that circulate in Cambodia by reconstructing a maximum likelihood tree using all the Cambodian N genes sequenced in this study and available on Genbank, as well as a representative sample of the genetic diversity of canine RABV across continents. The accession numbers of all sequences are listed in the `All_accession_IDs.csv` table. A dictionary of the variables of this table is available in `All_accession_IDs_dictionary.csv`. Finally, the nexus file `Global_RABV_ML.tree` contains the maximum likelihood tree reconstructed with IQ-TREE v2.0.6.     

### Discrete phylogeography

We then identified the probable movements between Cambodia and the other Southeast Asian countries using a discrete phylogeographic model implemented in BEAST v1.10 and BEAGLE v3. 

The input data and the parameterization of the discrete phylogeographic model are available in the following BEAST XML file: `BEAST_DTA_SEA.xml`.   

The marginal posterior distributions are available in `BEAST_DTA_SEA.log` and the maximum clade credibility tree is available in `BEAST_DTA_SEA.tree`.

### Continuous phylogeography

We further investigated the circulation of RABV among dog populations within Cambodia in a continuous phylogeographic analysis. 

The input data and model parameterization of the continous phylogeographic models for the N genes and the genomes are available in the `BEAST_RRW_analyses` folder. In the R code `R_scripts_RABV_CA.r`, one can find the analysis procedure to estimate the dispersal statistics made available in the `All_dispersal_statistics` folder.

### Landscape phylogeography

Finally, we estimated the association between environmental distances calculated for a selection of environmental rasters (see `Environmental_rasters`) and both lineage dispersal velocity and lineage dispersal location. Due to storage limitations, we do not provide here all the data extracted from the tree posterior distributions to perform the landscape analysis but a subset of them. Nevertheless, we provide the quantitative results of this analysis in the `All_seraphim_results` folder.
