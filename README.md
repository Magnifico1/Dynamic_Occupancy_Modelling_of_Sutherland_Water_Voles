# Dynamic_Occupancy_Modelling_of_Sutherland_Water_Voles

The data we used here is a subset of that described in Sutherland et al. (2014), a study of a population of water voles (*Arvicola amphibius*), a riparian specialist species, in northwest Scotland. Water vole colonies occupy discrete patches of lush riparian habitat embedded within a matrix of unsuitable habitat and exhibit high levels of patch turnover - i.e., they represent a rare example of a classically functioning metapopulation.

The patch network consists of $R=114$ patches, each surveyed between $J=2$
 and $J=4$
 times between July and August (breeding season), from 2009 to 2012 ($T=4$
). Water voles use latrines as territory marking so the detection-nondetection data are the number of visits in which at least one latrine was detected at a site - the data ($y$
) are counts (the number of detections across the $J$
 within-year sampling occasions). 
 
 Our objective was to use the water vole data to empirically test ‘area-isolation’ predictions. To test the ‘area’ prediction we modelled the extinction probability using a Logit Linear model with patch length as a covariate, and to test the ‘isolation’ prediction, we modelled the colonisation probability using a Logit Linear model with connectivity as a covariate. From the description of the metapopulation system, we developed a State-Space formulation of the Dynamic Occupancy model,
 fitted the model to the data using JAGS language while assessing convergence and fit, and conducted Posterior Predictive Checking by including a Discrepancy Metric.
