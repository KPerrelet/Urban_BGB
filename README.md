# Urban_BGB
Code for the article titled "Aquatic and terrestrial environmental DNA signals reveal decoupling of blue-green communities along an urbanization gradient"


Data is available on Zenodo (10.5281/zenodo.15210074) 


Raw sequences for the data generated are available on the European Nucleotide Archive under accession number PRJEB88517.



--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 

Urban_BGB_OTU.csv - OTU table. 

Rows = sites, columns = genus. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction).

 

Urban_BGB_taxo.csv - Taxonomy of all retained OTUs.

 

Urban_BGB_covariates.csv - All covariate values extracted for all sites. As many sites are private, their location cannot be disclosed publically.

 

Urban_BGB_life_history.csv - Life history (i.e., terrestrial, aquatic, or semi-terrestrial/semi-aquatic life cycle) for all shared OTUs.

 

metabarcoding_metadata.csv - metadata to process the raw sequences. 

 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Data dictionnary: 


Site_No: site ID

type: type of sample (water/soil)

Uniq_cd: alternative side ID

Smplng_d: sampling date

Connectivity_terrestrial: patch cohesion index, calculated using landscapemetrics

Connectivity_aquatic: mean distance between ponds in a 500 m buffer, calculated using landscapemetrics

Pond_Nb_500: number of ponds in a 500 m buffer (# ponds)

area_pond: size of the pond (m2)

area_site: size of the urban green space (m2)

Distance_Water: distance to the closest water body (m)

Distance_Forest: distance to the forest (m)

Distance_Pond: distance to the closest pond (m)

overwarming: mean deviation of nighttime temperature compare to the regional average (K)

no2: annual mean air NO2 concentration

population: population density in a 500 m buffer (# inhabitants)

employment: employment density in a 500 m buffer (# employees)

fish_PCR: fish presence, detected via PCR (T/F)

Walls: presence of walls at the pond banks (T/F)

tot_volume_filtered: total volume of water filtered across three water samples from the same site

Blue50_frac: fraction of blue surface in a 50 m buffer

Green50_frac: fraction of green surface in a 50 m buffer

Grey50_frac: fraction of grey surface in a 50 m buffer

LowManagement50_frac: fraction of low-maintained green surface in a 50 m buffer, based on habitat classification from Casanelles-Abella, J., Chauvier, Y., Zellweger, F., Villiger, P., Frey, D., Ginzler, C., ... & Pellissier, L. (2021). Applying predictive models to study the ecological properties of urban ecosystems: A case study in Zürich, Switzerland. Landscape and Urban Planning, 214, 104137.

LowRandom50_frac: fraction of low-maintained green surface in a 50 m buffer, based on random habitat classification

Blue500_frac: fraction of blue surface in a 500 m buffer

Green500_frac: fraction of green surface in a 500 m buffer

Grey500_frac: fraction of grey surface in a 500 m buffer

LowManagement500_frac: fraction of low-maintained green surface in a 500 m buffer, based on habitat classification from Casanelles-Abella, J., Chauvier, Y., Zellweger, F., Villiger, P., Frey, D., Ginzler, C., ... & Pellissier, L. (2021). Applying predictive models to study the ecological properties of urban ecosystems: A case study in Zürich, Switzerland. Landscape and Urban Planning, 214, 104137.

LowRandom500_frac: fraction of low-maintained green surface in a 500 m buffer, based on random habitat classification
