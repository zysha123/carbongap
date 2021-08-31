This repo contains table sheets from our recent study on carbon gap potential from vegetation

Some key imagees can be accessed at
Mean carbon gap flux for global vegetated area during 2001-2018
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/CarbongapMean
Mean NPP during 2001-2018
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/NPPMean
Land (vegetation) cover for the year 2018
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/landcoverMOD12Q1_2018
Relative NPP (relative contribution) for the year 2018
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/relativeNpp2018
Region segmentation for the year 2018
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/UniqueUnits2018
World land area
https://code.earthengine.google.com/?asset=users/zongyaosha/imgs2018/world_land_area
All other images are available if requested.

Common columns explanation
PCT_n%: corresponding vegetated area is first sliced into sub-areas using the n percentile of carbon gap flux (order by carbon gap flux from low to high) at an interval of 5%, i.e., 0-5%, 510%, …, 95-100%. All items related to PCT_n% is computed for each corresponding sub-area.
modis_landn: indicates the land (vegetation) cover type from IGBP classification schema of MODIS 12Q1 dataset, where n=1,2, ...17. Note that 11.(Permanent wetlands), 13 (Urban and built-up lands), 15 (Snow and ice), 16 (Barren), and 17 (Water bodies) were not included for the statistics.

Table sheets

Carbongap_area_distribution_wld.csv (for Fig. 4)
Accumulative total carbon gap against accumulative total vegetated area. Area unit: m2, carbon gap unit: gC
Accumu. Area_this_PCT: total area (m2) within this PCT
Accumu. Area%_this_PCT: total area % within this PCT
Accumu. Carbongap_this_PCT: total accumunated carbon gap within this PCT

carbongap_npp_population_relationship.csv
Carbongap_Max: the maximum carbon gap flux within this PCT
Area_this_PCT: total area (m2) within this PCT
Popu_Mean: mean polulation density within this PCT
Carbongap_Mean: mean carbon gap flux within this PCT
NPP_Mean: mean NPP within this PCT

code_wld_nature.py: PYTHON code for computation

continent_veg_type_npp_flux.csv
Vegetation NPP for each IGBP land (vegetation) cover type (gC/m2/year)  in each continent/region

continent_vegtype_carbongap_flux.csv
Carbon gap flux for each IGBP land (vegetation) cover type (gC/m2/year) in each continent/region

continent_vegtype_carbongap_total.csv
Total carbon gap (gC) for each continent/region in each year

continent_vegtype_npp_total.csv
Total NPP (gC) for each continent/region in each year

deciding_optimal_window_size.csv
Data used to decide optimal window size, which is ~20km. 

vegtype_area_in_continent.csv
Vegetated area (m2) based on IGBP land (vegetation) cover type in each continent/region

git_share_nm_samples_result.csv
Carbon gap and field data validation from Inner Mongolia, China. 
Detail information about headers is as follows:
ID: row id,	year: the year for the record collected,	lon: longitude,	lat: latitude,	u_id: unique id for the record,	grazing: grazing intensity (unit: sheep unit*10^4/km^2) for the administrative county where the sample locates,	biomass: above ground biomass inside the protection site,	o_biomass: biomass outside the site,	carbongap: carbon gap at the sample location,	NPP: npp at the sample location,	npp_ratio:carbon gap/npp*100(%),	biomass_ratio: (biomass-o_biomass)/o_biomass*100(%), 	DAY1: greenupdate,	DAY2: sample collection date,	DAYS: days between DAY1 and DAY2,	grazing_2: biomass consumption by grazing (1.8kg/sheep unit/per day, unit: g/m^2),	b_a_ratio: below-above ground biomass ratio,	biomass_t: total biomass inside the protection area,	o_biomass_t:total biomass outside the protection area,	o_biomass_t_c: grazed biomass+o_biomass_t, i.e., the total biomass outside the protection area considering grazing,	biomass_gap: Max((biomass_t-o_biomass_t_c)/o_biomass_t_c*100(%),0)
