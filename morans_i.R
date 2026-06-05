#Script Summary:
# 1. Run the (global) moran's i test
# 2. Run a localized moran's i test (LISA)


#NOTE: Moran's i is not appropriate test for the county level. This test calculates spatial autocorrelation---meaning
#it looks at a regions geographic neighbors to detect whether the spread of a variable is clustered or dispersed. There are
#only two sets of neighbors at the county level--Manhattan/the Bronx and Queens/Brooklyn. Staten Island would be disregareded
#and there isn't much statistical power in drawing conclusions about clustering for only 2 regions with 2 neighbors in each of
#them. For these reasons, we focus solely on the neighborhood level.

#To read more about (global) moran's i: https://www.paulamoraga.com/book-spatial/spatial-autocorrelation.html

#To read more about localized moran's i: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-how-cluster-and-outlier-analysis-anselin-local-m.html

#To read more about the spdep library: https://cran.r-project.org/web/packages/spdep/spdep.pdf


#import necessary libraries
library(sf)
library(spdep)

#main data with all necessary variables
data_neigh <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")

#mss only --- used for genomic testing
data_neigh <- read.csv("~\\clinical_and_genomic\\genomic_nyc_mss_colorec_20.csv")

#for stage progression and time_msk_90
data_neigh <- read.csv("~\\stages\\stages_nyc_nonmss_colorec_20_real.csv")

#this function appends the geometry data (spatial data) to the data set of variables. It also cleans up the spatial data 
#and maps it to the neighborhood level. Swap out the word "Neighborhood" for "County" or your zip code variable name to map the
#data to those region levels

#this is the same function that is in the geospaital_visualizations script
neigh_level <- function(df){
  lookup <- unique(df[, c("Neighborhood", "PT_ZIP_CD")])
  lookup$PT_ZIP_CD <- as.character(lookup$PT_ZIP_CD)
  
  geojson <- read_sf("~\\for_geospatial_vis\\NYC_zipcodes.geojson")
  
  names(geojson)[3] <- "PT_ZIP_CD"
  
  new_data = data.frame()
  
  for (i in 1:nrow(geojson)){
    zipcodes <- strsplit(geojson$PT_ZIP_CD[i], ",")[[1]]
    for (zipcode in zipcodes){
      new_row <- geojson[i, ]
      new_row$PT_ZIP_CD <- zipcode
      new_data <- rbind(new_data, new_row)
    }
  }
  
  new_data$PT_ZIP_CD <- trimws(new_data$PT_ZIP_CD)
  
  new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  new_data <- new_data %>%
    left_join(lookup, by = "PT_ZIP_CD")
  
  
  subset <- c("Neighborhood", "geometry")
  #subset <- c("County", "geometry")
  new_data <- new_data[subset]

  #names(main_data)[12] <- "Yost_Index"
  names(data_neigh)[13] <<- "Yost_Index"
  new_data <- na.omit(new_data)
  new_data <<- new_data %>%
    group_by(Neighborhood) %>%
    summarize(geometry = st_union(geometry), .groups = "drop")
  #st_cast(new_data, "MULTIPLOYGON")
  
}
neigh_level(data_neigh)

#merge the geometry data and the variables
new_data <- new_data %>%
  left_join(data_neigh, by = "Neighborhood")

#Returns strange errors when placed inside a function--just run manually
#this cleans up the spatial data--turns any "POLYGON" data to "MULTIPOLYGON" data. 
new_data$geometry <- st_sfc(lapply(new_data$geometry, function(g){
  if (st_geometry_type(g) == "POLYGON") st_cast(g, "MULTIPOLYGON") else g
}), crs = st_crs(new_data))







#perform moran's i on yost index
yost <- function(){
  #calculate median yost index for each neigh
  median_yost <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(median_yost = median(Yost_Index, na.rm = TRUE))
  
  #Moran's i takes one geometry per region---reduce data set to one entry per neighborhood
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  #merge median yost and subset of neighs
  new_data <- left_join(new_data, median_yost,by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  
  new_data$Neighborhood <- as.character(new_data$Neighborhood)
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  
  #PS: The excluded regions are upper west side, rockaway, and pelham - throgs neck
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally to the final test statistic. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)

  #ensure lw and new_data$Neighborhood are in the same order. This is crucial!
  lw_order <- attr(lw, "Neighborhood")
  stopifnot(all(new_data$Neighborhood == lw_order))
  
  #run moran's i!
  moran_yost <<- moran.test(new_data$median_yost, lw, zero.policy = TRUE)
  print(moran_yost)
  #run localized moran's i
  lisa_yost <<- localmoran_perm(new_data$median_yost, lw, nsim = 9999, zero.policy = TRUE)
  #*****to view the lisa results, click on the variable in the global environment pane!!!!!!!!
}
yost()

bmi <- function(){
  median_bmi <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(median_bmi = median(First_BMI_Measurement, na.rm = TRUE))
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, median_bmi, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  lw_order <- attr(lw, "Neighborhood")
  stopifnot(all(new_data$Neighborhood == lw_order))
  
  moran_bmi <<- moran.test(new_data$median_bmi, lw, zero.policy = TRUE)
  print(moran_bmi)
}
bmi()

age <- function(){
  median_age <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(median_age = median(age_at_diag, na.rm = TRUE))
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, median_age, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  lw_order <- attr(lw, "Neighborhood")
  stopifnot(all(new_data$Neighborhood == lw_order))
  
  moran_age <<- moran.test(new_data$median_age, lw, zero.policy = TRUE)
  print(moran_age)
  lisa_age <<- localmoran_perm(new_data$median_age, lw, nsim = 9999, zero.policy = TRUE)
}
age()

sex <- function(){
  perc_sex_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_sex = sum(Sex == "Female"),
      perc_sex = num_sex / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_sex_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_sex <<- moran.test(new_data$perc_sex, lw, zero.policy = TRUE)
  print(moran_sex)
}
sex()

right <- function(){
  perc_right_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_right = sum(Primary_Tumor_Location_Final == "Right"),
      perc_right = num_right / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_right_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_right <<- moran.test(new_data$perc_right, lw, zero.policy = TRUE)
  print(moran_right)
}
right()

left <- function(){
  perc_left_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_left = sum(Primary_Tumor_Location_Final == "Left"),
      perc_left = num_left / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_left_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_left <<- moran.test(new_data$perc_left, lw, zero.policy = TRUE)
  print(moran_left)
}
left()

rectum <- function(){
  perc_rectum_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_rectum = sum(Primary_Tumor_Location_Final == "Rectum"),
      perc_rectum = num_rectum / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_rectum_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_rectum <<- moran.test(new_data$perc_rectum, lw, zero.policy = TRUE)
  print(moran_rectum)
}
rectum()

mss <- function(){
  perc_mss_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_mss = sum(cohort == "MSS"),
      perc_mss = num_mss / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_mss_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_mss <- moran.test(new_data$perc_mss, lw, zero.policy = TRUE)
  print(moran_mss)
}
mss()

eur <- function(){
  perc_eur_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_eur = sum(Ancestry_Label == "EUR"),
      perc_eur = num_eur / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_eur_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  

  moran_eur <<- moran.test(new_data$perc_eur, lw, zero.policy = TRUE)
  print(moran_eur)
}
eur()

afr <- function(){
  perc_afr_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_afr = sum(Ancestry_Label == "AFR"),
      perc_afr = num_afr / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_afr_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_afr <<- moran.test(new_data$perc_afr, lw, zero.policy = TRUE)
  print(moran_afr)
}
afr()

eas <- function(){
  perc_eas_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_eas = sum(Ancestry_Label == "EAS"),
      perc_eas = num_eas / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_eas_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_eas <<- moran.test(new_data$perc_eas, lw, zero.policy = TRUE)
  print(moran_eas)
}
eas()

stage_IV_diag <- function(){
  perc_stage_diag_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_stage_diag = sum(Stage_at_Diagnosis_Final == "Stage_IV"),
      perc_stage_diag = num_stage_diag / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_stage_diag_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_stage_IV <<- moran.test(new_data$perc_stage_diag, lw, zero.policy = TRUE)
  print(moran_stage_IV)
}
stage_IV_diag()

stage_III_diag <- function(){
  perc_stage_diag_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_stage_diag = sum(Stage_at_Diagnosis_Final == "Stage_III"),
      perc_stage_diag = num_stage_diag / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_stage_diag_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  morane_stage_III <<- moran.test(new_data$perc_stage_diag, lw, zero.policy = TRUE)
  print(moran_stage_III)
}
stage_III_diag()

stage_I_or_II_diag <- function(){
  perc_stage_diag_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_stage_diag = sum(Stage_at_Diagnosis_Final == "Stage_I_or_II"),
      perc_stage_diag = num_stage_diag / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_stage_diag_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_I_or_II <<- moran.test(new_data$perc_stage_diag, lw, zero.policy = TRUE)
  print(moran_I_or_II)
}
stage_I_or_II_diag()

frac <- function(){
  median_frac <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(median_frac = median(Fraction_Genome_Altered_Portal, na.rm = TRUE))
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, median_frac, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_frac <<- moran.test(new_data$median_frac, lw, zero.policy = TRUE)
  print(moran_frac)
}
frac()

kras <- function(){
  data_neigh$KRAS <- as.character(data_neigh$KRAS)
  perc_kras_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_kras = sum(KRAS == "1"),
      perc_kras = num_kras / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_kras_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  lw_order <- attr(lw, "Neighborhood")
  stopifnot(all(new_data$Neighborhood == lw_order))
  
  moran.test(new_data$perc_kras, lw, zero.policy = TRUE)
  print(moran_kras)
  
  lisa_kras <<- localmoran_perm(new_data$perc_kras, lw, nsim = 9999, zero.policy = TRUE)
}
kras()

pik <- function(){
  data_neigh$PIK3CA <- as.character(data_neigh$PIK3CA)
  perc_pik_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_pik = sum(PIK3CA == "1"),
      perc_pik = num_pik / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_pik_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_pik <<- moran.test(new_data$perc_pik, lw, zero.policy = TRUE)
  print(moran_pik)
}
pik()

smad <- function(){
  data_neigh$SMAD4 <- as.character(data_neigh$SMAD4)
  perc_smad_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_smad = sum(SMAD4 == "1"),
      perc_smad = num_smad / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_smad_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_smad <<- moran.test(new_data$perc_smad, lw, zero.policy = TRUE)
  print(moran_smad)
}
smad()

stage_prog <- function(){

  perc_stage_prog_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_stage_prog = sum(Stage_Progression == "Already Metastatic"),
      perc_stage_prog = num_stage_prog / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_stage_prog_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_stage_prog <<- moran.test(new_data$perc_stage_prog, lw, zero.policy = TRUE)
  print(moran_stage_prog)
}
stage_prog()

time_diag <- function(){
  
  perc_time_diag_grouping <- data_neigh %>%
    group_by(Neighborhood) %>%
    summarize(
      total = n(),
      num_time_diag = sum(Time_Diag_to_MSK_90 == "Before 90 days"),
      perc_time_diag = num_time_diag / total)
  
  subset = c("Neighborhood", "geometry")
  new_data <- new_data[subset]
  new_data <- unique(new_data)
  
  new_data <- left_join(new_data, perc_time_diag_grouping, by = "Neighborhood")
  new_data <- st_as_sf(new_data)
  
  #this will print out 2 warning messages. It's fine. The warnings mean that our regions don't all touch--so it made n subgraphs
  #of the regions that do touch. This is because of our >=20 frequency count--some regions are excluded, disallowing for a 
  #continuous graph. Proceed as normal!
  neighbor_list <- poly2nb(new_data, row.names = new_data$Neighborhood, snap = 1e-6)
  #this basically makes an adjacency matrix of neighbors
  #style = "W" or row-standardized spatial weights makes each region contribute equally. 
  lw <- nb2listw(neighbor_list, style = "W", zero.policy = TRUE)
  
  
  
  moran_time <<- moran.test(new_data$perc_time_diag, lw, zero.policy = TRUE)
  print(moran_time)
}
time_diag()
