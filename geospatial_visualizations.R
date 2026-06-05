#Script Summary:
# 1.) Prepare and append geospatial mapping data to data frame of variables of interest
# 2.) Create geospatial visualizations with shading respective to the variables of interest


# *******The graphs corresponding to tests that yielded a statistically significant results were saved in the
# Graphs -> Geospatial Visualizations.********

#The main library used to make these visualizations--and do any geospatial data analysis--is "sf". You can read
#about it here: https://r-spatial.github.io/sf/

#load all necessary libraries
library(ggplot2)
library(sf)
library(dplyr)
library(scales)

#load in the necessary data
#NOTE: this is the data without the 20 frequency count threshold. This is on purpose. We still need the names of the regions that
#do not meet the threshold to shade them yellow (yellow = was not included in analysis)
main_data <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_0.csv")
main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)

#this is the data without the nonmss patients. We use this data when working with the genomic variables (KRAS, APC, etc.)
genomic_data <- read.csv("~\\clinical_and_genomic\\genomic_nyc_mss_colorec_0.csv")
genomic_data$PT_ZIP_CD <- as.character(genomic_data$PT_ZIP_CD)






#Start at the county level:






#Assign each county a geospatial datum 
county_level <- function(df){
  lookup <- unique(df[, c("County", "PT_ZIP_CD")]) #  reference table---just a table of all counties and corresponding zips
  
  geojson <- read_sf("~\\for_geospatial_vis\\NYC_zipcodes.geojson") #read in geojson file
  
  names(geojson)[3] <- "PT_ZIP_CD" #update the name of the zip col in the geojson so we can merge on zips
  
  new_data = data.frame() #initialize a df
  
  #the raw geojson is a bit wacky--the zips are entered strangely. Here we clean that up
  for (i in 1:nrow(geojson)){
    zipcodes <- strsplit(geojson$PT_ZIP_CD[i], ",")[[1]]
    for (zipcode in zipcodes){
      new_row <- geojson[i, ]
      new_row$PT_ZIP_CD <- zipcode
      new_data <- rbind(new_data, new_row)
    }
  }
  
  #remove any trailing spaces from the zip col in the new geojson (new_data)
  new_data$PT_ZIP_CD <- trimws(new_data$PT_ZIP_CD)
  
  #make sure the zips in the reference df and new_data are the same dtype so we can left_join
  lookup$PT_ZIP_CD <- as.character(lookup$PT_ZIP_CD)
  new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  
  #left_join !!!!!
  new_data <- new_data %>%
    left_join(lookup, by = "PT_ZIP_CD")
  
  #we can disregard the zip code data + all of the other random columns (like population) from the original geojson
  subset <- c("County", "geometry")
  new_data <- new_data[subset]
  
  #make sure there are no na values (there shouldn't be!!!)
  new_data <- na.omit(new_data)
  
  #in order to visualize the data using the "geometry" column of the geojson, we have to change it's "data type"
  #Here, we give the "geometry" column the "geometry" data type (that isn't exactly what's going on here, but you
  #get the jist). We also make sure there are no repeating geometries by grouping by county
  new_data <<- new_data %>%
    group_by(County) %>%
    summarize(geometry = st_union(geometry), .groups = "drop")
}

#Pass main_data first since we're going to start with socioeconomic, demographic, and clinical variables
county_level(main_data)

#Baseline map of county level
baseline_county <- function(){
  #We're manually assinging the colors we want to each county. They will be shaded with these colors
  final <- new_data %>%
    mutate(fill_color = case_when(
      County == "Bronx" ~ "#FF0000",
      County == "Kings" ~ "#0000FF",
      County == "Manhattan" ~ "#00FF00",
      County == "Queens" ~ "#c51b8a",
      County == "Staten Island" ~ "#34C9BD"
    ))
  
  #make plot!
  ggplot(final) +
    geom_sf(aes(fill = fill_color), color = "black") +
    theme(legend.position = "none")
}
baseline_county()

#yost and county
yost_county <- function(){
  
  median_yost <- main_data %>%
    group_by(County) %>%
    summarize(median_yost = median(Yost_Index, na.rm = TRUE))
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total <- total %>%
    left_join(median_yost, by = "County")
  
  final <- na.omit(total)
  
  
  ggplot(final) +
    geom_sf(aes(fill = median_yost), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of the Median Yost Index per NYC County of Patients with\nColorectal Cancer from MSK-IMPACT",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
  
  
}
yost_county()

#bmi and county
bmi_county <- function(){
  #get median bmi for each county
  median_bmi <- main_data %>%
    group_by(County) %>%
    summarize(median_bmi = median(First_BMI_Measurement, na.rm = TRUE))
  
  #combine variables of interest and data with geojson info
  total <- new_data %>%
    left_join(main_data, by = "County")
  #combine median_bmi (for color shading) and geojson/variables df
  total <- total %>%
    left_join(median_bmi, by = "County")
  
  final <- na.omit(total)
  
  #plot it!
  ggplot(final) +
    geom_sf(aes(fill = median_bmi), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of the Median First BMI Measurement per NYC County of\nPatients with Colorectal Cancer from MSK-IMPACT",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
bmi_county()

#ancestry labels
adm_perc_county <- function(df){
  adm_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_adm = sum(Ancestry_Label == "ADM", na.rm = TRUE),
              adm_perc = num_adm/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_adm <- total %>%
    left_join(adm_perc, by = "County")
  
  final_adm <- na.omit(total_adm)
  
  
  freq <- final_adm %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_adm <- merge(final_adm, freq, by = "County", all.x = TRUE)
  
  final_adm$fill_color <- ifelse(final_adm$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_adm %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_adm %>% filter(fill_color == "high"),
            aes(fill = adm_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Admixture Patients with Colorectal Cancer from\nMSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
adm_perc_county()

afr_perc_county <- function(df){
  afr_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_afr = sum(Ancestry_Label == "AFR", na.rm = TRUE),
              afr_perc = num_afr/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_afr <- total %>%
    left_join(afr_perc, by = "County")
  
  final_afr <- na.omit(total_afr)
  
  
  freq <- final_afr %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_afr <- merge(final_afr, freq, by = "County", all.x = TRUE)
  
  final_afr$fill_color <- ifelse(final_afr$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_afr %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_afr %>% filter(fill_color == "high"),
            aes(fill = afr_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of African Patients with Colorectal Cancer from\nMSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
afr_perc_county()

asj_perc_county <- function(df){
  asj_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_asj = sum(Ancestry_Label == "ASJ", na.rm = TRUE),
              asj_perc = num_asj/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_asj <- total %>%
    left_join(asj_perc, by = "County")
  
  final_asj <- na.omit(total_asj)
  
  
  freq <- final_asj %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_asj <- merge(final_asj, freq, by = "County", all.x = TRUE)
  
  final_asj$fill_color <- ifelse(final_asj$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_asj %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_asj %>% filter(fill_color == "high"),
            aes(fill = asj_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Ashkenazi Jewish Patients with Colorectal\nCancer from MSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
asj_perc_county()

eas_perc_county <- function(df){
  eas_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_eas = sum(Ancestry_Label == "EAS", na.rm = TRUE),
              eas_perc = num_eas/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_eas <- total %>%
    left_join(eas_perc, by = "County")
  
  final_eas <- na.omit(total_eas)
  
  
  freq <- final_eas %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_eas <- merge(final_eas, freq, by = "County", all.x = TRUE)
  
  final_eas$fill_color <- ifelse(final_eas$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_eas %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_eas %>% filter(fill_color == "high"),
            aes(fill = eas_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of East Asian Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
eas_perc_county()

eur_perc_county <- function(df){
  eur_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_eur = sum(Ancestry_Label == "EUR", na.rm = TRUE),
              eur_perc = num_eur/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_eur <- total %>%
    left_join(eur_perc, by = "County")
  
  final_eur <- na.omit(total_eur)
  
  
  freq <- final_eur %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_eur <- merge(final_eur, freq, by = "County", all.x = TRUE)
  
  final_eur$fill_color <- ifelse(final_eur$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_eur %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_eur %>% filter(fill_color == "high"),
            aes(fill = eur_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of European Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
eur_perc_county()

nam_perc_county <- function(df){
  nam_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_nam = sum(Ancestry_Label == "NAM", na.rm = TRUE),
              nam_perc = num_nam/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_nam <- total %>%
    left_join(nam_perc, by = "County")
  
  final_nam <- na.omit(total_nam)
  
  
  freq <- final_nam %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_nam <- merge(final_nam, freq, by = "County", all.x = TRUE)
  
  final_nam$fill_color <- ifelse(final_nam$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_nam %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_nam %>% filter(fill_color == "high"),
            aes(fill = nam_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Native American Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
nam_perc_county()

sas_perc_county <- function(df){
  sas_perc <- main_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_sas = sum(Ancestry_Label == "SAS", na.rm = TRUE),
              sas_perc = num_sas/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "County")
  total_sas <- total %>%
    left_join(sas_perc, by = "County")
  
  final_sas <- na.omit(total_sas)
  
  
  freq <- final_sas %>%
    count(County, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_sas <- merge(final_sas, freq, by = "County", all.x = TRUE)
  
  final_sas$fill_color <- ifelse(final_sas$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_sas %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_sas %>% filter(fill_color == "high"),
            aes(fill = sas_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of South Asian Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC County",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
sas_perc_county()

#for handling genomic data!!!!!!
county_level(genomic_data)

kras_perc_county <- function(){
  kras_perc <- genomic_data %>%
    group_by(County) %>%
    summarize(total_pop = n(),
              num_kras = sum(KRAS == "1", na.rm = TRUE),
              kras_perc = num_kras/total_pop
    )
  
  total <- new_data %>%
    left_join(genomic_data, by = "County")
  total_kras <- total %>%
    left_join(kras_perc, by = "County")
  
  final_kras <- na.omit(total_kras)
  st_geometry(final_kras) <- st_geometry(final_kras)
  
  ggplot(final_kras, aes(fill = kras_perc)) +
    geom_sf(color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of the Patients with KRAS Mutation and Colorectal\nCancer from MSK-IMPACT, per NYC County.",
      x = "Longitude of Counties in NYC",
      y = "Latitude of Counties in NYC"
    )
}
kras_perc_county()






#Now the neighborhood level:







#Assign each neighborhood a geospatial datum 
neigh_level <- function(df){
  lookup <- unique(df[, c("Neighborhood", "PT_ZIP_CD")])

  geojson <- read_sf("C:\\Users\\skrastm\\OneDrive - Memorial Sloan Kettering Cancer Center\\McKenzie\\Main_Datasets\\for_geospatial_vis\\NYC_zipcodes.geojson")

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
  lookup$PT_ZIP_CD <- as.character(lookup$PT_ZIP_CD)
  
  new_data <- new_data %>%
    left_join(lookup, by = "PT_ZIP_CD")


  subset <- c("Neighborhood", "geometry")
  #subset <- c("County", "geometry")
  new_data <- new_data[subset]

  #names(main_data)[12] <- "Yost_Index"
  names(df)[11] <- "Yost_Index"
  new_data <- na.omit(new_data)
  new_data <<- new_data %>%
    group_by(Neighborhood) %>%
    summarize(geometry = st_union(geometry), .groups = "drop")
}
neigh_level(main_data)

#basline maps
baseline_neigh <- function(){
  
  geom <- new_data
  
  subset <- c("County", "Neighborhood")
  data_subset <- main_mss_data[subset]
  
  new_data <- new_data %>%
    left_join(data_subset, by = "Neighborhood")
  
  new_data <- unique(new_data)
  
  final <- new_data %>%
    group_by(County) %>%
    mutate(dummy_perc = row_number() / n()) %>%
    ungroup()
  
  count <- main_mss_data %>%
    count(Neighborhood)
  
  final <- final %>%
    st_drop_geometry() %>%
    left_join(count, by = c("Neighborhood"))
    
  final <- final %>%
    mutate(fill_color = case_when(
      n < 20 ~ "gray80",
      County == "Bronx" ~ gradient_n_pal(c("#D485CF", "#E512D0"))(dummy_perc),
      County == "Kings" ~ gradient_n_pal(c("#F57878", "#BF1D1D"))(dummy_perc),
      County == "Manhattan" ~ gradient_n_pal(c("#55DB53", "#007D00"))(dummy_perc),
      County == "Queens" ~ gradient_n_pal(c("#2F50F5", "#02179E"))(dummy_perc),
      County == "Staten Island" ~ gradient_n_pal(c("#A8E3D7", "#34C9BD"))(dummy_perc),
    ))
  
  subset <- c("geometry", "Neighborhood")
  data_subset <- new_data[subset]
  
  final <- final %>%
    left_join(data_subset, by = "Neighborhood")
  
  final_truly <- data.frame()
  
  final <- st_sf(final, geometry = final$geometry)
  
  #num_neighs <- length(unique(new_data$Neighborhood))
  #colors <- colorRampPalette(c("#fde0dd", "#fa9fb5", "#c51b8a", "#7a0177"))(num_neighs)
  ggplot(final) +
    geom_sf(aes(fill = fill_color), color = "black") +
    scale_fill_identity() +
    #scale_fill_manual(values = colors) +
    theme(legend.position = "none")
}
baseline_neigh()

#yost
yost_neigh <- function(){
  median_yost <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(median_yost = median(Yost_Index, na.rm = TRUE))
  
  
  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #median_yost$PT_ZIP_CD <- as.character(median_yost$PT_ZIP_CD)
  #main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total <- total %>%
    left_join(median_yost, by = "Neighborhood")
  
  final <- na.omit(total)
  
  
  freq <- final %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  
  final <- merge(final, freq, by = "Neighborhood", all.x = TRUE)
  
  final$fill_color <- ifelse(final$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final %>% filter(fill_color == "high"),
            aes(fill = median_yost), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of the Median Yost Index per NYC Neighborhood of\nPatients with Colorectal Cancer from MSK-IMPACT",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    )
  
}
yost_neigh()

bmi_neigh <- function(){
  median_first_bmi_measurement <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(median_first_bmi_measurement = median(First_BMI_Measurement, na.rm = TRUE))
  
  
  #new_data$Neighborhood <- as.character(new_data$Neighborhood)
  #median_age_at_diag$Neighborhood <- as.character(median_age_at_diag$Neighborhood)
  #main_data$Neighborhood <- as.character(main_data$Neighborhood)
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total <- total %>%
    left_join(median_first_bmi_measurement, by = "Neighborhood")
  
  final <- na.omit(total)
  
  
  freq <- final %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  
  final <- merge(final, freq, by = "Neighborhood", all.x = TRUE)
  
  final$fill_color <- ifelse(final$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final %>% filter(fill_color == "high"),
            aes(fill = median_first_bmi_measurement), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
      name = "first_bmi\n_measurement"
    ) +
    labs(
      title = "Distribution of the Median First BMI Measurement per NYC\nNeighborhood of Patients with Colorectal Cancer from MSK-IMPACT",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
    
}
bmi_neigh()

age_neigh <- function(){
  median_age_at_diag <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(median_age_at_diag = median(age_at_diag, na.rm = TRUE))


  new_data$Neighborhood <- as.character(new_data$Neighborhood)
  median_age_at_diag$Neighborhood <- as.character(median_age_at_diag$Neighborhood)
  main_data$Neighborhood <- as.character(main_data$Neighborhood)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total <- total %>%
    left_join(median_age_at_diag, by = "Neighborhood")

  final <- na.omit(total)


  freq <- final %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()


  final <- merge(final, freq, by = "Neighborhood", all.x = TRUE)

  final$fill_color <- ifelse(final$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final %>% filter(fill_color == "high"),
            aes(fill = median_age_at_diag), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
      name = "median_age\n_at_diag"
    ) +
    labs(
      title = "Distribution of Age at Diagnosis of Patients with Colorectal\nCancer from MSK-IMPACT, per NYC Neighborhood ",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
age_neigh()

#Sex - Female !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fem_perc <- function(){
  fem_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_fem = sum(Sex == "Female", na.rm = TRUE),
              fem_perc = num_fem/total_pop
    )


  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)
  #fem_perc$PT_ZIP_CD <- as.character(fem_perc$PT_ZIP_CD)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_fem <- total %>%
    left_join(fem_perc, by = "Neighborhood")

  final_fem <- na.omit(total_fem)


  freq <- final_fem %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()

  final_fem <- merge(final_fem, freq, by = "Neighborhood", all.x = TRUE)

  final_fem$fill_color <- ifelse(final_fem$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final_fem %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_fem %>% filter(fill_color == "high"),
            aes(fill = fem_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
fem_perc()

#Primary Site Left !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
left_perc <- function(df){
  left_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_left = sum(Primary_Tumor_Location_Final == "Left", na.rm = TRUE),
              left_perc = num_left/total_pop
    )


  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)
  #left_perc$PT_ZIP_CD <- as.character(left_perc$PT_ZIP_CD)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_left <- total %>%
    left_join(left_perc, by = "Neighborhood")

  final_left <- na.omit(total_left)


  freq <- final_left %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()

  final_left <- merge(final_left, freq, by = "Neighborhood", all.x = TRUE)

  final_left$fill_color <- ifelse(final_left$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final_left %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_left %>% filter(fill_color == "high"),
            aes(fill = left_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Distribution of the Patients with Left-sided Colorectal\nCancer from MSK-IMPACT, per NYC Neighborhood",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
left_perc()

#Primary Site Right !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
right_perc <- function(df){
  right_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_right = sum(Primary_Tumor_Location_Final == "Right", na.rm = TRUE),
              right_perc = num_right/total_pop
    )


  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)
  #right_perc$PT_ZIP_CD <- as.character(right_perc$PT_ZIP_CD)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_right <- total %>%
    left_join(right_perc, by = "Neighborhood")

  final_right <- na.omit(total_right)

  freq <- final_right %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()

  final_right <- merge(final_right, freq, by = "Neighborhood", all.x = TRUE)

  final_right$fill_color <- ifelse(final_right$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final_right %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_right %>% filter(fill_color == "high"),
            aes(fill = right_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Distribution of the Patients with Right-sided Colorectal\nCancer from MSK-IMPACT, per NYC Neighborhood",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
right_perc()

#Primary Site Rectum !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rectum_perc <- function(){
  rectum_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_rectum = sum(Primary_Tumor_Location_Final == "Rectum", na.rm = TRUE),
              rectum_perc = num_rectum/total_pop
    )


  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #main_data$PT_ZIP_CD <- as.character(main_data$PT_ZIP_CD)
  #rectum_perc$PT_ZIP_CD <- as.character(rectum_perc$PT_ZIP_CD)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_rectum <- total %>%
    left_join(rectum_perc, by = "Neighborhood")

  final_rectum <- na.omit(total_rectum)

  freq <- final_rectum %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()

  final_rectum <- merge(final_rectum, freq, by = "Neighborhood", all.x = TRUE)

  final_rectum$fill_color <- ifelse(final_rectum$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final_rectum %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_rectum %>% filter(fill_color == "high"),
            aes(fill = rectum_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Distribution of the Patients with Rectum-sided Colorectal\nCancer from MSK-IMPACT, per NYC Neighborhood",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    )
}
rectum_perc()

nonmss_perc <- function(df){
  nonmss_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_nonmss = sum(cohort == "Non-MSS", na.rm = TRUE),
              nonmss_perc = num_nonmss/total_pop
    )


  #new_data$PT_ZIP_CD <- as.character(new_data$PT_ZIP_CD)
  #main_mss_data$PT_ZIP_CD <- as.character(main_mss_data$PT_ZIP_CD)

  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")

  #total$PT_ZIP_CD <- as.character(total$PT_ZIP_CD)
  #nonmss_perc$PT_ZIP_CD <- as.character(nonmss_perc$PT_ZIP_CD)

  total_nonmss <- total %>%
    left_join(nonmss_perc, by = "Neighborhood")

  final_nonmss <- na.omit(total_nonmss)

  freq <- final_nonmss %>%
    count(Neighborhood, name = "frequency")

  freq <- freq %>% st_drop_geometry()

  final_nonmss <- merge(final_nonmss, freq, by = "Neighborhood", all.x = TRUE)
  
  final_nonmss$fill_color <- ifelse(final_nonmss$frequency < 20, "low", "high")


  ggplot()+
    geom_sf(data = final_nonmss %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_nonmss %>% filter(fill_color == "high"),
            aes(fill = nonmss_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Microstallite Stability per NYC Neighborhood of\nPatients with Colorectal Cancer from MSK-IMPACT",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    )
}
nonmss_perc()

prior_perc <- function(){
  prior_meds_to_msk_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_prior = sum(Prior.Treatment.to.MSK..NLP. == "Prior medications to MSK", na.rm = TRUE),
              prior_meds_to_msk_perc = num_prior/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_prior <- total %>%
    left_join(prior_meds_to_msk_perc, by = "Neighborhood")
  
  final_prior <- na.omit(total_prior)
  
  freq <- final_prior %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_prior <- merge(final_prior, freq, by = "Neighborhood", all.x = TRUE)
  
  final_prior$fill_color <- ifelse(final_prior$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_prior %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_prior %>% filter(fill_color == "high"),
            aes(fill = prior_meds_to_msk_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
prior_perc()

smoke_perc <- function(){
  smoker_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_prior = sum(Smoking.History..NLP. == "Former/Current Smoker", na.rm = TRUE),
              smoker_perc = num_prior/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_smoke <- total %>%
    left_join(smoker_perc, by = "Neighborhood")
  
  final_smoke <- na.omit(total_smoke)
  
  freq <- final_smoke %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_smoke <- merge(final_smoke, freq, by = "Neighborhood", all.x = TRUE)
  
  final_smoke$fill_color <- ifelse(final_smoke$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_smoke %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_smoke %>% filter(fill_color == "high"),
            aes(fill = smoker_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
smoke_perc()

apc_perc <- function(){
  apc_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_apc = sum(APC == "1", na.rm = TRUE),
              apc_perc = num_apc/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_apc <- total %>%
    left_join(apc_perc, by = "Neighborhood")
  
  final_apc <- na.omit(total_apc)
  
  freq <- final_apc %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_apc <- merge(final_apc, freq, by = "Neighborhood", all.x = TRUE)
  
  final_apc$fill_color <- ifelse(final_apc$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_apc %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_apc %>% filter(fill_color == "high"),
            aes(fill = apc_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
apc_perc()

tp53_perc <- function(){
  tp53_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_tp53 = sum(TP53 == "1", na.rm = TRUE),
              tp53_perc = num_tp53/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_tp53 <- total %>%
    left_join(tp53_perc, by = "Neighborhood")
  
  final_tp53 <- na.omit(total_tp53)
  
  freq <- final_tp53 %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_tp53 <- merge(final_tp53, freq, by = "Neighborhood", all.x = TRUE)
  
  final_tp53$fill_color <- ifelse(final_tp53$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_tp53 %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_tp53 %>% filter(fill_color == "high"),
            aes(fill = tp53_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
tp53_perc()

braf_perc <- function(){
  braf_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_braf = sum(BRAF == "1", na.rm = TRUE),
              braf_perc = num_braf/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_braf <- total %>%
    left_join(braf_perc, by = "Neighborhood")
  
  final_braf <- na.omit(total_braf)
  
  freq <- final_braf %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_braf <- merge(final_braf, freq, by = "Neighborhood", all.x = TRUE)
  
  final_braf$fill_color <- ifelse(final_braf$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_braf %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_braf %>% filter(fill_color == "high"),
            aes(fill = braf_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
braf_perc()

smad4_perc <- function(){
  smad4_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_smad4 = sum(SMAD4 == "1", na.rm = TRUE),
              smad4_perc = num_smad4/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_smad4 <- total %>%
    left_join(smad4_perc, by = "Neighborhood")
  
  final_smad4 <- na.omit(total_smad4)
  
  freq <- final_smad4 %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_smad4 <- merge(final_smad4, freq, by = "Neighborhood", all.x = TRUE)
  
  final_smad4$fill_color <- ifelse(final_smad4$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_smad4 %>% filter(fill_color == "low"),
            fill = "gray80") +
    geom_sf(data = final_smad4 %>% filter(fill_color == "high"),
            aes(fill = smad4_perc), color = "black") +
    scale_fill_gradient(
      low = "#f1defa", high = "#0f0376"
    )
}
smad4_perc()

pik3ca_perc <- function(){
  pik3ca_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_pik3ca = sum(PIK3CA == "1", na.rm = TRUE),
              pik3ca_perc = num_pik3ca/total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_pik3ca <- total %>%
    left_join(pik3ca_perc, by = "Neighborhood")
  
  final_pik3ca <- na.omit(total_pik3ca)
  
  freq <- final_pik3ca %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_pik3ca <- merge(final_pik3ca, freq, by = "Neighborhood", all.x = TRUE)
  
  final_pik3ca$fill_color <- ifelse(final_pik3ca$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_pik3ca %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_pik3ca %>% filter(fill_color == "high"),
            aes(fill = pik3ca_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of the Patients with PIK3CA Mutation and Colorectal\nCancer from MSK-IMPACT, per NYC Neighborhood.",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    )
}
pik3ca_perc()

adm_perc <- function(df){
  adm_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_adm = sum(Ancestry_Label == "ADM", na.rm = TRUE),
              adm_perc = num_adm/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_adm <- total %>%
    left_join(adm_perc, by = "Neighborhood")
  
  final_adm <- na.omit(total_adm)
  
  
  freq <- final_adm %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_adm <- merge(final_adm, freq, by = "Neighborhood", all.x = TRUE)
  
  final_adm$fill_color <- ifelse(final_adm$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_adm %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_adm %>% filter(fill_color == "high"),
            aes(fill = adm_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Admixture Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )
}
adm_perc()

afr_perc <- function(df){
  afr_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_afr = sum(Ancestry_Label == "AFR", na.rm = TRUE),
              afr_perc = num_afr/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_afr <- total %>%
    left_join(afr_perc, by = "Neighborhood")
  
  final_afr <- na.omit(total_afr)
  
  
  freq <- final_afr %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_afr <- merge(final_afr, freq, by = "Neighborhood", all.x = TRUE)
  
  final_afr$fill_color <- ifelse(final_afr$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_afr %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_afr %>% filter(fill_color == "high"),
            aes(fill = afr_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of African Patients with Colorectal Cancer\nfrom MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )
}
afr_perc()

asj_perc <- function(df){
  asj_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_asj = sum(Ancestry_Label == "ASJ", na.rm = TRUE),
              asj_perc = num_asj/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_asj <- total %>%
    left_join(asj_perc, by = "Neighborhood")
  
  final_asj <- na.omit(total_asj)
  
  
  freq <- final_asj %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_asj <- merge(final_asj, freq, by = "Neighborhood", all.x = TRUE)
  
  final_asj$fill_color <- ifelse(final_asj$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_asj %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_asj %>% filter(fill_color == "high"),
            aes(fill = asj_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Ashkenazi Jewish Patients with Colorectal\nCancer from MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )
}
asj_perc()

eas_perc <- function(df){
  eas_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_eas = sum(Ancestry_Label == "EAS", na.rm = TRUE),
              eas_perc = num_eas/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_eas <- total %>%
    left_join(eas_perc, by = "Neighborhood")
  
  final_eas <- na.omit(total_eas)
  
  
  freq <- final_eas %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_eas <- merge(final_eas, freq, by = "Neighborhood", all.x = TRUE)
  
  final_eas$fill_color <- ifelse(final_eas$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_eas %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_eas %>% filter(fill_color == "high"),
            aes(fill = eas_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of East Asian Patients with Colorectal\nCancer from MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )
}
eas_perc()

eur_perc <- function(df){
  eur_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_eur = sum(Ancestry_Label == "EUR", na.rm = TRUE),
              eur_perc = num_eur/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_eur <- total %>%
    left_join(eur_perc, by = "Neighborhood")
  
  final_eur <- na.omit(total_eur)
  
  
  freq <- final_eur %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_eur <- merge(final_eur, freq, by = "Neighborhood", all.x = TRUE)
  
  final_eur$fill_color <- ifelse(final_eur$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_eur %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_eur %>% filter(fill_color == "high"),
            aes(fill = eur_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of European Patients with Colorectal\nCancer from MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    ) 
}
eur_perc()

nam_perc <- function(df){
  nam_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_nam = sum(Ancestry_Label == "NAM", na.rm = TRUE),
              nam_perc = num_nam/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_nam <- total %>%
    left_join(nam_perc, by = "Neighborhood")
  
  final_nam <- na.omit(total_nam)
  
  
  freq <- final_nam %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_nam <- merge(final_nam, freq, by = "Neighborhood", all.x = TRUE)
  
  final_nam$fill_color <- ifelse(final_nam$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_nam %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_nam %>% filter(fill_color == "high"),
            aes(fill = nam_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of Native American Patients with Colorectal\nCancer from MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )  
}
nam_perc()

sas_perc <- function(df){
  sas_perc <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_sas = sum(Ancestry_Label == "SAS", na.rm = TRUE),
              sas_perc = num_sas/total_pop
    )
  
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_sas <- total %>%
    left_join(sas_perc, by = "Neighborhood")
  
  final_sas <- na.omit(total_sas)
  
  
  freq <- final_sas %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_sas <- merge(final_sas, freq, by = "Neighborhood", all.x = TRUE)
  
  final_sas$fill_color <- ifelse(final_sas$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_sas %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_sas %>% filter(fill_color == "high"),
            aes(fill = sas_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF")
    )+
    labs(
      title = "Distribution of South Asian Patients with Colorectal\nCancer from MSK-IMPACT per NYC Neighborhood",
      x = "Longitude of Neighborhood in NYC",
      y = "Latitude of Neighborhood in NYC"
    )  
}
sas_perc()

#stage at diag
stage_diag_I_II_perc <- function(df){
  stage_I_II_diag_perc  <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_stage_I_II_diag_perc  = sum(Stage_at_Diagnosis_Final == "Stage_I_or_II", na.rm = TRUE),
              stage_I_II_diag_perc = num_stage_I_II_diag_perc /total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_stage_I_II_diag <- total %>%
    left_join(stage_I_II_diag_perc, by = "Neighborhood")
  
  final_stage_I_II_diag <- na.omit(total_stage_I_II_diag)
  
  
  freq <- final_stage_I_II_diag %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_stage_I_II_diag <- merge(final_stage_I_II_diag, freq, by = "Neighborhood", all.x = TRUE)
  
  final_stage_I_II_diag$fill_color <- ifelse(final_stage_I_II_diag$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_stage_I_II_diag %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_stage_I_II_diag %>% filter(fill_color == "high"),
            aes(fill = stage_I_II_diag_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Percentage of Patients from MSK-IMPACT by New York City\nNeighborhood with Stage I or II Colorectal Cancer at Diagnosis ",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
stage_diag_I_II_perc()

stage_diag_III_perc <- function(df){
  stage_III_diag_perc  <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_stage_III_diag_perc  = sum(Stage_at_Diagnosis_Final == "Stage_III", na.rm = TRUE),
              stage_III_diag_perc = num_stage_III_diag_perc /total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_stage_III_diag <- total %>%
    left_join(stage_III_diag_perc, by = "Neighborhood")
  
  final_stage_III_diag <- na.omit(total_stage_III_diag)
  
  
  freq <- final_stage_III_diag %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_stage_III_diag <- merge(final_stage_III_diag, freq, by = "Neighborhood", all.x = TRUE)
  
  final_stage_III_diag$fill_color <- ifelse(final_stage_III_diag$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_stage_III_diag %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_stage_III_diag %>% filter(fill_color == "high"),
            aes(fill = stage_III_diag_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Percentage of Patients from MSK-IMPACT by New York City\nNeighborhood with Stage III Colorectal Cancer at Diagnosis",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
stage_diag_III_perc()

stage_diag_IV_perc <- function(df){
  stage_IV_diag_perc  <- main_data %>%
    group_by(Neighborhood) %>%
    summarize(total_pop = n(),
              num_stage_IV_diag_perc  = sum(Stage_at_Diagnosis_Final == "Stage_IV", na.rm = TRUE),
              stage_IV_diag_perc = num_stage_IV_diag_perc /total_pop
    )
  
  total <- new_data %>%
    left_join(main_data, by = "Neighborhood")
  total_stage_IV_diag <- total %>%
    left_join(stage_IV_diag_perc, by = "Neighborhood")
  
  final_stage_IV_diag <- na.omit(total_stage_IV_diag)
  
  
  freq <- final_stage_IV_diag %>%
    count(Neighborhood, name = "frequency")
  
  freq <- freq %>% st_drop_geometry()
  
  final_stage_IV_diag <- merge(final_stage_IV_diag, freq, by = "Neighborhood", all.x = TRUE)
  
  final_stage_IV_diag$fill_color <- ifelse(final_stage_IV_diag$frequency < 20, "low", "high")
  
  
  ggplot()+
    geom_sf(data = final_stage_IV_diag %>% filter(fill_color == "low"),
            fill = "#F5F39D") +
    geom_sf(data = final_stage_IV_diag %>% filter(fill_color == "high"),
            aes(fill = stage_IV_diag_perc), color = "black") +
    scale_fill_gradientn( 
      colors = c("#FF0000", "#FFFFFF", "#0000FF"),
    ) +
    labs(
      title = "Percentage of Patients from MSK-IMPACT by New York City\nNeighborhood with Stage IV Colorectal Cancer at Diagnosis",
      x = "Longitude of Neighborhoods in NYC",
      y = "Latitude of Neighborhoods in NYC"
    ) 
}
stage_diag_IV_perc()
