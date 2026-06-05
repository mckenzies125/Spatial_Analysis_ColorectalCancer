#File Summary
# 1. Perform a baseline analysis of the data
# 2. Perform a univariate analysis: variable vs neighborhood
# 3. Visualize these relationships

library(forcats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(coin)


#Step 1: Perform baseline analysis of the data

#To see these visualizations without having to run the code, go to Graphs-> Preliminary Vis



#this data includes non-mss patients and has no threshold cutoff
nonmss_0 <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_0.csv")

#this vis shows the freq count of nyc neighborhoods in our data set, as well as the 20 freq cut off 
freq_count <- function(df){
  
  data <- df %>%
    count(Neighborhood) %>%
    mutate(
      Neighborhood = fct_infreq(Neighborhood, n),
      threshold = ifelse(n < 20, "Freq < 20", "Freq >= 20"))
  
  ggplot(data, aes(x = Neighborhood, y = n, fill = threshold)) +
    geom_col() +
    scale_fill_manual(values = c("Freq < 20" = "red", "Freq >= 20" = "blue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    labs(
      x = "New York City Neighborhoods of Patients with Colorectal Cancer from MSK-IMPACT",
      y = "Frequency Count of New York City Neighborhoods",
      title = "Frequency Count of New York City Neighborhoods of Patients\nwith Colorectal Cancer from MSK-IMPACT"
    )
  
}
freq_count(nonmss_0)


#this data is patient ids, sample ids, and curated subtype label of ALL cancer types
all_cur_subs <- read.csv("~\\raw_data\\MSK_impact_w_cur_sub.csv")

#this vis shows the curated subtypes freq count in our data
cur_subs_vis <- function(df){
  
  data <- df %>%
    mutate(Curated_Subtypes = fct_infreq(Curated_Subtypes))
  
  ggplot(data) +
    aes(x = Curated_Subtypes)+
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    labs(
      x = "Curated Subtypes in MSK-IMPACT Identified by MSK-MET",
      y = "Frequency Count of Curated Subtypes in MSK-IMPACT",
      title = "Frequency Count of MSK-MET Curated Subtypes Identified\nin MSK-IMPACT"
    )
  
}
cur_subs_vis(all_cur_subs)



#Step 2: Univariate analysis


#to see the resuls without running the code, go to Calculate_p_q_vals -> Univariate


help <- data20 %>%
  group_by(County) %>%
  summarise(lol = median(First_BMI_Measurement, na.rm  = TRUE))


#includes non-mss patients and a 20 freq count cut off ---- USE THIS DATA FOR CLINICAL TESTING (NOT GENOMIC)
data20 <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")
names(data20)[11] <- "Smoking_History"
names(data20)[13] <- "Yost_Index"
names(data20)[21] <- "Prior_Treatment"

#includes ONLY mss patients and a 20 freq count cut off ---- USE THIS DATA FOR GENOMIC TESTING (KRAS, TP53, ETC.)
data20_mss_cohort_only <- read.csv("~\\clinical_and_genomic\\genomic_nyc_mss_colorec_20.csv")
names(data20_mss_cohort_only)[11] <- "Smoking_History"
names(data20_mss_cohort_only)[13] <- "Yost_Index"
names(data20_mss_cohort_only)[21] <- "Prior_Treatment"

#this is a special data set for stage_progression and time_diag_to_msk_90 ---- there are A LOT (~400) null entries for these two variables, 
#and so we use a special data set without these null entries 
data20_stage <- read.csv("~\\stages\\stages_nyc_nonmss_colorec_20.csv")
names(data20_stage)[12] <- "Smoking_History"
names(data20_stage)[14] <- "Yost_Index"
names(data20_stage)[22] <- "Prior_Treatment"

#initialize the p value county level matrix
pvals_county <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA", "Fraction_Genome_Altered_Portal", "Stage_Progression", "Time_Diag_to_MSK_90", "Stage_at_Diagnosis_Final"), 
                          County = '')
#initialize the p value neighborhood level matrix
pvals_neigh <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA", "Fraction_Genome_Altered_Portal", "Stage_Progression", "Time_Diag_to_MSK_90", "Stage_at_Diagnosis_Final"), 
                          Neighborhood = '')
#initialize the q value county level matrix
qvals_county <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA", "Fraction_Genome_Altered_Portal", "Stage_Progression", "Time_Diag_to_MSK_90", "Stage_at_Diagnosis_Final"), 
                           County = '')
#initialize the q value neighborhood level matrix
qvals_neigh <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA", "Fraction_Genome_Altered_Portal", "Stage_Progression", "Time_Diag_to_MSK_90", "Stage_at_Diagnosis_Final"), 
                          Neighborhood = '')

#function to test all counties against each variable
allpvals_county <- function(df_mss_only, df, df_stage){
  #set a seed for reproducibility
  set.seed(42)
  
  #these must all be factors to work in the coin library
  df$County <- as.factor(df$County)
  df_stage$County <- as.factor(df_stage$County)
  df_mss_only$County <- as.factor(df_mss_only$County)
  
  df$Sex <- as.factor(df$Sex)
  df$Prior_Treatment <- as.factor(df$Prior_Treatment)
  df$Smoking_History <- as.factor(df$Smoking_History)
  df$Primary_Tumor_Location_Final <- as.factor(df$Primary_Tumor_Location_Final)
  df$Ancestry_Label <- as.factor(df$Ancestry_Label)
  df$cohort <- as.factor(df$cohort)
  df$Stage_at_Diagnosis_Final <- as.factor(df$Stage_at_Diagnosis_Final)
  
  df_mss_only$KRAS <- as.factor(df_mss_only$KRAS)
  df_mss_only$APC <- as.factor(df_mss_only$APC)
  df_mss_only$TP53 <- as.factor(df_mss_only$TP53)
  df_mss_only$BRAF <- as.factor(df_mss_only$BRAF)
  df_mss_only$SMAD4 <- as.factor(df_mss_only$SMAD4)
  df_mss_only$PIK3CA <- as.factor(df_mss_only$PIK3CA)
  
  df_stage$Stage_Progression <- as.factor(df_stage$Stage_Progression)
  df_stage$Time_Diag_to_MSK_90 <- as.factor(df_stage$Time_Diag_to_MSK_90)
  
  
  #county vs yost
  yost <- kruskal_test(
    Yost_Index ~ County,
    data = df
  )
  pvals_county[1,2] <<- pvalue(yost)
  
  
  #county vs age
  age <- kruskal_test(
    age_at_diag ~ County,
    data = df
  )
  pvals_county[2,2] <<- pvalue(age)
  
  
  #county vs location
  loc <- fisher.test(
    Primary_Tumor_Location_Final ~ County, data20
  )
  pvals_county[3,2] <<- pvalue(loc)
  
  
  #county vs sex
  sex <- chisq_test(
    Sex ~ County,
    data = df
  )
  pvals_county[4,2] <<- pvalue(sex)
  
  
  #county vs cohort
  mss <- chisq_test(
    cohort ~ County,
    data = df
  )
  pvals_county[5,2] <<- pvalue(mss)
  
  
  #county vs acest
  ancest <- chisq_test(
    Ancestry_Label ~ County,
    data = df
  )
  pvals_county[6,2] <<- pvalue(ancest)
  
  
  #county vs bmi
  bmi <- kruskal_test(
    First_BMI_Measurement ~ County,
    data = df
  )
  pvals_county[7,2] <<- pvalue(bmi)
  
  
  #county vs prior treatment
  prior <- chisq_test(
    Prior_Treatment ~ County,
    data = df
  )
  pvals_county[8,2] <<- pvalue(prior)
  
  
  #county vs smoking
  smoke <- chisq_test(
    Smoking_History ~ County,
    data = df
  )
  pvals_county[9,2] <<- pvalue(smoke)
  
  
  #county vs kras
  kras <- chisq_test(
    KRAS ~ County,
    data = df_mss_only
  )
  pvals_county[10,2] <<- pvalue(kras)
  
  
  #county vs apc
  apc <- chisq_test(
    APC ~ County,
    data = df_mss_only
  )
  pvals_county[11,2] <<- pvalue(apc)
  
  
  #county vs tp53
  tp <- chisq_test(
    TP53 ~ County,
    data = df_mss_only
  )
  pvals_county[12,2] <<- pvalue(tp)
  
  
  #county vs braf
  braf <- chisq_test(
    BRAF ~ County,
    data = df_mss_only
  )
  pvals_county[13,2] <<- pvalue(braf)
  
  
  #county vs smad4
  smad <- chisq_test(
    SMAD4 ~ County,
    data = df_mss_only
  )
  pvals_county[14,2] <<- pvalue(smad)
  
  
  #county vs pik3ca
  pik <- chisq_test(
    PIK3CA ~ County,
    data = df_mss_only
  )
  pvals_county[15,2] <<- pvalue(pik)
  
  
  #county vs fraction
  frac <- kruskal_test(
    Fraction_Genome_Altered_Portal ~ County,
    data = df
  )
  pvals_county[16,2] <<- pvalue(frac)

  
  #county vs stage_progression
  stage_prog <- chisq_test(
    Stage_Progression ~ County,
    data = df_stage
  )
  pvals_county[17,2] <<- pvalue(stage_prog)
  
  
  #county vs time
  time <- chisq_test(
    Time_Diag_to_MSK_90 ~ County,
    data = df_stage
  )
  pvals_county[18,2] <<- pvalue(time)
  
  
  #county vs stage at diag
  stage_diag <- chisq_test(
    Stage_at_Diagnosis_Final ~ County,
    data = df
  )
  pvals_county[19,2] <<- pvalue(stage_diag)
  
}
allpvals_neigh <- function(df_mss_only, df, df_stage){
  set.seed(42)
  
  df$Neighborhood <- as.factor(df$Neighborhood)
  df_stage$Neighborhood <- as.factor(df_stage$Neighborhood)
  df_mss_only$Neighborhood <- as.factor(df_mss_only$Neighborhood)
  
  df$Sex <- as.factor(df$Sex)
  df$Prior_Treatment <- as.factor(df$Prior_Treatment)
  df$Smoking_History <- as.factor(df$Smoking_History)
  df$Primary_Tumor_Location_Final <- as.factor(df$Primary_Tumor_Location_Final)
  df$Ancestry_Label <- as.factor(df$Ancestry_Label)
  df$cohort <- as.factor(df$cohort)
  df$Stage_at_Diagnosis_Final <- as.factor(df$Stage_at_Diagnosis_Final)
  
  df_mss_only$KRAS <- as.factor(df_mss_only$KRAS)
  df_mss_only$APC <- as.factor(df_mss_only$APC)
  df_mss_only$TP53 <- as.factor(df_mss_only$TP53)
  df_mss_only$BRAF <- as.factor(df_mss_only$BRAF)
  df_mss_only$SMAD4 <- as.factor(df_mss_only$SMAD4)
  df_mss_only$PIK3CA <- as.factor(df_mss_only$PIK3CA)
  
  df_stage$Stage_Progression <- as.factor(df_stage$Stage_Progression)
  df_stage$Time_Diag_to_MSK_90 <- as.factor(df_stage$Time_Diag_to_MSK_90)
  
  
  
  #county vs yost
  yost <- kruskal_test(
    Yost_Index ~ Neighborhood,
    data = df
  )
  pvals_neigh[1,2] <<- pvalue(yost)
  
  
  #county vs age
  age <- kruskal_test(
    age_at_diag ~ Neighborhood,
    data = df
  )
  pvals_neigh[2,2] <<- pvalue(age)
  
  
  #county vs location
  loc <- chisq_test(
    Primary_Tumor_Location_Final ~ Neighborhood,
    data = df
  )
  pvals_neigh[3,2] <<- pvalue(loc)
  
  
  #county vs sex
  sex <- chisq_test(
    Sex ~ Neighborhood,
    data = df
  )
  pvals_neigh[4,2] <<- pvalue(sex)
  
  
  #county vs cohort
  mss <- chisq_test(
    cohort ~ Neighborhood,
    data = df
  )
  pvals_neigh[5,2] <<- pvalue(mss)
  
  
  #county vs acest
  ancest <- chisq_test(
    Ancestry_Label ~ Neighborhood,
    data = df
  )
  pvals_neigh[6,2] <<- pvalue(ancest)
  
  
  #county vs bmi
  bmi <- kruskal_test(
    First_BMI_Measurement ~ Neighborhood,
    data = df
  )
  pvals_neigh[7,2] <<- pvalue(bmi)
  
  
  #county vs prior treatment
  prior <- chisq_test(
    Prior_Treatment ~ Neighborhood,
    data = df
  )
  pvals_neigh[8,2] <<- pvalue(prior)
  
  
  #county vs smoking
  smoke <- chisq_test(
    Smoking_History ~ Neighborhood,
    data = df
  )
  pvals_neigh[9,2] <<- pvalue(smoke)
  
  
  #county vs kras
  kras <- chisq_test(
    KRAS ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[10,2] <<- pvalue(kras)
  
  
  #county vs apc
  apc <- chisq_test(
    APC ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[11,2] <<- pvalue(apc)
  
  
  #county vs tp53
  tp <- chisq_test(
    TP53 ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[12,2] <<- pvalue(tp)
  
  
  #county vs braf
  braf <- chisq_test(
    BRAF ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[13,2] <<- pvalue(braf)
  
  
  #county vs smad4
  smad <- chisq_test(
    SMAD4 ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[14,2] <<- pvalue(smad)
  
  
  #county vs pik3ca
  pik <- chisq_test(
    PIK3CA ~ Neighborhood,
    data = df_mss_only
  )
  pvals_neigh[15,2] <<- pvalue(pik)
  
  
  #county vs fraction
  frac <- kruskal_test(
    Fraction_Genome_Altered_Portal ~ Neighborhood,
    data = df
  )
  pvals_neigh[16,2] <<- pvalue(frac)
  
  
  #county vs stage_progression
  stage_prog <- chisq_test(
    Stage_Progression ~ Neighborhood,
    data = df_stage
  )
  pvals_neigh[17,2] <<- pvalue(stage_prog)
  
  
  #county vs time
  time <- chisq_test(
    Time_Diag_to_MSK_90 ~ Neighborhood,
    data = df_stage
  )
  pvals_neigh[18,2] <<- pvalue(time)
  
  
  #county vs stage at diag
  stage_diag <- chisq_test(
    Stage_at_Diagnosis_Final ~ Neighborhood,
    data = df
  )
  pvals_neigh[19,2] <<- pvalue(stage_diag)
  
}
allqvals_county <- function(df){
  
  #turn the matrix into a vector that maintains its order
  pval_vector <- df %>% 
    select(-1) %>% 
    pivot_longer(everything(), values_to = "val") %>%
    filter(!is.na(val), val != "") %>%
    pull(val)
  
  #ensure the values are numbers, not strings
  pval_vector <- as.numeric(pval_vector)
  
  #adjust!
  qval_vector <- p.adjust(pval_vector, method = "BH")

  #place the qvals back in their respective places
  qvals_county[1,2] <<- qval_vector[1]
  qvals_county[2,2] <<- qval_vector[2]
  qvals_county[3,2] <<- qval_vector[3]
  qvals_county[4,2] <<- qval_vector[4]
  qvals_county[5,2] <<- qval_vector[5]
  qvals_county[6,2] <<- qval_vector[6]
  qvals_county[7,2] <<- qval_vector[7]
  qvals_county[8,2] <<- qval_vector[8]
  qvals_county[9,2] <<- qval_vector[9]
  qvals_county[10,2] <<- qval_vector[10]
  qvals_county[11,2] <<- qval_vector[11]
  qvals_county[12,2] <<- qval_vector[12]
  qvals_county[13,2] <<- qval_vector[13]
  qvals_county[14,2] <<- qval_vector[14]
  qvals_county[15,2] <<- qval_vector[15]
  qvals_county[16,2] <<- qval_vector[16]
  qvals_county[17,2] <<- qval_vector[17]
  qvals_county[18,2] <<- qval_vector[18]
  qvals_county[19,2] <<- qval_vector[19]
}
allqvals_neigh <- function(df){
  
  #turn the matrix into a vector that maintains its order
  pval_vector <- df %>%
    select(-1) %>%
    pivot_longer(everything(), values_to = "val") %>%
    filter(!is.na(val), val != "") %>%
    pull(val)
  
  #ensure the values are numbers, not strings
  pval_vector <- as.numeric(pval_vector)

  #adjust!
  qval_vector <- p.adjust(pval_vector, method = "BH")

  #place the qvals back in their respective places
  qvals_neigh[1,2] <<- qval_vector[1]
  qvals_neigh[2,2] <<- qval_vector[2]
  qvals_neigh[3,2] <<- qval_vector[3]
  qvals_neigh[4,2] <<- qval_vector[4]
  qvals_neigh[5,2] <<- qval_vector[5]
  qvals_neigh[6,2] <<- qval_vector[6]
  qvals_neigh[7,2] <<- qval_vector[7]
  qvals_neigh[8,2] <<- qval_vector[8]
  qvals_neigh[9,2] <<- qval_vector[9]
  qvals_neigh[10,2] <<- qval_vector[10]
  qvals_neigh[11,2] <<- qval_vector[11]
  qvals_neigh[12,2] <<- qval_vector[12]
  qvals_neigh[13,2] <<- qval_vector[13]
  qvals_neigh[14,2] <<- qval_vector[14]
  qvals_neigh[15,2] <<- qval_vector[15]
  qvals_neigh[16,2] <<- qval_vector[16]
  qvals_neigh[17,2] <<- qval_vector[17]
  qvals_neigh[18,2] <<- qval_vector[18]
  qvals_neigh[19,2] <<- qval_vector[19]
}

#run the functions
allpvals_neigh(data20_mss_cohort_only, data20, data20_stage)
allpvals_county(data20_mss_cohort_only, data20, data20_stage)
allqvals_county(pvals_county)
allqvals_neigh(pvals_neigh)

#save to csv
save_values <- function(){
  write.csv(pvals_county, "pvals_county.csv", row.names = FALSE)
  write.csv(qvals_county, "qvals_county.csv", row.names = FALSE)
  write.csv(pvals_neigh, "pvals_neigh.csv", row.names = FALSE)
  write.csv(qvals_neigh, "qvals_neigh.csv", row.names = FALSE)
}



#Step 3: visualize the significant relationships between county/neighborhood and these variables



#County first:



#yost and county
yost_county_graph <- function(df){
  ggplot(df, aes(reorder(County, Yost_Index, FUN = median, decreasing = TRUE), Yost_Index, fill = County)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Bronx" = "#D936C6", "Kings" = "#961818", "Manhattan" = "#109126", "Queens" = "#252FD9", "Staten Island" = "#56B4E9")) +
    labs(title = "Association Between Yost Index and County of Origin of Patients with\nColorectal Cancer from the MSK-Impact Data Set",
         x = "New York City Counties of Origin of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "Yost Index of Patients with Colorectal Cancer from MSK-IMPACT") +
    theme_minimal() +
    annotate("text", x = 5, y = 95, label = "p-value < 1e-10") + #p val copied from qvals_county
    theme(legend.position = "none")
  
}
yost_county_graph(data20)



#bmi and county
bmi_county_graph <- function(df){
  ggplot(df, aes(reorder(County, First_BMI_Measurement, FUN = median, decreasing = TRUE), First_BMI_Measurement, fill = County)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Bronx" = "#D936C6", "Kings" = "#961818", "Manhattan" = "#109126", "Queens" = "#252FD9", "Staten Island" = "#56B4E9")) +
    labs(title = "Association Between First BMI Measurement and County of Origin of\nPatients with Colorectal Cancer from the MSK-Impact Data Set",
         x = "New York City Counties of Origin of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "First BMI Measurement of Patients with Colorectal Cancer from MSK-IMPACT") +
    theme_minimal() +
    annotate("text", x = 1, y = 50, label = "p-value = 	9.133136e-08") + #p val copied from qvals_county
    theme(legend.position = "none") +
    ylim(15, 50)
  
}
bmi_county_graph(data20)


#ancestry vs county
ancest_county_graph <- function(df){
  
  jitter <- df %>%
    count(County, Ancestry_Label, name = "count") %>%
    group_by(County) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()

  jitter$Ancestry_Label <- factor(jitter$Ancestry_Label, levels = rev(c("EUR", "ADM", "AFR", "ASJ", "EAS", "NAM", "SAS")))
  
  eur <- jitter %>%
    filter(Ancestry_Label == "EUR") %>%
    select(County, eur_cnt = prop)
  
  jitter <- jitter %>%
    left_join(eur, by = "County") %>%
    mutate(
      County = fct_reorder(County, eur_cnt, .desc = TRUE)
      )
  
  ggplot(jitter, aes(x = County, y = prop, fill = Ancestry_Label)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Ancestry Label Distribution of Patients with Colorectal Cancer From\nMSK-IMPACT by New York City County",
         x = "County of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Ancestry Label Distribution of Patients with Colorectal Cancer")
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #theme(legend.position = "none")
}
ancest_county_graph(data20)


kras_county_graph <- function(df){
  
  df$KRAS <- as.factor(df$KRAS)
  
  jitter <- data20 %>%
    count(County, KRAS, name = "count") %>%
    group_by(County) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$KRAS <- factor(jitter$KRAS, levels = rev(c("1", "0")))
  
  KRASS <- jitter %>%
    filter(KRAS == "1") %>%
    select(County, KRAS_cnt = prop)
  
  jitter <- jitter %>%
    left_join(KRASS, by = "County") %>%
    mutate(
      County = fct_reorder(County, KRAS_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = County, y = prop, fill = KRAS)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No KRAS Mutation", "KRAS Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with KRAS Mutation and Colorectal\nCancer From MSK-IMPACT by New York City County",
         x = "County of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with KRAS Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
kras_county_graph(data20_mss_cohort_only)


lollipop_kras <- function(df){
  help <- df %>%
    mutate(fill_color = case_when(
      County == "Bronx" ~ gradient_n_pal(c("#E512D0"))(0.2),
      County == "Kings" ~ gradient_n_pal(c("#BF1D1D"))(0.4),
      County == "Manhattan" ~ gradient_n_pal(c("#007D00"))(0.6),
      County == "Queens" ~ gradient_n_pal(c("#02179E"))(0.8),
      County == "Staten Island" ~ gradient_n_pal(c("#34C9BD"))(1),
    ))
  
  lolli <- df %>%
    count(Neighborhood, KRAS, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    filter(KRAS == "1") %>%
    arrange(desc(prop)) %>%
    mutate(Neighborhood = factor(Neighborhood, levels = Neighborhood))
  
  # lolli <- df %>%
  #   filter(KRAS == 1) %>%
  #   count(Neighborhood, name = "freq") %>%
  #   arrange(freq) %>%
  #   mutate(Neighborhood = factor(Neighborhood, levels = Neighborhood))
  # 
  lolli <- lolli %>%
    left_join(help, by = "Neighborhood")
  
  lolli <- lolli %>%
    mutate(County = factor(fill_color,
                           levels = c("#E512D0", "#BF1D1D", "#007D00", "#02179E", "#34C9BD"),
                           labels = c("Bronx", "Kings", "Manhattan", "Queens", "Staten Island")))
  
  lolli <- lolli %>%
    mutate(Neighborhood = fct_reorder(Neighborhood, prop, .desc = TRUE))
  
  ggplot(lolli, aes(x = prop, y = Neighborhood, color = County)) +
    geom_segment(aes(x = 0, xend = prop, y = Neighborhood, yend = Neighborhood),
                 color = "grey60", show.legend = FALSE) +
    geom_point(size = 4, show.legend = TRUE) +
    scale_fill_identity(
      name = "County",
      guide = "legend"
    ) +
    labs(
      title = "Percentage of KRAS Mutations in each Neighborhood of Origin of\nPatients with Colorectal Cancer from MSK-IMPACT Data Set",
      x = "Percentage of KRAS Mutations of Patients with Colorectal Cancer from MSK-IMPACT",
      y = "Neighborhoods of Origin of Patients with Colorectal Cancer from MSK-IMPACT"
    )
}
lollipop_kras(data20_mss_cohort_only)




#Neighborhood:




#yost and neighborhood
yost_neigh_graph <- function(df){

  df <- df %>%
    mutate(fill_color = case_when(
      County == "Bronx" ~ gradient_n_pal(c("#D485CF", "#E512D0"))(0.2),
      County == "Kings" ~ gradient_n_pal(c("#F57878", "#BF1D1D"))(0.4),
      County == "Manhattan" ~ gradient_n_pal(c("#55DB53", "#007D00"))(0.6),
      County == "Queens" ~ gradient_n_pal(c("#2F50F5", "#02179E"))(0.8),
      County == "Staten Island" ~ gradient_n_pal(c("#A8E3D7", "#34C9BD"))(1),
    ))
  
  df$County <- trimws(df$County)
  
  #df$Neighborhood <- reorder(df$Neighborhood, df$Yost_Index, .fun = median, .desc = TRUE)
  ggplot(df, aes(reorder(Neighborhood, Yost_Index, FUN = median, decreasing = TRUE), Yost_Index, fill = fill_color)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    #scale_fill_manual(values = c("Bronx" = "#E69F00", "Kings" = "#56B4E9", "Manhattan" = "#009E73", "Queens" = "#8b1ff0", "Staten Island" = "#e33232")) +
    scale_fill_identity() +
    labs(title = "Association Between Yost Index and Neighborhood of Origin of Patients with\nColorectal Cancer from the MSK-Impact Data Set",
         x = "New York City Neighborhoods of Origin of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "Yost Index of Patients with Colorectal Cancer from MSK-IMPACT") +
    theme_minimal() +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 78, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 77.5, as.character(Neighborhood),'')), hjust = 0.5,vjust = 1.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 79, as.character(Neighborhood),'')), hjust = 1,vjust = -0.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 95, as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 23, y = 95, label = "p-value < 1e-10") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
  
}
yost_neigh_graph(data20)


#age and neighborhood
age_neigh_graph <- function(df){
  
  df <- df %>%
    mutate(fill_color = case_when(
      County == "Bronx" ~ gradient_n_pal(c("#D485CF", "#E512D0"))(0.2),
      County == "Kings" ~ gradient_n_pal(c("#F57878", "#BF1D1D"))(0.4),
      County == "Manhattan" ~ gradient_n_pal(c("#55DB53", "#007D00"))(0.6),
      County == "Queens" ~ gradient_n_pal(c("#2F50F5", "#02179E"))(0.8),
      County == "Staten Island" ~ gradient_n_pal(c("#A8E3D7", "#34C9BD"))(1),
    ))
  
  df$County <- trimws(df$County)
  
  #df$Neighborhood <- reorder(df$Neighborhood, df$Yost_Index, .fun = median, .desc = TRUE)
  ggplot(df, aes(reorder(Neighborhood, age_at_diag, FUN = median, decreasing = TRUE), age_at_diag, fill = fill_color)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    #scale_fill_manual(values = c("Bronx" = "#E69F00", "Kings" = "#56B4E9", "Manhattan" = "#009E73", "Queens" = "#8b1ff0", "Staten Island" = "#e33232")) +
    scale_fill_identity() +
    labs(title = "Association Between Age at Diagnosis and Neighborhood of Origin of Patients with\nColorectal Cancer from the MSK-Impact Data Set",
         x = "New York City Neighborhoods of Origin of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "Age at Diagnosis of Patients with Colorectal Cancer from MSK-IMPACT") +
    theme_minimal() +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 78, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 77.5, as.character(Neighborhood),'')), hjust = 0.5,vjust = 1.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 79, as.character(Neighborhood),'')), hjust = 1,vjust = -0.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 95, as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 23, y = 95, label = "p-value < 	0.007171") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
  
}
age_neigh_graph(data20)


#bmi vs neigh
bmi_neigh_graph <- function(df){
  df <- df %>%
    mutate(fill_color = case_when(
      County == "Bronx" ~ gradient_n_pal(c("#D485CF", "#E512D0"))(0.2),
      County == "Kings" ~ gradient_n_pal(c("#F57878", "#BF1D1D"))(0.4),
      County == "Manhattan" ~ gradient_n_pal(c("#55DB53", "#007D00"))(0.6),
      County == "Queens" ~ gradient_n_pal(c("#2F50F5", "#02179E"))(0.8),
      County == "Staten Island" ~ gradient_n_pal(c("#A8E3D7", "#34C9BD"))(1),
    ))
  
  ggplot(df, aes(reorder(Neighborhood, First_BMI_Measurement, FUN = median, decreasing = TRUE), First_BMI_Measurement, fill = fill_color)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    #scale_fill_manual(values = c("Bronx" = "#E69F00", "Kings" = "#56B4E9", "Manhattan" = "#009E73", "Queens" = "#8b1ff0", "Staten Island" = "#e33232")) +
    labs(title = "Association Between First BMI Measurement and Neighborhood of Origin of\nPatients with Colorectal Cancer from the MSK-Impact Data Set",
         x = "New York City Neighborhoods of Origin of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "First BMI Measurement of Patients with Colorectal Cancer") +
    theme_minimal() +
    scale_fill_identity() +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 78, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 77.5, as.character(Neighborhood),'')), hjust = 0.5,vjust = 1.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 79, as.character(Neighborhood),'')), hjust = 1,vjust = -0.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 95, as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 20, y = 52, label = "p-value = 	2.366422e-05") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1)) +
    ylim(14, 53)
}
bmi_neigh_graph(data20)



#ancestry and neigh
ancest_neigh_graph <- function(df){
  
  jitter <- df %>%
    count(Neighborhood, Ancestry_Label, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$Ancestry_Label <- factor(jitter$Ancestry_Label, levels = rev(c("EUR", "ADM", "AFR", "ASJ", "EAS", "NAM", "SAS")))
  
  eur <- jitter %>%
    filter(Ancestry_Label == "EUR") %>%
    select(Neighborhood, eur_cnt = prop)
  
  jitter <- jitter %>%
    left_join(eur, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, eur_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = Ancestry_Label)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Ancestry Label Distribution of Patients with Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Ancestry Label Distribution of Patients with Colorectal Cancer") +
  #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
  #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
  #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
  #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
ancest_neigh_graph(data20)





#The following were not found to be significant, but were requested to still be made:


loc_neigh_graph <- function(df){
  
  jitter <- df %>%
    count(Neighborhood, Primary_Tumor_Location_Final, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$Primary_Tumor_Location_Final <- factor(jitter$Primary_Tumor_Location_Final, levels = rev(c("Right", "Left", "Rectum")))
  
  right <- jitter %>%
    filter(Primary_Tumor_Location_Final == "Right") %>%
    select(Neighborhood, right_cnt = prop)
  
   jitter <- jitter %>%
     left_join(right, by = "Neighborhood") %>%
     mutate(
       Neighborhood = fct_reorder(Neighborhood, right_cnt, .desc = TRUE)
     )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = Primary_Tumor_Location_Final)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Rectum" = "#E69F00", "Left" = "#56B4E9", "Right" = "#009E73")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Primary Tumor Location of Patients with Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Primary Tumor Location of Patients with Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
loc_neigh_graph(data20)



mss_neigh_graph <- function(df){
  
  jitter <- df %>%
    count(Neighborhood, cohort, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$cohort <- factor(jitter$cohort, levels = rev(c("MSS", "Non-MSS")))
  
  mss <- jitter %>%
    filter(cohort == "MSS") %>%
    select(Neighborhood, mss_cnt = prop)
  
  jitter <- jitter %>%
    left_join(mss, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, mss_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = cohort)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#E69F00", "Non-MSS" = "#56B4E9")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Microsatellite Stability of Patients with Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Microsatellite Stability of Patients with Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
mss_neigh_graph(data20)



apc_neigh_graph <- function(df){
  
  df$APC <- as.factor(df$APC)
  
  jitter <- df %>%
    count(Neighborhood, APC, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$APC <- factor(jitter$APC, levels = rev(c("1", "0")))
  
  APCC <- jitter %>%
    filter(APC == "1") %>%
    select(Neighborhood, APC_cnt = prop)
  
  jitter <- jitter %>%
    left_join(APCC, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, APC_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = APC)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No APC Mutation", "APC Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with APC Mutation and Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with APC Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
apc_neigh_graph(data20_mss_cohort_only)



tp_neigh_graph <- function(df){
  
  df$TP53 <- as.factor(df$TP53)
  
  jitter <- df %>%
    count(Neighborhood, TP53, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$TP53 <- factor(jitter$TP53, levels = rev(c("1", "0")))
  
  TP533 <- jitter %>%
    filter(TP53 == "1") %>%
    select(Neighborhood, TP53_cnt = prop)
  
  jitter <- jitter %>%
    left_join(TP533, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, TP53_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = TP53)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No TP53 Mutation", "TP53 Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with TP53 Mutation and Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with TP53 Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
tp_neigh_graph(data20_mss_cohort_only)



braf_neigh_graph <- function(df){
  
  jitter <- data20_nonmss %>%
    count(Neighborhood, BRAF, name = "count") %>%
    complete(Neighborhood, BRAF = c(0,1), fill = list(count = 0)) %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$BRAF <- as.factor(jitter$BRAF)
  
  jitter$BRAF <- factor(jitter$BRAF, levels = rev(c("1", "0")))
  
  BRAFs <- jitter %>%
    filter(BRAF == "1") %>%
    select(Neighborhood, BRAF_cnt = prop)
  
  #View(BRAFs)
  
  jitter <- jitter %>%
    left_join(BRAFs, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, BRAF_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = BRAF)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No BRAF Mutation", "BRAF Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with BRAF Mutation and Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with BRAF Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
braf_neigh_graph(data20_mss_cohort_only)



smad_neigh_graph <- function(df){
  
  df$SMAD4 <- as.factor(df$SMAD4)
  
  jitter <- df %>%
    count(Neighborhood, SMAD4, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$SMAD4 <- factor(jitter$SMAD4, levels = rev(c("1", "0")))
  
  SMAD4S <- jitter %>%
    filter(SMAD4 == "1") %>%
    select(Neighborhood, SMAD4_cnt = prop)
  
  jitter <- jitter %>%
    left_join(SMAD4S, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, SMAD4_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = SMAD4)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No SMAD4 Mutation", "SMAD4 Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with SMAD4 Mutation and Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with SMAD4 Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
smad_neigh_graph(data20_mss_cohort_only)



pik_neigh_graph <- function(df){
  
  df$PIK3CA <- as.factor(df$PIK3CA)
  
  jitter <- df %>%
    count(Neighborhood, PIK3CA, name = "count") %>%
    group_by(Neighborhood) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  
  jitter$PIK3CA <- factor(jitter$PIK3CA, levels = rev(c("1", "0")))
  
  PIK3CAs <- jitter %>%
    filter(PIK3CA == "1") %>%
    select(Neighborhood, PIK3CA_cnt = prop)
  
  jitter <- jitter %>%
    left_join(PIK3CAs, by = "Neighborhood") %>%
    mutate(
      Neighborhood = fct_reorder(Neighborhood, PIK3CA_cnt, .desc = TRUE)
    )
  
  ggplot(jitter, aes(x = Neighborhood, y = prop, fill = PIK3CA)) +
    geom_col() +
    #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("1" = "#E69F00", "0" = "#56B4E9"), labels = c("No PIK3CA Mutation", "PIK3CA Mutation")) +
    #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Distribution of Patients with PIK3CA Mutation and Colorectal Cancer From\nMSK-IMPACT by New York City Neighborhood",
         x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
         y = "Percentage of Patients with PIK3CA Mutation and Colorectal Cancer") +
    #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
    #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
pik_neigh_graph(data20_mss_cohort_only)



stage_diag_neigh_graph <- function(df){
  
    jitter <- df %>%
      count(Neighborhood, Stage_at_Diagnosis_Final, name = "count") %>%
      group_by(Neighborhood) %>%
      mutate(prop = count / sum(count)) %>%
      ungroup()
    
    jitter$Stage_at_Diagnosis_Final <- factor(jitter$Stage_at_Diagnosis_Final, levels = rev(c("Stage_IV", "Stage_III", "Stage_I_or_II")))
    
    stage_diag <- jitter %>%
      filter(Stage_at_Diagnosis_Final == "Stage_IV") %>%
      select(Neighborhood, stage_cnt = prop)
    
    jitter <- jitter %>%
      left_join(stage_diag, by = "Neighborhood") %>%
      mutate(
        Neighborhood = fct_reorder(Neighborhood, stage_cnt, .desc = TRUE)
      )
    
    ggplot(jitter, aes(x = Neighborhood, y = prop, fill = Stage_at_Diagnosis_Final)) +
      geom_col() +
      #geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = c("Stage_IV" = "#E69F00", "Stage_III" = "#56B4E9", "Stage_I_or_II" = "#009E73")) +
      #geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
      theme_minimal() +
      labs(title = "Distribution of Stage at Diagnosis of Patients with Colorectal\nCancer from MSK-IMPACT by New York City Neighborhood",
           x = "Neighborhood of Origin of Patients with Colorectal Cancer from the MSK-IMPACT Cohort",
           y = "Colorectal Cancer Stage at Diagnosis of Patients from MSK-IMPACT") +
      #annotate("text", x = 6, y = 0.95, label = "p-value < 1e-07") +
      #annotate("text", x = 6, y = 0.91, label = "q-value < 1e-07") 
      #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
      #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 90, hjust = 1, vjust = 1))
}
stage_diag_neigh_graph(data20)


