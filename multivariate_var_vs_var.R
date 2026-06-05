#Script Summary:
# 1. Compile a matrix in which each cell corresponds to a p value resulted from testing its respective variables against each other
# 2. Compile a matrix in which each cell corresponds to a q value resulted from testing its respective variables against each other
# 3. Visualize these relationships 


#THIS is the multivariate analysis that Dr. Sanchez Vega requested


#coin library documentation: https://cran.r-project.org/web/packages/coin/vignettes/coin.pdf


#To see the values without running all the code, go to Calculate_p_q_vals -> Multivariate -> var_vs_var_str_by_neigh

library(coin)
library(dplyr)
library(tidyr)
library(ggplot2)


#Import the data
data20 <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")
names(data20)[11] <- "Smoking_History"
names(data20)[13] <- "Yost_Index"
names(data20)[21] <- "Prior_Treatment"

data20_mss_only <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")
names(data20_mss_only)[11] <- "Smoking_History"
names(data20_mss_only)[13] <- "Yost_Index"
names(data20_mss_only)[21] <- "Prior_Treatment"




#make the framework for the p-value matrix

# ------NOTE----- : variable names must match EXACTLY! You can adjust the names of the columns with the following line:
#names(dataframe name)[n] <- "New Column Name"
#where n is the column number
final_pvals <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA"), 
                                     Yost_Index = '', 
                                     Age_At_Diagnosis = '',
                                     Primary_Tumor_Location_Final = '',
                                     Sex = '',
                                     Cohort = '',
                                     Ancestry_Label = '',
                                     First_BMI_Measurement = '',
                                     Prior_Treatment = '',
                                     Smoking_History = '',
                                     KRAS = '',
                                     APC = '',
                                     TP53 = '',
                                     BRAF = '',
                                     SMAD4 = '',
                                     PIK3CA = '')



#function to calculate p-values
lots_of_pvals <- function(df, dfmss){

  #set a seed
  set.seed(42)
  temp <- dfmss
  temp_mss <- dfmss
  
  #turn all categorical variables into factors---the coin library only takes factors
  temp$Neighborhood <- as.factor(temp$Neighborhood)
  temp$Sex <- as.factor(temp$Sex)
  temp$Prior_Treatment <- as.factor(temp$Prior_Treatment)
  temp$Smoking_History <- as.factor(temp$Smoking_History)
  temp$Primary_Tumor_Location_Final <- as.factor(temp$Primary_Tumor_Location_Final)
  temp$Ancestry_Label <- as.factor(temp$Ancestry_Label)
  temp$cohort <- as.factor(temp$cohort)
  
  temp_mss$KRAS <- as.factor(temp_mss$KRAS)
  temp_mss$APC <- as.factor(temp_mss$APC)
  temp_mss$TP53 <- as.factor(temp_mss$TP53)
  temp_mss$BRAF <- as.factor(temp_mss$BRAF)
  temp_mss$SMAD4 <- as.factor(temp_mss$SMAD4)
  temp_mss$PIK3CA <- as.factor(temp_mss$PIK3CA)
  temp_mss$Neighborhood <- as.factor(temp$Neighborhood)
  temp_mss$Sex <- as.factor(temp$Sex)
  temp_mss$Prior_Treatment <- as.factor(temp_mss$Prior_Treatment)
  temp_mss$Smoking_History <- as.factor(temp_mss$Smoking_History)
  temp_mss$Primary_Tumor_Location_Final <- as.factor(temp_mss$Primary_Tumor_Location_Final)
  temp_mss$Ancestry_Label <- as.factor(temp_mss$Ancestry_Label)
  temp_mss$cohort <- as.factor(temp_mss$cohort)
  
  #filter out neighborhoods in list_neighs that did not make the threshold
  
  # ******SPEARMAN******

  #EXAMPLE: We test age at diagnosis against yost index and stratify by neighborhood
  age_yost <- spearman_test(
    age_at_diag ~ Yost_Index | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,2] <<- pvalue(age_yost)
  
  #Now we test yost against age and stratify by neighborhood---yes these values should be the same but for we compute this 
  #for completeness
  yost_age <- spearman_test(
    Yost_Index ~ age_at_diag | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,3] <<- pvalue(yost_age)
  
  yost_bmi <- spearman_test(
    First_BMI_Measurement ~ Yost_Index | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,8] <<- pvalue(yost_bmi)
  
  
  bmi_yost <- spearman_test(
    Yost_Index ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  print(bmi_yost)
  final_pvals[7,2] <<- pvalue(bmi_yost)

  
  age_bmi <- spearman_test(
    age_at_diag ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,8] <<- pvalue(age_bmi)
  
  bmi_age <- spearman_test(
    First_BMI_Measurement ~ age_at_diag | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,3] <<- pvalue(bmi_age)
  
  
  # *******KRUSKAL*******
  
  
  
  
  #yost vs location
  yost_loc <- kruskal_test(
    Yost_Index ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,4] <<- pvalue(yost_loc)
  
  loc_yost <- independence_test(
    Primary_Tumor_Location_Final ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Primary_Tumor_Location_Final" = c("Right" = 1, "Left" = 2, "Rectum" = 3)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,2] <<- pvalue(loc_yost)
    

  
  #yost vs sex
  yost_sex <- kruskal_test(
    Yost_Index ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,5] <<- pvalue(yost_sex)
    
  
  sex_yost <- independence_test(
    Sex ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,2] <<- pvalue(sex_yost)
  

  
  #yost vs ancestry
  yost_ancest <- kruskal_test(
    Yost_Index ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,7] <<- pvalue(yost_ancest)
    
  
  ancest_yost <- independence_test(
    Ancestry_Label ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Ancestry_Label" = c("ADM" = 1, "AFR" = 2, "ASJ" = 3, "EAS" = 4, "EUR" = 5, "NAM" = 6, "SAS" = 7)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,2] <<- pvalue(ancest_yost)

  
  
  #yost vs prior treatment
  yost_prior <- kruskal_test(
    Yost_Index ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,9] <<- pvalue(yost_prior)
  
  
  prior_yost <- independence_test(
    Prior_Treatment ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,2] <<- pvalue(prior_yost)
    
  
  
  #yost vs smoking
  yost_smoking <- kruskal_test(
    Yost_Index ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,10] <<- pvalue(yost_smoking)
  
  
  smoke_yost <- independence_test(
    Smoking_History ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,2] <<- pvalue(smoke_yost)
    

  
  #yost vs kras
  yost_kras <- kruskal_test(
    Yost_Index ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,11] <<- pvalue(yost_kras)
  
  
  kras_yost <- independence_test(
    KRAS ~ Yost_Index | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,2] <<- pvalue(kras_yost)
  
  
  #yost vs apc
  yost_apc <- kruskal_test(
    Yost_Index ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,12] <<- pvalue(yost_apc)
  
  
  apc_yost <- independence_test(
    APC ~ Yost_Index | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,2] <<- pvalue(apc_yost)
  
  
  #yost vs tp53
  yost_tp <- kruskal_test(
    Yost_Index ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,13] <<- pvalue(yost_tp)
  
  
  tp_yost <- independence_test(
    TP53 ~ Yost_Index | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,2] <<- pvalue(tp_yost)
  
  
  #yost vs braf
  yost_braf <- kruskal_test(
    Yost_Index ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,14] <<- pvalue(yost_braf)
  
  
  braf_yost <- independence_test(
    BRAF ~ Yost_Index | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,2] <<- pvalue(braf_yost)
  
  #yost vs smad4
  yost_smad4 <- kruskal_test(
    Yost_Index ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,15] <<- pvalue(yost_smad4)
  
  
  smad4_yost <- independence_test(
    SMAD4 ~ Yost_Index | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,2] <<- pvalue(smad4_yost)
  
  #yost vs pik3ca
  yost_pik <- kruskal_test(
    Yost_Index ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,16] <<- pvalue(yost_pik)
  
  
  pik_yost <- independence_test(
    PIK3CA ~ Yost_Index | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,2] <<- pvalue(pik_yost)
  
  #age vs location
  age_loc <- kruskal_test(
    age_at_diag ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,4] <<- pvalue(age_loc)
  
  
  loc_age <- independence_test(
    Primary_Tumor_Location_Final ~ age_at_diag | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Primary_Tumor_Location_Final" = c("Right" = 1, "Left" = 2, "Rectum" = 3)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,3] <<- pvalue(loc_age)
    
  
  #age vs sex
  age_sex <- kruskal_test(
    age_at_diag ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,5] <<- pvalue(age_sex)
    
  
  sex_age <- independence_test(
    Sex ~ age_at_diag | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,3] <<- pvalue(sex_age)
    
  
  
  #age vs ancestry
  age_ancest <- kruskal_test(
    age_at_diag ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,7] <<- pvalue(age_ancest)
    
  
  ancest_age <- independence_test(
    Ancestry_Label ~ age_at_diag | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Ancestry_Label" = c("ADM" = 1, "AFR" = 2, "ASJ" = 3, "EAS" = 4, "EUR" = 5, "NAM" = 6, "SAS" = 7)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,3] <<- pvalue(ancest_age)
    
  
  #age vs prior treatment
  age_prior <- kruskal_test(
    age_at_diag ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,9] <<- pvalue(age_prior)
    
  prior_age <- independence_test(
    Prior_Treatment ~ age_at_diag | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,3] <<- pvalue(prior_age)
    
  
  
  #age vs smoking
  age_smoking <- kruskal_test(
    age_at_diag ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,10] <<- pvalue(age_smoking)
  
  smoking_age <- independence_test(
    Smoking_History ~ age_at_diag | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,3] <<- pvalue(smoking_age)
  
  
  #age vs kras
  age_kras <- kruskal_test(
    age_at_diag ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,11] <<- pvalue(age_kras)
  
  
  kras_age <- independence_test(
    KRAS ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,3] <<- pvalue(kras_age)
  
  #age vs apc
  age_apc <- kruskal_test(
    age_at_diag ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,12] <<- pvalue(age_apc)
  
  
  apc_age <- independence_test(
    APC ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,3] <<- pvalue(apc_age)
  
  
  #age vs tp53
  age_tp <- kruskal_test(
    age_at_diag ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,13] <<- pvalue(age_tp)
  
  
  tp_age <- independence_test(
    TP53 ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,3] <<- pvalue(tp_age)
  
  
  #age vs braf
  age_braf <- kruskal_test(
    age_at_diag ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,14] <<- pvalue(age_braf)
  
  
  braf_age <- independence_test(
    BRAF ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,3] <<- pvalue(braf_age)
  
  
  #age vs smad4
  age_smad <- kruskal_test(
    age_at_diag ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,15] <<- pvalue(age_smad)
  
  
  smad_age <- independence_test(
    SMAD4 ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,3] <<- pvalue(smad_age)
  
  #age vs pik3ca
  age_pik <- kruskal_test(
    age_at_diag ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,16] <<- pvalue(age_pik)
  
  
  pik_age <- independence_test(
    PIK3CA ~ age_at_diag | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,3] <<- pvalue(pik_age)
    

  
  #bmi vs location
  bmi_loc <- kruskal_test(
    First_BMI_Measurement ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,4] <<- pvalue(bmi_loc)
    
  
  loc_bmi <- independence_test(
    Primary_Tumor_Location_Final ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Primary_Tumor_Location_Final" = c("Right" = 1, "Left" = 2, "Rectum" = 3)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,8] <<- pvalue(loc_bmi)
  
  
  #bmi vs sex
  bmi_sex <- kruskal_test(
    First_BMI_Measurement ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,5] <<- pvalue(bmi_sex)
  
  sex_bmi <- independence_test(
    Sex ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,8] <<- pvalue(sex_bmi)
    
  
  
  #bmi vs ancestry
  bmi_ancest <- kruskal_test(
    First_BMI_Measurement ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,7] <<- pvalue(bmi_ancest)
  
  
  ancest_bmi <- independence_test(
    Ancestry_Label ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    teststat = "scalar",
    scores = list("Ancestry_Label" = c("ADM" = 1, "AFR" = 2, "ASJ" = 3, "EAS" = 4, "EUR" = 5, "NAM" = 6, "SAS" = 7)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,8] <<- pvalue(ancest_bmi)
    

  
  #bmi vs prior treatment
  bmi_prior <- kruskal_test(
    First_BMI_Measurement ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,9] <<- pvalue(bmi_prior)
    
  
  prior_bmi <- independence_test(
    Prior_Treatment ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,8] <<- pvalue(prior_bmi)
    
  
  #bmi vs smoking
  bmi_smoking <- kruskal_test(
    First_BMI_Measurement ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,10] <<- pvalue(bmi_smoking)
  
  
  smoking_bmi <- independence_test(
    Smoking_History ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,8] <<- pvalue(smoking_bmi)
    

  
  #bmi vs kras
  bmi_kras <- kruskal_test(
    First_BMI_Measurement ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,11] <<- pvalue(bmi_kras)
  
  
  kras_bmi <- independence_test(
    KRAS ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,8] <<- pvalue(kras_bmi)
  
  #bmi vs apc
  bmi_apc <- kruskal_test(
    First_BMI_Measurement ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,12] <<- pvalue(bmi_apc)
  
  
  apc_bmi <- independence_test(
    APC ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,8] <<- pvalue(apc_bmi)
  
  
  #bmi vs tp53
  bmi_tp <- kruskal_test(
    First_BMI_Measurement ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,13] <<- pvalue(bmi_tp)
  
  
  tp_bmi <- independence_test(
    TP53 ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,8] <<- pvalue(tp_bmi)
  
  #bmi vs braf
  bmi_braf <- kruskal_test(
    First_BMI_Measurement ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,14] <<- pvalue(bmi_braf)
  
  
  braf_bmi <- independence_test(
    BRAF ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,8] <<- pvalue(braf_bmi)
  
  
  #bmi vs smad4
  bmi_smad <- kruskal_test(
    First_BMI_Measurement ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,15] <<- pvalue(bmi_smad)
  
  
  smad_bmi <- independence_test(
    SMAD4 ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,8] <<- pvalue(smad_bmi)
  
  #bmi vs pik3ca
  bmi_pik <- kruskal_test(
    First_BMI_Measurement ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,16] <<- pvalue(bmi_pik)
  
  
  pik_bmi <- independence_test(
    PIK3CA ~ First_BMI_Measurement | Neighborhood,
    data = temp_mss,
    teststat = "scalar",
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,8] <<- pvalue(pik_bmi)  
  
  # ********CHI^2********
  
  
  results_loc_sex <- independence_test(
    Primary_Tumor_Location_Final ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,5] <<- pvalue(results_loc_sex)
  
  
  results_sex_loc <- independence_test(
    Sex ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,4] <<- pvalue(results_sex_loc)
  
  
  
  ###ASK ABOUT!!!!!!
  #location vs ancestry
  results_loc_ancest <- independence_test(
    Primary_Tumor_Location_Final ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,7] <<- pvalue(results_loc_ancest)
  
  
  results_ancest_loc <- independence_test(
    Ancestry_Label ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,4] <<- pvalue(results_ancest_loc)
  
  
  #location vs prior treatment
  results_loc_priortreat <- independence_test(
    Primary_Tumor_Location_Final ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,9] <<- pvalue(results_loc_priortreat)
    
  
  results_prior_loc <- independence_test(
    Prior_Treatment ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,4] <<- pvalue(results_prior_loc)
  
  
  #location vs smoking
  results_loc_smoke <- independence_test(
    Primary_Tumor_Location_Final ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,10] <<- pvalue(results_loc_smoke)
    
  
  results_smoke_loc <- independence_test(
    Smoking_History ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,4] <<- pvalue(results_smoke_loc)
  
  #location vs kras
  results_loc_kras <- independence_test(
    Primary_Tumor_Location_Final ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,11] <<- pvalue(results_loc_kras)
  
  
  results_kras_loc <- independence_test(
    KRAS ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,4] <<- pvalue(results_kras_loc)
  
  
  #location vs apc
  results_loc_apc <- independence_test(
    Primary_Tumor_Location_Final ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,12] <<- pvalue(results_loc_apc)
  
  
  results_apc_loc <- independence_test(
    APC ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,4] <<- pvalue(results_apc_loc)
  
  #location vs tp53
  results_loc_tp <- independence_test(
    Primary_Tumor_Location_Final ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,13] <<- pvalue(results_loc_tp)
  
  
  results_tp_loc <- independence_test(
    TP53 ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,4] <<- pvalue(results_tp_loc)
  
  #location vs braf
  results_loc_braf <- independence_test(
    Primary_Tumor_Location_Final ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,14] <<- pvalue(results_loc_braf)
  
  
  results_braf_loc <- independence_test(
    BRAF ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,4] <<- pvalue(results_braf_loc)
  
  
  #location vs smad4
  results_loc_smad <- independence_test(
    Primary_Tumor_Location_Final ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,15] <<- pvalue(results_loc_smad)
  
  
  results_smad_loc <- independence_test(
    SMAD4 ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,4] <<- pvalue(results_smad_loc)
  
  
  #location vs pik3ca
  results_loc_pik <- independence_test(
    Primary_Tumor_Location_Final ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,16] <<- pvalue(results_loc_pik)
  
  results_pik_loc <- independence_test(
    PIK3CA ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,4] <<- pvalue(results_pik_loc)
  
  
  #Sex vs Ancestry
  results_sex_ancest <- independence_test(
    Sex ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,7] <<- pvalue(results_sex_ancest)
  
  
  results_ancest_sex <- independence_test(
    Ancestry_Label ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,5] <<- pvalue(results_ancest_sex)
  
  
  
  #sex vs prior treatment
  results_sex_prior <- independence_test(
    Sex ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,9] <<- pvalue(results_sex_prior)
    
  
  
  results_prior_sex <- independence_test(
    Prior_Treatment ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,5] <<- pvalue(results_prior_sex)
  
  
  
  #sex vs smoking
  results_sex_smoke <- independence_test(
    Sex ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,10] <<- pvalue(results_sex_smoke)
  
  
  results_smoke_sex <- independence_test(
    Smoking_History ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,5] <<- pvalue(results_smoke_sex)
  
  
  #sex vs kras
  results_sex_kras <- independence_test(
    Sex ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,11] <<- pvalue(results_sex_kras)
  
  
  results_kras_sex <- independence_test(
    KRAS ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,5] <<- pvalue(results_kras_sex)
    
  
  #sex vs apc
  results_sex_apc <- independence_test(
    Sex ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,12] <<- pvalue(results_sex_apc)
  
  
  results_apc_sex <- independence_test(
    APC ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,5] <<- pvalue(results_apc_sex)
  
  
  #sex vs tp53
  results_sex_tp <- independence_test(
    Sex ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,13] <<- pvalue(results_sex_tp)
  
  
  results_tp_sex <- independence_test(
    TP53 ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,5] <<- pvalue(results_tp_sex)
  
  #sex vs braf
  results_sex_braf <- independence_test(
    Sex ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,14] <<- pvalue(results_sex_braf)
  
  
  results_braf_sex <- independence_test(
    BRAF ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,5] <<- pvalue(results_braf_sex)
  
  
  #sex vs smad4
  results_sex_smad <- independence_test(
    Sex ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,15] <<- pvalue(results_sex_smad)
  
  
  results_smad_sex <- independence_test(
    SMAD4 ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,5] <<- pvalue(results_smad_sex)
  
  
  #sex vs pik
  results_sex_pik <- independence_test(
    Sex ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,16] <<- pvalue(results_sex_pik)
  
  
  results_pik_sex <- independence_test(
    PIK3CA ~ Sex | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,5] <<- pvalue(results_pik_sex)

  
  #prior treatment vs smoking
  results_prior_smoke <- independence_test(
    Prior_Treatment ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,10] <<- pvalue(results_prior_smoke)
  
  
  results_smoke_prior <- independence_test(
    Smoking_History ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,9] <<- pvalue(results_smoke_prior)
    
  
  #prior treatment vs kras
  results_prior_kras <- independence_test(
    Prior_Treatment ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,11] <<- pvalue(results_prior_kras)
  
  
  results_kras_prior <- independence_test(
    KRAS ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,9] <<- pvalue(results_kras_prior)
  
  
  #prior treatment vs apc
  results_prior_apc <- independence_test(
    Prior_Treatment ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,12] <<- pvalue(results_prior_apc)
  
  
  results_apc_prior <- independence_test(
    APC ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,9] <<- pvalue(results_apc_prior)
  
  
  #prior treatment vs tp53
  results_prior_tp <- independence_test(
    Prior_Treatment ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,13] <<- pvalue(results_prior_tp)
  
  
  results_tp_prior <- independence_test(
    TP53 ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,9] <<- pvalue(results_tp_prior)
  
  
  #prior treatment vs braf
  results_prior_braf <- independence_test(
    Prior_Treatment ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,14] <<- pvalue(results_prior_braf)
  
  
  results_braf_prior <- independence_test(
    BRAF ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,9] <<- pvalue(results_braf_prior)
  
  
  #prior treatment vs smad4
  results_prior_smad <- independence_test(
    Prior_Treatment ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,15] <<- pvalue(results_prior_smad)
  
  
  results_smad_prior <- independence_test(
    SMAD4 ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,9] <<- pvalue(results_smad_prior)
  
  
  #prior treatment vs PIK3CA
  results_prior_pik <- independence_test(
    Prior_Treatment ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,16] <<- pvalue(results_prior_pik)
  
  
  results_pik_prior <- independence_test(
    PIK3CA ~ Prior_Treatment | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,9] <<- pvalue(results_pik_prior)
  
  
  #ancestry vs prior treatment
  results_ancest_prior <- independence_test(
    Ancestry_Label ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,9] <<- pvalue(results_ancest_prior)
  
  
  results_prior_ancest <- independence_test(
    Prior_Treatment ~ Ancestry_Label| Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,7] <<- pvalue(results_prior_ancest)
  
  
  
  #ancestry vs smoking
  results_ancest_smoke <- independence_test(
    Ancestry_Label ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,10] <<- pvalue(results_ancest_smoke)
  
  
  results_smoke_ancest <- independence_test(
    Smoking_History ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,7] <<- pvalue(results_smoke_ancest)
  
  
  #ancestry vs kras
  results_ancest_kras <- independence_test(
    Ancestry_Label ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,11] <<- pvalue(results_ancest_kras)
  
  
  results_kras_ancest <- independence_test(
    KRAS ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,7] <<- pvalue(results_kras_ancest)
  
  
  
  #ancestry vs apc
  results_ancest_apc <- independence_test(
    Ancestry_Label ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,12] <<- pvalue(results_ancest_apc)
  
  
  results_apc_ancest <- independence_test(
    APC ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,7] <<- pvalue(results_apc_ancest)
  
  
  #ancestry vs tp53
  results_ancest_tp <- independence_test(
    Ancestry_Label ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,13] <<- pvalue(results_ancest_tp)
  
  
  results_tp_ancest <- independence_test(
    TP53 ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,7] <<- pvalue(results_tp_ancest)
  
  
  #ancestry vs braf
  results_ancest_braf <- independence_test(
    Ancestry_Label ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,14] <<- pvalue(results_ancest_braf)
  
  
  results_braf_ancest <- independence_test(
    BRAF ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,7] <<- pvalue(results_braf_ancest)
  
  
  #ancestry vs smad
  results_ancest_smad <- independence_test(
    Ancestry_Label ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,15] <<- pvalue(results_ancest_smad)
  
  
  results_smad_ancest <- independence_test(
    SMAD4 ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,7] <<- pvalue(results_smad_ancest)
  
  
  #ancestry vs pik3ca
  results_ancest_pik <- independence_test(
    Ancestry_Label ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,16] <<- pvalue(results_ancest_pik)
  
  
  results_pik_ancest <- independence_test(
    PIK3CA ~ Ancestry_Label | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,7] <<- pvalue(results_pik_ancest)
    
  
  #smoking history vs kras
  results_smoke_kras <- independence_test(
    Smoking_History ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,11] <<- pvalue(results_smoke_kras)
  
  
  results_kras_smoke <- independence_test(
    KRAS ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,10] <<- pvalue(results_kras_smoke)
  
  
  #smoking history vs apc
  results_smoke_apc <- independence_test(
    Smoking_History ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,12] <<- pvalue(results_smoke_apc)
  
  
  results_apc_smoke <- independence_test(
    APC ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,10] <<- pvalue(results_apc_smoke)
  
  #smoking history vs tp53
  results_smoke_tp <- independence_test(
    Smoking_History ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,13] <<- pvalue(results_smoke_tp)
  
  
  results_tp_smoke <- independence_test(
    TP53 ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,10] <<- pvalue(results_tp_smoke)
  
  
  #smoking history vs braf
  results_smoke_braf <- independence_test(
    Smoking_History ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,14] <<- pvalue(results_smoke_braf)
  
  
  results_braf_smoke <- independence_test(
    BRAF ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,10] <<- pvalue(results_braf_smoke)
  
  
  #smoking history vs smad4
  results_smoke_smad <- independence_test(
    Smoking_History ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,15] <<- pvalue(results_smoke_smad)
  
  
  results_smad_smoke <- independence_test(
    SMAD4 ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,10] <<- pvalue(results_smad_smoke)
  
  
  #smoking history vs pik
  results_smoke_pik <- independence_test(
    Smoking_History ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,16] <<- pvalue(results_smoke_pik)
  
  
  #smoking history vs pik
  results_pik_smoke <- independence_test(
    PIK3CA ~ Smoking_History | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,10] <<- pvalue(results_pik_smoke)
  
  
  #kras vs apc
  results_kras_apc <- independence_test(
    KRAS ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,12] <<- pvalue(results_kras_apc)
  
  
  results_apc_kras <- independence_test(
    APC ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,11] <<- pvalue(results_apc_kras)
  
  
  #kras vs tp53
  results_kras_tp <- independence_test(
    KRAS ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,13] <<- pvalue(results_kras_tp)
  
  results_tp_kras <- independence_test(
    TP53 ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,11] <<- pvalue(results_tp_kras)
  
  
  #kras vs BRAF
  results_kras_braf <- independence_test(
    KRAS ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,14] <<- pvalue(results_kras_braf)
  
  
  results_braf_kras <- independence_test(
    BRAF ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,11] <<- pvalue(results_braf_kras)
  
  
  #kras vs smad4
  results_kras_smad <- independence_test(
    KRAS ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,15] <<- pvalue(results_kras_smad)
  
  
  results_smad_kras <- independence_test(
    SMAD4 ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,11] <<- pvalue(results_smad_kras)
  
  
  
  #kras vs pik3c4
  results_kras_pik <- independence_test(
    KRAS ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,16] <<- pvalue(results_kras_pik)
  
  
  results_pik_kras <- independence_test(
    PIK3CA ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,11] <<- pvalue(results_pik_kras)
  
  
  #apc vs tp53
  results_apc_tp <- independence_test(
    APC ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,13] <<- pvalue(results_apc_tp)
  
  
  results_tp_apc <- independence_test(
    TP53 ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,12] <<- pvalue(results_tp_apc)
  
  
  #apc vs braf
  results_apc_braf <- independence_test(
    APC ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,14] <<- pvalue(results_apc_braf)
  
  
  results_braf_apc <- independence_test(
    BRAF ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,12] <<- pvalue(results_braf_apc)
  
  
  #apc vs smad4
  results_apc_smad <- independence_test(
    APC ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,15] <<- pvalue(results_apc_smad)
  
  
  results_smad_apc <- independence_test(
    SMAD4 ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,12] <<- pvalue(results_smad_apc)
  
  
  
  #apc vs pik3ca
  results_apc_pik <- independence_test(
    APC ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,16] <<- pvalue(results_apc_pik)
  
  
  results_pik_apc <- independence_test(
    APC ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,12] <<- pvalue(results_pik_apc)
  
  
  #tp53 vs braf
  results_tp_braf <- independence_test(
    TP53 ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,14] <<- pvalue(results_tp_braf)
  
  
  results_braf_tp <- independence_test(
    BRAF ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,13] <<- pvalue(results_braf_tp)
  
  
  #tp53 vs smad4
  results_tp_smad <- independence_test(
    TP53 ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,15] <<- pvalue(results_tp_smad)
  
  
  results_smad_tp <- independence_test(
    SMAD4 ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,13] <<- pvalue(results_smad_tp)
  
  
  
  #tp53 vs pik3ca
  results_tp_pik <- independence_test(
    TP53 ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,16] <<- pvalue(results_tp_pik)
  
  
  results_pik_tp <- independence_test(
    PIK3CA ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,13] <<- pvalue(results_pik_tp)
  
  
  #barf vs smad4
  results_braf_smad <- independence_test(
    BRAF ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,15] <<- pvalue(results_braf_smad)
  
  
  results_smad_braf <- independence_test(
    SMAD4 ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,14] <<- pvalue(results_smad_braf)
  
  
  #barf vs pik3ca
  results_braf_pik <- independence_test(
    BRAF ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,16] <<- pvalue(results_braf_pik)
  
  results_pik_braf <- independence_test(
    PIK3CA ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,14] <<- pvalue(results_pik_braf)
  
  
  #smad4 vs pik3ca
  results_pik_smad <- independence_test(
    PIK3CA ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,15] <<- pvalue(results_pik_smad)
  
  
  results_smad_pik <- independence_test(
    SMAD4 ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,16] <<- pvalue(results_smad_pik)
  
  
  #mss vs yost
  yost_mss <- kruskal_test(
    Yost_Index ~ cohort| Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,6] <<- pvalue(yost_mss)
    
  
  mss_yost <- independence_test(
    cohort ~ Yost_Index | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,2] <<- pvalue(mss_yost)
  

  
  
  #mss vs age
  age_mss <- kruskal_test(
    age_at_diag ~ cohort| Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,6] <<- pvalue(age_mss)
    
  
  mss_age <- independence_test(
    cohort ~ age_at_diag | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,3] <<- pvalue(mss_age)
  
  
  
  #mss vs bmi
  bmi_mss <- kruskal_test(
    First_BMI_Measurement ~ cohort| Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,6] <<- pvalue(bmi_mss)
  
  mss_bmi <- independence_test(
    cohort ~ First_BMI_Measurement | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,8] <<- pvalue(mss_bmi)
  
  
  
  #mss vs location
  cont_loc_mss <- independence_test(
    Primary_Tumor_Location_Final ~ cohort | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,6] <<- pvalue(cont_loc_mss)
  
  
  cont_mss_loc <- independence_test(
    cohort ~ Primary_Tumor_Location_Final | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,4] <<- pvalue(cont_mss_loc)
    
  
  #mss vs sex
  cont_sex_mss <- independence_test(
      Sex ~ cohort | Neighborhood,
      data = temp,
      distribution = approximate(nresample = 1000000)
    )
  final_pvals[4,6] <<- pvalue(cont_sex_mss)
  
  
  cont_mss_sex <- independence_test(
    cohort ~ Sex | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,5] <<- pvalue(cont_mss_sex)
    
  
  #mss vs ancestry
  cont_ancest_mss <- independence_test(
    Ancestry_Label ~ cohort | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,6] <<- pvalue(cont_ancest_mss)
  
  cont_mss_ancest <- independence_test(
    cohort ~ Ancestry_Label | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,7] <<- pvalue(cont_mss_ancest)
  
  
  #mss vs prior treatment
  cont_prior_mss <- independence_test(
    Prior_Treatment ~ cohort | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,6] <<- pvalue(cont_prior_mss)
  
  
  cont_mss_prior <- independence_test(
    cohort ~ Prior_Treatment | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,9] <<- pvalue(cont_mss_prior)
    
  
  #mss vs smoking
  cont_smoke_mss <- independence_test(
    Smoking_History ~ cohort | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,6] <<- pvalue(cont_smoke_mss)
  
  
  cont_mss_smoke <- independence_test(
    cohort ~ Smoking_History | Neighborhood,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,10] <<- pvalue(cont_mss_smoke)
    
  
  
  #mss vs kras
  cont_kras_mss <- independence_test(
    KRAS ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,6] <<- pvalue(cont_kras_mss)
  
  cont_mss_kras <- independence_test(
    cohort ~ KRAS | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,11] <<- pvalue(cont_mss_kras)
  
  
  #mss vs apc
  cont_apc_mss <- independence_test(
    APC ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,6] <<- pvalue(cont_apc_mss)
  
  
  cont_mss_apc <- independence_test(
    cohort ~ APC | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,12] <<- pvalue(cont_mss_apc)
  
  
  #mss vs tp53
  cont_tp_mss <- independence_test(
    TP53 ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,6] <<- pvalue(cont_tp_mss)
  
  cont_mss_tp <- independence_test(
    cohort ~ TP53 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,13] <<- pvalue(cont_mss_tp)
  
  
  #mss vs braf
  cont_braf_mss <- independence_test(
    BRAF ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,6] <<- pvalue(cont_braf_mss)
  
  
  cont_mss_braf <- independence_test(
    cohort ~ BRAF | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,14] <<- pvalue(cont_mss_braf)
  
  #mss vs smad4
  cont_smad_mss <- independence_test(
    SMAD4 ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,6] <<- pvalue(cont_smad_mss)
  
  
  cont_mss_smad <- independence_test(
    cohort ~ SMAD4 | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,15] <<- pvalue(cont_mss_smad)
  
  
  #mss vs pik3ca
  cont_pik_mss <- independence_test(
    PIK3CA ~ cohort | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,6] <<- pvalue(cont_pik_mss)
  
  
  cont_mss_pik <- independence_test(
    cohort ~ PIK3CA | Neighborhood,
    data = temp_mss,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,16] <<- pvalue(cont_mss_pik)
  
}

lots_of_pvals(data20, data20_mss_only)
write.csv(final_pvals, "var_vs_var_pvals.csv", row.names = FALSE)

#make the framework for the q-value matrix
final_qvals <- data.frame(variable_1 = c("Yost_Index", "Age_At_Diagnosis", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "First_BMI_Measurement", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA"), 
                          Yost_Index = '', 
                          Age_At_Diagnosis = '',
                          Primary_Tumor_Location_Final = '',
                          Sex = '',
                          Cohort = '',
                          Ancestry_Label = '',
                          First_BMI_Measurement = '',
                          Prior_Treatment = '',
                          Smoking_History = '',
                          KRAS = '',
                          APC = '',
                          TP53 = '',
                          BRAF = '',
                          SMAD4 = '',
                          PIK3CA = '')



#function to adjust p vals
qvals <- function(df){
  
  #squish matrix into vector while maintaining order
  pval_vector <- df %>%
    select(-1) %>%
    pivot_longer(everything(), values_to = "val") %>%
    filter(!is.na(val), val != "") %>%
    pull(val)
  
  #ensure all pvals are numbers, not strings
  pval_vector <- as.numeric(pval_vector)

  #adjust!
  qval_vector <- p.adjust(pval_vector, method = "BH")

  #put in q val matrix while maintaining order
  final_qvals[1,3] <<- qval_vector[1]
  final_qvals[1,4] <<- qval_vector[2]
  final_qvals[1,5] <<- qval_vector[3]
  final_qvals[1,6] <<- qval_vector[4]
  final_qvals[1,7] <<- qval_vector[5]
  final_qvals[1,8] <<- qval_vector[6]
  final_qvals[1,9] <<- qval_vector[7]
  final_qvals[1,10] <<- qval_vector[8]
  final_qvals[1,11] <<- qval_vector[9]
  final_qvals[1,12] <<- qval_vector[10]
  final_qvals[1,13] <<- qval_vector[11]
  final_qvals[1,14] <<- qval_vector[12]
  final_qvals[1,15] <<- qval_vector[13]
  final_qvals[1,16] <<- qval_vector[14]
  final_qvals[2,2] <<- qval_vector[15]
  final_qvals[2,4] <<- qval_vector[16]
  final_qvals[2,5] <<- qval_vector[17]
  final_qvals[2,6] <<- qval_vector[18]
  final_qvals[2,7] <<- qval_vector[19]
  final_qvals[2,8] <<- qval_vector[20]
  final_qvals[2,9] <<- qval_vector[21]
  final_qvals[2,10] <<- qval_vector[22]
  final_qvals[2,11] <<- qval_vector[23]
  final_qvals[2,12] <<- qval_vector[24]
  final_qvals[2,13] <<- qval_vector[25]
  final_qvals[2,14] <<- qval_vector[26]
  final_qvals[2,15] <<- qval_vector[27]
  final_qvals[2,16] <<- qval_vector[28]
  final_qvals[3,2] <<- qval_vector[29]
  final_qvals[3,3] <<- qval_vector[30]
  final_qvals[3,5] <<- qval_vector[31]
  final_qvals[3,6] <<- qval_vector[32]
  final_qvals[3,7] <<- qval_vector[33]
  final_qvals[3,8] <<- qval_vector[34]
  final_qvals[3,9] <<- qval_vector[35]
  final_qvals[3,10] <<- qval_vector[36]
  final_qvals[3,11] <<- qval_vector[37]
  final_qvals[3,12] <<- qval_vector[38]
  final_qvals[3,13] <<- qval_vector[39]
  final_qvals[3,14] <<- qval_vector[40]
  final_qvals[3,15] <<- qval_vector[41]
  final_qvals[3,16] <<- qval_vector[42]
  final_qvals[4,2] <<- qval_vector[43]
  final_qvals[4,3] <<- qval_vector[44]
  final_qvals[4,4] <<- qval_vector[45]
  final_qvals[4,6] <<- qval_vector[46]
  final_qvals[4,7] <<- qval_vector[47]
  final_qvals[4,8] <<- qval_vector[48]
  final_qvals[4,9] <<- qval_vector[49]
  final_qvals[4,10] <<- qval_vector[50]
  final_qvals[5,2] <<- qval_vector[51]
  final_qvals[5,3] <<- qval_vector[52]
  final_qvals[5,4] <<- qval_vector[53]
  final_qvals[5,5] <<- qval_vector[54]
  final_qvals[5,7] <<- qval_vector[55]
  final_qvals[5,8] <<- qval_vector[56]
  final_qvals[5,9] <<- qval_vector[57]
  final_qvals[5,10] <<- qval_vector[58]
  final_qvals[5,11] <<- qval_vector[59]
  final_qvals[5,12] <<- qval_vector[60]
  final_qvals[5,13] <<- qval_vector[61]
  final_qvals[5,14] <<- qval_vector[62]
  final_qvals[5,15] <<- qval_vector[63]
  final_qvals[5,16] <<- qval_vector[64]
  final_qvals[6,2] <<- qval_vector[65]
  final_qvals[6,3] <<- qval_vector[66]
  final_qvals[6,4] <<- qval_vector[67]
  final_qvals[6,5] <<- qval_vector[68]
  final_qvals[6,6] <<- qval_vector[69]
  final_qvals[6,8] <<- qval_vector[70]
  final_qvals[6,9] <<- qval_vector[71]
  final_qvals[6,10] <<- qval_vector[72]
  final_qvals[6,11] <<- qval_vector[73]
  final_qvals[6,12] <<- qval_vector[74]
  final_qvals[6,13] <<- qval_vector[75]
  final_qvals[6,14] <<- qval_vector[76]
  final_qvals[6,15] <<- qval_vector[77]
  final_qvals[6,16] <<- qval_vector[78]
  final_qvals[7,2] <<- qval_vector[79]
  final_qvals[7,3] <<- qval_vector[80]
  final_qvals[7,4] <<- qval_vector[81]
  final_qvals[7,5] <<- qval_vector[82]
  final_qvals[7,6] <<- qval_vector[83]
  final_qvals[7,7] <<- qval_vector[84]
  final_qvals[7,9] <<- qval_vector[85]
  final_qvals[7,10] <<- qval_vector[86]
  final_qvals[7,11] <<- qval_vector[87]
  final_qvals[7,12] <<- qval_vector[88]
  final_qvals[7,13] <<- qval_vector[89]
  final_qvals[7,14] <<- qval_vector[90]
  final_qvals[7,15] <<- qval_vector[91]
  final_qvals[7,16] <<- qval_vector[92]
  final_qvals[8,2] <<- qval_vector[93]
  final_qvals[8,3] <<- qval_vector[94]
  final_qvals[8,4] <<- qval_vector[95]
  final_qvals[8,5] <<- qval_vector[96]
  final_qvals[8,6] <<- qval_vector[97]
  final_qvals[8,7] <<- qval_vector[98]
  final_qvals[8,8] <<- qval_vector[99]
  final_qvals[8,10] <<- qval_vector[100]
  final_qvals[8,11] <<- qval_vector[101]
  final_qvals[8,12] <<- qval_vector[102]
  final_qvals[8,13] <<- qval_vector[103]
  final_qvals[8,14] <<- qval_vector[104]
  final_qvals[8,15] <<- qval_vector[105]
  final_qvals[8,16] <<- qval_vector[106]
  final_qvals[9,2] <<- qval_vector[107]
  final_qvals[9,3] <<- qval_vector[108]
  final_qvals[9,4] <<- qval_vector[109]
  final_qvals[9,5] <<- qval_vector[110]
  final_qvals[9,6] <<- qval_vector[111]
  final_qvals[9,7] <<- qval_vector[112]
  final_qvals[9,8] <<- qval_vector[113]
  final_qvals[9,9] <<- qval_vector[114]
  final_qvals[9,11] <<- qval_vector[115]
  final_qvals[9,12] <<- qval_vector[116]
  final_qvals[9,13] <<- qval_vector[117]
  final_qvals[9,14] <<- qval_vector[118]
  final_qvals[9,15] <<- qval_vector[119]
  final_qvals[9,16] <<- qval_vector[120]
  final_qvals[10,2] <<- qval_vector[121]
  final_qvals[10,3] <<- qval_vector[122]
  final_qvals[10,4] <<- qval_vector[123]
  final_qvals[10,5] <<- qval_vector[124]
  final_qvals[10,7] <<- qval_vector[125]
  final_qvals[10,8] <<- qval_vector[126]
  final_qvals[10,9] <<- qval_vector[127]
  final_qvals[10,10] <<- qval_vector[128]
  final_qvals[10,12] <<- qval_vector[129]
  final_qvals[10,13] <<- qval_vector[130]
  final_qvals[10,14] <<- qval_vector[131]
  final_qvals[10,15] <<- qval_vector[132]
  final_qvals[10,16] <<- qval_vector[133]
  final_qvals[11,2] <<- qval_vector[134]
  final_qvals[11,3] <<- qval_vector[135]
  final_qvals[11,4] <<- qval_vector[136]
  final_qvals[11,5] <<- qval_vector[137]
  final_qvals[11,7] <<- qval_vector[138]
  final_qvals[11,8] <<- qval_vector[139]
  final_qvals[11,9] <<- qval_vector[140]
  final_qvals[11,10] <<- qval_vector[141]
  final_qvals[11,11] <<- qval_vector[142]
  final_qvals[11,13] <<- qval_vector[143]
  final_qvals[11,14] <<- qval_vector[144]
  final_qvals[11,15] <<- qval_vector[145]
  final_qvals[11,16] <<- qval_vector[146]
  final_qvals[12,2] <<- qval_vector[147]
  final_qvals[12,3] <<- qval_vector[148]
  final_qvals[12,4] <<- qval_vector[149]
  final_qvals[12,5] <<- qval_vector[150]
  final_qvals[12,7] <<- qval_vector[151]
  final_qvals[12,8] <<- qval_vector[152]
  final_qvals[12,9] <<- qval_vector[153]
  final_qvals[12,10] <<- qval_vector[154]
  final_qvals[12,11] <<- qval_vector[155]
  final_qvals[12,12] <<- qval_vector[156]
  final_qvals[12,14] <<- qval_vector[157]
  final_qvals[12,15] <<- qval_vector[158]
  final_qvals[12,16] <<- qval_vector[159]
  final_qvals[13,2] <<- qval_vector[160]
  final_qvals[13,3] <<- qval_vector[161]
  final_qvals[13,4] <<- qval_vector[162]
  final_qvals[13,5] <<- qval_vector[163]
  final_qvals[13,7] <<- qval_vector[164]
  final_qvals[13,8] <<- qval_vector[165]
  final_qvals[13,9] <<- qval_vector[166]
  final_qvals[13,10] <<- qval_vector[167]
  final_qvals[13,11] <<- qval_vector[168]
  final_qvals[13,12] <<- qval_vector[169]
  final_qvals[13,13] <<- qval_vector[170]
  final_qvals[13,15] <<- qval_vector[171]
  final_qvals[13,16] <<- qval_vector[172]
  final_qvals[14,2] <<- qval_vector[173]
  final_qvals[14,3] <<- qval_vector[174]
  final_qvals[14,4] <<- qval_vector[175]
  final_qvals[14,5] <<- qval_vector[176]
  final_qvals[14,7] <<- qval_vector[177]
  final_qvals[14,8] <<- qval_vector[178]
  final_qvals[14,9] <<- qval_vector[179]
  final_qvals[14,10] <<- qval_vector[180]
  final_qvals[14,11] <<- qval_vector[181]
  final_qvals[14,12] <<- qval_vector[182]
  final_qvals[14,13] <<- qval_vector[183]
  final_qvals[14,14] <<- qval_vector[184]
  final_qvals[14,16] <<- qval_vector[185]
  final_qvals[15,2] <<- qval_vector[186]
  final_qvals[15,3] <<- qval_vector[187]
  final_qvals[15,4] <<- qval_vector[188]
  final_qvals[15,5] <<- qval_vector[189]
  final_qvals[15,7] <<- qval_vector[190]
  final_qvals[15,8] <<- qval_vector[191]
  final_qvals[15,9] <<- qval_vector[192]
  final_qvals[15,10] <<- qval_vector[193]
  final_qvals[15,11] <<- qval_vector[194]
  final_qvals[15,12] <<- qval_vector[195]
  final_qvals[15,13] <<- qval_vector[196]
  final_qvals[15,14] <<- qval_vector[197]
  final_qvals[15,15] <<- qval_vector[198]
}

qvals(pvals_20)
write.csv(final_qvals, "var_vs_var_qvals.csv", row.names = FALSE)








#VISUALIZATIONS






#these are basic scatterplots:


bmi_age_graph <- function(df){

  age_vector <- df %>%
    group_by(Neighborhood) %>%
    summarise(median_age = median(age_at_diag, na.rm = TRUE))
  
  bmi_vector <- df %>%
    group_by(Neighborhood) %>%
    summarise(median_bmi = median(First_BMI_Measurement, na.rm = TRUE)) 

  merge <- merge(bmi_vector, age_vector, by = "Neighborhood")
  
  ggplot(merge) +
    aes(x = median_bmi, y = median_age, label = Neighborhood) +
    geom_point(color = "#e33232") +
    labs(
      x = "Median First Body Mass Index Measurement",
      y = "Median Age at Diagnosis of Patients at MSK",
      title = "Median Age at Diagnosis and Median First Body Mass Index Measurement of Patients with\nColorectal Cancer from the MSK-Impact Data Set, Grouped by New York City\nNeighborhood"
    ) +
    geom_point() + geom_text(aes(label=ifelse(median_bmi<24,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    geom_point() + geom_text(aes(label=ifelse(median_bmi>28.5,as.character(Neighborhood),'')), hjust = 0.8,vjust = -0.8) +
    geom_point() + geom_text(aes(label=ifelse(median_age>62,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    geom_point() + geom_text(aes(label=ifelse(median_age<49,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 24.25, y = 51, label = "p-value = 0.004077") +
    annotate("text", x = 24.25, y = 50.5, label = "q-value = 0.0171234") +
    theme_minimal()
}
bmi_age_graph(data20)


#bmi vs yost index
bmi_yost_graph <- function(df){
  yost_vector <- df %>%
    group_by(Neighborhood) %>%
    summarise(median_yost = median(Yost_Index, na.rm = TRUE))
  
  bmi_vector <- df %>%
    group_by(Neighborhood) %>%
    summarise(median_bmi = median(First_BMI_Measurement, na.rm = TRUE)) 
  
  merge <- merge(bmi_vector, yost_vector, by = "Neighborhood")
  
  ggplot(merge) +
    aes(x = median_bmi, y = median_yost, label = Neighborhood) +
    geom_point(color = "#e33232") +
    labs(
      x = "Median First Body Mass Index Measurement",
      y = "Median Yost Index of Patients at MSK by Neighborhood",
      title = "Median Age at Diagnosis and Median First Body Mass Index Measurement of Patients with\nColorectal Cancer from the MSK-Impact Data Set, Grouped by New York City\nNeighborhood"
    ) +
    geom_point() + geom_text(aes(label=ifelse(median_bmi==23.9,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    geom_point() + geom_text(aes(label=ifelse(median_bmi==24.05,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    geom_point() + geom_text(aes(label=ifelse(median_bmi==24.7,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #	Gramercy Park - Murray Hill , Lower Manhattan , Upper East Side
    annotate("text", x = 24.75, y = 51, label = "p-value = 0.004564") +
    annotate("text", x = 24.75, y = 49, label = "q-value = 0.018432") +
    theme_minimal()
}
bmi_yost_graph(data20)



#bmi vs yost index: this one isn't grouped by medians
bmi_yost_graph_all <- function(df){

  ggplot(df) +
    aes(x = First_BMI_Measurement, y = Yost_Index) +
    geom_point(color = "#e33232") +
    labs(
      x = "First Body Mass Index Measurement",
      y = "Yost Index of Patients at MSK by Neighborhood",
      title = "Age at Diagnosis and First Body Mass Index Measurement of Patients with\nColorectal Cancer from the MSK-Impact Data Set, Grouped by New York City\nNeighborhood"
    ) +
    #geom_point() + geom_text(aes(label=ifelse(median_bmi==23.9,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_bmi==24.05,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_bmi==24.7,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #	Gramercy Park - Murray Hill , Lower Manhattan , Upper East Side
    annotate("text", x = 16, y = 93, label = "Rho = 0.01105546") +
    annotate("text", x = 16, y = 87, label = "p-value = 0.7263") +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    
    theme_minimal()
}
bmi_yost_graph_all(data20)




#Categorical bocx plots


#For help understanding the interpretation of these plots, open *this* function:
age_loc_graph <- function(df){
  #x axis has 3 cats, each being from Primary Tumor Loc. Each point is one neighborhood. The height of each point is the
  #median age at diagnosis of patients within that neighborhood AND with left, right, or rectum primary locs. 
  
  #the line in the box plot is the median age for patients in the left, rectum, or right categories. 
  
  
  #calculate the median age per primary tumor loc. ---this dictates the height of each point (a neighborhood) on the plot
  age_loc_vector <- df %>%
    group_by(Neighborhood, Primary_Tumor_Location_Final) %>%
    summarise(median_age = median(age_at_diag, na.rm = TRUE))
  
  ggplot(age_loc_vector, aes(Primary_Tumor_Location_Final, median_age, fill = Primary_Tumor_Location_Final)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9", "Rectum" = "#009E73")) +
    labs(title = "Median Age at Diagnosis by Neighborhood, grouped by Primary Tumor\nLocation of Patients with Colorectal Cancer from the MSK-Impact Data Set",
         x = "Primary Tumor Location",
         y = "Median Age at Diagnosis of Patients at MSK by Neighborhood") +
    theme_minimal() +
    #geom_point() + geom_text(aes(label=ifelse(median_age>68,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    geom_point() + geom_text(aes(label=ifelse(median_age<45,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 1, y = 70, label = "p-value = 9e-06") +
    annotate("text", x = 1, y = 69, label = "q-value = 7.269231e-05") +
    theme(legend.position = "none")
  
}
age_loc_graph(data20)


yost_ancest_graph <- function(df){
  yost_ancest_vector <- df %>%
    group_by(Neighborhood, Ancestry_Label) %>%
    summarise(median_yost = median(Yost_Index, na.rm = TRUE))
  
  ggplot(yost_ancest_vector, aes(Ancestry_Label, median_yost, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    labs(title = "Median Yost Index by Neighborhood, grouped by Ancestry of Patients with\nColorectal Cancer from the MSK-Impact Data Set",
         x = "Ancestry Labels of Patients at MSK",
         y = "Median Yost Index") +
    theme_minimal() +
    geom_point() + geom_text(aes(label=ifelse(median_yost == 78, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4) +
    geom_point() + geom_text(aes(label=ifelse(median_yost == 77.5, as.character(Neighborhood),'')), hjust = 0.5,vjust = 1.5) +
    geom_point() + geom_text(aes(label=ifelse(median_yost == 79, as.character(Neighborhood),'')), hjust = 1,vjust = -0.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 95, as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 1, y = 86, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 83, label = "q-value < 1e-07") +
    theme(legend.position = "none")
  
}
yost_ancest_graph(data20)


age_smoke_graph <- function(df){
  age_smoke_vector <- df %>%
    group_by(Neighborhood, Smoking_History) %>%
    summarise(median_age = median(age_at_diag, na.rm = TRUE))
  
  ggplot(age_smoke_vector, aes(Smoking_History, median_age, fill = Smoking_History)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Former/Current Smoker" = "#8b1ff0", "Never" = "#009E73")) +
    labs(title = "Median Age at Diagnosis by Neighborhood, grouped by Sex of Patients with\nColorectal Cancer from the MSK-Impact Data Set",
         x = "Former/Current Smoker Patients at MSK",
         y = "Median Age at Diagnosis") +
    theme_minimal() +
    geom_point() + geom_text(aes(label=ifelse(median_age<45,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_age<45,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 1, y = 70, label = "p-value = 5.8e-05") + 
    annotate("text", x = 1, y = 69, label = "q-value = 0.000358") +
    theme(legend.position = "none")
  
}
age_smoke_graph(data20)


ancest_locLeft_graph <- function(df){
    
  jitter <- df %>%
    group_by(Neighborhood, Ancestry_Label) %>% 
    summarize(
      mean_left = mean(Primary_Tumor_Location_Final == "Left"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Ancestry_Label, y = mean_left, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with Left-Sided Colorectal\nCancer and Ancestry of Patients with Colorectal Cancer, Layered Underneath\nthe Percentage of Patients with Left-Sided Colorectal Cancer per Ancestry\nLabel and Neighborhood",
         x = "Ancestry Label of Patients at MSK",
         y = "Percentage of Patients with Left-Sided Colorectal Cancer") +
    annotate("text", x = 1.25, y = 0.9, label = "p-value = 0.030567") +
    annotate("text", x = 1.25, y = 0.85, label = "q-value = 0.100298") +
    theme(legend.position = "none")
    
}
ancest_locLeft_graph(data20)


ancest_locRight_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Ancestry_Label) %>% 
    summarize(
      mean_left = mean(Primary_Tumor_Location_Final == "Right"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Ancestry_Label, y = mean_left, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with Right-Sided Colorectal\nCancer and Ancestry of Patients with Colorectal Cancer, Layered Underneath\nthe Percentage of Patients with Right-Sided Colorectal Cancer per Ancestry\nLabel and Neighborhood",
         x = "Ancestry Label of Patients at MSK",
         y = "Percentage of Patients with Right-Sided Colorectal Cancer") +
    theme(legend.position = "none")
  
}
ancest_locRight_graph(data20)



ancest_locRectum_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Ancestry_Label) %>% 
    summarize(
      mean_left = mean(Primary_Tumor_Location_Final == "Rectum"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Ancestry_Label, y = mean_left, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with Rectum-Sided Colorectal\nCancer and Ancestry of Patients with Colorectal Cancer, Layered Underneath\nthe Percentage of Patients with Rectum-Sided Colorectal Cancer per Ancestry\nLabel and Neighborhood",
         x = "Ancestry Label of Patients at MSK",
         y = "Percentage of Patients with Rectum-Sided Colorectal Cancer") +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    theme(legend.position = "none")
  
}
ancest_locRectum_graph(data20)



loc_MSS_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Primary_Tumor_Location_Final) %>% 
    summarize(
      mean_mss = mean(cohort == "MSS"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Primary_Tumor_Location_Final, y = mean_mss, fill = Primary_Tumor_Location_Final)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9", "Rectum" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with Microsatellite Stable\nColorectal Cancer and Primary Tumor Location of Patients with Colorectal\nCancer, Layered Underneath the Percentage of Patients with Microsatellite\nStable Colorectal Cancer per Primary Tumor Location and Neighborhood",
         x = "Primary Tumor Location of Patients at MSK",
         y = "Percentage of Patients with Microsatellite Stable Colorectal Cancer") +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    annotate("text", x = 1, y = 0.7, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.68, label = "q-value < 1e-07") +
    theme(legend.position = "none") 
  
}
loc_MSS_graph(data20)


ancest_smoke_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Ancestry_Label) %>% 
    summarize(
      mean_smoke = mean(Smoking_History == "Former/Current Smoker"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Ancestry_Label, y = mean_smoke, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Former/Current Smoker Colorectal Cancer\nPatients and Ancestry of Patients, Layered Underneath the Percentage of\nFormer/Current Smoker Patients with Colorectal Cancer per Ancestry Label\nand Neighborhood",
         x = "Ancestry Label of Patients at MSK",
         y = "Percentage of Former/Current Smoker Patients with Colorectal Cancer") +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_left == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    annotate("text", x = 1.25, y = 0.95, label = "p-value = 0.000543") +
    annotate("text", x = 1.25, y = 0.9, label = "q-value = 0.002715") +
    theme(legend.position = "none")
  
}
ancest_smoke_graph(data20)



sex_smoke_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Sex) %>% 
    summarize(
      mean_smoke = mean(Smoking_History == "Former/Current Smoker"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Sex, y = mean_smoke, fill = Sex)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Male" = "#E69F00", "Female" = "#56B4E9")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Former/Current Smoker Colorectal Cancer\nPatients and Sex of Patients, Layered Underneath the Percentage of\nFormer/Current Smoker Patients with Colorectal Cancer per Sex and\nNeighborhood",
         x = "Sex of Patients at MSK",
         y = "Percentage of Former/Current Smoker Patients with Colorectal Cancer") +
    geom_point() + geom_text(aes(label=ifelse(mean_smoke > 0.7, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    geom_point() + geom_text(aes(label=ifelse(mean_smoke < 0.05, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 3) +
    annotate("text", x = 1, y = 0.71, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.68, label = "q-value < 1e-07") +
    theme(legend.position = "none")
  
}
sex_smoke_graph(data20)



sex_mss_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Sex) %>% 
    summarize(
      mean_mss = mean(cohort == "MSS"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Sex, y = mean_mss, fill = Sex)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Male" = "#E69F00", "Female" = "#56B4E9")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with Microsatellite Stable\nColorectal Cancer and Sex of Patients with Colorectal Cancer, Layered\nUnderneath the Percentage of Patients with Microsatellite Stable Colorectal\nCancer per Sex and Neighborhood",
         x = "Sex of Patients at MSK",
         y = "Percentage of Patients with Microsatellite Stable Colorectal Cancer") +
    geom_point() + geom_text(aes(label=ifelse(mean_mss < 0.7, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_smoke < 0.05, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 3) +
    annotate("text", x = 2, y = 0.74, label = "p-value = 0.043368") +
    annotate("text", x = 2, y = 0.725, label = "q-value: 0.130578") +
    theme(legend.position = "none")
  
}
sex_mss_graph(data20)



loc_kras_graph <- function(df){
  
  jitter <- data20 %>%
    group_by(Neighborhood, Primary_Tumor_Location_Final) %>% 
    summarize(
      mean_kras = mean(KRAS == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Primary_Tumor_Location_Final, y = mean_kras, fill = Primary_Tumor_Location_Final)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9", "Rectum" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with KRAS Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with KRAS Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Primary Tumor Location of Patients at MSK",
         y = "Percentage of Patients with KRAS Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_kras == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 4) +
    geom_point() + geom_text(aes(label=ifelse(mean_kras == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    annotate("text", x = 1, y = 0.9, label = "p-value = 0.000114") +
    annotate("text", x = 1, y = 0.85, label = "q-value = 0.000665") +
    theme(legend.position = "none") 
  
}
loc_kras_graph(data_mss_only)



mss_kras_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_kras = mean(KRAS == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_kras, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with KRAS Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with KRAS Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with KRAS Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_kras == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    geom_point() + geom_text(aes(label=ifelse(mean_kras > 0.625, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_kras < 0.25, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #mss lowest point is Gramcery Park - Murray Hill
    annotate("text", x = 1, y = 0.9, label = "p-value = 0.000383") +
    annotate("text", x = 1, y = 0.85, label = "q-value = 0.002011") +
    theme(legend.position = "none") 
  
}
mss_kras_graph(data_mss_only)



yost_apc_graph <- function(df){
  labelss <- c("No APC Mutation", "APC Mutation")
  df$APC <- as.character(df$APC)
  yost_apc_vector <- df %>%
    group_by(Neighborhood, APC) %>%
    summarise(median_yost = median(Yost_Index, na.rm = TRUE))
  
  ggplot(yost_apc_vector, aes(APC, median_yost, fill = APC)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("1" = "#8b1ff0", "0" = "#009E73")) +
    labs(title = "Median Yost Index by Neighborhood, grouped by APC Mutation\nClassification of Patients with Colorectal Cancer from the MSK-Impact\nData Set",
         x = "APC Mutation Classification of Patients with Colorectal Cancer from MSK-IMPACT",
         y = "Median Yost Index of Patients with Colorectal Cancer from MSK-IMPACT") +
    theme_minimal() +
    scale_x_discrete(labels = labelss) +
    geom_point() + geom_text(aes(label=ifelse(median_yost > 68, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 77.5, as.character(Neighborhood),'')), hjust = 0.5,vjust = 1.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 79, as.character(Neighborhood),'')), hjust = 1,vjust = -0.5) +
    #geom_point() + geom_text(aes(label=ifelse(median_yost == 95, as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    annotate("text", x = 2, y = 68, label = "p-value = 0.009079") +
    annotate("text", x = 2, y = 65, label = "q-value = 0.035307") +
    theme(legend.position = "none")
  
}
yost_apc_graph(data_mss_only)



mss_apc_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_apc = mean(APC == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_apc, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with APC Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with APC Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with APC Mutation") +
    geom_point() + geom_text(aes(label=ifelse(mean_apc == 1, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 4) +
    #geom_point() + geom_text(aes(label=ifelse(mean_apc == 0, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    #non mss smallest point is bed stuy - crown heights AND Southwest Queens
    #mss lowest point is Gramcery Park - Murray Hill
    annotate("text", x = 1, y = 0.25, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.2, label = "q-value < 1e-07") +
    theme(legend.position = "none") 
  
}
mss_apc_graph(data_mss_only)



age_tp_graph <- function(df){
  labelss <- c("No TP53 Mutation", "TP53 Mutation")
  df$TP53 <- as.character(df$TP53)
  age_tp_vector <- df %>%
    group_by(Neighborhood, TP53) %>%
    summarise(median_age = median(age_at_diag, na.rm = TRUE))
  
  ggplot(age_tp_vector, aes(TP53, median_age, fill = TP53)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("1" = "#8b1ff0", "0" = "#009E73")) +
    labs(title = "Median Age at Diagnosis by Neighborhood, grouped by TP53 Mutation\nClassification of Patients with Colorectal Cancer from the MSK-Impact Data\nSet",
         x = "Percentage of Patients with TP53 Mutation",
         y = "Median Age at Diagnosis") +
    scale_x_discrete(labels = labelss) +
    theme_minimal() +
    geom_point() + geom_text(aes(label=ifelse(median_age>70,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_age>65,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #tp53 mut highest point : union square - lower east side
    annotate("text", x = 2, y = 75, label = "p-value = 0.012445") + 
    annotate("text", x = 2, y = 73, label = "q-value = 0.046669") +
    theme(legend.position = "none")
}
age_tp_graph(data_mss_only)



age_tp_graph <- function(df){
  labelss <- c("No TP53 Mutation", "TP53 Mutation")
  df$TP53 <- as.character(df$TP53)
  age_tp_vector <- df %>%
    group_by(Neighborhood, TP53) %>%
    summarise(median_age = median(age_at_diag, na.rm = TRUE))
  
  ggplot(age_tp_vector, aes(TP53, median_age, fill = TP53)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("1" = "#8b1ff0", "0" = "#009E73")) +
    labs(title = "Median Age at Diagnosis by Neighborhood, grouped by TP53 Mutation\nClassification of Patients with Colorectal Cancer from the MSK-Impact Data\nSet",
         x = "Percentage of Patients with TP53 Mutation",
         y = "Median Age at Diagnosis") +
    scale_x_discrete(labels = labelss) +
    theme_minimal() +
    geom_point() + geom_text(aes(label=ifelse(median_age>70,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_age>65,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #tp53 mut highest point : union square - lower east side
    annotate("text", x = 2, y = 75, label = "p-value = 0.012445") + 
    annotate("text", x = 2, y = 73, label = "q-value = 0.046669") +
    theme(legend.position = "none")
}
age_tp_graph(data_mss_only)



loc_tp_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Primary_Tumor_Location_Final) %>% 
    summarize(
      mean_tp = mean(TP53 == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Primary_Tumor_Location_Final, y = mean_tp, fill = Primary_Tumor_Location_Final)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9", "Rectum" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with TP53 Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with TP53 Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Primary Tumor Location of Patients at MSK",
         y = "Percentage of Patients with TP53 Mutation") +
    geom_point() + geom_text(aes(label=ifelse(mean_tp < 0.45, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 4) +
    #minimum right vals = Rockaway and willowbrook
    annotate("text", x = 1, y = 0.5, label = "p-value = 0.003034") +
    annotate("text", x = 1, y = 0.47, label = "q-value = 0.013851") +
    theme(legend.position = "none") 
  
}
loc_tp_graph(data_mss_only)



mss_tp_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_tp = mean(TP53 == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_tp, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with TP53 Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with TP53 Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with TP53 Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_tp == 1, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 4) +
    #highest points for non-mss are bed stuy -crown heights and rockaway
    #mss lowest point is Gramcery Park - Murray Hill
    annotate("text", x = 1, y = 0.25, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.2, label = "q-value < 1e-07") +
    theme(legend.position = "none") 
  
}
mss_tp_graph(data_mss_only)




loc_braf_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, Primary_Tumor_Location_Final) %>% 
    summarize(
      mean_braf = mean(BRAF == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Primary_Tumor_Location_Final, y = mean_braf, fill = Primary_Tumor_Location_Final)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("Left" = "#E69F00", "Right" = "#56B4E9", "Rectum" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with BRAF Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with BRAF Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Primary Tumor Location of Patients at MSK",
         y = "Percentage of Patients with BRAF Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_braf > 0.3, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 4) +
    #maximum right vals = stapleton st g and upper east side
    annotate("text", x = 1, y = 0.5, label = "p-value = 1.7e-05") +
    annotate("text", x = 1, y = 0.47, label = "q-value = 0.000119") +
    theme(legend.position = "none") 
  
}
loc_braf_graph(data_mss_only)


mss_braf_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_braf= mean(BRAF == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_braf, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with BRAF Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with BRAF Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with BRAF Mutation") +
    geom_point() + geom_text(aes(label=ifelse(mean_braf > 0.8, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    annotate("text", x = 1, y = 0.9, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.85, label = "q-value < 1e-07") +
    theme(legend.position = "none") 
  
}
mss_braf_graph(data_mss_only)




mss_smad_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_smad= mean(SMAD4 == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_smad, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with SMAD4 Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with SMAD4 Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with SMAD4 Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_smad > 0.05, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    geom_point() + geom_text(aes(label=ifelse(mean_smad > 0.45, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    annotate("text", x = 1, y = 0.49, label = "p-value = 0.013573") +
    annotate("text", x = 1, y = 0.47, label = "q-value = 0.049144") +
    theme(legend.position = "none") 
  
}
mss_smad_graph(data_mss_only)



mss_pik_graph <- function(df){
  
  jitter <- df %>%
    group_by(Neighborhood, cohort) %>% 
    summarize(
      mean_pik= mean(PIK3CA == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = cohort, y = mean_pik, fill = cohort)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("MSS" = "#8b1ff0", "Non-MSS" = "#009E73")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with PIK3CA Mutation Colorectal\nCancer and Primary Tumor Location of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with PIK3CA Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Microsatellite Stability of Patients with Colorectal Cancer at MSK",
         y = "Percentage of Patients with PIK3CA Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_smad > 0.05, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    #geom_point() + geom_text(aes(label=ifelse(mean_smad > 0.45, as.character(Neighborhood),'')), hjust = 0,vjust = 0, size = 2) +
    annotate("text", x = 1, y = 0.52, label = "p-value < 1e-07") +
    annotate("text", x = 1, y = 0.47, label = "q-value < 1e-07") +
    theme(legend.position = "none") 
  
}
mss_pik_graph(data_mss_only)



bmi_pik_graph <- function(df){
  labelss <- c("No PIK3CA Mutation", "PIK3CA Mutation")
  df$PIK3CA <- as.character(df$PIK3CA)
  bmi_pik_vector <- df %>%
    group_by(Neighborhood, PIK3CA) %>%
    summarise(median_bmi = median(First_BMI_Measurement, na.rm = TRUE))
  
  ggplot(bmi_pik_vector, aes(PIK3CA, median_bmi, fill = PIK3CA)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("1" = "#8b1ff0", "0" = "#009E73")) +
    labs(title = "Median BMI at Diagnosis by Neighborhood, grouped by PIK3CA Mutation\nClassification of Patients with Colorectal Cancer from the MSK-Impact Data\nSet",
         x = "Percentage of Patients with PIK3CA Mutation",
         y = "Median BMI at Diagnosis") +
    scale_x_discrete(labels = labelss) +
    theme_minimal() +
    geom_point() + geom_text(aes(label=ifelse(median_bmi>30,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #geom_point() + geom_text(aes(label=ifelse(median_age>65,as.character(Neighborhood),'')), hjust = 0,vjust = 0) +
    #tp53 mut highest point : union square - lower east side
    annotate("text", x = 1, y = 22, label = "p-value = 0.003633") + 
    annotate("text", x = 1, y = 21.5, label = "q-value = 0.015894") +
    theme(legend.position = "none")
}
bmi_pik_graph(data_mss_only)


ancest_kras_graph <- function(df){
  
  jitter <- data20 %>%
    group_by(Neighborhood, Ancestry_Label) %>% 
    summarize(
      mean_kras = mean(KRAS == "1"),
      .groups = 'drop'
    )
  
  ggplot(jitter, aes(x = Ancestry_Label, y = mean_kras, fill = Ancestry_Label)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("ADM" = "#E69F00", "AFR" = "#56B4E9", "ASJ" = "#009E73", "EAS" = "#8b1ff0", "EUR" = "#e33232", "NAM" = "#271f9c", "SAS" = "#e956e2")) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(title = "Association Between Percentage of Patients with KRAS Mutation Colorectal\nCancer and Ancestry Label of Patients with Colorectal Cancer,\nLayered Underneath the Percentage of Patients with KRAS Mutation per\nPrimary Tumor Location and Neighborhood",
         x = "Primary Tumor Location of Patients at MSK",
         y = "Percentage of Patients with KRAS Mutation") +
    #geom_point() + geom_text(aes(label=ifelse(mean_kras == 1, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 4) +
    #geom_point() + geom_text(aes(label=ifelse(mean_kras == 0, as.character(Neighborhood),'')), hjust = 0,vjust = -0.4, size = 2) +
    #annotate("text", x = 1, y = 0.9, label = "p-value = 0.000114") +
    #annotate("text", x = 1, y = 0.85, label = "q-value = 0.000665") +
    theme(legend.position = "none") 
  
}
ancest_kras_graph(data_mss_only)
