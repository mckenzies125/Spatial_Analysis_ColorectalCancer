#Script Summary:
# 1. Compile a matrix in which each cell corresponds to a p value resulted from testing its respective variables against neighborhood
# 2. Compile a matrix in which each cell corresponds to a q value resulted from testing its respective variables against neighborhood


#THIS is the multivariate analysis that Dr. Sharafudeen Abubakar and I believed could also support our analysis


#In this script, we test the independence of each variable and neighborhood, and stratify by
#a second variable, treating it as a confounder. We use coin's independence test for the following reasons:
# 1. all variable must be categorical if we want to stratify by then---meaning we no longer have
#    use for Kruskal-Wallis.
# 2. by design, chi^2/fisher's exact tests cannot include stratification 
# 3. another option is Cochran-Mantel-Haenszel test, but one of the assumptions of this test is that
#    the odds ratio across all strata is the same---our data almost definitely violates this assumption
# Therefore, we opt for a customizable independence test like the one in the coin library

#coin library documentation: https://cran.r-project.org/web/packages/coin/vignettes/coin.pdf


#This script is written so that it doesn't HAVE to be run while on the cluster. However, for 
#a significantly faster outcome, run on the cluster


#To see the values without running all the code, go to Calculate_p_q_vals -> Multivariate -> var_vs_neigh_str_by_var


#libraries
library(coin)
library(dplyr)

data20 <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")
names(data20)[11] <- "Smoking_History"
names(data20)[13] <- "Yost_Index"
names(data20)[21] <- "Prior_Treatment"

data20_mss_only <- read.csv("~\\clinical_and_genomic\\genomic_nyc_mss_colorec_20.csv")
names(data20_mss_only)[11] <- "Smoking_History"
names(data20_mss_only)[13] <- "Yost_Index"
names(data20_mss_only)[21] <- "Prior_Treatment"


#initialize pvals matrix
final_pvals <- data.frame(variable_1 = c("Yost_Index_Group", "Age_At_Diagnosis_Group", "Primary_Tumor_Location_Final", "Sex", "Cohort", "Ancestry_Label", "BMI_Groups", "Prior_Treatment", "Smoking_History", "KRAS", "APC", "TP53", "BRAF", "SMAD4", "PIK3CA"), 
                          Yost_Index_Group = '', 
                          Age_At_Diagnosis_Group = '',
                          Primary_Tumor_Location_Final = '',
                          Sex = '',
                          Cohort = '',
                          Ancestry_Label = '',
                          BMI_Groups = '',
                          Prior_Treatment = '',
                          Smoking_History = '',
                          KRAS = '',
                          APC = '',
                          TP53 = '',
                          BRAF = '',
                          SMAD4 = '',
                          PIK3CA = '')

#make yost index into a categorical variable---based on categories provided by Dr. Sharafudeen Abubakar
yost_grouping_data20 <- function(df){
  df["Yost_Index_Group"] = ''
  for (i in 1:nrow(df)){
    if (df[i, "Yost_Index"] <= 20.5){
      df[i, "Yost_Index_Group"] = "High_affluence"
    }else if (df[i, "Yost_Index"] > 20.5 && df[i, "Yost_Index"] <= 49){
      df[i, "Yost_Index_Group"] = "Middle_affluence"
    }else{
      df[i, "Yost_Index_Group"] = "Low_affluence"
    }
  }
  data20 <<- df
}
yost_grouping_data20(data20)

#make yost index into a categorical variable---based on categories provided by Dr. Sharafudeen Abubakar
yost_grouping_data20_mss_only <- function(df){
  df["Yost_Index_Group"] = ''
  for (i in 1:nrow(df)){
    if (df[i, "Yost_Index"] <= 20.5){
      df[i, "Yost_Index_Group"] = "High_affluence"
    }else if (df[i, "Yost_Index"] > 20.5 && df[i, "Yost_Index"] <= 49){
      df[i, "Yost_Index_Group"] = "Middle_affluence"
    }else{
      df[i, "Yost_Index_Group"] = "Low_affluence"
    }
  }
  data20_mss_only <<- df
}
yost_grouping_data20_mss_only(data20_mss_only)

#calculate pvals
neigh_vs_var_control_var <- function(df, dfmss){
  set.seed(42)
  temp <- data20
  temp_mss_only <- data20_mss_only
  
  temp$Neighborhood <- as.factor(temp$Neighborhood)
  temp$Sex <- as.factor(temp$Sex)
  temp$Prior_Treatment <- as.factor(temp$Prior_Treatment)
  temp$Smoking_History <- as.factor(temp$Smoking_History)
  temp$Primary_Tumor_Location_Final <- as.factor(temp$Primary_Tumor_Location_Final)
  temp$Ancestry_Label <- as.factor(temp$Ancestry_Label)
  temp$Yost_Index_Group <- as.factor(temp$Yost_Index_Group)
  temp$age_dx_final_group <- as.factor(temp$age_dx_final_group)
  temp$BMI_Groups <- as.factor(temp$BMI_Groups)
  temp$cohort <- as.factor(temp$cohort)
  
  temp_mss_only$KRAS <- as.factor(temp_mss_only$KRAS)
  temp_mss_only$APC <- as.factor(temp_mss_only$APC)
  temp_mss_only$TP53 <- as.factor(temp_mss_only$TP53)
  temp_mss_only$BRAF <- as.factor(temp_mss_only$BRAF)
  temp_mss_only$SMAD4 <- as.factor(temp_mss_only$SMAD4)
  temp_mss_only$PIK3CA <- as.factor(temp_mss_only$PIK3CA)
  temp_mss_only$Neighborhood <- as.factor(temp_mss_only$Neighborhood)
  temp_mss_only$Sex <- as.factor(temp_mss_only$Sex)
  temp_mss_only$Prior_Treatment <- as.factor(temp_mss_only$Prior_Treatment)
  temp_mss_only$Smoking_History <- as.factor(temp_mss_only$Smoking_History)
  temp_mss_only$Primary_Tumor_Location_Final <- as.factor(temp_mss_only$Primary_Tumor_Location_Final)
  temp_mss_only$Ancestry_Label <- as.factor(temp_mss_only$Ancestry_Label)
  temp_mss_only$Yost_Index_Group <- as.factor(temp_mss_only$Yost_Index_Group)
  temp_mss_only$age_dx_final_group <- as.factor(temp_mss_only$age_dx_final_group)
  temp_mss_only$BMI_Groups <- as.factor(temp_mss_only$BMI_Groups)
  temp_mss_only$cohort <- as.factor(temp_mss_only$cohort)
  
  
  age_yost <- independence_test(
    age_dx_final_group ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,2] <<- pvalue(age_yost)
  rm(age_yost)
  
  yost_age <- independence_test(
    Yost_Index_Group ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,3] <<- pvalue(yost_age)
  rm(yost_age)
  
  yost_bmi <- independence_test(
    Yost_Index_Group ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,8] <<- pvalue(yost_bmi)
  rm(yost_bmi)
  gc()
  
  bmi_yost <- independence_test(
    BMI_Groups ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,2] <<- pvalue(bmi_yost)
  rm(bmi_yost)
  
  age_bmi <- independence_test(
    age_dx_final_group ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,8] <<- pvalue(age_bmi)
  rm(age_bmi)
  
  bmi_age <- independence_test(
    BMI_Groups ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,3] <<- pvalue(bmi_age)
  rm(bmi_age)
  gc()
  
  
  
  
  #yost vs location
  yost_loc <- independence_test(
    Yost_Index_Group ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,4] <<- pvalue(yost_loc)
  rm(yost_loc)
  
  loc_yost <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,2] <<- pvalue(loc_yost)
  rm(loc_yost)
  
  #yost vs sex
  yost_sex <- independence_test(
    Yost_Index_Group ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,5] <<- pvalue(yost_sex)
  rm(yost_sex)
  gc()
  
  sex_yost <- independence_test(
    Sex ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,2] <<- pvalue(sex_yost)
  rm(sex_yost)
  
  #yost vs ancestry
  yost_ancest <- independence_test(
    Yost_Index_Group ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,7] <<- pvalue(yost_ancest)
  rm(yost_ancest)
  
  ancest_yost <- independence_test(
    Ancestry_Label ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,2] <<- pvalue(ancest_yost)
  rm(ancest_yost)
  gc()
  
  #yost vs prior treatment
  yost_prior <- independence_test(
    Yost_Index_Group ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,9] <<- pvalue(yost_prior)
  rm(yost_prior)
  
  prior_yost <- independence_test(
    Prior_Treatment ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,2] <<- pvalue(prior_yost)
  rm(prior_yost)
  
  #yost vs smoking
  yost_smoking <- independence_test(
    Yost_Index_Group ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,10] <<- pvalue(yost_smoking)
  rm(yost_smoking)
  gc()
  
  smoke_yost <- independence_test(
    Smoking_History ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,2] <<- pvalue(smoke_yost)
  rm(smoke_yost)
  
  
  
  
  #yost vs kras
  yost_kras <- independence_test(
    Yost_Index_Group ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,11] <<- pvalue(yost_kras)
  rm(yost_kras)
  
  kras_yost <- independence_test(
    KRAS ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,2] <<- pvalue(kras_yost)
  rm(kras_yost)
  gc()
  
  #yost vs apc
  yost_apc <- independence_test(
    Yost_Index_Group ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,12] <<- pvalue(yost_apc)
  rm(yost_apc)
  
  apc_yost <- independence_test(
    APC ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,2] <<- pvalue(apc_yost)
  rm(apc_yost)
  
  #yost vs tp53
  yost_tp <- independence_test(
    Yost_Index_Group ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,13] <<- pvalue(yost_tp)
  rm(yost_tp)
  gc()
  
  tp_yost <- independence_test(
    TP53 ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,2] <<- pvalue(tp_yost)
  rm(tp_yost)
  
  #yost vs braf
  yost_braf <- independence_test(
    Yost_Index_Group ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,14] <<- pvalue(yost_braf)
  rm(yost_braf)
  
  braf_yost <- independence_test(
    BRAF ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,2] <<- pvalue(braf_yost)
  rm(braf_yost)
  gc()
  
  #yost vs smad4
  yost_smad4 <- independence_test(
    Yost_Index_Group ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,15] <<- pvalue(yost_smad4)
  rm(yost_smad4)
  
  smad4_yost <- independence_test(
    SMAD4 ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,2] <<- pvalue(smad4_yost)
  rm(smad4_yost)
  
  #yost vs pik3ca
  yost_pik <- independence_test(
    Yost_Index_Group ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,16] <<- pvalue(yost_pik)
  rm(yost_pik)
  gc()
  
  pik_yost <- independence_test(
    PIK3CA ~ Neighborhood | Yost_Index_Group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,2] <<- pvalue(pik_yost)
  rm(pik_yost)
  
  #age vs location
  age_loc <- independence_test(
    age_dx_final_group ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,4] <<- pvalue(age_loc)
  rm(age_loc)
  
  loc_age <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,3] <<- pvalue(loc_age)
  rm(loc_age)
  gc()
  
  #age vs sex
  age_sex <- independence_test(
    age_dx_final_group ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,5] <<- pvalue(age_sex)
  rm(age_sex)
  
  
  sex_age <- independence_test(
    Sex ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,3] <<- pvalue(sex_age)
  rm(sex_age)
  
  #age vs ancestry
  age_ancest <- independence_test(
    age_dx_final_group ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,7] <<- pvalue(age_ancest)
  rm(age_ancest)
  gc()
  
  ancest_age <- independence_test(
    Ancestry_Label ~ Neighborhood | age_dx_final_group, 
    data = temp,
    #scores = list("Ancestry_Label" = c("ADM" = 1, "AFR" = 2, "ASJ" = 3, "EAS" = 4, "EUR" = 5, "NAM" = 6, "SAS" = 7)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,3] <<- pvalue(ancest_age)
  rm(ancest_age)
  
  #age vs prior treatment
  age_prior <- independence_test(
    age_dx_final_group ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,9] <<- pvalue(age_prior)
  rm(age_prior)
  
  
  prior_age <- independence_test(
    Prior_Treatment ~ Neighborhood | age_dx_final_group, 
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,3] <<- pvalue(prior_age)
  rm(prior_age)
  gc()
  
  #age vs smoking
  age_smoking <- independence_test(
    age_dx_final_group ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,10] <<- pvalue(age_smoking)
  rm(age_smoking)
  
  
  smoking_age <- independence_test(
    Smoking_History ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,3] <<- pvalue(smoking_age)
  rm(smoking_age)
  
  
  #age vs kras
  age_kras <- independence_test(
    age_dx_final_group ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,11] <<- pvalue(age_kras)
  rm(age_kras)
  gc()
  
  kras_age <- independence_test(
    KRAS ~ Neighborhood | age_dx_final_group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,3] <<- pvalue(kras_age)
  rm(kras_age)
  
  #age vs apc
  age_apc <- independence_test(
    age_dx_final_group ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,12] <<- pvalue(age_apc)
  rm(age_apc)
  
  apc_age <- independence_test(
    APC ~ Neighborhood | age_dx_final_group, 
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,3] <<- pvalue(apc_age)
  rm(apc_age)
  gc()
  
  #age vs tp53
  age_tp <- independence_test(
    age_dx_final_group ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,13] <<- pvalue(age_tp)
  rm(age_tp)
  
  tp_age <- independence_test(
    TP53 ~ Neighborhood | age_dx_final_group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,3] <<- pvalue(tp_age)
  rm(tp_age)
  
  #age vs braf
  age_braf <- independence_test(
    age_dx_final_group ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,14] <<- pvalue(age_braf)
  rm(age_braf)
  gc()
  
  braf_age <- independence_test(
    BRAF ~ Neighborhood | age_dx_final_group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,3] <<- pvalue(braf_age)
  rm(braf_age)
  
  #age vs smad4
  age_smad <- independence_test(
    age_dx_final_group ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,15] <<- pvalue(age_smad)
  rm(age_smad)
  
  smad_age <- independence_test(
    SMAD4 ~ Neighborhood | age_dx_final_group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,3] <<- pvalue(smad_age)
  rm(smad_age)
  gc()
  
  #age vs pik3ca
  age_pik <- independence_test(
    age_dx_final_group ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,16] <<- pvalue(age_pik)
  rm(age_pik)
  
  pik_age <- independence_test(
    PIK3CA ~ Neighborhood | age_dx_final_group,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,3] <<- pvalue(pik_age)
  rm(pik_age)
  
  #bmi vs location
  bmi_loc <- independence_test(
    BMI_Groups ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,4] <<- pvalue(bmi_loc)
  rm(bmi_loc)
  gc()
  
  
  loc_bmi <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | BMI_Groups,
    data = temp,
    #scores = list("Primary_Tumor_Location_Final" = c("Right" = 1, "Left" = 2, "Rectum" = 3)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,8] <<- pvalue(loc_bmi)
  rm(loc_bmi)
  
  
  
  
  
  
  #bmi vs sex
  bmi_sex <- independence_test(
    BMI_Groups ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,5] <<- pvalue(bmi_sex)
  rm(bmi_sex)
  
  
  sex_bmi <- independence_test(
    Sex ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,8] <<- pvalue(sex_bmi)
  rm(sex_bmi)
  gc()
  
  #bmi vs ancestry
  bmi_ancest <- independence_test(
    BMI_Groups ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,7] <<- pvalue(bmi_ancest)
  rm(bmi_ancest)
  
  ancest_bmi <- independence_test(
    Ancestry_Label ~ Neighborhood | BMI_Groups,
    data = temp,
    #scores = list("Ancestry_Label" = c("ADM" = 1, "AFR" = 2, "ASJ" = 3, "EAS" = 4, "EUR" = 5, "NAM" = 6, "SAS" = 7)),
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,8] <<- pvalue(ancest_bmi)
  rm(ancest_bmi)
  
  
  #bmi vs prior treatment
  bmi_prior <- independence_test(
    BMI_Groups ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,9] <<- pvalue(bmi_prior)
  rm(bmi_prior)
  gc()
  
  prior_bmi <- independence_test(
    Prior_Treatment ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,8] <<- pvalue(prior_bmi)
  rm(prior_bmi)
  
  
  #bmi vs smoking
  bmi_smoking <- independence_test(
    BMI_Groups ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,10] <<- pvalue(bmi_smoking)
  rm(bmi_smoking)
  
  
  smoking_bmi <- independence_test(
    Smoking_History ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,8] <<- pvalue(smoking_bmi)
  rm(smoking_bmi)
  gc()
  
  
  
  #bmi vs kras
  bmi_kras <- independence_test(
    BMI_Groups ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,11] <<- pvalue(bmi_kras)
  rm(bmi_kras)
  
  kras_bmi <- independence_test(
    KRAS ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,8] <<- pvalue(kras_bmi)
  rm(kras_bmi)
  
  #bmi vs apc
  bmi_apc <- independence_test(
    BMI_Groups ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,12] <<- pvalue(bmi_apc)
  rm(bmi_apc)
  gc()
  
  apc_bmi <- independence_test(
    APC ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,8] <<- pvalue(apc_bmi)
  rm(apc_bmi)
  
  #bmi vs tp53
  bmi_tp <- independence_test(
    BMI_Groups ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,13] <<- pvalue(bmi_tp)
  rm(bmi_tp)
  
  tp_bmi <- independence_test(
    TP53 ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,8] <<- pvalue(tp_bmi)
  rm(tp_bmi)
  gc()
  
  #bmi vs braf
  bmi_braf <- independence_test(
    BMI_Groups ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,14] <<- pvalue(bmi_braf)
  rm(bmi_braf)
  
  braf_bmi <- independence_test(
    BRAF ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,8] <<- pvalue(braf_bmi)
  rm(braf_bmi)
  
  #bmi vs smad4
  bmi_smad <- independence_test(
    BMI_Groups ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,15] <<- pvalue(bmi_smad)
  rm(bmi_smad)
  gc()
  
  smad_bmi <- independence_test(
    SMAD4 ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,8] <<- pvalue(smad_bmi)
  rm(smad_bmi)
  
  #bmi vs pik3ca
  bmi_pik <- independence_test(
    BMI_Groups ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,16] <<- pvalue(bmi_pik)
  rm(bmi_pik)
  
  pik_bmi <- independence_test(
    PIK3CA ~ Neighborhood | BMI_Groups,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,8] <<- pvalue(pik_bmi)  
  rm(pik_bmi)
  gc()
  
  
  
  
  results_loc_sex <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,5] <<- pvalue(results_loc_sex)
  rm(results_loc_sex)
  
  results_sex_loc <- independence_test(
    Sex ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,4] <<- pvalue(results_sex_loc)
  rm(results_sex_loc)
  
  
  #location vs ancestry
  results_loc_ancest <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,7] <<- pvalue(results_loc_ancest)
  rm(results_loc_ancest)
  gc()
  
  results_ancest_loc <- independence_test(
    Ancestry_Label ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,4] <<- pvalue(results_ancest_loc)
  rm(results_ancest_loc)
  
  #location vs prior treatment
  results_loc_priortreat <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,9] <<- pvalue(results_loc_priortreat)
  rm(results_loc_priortreat)
  
  results_prior_loc <- independence_test(
    Prior_Treatment ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,4] <<- pvalue(results_prior_loc)
  rm(results_prior_loc)
  gc()
  
  #location vs smoking
  results_loc_smoke <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,10] <<- pvalue(results_loc_smoke)
  rm(results_loc_smoke)
  
  results_smoke_loc <- independence_test(
    Smoking_History ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,4] <<- pvalue(results_smoke_loc)
  rm(results_smoke_loc)
  
  #location vs kras
  results_loc_kras <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,11] <<- pvalue(results_loc_kras)
  rm(results_loc_kras)
  gc()
  
  results_kras_loc <- independence_test(
    KRAS ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,4] <<- pvalue(results_kras_loc)
  rm(results_kras_loc)
  
  #location vs apc
  results_loc_apc <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,12] <<- pvalue(results_loc_apc)
  rm(results_loc_apc)
  
  results_apc_loc <- independence_test(
    APC ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,4] <<- pvalue(results_apc_loc)
  rm(results_apc_loc)
  gc()
  
  #location vs tp53
  results_loc_tp <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,13] <<- pvalue(results_loc_tp)
  rm(results_loc_tp)
  
  results_tp_loc <- independence_test(
    TP53 ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,4] <<- pvalue(results_tp_loc)
  rm(results_tp_loc)
  
  #location vs braf
  results_loc_braf <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,14] <<- pvalue(results_loc_braf)
  rm(results_loc_braf)
  gc()
  
  results_braf_loc <- independence_test(
    BRAF ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,4] <<- pvalue(results_braf_loc)
  rm(results_braf_loc)
  
  #location vs smad4
  results_loc_smad <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,15] <<- pvalue(results_loc_smad)
  rm(results_loc_smad)
  
  results_smad_loc <- independence_test(
    SMAD4 ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,4] <<- pvalue(results_smad_loc)
  rm(results_smad_loc)
  gc()
  
  #location vs pik3ca
  results_loc_pik <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,16] <<- pvalue(results_loc_pik)
  rm(results_loc_pik)
  
  results_pik_loc <- independence_test(
    PIK3CA ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,4] <<- pvalue(results_pik_loc)
  rm(results_pik_loc)
  
  
  #Sex vs Ancestry
  results_sex_ancest <- independence_test(
    Sex ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,7] <<- pvalue(results_sex_ancest)
  rm(results_sex_ancest)
  gc()
  
  results_ancest_sex <- independence_test(
    Ancestry_Label ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,5] <<- pvalue(results_ancest_sex)
  rm(results_ancest_sex)
  
  
  #sex vs prior treatment
  results_sex_prior <- independence_test(
    Sex ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,9] <<- pvalue(results_sex_prior)
  rm(results_sex_prior)
  
  
  results_prior_sex <- independence_test(
    Prior_Treatment ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,5] <<- pvalue(results_prior_sex)
  rm(results_prior_sex)
  gc()
  
  #sex vs smoking
  results_sex_smoke <- independence_test(
    Sex ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,10] <<- pvalue(results_sex_smoke)
  rm(results_sex_smoke)
  
  results_smoke_sex <- independence_test(
    Smoking_History ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,5] <<- pvalue(results_smoke_sex)
  rm(results_smoke_sex)
  
  #sex vs kras
  results_sex_kras <- independence_test(
    Sex ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,11] <<- pvalue(results_sex_kras)
  rm(results_sex_kras)
  gc()
  
  results_kras_sex <- independence_test(
    KRAS ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,5] <<- pvalue(results_kras_sex)
  rm(results_kras_sex)
  
  #sex vs apc
  results_sex_apc <- independence_test(
    Sex ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,12] <<- pvalue(results_sex_apc)
  rm(results_sex_apc)
  
  results_apc_sex <- independence_test(
    APC ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,5] <<- pvalue(results_apc_sex)
  rm(results_apc_sex)
  gc()
  
  #sex vs tp53
  results_sex_tp <- independence_test(
    Sex ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,13] <<- pvalue(results_sex_tp)
  rm(results_sex_tp)
  
  results_tp_sex <- independence_test(
    TP53 ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,5] <<- pvalue(results_tp_sex)
  rm(results_tp_sex)
  
  #sex vs braf
  results_sex_braf <- independence_test(
    Sex ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,14] <<- pvalue(results_sex_braf)
  rm(results_sex_braf)
  gc()
  
  results_braf_sex <- independence_test(
    BRAF ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,5] <<- pvalue(results_braf_sex)
  rm(results_braf_sex)
  
  
  #sex vs smad4
  results_sex_smad <- independence_test(
    Sex ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,15] <<- pvalue(results_sex_smad)
  rm(results_sex_smad)
  
  results_smad_sex <- independence_test(
    SMAD4 ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,5] <<- pvalue(results_smad_sex)
  rm(results_smad_sex)
  gc()
  
  #sex vs pik
  results_sex_pik <- independence_test(
    Sex ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,16] <<- pvalue(results_sex_pik)
  rm(results_sex_pik)
  
  results_pik_sex <- independence_test(
    PIK3CA ~ Neighborhood | Sex,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,5] <<- pvalue(results_pik_sex)
  rm(results_pik_sex)
  
  #prior treatment vs smoking
  results_prior_smoke <- independence_test(
    Prior_Treatment ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,10] <<- pvalue(results_prior_smoke)
  rm(results_prior_smoke)
  gc()
  
  results_smoke_prior <- independence_test(
    Smoking_History ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,9] <<- pvalue(results_smoke_prior)
  rm(results_smoke_prior)
  
  #prior treatment vs kras
  results_prior_kras <- independence_test(
    Prior_Treatment ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,11] <<- pvalue(results_prior_kras)
  rm(results_prior_kras)
  
  results_kras_prior <- independence_test(
    KRAS ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,9] <<- pvalue(results_kras_prior)
  rm(results_kras_prior)
  gc()
  
  #prior treatment vs apc
  results_prior_apc <- independence_test(
    Prior_Treatment ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,12] <<- pvalue(results_prior_apc)
  rm(results_prior_apc)
  
  results_apc_prior <- independence_test(
    APC ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,9] <<- pvalue(results_apc_prior)
  rm(results_apc_prior)
  
  #prior treatment vs tp53
  results_prior_tp <- independence_test(
    Prior_Treatment ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,13] <<- pvalue(results_prior_tp)
  rm(results_prior_tp)
  gc()
  
  results_tp_prior <- independence_test(
    TP53 ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,9] <<- pvalue(results_tp_prior)
  rm(results_tp_prior)
  
  #prior treatment vs braf
  results_prior_braf <- independence_test(
    Prior_Treatment ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,14] <<- pvalue(results_prior_braf)
  rm(results_prior_braf)
  
  results_braf_prior <- independence_test(
    BRAF ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,9] <<- pvalue(results_braf_prior)
  rm(results_braf_prior)
  gc()
  
  #prior treatment vs smad4
  results_prior_smad <- independence_test(
    Prior_Treatment ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,15] <<- pvalue(results_prior_smad)
  rm(results_prior_smad)
  
  results_smad_prior <- independence_test(
    SMAD4 ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,9] <<- pvalue(results_smad_prior)
  rm(results_smad_prior)
  
  #prior treatment vs PIK3CA
  results_prior_pik <- independence_test(
    Prior_Treatment ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,16] <<- pvalue(results_prior_pik)
  rm(results_prior_pik)
  gc()
  
  results_pik_prior <- independence_test(
    PIK3CA ~ Neighborhood | Prior_Treatment,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,9] <<- pvalue(results_pik_prior)
  rm(results_pik_prior)
  
  #smoking history vs kras
  results_smoke_kras <- independence_test(
    Smoking_History ~ KRAS | Neighborhood,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,11] <<- pvalue(results_smoke_kras)
  rm(results_smoke_kras)
  
  results_kras_smoke <- independence_test(
    KRAS ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,10] <<- pvalue(results_kras_smoke)
  rm(results_kras_smoke)
  gc()
  
  #smoking history vs apc
  results_smoke_apc <- independence_test(
    Smoking_History ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,12] <<- pvalue(results_smoke_apc)
  rm(results_smoke_apc)
  
  results_apc_smoke <- independence_test(
    APC ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,10] <<- pvalue(results_apc_smoke)
  rm(results_apc_smoke)
  
  #smoking history vs tp53
  results_smoke_tp <- independence_test(
    Smoking_History ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,13] <<- pvalue(results_smoke_tp)
  rm(results_smoke_tp)
  gc()
  
  results_tp_smoke <- independence_test(
    TP53 ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,10] <<- pvalue(results_tp_smoke)
  rm(results_tp_smoke)
  
  #smoking history vs braf
  results_smoke_braf <- independence_test(
    Smoking_History ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,14] <<- pvalue(results_smoke_braf)
  rm(results_smoke_braf)
  
  results_braf_smoke <- independence_test(
    BRAF ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,10] <<- pvalue(results_braf_smoke)
  rm(results_braf_smoke)
  gc()
  
  #smoking history vs smad4
  results_smoke_smad <- independence_test(
    Smoking_History ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,15] <<- pvalue(results_smoke_smad)
  rm(results_smoke_smad)
  
  results_smad_smoke <- independence_test(
    SMAD4 ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,10] <<- pvalue(results_smad_smoke)
  rm(results_smad_smoke)
  
  #smoking history vs pik
  results_smoke_pik <- independence_test(
    Smoking_History ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,16] <<- pvalue(results_smoke_pik)
  rm(results_smoke_pik)
  gc()
  
  
  #kras vs apc
  results_kras_apc <- independence_test(
    KRAS ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,12] <<- pvalue(results_kras_apc)
  rm(results_kras_apc)
  
  results_apc_kras <- independence_test(
    APC ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,11] <<- pvalue(results_apc_kras)
  rm(results_apc_kras)
  
  #kras vs tp53
  results_kras_tp <- independence_test(
    KRAS ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,13] <<- pvalue(results_kras_tp)
  rm(results_kras_tp)
  gc()
  
  
  results_tp_kras <- independence_test(
    TP53 ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,11] <<- pvalue(results_tp_kras)
  rm(results_tp_kras)
  
  #kras vs BRAF
  results_kras_braf <- independence_test(
    KRAS ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,14] <<- pvalue(results_kras_braf)
  rm(results_kras_braf)
  
  results_braf_kras <- independence_test(
    BRAF ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,11] <<- pvalue(results_braf_kras)
  rm(results_braf_kras)
  gc()
  
  #kras vs smad4
  results_kras_smad <- independence_test(
    KRAS ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,15] <<- pvalue(results_kras_smad)
  rm(results_kras_smad)
  
  results_smad_kras <- independence_test(
    SMAD4 ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,11] <<- pvalue(results_smad_kras)
  rm(results_smad_kras)
  
  
  #kras vs pik3c4
  results_kras_pik <- independence_test(
    KRAS ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,16] <<- pvalue(results_kras_pik)
  rm(results_kras_pik)
  gc()
  
  results_pik_kras <- independence_test(
    PIK3CA ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,11] <<- pvalue(results_pik_kras)
  rm(results_pik_kras)
  
  #apc vs tp53
  results_apc_tp <- independence_test(
    APC ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,13] <<- pvalue(results_apc_tp)
  rm(results_apc_tp)
  
  results_tp_apc <- independence_test(
    TP53 ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,12] <<- pvalue(results_tp_apc)
  rm(results_tp_apc)
  gc()
  
  #apc vs braf
  results_apc_braf <- independence_test(
    APC ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,14] <<- pvalue(results_apc_braf)
  rm(results_apc_braf)
  
  results_braf_apc <- independence_test(
    BRAF ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,12] <<- pvalue(results_braf_apc)
  rm(results_braf_apc)
  
  #apc vs smad4
  results_apc_smad <- independence_test(
    APC ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,15] <<- pvalue(results_apc_smad)
  rm(results_apc_smad)
  gc()
  
  results_smad_apc <- independence_test(
    SMAD4 ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,12] <<- pvalue(results_smad_apc)
  rm(results_smad_apc)
  
  
  #apc vs pik3ca
  results_apc_pik <- independence_test(
    APC ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,16] <<- pvalue(results_apc_pik)
  rm(results_apc_pik)
  
  results_pik_apc <- independence_test(
    APC ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,12] <<- pvalue(results_pik_apc)
  rm(results_pik_apc)
  gc()
  
  #tp53 vs braf
  results_tp_braf <- independence_test(
    TP53 ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,14] <<- pvalue(results_tp_braf)
  rm(results_tp_braf)
  
  results_braf_tp <- independence_test(
    BRAF ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,13] <<- pvalue(results_braf_tp)
  rm(results_braf_tp)
  
  #tp53 vs smad4
  results_tp_smad <- independence_test(
    TP53 ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,15] <<- pvalue(results_tp_smad)
  rm(results_tp_smad)
  gc()
  
  results_smad_tp <- independence_test(
    SMAD4 ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,13] <<- pvalue(results_smad_tp)
  rm(results_smad_tp)
  
  
  #tp53 vs pik3ca
  results_tp_pik <- independence_test(
    TP53 ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,16] <<- pvalue(results_tp_pik)
  rm(results_tp_pik)
  
  results_pik_tp <- independence_test(
    PIK3CA ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,13] <<- pvalue(results_pik_tp)
  rm(results_pik_tp)
  gc()
  
  #barf vs smad4
  results_braf_smad <- independence_test(
    BRAF ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,15] <<- pvalue(results_braf_smad)
  rm(results_braf_smad)
  
  results_smad_braf <- independence_test(
    SMAD4 ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,14] <<- pvalue(results_smad_braf)
  rm(results_smad_braf)
  
  #barf vs pik3ca
  results_braf_pik <- independence_test(
    BRAF ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,16] <<- pvalue(results_braf_pik)
  rm(results_braf_pik)
  gc()
  
  
  results_pik_braf <- independence_test(
    PIK3CA ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,14] <<- pvalue(results_pik_braf)
  rm(results_pik_braf)
  
  #smad4 vs pik3ca
  results_pik_smad <- independence_test(
    PIK3CA ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,15] <<- pvalue(results_pik_smad)
  rm(results_pik_smad)
  
  results_smad_pik <- independence_test(
    SMAD4 ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,16] <<- pvalue(results_smad_pik)
  rm(results_smad_pik)
  gc()
  
  
  #mss vs yost
  yost_mss <- independence_test(
    Yost_Index_Group ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[1,6] <<- pvalue(yost_mss)
  rm(yost_mss)
  
  mss_yost <- independence_test(
    cohort ~ Neighborhood | Yost_Index_Group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,2] <<- pvalue(mss_yost)
  rm(mss_yost)
  
  #mss vs age
  age_mss <- independence_test(
    age_dx_final_group ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[2,6] <<- pvalue(age_mss)
  rm(age_mss)
  gc()
  
  mss_age <- independence_test(
    cohort ~ Neighborhood | age_dx_final_group,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,3] <<- pvalue(mss_age)
  rm(mss_age)
  
  
  #mss vs bmi
  bmi_mss <- independence_test(
    BMI_Groups ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[7,6] <<- pvalue(bmi_mss)
  rm(bmi_mss)
  
  mss_bmi <- independence_test(
    cohort ~ Neighborhood | BMI_Groups,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,8] <<- pvalue(mss_bmi)
  rm(mss_bmi)
  gc()
  
  
  #mss vs location
  cont_loc_mss <- independence_test(
    Primary_Tumor_Location_Final ~ Neighborhood | cohort ,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[3,6] <<- pvalue(cont_loc_mss)
  rm(cont_loc_mss)
  
  cont_mss_loc <- independence_test(
    cohort ~ Neighborhood | Primary_Tumor_Location_Final,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,4] <<- pvalue(cont_mss_loc)
  rm(cont_mss_loc)
  
  #mss vs sex
  cont_sex_mss <- independence_test(
    Sex ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[4,6] <<- pvalue(cont_sex_mss)
  rm(cont_sex_mss)
  gc()
  
  cont_mss_sex <- independence_test(
    cohort ~ Neighborhood | Sex,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,5] <<- pvalue(cont_mss_sex)
  rm(cont_mss_sex)
  
  #mss vs ancestry
  cont_ancest_mss <- independence_test(
    Ancestry_Label ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,6] <<- pvalue(cont_ancest_mss)
  rm(cont_ancest_mss)
  
  cont_mss_ancest <- independence_test(
    cohort ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,7] <<- pvalue(cont_mss_ancest)
  rm(cont_mss_ancest)
  gc()
  
  #mss vs prior treatment
  cont_prior_mss <- independence_test(
    Prior_Treatment ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,6] <<- pvalue(cont_prior_mss)
  rm(cont_prior_mss)
  
  cont_mss_prior <- independence_test(
    cohort ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,9] <<- pvalue(cont_mss_prior)
  rm(cont_mss_prior)
  
  #mss vs smoking
  cont_smoke_mss <- independence_test(
    Smoking_History ~ Neighborhood | cohort,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,6] <<- pvalue(cont_smoke_mss)
  rm(cont_smoke_mss)
  gc()
  
  cont_mss_smoke <- independence_test(
    cohort ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,10] <<- pvalue(cont_mss_smoke)
  rm(cont_mss_smoke)
  
  
  #ancestry vs prior treatment
  results_ancest_prior <- independence_test(
    Ancestry_Label ~ Neighborhood | Prior_Treatment,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,9] <<- pvalue(results_ancest_prior)
  rm(results_ancest_prior)
  
  results_prior_ancest <- independence_test(
    Prior_Treatment ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[8,7] <<- pvalue(results_prior_ancest)
  rm(results_prior_ancest)
  
  #ancestry vs smoking
  results_ancest_smoke <- independence_test(
    Ancestry_Label ~ Neighborhood | Smoking_History,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,10] <<- pvalue(results_ancest_smoke)
  rm(results_ancest_smoke)
  gc()
  
  results_smoke_ancest <- independence_test(
    Smoking_History ~ Neighborhood | Ancestry_Label,
    data = temp,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[9,7] <<- pvalue(results_smoke_ancest)
  rm(results_smoke_ancest)
  
  
  #ancestry vs kras
  results_ancest_kras <- independence_test(
    Ancestry_Label ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,11] <<- pvalue(results_ancest_kras)
  rm(results_ancest_kras)
  
  results_kras_ancest <- independence_test(
    KRAS ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,7] <<- pvalue(results_kras_ancest)
  rm(results_kras_ancest)
  gc()
  
  
  #ancestry vs apc
  results_ancest_apc <- independence_test(
    Ancestry_Label ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,12] <<- pvalue(results_ancest_apc)
  rm(results_ancest_apc)
  
  
  
  results_apc_ancest <- independence_test(
    APC ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,7] <<- pvalue(results_apc_ancest)
  rm(results_apc_ancest)
  
  #ancestry vs tp53
  results_ancest_tp <- independence_test(
    Ancestry_Label ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,13] <<- pvalue(results_ancest_tp)
  rm(results_ancest_tp)
  gc()
  
  results_tp_ancest <- independence_test(
    TP53 ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,7] <<- pvalue(results_tp_ancest)
  rm(results_tp_ancest)
  
  #ancestry vs braf
  results_ancest_braf <- independence_test(
    Ancestry_Label ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,14] <<- pvalue(results_ancest_braf)
  rm(results_ancest_braf)
  
  results_braf_ancest <- independence_test(
    BRAF ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,7] <<- pvalue(results_braf_ancest)
  rm(results_braf_ancest)
  gc()
  
  
  #ERROR
  #ancestry vs smad
  results_ancest_smad <- independence_test(
    Ancestry_Label ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,15] <<- pvalue(results_ancest_smad)
  rm(results_ancest_smad)
  
  results_smad_ancest <- independence_test(
    SMAD4 ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,7] <<- pvalue(results_smad_ancest)
  rm(results_smad_ancest)
  
  
  #ancestry vs pik3ca
  results_ancest_pik <- independence_test(
    Ancestry_Label ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[6,16] <<- pvalue(results_ancest_pik)
  rm(results_ancest_pik)
  gc()
  
  
  results_pik_ancest <- independence_test(
    PIK3CA ~ Neighborhood | Ancestry_Label,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,7] <<- pvalue(results_pik_ancest)
  rm(results_pik_ancest)

  #smoking history vs pik
  results_pik_smoke <- independence_test(
    PIK3CA ~ Neighborhood | Smoking_History,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,10] <<- pvalue(results_pik_smoke)
  rm(results_pik_smoke)
  
  #mss vs kras
  cont_kras_mss <- independence_test(
    KRAS ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[10,6] <<- pvalue(cont_kras_mss)
  rm(cont_kras_mss)
  
  cont_mss_kras <- independence_test(
    cohort ~ Neighborhood | KRAS,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,11] <<- pvalue(cont_mss_kras)
  rm(cont_mss_kras)
  gc()
  
  # #mss vs apc
  cont_apc_mss <- independence_test(
    APC ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[11,6] <<- pvalue(cont_apc_mss)
  rm(cont_apc_mss)
  
  cont_mss_apc <- independence_test(
    cohort ~ Neighborhood | APC,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,12] <<- pvalue(cont_mss_apc)
  rm(cont_mss_apc)
  
  #mss vs tp53
  cont_tp_mss <- independence_test(
    TP53 ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[12,6] <<- pvalue(cont_tp_mss)
  rm(cont_tp_mss)
  gc()
  
  cont_mss_tp <- independence_test(
    cohort ~ Neighborhood | TP53,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,13] <<- pvalue(cont_mss_tp)
  rm(cont_mss_tp)
  
  #mss vs braf
  cont_braf_mss <- independence_test(
    BRAF ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[13,6] <<- pvalue(cont_braf_mss)
  rm(cont_braf_mss)
  
  cont_mss_braf <- independence_test(
    cohort ~ Neighborhood | BRAF,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,14] <<- pvalue(cont_mss_braf)
  rm(cont_mss_braf)
  gc()
  
  #mss vs smad4
  cont_smad_mss <- independence_test(
    SMAD4 ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[14,6] <<- pvalue(cont_smad_mss)
  rm(cont_smad_mss)
  
  cont_mss_smad <- independence_test(
    cohort ~ Neighborhood | SMAD4,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,15] <<- pvalue(cont_mss_smad)
  rm(cont_mss_smad)
  
  #mss vs pik3ca
  cont_pik_mss <- independence_test(
    PIK3CA ~ Neighborhood | cohort,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[15,6] <<- pvalue(cont_pik_mss)
  rm(cont_pik_mss)
  gc()
  
  cont_mss_pik <- independence_test(
    cohort ~ Neighborhood | PIK3CA,
    data = temp_mss_only,
    distribution = approximate(nresample = 1000000)
  )
  final_pvals[5,16] <<- pvalue(cont_mss_pik)
  rm(cont_mss_pik)
  gc()
  
}
neigh_vs_var_control_var(data20, data20_mss_only)

#initialize qvals matrix
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

#calculate qvals
qvals <- function(df){
  
  pval_vector <- df %>%
    select(-1) %>%
    pivot_longer(everything(), values_to = "val") %>%
    filter(!is.na(val), val != "") %>%
    pull(val)
  
  pval_vector <- as.numeric(pval_vector)
  #print(pval_vector)
  
  qval_vector <- p.adjust(pval_vector, method = "BH")
  #print(qval_vector)
  
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
  final_qvals[4,10] <<- qval_vector[51]
  final_qvals[4,11] <<- qval_vector[52]
  final_qvals[4,12] <<- qval_vector[53]
  final_qvals[4,13] <<- qval_vector[54]
  final_qvals[4,14] <<- qval_vector[55]
  final_qvals[4,15] <<- qval_vector[56]
  final_qvals[4,16] <<- qval_vector[57]
  final_qvals[5,2] <<- qval_vector[58]
  final_qvals[5,3] <<- qval_vector[59]
  final_qvals[5,4] <<- qval_vector[60]
  final_qvals[5,5] <<- qval_vector[61]
  final_qvals[5,7] <<- qval_vector[62]
  final_qvals[5,8] <<- qval_vector[63]
  final_qvals[5,9] <<- qval_vector[64]
  final_qvals[5,10] <<- qval_vector[65]
  final_qvals[6,2] <<- qval_vector[66]
  final_qvals[6,3] <<- qval_vector[67]
  final_qvals[6,4] <<- qval_vector[68]
  final_qvals[6,5] <<- qval_vector[69]
  final_qvals[6,6] <<- qval_vector[70]
  final_qvals[6,8] <<- qval_vector[71]
  final_qvals[6,9] <<- qval_vector[72]
  final_qvals[6,10] <<- qval_vector[73]
  final_qvals[6,11] <<- qval_vector[74]
  final_qvals[6,12] <<- qval_vector[75]
  final_qvals[6,13] <<- qval_vector[76]
  final_qvals[6,14] <<- qval_vector[77]
  final_qvals[6,15] <<- qval_vector[78]
  final_qvals[6,16] <<- qval_vector[79]
  final_qvals[7,2] <<- qval_vector[80]
  final_qvals[7,3] <<- qval_vector[81]
  final_qvals[7,4] <<- qval_vector[82]
  final_qvals[7,5] <<- qval_vector[83]
  final_qvals[7,6] <<- qval_vector[84]
  final_qvals[7,7] <<- qval_vector[85]
  final_qvals[7,9] <<- qval_vector[86]
  final_qvals[7,10] <<- qval_vector[87]
  final_qvals[7,11] <<- qval_vector[88]
  final_qvals[7,12] <<- qval_vector[89]
  final_qvals[7,13] <<- qval_vector[90]
  final_qvals[7,14] <<- qval_vector[91]
  final_qvals[7,15] <<- qval_vector[92]
  final_qvals[7,16] <<- qval_vector[93]
  final_qvals[8,2] <<- qval_vector[94]
  final_qvals[8,3] <<- qval_vector[95]
  final_qvals[8,4] <<- qval_vector[96]
  final_qvals[8,5] <<- qval_vector[97]
  final_qvals[8,6] <<- qval_vector[98]
  final_qvals[8,7] <<- qval_vector[99]
  final_qvals[8,8] <<- qval_vector[100]
  final_qvals[8,10] <<- qval_vector[101]
  final_qvals[8,11] <<- qval_vector[102]
  final_qvals[8,12] <<- qval_vector[103]
  final_qvals[8,13] <<- qval_vector[104]
  final_qvals[8,14] <<- qval_vector[105]
  final_qvals[8,15] <<- qval_vector[106]
  final_qvals[8,16] <<- qval_vector[107]
  final_qvals[9,2] <<- qval_vector[108]
  final_qvals[9,3] <<- qval_vector[109]
  final_qvals[9,4] <<- qval_vector[110]
  final_qvals[9,5] <<- qval_vector[111]
  final_qvals[9,6] <<- qval_vector[112]
  final_qvals[9,7] <<- qval_vector[113]
  final_qvals[9,8] <<- qval_vector[114]
  final_qvals[9,9] <<- qval_vector[115]
  final_qvals[9,11] <<- qval_vector[116]
  final_qvals[9,12] <<- qval_vector[117]
  final_qvals[9,13] <<- qval_vector[118]
  final_qvals[9,14] <<- qval_vector[119]
  final_qvals[9,15] <<- qval_vector[120]
  final_qvals[9,16] <<- qval_vector[121]
  final_qvals[10,2] <<- qval_vector[122]
  final_qvals[10,3] <<- qval_vector[123]
  final_qvals[10,4] <<- qval_vector[124]
  final_qvals[10,5] <<- qval_vector[125]
  final_qvals[10,7] <<- qval_vector[126]
  final_qvals[10,8] <<- qval_vector[127]
  final_qvals[10,9] <<- qval_vector[128]
  final_qvals[10,10] <<- qval_vector[129]
  final_qvals[10,12] <<- qval_vector[130]
  final_qvals[10,13] <<- qval_vector[131]
  final_qvals[10,14] <<- qval_vector[132]
  final_qvals[10,15] <<- qval_vector[133]
  final_qvals[10,16] <<- qval_vector[134]
  final_qvals[11,2] <<- qval_vector[135]
  final_qvals[11,3] <<- qval_vector[136]
  final_qvals[11,4] <<- qval_vector[137]
  final_qvals[11,5] <<- qval_vector[138]
  final_qvals[11,7] <<- qval_vector[139]
  final_qvals[11,8] <<- qval_vector[140]
  final_qvals[11,9] <<- qval_vector[141]
  final_qvals[11,10] <<- qval_vector[142]
  final_qvals[11,11] <<- qval_vector[143]
  final_qvals[11,13] <<- qval_vector[144]
  final_qvals[11,14] <<- qval_vector[145]
  final_qvals[11,15] <<- qval_vector[146]
  final_qvals[11,16] <<- qval_vector[147]
  final_qvals[12,2] <<- qval_vector[148]
  final_qvals[12,3] <<- qval_vector[149]
  final_qvals[12,4] <<- qval_vector[150]
  final_qvals[12,5] <<- qval_vector[151]
  final_qvals[12,7] <<- qval_vector[152]
  final_qvals[12,8] <<- qval_vector[153]
  final_qvals[12,9] <<- qval_vector[154]
  final_qvals[12,10] <<- qval_vector[155]
  final_qvals[12,11] <<- qval_vector[156]
  final_qvals[12,12] <<- qval_vector[157]
  final_qvals[12,14] <<- qval_vector[158]
  final_qvals[12,15] <<- qval_vector[159]
  final_qvals[12,16] <<- qval_vector[160]
  final_qvals[13,2] <<- qval_vector[161]
  final_qvals[13,3] <<- qval_vector[162]
  final_qvals[13,4] <<- qval_vector[163]
  final_qvals[13,5] <<- qval_vector[164]
  final_qvals[13,7] <<- qval_vector[165]
  final_qvals[13,8] <<- qval_vector[166]
  final_qvals[13,9] <<- qval_vector[167]
  final_qvals[13,10] <<- qval_vector[168]
  final_qvals[13,11] <<- qval_vector[169]
  final_qvals[13,12] <<- qval_vector[170]
  final_qvals[13,13] <<- qval_vector[171]
  final_qvals[13,15] <<- qval_vector[172]
  final_qvals[13,16] <<- qval_vector[173]
  final_qvals[14,2] <<- qval_vector[174]
  final_qvals[14,3] <<- qval_vector[175]
  final_qvals[14,4] <<- qval_vector[176]
  final_qvals[14,5] <<- qval_vector[177]
  final_qvals[14,7] <<- qval_vector[178]
  final_qvals[14,8] <<- qval_vector[179]
  final_qvals[14,9] <<- qval_vector[180]
  final_qvals[14,10] <<- qval_vector[181]
  final_qvals[14,11] <<- qval_vector[182]
  final_qvals[14,12] <<- qval_vector[183]
  final_qvals[14,13] <<- qval_vector[184]
  final_qvals[14,14] <<- qval_vector[185]
  final_qvals[14,16] <<- qval_vector[186]
  final_qvals[15,2] <<- qval_vector[187]
  final_qvals[15,3] <<- qval_vector[188]
  final_qvals[15,4] <<- qval_vector[189]
  final_qvals[15,5] <<- qval_vector[190]
  final_qvals[15,7] <<- qval_vector[191]
  final_qvals[15,8] <<- qval_vector[192]
  final_qvals[15,9] <<- qval_vector[193]
  final_qvals[15,10] <<- qval_vector[194]
  final_qvals[15,11] <<- qval_vector[195]
  final_qvals[15,12] <<- qval_vector[196]
  final_qvals[15,13] <<- qval_vector[197]
  final_qvals[15,14] <<- qval_vector[198]
  final_qvals[15,15] <<- qval_vector[199]
}
qvals(final_pvals)

#save to csv if so desired
write.csv(final_qvals, "confound_qvals.csv", row.names = FALSE)
