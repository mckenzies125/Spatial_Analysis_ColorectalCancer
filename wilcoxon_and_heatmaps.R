#Script Summary:
# 1. Perform Wilcoxon-Ranked Sum Testing to quantify which neighborhoods or counties had the most/if any dissimilarity
# 2. Visualize this in the form of a heat map


#The raw pvals/qvals from this analysis were not saved. However, the heatmaps were. They can be found at 
# Graphs -> Heatmaps

library(plyr)
library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(coin)

#data!
data <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_0.csv")
data_neigh_mss_only <- read.csv("~\\clinical_and_genomic\\genomic_nyc_mss_colorec_20.csv")
data_neigh <- read.csv("~\\clinical_and_genomic\\genomic_nyc_nonmss_colorec_20.csv")

names(data_neigh_mss_only)[13] <- "Yost_Index"
names(data_neigh)[13] <- "Yost_Index"
names(data)[11] <- "Yost_Index"


#everything must be a factor
data_neigh$Ancestry_Label = as.factor(data_neigh$Ancestry_Label)
data_neigh$Neighborhood = as.factor(data_neigh$Neighborhood)
data_neigh_mss_only$KRAS = as.factor(data_neigh_mss_only$KRAS)
data_neigh_mss_only$APC = as.factor(data_neigh_mss_only$APC)
data_neigh_mss_only$TP53 = as.factor(data_neigh_mss_only$TP53)
data_neigh_mss_only$BRAF = as.factor(data_neigh_mss_only$BRAF)
data_neigh_mss_only$SMAD4 = as.factor(data_neigh_mss_only$SMAD4)
data_neigh_mss_only$PIK3CA = as.factor(data_neigh_mss_only$PIK3CA)
data_neigh_mss_only$Primary_Tumor_Location_Final = as.factor(data_neigh_mss_only$Primary_Tumor_Location_Final)
data_neigh$cohort = as.factor(data_neigh$cohort)
data$Ancestry_Label = as.factor(data$Ancestry_Label)

#Descriptions of how this works can be found in this first function:
yost_county <- function(df){
  data_unique <- unique(df$County) #extract all county names
  combos <- combn(data_unique, 2, simplify = FALSE) #get all unique combinations of counties (no double pairings ie Bronx - Bronx)
  wilcoxon <- function(county1, county2, df2){ #this function performs wilcoxon pairwise tests on the counties
    subset_df_counties <- df2[df2$County %in% c(county1, county2), ] #we pass in 2 county names. Extract only the patients from those counties
    subset_df_counties$County <- droplevels(factor(subset_df_counties$County)) #continuation of previous line
    
    testing <- wilcox_test(Yost_Index ~ County, data = subset_df_counties, distribution = "exact") #perform actual wilcoxon test
    pval <- as.numeric(pvalue(testing)) #ensure pval is a number
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval)) #return a df of pvals per county
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) wilcoxon(x[1], x[2], df))) #call the wilcoxon function for EVERY county pairing we made
  
  #mirror the heatmap
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      dplyr::rename(County_1 = County_2, County_2 = County_1) #specify namespace in case the rename function is hidden from plyr
  )
  
  #adjust p vals
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  #construct qval matrix WHILE MAINTAINING ORDER
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  #construct melt (heatmap) plot
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #trans = "log10", #log if desired...only recommended for yost index
      low = "#b30000",
      high = "#fef0d9",
      limits = c(1e-15, 1e-5),
      oob = scales::squish,
      name = "q_value") + #if log is taken, change this legend title to "log_q_value
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), ) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of County-Level Differences\nin Yost Index of Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "County 1",
         y = "County 2")
}
yost_county(data_neigh)


yost_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  wilcoxon <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- wilcox_test(Yost_Index ~ Neighborhood, data = subset_df_counties, distribution = "exact")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) wilcoxon(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      dplyr::rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  #order = c("Pelham - Throgs Neck", "Bedford Stuyvesant - Crown Heights", "Bensonhurst - Bay Ridge", "Borough Park", "Canarsie - Flatlands", "Coney Island - Sheepshead Bay", "Downtown - Heights - Park Slope", "East Flatbush - Flatbush", "Gramercy Park - Murray Hill", "Union Square - Lower East Side", "Upper East Side", "Upper West Side", "Bayside - Little Neck", "Flushing - Clearview", "Jamaica", "Long Island City - Astoria", "Ridgewood - Forest Hills", "Rockaway", "Southeast Queens", "Southwest Queens", "South Beach - Tottenville", "Stapleton - St. George", "Willowbrook")
  #q_matrix <- q_matrix[order, order]
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.01, 0.06),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), ) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin Yost Index of Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2") 
}
yost_neigh(data_neigh)



age_neigh <- function(df){
  data_unique <- unique(data_neigh$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  wilcoxon <- function(county1, county2, df2){
    subset_df_counties <- df2[df2$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- wilcox_test(age_at_diag ~ Neighborhood, data = subset_df_counties, distribution = "exact")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) wilcoxon(x[1], x[2], data_neigh)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  print(pairs$p_value)
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  print(pairs$q_value)
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.01, 0.5),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin Age of Diagnosis of Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
age_neigh(data_neigh)


bmi_county <- function(df){
  data_unique <- unique(data_neigh$County)
  combos <- combn(data_unique, 2, simplify = FALSE)
  wilcoxon <- function(county1, county2, df2){
    subset_df_counties <- df2[df2$County %in% c(county1, county2), ]
    subset_df_counties$County <- droplevels(factor(subset_df_counties$County))
    
    testing <- wilcox_test(First_BMI_Measurement ~ County, data = subset_df_counties, distribution = "exact")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) wilcoxon(x[1], x[2], data_neigh)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      limits = c(1e-7, 0.1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of County-Level Differences\nin First Body Mass Index Measurement of Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "County 1",
         y = "County 2")
}
bmi_county(data)


bmi_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  wilcoxon <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- wilcox_test(First_BMI_Measurement ~ Neighborhood, data = subset_df_counties, distribution = "exact")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) wilcoxon(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      dplyr::rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.01, 0.06),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), ) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin First Body Mass Index Measurements of Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
bmi_neigh(data_neigh)


ancest_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(Ancestry_Label ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      dplyr::rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.01, 0.05),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), ) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin Ancestry Label Distribution, from Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
ancest_neigh(data_neigh)


kras_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(KRAS ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "log_q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin KRAS Mutation Distribution, from Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
kras_neigh(data_neigh_mss_only)


kras_county <- function(df){
  data_unique <- unique(data_neigh_mss_only$County)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df2[df2$County %in% c(county1, county2), ]
    subset_df_counties$County <- droplevels(factor(subset_df_counties$County))
    
    testing <- chisq_test(KRAS ~ County, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], data_neigh_mss_only)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "log_q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin KRAS Mutation Distribution, from Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
kras_county(data_neigh_mss_only)


apc_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(APC ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10",
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin APC Mutation Distribution, from Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
apc_neigh(data_neigh_mss_only)



tp53_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(TP53 ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #trans = "log10",
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin TP53 Mutation Distribution, from Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
tp53_neigh(data_neigh_mss_only)



braf_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(BRAF ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin BRAF Mutation Distribution in Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
braf_neigh(data_neigh_mss_only)



smad_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(SMAD4 ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin SMAD4 Mutation Distribution in Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
smad_neigh(data_neigh_mss_only)



pik_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(PIK3CA ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      #colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      low = "#b30000",
      high = "#fef0d9",
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences\nin PIK3CA Mutation Distribution in Patients with Colorectal Cancer in the\nMSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
pik_neigh(data_neigh_mss_only)



mss_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(cohort ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences in\nMicrostatellite and Nonmicrosatellite Stable Distribution in Patients with Colorectal\nCancer in the MSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
mss_neigh(data_neigh_mss_only)



loc_neigh <- function(df){
  data_unique <- unique(df$Neighborhood)
  combos <- combn(data_unique, 2, simplify = FALSE)
  chi2 <- function(county1, county2, df2){
    subset_df_counties <- df[df$Neighborhood %in% c(county1, county2), ]
    subset_df_counties$Neighborhood <- droplevels(factor(subset_df_counties$Neighborhood))
    
    testing <- chisq_test(Primary_Tumor_Location_Final ~ Neighborhood, data = subset_df_counties, distribution = "asymptotic")
    pval <- as.numeric(pvalue(testing))
    return(data.frame(County_1 = county1, County_2 = county2, p_value = pval))
  }
  pairs_results <- do.call(rbind, lapply(combos, function(x) chi2(x[1], x[2], df)))
  
  pairs <- bind_rows(
    pairs_results,
    pairs_results %>%
      rename(County_1 = County_2, County_2 = County_1)
  )
  
  pairs$q_value <- p.adjust(pairs$p_value, method = "BH")
  
  q_matrix <- acast(pairs, County_1 ~ County_2, value.var = "q_value")
  
  melt(q_matrix) |>
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
      #trans = "log10", 
      limits = c(0.05, 1),
      oob = scales::squish,
      name = "q_value") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q-Values for the Wilcoxon-Pairwise Comparisons of Neighborhood-Level Differences in\nPrimary Tumor Location of in Patients with Colorectal\nCancer in the MSK-IMPACT Data Set",
         x = "Neighborhood 1",
         y = "Neighborhood 2")
}
loc_neigh(data_neigh)




