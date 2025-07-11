# Selection of terrain-based variables
# Following the DATABASE R.script, this script aims at selecting the terrain-based variables to later reduce the dimensionality of the dataset using a PCA.

# Library
library(FactoMineR)
library(dplyr)
library(tidyr)
library(tidyverse)
library(corrplot)
library(RColorBrewer)
library(missMDA)
library(factoextra)
library(sf)
library(data.table)
library(ggplot2)
library(reshape2)



# STEP 1 : MERGE CATEGORICAL (vegetation) AND NUMERICAL VARIABLES (terrain-based)



files <- list.files("/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ANALYSIS/Data_sampling_Dispersal_period/Modified_date_hester", pattern = ".rds", full.names = TRUE)
PCA_eagle <- data.frame()

# Loop through the files
for (ind in seq_along(files)) {
      print(ind)
      
      # Step 1.1 : load and read the data in R 
      eagle <- readRDS(files[ind])
      
      # Step 1.2 : Transform categorical variables into continuous variables 
      # Filter the dataset to keep the individuals with more than 30 locations
      eagle_subset <- eagle %>%
            group_by(individual.local.identifier) %>%
            filter(n() >= 30) %>%
            dplyr::select(CLC, individual.local.identifier.per.week)
      # Count the location per type of land cover and per individuals
      eagle_counts_by_ind <- eagle_subset %>% 
            group_by(individual.local.identifier.per.week, CLC) %>% 
            summarise(count = n(), .groups = "drop") %>%
            ungroup()
      # Remove the NA row
      eagle_na_omit <- na.omit(eagle_counts_by_ind)
      # Group by individual and add a column with the total count per individual
      eagle0 <- eagle_na_omit %>%
            group_by(individual.local.identifier.per.week) %>%
            mutate(total_count = sum(count)) %>%
            ungroup()
      # Add a column with the percentage of points per category per individual
      eagle1 <- eagle0 %>%
            mutate(percentage = count / total_count * 100)
      # Drop the geometry column and the undesired columns 
      eagle2 <- st_drop_geometry(eagle1)
      eagle3 <- eagle2[, !names(eagle2) %in% c("total_count", "count")]
      # Create a new data.frame with for each columns the name of Land cover type 
      eagle4 <- eagle3 %>% 
            group_by(CLC) %>%
            mutate(percentage = as.character(percentage)) %>%
            pivot_longer(-c("CLC", "individual.local.identifier.per.week")) %>%
            pivot_wider(id_cols = individual.local.identifier.per.week, names_from = c(CLC)) %>%
            type.convert(as.is = TRUE)
      
      # Impute 0 value to the NA (nb : it means that there is 0 percent of GPS location in the category and not that we don't have any data)
      eagle4[is.na(eagle4)] <- 0
      # Convert to numerical data.frame
      eagle5 <- as.data.frame(eagle4)
      eagle5 <- eagle5 %>% 
            mutate(across(-individual.local.identifier.per.week, as.numeric))
      eagle5$individual.local.identifier.per.week <- eagle4$individual.local.identifier.per.week
      
      # Step 1.3 : Merge categories 
      # Low vegetation 
      low_vegetation <- c("Permanent herbaceous", "Periodically herbaceous", "Non and sparsely vegetated")
      columns <- intersect(low_vegetation, colnames(eagle5))
      if (length(columns) > 0) {
            eagle5$Low_vegetation <- rowSums(eagle5[, columns], na.rm = TRUE)
      }
      # Forest
      forest <- c("Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees")
      forest_columns <- intersect(forest, colnames(eagle5))
      if (length(forest_columns) > 0) {
            eagle5$Woody_broadleaved_tree <- rowSums(eagle5[, forest_columns], na.rm = TRUE)
      }
      # Remove the previous category that have been merged
      vegetation <- eagle5[, !(names(eagle5) %in% c("Periodically herbaceous", "Permanent herbaceous", "Non and sparsely vegetated", "Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees"))]
      
      # Step 1.4 : Derivative variables of the landscape (min, max, means, var) and reduction of the data.frame to one row per individual
      # Filter the dataset to keep the individuals with more than 30 locations
      sample0 <- eagle %>%
            group_by(individual.local.identifier.per.week) %>%
            filter(n() >= 30) %>%
            st_drop_geometry()
      # Remove unnecessary columns and drop geometry 
      sample1 <- sample0[, !(names(sample0) %in% c("timestamp", "location.long", "location.lat", "height.above.ellipsoid", "days_since_emig", "stage", "CLC"))]
      
      # Step 1.5 : give a weight per individuals 
      ind_list1 <- split(sample1, sample1$individual.local.identifier.per.week)
      gpsPoints_Ind <- as.data.frame(sapply(ind_list1, nrow)) # n. gps point per ind
      gpsPoints_Ind$individual.local.identifier.per.week <- rownames(gpsPoints_Ind)
      gpsPoints_Ind$percentage <- as.numeric(gpsPoints_Ind[,1]) / nrow(eagle_subset) # n. gps points per individual divided by tot n of gps points
      gpsPoints_Ind$weight <- abs(gpsPoints_Ind$percentage - 1) # absolute value (reversed weight of the individuals)
      #sum(gpsPoints_perInd$percentage)
      
      # Step 1.6 : calculate the mean, var, min, max for each variables 
      summary <- rbindlist(lapply(ind_list1, function(indi_data) {
            indi_summary <- indi_data %>%
                  select_if(is.numeric) %>%
                  summarise(across(everything(), .fns = list(
                        mean = ~ mean(., na.rm = TRUE),
                        var = ~ var(., na.rm = TRUE),
                        max = ~ max(., na.rm = TRUE),
                        min = ~ min(., na.rm = TRUE),
                        q1 = ~ quantile(., probs = 0.25, na.rm = TRUE)
                  )))
            return(indi_summary)
      }))
      summary <- merge(summary, gpsPoints_Ind[,2:4], by="individual.local.identifier.per.week", all.x=T)
      
      # Step 1.7 : merge the dataset with the Corine Land Cover dataset 
      tot_eagle <- cbind(summary, vegetation)
      
      # Step 1.8 : suppress the variables pre-selected
      # Remove column with name in double
      ge <- tot_eagle[,-1]
      
      # Scale the continuous variables
      scale_continuous <- ge%>%
            select_if(is.numeric) %>%
            scale() %>%
            as.data.frame()
      
      scale_continuous$weight <- NULL
      scale_continuous_df <- cbind(scale_continuous, data.frame(weight=ge$weight))
      
      # Add the column name
      scale_continuous_df$individual.local.identifier.per.week <- tot_eagle$individual.local.identifier.per.week
      
      # Rbind the results
      PCA_eagle <- bind_rows(PCA_eagle, scale_continuous_df)
}



# STEP 2 : CORRELATION MATRIX --------------------------------------------------


# replace NA by 0
correlation <- PCA_eagle
correlation[is.na(correlation)] <- 0  
df_corr <- correlation %>%
      select(-individual.local.identifier.per.week, -weight)

# identify and drop constant columns
zero_var_cols <- names(which(apply(df_corr, 2, sd) == 0))
if (length(zero_var_cols) > 0) {
      message("Dropping constant columns: ", paste(zero_var_cols, collapse = ", "))
      df_corr <- df_corr %>% select(-all_of(zero_var_cols))}

# compute the correlation matrix
cor.mat <- round(cor(df_corr, use = "pairwise.complete.obs"), 2)
cor.mat.long <- reshape2::melt(cor.mat)
thr <- 0.60
filt <- cor.mat.long %>%
      filter(Var1 != Var2, abs(value) >= thr)

interesting_vars <- unique(c(filt$Var1, filt$Var2))

# subset the wide matrix to only those variables
cor.mat.sub <- cor.mat[interesting_vars, interesting_vars]

# melt the reduced matrix
cor.mat.sub.long <- reshape2::melt(cor.mat.sub)
colnames(cor.mat.sub.long) <- c("Var1","Var2","value")

# visualisation
ggplot(cor.mat.sub.long, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(
            colors    = rev(brewer.pal(11, "RdBu")), 
            values    = scales::rescale(seq(-1, 1, length.out = 11)), 
            limits    = c(-1, 1),
            na.value  = "grey90"
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_fixed() +
      theme_minimal(base_family = "Avenir") +
      theme(
            axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y     = element_text(size = 8),
            axis.title      = element_blank(),
            panel.grid      = element_blank(),
            legend.position = "bottom",
            legend.text     = element_text(size = 8),
            legend.title    = element_text(size = 10)
      )



# STEP 3 : DENDOGRAM OF CORRELATED VARIABLES -----------------------------------


dist_mat <- as.dist(1 - abs(cor.mat.sub))
hc <- hclust(dist_mat, method = "average")  

# visualisation
plot(hc, main = "Cluster dendrogram of selected variables", 
     xlab = "", sub = "", cex = 0.8)

# Using both the correlation matrix and the dendogram compare the groups of correlated variables and select one of the variation of the variable between the mean, maximum, minimum, median and first quartile.
to_remove <- c("sealed", "water", "percentage", "Low-growing woody plants", "height_ground_max", "height_ground_min", "height_ground_var", "height_ground_mean", "height_ground_q1",
      "elevation_var", "elevation_mean", "elevation_q1",
      "NEAREST_WIND_TURBINE_max", "NEAREST_WIND_TURBINE_mean", "NEAREST_WIND_TURBINE_var","NEAREST_WIND_TURBINE_q1",
      "DIST_AERIALWAYS_max", "DIST_AERIALWAYS_var", "DIST_AERIALWAYS_mean", "DIST_AERIALWAYS_q1",
      "DIST_POWER_LINE_max", "DIST_POWER_LINE_mean","DIST_POWER_LINE_var","DIST_POWER_LINE_q1",
      "DIST_ROADS_var", "DIST_ROADS_mean", "DIST_ROADS_max", "DIST_ROADS_q1",
      "TRI_class_mean", "TRI_class_var", "TRI_class_max", "TRI_class_min",
      "TPI_intermediate_class_max", "TPI_intermediate_class_mean", "TPI_intermediate_class_min", "TPI_intermediate_class_var",
      "TPI_small_class_max", "TPI_small_class_mean", "TPI_small_class_min", "TPI_small_class_var",
      "TPI_intermediate_brut_mean", "TPI_intermediate_brut_var", "TPI_intermediate_brut_min", "TPI_intermediate_brut_max", 
      "TPI_small_brut_max", "TPI_small_brut_var", "TPI_small_brut_min", "TPI_small_brut_mean",
      "DIST_TO_RIDGE_min", "DIST_TO_RIDGE_var", "DIST_TO_RIDGE_mean",
      "NEAREST_BUILT_SURFACE_max","NEAREST_BUILT_SURFACE_var","NEAREST_BUILT_SURFACE_mean", "NEAREST_BUILT_SURFACE_q1",
      "TRI_brut_var", "TRI_brut_max", "TRI_brut_q1",
      "slope_mean", "slope_var","slope_q1",
      "eastness_var", "eastness_min", "eastness_max", "eastness_mean",
      "northness_mean", "northness_q1", "northness_var", "northness_min",
      "DIST_TO_TALL_BUILDING_var", "DIST_TO_TALL_BUILDING_min", "DIST_TO_TALL_BUILDING_mean","DIST_TO_TALL_BUILDING_max",
      "BUILDING_HEIGHT_mean", "BUILDING_HEIGHT_var", "BUILDING_HEIGHT_max",
      "DENSITY_var", "DENSITY_max", "height_ground_max", "height_ground_min", "height_ground_var", "height_ground_mean", "height_ground_q1", "Snow and ice")

df_uncorr2 <- df_corr %>% 
      select(-any_of(to_remove))

eagle_weight0 <- PCA_eagle[, !(colnames(PCA_eagle) %in% "weight")]

# Turn bird IDs into row names
eagle_weight0$individual.local.identifier.per.week <-
      make.unique(as.character(eagle_weight0$individual.local.identifier.per.week))

PCA_data <- eagle_weight0 %>%
      column_to_rownames(var = "individual.local.identifier.per.week")

# Keep only your de‐correlated variables
PCA_data <- PCA_data[, colnames(df_uncorr2)]

# Run the PCA on that reduced set
PCA_1 <- FactoMineR::PCA(
      PCA_data,
      row.w = PCA_eagle$weight,
      ncp   = ncol(PCA_data),
      graph = FALSE)

plot(
      PCA_1,
      choix = "var",       # <-- switch to variables
      axes  = c(1, 2),     # PC1 vs PC2
      cex   = 0.8,         # shrink labels if you like
      main  = "Variables – PCA axes 1 & 2")

# Remove additionnal variables
to_remove2 <- c("Water", "Sealed",
      "eastness_var",  "eastness_min", "eastness_max",  "eastness_mean", "eastness_q1",
      "northness_mean","northness_max","northness_var",
      "DENSITY_mean",
      "slope_max",
      "DIST_TO_TALL_BUILDING_q1",
      "slope_min",     "TPI_small_brut_q1")

PCA_data2 <- PCA_data[, !(colnames(PCA_data) %in% to_remove2)]

# PCA and visualisation
PCA_2 <- FactoMineR::PCA(
      PCA_data2,
      row.w      = PCA_eagle$weight,
      ncp        = ncol(PCA_data2),
      scale.unit = TRUE,
      graph      = FALSE
)

plot(
      PCA_2,
      choix = "var",       
      axes  = c(1, 2),     
      cex   = 0.8,        
      main  = "Variables – PCA axes 1 & 2")




# STEP 4 : LOOP THROUGHT ALL THE FILES TO PREPARE THE DATASET USING PRESELECTED VARIABLES



list_file <- list.files("/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ANALYSIS/Data_sampling_Dispersal_period/Modified_date_hester", pattern = ".rds", full.names = TRUE)
PCA_eagle_final <- data.frame()

# Loop through the files
for (ind in seq_along(list_file)) {
      print(ind)
      
      # Step 4.1 : load and read the data in R 
      data <- readRDS(list_file[ind])
      
      # Step 4.2 : Transform categorical variables into continuous variables 
      # Filter the dataset to keep the individuals with more than 30 locations
      data_subset <- data %>%
            group_by(individual.local.identifier) %>%
            filter(n() >= 30) %>%
            dplyr::select(CLC, individual.local.identifier.per.week)
      # Count the location per type of land cover and per individuals
      data_counts_by_ind <- data_subset %>% 
            group_by(individual.local.identifier.per.week, CLC) %>% 
            summarise(count = n(), .groups = "drop") %>%
            ungroup()
      # Remove the NA row
      data_na_omit <- na.omit(data_counts_by_ind)
      # Group by individual and add a column with the total count per individual
      data_0 <- data_na_omit %>%
            group_by(individual.local.identifier.per.week) %>%
            mutate(total_count = sum(count)) %>%
            ungroup()
      # Add a column with the percentage of points per category per individual
      data_1 <- data_0 %>%
            mutate(percentage = count / total_count * 100)
      # Drop the geometry column and the undesired columns 
      data_2 <- st_drop_geometry(data_1)
      data_3 <- data_2[, !names(data_2) %in% c("total_count", "count")]
      # Create a new data.frame with for each columns the name of Land cover type 
      data_4 <- data_3 %>% 
            group_by(CLC) %>%
            mutate(percentage = as.character(percentage)) %>%
            pivot_longer(-c("CLC", "individual.local.identifier.per.week")) %>%
            pivot_wider(id_cols = individual.local.identifier.per.week, names_from = c(CLC)) %>%
            type.convert(as.is = TRUE)
      
      # Impute 0 value to the NA (nb : it means that there is 0 percent of GPS location in the category and not that we don't have any data)
      data_4[is.na(data_4)] <- 0
      # Convert to numerical data.frame
      data_5 <- as.data.frame(data_4)
      data_5 <- data_5 %>% 
            mutate(across(-individual.local.identifier.per.week, as.numeric))
      data_5$individual.local.identifier.per.week <- data_4$individual.local.identifier.per.week
      
      # Step 4.3 : Merge categories 
      # Low vegetation 
      columns_low_vegetation <- c("Permanent herbaceous", "Periodically herbaceous", "Non and sparsely vegetated")
      existing_columns <- intersect(columns_low_vegetation, colnames(data_5))
      if (length(existing_columns) > 0) {
            data_5$Low_vegetation <- rowSums(data_5[, existing_columns], na.rm = TRUE)
      }
      # Forest
      columns_forest <- c("Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees")
      existing_forest_columns <- intersect(columns_forest, colnames(data_5))
      if (length(existing_forest_columns) > 0) {
            data_5$Woody_broadleaved_tree <- rowSums(data_5[, existing_forest_columns], na.rm = TRUE)
      }
      # Remove the previous category that have been merged
      eagle_vegetation <- data_5[, !(names(data_5) %in% c("Periodically herbaceous", "Permanent herbaceous", "Non and sparsely vegetated", "Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees"))]
      
      # Step 2.4 : Derivative variables of the landscape (min, max, means, var) and reduction of the data.frame to one row per individual
      # Filter the dataset to keep the individuals with more than 30 locations
      sample <- data %>%
            group_by(individual.local.identifier.per.week) %>%
            filter(n() >= 30) %>%
            st_drop_geometry()
      # Remove unnecessary columns and drop geometry 
      sample_1 <- sample[, !(names(sample) %in% c("timestamp", "location.long", "location.lat", "height.above.ellipsoid", "days_since_emig", "stage", "CLC"))]
      
      # Step 4.5 : give a weight per individuals 
      ind_list <- split(sample_1, sample_1$individual.local.identifier.per.week)
      gpsPoints_perInd <- as.data.frame(sapply(ind_list, nrow)) # n. gps point per ind
      gpsPoints_perInd$individual.local.identifier.per.week <- rownames(gpsPoints_perInd)
      gpsPoints_perInd$percentage <- as.numeric(gpsPoints_perInd[,1]) / nrow(data_subset) # n. gps points per individual divided by tot n of gps points
      gpsPoints_perInd$weight <- abs(gpsPoints_perInd$percentage - 1) # absolute value (reversed weight of the individuals)
      #sum(gpsPoints_perInd$percentage)
      
      
      # Step 4.6 : calculate the mean, var, min, max for each variables 
      individual_summary <- rbindlist(lapply(ind_list, function(ind_data) {
            ind_summary <- ind_data %>%
                  select_if(is.numeric) %>%
                  summarise(across(everything(), .fns = list(
                        mean = ~ mean(., na.rm = TRUE),
                        var = ~ var(., na.rm = TRUE),
                        max = ~ max(., na.rm = TRUE),
                        min = ~ min(., na.rm = TRUE),
                        q1 = ~ quantile(., probs = 0.25, na.rm = TRUE)
                  )))
            return(ind_summary)
      }))
      individual_summary <- merge(individual_summary, gpsPoints_perInd[,2:4], by="individual.local.identifier.per.week", all.x=T)
      
      # Step 4.7 : merge the dataset with the Corine Land Cover dataset 
      golden_eagle_tot <- cbind(individual_summary, eagle_vegetation)
      
      # Step 4.8 : suppress the variables pre-selected
      # Remove column with name in double
      eagle_1 <- golden_eagle_tot[,-1]
      
      # Scale the continuous variables
      eagle_scale_continuous_var <- eagle_1%>%
            select_if(is.numeric) %>%
            scale() %>%
            as.data.frame()
      
      eagle_scale_continuous_var$weight <- NULL
      eagle_scale_continuous_var_df <- cbind(eagle_scale_continuous_var, data.frame(weight=eagle_1$weight))
      
      # Remove the columns with NA values : building_height_min and density_min
      df <- eagle_scale_continuous_var_df[, !(colnames(eagle_scale_continuous_var_df) %in% c(
            "BUILDING_HEIGHT_min", "DIST_TO_RIDGE_min", 
            "step_dist_max", "step_dist_var", "step_dist_mean", "step_dist_min", "step_dist_q1", 
            "height_ground_max", "height_ground_min", "height_ground_var", "height_ground_mean", "height_ground_q1",
            "elevation_var", "elevation_mean", "elevation_q1",
            "slope_mean", "slope_var", "slope_max", "slope_min", "slope_q1",
            "TRI_brut_var", "TRI_brut_q1", "TRI_brut_max",
            "DIST_TO_RIDGE_mean", "DIST_TO_RIDGE_var",
            "northness_mean", "northness_var", "northness_q1", "northness_max", "northness_min",
            "eastness_max", "eastness_mean", "eastness_min", "eastness_var", "eastness_q1",
            "TPI_small_brut_mean", "TPI_small_brut_var", "TPI_small_brut_max", "TPI_small_brut_q1", "TPI_small_brut_min",
            "DIST_AERIALWAYS_var", "DIST_AERIALWAYS_max", "DIST_AERIALWAYS_mean", "DIST_AERIALWAYS_mean", "DIST_AERIALWAYS_q1",
            "NEAREST_BUILT_SURFACE_max","NEAREST_BUILT_SURFACE_var", "NEAREST_BUILT_SURFACE_mean", "NEAREST_BUILT_SURFACE_q1",
            "NEAREST_WIND_TURBINE_mean", "NEAREST_WIND_TURBINE_var", "NEAREST_WIND_TURBINE_q1",
            "Snow and ice", "Water", "Sealed", "Low-growing woody plants",
            "BUILDING_HEIGHT_var", "BUILDING_HEIGHT_max", "BUILDING_HEIGHT_mean", "BUILDING_HEIGHT_q1", 
            "DIST_TO_TALL_BUILDING_max", "DIST_TO_TALL_BUILDING_mean", "DIST_TO_TALL_BUILDING_min", "DIST_TO_TALL_BUILDING_var", "DIST_TO_TALL_BUILDING_q1",
            "DIST_POWER_LINE_max", "DIST_POWER_LINE_mean","DIST_POWER_LINE_var", "DIST_POWER_LINE_q1",
            "DIST_ROADS_mean", "DIST_ROADS_var", "DIST_ROADS_q1", "DIST_ROADS_q1", "DIST_ROADS_max",
            "DENSITY_mean", "DENSITY_min", "DENSITY_var", "DENSITY_max", "DENSITY_q1",
            "NEAREST_WIND_TURBINE_max", "percentage"))]
      
      # Add the column name
      df$individual.local.identifier.per.week <- golden_eagle_tot$individual.local.identifier.per.week
      
      # Rbind the results
      PCA_eagle_final <- rbind(PCA_eagle_final, df)
}

# Create a copy of the first dataset without the weight column as so that the column is not in the PCA
eagle_without_weight <- PCA_eagle_final[, !(colnames(PCA_eagle_final) %in% c("weight"))]

# Add row names to display the names of the birds in the PCA
PCA_data_3 <- eagle_without_weight %>%
      column_to_rownames(var = "individual.local.identifier.per.week")

# PCA plot 
PCA_3 <- FactoMineR::PCA(PCA_data_3[,1:13], row.w = PCA_eagle_final$weight, ncp = 13)
summary(PCA)
plot(
      PCA_3,
      choix = "var",       
      axes  = c(1, 2),     
      cex   = 0.8,        
      main  = "Variables – PCA axes 1 & 2"
)


#write_rds(PCA_eagle_final, "your_path_to_repository.rds")

