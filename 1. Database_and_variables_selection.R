# JUVENILE GOLDEN EAGLE DATABASE FOR TERRESTRIAL ANALYSIS OF THE LANDSCAPE 
# This R script is made to ; I. Select the eagle GPS location by temporal sequence period (fledging / emigration) (step 0 to 2), 
# II. extract from each previously constructed raster layers, the values of each pixels under the localisation of a golden eagle (step 3), 
# III. summarize the collected data per individual using several statistical tools (mean, maximum, minimum, etc ...) (step 4),
# IV. operating a data-reduction by pre-selecting variables that will be used in the PCA (step 5).
# Louise FAURE, MPI, 28.07.2025


#library
library(terra)
library(lubridate)
library(dplyr)
library(purrr)
library(sf)
library(move)
library(tidyr)
library(FactoMineR)
library(tidyverse)
library(corrplot)
library(RColorBrewer)
library(missMDA)
library(factoextra)
library(data.table)
library(ggplot2)
library(reshape2)


# STEP 0 : load the required data -------------------------------------------------------------------------------------------------------------------------
# Geoid and raster
geo <- terra::rast("us_nga_egm96_15.tif")
utm_crs <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"

# Topographic variables
elevation_25 <- rast("region-alpes-dem.tif")
slope_25 <- rast("slope_25.tif")
TRI_25 <- rast("Terrain Ruggedness Index (TRI)Corrected.tif")
TPI_25_small_scale <- rast("Topographic Position Index.tif")
INDEX_DIST_TO_RIDGE <- rast("distance_to_ridge_line_complete_version.tif")

# Human variables
CLC_ALPS <- terra::rast("CLC_longlat.tif")
BUILDING_HEIGHT <- terra::rast("building_long_lat.tif")
DIST_TO_10M_BUILDING <- rast("DIST_TO_10MBUILDING_longlat.tif")
DIST_TO_POWER_LINE <- rast("dist_powerline_longlat.tif") 
DIST_TO_ROADS <- rast("dist_road_longlat.tif") 
DIST_TO_AERIALWAYS <- rast("DIST_AERIALWAYS_longlat.tif") 
DIST_TO_NEAREST_BUILT_SURFACE <- rast("dist_to_settlement_res50.tif")
DIST_TO_WIND_TURBINE <- rast("raster_distance_wind_turbine_CRS_CLC.tif")

# set output directory 
output_dir <- "enter your output directory"

# Emigration date
emig_dates <- readRDS("fledging_emigration.rds")

# Update the date for Matsch19 (Swiss Ornithological Institute estimation)
emig_dates$emigration_dt <- as.POSIXct(emig_dates$emigration_dt)
emig_dates$emigration_dt[emig_dates$individual.local.identifier == "Matsch19 (eobs 7035)"] <- as.POSIXct("2020-03-08  09:00:08")
# Update the date for Kastelbell19 (Hester Bronnvik)
emig_dates$emigration_dt[emig_dates$individual.local.identifier == "Kastelbell19 (eobs 7034)"] <- as.POSIXct("2020-02-04 15:40:14")
emig_dates <- emig_dates %>%
  mutate(
    fledging_dt = as.POSIXct(fledging_dt, tz = "UTC"),
    emigration_dt = as.POSIXct(emigration_dt, tz = "UTC")
  )  %>%
  filter(!is.na(fledging_dt) & !is.na(emigration_dt))




# STEP 1 : list all the file ------------------------------------------------------------------------------------------------------------------------------
  # (loop them all and skip inside the loop)
  list_file_all <- list.files(
      "gps_for_Louise",
      full.names = TRUE, pattern = "\\.rds$"
    )
valid_ids <- emig_dates$individual.local.identifier





# STEP 2 : Define "w" as the temporal segmentation of the data --------------------------------------------------------------------------------------------
week_definitions <- emig_dates %>%
  mutate(week_list = map2(fledging_dt, emigration_dt, ~{
    # Week 0 = pre-dispersal
    weeks <- list(data.frame(week = 0, start_day = .x, end_day = .y))
    # Week 1 to 15 (7 days-sequences after emigration_dt)
    for (w in 1:15) {
      weeks[[w + 1]] <- data.frame(
        week = w,
        start_day = .y + days(7 * (w - 1)),
        end_day = .y + days(7 * w)
      )
    }
    bind_rows(weeks)
  })) %>%
  tidyr::unnest(week_list) %>%
  dplyr::select(individual = individual.local.identifier, week, start_day, end_day)




# STEP 3 : week x individual ----------------------------------------------------------------------------------------------------------------------------
for (w in 0:15) {
  data_week <- data.frame()
  # Loop through the files
  
  for (ind in seq_along(list_file_all)) {
    
    print(ind)
    # Load and read the data in R
    data <- readRDS(list_file_all[ind])
    id   <- unique(data$individual.local.identifier)
    # SKIP any file whose id is not in our emigration table and not “Johnsbach”
    if (!(id %in% valid_ids || grepl("Johnsbach", id))) next
    # SAMPLING BY WEEK
    data <- data[, c("timestamp", "location.long", "location.lat", "height.above.ellipsoid","individual.local.identifier" )]
    
    id <- unique(data$individual.local.identifier)
    
    week_info <- week_definitions %>%
      filter(individual == id, week == w)  # with "w" from the external loop
    
    if (nrow(week_info) == 0) next  # skip if there is no information for this individual
    
    start_day <- as.POSIXct(week_info$start_day, tz = "UTC")
    end_day   <- as.POSIXct(week_info$end_day,   tz = "UTC")
    
    data_sample <- data %>%
      filter(timestamp >= start_day & timestamp < end_day)
    if (nrow(data_sample) <= 30) {next}  # keep only individuals with more than 30 locations
    
    # Remove all double timestamp
    rows_to_delete <- unlist(sapply(move::getDuplicatedTimestamps(x = as.factor(data_sample$individual.local.identifier), timestamps = as.POSIXct(data_sample$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")), "[", -1)) #get all but the first row of each set of duplicate rows
    
    if(length(rows_to_delete) > 0){
      data_timestamp <- data_sample[-rows_to_delete,]
    } else {
      data_timestamp <- data_sample
    }
    
    # SAMPLE BY EVERY 15 MINUTES
    eagle_15m <- data_timestamp %>%
      mutate(timestamp_rounded = floor_date(timestamp, "15 minutes")) %>% 
      group_by(timestamp_rounded) %>%
      slice(1)
    
    eagle_15m <- as.data.frame(eagle_15m)
    
    # SAMPLING BY 150 METERS DISTANCE BETWEEN LOCATION POINTS
    data_move <- move(x = eagle_15m$location.long, y = eagle_15m$location.lat,
                      time = eagle_15m$timestamp,
                      proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                      data = eagle_15m, animal = eagle_15m$individual.local.identifier)
    
    data_move@data$step_dist <- c(NA, distance(data_move)) # to calculate step distances
    
    while(any(data_move$step_dist[!is.na(data_move$step_dist)] < 150)){
      toRemove <- min(which(data_move$step_dist < 150))
      data_move <- data_move[-toRemove,]
      data_move$step_dist <- c(NA, distance(data_move))
    } # to select GPS location points evry 150 meters
    
    
    # EXTRACTION OF VALUES FOR TOPOGRAPHIC & SETTLEMENT RASTERS 
    # Transform data_df into a spatial object
    data_df <- as.data.frame(data_move)
    data_df$individual.local.identifier <- rep(move::idData(data_move), nrow(data_df))
    individual_df_sf <- st_as_sf(data_df, 
                                 coords = c("location.long", "location.lat"),
                                 crs=crs("+proj=longlat +datum=WGS84 +no_defs"),
                                 remove = F)
    
    # Adjust the CRS of the individual points before extracting pixel values
    individual_df_sf <- st_transform(individual_df_sf, crs(slope_25)) 
    
    # Extract for each layers the values of the pixel under the GPS location
    elevation_values <- terra::extract(elevation_25, individual_df_sf)
    slope_values <- terra::extract(slope_25, individual_df_sf)
    TPI_small_brut_values <- terra::extract(TPI_25_small_scale, individual_df_sf)
    TRI_brut_values <- terra::extract(TRI_25, individual_df_sf)
    DIST_TO_RIDGE_values <- terra::extract(INDEX_DIST_TO_RIDGE, individual_df_sf)
    DIST_TO_NEAREST_BUILT_SURFACE_values <- terra::extract(DIST_TO_NEAREST_BUILT_SURFACE, individual_df_sf)
    
    # Create new columns for each variables
    data_topo <- cbind(individual_df_sf, data.frame(
      elevation = elevation_values$eu_dem_v11_E40N20, 
      slope = slope_values$slope_25, 
      TPI_small_brut = TPI_small_brut_values$`Topographic Position Index`, 
      TRI_brut = TRI_brut_values$`Terrain Ruggedness Index (TRI)Corrected`, 
      DIST_TO_RIDGE = DIST_TO_RIDGE_values$distance_to_ridge_line_mask,
      NEAREST_BUILT_SURFACE = DIST_TO_NEAREST_BUILT_SURFACE_values$dist_to_settlement_res50))
    
    
    # EXTRACTION OF VALUES FOR HUMAN-MADE INFRASTRUCTURES RASTERS 
    # Adjust the CRS of individual point before extracting pixel values
    individual_layer <- st_transform(data_topo, crs(CLC_ALPS))
    
    # Extract for each layers the values of the pixel under the GPS location
    CLC_values <- terra::extract(CLC_ALPS, individual_layer)
    BUILDING_HEIGHT_values <- terra::extract(BUILDING_HEIGHT, individual_layer)
    DIST_TO_10M_BUILDING_values <- terra::extract(DIST_TO_10M_BUILDING, individual_layer)
    DIST_TO_POWER_LINE_values <- terra::extract(DIST_TO_POWER_LINE, individual_layer)
    DIST_TO_ROADS_values <- terra::extract(DIST_TO_ROADS, individual_layer)
    DIST_TO_AERIALWAYS_values <- terra::extract(DIST_TO_AERIALWAYS, individual_layer)
    DIST_TO_WIND_TURBINE_values <- terra::extract(DIST_TO_WIND_TURBINE, individual_layer)
    
    # Create new columns for each variables
    eagle_human <- cbind(individual_layer, data.frame(
      CLC = CLC_values$Class_name,
      BUILDING_HEIGHT = BUILDING_HEIGHT_values$GHS_BUILT_H_ANBH_E2018_GLOBE_R2022A_54009_100_V1_0_R4_C20, 
      DIST_TO_TALL_BUILDING = DIST_TO_10M_BUILDING_values$GHS_BUILT_H_ANBH_E2018_GLOBE_R2022A_54009_100_V1_0_R4_C20, 
      DIST_POWER_LINE = DIST_TO_POWER_LINE_values$map_power_line, 
      DIST_ROADS = DIST_TO_ROADS_values$raster_distance_to_roads2, 
      DIST_AERIALWAYS = DIST_TO_AERIALWAYS_values$raster_distance_to_aerialways, 
      NEAREST_WIND_TURBINE = DIST_TO_WIND_TURBINE_values$raster_distance_wind_turbine))
    
    # Rbind the results to add rows for each individuals 
    data_week <- rbind(data_week, eagle_human)
    }
  if (nrow(data_week) > 0) {
    # suppress useless columns
    data_week_clean <- data_week[, !(names(data_week) %in% c(
      "coords.x1", "coords.x2", "sensor", "timestamps", "optional", "dem", "geoid", 
      "height_msl", "step_dist", "height.above.ellipsoid", "days_since_emig", 
      "timestamp_rounded", "stage", "from.idData.rep.1..n.locs.from....."))]
    
    # Safeguard within a file specific to a given time period (F = fledgling / D = dispersal)
    fname <- if (w == 0) "F.rds" else paste0("D", w, ".rds")
    saveRDS(data_week_clean, file = file.path(output_dir, fname))
  }
}





# STEP 4 : loop to merge categorical (vegetation) and numerical variables (terrain-based) ------------------------------------------------------------------
files <- list.files("your path to your directory", pattern = ".rds", full.names = TRUE)
PCA_eagle <- data.frame()

for (ind in seq_along(files)) {
      print(ind)
  
  # create a new column and rename individuals based on their time period
  cat("→ files", ind, ":", basename(files[ind]), "\n")
  
  eagle <- readRDS(files[ind])
  
  week_tag <- tools::file_path_sans_ext(basename(files[ind]))  # to retain the name of the file (F, D1, D2, etc ...)
  
  eagle$individual.local.identifier.per.week <- sub(
      "\\s*\\(.*\\)$", 
      paste0(" ", week_tag), 
      trimws(as.character(eagle$individual.local.identifier))
    )
  
  # TRANSFORMATION OF CATEGORICAL VALUES TO NUMERICAL VARIABLES
  # Filter the data to keep the individuals with more than 30 locations
  eagle_subset <- eagle %>%
            group_by(individual.local.identifier.per.week) %>%
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
  
  # CReate a list of all the Land Cover categories    
  all_clc_categories <- c(
        "Sealed",
        "Woody needle leaved trees",
        "Woody Broadleaved deciduous trees",
        "Woody Broadleaved evergreen trees",
        "Low-growing woody plants",
        "Permanent herbaceous",
        "Periodically herbaceous",
        "Lichens and mosses",
        "Non and sparsely vegetated",
        "Water",
        "Snow and ice")
  
   # Add missing columns filled with 0
   missing_cols <- setdiff(all_clc_categories, colnames(eagle4))
      for (col in missing_cols) {
        eagle4[[col]] <- 0}
   
   # Impute 0 value to the NA (nb : it means that there is 0 percent of GPS location in the land cover category)
   eagle4[is.na(eagle4)] <- 0
   
   # Convert to numerical data.frame
   eagle5 <- as.data.frame(eagle4)
   
   eagle5 <- eagle5 %>% 
            mutate(across(-individual.local.identifier.per.week, as.numeric))
   
   eagle5$individual.local.identifier.per.week <- eagle4$individual.local.identifier.per.week
   
   # Merge categories 
   low_vegetation <- c("Permanent herbaceous", "Periodically herbaceous") # grassland, annual crops 
   columns <- intersect(low_vegetation, colnames(eagle5))
   
   if (length(columns) > 0) {
            eagle5$Low_vegetation <- rowSums(eagle5[, columns], na.rm = TRUE)
   }
   
   forest <- c("Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees") # Forest
   forest_columns <- intersect(forest, colnames(eagle5))
   
   if (length(forest_columns) > 0) {
            eagle5$Woody_broadleaved_tree <- rowSums(eagle5[, forest_columns], na.rm = TRUE)
   }
   
   # Remove the previous category that have been merged
   vegetation <- eagle5[, !(names(eagle5) %in% c("Periodically herbaceous", "Permanent herbaceous", "Woody needle leaved trees", "Woody Broadleaved deciduous trees", "Woody Broadleaved evergreen trees"))]
      
   # SUMMARISE THE LANDSCAPE VARIABLES PER INDIVIDUALS USING THE MIN, MAX, MEAN, VAR, ETC ...
   # Filter the dataset to keep the individuals with more than 30 locations
   sample0 <- eagle %>%
            group_by(individual.local.identifier.per.week) %>%
            filter(n() >= 30) %>%
            st_drop_geometry()
   
   # Remove unnecessary columns and drop geometry 
   sample1 <- sample0[, !(names(sample0) %in% c("timestamp", "location.long", "location.lat", "height.above.ellipsoid", "days_since_emig", "stage", "CLC"))]
      
   # Attribute a weight per individuals 
   ind_list1 <- split(sample1, sample1$individual.local.identifier.per.week)
   gpsPoints_Ind <- as.data.frame(sapply(ind_list1, nrow)) # number of gps point per ind
   gpsPoints_Ind$individual.local.identifier.per.week <- rownames(gpsPoints_Ind)
   gpsPoints_Ind$percentage <- as.numeric(gpsPoints_Ind[,1]) / nrow(eagle_subset) # number gps points per individual divided by tot number of gps points
   gpsPoints_Ind$weight <- abs(gpsPoints_Ind$percentage - 1) # absolute value (reversed weight of the individuals)
      
    # Calculate the mean, var, min, max for each variables 
    summary <- rbindlist(lapply(ind_list1, function(indi_data) {
            indi_summary <- indi_data %>%
                  select_if(is.numeric) %>%
                  summarise(across(everything(), .fns = list(
                        mean = ~ mean(., na.rm = TRUE),
                        var = ~ var(., na.rm = TRUE),
                        max = ~ max(., na.rm = TRUE),
                        min = ~ min(., na.rm = TRUE),
                        q1 = ~ quantile(., probs = 0.10, na.rm = TRUE)
                  )))
            return(indi_summary)
      }))
    
     summary <- merge(summary, gpsPoints_Ind[,2:4], by="individual.local.identifier.per.week", all.x=T)
      
    # Merge the resulting dataset with the vegetation dataset by individual.local.identifier.per.week 
      merged_data <- merge(summary, vegetation, by = "individual.local.identifier.per.week", all.x = TRUE)
      
    # Rbind the different files 
      PCA_eagle <- rbind(PCA_eagle, merged_data)
}




# STEP 5 : Dataset reduction ---------------------------------------------------------------------------------------------------------------------------- 
df_corr <- PCA_eagle

# Drop the 'percentage' column (already used for weighting, not needed in PCA)
df_corr <- as.data.frame(df_corr) %>% 
  dplyr::select(-percentage)

# Identify variables to scale (exclude ID and weight)
vars_to_scale <- names(df_corr)[
  sapply(df_corr, is.numeric) & 
    names(df_corr) != "weight"]


# FIRST FILTER : DROP CONSTANT COLUMNS
# Identify and drop null and columns with mostly 0, e.i., columns with 90% of values between -0.0001 & 0.0001
zero_only_vars <- vars_to_scale[
  sapply(df_corr[vars_to_scale], function(x) all(x == 0))]

sparse_prop_vars <- vars_to_scale[
  sapply(df_corr[vars_to_scale], function(x) {
    prop_near_0 <- mean(x > -0.0001 & x < 0.0001, na.rm = TRUE)
    prop_near_0 > 0.80
  })]

message("Null columns : ", paste(zero_only_vars, collapse = ", "))
message("Binary columns  : ", paste(sparse_prop_vars, collapse = ", "))
vars_to_remove <- union(zero_only_vars, sparse_prop_vars)

df_corr_clean <- df_corr %>%
  dplyr::select(-all_of(vars_to_remove))


# Scale numeric columns (except weight & ID), replace NA by 0
vars_to_scale_clean <- setdiff(
  names(df_corr_clean),
  c("individual.local.identifier.per.week", "weight")
) 

df_corr_clean[vars_to_scale_clean] <- scale(df_corr_clean[vars_to_scale_clean])
df_corr_clean[is.na(df_corr_clean)] <- 0 



# SECOND FILTER : WITHIN GROUP SELECTION 
# Define variable groups 
variable_groups <- list(
  elevation   = grep("^elevation_", names(df_corr_clean), value = TRUE),
  slope   = grep("^slope_", names(df_corr_clean), value = TRUE),
  TRI     = grep("^TRI_", names(df_corr_clean), value = TRUE),
  TPI     = grep("^TPI_", names(df_corr_clean), value = TRUE),
  ridge   = grep("^DIST_TO_RIDGE_", names(df_corr_clean), value = TRUE),
  height  = grep("^BUILDING_HEIGHT_", names(df_corr_clean), value = TRUE),
  heightdist  = grep("^DIST_TO_TALL_BUILDING_", names(df_corr_clean), value = TRUE),
  power   = grep("^DIST_POWER_LINE_", names(df_corr_clean), value = TRUE),
  roads   = grep("^DIST_ROADS_", names(df_corr_clean), value = TRUE),
  aerial  = grep("^DIST_AERIALWAYS_", names(df_corr_clean), value = TRUE),
  built   = grep("^NEAREST_BUILT_SURFACE_", names(df_corr_clean), value = TRUE),
  wind    = grep("^NEAREST_WIND_TURBINE_", names(df_corr_clean), value = TRUE)
  # add more groups as needed
)

# Compute correlation matrices within each group
group_cor_matrices <- lapply(variable_groups, function(var_names) {
  vars <- df_corr_clean[, var_names, drop = FALSE]
  if (ncol(vars) > 1) {
    cor(vars, use = "complete.obs")
  } else {
    NULL
  }
})

# Display the correlation matrices for inspection
for (group in names(group_cor_matrices)) {
  cat("\n=== Correlation matrix for:", group, "===\n")
  print(round(group_cor_matrices[[group]], 2))
} # remove one to three variables from each groups. 

# variables we decided to remove based on the correlation matrices
vars_to_remove_2 <- c(
  "elevation_mean", "elevation_var", "elevation_q1",
  "slope_q1", "slope_var",
  "TRI_brut_q1", "TRI_brut_var", 
  "TPI_small_brut_q1", "TPI_small_brut_var", 
  "DIST_TO_RIDGE_var", "DIST_TO_RIDGE_mean",
  "BUILDING_HEIGHT_var", "BUILDING_HEIGHT_max", 
  "DIST_TO_TALL_BUILDING_mean", "DIST_TO_TALL_BUILDING_q10", "DIST_TO_TALL_BUILDING_max",
  "DIST_POWER_LINE_q1", "DIST_POWER_LINE_var", "DIST_POWER_LINE_mean",
  "DIST_ROADS_q1", "DIST_ROADS_mean","DIST_ROADS_var",
  "DIST_AERIALWAYS_mean", "DIST_AERIALWAYS_q1", "DIST_AERIALWAYS_var", 
  "NEAREST_BUILT_SURFACE_var", "NEAREST_BUILT_SURFACE_q1", "NEAREST_BUILT_SURFACE_mean",
  "NEAREST_WIND_TURBINE_max", "NEAREST_WIND_TURBINE_q1", "NEAREST_WIND_TURBINE_mean", "NEAREST_WIND_TURBINE_var"
)

# Ensure variables exist before removing
vars_to_remove_2 <- intersect(vars_to_remove_2, names(df_corr_clean))
df_corr_sd_filter <- df_corr_clean %>% dplyr::select(-all_of(vars_to_remove_2))



# THIRD FILTER: BETWEEN GROUP OF VARIABLE COMPARISON AND SELECTION
# Prepare dataset by removing ID and weight
df_corr_rd_filter <- df_corr_sd_filter %>% 
  dplyr::select(-c(individual.local.identifier.per.week, weight))

# Compute the correlation matrix
cor.mat <- round(cor(df_corr_rd_filter, use = "pairwise.complete.obs"), 2)
cor.mat.long <- reshape2::melt(cor.mat)

# Filter for highly correlated pairs (above threshold)
thr <- 0.50
filt <- cor.mat.long %>%
  filter(Var1 != Var2, abs(value) >= thr)

# Select variables involved in high correlation
interesting_vars <- unique(c(filt$Var1, filt$Var2))

# Subset the wide matrix to only those variables
cor.mat.sub <- cor.mat[interesting_vars, interesting_vars]


         
# ------------------------------------------------------------------------------- Supplementary material : Figure 1. 

orig_order <- colnames(df_corr_clean)[colnames(df_corr_clean) %in% interesting_vars]
cor.mat.sub.long <- reshape2::melt(cor.mat.sub)
colnames(cor.mat.sub.long) <- c("Var1","Var2","value")
cor.mat.sub.long <- cor.mat.sub.long %>%
  mutate(
    Var1 = factor(Var1, levels = orig_order),
    Var2 = factor(Var2, levels = orig_order))

p <- ggplot(cor.mat.sub.long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) +  
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    values = scales::rescale(seq(-1, 1, length.out = 11)),
    limits = c(-1, 1),
    na.value = "grey90",
    guide = guide_colorbar(
      direction = "vertical",
      barwidth = unit(0.5, "cm"),
      barheight = unit(8, "cm"),
      title.position = "top"
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal(base_family = "Avenir") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
print(p)
# ------------------------------------------------------------------------------- (end) Supplementary material : Figure 1


         
# Remove the selected variables (correlation > 0.60 ± 0.03)
vars_to_remove_pca_pre = c(
  "DIST_TO_RIDGE_q1",  "slope_mean",
  "TPI_small_brut_mean", "slope_min", "slope_max", "DIST_TO_TALL_BUILDING_q1")
vars_to_remove_pca_pre <- intersect(vars_to_remove_pca_pre, names(df_corr_sd_filter))
df_var_pc <- df_corr_sd_filter %>% dplyr::select(-all_of(vars_to_remove_pca_pre))



# FOURTH FILTER : SUPPRESSION OF VARIABLES WITH VERY SMALL CONTRIBUTION TO THE PCs
# Keep the identifier and weights for PCA
identifiers <- df_var_pc$individual.local.identifier.per.week
weights     <- df_var_pc$weight

# Subset numeric-only variables for PCA
PCA_data2 <- df_var_pc %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-weight)  # weight is passed separately

# Run PCA using FactoMineR and identify variables with weighted contribution and communalities (cos2).
PCA_2 <- FactoMineR::PCA(
  PCA_data2,
  row.w      = weights,
  ncp        = ncol(PCA_data2),
  scale.unit = TRUE,
  graph      = FALSE
)

# Select the dimension that explain at least 50% of the variance of the dataset
eigvals <- PCA_2$eig[, "percentage of variance"] / 100
cumvar  <- cumsum(eigvals)
k <- which(cumvar >= 0.50)[1]  # with k representing the number of dimension that explain 50% of the variance

# Weighted contribution per dimensions
contrib_mat <- PCA_2$var$contrib[, 1:k]  
weighted_contrib <- rowSums( sweep(contrib_mat, 2, eigvals[1:k], `*`) )

# communalities
cos2_mat   <- PCA_2$var$cos2[, 1:k]
communalities <- rowSums(cos2_mat)

# Threshold
th_contrib <- mean(weighted_contrib)
th_comm    <- 0.50

low_contrib_vars <- names(weighted_contrib)[weighted_contrib < th_contrib]
low_comm_vars    <- names(communalities)[ communalities  < th_comm]
vars_to_drop <- intersect(low_contrib_vars, low_comm_vars)

# print
cat("Retained k axes :", k, " (", round(cumvar[k]*100,1), "% de variance )\n")
cat("Threshold weighted contribution :", round(th_contrib,4), "\n")
cat("Threshold communalities   :", th_comm, "\n\n")
cat("Variable with low weighted contribution :\n", 
    paste0("  • ", low_contrib_vars, collapse = "\n"), "\n\n")
cat("Variables with low communalities :\n", 
    paste0("  • ", low_comm_vars, collapse = "\n"), "\n\n")
cat("Variables that meet the two criteria :\n", 
    paste0("  • ", vars_to_drop, collapse = "\n"), "\n")

# Adapt your selection of the variable to remove based on both these statistical indicators and the ecological importance of the variables
manual_remove <- c(
  "TPI_small_brut_min",
  "TRI_brut_min",
  "DIST_TO_RIDGE_max",
  "DIST_TO_TALL_BUILDING_var","DIST_TO_TALL_BUILDING_min",
  "DIST_POWER_LINE_max",
  "DIST_ROADS_max",        
  "BUILDING_HEIGHT_mean", 
  "DIST_AERIALWAYS_max",
  "TRI_brut_max", 
  "NEAREST_BUILT_SURFACE_max", "Low-growing woody plants"
)

# Remove the variables
manual_remove <- intersect(manual_remove, names(df_var_pc))

df_pca_final <- df_var_pc %>%
  dplyr::select(-all_of(manual_remove))

# Set aside individual identifier and weight, ensure all variables are numerical for the PCA
identifiers2 <- df_pca_final$individual.local.identifier.per.week

weights2     <- df_pca_final$weight

PCA_data3 <- df_pca_final %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-weight)

# Run the PCA and plot the results
PCA_final <- FactoMineR::PCA(
  PCA_data3,
  row.w      = weights2,
  ncp        = ncol(PCA_data3),
  scale.unit = TRUE,
  graph      = FALSE
)

factoextra::fviz_pca_var(
  PCA_final,
  col.var       = "contrib",
  gradient.cols = c("#00AFBB","#E7B800","#FC4E07"),
  repel         = TRUE
)

#saveRDS(df_pca_final, "your_path/dataset.rds")

