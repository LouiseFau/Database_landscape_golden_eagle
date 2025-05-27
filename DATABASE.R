# JUVENILE GOLDEN EAGLE DATABASE FOR TERRESTRIAL ANALYSIS OF THE LANDSCAPE - LOCATION POINT SCALE
# This R script is made to ; I. Select the eagle GPS location after emigration date (sampling part), II. calculate the distance between each GPS location points and flight altitude, III. extract from each previously constructed layers (see script 1 and 2), the values of each pixels under the localisation of a golden eagle


#library
library(terra)
library(lubridate)
library(dplyr)
library(purrr)
library(sf)
library(move)

# STEP 0 : load the required data ----------------------------------------------

# Emigration date
emig_dates <- readRDS("emigration_dates.rds")

# Update the date for Matsch19 (Swiss Ornithological Institute estimation)
emig_dates$emigration_dt <- as.POSIXct(emig_dates$emigration_dt)
emig_dates$emigration_dt[emig_dates$individual.local.identifier == "Matsch19 (eobs 7035)"] <- as.POSIXct("2020-03-08  09:00:08")
# Update the date for Kastelbell19 (based on Hester's work)
emig_dates$emigration_dt[emig_dates$individual.local.identifier == "Kastelbell19 (eobs 7034)"] <- as.POSIXct("2020-02-04 15:40:14")


# Geoid and raster
geo <- terra::rast("/home/louise/Desktop/Topography/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("/home/louise/Desktop/Topography/data/pretraitements/Region-Alpes-Dem/DEM_CRS_Lonlat/dem_lonlat.tif") 
utm_crs <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"

# Topographic variables
elevation_25 <- rast("/home/louise/Desktop/Topography/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
slope_25 <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/slope/slope.tif")
northness_25 <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/northness/northness.tif")
eastness_25 <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/eastness/eastness.tif")
TRI_25 <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/TRI/TRI.sdat")
TPI_25_small_scale <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/TPI/TPI.sdat")
INDEX_DIST_TO_RIDGE <- rast("/home/louise/Desktop/Topography/data/pretraitements/DEM_25m/Distance_to_Ridge/map_distance_ridge.tif")

# Human variables
CLC_ALPS <- terra::rast("/home/louise/Desktop/Topography/data_landuse/Data_CLC/CLC_Alps/CLC_longlat/CLC_longlat.tif")
BUILDING_HEIGHT <- terra::rast("/home/louise/Desktop/Topography/data_landuse/Building/Building_height/BUILDING_HIGHT/BUILDING_longlat/building_long_lat.tif")
DIST_TO_10M_BUILDING <- rast("/home/louise/Desktop/Topography/data_landuse/Building/Distance_to_hight_building/10M_building/10MBUILDING_DIST_Longlat/DIST_TO_10MBUILDING_longlat.tif")
DIST_TO_POWER_LINE <- rast("/home/louise/Desktop/Topography/data_landuse/raster_map_distance/raster_distance_power_line/DIST_POWERLINE_longlat/dist_powerline_longlat.tif") 
DIST_TO_ROADS <- rast("/home/louise/Desktop/Topography/data_landuse/raster_map_distance/raster_distance_roads/dist_road_longlat/dist_road_longlat.tif") 
DIST_TO_AERIALWAYS <- rast("/home/louise/Desktop/Topography/data_landuse/raster_map_distance/raster_distance_aerialways/DIST_AERIALWAYS_longlat/DIST_AERIALWAYS_longlat.tif") 
SETTLEMENT_DENSITY <- rast("/home/louise/Desktop/Topography/data_landuse/settlment_density/settlment_density/density_of_built_surface_3.tif")
DIST_TO_NEAREST_BUILT_SURFACE <- rast("/home/louise/Desktop/Topography/data_landuse/raster_map_distance/dist_settlment/distance_moins_NA_Zone.tif")
DIST_TO_WIND_TURBINE <- rast("/home/louise/Desktop/Topography/data_landuse/raster_map_distance/distance_wind_turbine/raster_distance_wind_turbine_CRS_CLC.tif")

# list all the file
list_file <- list.files("/home/louise/Desktop/Topography/gps_for_Louise", pattern = ".rds", full.names = TRUE)

# Create an empty data.frame that will stored all the extracted values from the golden eagle locations
data_eagle <- data.frame()


# Loop through the files
for (ind in seq_along(list_file)) {
  
  print(ind)
  # STEP 1 : load and read the data in R ---------------------------------------
  data <- readRDS(list_file[ind])
  
  # STEP 2 : Sampling : assign life stages: filter post-emigration -------------
  # Extract the individual identifier columns, lat, long and height above ellipsoid and suppress all the other columns
  data <- data[, c("timestamp", "location.long", "location.lat", "height.above.ellipsoid","individual.local.identifier" )]
  
  # Filter GPS location after the emigration date and conserve only the location after emigration date
  d <- emig_dates %>%
    filter(emig_dates$individual.local.identifier == unique(data$individual.local.identifier))
  
  d$emigration_dt <- as.POSIXct(d$emigration_dt, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Check if there is a match in emig_dates
  if (nrow(d) > 0) {
    data_emig <- data %>%
    mutate(stage = as.logical(data$timestamp >= d$emigration_dt)) # create a new columns with the stage, i.e, dispersed or not dispersed
    data_emig <- data_emig[data_emig$stage == "TRUE", ] # remove all the GPS location before emigration date
    
    # Create a new columns with a numeric values for the number of days after emigration
    data_emig$days_since_emig <- as.numeric(difftime(data_emig$timestamp, d$emigration_dt, units = "days"))
    
    # Filter for one week after emigration date
    data_sample <- data_emig %>%
      filter(days_since_emig > 105 & days_since_emig <= 112)}
      #filter(days_since_emig <= 7)}
  
  # As some individual have a dispersal date but no GPS location after dispersal, we should not take into account these individuals
  if (nrow(data_sample) > 30) {
    
    # Step 3a : Sampling of the GPS location points ----------------------------
    # Remove all double timestamp
    rows_to_delete <- unlist(sapply(move::getDuplicatedTimestamps(x = as.factor(data_sample$individual.local.identifier), timestamps = as.POSIXct(data_sample$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")), "[", -1)) #get all but the first row of each set of duplicate rows

    if(length(rows_to_delete) > 0){
      data_timestamp <- data_sample[-rows_to_delete,]
      } else {
        data_timestamp <- data_sample
        }
    
    # Sampling every 15 minutes
    eagle_15m <- data_timestamp %>%
      mutate(timestamp_rounded = floor_date(timestamp, "15 minutes")) %>% 
      group_by(timestamp_rounded) %>%
      slice(1)
    
    eagle_15m <- as.data.frame(eagle_15m)
    
    # Create a move object for the individual
    data_move <- move(x = eagle_15m$location.long, y = eagle_15m$location.lat,
                    time = eagle_15m$timestamp,
                    proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                    data = eagle_15m, animal = eagle_15m$individual.local.identifier)
    
    # Calculate step distances
    data_move@data$step_dist <- c(NA, distance(data_move))
    #data_move@data$timeLag_m <- c(NA, timeLag(data_move, units="mins"))
    
    # Select GPS location points every 150 meters
    while(any(data_move$step_dist[!is.na(data_move$step_dist)] < 150)){
      toRemove <- min(which(data_move$step_dist < 150))
      data_move <- data_move[-toRemove,]
      data_move$step_dist <- c(NA, distance(data_move))
    }
    
    # Select GPS location points every 5 minutes
    #move_5min <- as.data.frame(data_move)
    #indAmt <- mk_track(tbl=move_5min, all_cols=T,
                       #.x=location.long, .y=location.lat, crs = st_crs(4326),
                       #.t=timestamp, order_by_ts = T, check_duplicates = T)
    #indAmt <- track_resample(indAmt, rate = minutes(15), tolerance = minutes(2), start = 1)
    #move_5min <- as_move(indAmt)


    # Step 3c : Calculate height above the ground ------------------------------
    # Convert move object into a data.frame 
    data_df <- as.data.frame(data_move)
    
    # Convert data.frame to spatial object
    data_df_spatial <- st_as_sf(data_df, 
                              coords = c("location.long", "location.lat"),
                              crs=crs("+proj=longlat +datum=WGS84 +no_defs"),
                              remove = F)
    
    # Adjust the CRS of the spatial object so that it match the CRS of the dem and geoid raster
    data_df_spatial_sf <- st_transform(data_df_spatial, crs(dem))
    
    # Extract elevation values
    data_df_spatial_sf$dem <- terra::extract(x = dem, y = data_df_spatial_sf[,c("location.long","location.lat")], method = "bilinear")[,2] 
    data_df_spatial_sf$geoid <- terra::extract(x = geo, y = data_df_spatial_sf[,c("location.long","location.lat")], method = "bilinear")[,2] 
    
    # Calculate flight height as height above mean sea level - dem. height above sea level is height above ellipsoid - geoid
    data_flight_height <- data_df_spatial_sf %>% 
      mutate(height_msl = height.above.ellipsoid - geoid) %>% 
      mutate(height_ground = height_msl - dem)
    
    # STEP 4 : extract values from alpines topographic layers ---------------------
    # Transform data_df into a spatial object
    individual_df_sf <- st_as_sf(data_flight_height, 
                               coords = c("location.long", "location.lat"),
                               crs=crs("+proj=longlat +datum=WGS84 +no_defs"),
                               remove = F)
    
    # Adjust the CRS
    individual_df_sf <- st_transform(individual_df_sf, crs(slope_25))
    
    
    # Extract for each layers the values of the pixel under the GPS location
    elevation_values <- terra::extract(elevation_25, individual_df_sf)
    slope_values <- terra::extract(slope_25, individual_df_sf)
    northness_values <- terra::extract(northness_25, individual_df_sf)
    eastness_values <- terra::extract(eastness_25, individual_df_sf)
    TPI_small_brut_values <- terra::extract(TPI_25_small_scale, individual_df_sf)
    TRI_brut_values <- terra::extract(TRI_25, individual_df_sf)
    DIST_TO_RIDGE_values <- terra::extract(INDEX_DIST_TO_RIDGE, individual_df_sf)
    
    # Create new columns for each variables 
    data_topo <- cbind(individual_df_sf, data.frame(
      elevation = elevation_values$eu_dem_v11_E40N20, 
      slope = slope_values$slope, 
      northness = northness_values$aspect,
      eastness = eastness_values$aspect, 
      TPI_small_brut = TPI_small_brut_values$TPI, 
      TRI_brut = TRI_brut_values$TRI, 
      DIST_TO_RIDGE = DIST_TO_RIDGE_values$distance_ridge))
    
    
    # Step 4 : Extract human based features of the golden eagle landscape --------
    # Adjust the CRS 
    individual_layer <- st_transform(data_topo, crs(CLC_ALPS))
    
    # Extract for each layers the values of the pixel under the GPS location
    CLC_values <- terra::extract(CLC_ALPS, individual_layer)
    BUILDING_HEIGHT_values <- terra::extract(BUILDING_HEIGHT, individual_layer)
    DIST_TO_10M_BUILDING_values <- terra::extract(DIST_TO_10M_BUILDING, individual_layer)
    DIST_TO_POWER_LINE_values <- terra::extract(DIST_TO_POWER_LINE, individual_layer)
    DIST_TO_ROADS_values <- terra::extract(DIST_TO_ROADS, individual_layer)
    DIST_TO_AERIALWAYS_values <- terra::extract(DIST_TO_AERIALWAYS, individual_layer)
    SETTLEMENT_DENSITY_values <- terra::extract(SETTLEMENT_DENSITY, individual_layer)
    DIST_TO_NEAREST_BUILT_SURFACE_values <- terra::extract(DIST_TO_NEAREST_BUILT_SURFACE, individual_layer)
    DIST_TO_WIND_TURBINE_values <- terra::extract(DIST_TO_WIND_TURBINE, individual_layer)
    
    # Create new columns for each variables
    eagle_human <- cbind(individual_layer, data.frame(
      CLC = CLC_values$Class_name,
      BUILDING_HEIGHT = BUILDING_HEIGHT_values$GHS_BUILT_H_ANBH_E2018_GLOBE_R2022A_54009_100_V1_0_R4_C20, 
      DIST_TO_TALL_BUILDING = DIST_TO_10M_BUILDING_values$GHS_BUILT_H_ANBH_E2018_GLOBE_R2022A_54009_100_V1_0_R4_C20, 
      DIST_POWER_LINE = DIST_TO_POWER_LINE_values$map_power_line, 
      DIST_ROADS = DIST_TO_ROADS_values$raster_distance_to_roads2, 
      DIST_AERIALWAYS = DIST_TO_AERIALWAYS_values$raster_distance_to_aerialways, 
      DENSITY = SETTLEMENT_DENSITY_values$focal_sum,
      NEAREST_BUILT_SURFACE = DIST_TO_NEAREST_BUILT_SURFACE_values$distance_moins_NA_Zone,
      NEAREST_WIND_TURBINE = DIST_TO_WIND_TURBINE_values$raster_distance_wind_turbine))
    
    # Rbind the results to add rows for each individuals 
    data_eagle <- rbind(data_eagle, eagle_human)
    
    # Clean the database to suppress the undesired columns
    data_eagle_landscape <- data_eagle[, !(names(data_eagle) %in% c("coords.x1", "coords.x2", "sensor", "timestamps", "optional", "dem", "geoid", "height_msl", "step_dist", "height.above.ellipsoid", "days_since_emig", "timestamp_rounded", "stage"))]
    
    # Rename the columns names
    names_new <- sub("\\s*\\(.*\\)", " D16", data_eagle_landscape$individual.local.identifier)
    data_eagle_landscape$individual.local.identifier.per.week <- names_new
    
    # Get ride of the eagle that fly outside the alpin area
    data_eagle_landscape <- data_eagle_landscape[data_eagle_landscape$individual.local.identifier!="MonteAdone22 (eobs 10540)",]
  }
}

# Save the data.frame
saveRDS(data_eagle_landscape, file = "name_of_your_file.rds")