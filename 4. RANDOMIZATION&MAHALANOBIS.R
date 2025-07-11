# RANDOMIZATION & MAHALANOBIS DISTANCE
# First aim of the script (step 1 to 4): calculating the Mahalanobis distance between each 7-day individual sequence of the dispersal period and evaluating individual consistency in the use of landscape features using a randomization test
# Second aim of the script (step 5 to 7): calculating the Mahalanobis distance between the pre-dispersal period and each of the 7-day sequence of the dispersal period.  
# Louise Faure, 09.07.2025.


# Library
library(FactoMineR)     
library(factoextra)     
library(dplyr)          
library(tidyr)          
library(ggplot2)        
library(reshape2)       
library(ggrepel)        
library(ggforce)        
library(tibble)
library(gridExtra)


# STEP 1 : LOAD PCA DATA & SELECTION OF NUMBER OF EIGENVALUES ------------------

# Load the dataset
eagle_PCA <- readRDS("/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ANALYSIS/database_for_PCA/results_F_d15.rds")

# Create a copy of the first dataset without the weight column as so that the column is not in the PCA
eagle_without_weight <- eagle_PCA[, !(colnames(eagle_PCA) %in% c("weight"))]

# Add row names to display the names of the birds in the PCA
PCA_data <- eagle_without_weight %>%
      column_to_rownames(var = "individual.local.identifier.per.week")

# Get the PCA (nb: the number of dimesions is equal to the number of variables)
PCA <- FactoMineR::PCA(PCA_data[,1:13], row.w = eagle_PCA$weight, ncp = 13)
summary(PCA)

# Selection of the variable based on the size (Jollife, 2002, Peres-Neto, Jackson ad Somers, 2005)
eigenvalues <- get_eigenvalue(PCA)

# Extract the PCA for the previously identified principal components 
PCA_dim_reduced <- PCA$ind$coord[, 1:7] # we keep only the 6 first dimension that explains 84,4 percent of the variation



#STEP 2 : CALCULATION OF THE MAHALANOBIS DISTANCE BETWEEN EACH WEEK OF THE DISPERSAL PERIOD
# Re order by date 
nameOrder <- sapply(strsplit(rownames(PCA_dim_reduced), " "), function(str){paste(str[-length(str)], collapse=" ")}) 
weekOrder <- sapply(strsplit(rownames(PCA_dim_reduced), " "), function(str){str[length(str)]}) 
weekOrder <- as.numeric(sub("D","",sub("F","0",weekOrder))) # change the D into numerical values and F to 0

PCA_dim_reduced_ordered <- PCA_dim_reduced[order(nameOrder, weekOrder),]

# suppress all the row with individual name that end by "F" to suppress the fledging period (pre-dispersal)
PCA_dim_reduced_ordered_withoutF <- PCA_dim_reduced_ordered[!grepl(" F$", rownames(PCA_dim_reduced_ordered)), ]

# reorder by names
nameSplit <- sapply(strsplit(rownames(PCA_dim_reduced_ordered_withoutF), " "), function(str){paste(str[-length(str)], collapse=" ")})
matList <- split(as.data.frame(PCA_dim_reduced_ordered_withoutF), nameSplit)

# Calculate distances between each week of the dispersal period (from week 1 to week 2, from week 2 to week 3, etc ...)
distMatList <- lapply(matList, function(ind) {
      if (nrow(ind) < 2) return(NULL)
      distances <- sapply(1:(nrow(ind) - 1), function(i) {
            sqrt(sum((ind[i, ] - ind[i + 1, ])^2))
      })
      names(distances) <- paste0("W", 1:(length(distances)), "_to_W", 2:(length(distances) + 1))
      return(distances)
})

# Convert list of distance vectors into a long-format dataframe
combinedDistMat <- do.call(rbind, lapply(names(distMatList), function(ind_name) {
      dist_vec <- distMatList[[ind_name]]
      if (is.null(dist_vec)) return(NULL)
      
      data.frame(
            individual = ind_name,
            week_pair = names(dist_vec),
            distance = as.numeric(dist_vec)
      )
}))


# calculate the mean, max and sum of distance per individuals 
all_obs_dists <- combinedDistMat$distance   

# Compute the three pooled statistics 
observed_mean_dist <- mean(all_obs_dists)  
observed_sum_dist  <- sum(all_obs_dists)    
observed_max_dist  <- max(all_obs_dists)



# STEP 3 : GENERATE THE NULL DISTRIBUTION --------------------------------------
set.seed(123)  # to ensure reproducibility by ensuring fix random seed
n_permutation <- 1000 # number of permutations for the Monte Carlo Permutation Procedure (MCPP)

# Flatten all PCA windows into a single matrix with individual and week info
df_all <- as.data.frame(PCA_dim_reduced_ordered)
df_all <- df_all[!grepl(" F\\s*$", rownames(df_all)), ]
df_all$individual <- sapply(strsplit(rownames(df_all), " "), function(str) {
      paste(str[-length(str)], collapse = " ")
})
df_all$week <- as.numeric(gsub("D", "", sapply(strsplit(rownames(df_all), " "), function(str) str[length(str)])))


# Keep only individuals with >= 2 weeks of data (because individuals with less than 2 weeks cannot contribute to distance)
counts <- table(df_all$individual)
valid_inds <- names(counts)[counts >= 2]
df_all <- df_all[df_all$individual %in% valid_inds, ]

# controle number of individuals
length(unique(df_all$individual))
table(df_all$individual)

# prepare fake individuals based on the number of 'real' individuals
ind_sizes   <- table(df_all$individual)       
ind_labels  <- names(ind_sizes)              
n_rows_tot  <- nrow(df_all)
print(ind_sizes)
print(n_rows_tot)

# compute mean distances for each permutation
null_stats <- data.frame(mean = numeric(n_permutation),
                         sum  = numeric(n_permutation),
                         max  = numeric(n_permutation)) # vector to store mean, max and sum of distances

# loop
for (perm in 1:n_permutation) {
      shuffled_idx <- sample.int(n_rows_tot)
      # Assign rows to fake individuals so that each fake bird gets the same number of rows as its real counterpart
      fake_id_vec <- rep(names(ind_sizes), times = as.integer(ind_sizes))   
      df_perm     <- df_all[shuffled_idx, ]               
      df_perm$fake_id <- fake_id_vec                      
      
      matList <- lapply(split(df_perm, df_perm$fake_id), function(sub) {
            as.matrix(sub[ , 1:7, drop = FALSE])      # keep 7 PC dims
      })
      # Compute week-to-week distances
      distMatList <- lapply(matList, function(ind) {
            if (nrow(ind) < 2) return(NULL)
            sapply(seq_len(nrow(ind) - 1), function(i) {
                  sqrt(sum((ind[i, ] - ind[i + 1, ])^2))
            })
      })
      # Pool all distances from this permutation
      all_dists <- unlist(distMatList, use.names = FALSE)
      
      # Store summary stats
      null_stats$mean[perm] <- mean(all_dists)
      null_stats$sum[perm]  <- sum(all_dists)
      null_stats$max[perm]  <- max(all_dists)
}
   

# STEP 4 : P-VALUE CALCULATION -------------------------------------------------
p_value_mean <- (sum(null_stats$mean >= observed_mean_dist) + 1) / (n_perm + 1)
p_value_sum  <- (sum(null_stats$sum  >= observed_sum_dist)  + 1) / (n_perm + 1)
p_value_max  <- (sum(null_stats$max  >= observed_max_dist)  + 1) / (n_perm + 1)

# print results
cat("P-values:\n")
cat("  Mean distance: ", format(p_value_mean, digits = 10), "\n")
cat("  Sum distance: ", format(p_value_mean, digits = 10), "\n")
cat("  Max distance: ", format(p_value_mean, digits = 10), "\n")

# display results under a table format
null_stats_long <- null_stats %>%
      pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

observed_df <- data.frame(
      Metric = c("mean", "sum", "max"),
      Observed = c(observed_mean_dist, observed_sum_dist, observed_max_dist)
)

summary_table <- data.frame(
      Metric = c("Mean", "Sum", "Max"),
      Observed = round(c(observed_mean_dist, observed_sum_dist, observed_max_dist), 5),
      Null_Mean = round(c(mean(null_stats$mean), mean(null_stats$sum), mean(null_stats$max)), 5),
      P_value = signif(c(p_value_mean, p_value_sum, p_value_max), 5)
)
grid.table(summary_table)

# three histogram for the sum, mean and median with the observed and randomely simulated data
combined_plot_df <- null_stats_long %>%
      mutate(Source = "Simulated") %>%
      bind_rows(
            observed_df %>% rename(Value = Observed) %>%
                  mutate(Source = "Observed"))
ggplot(combined_plot_df, aes(x = Value, fill = Source)) +
      geom_histogram(data = subset(combined_plot_df, Source == "Simulated"),
                     bins = 40, colour = "white", alpha = 0.8) +
      geom_vline(data = subset(combined_plot_df, Source == "Observed"),
                 aes(xintercept = Value),
                 linetype = "dashed", colour = "red", linewidth = 1) +
      facet_wrap(~Metric, scales = "free", ncol = 3) +
      labs(
            title = "Comparison of Simulated and Observed Distributions",
            x = "Value", y = "Count"
      ) +
      scale_fill_manual(values = c("Simulated" = "steelblue", "Observed" = "red")) +
      theme_minimal(base_size = 14)




# STEP 5 : CALCULATION OF THE MAHALANOBIS DISTANCE BETWEEN PRE-DISPERSAL AND 7-DAYS SEQUENCE OF THE DISPERSAL PERIOD
# reorder by names
name <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "), function(str){paste(str[-length(str)], collapse=" ")})
List <- split(as.data.frame(PCA_dim_reduced_ordered), name)

# Calculate distances 
MahalanobisList <- lapply(List, function(ind){
      MahalanobisDist <- dist(ind, method = "euclidean", diag=T, upper=T)
      # MahalanobisDist <- ind[-1,]
      # ind[-nrow(ind),]
      return(MahalanobisDist)})

# Convert the title of the colomn to 0, 1, 2, 3 ...
MahalanobisList <- lapply(MahalanobisList, function(distance) {
      distance <- as.matrix(distance)  # Convert distance object to matrix
      colnames(distance) <- seq_len(ncol(distance)) - 1
      return(distance)})


# Keep the first row of each matrix
MahalanobisList <- lapply(MahalanobisList, function(distance) {
      First_row <- distance[1, , drop = FALSE]
      return(First_row)})

# Find the maximum number of columns
max_col <- max(sapply(MahalanobisList, ncol))

# Add missing columns with NA values to align the matrices
alignedMahalanobisList <- lapply(MahalanobisList, function(distance) {
      Num_missing_cols <- max_col - ncol(distance)
      if (Num_missing_cols > 0) {
            distance <- cbind(distance, matrix(NA, nrow = nrow(distance), ncol = Num_missing_cols))
      }
      return(distance)
})

# Bind each row of alignedDistances
combinedDistances <- do.call(rbind, alignedMahalanobisList)



# Step 6 : CLEAN DATA BEFORE VISUALISATION -------------------------------------
# Remove the rows for which there are NA in the column called 1
combinedDistances_df <- as.data.frame(combinedDistances)
combinedDistances_df <- combinedDistances_df[!is.na(combinedDistances_df$`1`), ]
combinedDistances_df <- combinedDistances_df %>%
      rownames_to_column(var = "Individuals")
combinedDistances_df$Individuals <- gsub(" F", "", combinedDistances_df$Individuals)
combinedDistances_df$Individuals <- trimws(combinedDistances_df$Individuals)

# Add death information
death <- read.csv("/Volumes/My Passport/Louise/memoire_aigle/wetransfer_death_document-csv_2023-06-30_1528/death_document.csv")
death$individual.local.identifier <- gsub("\\(.*\\)", "", death$individual.local.identifier)
death$individual.local.identifier <- trimws(death$individual.local.identifier)
death$individual.local.identifier <- trimws(death$individual.local.identifier)

# Filter death file to keep rows with matching names to combinedDistMat_df
death_filtered <- death[death$individual.local.identifier %in% combinedDistances_df$Individuals, ]

# Controle the matching names
#missing_names <- setdiff(combinedDistances_df$Individuals, death_filtered$individual.local.identifier)
#print(missing_names)

# Merge the dataframe based on matching names
merged_df <- merge(combinedDistances_df, death, by.x = "Individuals", by.y = "individual.local.identifier", all.x = TRUE)

# Suppress rows with lack of information
na_rows <- apply(merged_df[, 2:16], 1, function(x) any(is.na(x)))

# Subset the dataframe based on the condition
df <- merged_df[!(na_rows & merged_df$loss_cause != "death"), ]

# Create rownames
df<- df[, !(names(subset_df) %in% c("cause_death", "loss_location", "loss_cause"))]



# STEP 7 : CLUSTERING & DATA VISUALISATION -------------------------------------
# Reshape the data from wide to long format
rownames(subset_df) <- NULL

# Create a new column with the variance 
subset_df_1 <- subset_df %>%
      mutate_at(vars(2:17), as.numeric)

subset_df_1$variance <- apply(subset_df[, 2:17], 1, var, na.rm = TRUE)

# Create a new column with the maximum value 
subset_df_1$maximum <- apply(subset_df[,  2:17], 1, max, na.rm = TRUE)

# Create three clusters for forest and mountain eagle
subset_df_excluded <- subset_df_1[subset_df_1$Individuals != "Almen18", ]

# Group individual using a k-means clustering on individuals Mahalanobis distance variance and maximum
subset_df_excluded <- subset_df_excluded %>%
      filter(!is.na(variance), !is.na(maximum), is.finite(variance), is.finite(maximum))
cluster_variance <- kmeans(subset_df_excluded[, c("variance", "maximum")], centers = 3, nstart = 25)
subset_df_excluded$cluster <- cluster_variance$cluster

# add landscape utilization strategy of the pre-dispersal period (F)
lus_fledging <- readRDS("/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ARTICLE/supplementary_data/LUS_fledging.rds")

# Plot the result with x = variance, y = maximum
colour_mahalanobis <- c("forest" = "darkcyan", mountain = "salmon", valley = "azure3")

ggplot(mahalanobis, aes(x = variance, y = maximum, color = type, shape = as.factor(cluster))) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = Individuals, col = "dark"), size = 3, box.padding = 0.5, segment.color = "transparent") +
      geom_mark_ellipse(aes(x = variance, y = maximum, color = as.factor(cluster)),
                        expand = unit(2, "mm"), radius = unit(2, "mm"), con.colour = "black",
                        con.size = 0.1, con.type = "elbow", con.linetype = 1,
                        con.border = "one", con.cap = unit(0.7, "mm"), label.buffer = unit(8, "mm")) +
      labs(x = "Variance of the Mahalanobis Distance", y = "Maximum of the Mahalanobis Distance") +
      theme(text = element_text(family = "Avenir"),
            axis.title = element_text(face = "plain"),
            axis.text = element_text(size = 10),
            axis.title.x = element_text(size = 12),
            legend.text = element_text(size = 7),
            legend.key.size = unit(1, "cm"),  # Adjust the size of the legend key
            legend.key.width = unit(1, "cm"),  # Adjust the width of the legend key
            legend.text.align = 0.4,
            plot.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "cm"),
            legend.margin = margin(4, -1, 0, 1, unit = "cm"),
            legend.title = " ") +  
      theme_gray()  


