# Between and within individual variation in spatial behavior
# This script is aims at: I. identifying variation in landscape feature use across individuals using a PCA (Step 1 to 2) and a clusterisation method (Step 3 to 4).
# II. Testing individual consistency in landscape use by comparing the observed mean distance between consecutive individual 7-day sequence point to randomely re-assigned coordinates.
# III. Testing the dependance between the behavior during pre-dispersal and dispersal by comparing the observed distance value between each week of the dispersal phase and the pre-dispersal point to the randomely
# re assigned location point.
# Louise FAURE, MPI, 23/06/2023

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
library(NbClust)
library(mclust)

eagle_PCA <- readRDS("path_to_your_data/resultsFtoD15.rds")

# STEP 1 : SELECTION OF PCA DIMENSION --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copy of the first dataset without the weight column
eagle_without_weight <- eagle_PCA[, !(colnames(eagle_PCA) %in% c("weight"))]

# Add row names to display the names of the birds in the PCA
PCA_data <- eagle_without_weight %>%
  column_to_rownames(var = "individual.local.identifier.per.week")

# PCA loadings (nb: the number of dimensions is equal to the number of variables)
PCA <- FactoMineR::PCA(PCA_data, row.w = eagle_PCA$weight, ncp = 12)
summary(PCA)

# ------------------------------------------------------------------------------ Supplementary material : F2. A and B
eig_df <- as.data.frame(get_eigenvalue(PCA))
eig_df$dim_num <- 1:nrow(eig_df)
eig_df$percentage <- round(eig_df$variance.percent, 1)


custom_blues <- rev(colorRampPalette(brewer.pal(9, "Blues"))(nrow(eig_df)))
eig <- ggplot(eig_df, aes(x = dim_num, y = variance.percent)) +
      geom_col(fill = custom_blues, color = "grey30", width = 0.8) +
      geom_line(color = "black", linewidth = 0.8) +
      geom_point(color = "black", size = 2) +
      geom_text(aes(label = percentage), vjust = -0.8, hjust = 0.1, size = 4, color = "black") +
      scale_x_continuous(breaks = eig_df$dim_num) +
      labs(x = "Principal components",
            y = "Variance (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
            axis.title = element_text(size = 16),
            panel.grid.major.y = element_line(color = "gray85"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
      ) # 4 first dimensions selected
print(eig) 

plot.new()
pal <- colorRampPalette(c("#F7FBFF", "#4292C6", "#08519C"))(100)
corrplot(
  PCA$var$contrib[, 1:4],
  is.corr     = FALSE,
  method      = "color",
  col         = pal,
  addgrid.col = "grey85",
  addCoef.col = "black",
  number.cex  = 0.7,  
  tl.col      = "black",
  tl.srt      = 45,
  tl.cex      = 0.8,
  cl.pos      = "n"
)
# ------------------------------------------------------------------------------ (end) Supplementary material : F2. A and B


# Selection of the variable and individuals points based on the size (Jollife, 2002, Peres-Neto, Jackson ad Somers, 2005)
PCA_ind_reduced <- PCA$ind$coord[, 1:4] # variables
PCA_var_reduced <- PCA$var$coord[, 1:4] # individual 7-day sequence point



# STEP 2 : PCA & LANDSCAPE UTILIZATION STRATEGIES---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Cluster identification
set.seed(123)
cluster_ind <- kmeans(PCA_ind_reduced, centers = 3, nstart = 25)
grp_ind    <- as.factor(cluster_ind$cluster) 

cluster_ind$size

cluster_var <- kmeans(PCA_var_reduced, centers = 3, nstart = 25)
grp_var <- as.factor(cluster_var$cluster)

# Visualization
cluster_palette <- c("1" = "darkcyan", "2" = "salmon", "3" = "azure3")
biplot <- fviz_pca_biplot(
      PCA,
      col.ind      = grp_ind,           
      col.var = grp_var,
      palette      = cluster_palette, 
      addEllipses  = TRUE,
      repel        = TRUE,              
      label        = "var",             
      ellipse.type = "t",
      ellipse.level= 0.95,
      arrow.col    = "black",
      arrow.size   = 0.3,
      legend.title = " "
) + theme_minimal()
print(biplot)

# Re-adjust the cluster palette and visualization
cluster_palette <- c("1" = "darkcyan", "2" = "azure3", "3" = "salmon")
biplot <- fviz_pca_biplot(
      PCA,
      col.ind      = grp_ind,           
      col.var      = "black",           
      palette      = cluster_palette, 
      addEllipses  = TRUE,
      repel        = TRUE,              
      label        = "var",             
      ellipse.type = "t",
      ellipse.level= 0.95,
      arrow.col    = "black",
      arrow.size   = 0.3,
      legend.title = " "
) + theme_minimal()
print(biplot)



# STEP 4 : CLUSTER STRENGHT INDICATORS -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Small WCSS and large BCSS
WCSS <- cluster_ind$tot.withinss # Within Cluster Sum of Squares (= total amount of variation inside the clusters once the data is partitioned)
BCSS <- cluster_ind$betweenss # Between Cluster Sum of Squares (= amount of variation explained by the separation between cluster centers)

# Set the range of cluster numbers to evaluate
min_clusters <- 2
max_clusters <- 10

# Initialize vectors to store WCSS and BCSS values for each cluster number
wcss <- numeric(max_clusters - min_clusters + 1)
bcss <- numeric(max_clusters - min_clusters + 1)

# Loop through each cluster number
for (k in min_clusters:max_clusters) {
  # Perform k-means clustering
  kmeans_result <- kmeans(PCA_ind_reduced, centers = k, nstart = 25)
  
  # Calculate WCSS
  wcss[k - min_clusters + 1] <- kmeans_result$tot.withinss
  
  # Calculate BCSS
  bcss[k - min_clusters + 1] <- kmeans_result$betweenss
}

# Print WCSS and BCSS for each cluster number
for (k in min_clusters:max_clusters) {
  cat("Cluster number:", k, "\tWCSS:", wcss[k - min_clusters + 1], "\tBCSS:", bcss[k - min_clusters + 1], "\n")
}

# Create a data frame with cluster numbers, WCSS, and BCSS
cluster_results <- data.frame(
  ClusterNumber = min_clusters:max_clusters,
  WCSS = wcss,
  BCSS = bcss
)

# standardize
cluster_results$WCSS <- scale(cluster_results$WCSS)[,1]
cluster_results$BCSS <- scale(cluster_results$BCSS)[,1]

cluster_results$diff <- abs(cluster_results$BCSS - cluster_results$WCSS)
cluster_results$alpha <- scales::rescale(max(cluster_results$diff) - cluster_results$diff, to = c(0.1, 0.8))

# ------------------------------------------------------------------------------ Supplementary material : Figure 3
ribbon_df <- bind_rows(
      cluster_results %>% dplyr::select(x = ClusterNumber, y = WCSS, alpha),
      cluster_results %>% arrange(desc(ClusterNumber)) %>%
        dplyr::select(x = ClusterNumber, y = BCSS, alpha))

ribbon_layers <- lapply(1:(nrow(cluster_results) - 1), function(i) {
      geom_ribbon(
            data = cluster_results[i:(i+1), ],
            aes(x = ClusterNumber, ymin = WCSS, ymax = BCSS),
            fill = "#e41a1c", alpha = cluster_results$alpha[i],
            inherit.aes = FALSE)})

ggplot(cluster_results, aes(x = ClusterNumber)) +
      c(ribbon_layers,
        list(geom_line(aes(y = WCSS, linetype = "WCSS", color = "Metric"), linewidth = 0.6),
              geom_line(aes(y = BCSS, linetype = "BCSS", color = "Metric"), linewidth = 0.6),
              geom_point(aes(y = WCSS), color = "black", size = 2),
              geom_point(aes(y = BCSS), color = "black", size = 2),
              scale_color_manual(values = c("Metric" = "black")),
              scale_linetype_manual(values = c("WCSS" = "dashed", "BCSS" = "solid")),
              labs(x = "Number of clusters", y = "Standardized sum of squares distances"),
              theme_minimal(base_size = 14) +
                    theme(
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 12),
                          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                          legend.title = element_text(size = 14),
                          legend.text = element_text(size = 12),
                    )
        )
      ) -> full_plot
full_plot # optimal cluster between k = 2 and k = 5
# ------------------------------------------------------------------------------ (end) Supplementary material : Figure 3

# 30 indices (Charrad et al., 2014)
df_scaled <- scale(PCA_ind_reduced)
NbClust(PCA_ind_reduced, diss = NULL, distance = "manhattan", min.nc = 2, max.nc = 6, method = "kmeans") # optimal number of cluster : 2, 6 or 3. 


# Bayesian Information Criterion
d_clust <- Mclust(PCA_ind_reduced, G = 1:6, 
                  modelNames = mclust.options("emModelNames"))
d_clust$BIC # optimal number of clusters : 5, 3 or 6. 




# STEP 5 : INDIVIDUAL LANDSCAPE UTILIZATION STRATEGY CONSISTENCY ----------------------------------------------------------------------------------------------------------------------------------------------------
# Re order by date 
nameOrder <- sapply(strsplit(rownames(PCA_ind_reduced), " "), function(str){paste(str[-length(str)], collapse=" ")}) # first part of the code divide the matrix per space and get rid of the last part with space
weekOrder <- sapply(strsplit(rownames(PCA_ind_reduced), " "), function(str){str[length(str)]}) # merge the names that have a space
weekOrder <- as.numeric(sub("D","",sub("F","0",weekOrder))) # change the D into numerical values and F to 0

PCA_ind_reduced_ordered <- PCA_ind_reduced[order(nameOrder, weekOrder),]


idx <- order(nameOrder, weekOrder)
PCA_ind_reduced_ordered <- PCA_ind_reduced[idx, ]
nameOrder_ordered <- nameOrder[idx]

matList <- split(as.data.frame(PCA_ind_reduced_ordered), nameOrder_ordered)


# Create long format matrix of cluster assignments
LS <- split(cluster_ind$cluster, nameOrder)

LS_long <- imap_dfr(LS, function(v, ind) {
  v <- unlist(v, use.names = TRUE)
  nm <- names(v)
  week_tag <- sub(".*\\s+", "", nm)                       # last token
  week_num <- ifelse(week_tag == "F", 0L,
                     as.integer(sub("^D", "", week_tag)))
  tibble(
    Individual = ind,
    WeekTag    = week_tag,
    WeekNum    = week_num,
    Cluster    = as.integer(v)
  )
}) %>%
  arrange(.data$Individual, .data$WeekNum, .data$WeekTag)


#  Build wide matrix Individuals × Weeks (0:15) from LS_long
weeks_all <- 0:15

mat <- LS_long %>%
  group_by(Individual, WeekNum) %>%
  summarise(Cluster = first(Cluster), .groups = "drop") %>%  # Or any aggregation
  complete(Individual, WeekNum = weeks_all) %>%
  pivot_wider(names_from = WeekNum, values_from = Cluster) %>%
  arrange(Individual) %>%
  rename(Individuals = Individual)


# remove individual with NA for all the weeks of the dispersal phase
mat <- mat %>%
  filter(rowSums(!is.na(select(., `1`:`15`))) > 0)

# Merge with death metadata
death <- read.csv("C:/Users/lfaure7/OneDrive/MEMOIRE M2/Sauvegarde - My Passport/memoire_aigle/wetransfer_death_document-csv_2023-06-30_1528/death_document.csv",
                  stringsAsFactors = FALSE)
death$individual.local.identifier <- gsub("\\(.*\\)", "", death$individual.local.identifier)
death$individual.local.identifier <- trimws(death$individual.local.identifier)

merged_df <- mat %>%
  left_join(death, by = c("Individuals" = "individual.local.identifier"))

# Replace NA weeks by 5 (death) or 20 (logger lost/stop/no info)
lc_col <- intersect(c("loss_cause", "cause_death"), names(merged_df))
if (length(lc_col) == 0) {
  stop("No loss cause column found. Expected one of: loss_cause, cause_death.")
}
lc_col <- lc_col[1]

week_cols <- intersect(colnames(merged_df), as.character(weeks_all))

# ensure week cols are numeric
merged_df[week_cols] <- lapply(merged_df[week_cols], as.numeric)

lc_vec <- merged_df[[lc_col]]

for (i in seq_len(nrow(merged_df))) {
  x  <- as.numeric(merged_df[i, week_cols])
  lc <- lc_vec[i]
  
  # rle on reversed vector to detect trailing NAs
  r <- rle(rev(is.na(x)))
  na_trail_length <- if (length(r$values) > 0 && r$values[1]) r$lengths[1] else 0
  n <- length(x)
  
  for (k in seq_along(x)) {
    if (is.na(x[k])) {
      # If NA is in the trailing block *and* loss cause is "death"
      if (!is.na(lc) && lc == "death" && k > (n - na_trail_length)) {
        x[k] <- 5L
      } else {
        x[k] <- 20L
      }
    }
  }
  
  merged_df[i, week_cols] <- as.list(x)
}


# We'll keep a working copy with only week columns and identifiers
merged_df_rownames <- merged_df %>%
  select(Individuals, all_of(week_cols), any_of(lc_col)) %>%
  column_to_rownames("Individuals")

# Majority cluster + frequency computations

threshold <- 8  # strict majority > 8 weeks
min_each  <- 2  # min count for "no clear majority" group

merged_df_rownames <- merged_df_rownames %>%
  rownames_to_column("Individuals") %>%
  rowwise() %>%
  mutate(
    Majority_cluster = {
      weeks <- c_across(all_of(week_cols))
      if (any(weeks == 5,  na.rm = TRUE)) {
        5L
      } else if (any(weeks == 20, na.rm = TRUE)) {
        20L
      } else {
        tab  <- sort(table(weeks), decreasing = TRUE)
        cats <- as.integer(names(tab))
        cnts <- as.integer(tab)
        
        top_cat    <- cats[1]
        top_cnt    <- cnts[1]
        second_cnt <- if (length(cnts) >= 2) cnts[2] else 0L
        third_cnt  <- if (length(cnts) >= 3) cnts[3] else NA_integer_
        
        num_behav <- length(cats)
        
        has_strict_majority <-
          top_cat %in% 1:3 &&
          top_cnt  >  threshold &&
          top_cnt  >  second_cnt &&
          (is.na(third_cnt) || second_cnt != third_cnt)
        
        all_three_big_enough <- num_behav >= 3 &&
          all(cnts[1:3] > min_each, na.rm = TRUE)
        
        if (num_behav >= 3 && all_three_big_enough && !has_strict_majority) {
          10L
        } else {
          top_cat
        }
      }
    },
    Freq_majority = {
      weeks <- c_across(all_of(week_cols))
      if (Majority_cluster %in% 1:4) sum(weeks == Majority_cluster) else NA_integer_
    }
  ) %>%
  ungroup()

# Example tweak you had:
merged_df_rownames$Majority_cluster[merged_df_rownames$Majority_cluster == 2] <- 4L

# count of 20
merged_df_rownames <- merged_df_rownames %>%
  mutate(n20 = rowSums(across(all_of(week_cols), ~ . == 20L), na.rm = TRUE))

# Ordering and plotting

custom_order <- c(20, 5, 1, 3, 4, 10)

ordered_with_majority <- merged_df_rownames %>%
  arrange(
    factor(Majority_cluster, levels = custom_order),
    dplyr::case_when(
      Majority_cluster == 20 ~ n20,
      TRUE                   ~ -Freq_majority
    )
  )

order1 <- ordered_with_majority %>%
  select(-Majority_cluster, -Freq_majority, -n20, -loss_cause)

melted_data_ordered <- reshape2::melt(order1, id.vars = "Individuals")
melted_data_ordered$value <- as.factor(melted_data_ordered$value)
melted_data_ordered <- melted_data_ordered %>%
  mutate(Individuals = factor(Individuals, levels = order1$Individuals))

# ------------------------------------------------------------------------------ plot Figure 2
LUS_accross_week <- ggplot(melted_data_ordered,
                           aes(x = variable, y = Individuals, fill = value)) +
  geom_tile(color = "white", size = 0.1) +
  labs(x = "Weeks after emigration date", y = "",
       fill = "Landscape utilization\nstrategies") +
  scale_fill_manual(values = c(
    "1" = "darkcyan", "2" = "azure3", "3" = "salmon",
    "5" = "gray7", "10" = "gray7", "20" = "white"
  )) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 9),
    legend.title = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.key.height = unit(1, "pt"),
    legend.key.width = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 7)))
print(LUS_accross_week)
# ------------------------------------------------------------------------------ plot Figure 2


# STEP 6 : CALCULATION OF THE EUCLIDEAN STEP DISTANCES BETWEEN EACH WEEK --------------------------------------------------------------------------------------------------------------------------------------------
PCA_dim_reduced <- PCA$ind$coord[, 1:7] # keep first 7 dimensions (explains ~87%)
# Reorder by date 
nameOrder <- sapply(strsplit(rownames(PCA_dim_reduced), " "), function(str){paste(str[-length(str)], collapse=" ")}) 
weekOrder <- sapply(strsplit(rownames(PCA_dim_reduced), " "), function(str){str[length(str)]}) 
weekOrder <- as.numeric(sub("D","",sub("F","0",weekOrder))) # change the D into numerical values and F to 0

PCA_dim_reduced_ordered <- PCA_dim_reduced[order(nameOrder, weekOrder), ]

# suppress all the rows with individual name that end by "F" to suppress the fledging period (pre-dispersal)
PCA_dim_reduced_ordered_withoutF <- PCA_dim_reduced_ordered[!grepl(" F$", rownames(PCA_dim_reduced_ordered)), ]

# reorder by names
nameSplit <- sapply(strsplit(rownames(PCA_dim_reduced_ordered_withoutF), " "), function(str){paste(str[-length(str)], collapse=" ")})
matList <- split(as.data.frame(PCA_dim_reduced_ordered_withoutF), nameSplit)

# Calculate distances between each week of the dispersal period (from week 1 to week 2, etc.)
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
    week_pair  = names(dist_vec),
    distance   = as.numeric(dist_vec)
  )
}))

dist_nonnull <- distMatList[!sapply(distMatList, is.null)] # filter to remove individuals that have less than 2 locations

obs_distances <- map_dfr(
  dist_nonnull,           
  ~ as_tibble_row(.x),    
  .id = "individual"      
)

# Individual name, sum and mean of distances
obs_stat <- obs_distances %>%
  rowwise() %>%
  mutate(
    obs_mean = mean(c_across(starts_with("W")), na.rm = TRUE),
    obs_sum  = sum (c_across(starts_with("W")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(individual, obs_mean, obs_sum)



# STEP 7 : GENERATE THE NULL DISTRIBUTION PER INDIVIDUALS ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Randomise within-week only: permute the Dim.1–Dim.6 vectors among individuals present in the SAME week.
set.seed(123)  # ensure reproducibility 

for_rdn <- matList[sapply(matList, nrow) >= 2]               # keep individuals with >=2 weeks
names(for_rdn) <- names(matList)[sapply(matList, nrow) >= 2]

# prepare data
for_rdn <- lapply(names(for_rdn), function(id) {
  df <- for_rdn[[id]]
  df$week       <- seq_len(nrow(df))                         
  df$individual <- id
  rownames(df) <- NULL
  df
})

combined_for_rdn <- bind_rows(for_rdn)


dims <- paste0("Dim.", 1:6)
n_permutation <- 1000

random_distances <- lapply(seq_len(n_permutation), function(p){
  
  # Split rows by week, and inside each week permute only the Dim columns.
  wk_splits <- split(combined_for_rdn, combined_for_rdn$week)
  
  wk_shuffled <- lapply(wk_splits, function(dfwk) {
    if (nrow(dfwk) <= 1L) return(dfwk)                 
    perm_idx <- sample.int(nrow(dfwk))
    dfwk[dims] <- dfwk[perm_idx, dims, drop = FALSE]   
    dfwk
  })
  
  # Recombine all weeks: each individual keeps its own weeks, but coordinates were reassigned within-week.
  shuffled_all <- bind_rows(wk_shuffled)
  
  # split by individual to compute step distances along increasing week index
  shuffled_list <- split(shuffled_all, shuffled_all$individual)
  
  dist_list2 <- lapply(shuffled_list, function(ind) {
    ind <- ind[order(ind$week), , drop = FALSE]        
    coords <- ind[, dims, drop = FALSE]
    if (nrow(coords) < 2L) return(NULL)
    
    distances <- sapply(1:(nrow(coords) - 1L), function(i) {
      sqrt(sum((as.numeric(coords[i, ]) - as.numeric(coords[i + 1L, ]))^2))
    })
    
    names(distances) <- paste0("W", 1:(length(distances)), "_to_W", 2:(length(distances) + 1L))
    t(distances) %>% as.data.frame()
  }) %>%
    bind_rows(.id = "individual") %>%
    mutate(permutation = p) %>%
    as.data.frame()
  
  dist_list2
}) %>%
  bind_rows() %>%
  as.data.frame()

# compute mean and sum per row from random_distances
random_stat <- random_distances %>%
  mutate(
    mean_dist_rdn = rowMeans(dplyr::select(., -individual, -permutation), na.rm = TRUE),
    sum_dist_rdn  = rowSums (dplyr::select(., -individual, -permutation), na.rm = TRUE)
  )

# extract observed and random values for each individual & compute p-values
p_vals <- lapply(unique(obs_stat$individual), function(bird) {
  
  obs <- obs_stat    %>% filter(individual == bird)
  rnd <- random_stat %>% filter(individual == bird)
  
  # one-sided p-values for mean
  obs$p_less_mean <- sum(rnd$mean_dist_rdn <= obs$obs_mean, na.rm = TRUE) / n_permutation
  obs$p_more_mean <- sum(rnd$mean_dist_rdn >= obs$obs_mean, na.rm = TRUE) / n_permutation
  
  # one-sided p-values for sum
  obs$p_less_sum  <- sum(rnd$sum_dist_rdn  <= obs$obs_sum,  na.rm = TRUE) / n_permutation
  obs$p_more_sum  <- sum(rnd$sum_dist_rdn  >= obs$obs_sum,  na.rm = TRUE) / n_permutation
  
  obs
}) %>%
  bind_rows()


# ------------------------------------------------------------------------------ plot Figure S2.B
ord_vec <- order1$Individuals 

p_plot <- tibble(individual = ord_vec) %>%
  left_join(p_vals, by = "individual") %>%
  left_join(
    ordered_with_majority %>% select(Individuals, Majority_cluster),
    by = c("individual" = "Individuals")
  ) %>%
  mutate(
    individual_f = factor(individual, levels = ord_vec),
    clus = factor(Majority_cluster)
  ) %>%
  arrange(individual_f)

cols <- c(`1` = "darkcyan",
          `3` = "salmon",
          `4` = "azure3",
          `5` = "gray6",
          `10` = "black",
          `20` = "black")  


p_plot <- p_plot %>%
  mutate(
    shape = case_when(
      clus %in% c(10, 20) ~ 21,  # circle with contour line
      TRUE                ~ 16   # plain circle
    ),
    fill = case_when(
      clus == 10 ~ p_less_mean,  
      clus == 20 ~ NA_real_,     
      TRUE       ~ NA_real_))


rdn <- ggplot(
  p_plot,
  aes(y = individual_f, x = p_less_mean,
      colour = clus, fill = fill, shape = shape)
) +
  geom_segment(aes(x = 0, xend = p_less_mean, yend = individual_f),
               linewidth = 2, colour = "grey80") +
  geom_point(size = 4, stroke = 0.7) +
  scale_colour_manual(values = cols, name = "Cluster") +
  scale_shape_identity() +
  scale_fill_gradientn(
    colours = c("salmon", "darkcyan"),
    na.value = "white", guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, .02))) +
  labs(x = "p-value (p_less_mean)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor     = element_blank())
print(rdn)

# ------------------------------------------------------------------------------ Supplementary material : Figure 4

rs <- random_stat %>%
  left_join(obs_stat, by = "individual")

rs_plot <- rs %>%
  semi_join(tibble(individual = ord_vec), by = "individual") %>%  # keep only ordered birds
  mutate(individual = factor(individual, levels = ord_vec)) %>%   # set facet order
  arrange(individual)

RandomVsObservedMean <- ggplot(rs_plot, aes(x = mean_dist_rdn)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "grey90", colour = "black") +
  geom_density(colour = "darkred") +
  geom_vline(aes(xintercept = obs_mean), color = "black", size = 0.9) +
  facet_wrap(~ individual, scales = "free_y", ncol = 5) +
  labs(
    title = "Null distribution of the mean Euclidean step distances (within-week shuffle)",
    x = "Mean distance",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(size = 8),
    axis.text.y      = element_text(size = 6),   
    axis.title.y     = element_text(size = 6),
    panel.spacing    = unit(0.5, "lines"),
    panel.grid.minor = element_blank()
  ) 

RandomVsObservedMean
# ------------------------------------------------------------------------------ Supplementary material : Figure 4


# STEP 8 : distance between pre-dispersal and 7-day sequence (observed and randomly generated) ----------------------------------------------------------------------------------------------------------------------
# reorder by names
name <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "), function(str){paste(str[-length(str)], collapse=" ")})
List <- split(as.data.frame(PCA_dim_reduced_ordered), name)

# Calculate distances between natal territory and each week of the dispersal period 
MahalanobisList <- lapply(List, function(ind){
  MahalanobisDist <- dist(ind, method = "euclidean", diag = TRUE, upper = TRUE)
  return(MahalanobisDist)
})

# Convert the titles of the columns to 0, 1, 2, 3 ...
MahalanobisList <- lapply(MahalanobisList, function(distance) {
  distance <- as.matrix(distance)  # Convert distance object to matrix
  colnames(distance) <- seq_len(ncol(distance)) - 1
  distance
})

# Keep the first row of each matrix (distance from natal territory, week 0)
MahalanobisList <- lapply(MahalanobisList, function(distance) {
  distance[1, , drop = FALSE]
})

# Find the maximum number of columns
max_col <- max(sapply(MahalanobisList, ncol))

# Add missing columns with NA values to align the matrices
alignedMahalanobisList <- lapply(MahalanobisList, function(distance) {
  Num_missing_cols <- max_col - ncol(distance)
  if (Num_missing_cols > 0) {
    distance <- cbind(distance, matrix(NA, nrow = nrow(distance), ncol = Num_missing_cols))
  }
  distance
})

# Bind each row and reorganise dataset
combinedDistances <- do.call(rbind, alignedMahalanobisList)

combinedDistances2 <- combinedDistances %>%
  as.data.frame() %>%                        
  rownames_to_column("individual") %>%      
  pivot_longer(
    cols      = -individual,
    names_to  = "week",
    values_to = "distance"
  ) 


# Compare distances to week 0: observed vs week-stratified randomisation
dims6 <- paste0("Dim.", 1:6)

# Build long table from PCA_dim_reduced_ordered (contains week 0 = "F")
name_vec <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "),
                   function(str){paste(str[-length(str)], collapse=" ")})
week_chr <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "),
                   function(str){str[length(str)]})
week_num <- as.numeric(sub("D","", sub("F","0", week_chr)))  # F->0, Dk->k

df_long6 <- as.data.frame(PCA_dim_reduced_ordered)[, dims6]
df_long6$individual <- name_vec
df_long6$week <- week_num
df_long6 <- df_long6 %>% arrange(individual, week)

# Keep only individuals that have week 0 (anchor)
inds_with_w0 <- df_long6 %>%
  filter(week == 0) %>%
  pull(individual) %>%
  unique()

df_long6 <- df_long6 %>% filter(individual %in% inds_with_w0)

# Map week 0 vectors
w0_by_ind <- df_long6 %>%
  filter(week == 0) %>%
  group_by(individual) %>%
  slice(1) %>%
  ungroup() %>%
  { setNames(split(.[, dims6], .$individual), .$individual) }

# Build a table of week-0 coordinates per individual
w0_df <- df_long6 %>%
  filter(week == 0) %>%
  select(individual, all_of(dims6)) %>%
  rename_with(~ paste0(., "_w0"), all_of(dims6))

# Observed distances to week 0 for each (individual, week > 0)
obs_dw0 <- df_long6 %>%
  filter(week > 0) %>%
  left_join(w0_df, by = "individual") %>%
  mutate(
    dist_obs = {
      cur  <- as.matrix(select(., all_of(dims6)))
      base <- as.matrix(select(., paste0(dims6, "_w0")))
      sqrt(rowSums((cur - base)^2))
    }
  ) %>%
  select(individual, week, dist_obs) %>%
  arrange(individual, week)

# Per-week observed summary
obs_week_stats <- obs_dw0 %>%
  group_by(week) %>%
  summarise(
    n_obs = n(),
    mean_obs = mean(dist_obs, na.rm = TRUE),
    median_obs = median(dist_obs, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare week-wise pools for randomisation
by_week <- split(df_long6, df_long6$week)

# Safety check: each week >0 should have at least 2 distinct individuals
check_weeks <- sapply(names(by_week), function(w){
  length(unique(by_week[[w]]$individual))
})


n_permutation <- 1000
set.seed(123)

# Randomisation: for each week w>0 and individual i present at w : take one other individual j in i present at week w; compute distance between w0(i) and vec(j, w).
random_dw0 <- lapply(seq_len(n_permutation), function(p){
  
  per_week <- lapply(names(by_week), function(wc){
    w <- as.integer(wc)
    if (w == 0) return(NULL)
    
    dfw <- by_week[[wc]]
    # Filter to individuals that have week0
    dfw <- dfw[dfw$individual %in% names(w0_by_ind), , drop = FALSE]
    if (nrow(dfw) == 0) return(NULL)
    
    # For each focal individual at this week
    rows <- lapply(seq_len(nrow(dfw)), function(r){
      id_i <- dfw$individual[r]
      # pool of others at same week
      pool <- dfw[dfw$individual != id_i, , drop = FALSE]
      if (nrow(pool) == 0) return(NULL)  
      
      pick <- pool[sample.int(nrow(pool), 1L), dims6, drop = FALSE]
      vj <- as.numeric(pick[1, ])
      
      w0_row <- w0_by_ind[[id_i]]
      if (is.null(w0_row)) return(NULL)
      w0 <- as.numeric(w0_row[1, dims6, drop = TRUE])
      
      dist_rdn <- sqrt(sum((vj - w0)^2))
      
      tibble(individual = id_i, week = w, permutation = p, dist_rdn = dist_rdn)
    })
    
    bind_rows(rows)
  })
  
  bind_rows(per_week)
}) %>%
  bind_rows()


# ------------------------------------------------------------------------------ Supplementary material : Per-(individual, week) null distributions
# Merge random and observed at (i, w)
iw_merged <- random_dw0 %>%
  left_join(obs_dw0, by = c("individual", "week"))

# Safety: remove pairs that lack observed distance
iw_merged <- iw_merged %>% filter(!is.na(dist_obs))


# 1) Density per (individual, week) with observed line

# Choose which individuals and weeks to show
if (exists("order1") && "Individuals" %in% names(order1)) {
  full_ord <- order1$Individuals
} else {
  full_ord <- sort(unique(iw_merged$individual))
}

inds_to_plot  <- full_ord                 # or e.g. head(full_ord, 12)
weeks_to_plot <- sort(unique(iw_merged$week))  # or, e.g., 1:10

dens_data <- iw_merged %>%
  filter(individual %in% inds_to_plot,
         week %in% weeks_to_plot) %>%
  mutate(individual = factor(individual, levels = full_ord),
         week = factor(week, levels = sort(unique(week))))

obs_lines <- obs_dw0 %>%
  filter(individual %in% inds_to_plot,
         week %in% weeks_to_plot) %>%
  mutate(individual = factor(individual, levels = full_ord),
         week = factor(week, levels = sort(unique(week))))

Plot_IW_Density <- ggplot(dens_data, aes(x = dist_rdn)) +
  geom_density(fill = "grey80", colour = "black", adjust = 1.1) +
  geom_vline(data = obs_lines, aes(xintercept = dist_obs),
             color = "black", linewidth = 0.7) +
  facet_grid(individual ~ week, scales = "free_y") +
  labs(
    title = "Per-(individual, week) null distributions of distance to week 0",
    x = "Randomised distance (other individual at same week)",
    y = "Density"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text.x    = element_text(size = 8),
    strip.text.y    = element_text(size = 8),
    axis.text.y     = element_text(size = 6),
    axis.title.y    = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

Plot_IW_Density


#WILCOXON TEST
# Compute null mean per (i,w):
null_mean_iw <- random_dw0 %>%
  group_by(individual, week) %>%
  summarise(dist_null_mean = mean(dist_rdn, na.rm = TRUE),
            .groups = "drop")

# Assemble paired differences Δ = obs – null_mean:
iw_diff <- obs_dw0 %>%
  inner_join(null_mean_iw, by = c("individual", "week")) %>%
  mutate(delta = dist_obs - dist_null_mean)

# Wilcoxon signed‑rank test (one‑sided: delta < 0)
wilcox.test(iw_diff$delta,
            alternative = "less",
            paired      = FALSE,   
            conf.level  = 0.95)

# Paired t‑test
t.test(iw_diff$delta,
       alternative = "less",
       paired      = FALSE)
