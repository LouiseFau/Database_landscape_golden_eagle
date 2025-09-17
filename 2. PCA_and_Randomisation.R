# R. Script to identify variations between and within individuals
# This script is following the second loop of the R.Script called "Database_and_variables_selection". This script groups individuals based on their landscape utilisation strategies defines during the Steps 1 to 2, 
# using a clustering method (Step 3 to 4). Then, it identifies the consistency in landscape features used by measuring the Mahalanobis distance between (I) consecutive individual-week points of the dispersal phase 
# compared to a random distribution (Step 4 to 5) and (II) pre-dispersal and each of the individual-week points of the dispersal phase (Step 6). 
# L. Faure, MPI, 17.09.2025

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
library(gt)

eagle_PCA <- readRDS("your previous repository path.rds")

# STEP 1 : SELECTION OF PCA DIMENSION ------------------------------------------
# Copy of the first dataset without the weight column
eagle_without_weight <- eagle_PCA[, !(colnames(eagle_PCA) %in% c("weight"))]

# Add row names to display the names of the birds in the PCA
PCA_data <- eagle_without_weight %>%
  column_to_rownames(var = "individual.local.identifier.per.week")

# PCA loadings
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
PCA_ind_reduced <- PCA$ind$coord[, 1:4] 
PCA_var_reduced <- PCA$var$coord[, 1:4] 



# STEP 2 : PCA & LANDSCAPE UTILIZATION STRATEGIES-------------------------------
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



# STEP 4 : CLUSTER STRENGHT INDICATORS -----------------------------------------
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

# Standardize
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
# ------------------------------------------------------------------------------(end) Supplementary material : Figure 3

# 30 indices (Charrad et al., 2014)
df_scaled <- scale(PCA_ind_reduced)
NbClust(PCA_ind_reduced, diss = NULL, distance = "manhattan", min.nc = 2, max.nc = 6, method = "kmeans") # optimal number of cluster : 2, 6 or 3. 


# Bayesian Information Criterion
d_clust <- Mclust(PCA_ind_reduced, G = 1:6, 
                  modelNames = mclust.options("emModelNames"))
d_clust$BIC # optimal number of clusters : 5, 3 or 6. 





# STEP 5 : INDIVIDUAL LANDSCAPE UTILIZATION STRATEGY CONSISTENCY ---------------
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
  week_tag <- sub(".*\\s+", "", nm)                       
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


#  Build wide matrix individuals/weeks 
weeks_all <- 0:15

mat <- LS_long %>%
  group_by(Individual, WeekNum) %>%
  summarise(Cluster = first(Cluster), .groups = "drop") %>%  # Or any aggregation
  complete(Individual, WeekNum = weeks_all) %>%
  pivot_wider(names_from = WeekNum, values_from = Cluster) %>%
  arrange(Individual) %>%
  rename(Individuals = Individual)


# Remove individual with NA for all the weeks of the dispersal phase
mat <- mat %>%
  filter(rowSums(!is.na(select(., `1`:`15`))) > 0)

# Merge with death metadata
death <- read.csv("death_document.csv",
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

# Ensure week cols are numeric
merged_df[week_cols] <- lapply(merged_df[week_cols], as.numeric)

lc_vec <- merged_df[[lc_col]]

for (i in seq_len(nrow(merged_df))) {
  x  <- as.numeric(merged_df[i, week_cols])
  lc <- lc_vec[i]
  
  r <- rle(rev(is.na(x)))
  na_trail_length <- if (length(r$values) > 0 && r$values[1]) r$lengths[1] else 0
  n <- length(x)
  
  for (k in seq_along(x)) {
    if (is.na(x[k])) {
      if (!is.na(lc) && lc == "death" && k > (n - na_trail_length)) {
        x[k] <- 5L
      } else {
        x[k] <- 20L
      }
    }
  }
  
  merged_df[i, week_cols] <- as.list(x)
}


# working copy with only week columns and identifiers
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

merged_df_rownames$Majority_cluster[merged_df_rownames$Majority_cluster == 2] <- 4L

# Count of 20
merged_df_rownames <- merged_df_rownames %>%
  mutate(n20 = rowSums(across(all_of(week_cols), ~ . == 20L), na.rm = TRUE))




# STEP 5 : CONSISTENCY IN LANDSCAPE USE ACROSS WEEKS ---------------------------
PCA_dim_reduced_rdn <- PCA$ind$coord[, 1:7] # keep first 7 dimensions (explains ~87%)

# Reorder by date 
nameOrder <- sapply(strsplit(rownames(PCA_dim_reduced_rdn), " "), function(str){paste(str[-length(str)], collapse=" ")}) 
weekOrder <- sapply(strsplit(rownames(PCA_dim_reduced_rdn), " "), function(str){str[length(str)]}) 
weekOrder <- as.numeric(sub("D","",sub("F","0",weekOrder))) # change the D into numerical values and F to 0

PCA_dim_reduced_ordered <- PCA_dim_reduced_rdn[order(nameOrder, weekOrder), ]

# Suppress all the rows with individual name that end by "F" to suppress the fledging period (pre-dispersal)
PCA_dim_reduced_ordered_withoutF <- PCA_dim_reduced_ordered[!grepl(" F$", rownames(PCA_dim_reduced_ordered)), ]

# Reorder by names
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


# Generate the null distribution per individuals
set.seed(123) 
dims <- paste0("Dim.", 1:6) # permutation on the first 6 PCs
n_permutation <- 1000 # number of permutation

for_rdn <- matList[sapply(matList, nrow) >= 2]               # keep individuals with >=2 weeks
names(for_rdn) <- names(matList)[sapply(matList, nrow) >= 2]

# Data preparation
for_rdn <- lapply(names(for_rdn), function(id) {
  df <- for_rdn[[id]]
  df$week       <- seq_len(nrow(df))                         # within-individual sequential index 1..15 (max)
  df$individual <- id
  rownames(df) <- NULL
  df
})

combined_for_rdn <- bind_rows(for_rdn)

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
  
# Split by individual to compute step distances along increasing week index
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

# Compute mean per row 
random_stat <- random_distances %>%
  mutate(mean_dist_rdn = rowMeans(dplyr::select(., -individual, -permutation), na.rm = TRUE))

# Extract observed and random values for each individual & compute p-values
p_vals <- lapply(unique(obs_stat$individual), function(bird) {
  
  obs <- obs_stat    %>% filter(individual == bird)
  rnd <- random_stat %>% filter(individual == bird)
  
  obs$p_less_mean <- sum(rnd$mean_dist_rdn <= obs$obs_mean, na.rm = TRUE) / n_permutation
  obs$p_more_mean <- sum(rnd$mean_dist_rdn >= obs$obs_mean, na.rm = TRUE) / n_permutation
  
  obs
}) %>%
  bind_rows()



#  STEP 6 : CONSISTENCY IN LANDSCAPE USE BETWEEN PRE-DISPERSAL AND 7-DAYS SEQUENCE
name_vec <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "),
                   function(str){paste(str[-length(str)], collapse=" ")})
week_chr <- sapply(strsplit(rownames(PCA_dim_reduced_ordered), " "),
                   function(str){str[length(str)]})
week_num <- as.numeric(sub("D","", sub("F","0", week_chr)))  # F->0, Dk->k

df_long6 <- as.data.frame(PCA_dim_reduced_ordered)[, dims]
df_long6$individual <- name_vec
df_long6$week <- week_num
df_long6 <- df_long6 %>% arrange(individual, week)

# Keep only individuals with pre-dispersal period
inds_with_w0 <- df_long6 %>%
  filter(week == 0) %>%
  pull(individual) %>%
  unique()

df_long6 <- df_long6 %>% filter(individual %in% inds_with_w0)

w0_by_ind <- df_long6 %>%
  filter(week == 0) %>%
  group_by(individual) %>%
  slice(1) %>%
  ungroup() %>%
  { setNames(split(.[, dims], .$individual), .$individual) }

# Table of coordinates for pre-dispersal
w0_df <- df_long6 %>%
  filter(week == 0) %>%
  select(individual, all_of(dims)) %>%
  rename_with(~ paste0(., "_w0"), all_of(dims))

# Compute distances from pre-dispersal to each week of the dispersal phase 
obs_dw0 <- df_long6 %>%
  filter(week > 0) %>%
  left_join(w0_df, by = "individual") %>%
  mutate(
    dist_obs = {
      cur  <- as.matrix(select(., all_of(dims)))
      base <- as.matrix(select(., paste0(dims, "_w0")))
      sqrt(rowSums((cur - base)^2))
    }
  ) %>%
  select(individual, week, dist_obs) %>%
  arrange(individual, week)


by_week <- split(df_long6, df_long6$week)

# Generate the null distribution
random_dw0 <- lapply(seq_len(n_permutation), function(p){
  per_week <- lapply(names(by_week), function(wc){
    w <- as.integer(wc)
    if (w == 0) return(NULL)
    dfw <- by_week[[wc]]
    dfw <- dfw[dfw$individual %in% names(w0_by_ind), , drop = FALSE]
    if (nrow(dfw) == 0) return(NULL)
    rows <- lapply(seq_len(nrow(dfw)), function(r){
      id_i <- dfw$individual[r]
      pool <- dfw[dfw$individual != id_i, , drop = FALSE]
      if (nrow(pool) == 0) return(NULL)
      pick <- pool[sample.int(nrow(pool), 1L), dims, drop = FALSE]
      vj <- as.numeric(pick[1, ])
      w0_row <- w0_by_ind[[id_i]]
      if (is.null(w0_row)) return(NULL)
      w0 <- as.numeric(w0_row[1, dims, drop = TRUE])
      dist_rdn <- sqrt(sum((vj - w0)^2))
      tibble(individual = id_i, week = w, permutation = p, dist_rdn = dist_rdn)
    })
    bind_rows(rows)
  })
  bind_rows(per_week)
}) %>% bind_rows()

# Compute differences RANDOM - OBSERVED for all randomizations to identify individual that stick to pre-dispersal behavior
delta_iwp <- random_dw0 %>%
  left_join(obs_dw0, by = c("individual", "week")) %>%
  filter(!is.na(dist_obs)) %>%
  mutate(delta = dist_rdn - dist_obs)

# One mean difference per (individual, permutation)
perm_means <- delta_iwp %>% 
  group_by(individual, permutation) %>%
  summarise(bar_delta = mean(delta, na.rm = TRUE), .groups = "drop")

ind_pvals <- perm_means %>% 
  group_by(individual) %>%
  summarise(
    p_val = (sum(bar_delta <= 0) + 1) / (n() + 1),   # +1 correction
    .groups = "drop"
  )

# ------------------------------------------------------------------------------ plot Figure 2
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

ord_vec <- levels(melted_data_ordered$Individuals)

pvals_joined <- order1 %>%
  dplyr::select(Individuals) %>%
  # p-values randomisation 1 
  left_join(
    p_vals %>% dplyr::select(Individuals = individual, p_less_mean),
    by = "Individuals"
  ) %>%
  # p-values randomisation 2 
  left_join(
    ind_pvals %>% dplyr::select(Individuals = individual, p_val),
    by = "Individuals"
  )

order1_ext <- order1 %>%
  mutate(
    `Consistency_7d` = 20L,   # 20 = white 
    `Consistency_w0` = 20L
  )


melted_data_ordered <- reshape2::melt(order1_ext, id.vars = "Individuals")
melted_data_ordered$value <- as.factor(melted_data_ordered$value)


melted_data_ordered <- melted_data_ordered %>%
  mutate(Individuals = factor(Individuals, levels = order1$Individuals))


stars_long <- pvals_joined %>%
  transmute(
    Individuals,
    Consistency_7d = ifelse(!is.na(p_less_mean) & p_less_mean < 0.05, "*", ""),
    Consistency_w0 = ifelse(!is.na(p_val)       & p_val       < 0.05, "*", "")
  ) %>%
  pivot_longer(-Individuals, names_to = "variable", values_to = "star") %>%
  filter(star == "*")  


x_levels <- c(as.character(0:15), "Consistency_7d", "Consistency_w0")
x_labels <- c(as.character(0:15), "", "")


melted_data_ordered$variable <- factor(melted_data_ordered$variable, levels = x_levels)
stars_long$variable          <- factor(stars_long$variable,          levels = x_levels)

LUS_accross_week <- ggplot(melted_data_ordered,
                           aes(x = variable, y = Individuals, fill = value)) +
  geom_tile(color = "white", size = 0.1) +
  labs(x = "Weeks after emigration date",
       y = "",
       fill = "Landscape utilization\nstrategies") +
  scale_x_discrete(limits = x_levels, labels = x_labels) +  
  scale_fill_manual(values = c(
    "1"  = "darkcyan",
    "2"  = "azure3",
    "3"  = "salmon",
    "5"  = "gray7",
    "10" = "gray7",
    "20" = "white"     
  )) +
  geom_text(
    data = stars_long,
    aes(x = variable, y = Individuals, label = star, colour = variable),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold"
  ) +
  scale_colour_manual(
    values = c("Consistency_7d" = "black", "Consistency_w0" = "darkorange"),
    guide = "none"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 9),
    legend.title = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.key.height = unit(1, "pt"),
    legend.key.width  = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 7)))

print(LUS_accross_week)

# ------------------------------------------------------------------------------ Supplementary material : Figure 4 and Figure 5
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
    title = "Null distribution of the mean Mahalanobis step distances",
    x = "Mean Mahalanobis distance",
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

tbl_consistency <- pvals_joined %>%
  dplyr::mutate(
    `Consistency across 7-days sequences` = p_less_mean,
    `Consistency between pre-dispersal and dispersal periods` = p_val
  ) %>%
  dplyr::select(
    Individuals,
    `Consistency across 7-days sequences`,
    `Consistency between pre-dispersal and dispersal periods`
  ) %>%
  dplyr::semi_join(tibble::tibble(Individuals = ord_vec), by = "Individuals") %>%
  dplyr::mutate(Individuals = factor(Individuals, levels = ord_vec)) %>%
  dplyr::arrange(Individuals) %>%
  dplyr::mutate(Individuals = as.character(Individuals))

headers <- c(
  "Consistency across 7-days sequences",
  "Consistency between pre-dispersal and dispersal periods"
)

col_px <- 140

gt_consistency <- tbl_consistency %>%
  gt::gt(rowname_col = "Individuals") %>%
  gt::cols_label(
    Individuals = gt::md("Individuals"),
    `Consistency across 7-days sequences` =
      gt::html("Consistency<br>across 7-days sequences"),
    `Consistency between pre-dispersal and dispersal periods` =
      gt::html("Consistency between<br>pre-dispersal and dispersal periods")
  ) %>%
  gt::cols_align(align = "center", columns = dplyr::everything()) %>%
  gt::cols_align(align = "left", columns = dplyr::all_of("Individuals")) %>%
  gt::fmt_number(columns = -Individuals, decimals = 3) %>%
  gt::fmt_missing(columns = dplyr::everything(), missing_text = "") %>%
  gt::cols_width(
    dplyr::everything() ~ gt::px(col_px)
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(12),
    data_row.padding = gt::px(4),
    column_labels.padding = gt::px(6)
  )

gt_consistency
