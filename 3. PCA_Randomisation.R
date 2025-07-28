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
