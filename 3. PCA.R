# R. Script to identify the different landscape utilization strategies through time comparing different individual
# This script is following the second loop of the R.Script called "Variables selection". The script include a PCA (Step 1 to 2) and a clustering method (Step 3 to 4).
# L. Faure, MPI, 23/06/2023

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

eagle_PCA <- readRDS("/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ANALYSIS/database_for_PCA/results_F_d15.rds")

# STEP 1 : SELECTION OF PCA DIMENSION ------------------------------------------
# Copy of the first dataset without the weight column
eagle_without_weight <- eagle_PCA[, !(colnames(eagle_PCA) %in% c("weight"))]

# Add row names to display the names of the birds in the PCA
PCA_data <- eagle_without_weight %>%
  column_to_rownames(var = "individual.local.identifier.per.week")

# PCA loadings (nb: the number of dimensions is equal to the number of variables)
PCA <- FactoMineR::PCA(PCA_data[,1:13], row.w = eagle_PCA$weight, ncp = 13)
summary(PCA)

# Scree plot of eigenvalues (scree plot of the eignevalues)
# Extract and prepare eigenvalue data
eig_df <- as.data.frame(get_eigenvalue(PCA))
eig_df$dim_num <- 1:nrow(eig_df)
eig_df$percentage <- round(eig_df$variance.percent, 1)

# Scree plot
custom_blues <- rev(colorRampPalette(brewer.pal(9, "Blues"))(nrow(eig_df)))
ggplot(eig_df, aes(x = dim_num, y = variance.percent)) +
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
      ) # 5 first dimensions selected

# Selection of the variable and individuals points based on the size (Jollife, 2002, Peres-Neto, Jackson ad Somers, 2005)
PCA_ind_reduced <- PCA$ind$coord[, 1:5] # variables
PCA_var_reduced <- PCA$var$coord[, 1:5] # individual 7-day sequence point



# STEP 2 : PCA & LANDSCAPE UTILIZATION STRATEGIES-------------------------------
# Cluster identification
set.seed(123)
cluster_ind <- kmeans(PCA_ind_reduced, centers = 3, nstart = 25)
grp_ind    <- as.factor(cluster_ind$cluster) 

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
cluster_palette <- c("1" = "salmon", "2" = "azure3", "3" = "darkcyan")
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

# standardize
cluster_results$WCSS <- scale(cluster_results$WCSS)[,1]
cluster_results$BCSS <- scale(cluster_results$BCSS)[,1]

cluster_results$diff <- abs(cluster_results$BCSS - cluster_results$WCSS)
cluster_results$alpha <- scales::rescale(max(cluster_results$diff) - cluster_results$diff, to = c(0.1, 0.8))


ribbon_df <- bind_rows(
      cluster_results %>% select(x = ClusterNumber, y = WCSS, alpha),
      cluster_results %>% arrange(desc(ClusterNumber)) %>%
            select(x = ClusterNumber, y = BCSS, alpha))

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

# 30 indices (Charrad et al., 2014)
df_scaled <- scale(PCA_ind_reduced)
NbClust(PCA_ind_reduced, diss = NULL, distance = "manhattan", min.nc = 2, max.nc = 6, method = "kmeans") # optimal number of cluster : 4, 2 or 3/6. 

# Bayesian Information Criterion
d_clust <- Mclust(PCA_ind_reduced, G = 1:6, 
                  modelNames = mclust.options("emModelNames"))
d_clust$BIC # optimal number of clusters : 5, 6 or 3.




# STEP 5 : INDIVIDUAL LANDSCAPE UTILIZATION STRATEGY CONSISTENCY ---------------
# Re order by date 
nameOrder <- sapply(strsplit(rownames(PCA_ind_reduced), " "), function(str){paste(str[-length(str)], collapse=" ")}) # first part of the code divide the matrix per space and get rid of the last part with space
weekOrder <- sapply(strsplit(rownames(PCA_ind_reduced), " "), function(str){str[length(str)]}) # merge the names that have a space
weekOrder <- as.numeric(sub("D","",sub("F","0",weekOrder))) # change the D into numerical values and F to 0

PCA_ind_reduced_ordered <- PCA_ind_reduced[order(nameOrder, weekOrder),]

# reorder by names
matList <- split(as.data.frame(PCA_ind_reduced_ordered), nameOrder)
combined_df <- do.call(rbind, lapply(seq_along(matList), function(i) {
      df <- matList[[i]]
      rownames(df) <- paste0(rownames(df), ".", i)  # Add a unique identifier
      df}))

# Reset the row names without the unique identifier
rownames(combined_df) <- gsub("\\.\\d+$", "", rownames(combined_df))

# Create long format matrix of cluster assignments
ls <- split(cluster_ind$cluster, nameOrder)

ls_T <- lapply(1:length(ls), function(i) {
  dfT <- t(as.matrix(ls[[i]]))
  rownames(dfT) <- names(ls)[i]
  colnames(dfT) <- 0:(ncol(dfT) - 1)  
  dfT})


df_list <- lapply(seq_along(ls_T), function(i) {
  df <- as.data.frame(ls_T[[i]])
  df})

mat <- bind_rows(df_list, .id = "Cluster")
mat <- subset(mat, select = -Cluster)

#Merge with death metadata and remove rows with missing values
mat <- mat[!is.na(mat$`1`), ]
mat <- mat %>%
  rownames_to_column(var = "Individuals")
mat$Individuals <- trimws(mat$Individuals)

death <- read.csv("/Volumes/My Passport/Louise/memoire_aigle/wetransfer_death_document-csv_2023-06-30_1528/death_document.csv")
death$individual.local.identifier <- gsub("\\(.*\\)", "", death$individual.local.identifier)
death$individual.local.identifier <- trimws(death$individual.local.identifier)
death$individual.local.identifier <- trimws(death$individual.local.identifier)

death_filtered <- death[death$individual.local.identifier %in% mat$Individuals, ]

# Controle the matching names
#missing_names <- setdiff(mat$Individuals, death_filtered$individual.local.identifier)
#print(missing_names)

# Merge the dataframe based on matching names
merged_df <- merge(mat, death, by.x = "Individuals", by.y = "individual.local.identifier", all.x = TRUE)

# Create rownames
merged_df_rownames <- merged_df %>%
  column_to_rownames(var = "Individuals")
merged_df_rownames <- merged_df_rownames[, !(names(merged_df_rownames) %in% c("cause_death", "loss_location"))]

# For the rows that have "death", change NA to number 4
merged_df_rownames <- merged_df_rownames %>%
      mutate(across(
            where(is.numeric),
            ~ if_else(loss_cause == "death" & is.na(.), as.numeric(4), as.numeric(.))
      ))
merged_df_rownames <- merged_df_rownames[complete.cases(merged_df_rownames), ]
  
# Rownames to columns
merged_df_rownames <- merged_df_rownames %>%
  rownames_to_column(var = "Individuals")
merged_df_rownames <- merged_df_rownames[, !(names(merged_df_rownames) %in% c("loss_cause"))]

# calculate majority and frequency of each strategy
merged_df_rownames <- merged_df_rownames %>%
      rowwise() %>%
      mutate(
            Majority_cluster = {
                  vals <- c_across(`0`:`15`)
                  tab <- sort(table(vals), decreasing = TRUE)
                  top_val <- as.integer(names(tab)[1])
                  if ((length(tab) > 1 && tab[1] == tab[2]) || !(top_val %in% c(1, 2, 3))) {
                        4  # no majority or majority is not 1/2/3
                  } else {
                        top_val
                  }
            },
            Freq_majority = {
                  vals <- c_across(`0`:`15`)
                  tab <- table(vals)
                  top_val <- as.integer(names(tab)[which.max(tab)])
                  if ((length(tab) > 1 && sort(tab, decreasing = TRUE)[1] == sort(tab, decreasing = TRUE)[2]) ||
                      !(top_val %in% c(1, 2, 3))) {
                        NA  # no clear majority or invalid category
                  } else {
                        max(tab)
                  }
            }
      ) %>%
      ungroup()

# make sure the category 4 appear at the bottom of the table and suppress the individual Almen18
merged_df_rownames <- merged_df_rownames %>% filter(Individuals != "Almen18") # individual deseased after 2 weeks
custom_order <- c("Flüela20", "Almen18", "Sils20", "Flüela1 21", "ValGrande19")

# Separate and sort individuals in two groups (with and without majority)
ordered_with_majority <- merged_df_rownames %>%
      filter(Majority_cluster != 4) %>%
      arrange(Majority_cluster, desc(Freq_majority))

ordered_no_majority <- merged_df_rownames %>%
      filter(Majority_cluster == 4) %>%
      filter(Individuals %in% custom_order) %>%
      mutate(Individuals = factor(Individuals, levels = custom_order)) %>%
      arrange(Individuals)

# Bind both datasets together
merged_df_rownames <- bind_rows(ordered_with_majority, ordered_no_majority)
merged_df_rownames$Individuals <- factor(merged_df_rownames$Individuals,
                                         levels = merged_df_rownames$Individuals)

# suppress temporary columns
merged_df_rownames <- merged_df_rownames %>%
      select(-Majority_cluster, -Freq_majority)

# change format for ggplot2
melted_data_ordered <- reshape2::melt(merged_df_rownames, id.vars = "Individuals")
melted_data_ordered$value <- factor(melted_data_ordered$value)

# Raster plot : change colours code based on the last colour code for the PCA
ggplot(melted_data_ordered, aes(x = variable, y = Individuals, fill = value)) +
      geom_tile(color = "white", size = 0.1) +
      labs(x = "Weeks after emigration date", y = "", fill = "Landscape utilization\nstrategies") +
      scale_fill_manual(values = c(
            "1" = "salmon", "2" = "azure3", "3" = "darkcyan", "4" = "gray7", "10" = "#999999"    # check colours based on the last PCA biplot colour code
      )) +
      scale_y_discrete(limits = rev(levels(melted_data_ordered$Individuals))) +  # <--- fix here
      theme_minimal() +
      theme(
            text = element_text(family = "Avenir"),
            axis.title = element_text(face = "bold", size = 9),
            legend.title = element_text(size = 10),
            axis.title.x = element_text(size = 12),
            legend.text = element_text(size = 8),
            legend.key.height = unit(1, "pt"),
            legend.key.width = unit(0.5, "cm")
      ) +
      guides(fill = guide_legend(override.aes = list(size = 7)))

# fledging landscape utilization strategy and export 
lus_fledging <- merged_df_rownames %>%
      select(Individuals, `0`) %>%
      mutate(`0` = recode(`0`,
                          `1` = "rangetop",
                          `2` = "valley",
                          `3` = "forest"))

#saveRDS(lus_fledging, file = "/Volumes/My Passport/Louise/memoire_aigle/Kami_ordi/wetransfer_basic-analysis_2023-07-25_1534/Basic-analysis/ARTICLE/supplementary_data/LUS_fledging.rds")
