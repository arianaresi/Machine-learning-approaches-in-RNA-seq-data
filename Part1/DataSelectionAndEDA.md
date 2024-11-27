# Part 0: Load libraries 

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#----01 package loading-----
library("DESeq2")
library("ComplexHeatmap")
library("ggplot2")
library('pheatmap')
library("RNAseqQC")
library("ensembldb")
library("dplyr")
library("purrr")
library("magrittr")
library("EnhancedVolcano")
library("RColorBrewer")
library("AnnotationDbi")
library("tidyr")
library("tibble")
library("rnaseqGene")
library("stringr")
library("ggplot2")
library("apeglm")
library("tibble")
library("ggalt")
library("corrplot") 
library("corrr")
library("tidymodels")
library("janitor")
library("DataExplorer")
library("dplyr")
library("ggplot2")
library("igraph")
library("ggrepel")
library("skimr")
library("ggrepel")
library("DataExplorer")
```

# Part 1: Dataset Selection and Initial Analysis

## 1.1 Dataset Documentation

This RNA-seq dataset is derived from a project by the Cancer Genetics & Bioinformatics Laboratory at the International Laboratory for Human Genome Research. The data from this project have not been published yet, so we will not share very specific details.  

The data from this study come from an experiment involving four groups of rats, each consisting of 20 rats, for a total of 80 rats, subjected to different feeding protocols. Two of the groups were administered a hepatotoxic to induce the development of liver carcinoma, while the other two groups were kept as healthy controls. The main objective was to assess changes in gene expression related to the diets and the administration of the hepatotoxic.

Final distribution of the groups: 

- Group 1: healthy (normal diet) 
- Group 2: hepatotoxic (normal diet) 
- Group 3: healthy (treatment diet) 
- Group 4: hepatotoxic (treatment diet)

As initial data, we have the raw counts for each of the 80 samples, covering a total of 15,309 genes (features). In this project, our goal is to classify the experimental groups based on the GSVA scores for each pathway across the samples, capturing the variation in pathway activity. To achieve this, we will implement multinomial logistic regression models, Decision Trees, and Random Forests, aiming to identify the model with the best performance. 

## 1.2 Exploratory Data Analysis 

```{r, echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, fig.show = 'hide'}
#---- 01 Load data ----
# We load the raw counts stored within the dds object, which also contains the experimental design and other relevant data for the analysis.
dds <- readRDS("/home/user/Documents/Proyectod_ML/dds_object.rds")

#---- 02 normalization ----
dds <- estimateSizeFactors(dds) # Normalization of sample counts
vds<- vst(dds,blind = FALSE) # A variance stabilizing transformation (VST) is applied

# head(assay(vds), 3)
vst_counts <- assay(vds)
colnames(vst_counts) <- substr(colnames(vst_counts), 1, 2)
```

```{r, message = FALSE,}
#---- 03 Exploratory Data Analysis ----
# create_report(vst_counts) 
```

#### Data Quality and Missing Values

We generated a report using the DataExplorer package to gain a general overview of the data. Our dataset contains 80 samples and 15,309 genes, consisting solely of continuous variables. Furthermore, there are no missing values in any column or row of the dataset.   

#### Histograms

The report shows the distribution of gene expression for each sample. Most samples exhibit similar distributions, slightly skewed to the left, which is expected since the data has already been normalized. However, samples 61 and 62, belonging to group 4, display noticeably different distributions. According to the MultiQC analysis, these samples are not of poor quality but have a significantly lower number of reads compared to the other samples. This discrepancy could explain the globally lower gene expression observed in the histograms for these samples.

#### Correlation Analysis

Our dataset has been normalized using VST (Variance Stabilizing Transformation). This method stabilizes the variance of gene counts, reducing variability among samples. In our case, the correlation plot shows that the samples are highly positively correlated, likely as a result of the VST normalization process.


```{r, echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, fig.show = 'hide'}
#---- 04 PCA ----
plotPCA(vds, intgroup= c("condition"))
pca_data <- plotPCA(vds, intgroup = "condition", returnData = TRUE)
```


```{r, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE, fig.show = 'hide'}
# We add a variable that indicates the type of dieta 
pca_data <- pca_data %>% 
  dplyr::mutate(Diet = case_when(  condition=="AL" | condition=="AL_DEN" ~ "normal",
                                    condition=="TCR" | condition=="TCR_DEN" ~ "treatment")) 

#We changed the name of the groups
pca_data <- pca_data %>%
  dplyr::mutate(condition = case_when(
    condition == "AL" ~ "Group1",
    condition == "AL_DEN" ~ "Group2",
    condition == "TCR" ~ "Group3",
    condition == "TCR_DEN" ~ "Group4"
  ))
pca_data$name <- substr(pca_data$name, 1, 2)
```

```{r, warning = FALSE}
#Plot
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(x = PC1, y = PC2, col = Diet)) +
  geom_point(aes(shape = condition), size = 3) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 10) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_bw() +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.03) +
  ggtitle("PCA") + 
  scale_colour_manual(values = c("#009ACD", "#CD2626"))

```

#### PCA of gene expression 

We can see that the first two principal components explain 76% and 8% of the variance, respectively, totaling 84%. In the PCA, two well-defined groups emerge based on the type of diet: the first group corresponds to the normal diet (Blue), and the second group corresponds to the treatment diet (Red). This suggests that diet plays an important role in explaining the variation in our data and is likely influencing gene expression.

Both groups show dispersion, indicating that each group has internal variability and is not completely homogeneous. PC1 captures most of the variance, likely reflecting differences due to the diet type, while PC2 seems to capture more subtle differences between the diets. Additionally, we observe outliers in the treatment diet group, specifically samples 62 and 64, which were also identified in the DataExplorer report.

```{r}
#---- 05 Outliers ----
explore <- read.csv("/home/user/Documents/Proyectod_ML/explore.csv") # load the data
explore <-  explore %>% # drop the sample column 
  select( -X )

df_long <- explore %>% # Adjust the data set
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Calculate the quartiles and the lower and upper limits for each variable in the data set.
df_outliers <- df_long %>%
  group_by(variable) %>%
  mutate(
    Q1 = quantile(value, 0.25),  # First quartile
    Q3 = quantile(value, 0.75),  # Third quartile
    IQR = Q3 - Q1,               # Interquartile range
    lower_limit = Q1 - 1.5 * IQR,  # Lower limit
    upper_limit = Q3 + 1.5 * IQR,  # Upper limit
    is_outlier = value < lower_limit | value > upper_limit  # Identify outliers
  ) 

#Outliers plot 
ggplot(df_outliers, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  geom_point(data = filter(df_outliers, is_outlier), 
             aes(x = variable, y = value), 
             color = "red", size = 3) +  
  labs(title = "Outliers") +
  theme_minimal()

```

#### Outliers

Additionally, we examined our metadata to identify any variables that could be affecting samples 62 and 64. While these samples are outliers in terms of RNA Integrity Number (RIN), they are not outliers in the amount of RNA present, and they do not show poor quality in the FastQC report. However, they do have a significantly lower number of sequences compared to other samples. This places them as outliers, but it does not compromise the overall quality of the dataset.

We decided to keep these samples in the analysis because the machine learning models showed good performance (77% accuracy across all models). This suggests that, although samples 62 and 64 are outliers, they do not negatively impact the models' ability to make correct classifications. Furthermore, retaining these samples helps increase the diversity and representativeness of the dataset, which may improve the robustness and generalizability of our models.
 
### Gene Set Variation Analysis (GSVA)

Now we will perform a Gene Set Variation Analysis for our dataset. 

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#| label: load-packages
#| message: false
#| 
library(recount3)
library(GSVA)
library(GSVAdata)
library(msigdbr)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(org.Rn.eg.db) 
``` 


```{r}
#--------- 01 Res object ----------
expr_matrix <- vst_counts

# Creating sample metadata
col_data <- data.frame(
  sample_id = colnames(expr_matrix),
  condition = rep(c("control", "treated"), each = ncol(expr_matrix) / 2), 
  row.names = colnames(expr_matrix)
)

# Creating the SummarizedExperiment object
rse <- SummarizedExperiment(
  assays = list(counts = as.matrix(expr_matrix)),
  colData = col_data
)

# SummarizedExperiment object
# rse
```

```{r}
#--------- 02 Pathway-data ---------
#| label: pathway-data
# Load C2 canonical pathways for Rattus norvegicus
c2_gene_sets <- msigdbr(species = "Rattus norvegicus", 
                        category = "C2") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
print(paste("Number of pathways loaded:", length(c2_gene_sets)))
```

```{r}
#--------- 03 Data-cleaning ---------
#| label: clean-data

# Create the expression matrix from the data in `rse`
expr_matrix <- assays(rse)$counts 
gene_symbols <- rownames(rse)  
rownames(expr_matrix) <- gene_symbols

# Remove genes with duplicate names or `NA`
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ] 
# Remove duplicates
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]  

# Get the rse metadata
metadata <- colData(rse) %>% 
  as.data.frame()

print(paste("Número de genes después de limpiar:", nrow(expr_matrix)))
print(paste("Número de muestras:", ncol(expr_matrix)))

```
```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#--------- 04 Match gene names ---------
gene_ids <- rownames(expr_matrix)

# Mapping Ensembl id's to gene symbols
gene_names <- mapIds(org.Rn.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Number of genes without a symbol
# sum(is.na(gene_names))

# Replace Ensembl identifiers with gene names
rownames(expr_matrix) <- gene_names

# Delete rows with names NA
expr_matrix <- expr_matrix[!is.na(gene_names), ]

# Verificar el resultado
# head(rownames(expr_matrix))

```

```{r}
#--------- 05 Filter the expression matrix --------- 
# Obtaining the unique genes of the pathways
genes_in_pathways <- unique(unlist(c2_gene_sets))

# Filter the expression matrix
expr_matrix <- expr_matrix[rownames(expr_matrix) %in% genes_in_pathways, ]

print(paste("Number of genes after filtering:", nrow(expr_matrix)))

```

```{r}
#--------- 06 Filter duplicate genes --------- 
# Delete or merge duplicate genes 
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ] 

# Merge duplicate genes (taking the average of the expression values)
expr_matrix <- aggregate(expr_matrix, by = list(rownames(expr_matrix)), FUN = mean)

# Reassign gene names
rownames(expr_matrix) <- expr_matrix$Group.1  
# Remove the extra column from the groups
expr_matrix <- expr_matrix[, -1]  

# Check for duplicates
anyDuplicated(rownames(expr_matrix))  
# This should return 0 if there are no duplicates

```

```{r,  message = FALSE}
#--------- 07 Run gsva --------- 
#| label: run-gsva

# Set up GSVA parameters
expr_matrix <- as.matrix(expr_matrix)
gsvaPar <- GSVA::gsvaParam(expr = expr_matrix,
                           geneSets = c2_gene_sets,
                           kcdf = "Gaussian")  

# GSVA
gsva_results <- gsva(gsvaPar, verbose = TRUE)

# Results 
print(paste("Dimensions of GSVA results:", 
            nrow(gsva_results), "pathways by", 
            ncol(gsva_results), "samples"))

```


```{r}
#--------- 08 Basic-analysis --------- 
# Calculate the average enrichment scores for each pathway
mean_enrichment <- rowMeans(gsva_results)  

# Sort pathways descendingly by enrichment score
sorted_pathways <- sort(mean_enrichment, decreasing = TRUE)

# Select the top 20 pathways
top_pathways <- head(sorted_pathways, 20)

# Clean up pathway names 
pathway_names <- names(top_pathways)

# Eliminate common prefixes
shortened_names <- gsub("REACTOME_|KEGG_|BIOCARTA_|PID_|WP_", "", pathway_names)
shortened_names <- gsub("_", " ", shortened_names)

# Assign the new names to the pathways
names(top_pathways) <- shortened_names

# Print the 5 most enriched pathways
print("Top 5 enriched pathways:")
head(top_pathways, 5)
```
 
```{r, echo = TRUE, results = 'hide', message = FALSE, warning = FALSE, fig.show = 'hide'}
#--------- 09 Subset pathways --------- 
# Calculate average enrichment scores
mean_enrichment <- rowMeans(gsva_results)

# Select the first 100 most enriched routes
top_pathways <- sort(mean_enrichment, decreasing = TRUE)[1:25]

# Subset of GSVA results with the first 100 routes
gsva_subset <- gsva_results[names(top_pathways), ]

# Remove NA, NaN, or Inf values
gsva_subset_clean <- gsva_subset
gsva_subset_clean <- gsva_subset_clean[complete.cases(gsva_subset_clean), ]
gsva_subset_clean <- gsva_subset_clean[, apply(gsva_subset_clean, 2, function(x) all(is.finite(x)))]

# Heatmap 
pheatmap(gsva_subset_clean, 
         main = "Top 100 Enriched C2 Pathways", 
         scale = "row", 
         show_colnames = FALSE, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         clustering_method = "complete", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 4  
)


```
 
 
 
```{r}
#--------- 10 Clean final data --------- 
# Remove NA, NaN, or Inf values from the entire dataset
gsva_clean <- gsva_results
gsva_clean <- gsva_clean[complete.cases(gsva_clean), ]
gsva_clean <- gsva_clean[, apply(gsva_clean, 2, function(x) all(is.finite(x)))]

# Dimensions
print(paste("Número de filas y columnas después de la limpieza:", dim(gsva_clean)))

# write.csv
# write.csv(gsva_clean, "/home/user/Documents/Proyectod_ML/gsva_clean.csv", row.names = TRUE)

```
#### Final result 

We obtained a data set with 6347 pathways and kept our 80 samples after GSVA 

## Required Dimensionality Reduction Analysis

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#----01 package loading-----
library(corrplot)
library(tidymodels)
library(corrr)
library(palmerpenguins)
library(janitor)
library(tidyverse)
library(DataExplorer)
library(dplyr)
library(ggplot2)
```

```{r}
#----02 Load the data set -----
df <-  read.csv("/home/user/Documents/Proyectod_ML/gsva_clean.csv")
df <- df %>% rename(pathway = X)
df <- drop_na(df)

# Transpose the data frame
df <- t(df)
df <- as.data.frame(df)
colnames(df) <- as.character(df[1,])
df <- df[-1,]
df <- tibble::rownames_to_column(df, var = "sampleID")
df <- df %>%
  mutate(across(-sampleID, as.numeric))
# str(df)
# View(df)

# Load metadata to add groups 
df$group <- pca_data$condition
df$condition <- pca_data$Diet

#Add hepatotoxic
df <- df %>%
  mutate(hepatotoxic = case_when(
    group %in% c("Group1", "Group3") ~ "healthy",   
    group %in% c("Group2", "Group4") ~ "cancer"  
  ))

df <- df %>%
  select(sampleID, group, condition, hepatotoxic, everything()) 
# write.csv(df, file = "/home/user/Documents/Proyectod_ML/df.csv")
```

```{r}
#---- 03 Recipe ----
data_recipe <-
  recipe(~., data = df) %>%
  update_role(sampleID, group, condition, hepatotoxic, new_role = "id") %>%
  step_naomit(all_predictors(), id = "na") %>%   # Delete NA
  step_normalize(all_predictors(), id= "normalize") %>% # Normalize predictors
  step_pca(all_predictors(), num_comp = 4,  id = "pca") %>%  # Add PCA
  prep()

# Extract the PCA
data_pca <- data_recipe %>%
  tidy(id = "pca")
# data_pca

# Extract NAs IDs
data_na <- data_recipe %>%
  tidy(id = "na")
# data_na

# Extract Normalize IDs
data_nor <- data_recipe %>%
  tidy(id = "normalize")
# data_nor

```

#### Variance ratio

```{r}
#---- 03 Scree Plot ----
#Extract the variance
variance <- data_recipe %>%
  tidy(id = "pca", type = "variance")

# Filter for percent variance
percent_variance <- variance %>%
  filter(terms == "percent variance")

# Plot
ggplot(percent_variance, aes(x = component, y = value)) +
  geom_col() +
  ylab("Percent Variance Explained") +
  xlab("Principal Component") +
  ggtitle("Scree Plot") 

# We obtain the percentage of variance explained by each component
variance_table <- data_recipe %>%
  tidy(id = "pca", type = "variance") %>%
  filter(terms == "percent variance") %>%
  mutate(percent_variance = round(value, 2)) %>% 
  select(component, percent_variance)

#  Table
head(variance_table,10)

```

####  Number of components selected

We have 80 principal components, of which the first 4 explain 62.04% of the total variance. Starting from the fifth component, each accounts for less than 4% of the variance. Therefore, we decided to consider only the first 4 components and discard the remaining ones. 


```{r, warning=FALSE}
#---- 04 PCA ----
#Prep
pca_estimates <- prep(data_recipe)

#Bake 
pca_results <- bake(pca_estimates, new_data = NULL)

#PCA
pca_results$sampleID <- substr(pca_results$sampleID , 1, 3)
ggplot(pca_results, aes(PC1, PC2, color = hepatotoxic, shape = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = sampleID), vjust = -0.5, size = 3) + 
  stat_ellipse(aes(group = hepatotoxic), level = 0.95) +
  labs(title = "PC1 vs PC2", x = "PC 1", y = "PC 2") +
  theme_minimal() +
  theme(legend.position = "right")

```

#### PCA 1 and 2 

In this first PCA, we can observe that PC1 explains the largest portion of the data variability, accounting for 36.29%, while PC2 explains 13.61%, resulting in a total of 49.9% of the variability captured by these two components. Two clusters are moderately separated along PC1: the first cluster includes the groups that developed cancer due to the hepatotoxic agent (red), and the second cluster includes the healthy groups (blue).

Some samples from Group 4 (treatment diet + hepatotoxic agent) overlap with the healthy samples. This makes biological sense, as the treatment diet used (Time Caloric-Restriction Protocol) has been shown to have a protective effect in models of liver damage. Time-caloric restriction not only reduces hepatomegaly but also prevents the increase in blood leukocytes induced by hepatotoxic agents such as diethylnitrosamine. Furthermore, this protocol preserves the liver's functional and histological characteristics, favoring fibrotic areas over cirrhotic ones and limiting the progression to advanced stages of liver damage.

Although an increase in collagen deposits was observed, the protocol demonstrated positive effects, such as preventing systemic inflammation, reducing carcinoembryonic antigen levels, and limiting neoplastic transformation. Time Caloric-Restriction is an effective strategy to mitigate the adverse effects of liver damage and preserve liver function (Molina-Aguilar et al., 2017). This explains the behavior of the samples from Group 4, as we can observe expression patterns and pathways resembling healthy conditions in this group, despite the presence of the hepatotoxic agent.

```{r, warning=FALSE}
# Prep
pca_estimates <- prep(data_recipe)

# Bake 
pca_results <- bake(pca_estimates, new_data = NULL)

# Biplot para PC3 y PC4
pca_results$sampleID <- substr(pca_results$sampleID , 1, 3)
ggplot(pca_results, aes(PC3, PC4, color = condition, shape = group)) +
  geom_point(size = 3) +  
   geom_text(aes(label = sampleID), vjust = -0.5, size = 3) +
  stat_ellipse(aes(group = condition), level = 0.95) +
  labs(title = "PC3 vs PC4", x = "PC 3", y = "PC 4") +
  theme_minimal() +
  theme(legend.position = "right")  # Coloca la leyenda a la derecha
```
#### PCA 3 and 4 

In this PCA, we are analyzing PC3, which captures 6.93% of the variability, and PC4, which captures 5.21%. Together, they account for 12.14% of the total variability. There is a clear differentiation between the samples based on diet—normal (red) and treatment (blue)—especially along PC4. This indicates that we are also observing differences in the pathways expressed as a result of the type of diet.

We have the following:

- PC1 and PC2 capture the global differences between healthy samples and cancer samples.
- PC3 and PC4 reflect diet-specific effects, distinguishing between treatment and normal diets.

```{r}
#---- 05 Importance plot ----

################################## PC1 ########################################
# Get PCA loadings from recipe 
loadings <- data_recipe %>%
  tidy(id = "pca", type = "coef")

# Filter for only one componente 
loadings_pc1 <- loadings %>%
  filter(component == "PC1") %>%
  arrange(desc(value))  

# Select the 10 variables with the greatest positive and negative contribution
top_contrib_positive <- loadings_pc1 %>%
  top_n(10, value)  

top_contrib_negative <- loadings_pc1 %>%
  top_n(-10, value) 

top_contrib <- bind_rows(top_contrib_positive, top_contrib_negative)

# Importance plot
ggplot(top_contrib, aes(x = reorder(terms, value), y = value)) +
  geom_col(fill = "steelblue") +
  coord_flip() +  
  xlab("Variables") +
  ylab("Contribution to PC1") +
  ggtitle("Top 10 contributing variables to PC1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),  
    axis.text.y = element_text(size = 6),  
    axis.title = element_text(size = 12),   
    plot.title = element_text(size = 14)   
  )

################################## PC2 ########################################

loadings <- data_recipe %>%
  tidy(id = "pca", type = "coef")

loadings_pc2 <- loadings %>%
  filter(component == "PC2") %>%
  arrange(desc(value)) 

top_contrib_positive <- loadings_pc2 %>%
  top_n(10, value)  

top_contrib_negative <- loadings_pc2 %>%
  top_n(-10, value) 

top_contrib <- bind_rows(top_contrib_positive, top_contrib_negative)

ggplot(top_contrib, aes(x = reorder(terms, value), y = value)) +
  geom_col(fill = "pink") +
  coord_flip() +  
  xlab("Variables") +
  ylab("Contribution to PC2") +
  ggtitle("Top 10 contributing variables to PC2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),  
    axis.text.y = element_text(size = 6),  
    axis.title = element_text(size = 8),   
    plot.title = element_text(size = 8)    
  )

################################## PC3 ########################################
loadings <- data_recipe %>%
  tidy(id = "pca", type = "coef")

loadings_pc3 <- loadings %>%
  filter(component == "PC3") %>%
  arrange(desc(value))  

top_contrib_positive <- loadings_pc3 %>%
  top_n(10, value) 

top_contrib_negative <- loadings_pc3 %>%
  top_n(-10, value)  

top_contrib <- bind_rows(top_contrib_positive, top_contrib_negative)

ggplot(top_contrib, aes(x = reorder(terms, value), y = value)) +
  geom_col(fill = "#B23AEE") +
  coord_flip() +  
  xlab("Variables") +
  ylab("Contribution to PC3") +
  ggtitle("Top 10 contributing variables to PC3") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),  
    axis.text.y = element_text(size = 6),  
    axis.title = element_text(size = 8),   
    plot.title = element_text(size = 8)    
  )

################################## PC4 ########################################
loadings <- data_recipe %>%
  tidy(id = "pca", type = "coef")

loadings_pc4 <- loadings %>%
  filter(component == "PC4") %>%
  arrange(desc(value))  

top_contrib_positive <- loadings_pc4 %>%
  top_n(10, value) 

top_contrib_negative <- loadings_pc4 %>%
  top_n(-10, value)  

top_contrib <- bind_rows(top_contrib_positive, top_contrib_negative)

ggplot(top_contrib, aes(x = reorder(terms, value), y = value)) +
  geom_col(fill = "#00CDCD") +
  coord_flip() + 
  xlab("Variables") +
  ylab("Contribution to PC4") +
  ggtitle("PC4 Top 10 contributing variables to PC4") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),  
    axis.text.y = element_text(size = 6),  
    axis.title = element_text(size = 8),   
    plot.title = element_text(size = 8)    
  )
```

#### PC1

This principal component captures:

- Positive Contributions: Processes related to cell differentiation, signaling pathways (e.g., EGFR signaling, focal adhesion), and cell cycle regulation (e.g., Rho GTPase cycle, CDC42 GTPase cycle).

- Negative Contributions: Processes related to mitochondrial function, polyamine metabolism, and degradation pathways, such as SCF-SKP2-mediated degradation.

PC1 separates the samples based on cancer status, distinguishing cancerous samples by processes related to cell proliferation, signaling, and differentiation, while healthy samples are characterized by mitochondrial and metabolic processes. This suggests that PC1 helps differentiate between healthy and cancer samples based on their gene expression profiles.

#### PC2

This principal component captures:

- Positive Contributions: Processes related to liver cancer, peripheral clock regulation, zinc homeostasis, and chylomicron clearance, linked to metabolism and liver function.

- Negative Contributions: Cancer resistance mechanisms, such as cisplatin resistance and various cancer subtypes, indicating downregulated processes in cancer.

PC2 separates samples based on cancer status, with cancer samples showing disrupted metabolism, DNA damage response, and treatment resistance, while healthy samples are associated with normal metabolic and cellular functions.

#### PC3 

This principal component captures:

-Positive Contributions: Processes related to obesity, myeloid cell development, inflammatory response, and signaling by interleukins, which are linked to immune response and cell differentiation.

- Negative Contributions: Processes involving negative regulation of signaling pathways, such as phospholipase C-mediated cascade and TCF-dependent signaling, indicating potential downregulation in certain cancer-related pathways.

PC3 primarily separates samples based on diet, with distinct clustering based on dietary factors rather than cancer status, highlighting differences in immune and metabolic responses.

#### PC4

This principal component captures:

- Positive Contributions: Processes related to cell proliferation, cancer, and cell cycle regulation.
- Negative Contributions: Metabolic processes, antioxidant responses, and cellular stress.

PC4 separates the samples based on diet type, with a distinction between healthy and treatment diet groups. This component helps to differentiate between samples based on the influence of the diet.

#### Variable selection 

We decided to select the 25 most relevant variables from each of the 4 principal components to build our final dataset, which we will use to develop our machine learning models. This will help us reduce the dimensionality of the data while retaining the most significant features that explain the variability within it. 

```{r}
#---- 06 Variable selection ----
# store the most important variables of each component
important_vars <- list()

# Filter contributions for the 4 PCs
for (pc in paste0("PC", 1:4)) {
  loadings_pc <- loadings %>%
    filter(component == pc) %>%
    arrange(desc(abs(value)))  # # Sort by magnitude of contribution 
  # Select the 25 variables with the greatest contribution
  top_vars_pc <- loadings_pc %>%
    slice_head(n = 25) %>%
    pull(terms)  
  
  # list of variables
  important_vars[[pc]] <- top_vars_pc
}

# Merge all selected variables into a single vector 
selected_vars <- unique(unlist(important_vars))
extra_vars <- c("sampleID",  "group", "condition", "hepatotoxic")

# Combine both vectors
selected_vars <- c(selected_vars, extra_vars)
selected_vars <- unique(selected_vars)
```

```{r}
# Subset of the original dataset with the selected variables
reduced_data <- df[, selected_vars]
reduced_data <- reduced_data %>%
  select(sampleID, group, condition, hepatotoxic, everything())

# write.csv(reduced_data, file = "/home/user/Documents/Proyectod_ML/reduced_data.csv")
```

## 1.3 Biological Relevance

This dataset and analysis aim to advance hepatocellular carcinoma (HCC) research. By using GSVA to classify samples based on gene set activity, we can identify patterns of gene expression associated with specific diets and key pathways involved in liver disease progression. Notably, the TCR protocol's ability to slow fibrosis progression (Molina-Aguilar et al., 2017) provides a valuable opportunity to explore early stages of liver disease and improve clinical outcomes.

The findings from this analysis could have substantial clinical and research implications. By applying machine learning models to classify patient samples based on pathway signatures, we aim to develop predictive tools that could expand into broader clinical applications, such as early diagnosis and personalized treatment strategies. However, some limitations should be considered. The sample size, while informative, may limit generalizability, and the GSVA method captures only part of the complexity of liver disease biology.
