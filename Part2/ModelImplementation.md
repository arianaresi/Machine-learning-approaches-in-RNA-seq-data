# **Part 2: Machine Leaning Approach**

## **2.1 Problem Formulation**

This project aims to classify rats based on their feeding schedules and the impact of hepatotoxin administration, using pathway enrichment scores derived from RNA sequencing (RNAseq) data. These scores were computed using Gene Set Variation Analysis (GSVA), which summarizes the activity of biological pathways instead of focusing on individual gene expression.

By leveraging *machine learning* models, the study seeks to identify patterns in pathway activity associated with feeding schedules, while accounting for the effects of hepatotoxicity.

This project employs a *supervised* learning approach, as the dataset includes labeled samples that indicate the feeding schedule and treatment condition for each rat. The objective is to develop a model capable of accurately classifying unknown rats into their respective feeding schedule categories based on pathway activity profiles.

The performance of the model will be *evaluated* using standard classification metrics such as accuracy, ROC-AUC, and Kappa score. By focusing on pathway-level data, this approach aims to provide a more functional and biologically relevant interpretation of the relationships between feeding schedules, hepatotoxin exposure, and their effects on cellular pathways.

## **2.2 Required Model Implementation**

## Multi-class classification

We implemented a multinomial logistic regression model to classify the four experimental groups in our dataset, using the 25 most relevant variables from each of the 4 principal components. This approach allowed us to reduce the dataset's dimensionality and focus the model on the most representative features, maximizing its predictive capacity while minimizing the noise from less informative variables.

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#----01 package loading-----
library(tidyverse)
library(tidymodels)
```

We loaded the dataset selected from the principal component analysis. Our dataset contains 80 samples and 105 pathways, with no missing data. Each group contains 20 samples. In terms of diet, two of these groups follow the treatment diet, while the other two follow a normal diet. Additionally, there are two groups with cancer and two healthy groups. Therefore, there is no bias in the number of samples across the different categories.

```{r}
#----01 Load the data -----
df <- read.csv("/home/user/Documents/Proyectod_ML/reduced_data.csv")
# dim(df)

#----02 Clean the data ----
rattus_clean <- df %>%
  drop_na()

rattus_clean <- rattus_clean %>%
  select(-X )
rattus_clean$group <- as.factor(rattus_clean$group)
rattus_clean$condition <- as.factor(rattus_clean$condition)
rattus_clean$hepatotoxic <- as.factor(rattus_clean$hepatotoxic)
rattus_clean <- rattus_clean %>%
  select(-condition, -hepatotoxic, -sampleID )
# str(rattus_clean)
# View(rattus_clean)
```


```{r}
#----03 Logistic Regression ----
# Split the data into training and testing
set.seed(123)  
rattus_split <- initial_split(rattus_clean, prop = 0.8, strata = group)
rattus_train <- training(rattus_split)
rattus_test <- testing(rattus_split)
```

```{r}
#----04 Specify the model ----
log_reg_model <- multinom_reg() %>%
  set_engine("nnet")  # We use the 'nnet' engine for multinomial logistic regression

```

```{r}
#----05 Recipe -----
rattus_recipe <- recipe(group ~ ., data = rattus_train) %>%
  step_normalize(all_numeric_predictors()) %>%  
  step_dummy(all_nominal_predictors())  # Creating dummy variables for categorical variables

```

```{r}
#----06 workflow -----
rattus_workflow <- workflow() %>%
  add_model(log_reg_model) %>%
  add_recipe(rattus_recipe)
```

```{r}
#----07 Fit the model -----
rattus_fit <- rattus_workflow %>%
  fit(data = rattus_train)

```

```{r}
#----08 Predict on test set -----
predictions <- rattus_fit %>%
  predict(new_data = rattus_test) %>%
  bind_cols(rattus_test)

```

```{r}
#----09 Performance -----
metrics <- predictions %>%
  metrics(truth = group, estimate = .pred_class)
# Metrics
metrics
```

These two metrics indicate very good model performance. The 87.5% accuracy means that the model correctly assigned the samples to their respective groups 87.5% of the time in the test set. On the other hand, the Kappa index of 0.833 shows that, in addition to accuracy, the model has high consistency in its predictions.

```{r, message = FALSE}
#----10 Confusion Matrix -----
conf_mat <- predictions %>%
  conf_mat(truth = group, estimate = .pred_class)

autoplot(conf_mat, type = "heatmap") + 
  scale_fill_gradient(low = "white", high = "pink")

```

We can observe that we have a high number of true positives and very few false positives and false negatives, suggesting that the model is correctly classifying most of the samples and making few errors.

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#----11 Extraction of parameters -----
# Coefficients for the multinomial model
model_fit <- extract_fit_parsnip(rattus_fit)

# Fitted values
model_fit$fit$fitted.values
coef(model_fit$fit)

# Extract the coefficients from the multinomial model
coefs <- coef(model_fit$fit)
coefs_df <- as.data.frame(coefs)
coefs_df$class <- rownames(coefs_df)  
coefs_long <- coefs_df %>%
  pivot_longer(cols = -class, names_to = "variable", values_to = "coefficient")

```

```{r}
#----12 Importance plot -----
# Sort the coefficients by absolute value and select the 15 most important variables
top_15_vars <- coefs_long %>%
  mutate(abs_coefficient = abs(coefficient)) %>%
  group_by(variable) %>%
  summarise(mean_abs_coefficient = mean(abs_coefficient)) %>%
  arrange(desc(mean_abs_coefficient)) %>%
  head(20)

# Filter the original dataframe to keep only the 15 most important variables
top_15_coefs <- coefs_long %>%
  filter(variable %in% top_15_vars$variable)

# Importance plot 
ggplot(top_15_coefs, aes(x = reorder(variable, abs(coefficient)), y = coefficient, fill = class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 15 Most Important Variables in the Multinomial Model",
       x = "Variables",
       y = "Coefficient",
       fill = "Clase") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =6 ),  
        axis.text.y = element_text(size = 6)) +  
  scale_x_discrete(labels = function(x) gsub("_", " ", x))  

```

The pathways that contribute the most to the classifier are mainly related to lipid metabolism, treatment response, cell cycle, antioxidant defense, and immune signaling. These are key areas in cancer development and progression. This suggests that the model is capturing relevant information across the different groups.

#### Regularization

We used Elastic Net as a regularization method, as it combines Lasso (L1) and Ridge (L2). Given that we have 100 pathways and 80 samples, this method is well-suited for handling high-dimensional situations, as it helps select relevant features while controlling for overfitting. By combining the properties of Lasso and Ridge, Elastic Net can improve model accuracy and manage multicollinearity between variables, which is crucial when working with a large number of features relative to the number of available samples.
```{r}
#----12 Regularization -----
elastic_net_spec <- multinom_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet")

```

```{r}
#----13 cross validation -----
rattus_folds <- vfold_cv(rattus_clean, v = 5, strata = group)
```

```{r}
#----14 hyperparameters -----

# Define the hyperparameter grid:
elastic_net_grid <- grid_regular(
  penalty(range = c(-4, 1)),  
  mixture(range = c(0, 1)),   # Mixture between Lasso (L1) and Ridge (L2)
  levels = 5                  
)

#  workflow:
rattus_recipe <- recipe(group ~ ., data = rattus_train) %>%
  step_normalize(all_numeric_predictors()) %>%  
  step_dummy(all_nominal_predictors())  

rattus_workflow <- workflow() %>%
  add_model(elastic_net_spec) %>%
  add_recipe(rattus_recipe)

#Tune_grid
elastic_net_tune <- tune_grid(
  rattus_workflow,  
  resamples = rattus_folds, 
  grid = elastic_net_grid  
)

```

```{r, echo = TRUE, results = 'hide', message = FALSE, warning=FALSE}
#----15 metrics -----
elastic_net_tune %>%
  collect_metrics()

```


```{r}
autoplot(elastic_net_tune) +
  ggtitle("Tuning Results for Elastic Net Model") +
  theme_minimal()
```

With low or moderate regularization, the performance metrics (accuracy, brier_class, roc_auc) are good, and the Lasso proportion has no significant impact. In contrast, with high regularization , all metrics deteriorate, especially with higher Lasso proportions, indicating an overly simplistic model. The best results are observed with intermediate Lasso proportions (between 0.25 and 0.75), achieving a good balance across metrics.


#### Best model

```{r}
# Select the best model based on the roc_auc metric
best_elastic_net <- select_best(elastic_net_tune, metric = "roc_auc")
print(best_elastic_net)
```
The model uses weak regularization (penalty = 1e-04) with 25% Lasso (L1) and 75% Ridge (L2). This means it focuses more on smoothing coefficients like Ridge but still allows some coefficients to be zero thanks to the Lasso component. 

#### Comparación de modelos
```{r}
#----16 Model comparison -----
# Tune the model with the best hyperparameters
final_workflow <- finalize_workflow(rattus_workflow, best_elastic_net)

# Train the model with the full training set
elastic_net_fit <- fit(final_workflow, data = rattus_train)

# Function to evaluate the model and calculate only general metrics
evaluate_model <- function(model_fit) {
  predictions <- augment(model_fit, new_data = rattus_test)
  metrics <- predictions %>%
    metrics(truth = group, estimate = .pred_class)
  
  return(metrics)
}

# Evaluate the Elastic Net model
elastic_net_results <- evaluate_model(elastic_net_fit)

# metrics 
elastic_net_results
```

The model with all variables outperforms the Elastic Net in accuracy (87.5% vs. 75%) and Kappa index (83.33% vs. 66.67%), indicating that the full model correctly predicts more cases and has greater agreement between the predictions and the true classes.

The Elastic Net model performs worse in both metrics compared to the model with all variables. This is likely due to Elastic Net introducing regularization, penalizing the coefficients, and reducing the model's complexity.





```{r}
#Importance plot 
coefs <- tidy(elastic_net_fit)

coefs$class <- rep(c("Group1", "Group2", "Group3", "Group4"), each = nrow(coefs) / 4)

top_15_vars <- coefs %>%
  mutate(abs_coefficient = abs(estimate)) %>%
  arrange(desc(abs_coefficient)) %>%
  head(15)


top_15_coefs <- coefs %>%
  filter(term %in% top_15_vars$term)


library(ggplot2)

ggplot(top_15_coefs, aes(x = reorder(term, abs(estimate)), y = estimate, fill = class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 15 Variables ",
       x = "Variables",
       y = "Coeficiente",
       fill = "Clase") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  
        axis.text.y = element_text(size = 6),  
        plot.title = element_text(size = 10),  
        axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8)) +  
  scale_x_discrete(labels = function(x) gsub("_", " ", x))
```

The biological pathways identified in the elastic net model are linked to key cellular and metabolic processes that affect health and treatment response. For example, the REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION pathway is involved in antigen presentation to T cells, crucial for adaptive immune response, suggesting a role in defense against infections and cancer. Other pathways, such as REACTOME_NR1H2_NR1H3_REGULATE_GENE_EXPRESSION_TO_LIMIT_CHOLESTEROL_UPTAKE, regulate cholesterol uptake, which has implications for metabolic and cardiovascular diseases. The catalytic cycle of flavin-containing monooxygenases pathway is associated with drug metabolism and cellular detoxification, while the steroidogenesis pathway with glucocorticoid and mineralocorticoid receptors is related to inflammation regulation and stress response. Together, these pathways reflect a network of biological processes that impact cellular homeostasis, immune response, metabolism, and hormonal regulation, all of which are crucial for various clinical conditions.

### Decision trees and Random forest models

For the project, we used decision tree and random forest models to predict the categories of rats evaluated in the study (group 1, group 2, group 3 and group 4), each with specific treatment and diet schedules. For further explanation, some explanations of each of the models will be included below.

A decision tree is a non-parametric supervised learning algorithm, which is utilized for both classification and regression tasks. It has a hierarchical, tree structure, which consists of a root node, branches, internal nodes and leaf nodes (Ibm, 2024).

Meanwhile, random forest is a machine learning algorithm that combines the output of multiple decision trees to reach a single result. It handles both classification and regression problems (Ibm, 2024).

```{r, message = FALSE, warning = FALSE}
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)
library(corrplot) 
library(corrr)
library(tidymodels)
library(janitor)
library(DataExplorer)
library(tidyverse)
library(tidymodels)
library(vip)
library(viridis)

df <- read.csv("/Users/arianasilva/Documents/ML Renan/Proyecto/reduced_data.csv")
dim(df)
```

We have included the dataset containing the pathways, groups, conditions, and treatments. Given that conditions and treatments are integral components of the experimental design for group formation, we hypothesize that they may be strongly correlated with our target variable. Nevertheless, we will conduct a Pearson's Chi-squared test to assess this relationship.

```{r}
# Create contigency table between condition and hepatotoxic, condition and group, hepatotoxic and group
table_condition_hepatotoxic <- table(df$condition, df$hepatotoxic)

table_condition_group <- table(df$condition, df$group)

table_hepatotoxic_group <- table(df$hepatotoxic, df$group)

# Chi-squeared test
chi_sq_condition_hepatotoxic <- chisq.test(table_condition_hepatotoxic)

chi_sq_condition_group <- chisq.test(table_condition_group)

chi_sq_hepatotoxic_group <- chisq.test(table_hepatotoxic_group)

# Print results
print(chi_sq_condition_hepatotoxic)

print(chi_sq_condition_group)

print(chi_sq_hepatotoxic_group)
```

**Condition vs. Hepatotoxic:** The result yields an X-squared value of 0, with a p-value of 1; this thereby indicates no significant association between the condition and hepatotoxic variable. This indicates that these two variables are independent and there is no meaningful relationship between these variables in the dataset.

**Condition vs. Group:** The test gives an X-squared value of 80 with a p-value smaller than 2.2e-16, which is highly significant. This therefore indicates very strong association between the condition and group variables; hence, these two variables are not independent and may be related.

**Hepatotoxic vs. Group:** In the test, X-squared value of 80 was given with a p-value smaller than 2.2e-16, also indicating strong and statistically significant association of hepatotoxic with group.

In other words, while there is no significant relationship between condition and hepatotoxic, both condition and hepatotoxic are strongly and significantly related to the group variable. Therefore, we will eliminate both the conditional variable and the hepatotoxic variable for the creation of the model.

Then, we performed data preprocessing on the `rattus_clean` dataset by cleaning the data (removing missing values and irrelevant columns) and normalizing predictors. The goal is to prepare the data for further analysis or modeling.

```{r}
rattus_clean <- df %>%
  drop_na()

# Check the number of remaining rows
nrow(rattus_clean)
rattus_clean <- rattus_clean %>%
  select(-c(condition, hepatotoxic, X, sampleID))
rattus_clean$group <- as.factor(rattus_clean$group)
View(rattus_clean)

# Recipe
data_recipe <- recipe(~., data = rattus_clean) %>%
  update_role(group, new_role = "id") %>%
  step_naomit(all_predictors(), id = "na") %>%   # Delete NA
  step_normalize(all_predictors(), id = "normalize") %>% # Normalize predictors
  prep()
```

### Splitting data

The dataset is split into training and testing subsets using initial_split, ensuring the stratification is done by the *group* variable, representing the different groups of rats based on their feeding schedule.

```{r}
# Data split
set.seed(123)
subset_split <- initial_split(rattus_clean, prop = 0.7, strata = group)
subset_train <- training(subset_split)
subset_test <- testing(subset_split)
```

### Recipe for classification

A preprocessing recipe was created, specifying how the data should be transformed before being used in the model. The target variable (group) is considered an identifier, with all other variables being predictors, involving dummy encoding for categorical variables and normalization for numerical variables. 

```{r}
# Recipe for the mode
subset_recipe <- recipe(group ~ ., data = subset_train) %>%
  step_dummy(all_nominal_predictors()) %>%  # Convert categorical features into dummy variables (binary)
  step_zv(all_predictors()) # Remove all predictors with zero variance
```

### Cross-validation

10-fold cross-validation is used to evaluate the model's performance, where the data is split into 10 subsets (folds), and the model is trained on 9 folds, with the remaining fold used as the test set. 

```{r}
set.seed(123)
subset_folds <- vfold_cv(subset_train, v = 10, repeats = 5, strata = group)
```

### Classification

This is the moment when the decision tree and random forest models begin to be built, with the subsequent creation of workflows where the recipe created previously is combined with the specific model.

```{r}
library(ranger)
# Specifying decision tree model
dt_spec <- decision_tree() %>%
  set_mode("classification") %>% # Set the mode to sorting
  set_engine("rpart") # "rpart" is tipically used in decision tree models

# Specifying random forest model
rf_spec <- rand_forest() %>%
  set_mode("classification") %>%  # Set the mode to sorting
  set_engine("ranger", importance = "impurity") # "ranger" is an eficient engine for random forest

dt_workflow <- workflow() %>%
  add_model(dt_spec) %>%
  add_recipe(subset_recipe)

rf_workflow <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(subset_recipe)
```

The Random Forest model is trained and evaluated using the cross-validation folds. The workflow combines the preprocessing recipe and the model, and the fit_resamples() function trains the model on each fold while collecting performance metrics.

### Model evaluation

Performance metrics from the cross-validation process are aggregated using collect_metrics.

```{r}
dt_res <- dt_workflow %>%
  fit_resamples(
    resamples = subset_folds, 
    metrics = metric_set(accuracy, roc_auc, kap),
    control = control_resamples(save_pred = TRUE)
  )

rf_res <- rf_workflow %>%
  fit_resamples(
    resamples = subset_folds,
    metrics = metric_set(accuracy, roc_auc, kap),
    control = control_resamples(save_pred = TRUE)
  )

# Decision tree metrics
dt_metrics <- dt_res %>%
  collect_metrics()
# Random forest metrics
rf_metrics <- rf_res %>%
  collect_metrics()

# Combine both metrics
combined_metrics <- bind_rows(
  dt_metrics %>% mutate(model = "Decision Tree"),
  rf_metrics %>% mutate(model = "Random Forest")
)

print(combined_metrics)

ggplot(combined_metrics, aes(x = model, y = mean, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  facet_wrap(~.metric, scales = "free_y") + 
  labs(title = "Metrics comparation between models", 
       x = "Model", 
       y = "Mean metric") +
  theme_minimal() +
  scale_fill_manual(values = c("Decision Tree" = "#CDAA7D", "Random Forest" = "lightblue3"))

```

**Accuracy:** In decision tree model, on average, it correctly classified 70% of the cases, in contrast with random forest model, which correctly classified 80% of the cases, better than decision tree model.

**Kappa score:** Kappa score of 0.6111 for the decision tree model indicates moderate agreement between predictions and actual values, while the score of 0.7444 for the random forest model suggests better agreement.

**ROC AUC:** Since the ROC AUC of the decision tree model is 0.8842 and that of the random forest is 0.9481, we can conclude that the random forest performs significantly better in classifying this data set.

Subsequently, we fitted a decision tree model to the training data, visualized the feature importance, and extracted the most influential variables. In a biological context, these key features correspond to the pathways most impacted by the model, helping to identify which biological processes play a crucial role in the observed outcomes.

```{r}
# Visualizar la importancia de características para Decision Tree
modelo_dt <- dt_workflow %>%
  fit(subset_train) %>%
  extract_fit_parsnip()

vip(modelo_dt, geom = "col", aesthetics = list(fill = "#CDAA7D")) +
  labs(title = "Feature Importance - Decision Tree",
       x = "Features",
       y = "Importance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # Reduce size and rotate labels
    axis.title = element_text(size = 8), # Adjust axis title size
    plot.title = element_text(size = 8), # Adjust plot title size
    plot.margin = margin(8, 8, 8, 8), # Adjust margins to save space
    axis.text.y = element_text(size = 6), # Reduce size of y-axis labels (importance values)
    axis.ticks = element_line(size = 0.5), # Reduce size of axis ticks
    legend.position = "none" # Hide legend if not necessary
  ) +
  coord_flip() # Flip coordinates if the variable names are still too large
```
In the decision tree model, we observe that the most influential pathways are associated with biological processes such as drug metabolism, gene silencing through methylation, pathways regulated by the WT1 gene, which plays key roles in development and cell survival, as well as immune system responses and the metabolism of specific substances.

Similarly, the same approach was applied to the random forest model. By fitting the model to the training data and extracting the feature importance, we identified the pathways most influential in the random forest model. 

```{r}
# Visualizar la importancia de características para Decision Tree
modelo_rf <- rf_workflow %>%
  fit(subset_train) %>%
  extract_fit_parsnip()

vip(modelo_rf, geom = "col", aesthetics = list(fill = "lightblue3")) +
  labs(title = "Feature Importance - RF",
       x = "Features",
       y = "Importance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # Reduce size and rotate labels
    axis.title = element_text(size = 8), # Adjust axis title size
    plot.title = element_text(size = 8), # Adjust plot title size
    plot.margin = margin(8, 8, 8, 8), # Adjust margins to save space
    axis.text.y = element_text(size = 6), # Reduce size of y-axis labels (importance values)
    axis.ticks = element_line(size = 0.5), # Reduce size of axis ticks
    legend.position = "none" # Hide legend if not necessary
  ) +
  coord_flip() # Flip coordinates if the variable names are still too large
```
In the random forest model, we observe pathways that regulate hormones and metabolic functions, adipogenesis, gene silencing, apoptosis, and liver toxicity. Additionally, pathways related to nutrient absorption and immune system functions are also present.

### **Model prediction**

After the model has been trained, it is utilized to predict outcomes on the test set. The accuracy of the model is assessed by comparing its predictions with the actual labels of the test set to determine its ability to work with new data. 

```{r}
# Fit the best performing model to the entire training set
dt_best_model <- dt_workflow
dt_final_fit <- dt_best_model %>%
  last_fit(split = subset_split)

rf_best_model <- rf_workflow
rf_final_fit <- rf_best_model %>%
  last_fit(split = subset_split)

# Collect the results
dt_final_results <- collect_metrics(dt_final_fit)
rf_final_results <- collect_metrics(rf_final_fit)

print(dt_final_results)
print(rf_final_results)
```

**Accuracy:** the decision tree model correctly classified, on average, 50% of the cases, while the random forest model demonstrated a higher accuracy, correctly classifying 79% of the cases.

**ROC AUC:** the decision tree model achieved a ROC AUC of 0.7106, compared to 0.9481 for the random forest model. This indicates that the random forest significantly outperforms the decision tree in terms of classification ability, as reflected by the higher ROC AUC value.

**Brier class:** the decision tree model has a score of 0.4197 indicates the model's predicted probabilities are fairly accurate, however, there is still potential for enhancement. Meanwhile, in random forest the brier class value is 0.1987, which shows good probability calibration for multiclass predictions.

Then, a confusion matrix was created to assess how well the classification model performed by comparing the predicted and actual class labels in the test set. This matrix shows how well the model categorizes each group, detecting true positives, false positives, true negatives, and false negatives.

```{r}
dt_predictions <- dt_final_fit %>%
  collect_predictions()
rf_predictions <- rf_final_fit %>%
  collect_predictions()

dt_conf_mat <- dt_predictions %>%
  conf_mat(truth = group, estimate = .pred_class)
rf_conf_mat <- rf_predictions %>%
  conf_mat(truth = group, estimate = .pred_class)

autoplot(dt_conf_mat, type = "heatmap")
```

```{r}
autoplot(rf_conf_mat, type = "heatmap")
```

While *decision tree* showed decent performance, especially in group 2, but struggled more with misclassifying instances from group 1 and group 4, *random forest* demonstrated better overall performance, particularly in group 4, with fewer misclassifications across all groups. It showed a balanced prediction approach, correctly identifying instances in group 1, group 3, and group 4.

### **Hyperparameter adjustment**

Also, we tuned the hyperparameters of a Random Forest model by fine-tuning parameters like mtry, number of trees (trees), and minimum data points for a split (min_n). It establishes a range of potential values for every parameter and employs cross-validation to assess the model's effectiveness. The optimal hyperparameter combination is chosen considering accuracy, and the test set is used to train and assess the final model. The model's performance is optimized by visualizing the results with a confusion matrix.

```{r}
# Especificar modelo de Random Forest con hiperparámetros ajustables
rf_spec_tuned <- rand_forest(
  mtry = tune(),       # Ajustar el número de predictores seleccionados en cada división
  trees = tune(),
  min_n = tune()       # Ajustar el número mínimo de datos para realizar una división
) %>%
  set_mode("classification") %>%
  set_engine("ranger", importance = "impurity")

# Actualizar el workflow con el modelo que tiene tune()
rf_workflow_tuned <- workflow() %>%
  add_model(rf_spec_tuned) %>%
  add_recipe(subset_recipe)

# Definir la cuadrícula de hiperparámetros
grid_rf <- grid_regular(
  mtry(range = c(1, 10)), # Rango para 'mtry'
  min_n(range = c(2, 20)), 
  trees(range = c(100, 800)),# Rango para 'min_n'
  levels = 5                    # Número de valores a probar en cada rango
)

# Ajuste del modelo usando tune_grid
set.seed(123) # Fijar semilla para reproducibilidad
rf_tuned_res <- rf_workflow_tuned %>%
  tune_grid(
    resamples = subset_folds,     # Validación cruzada definida previamente
    grid = grid_rf,               # Cuadrícula de hiperparámetros
    metrics = metric_set(accuracy, roc_auc, kap),  # Métricas a evaluar
    control = control_grid(save_pred = TRUE)  # Guardar predicciones
  )

# Seleccionar los mejores hiperparámetros según la métrica deseada (por ejemplo, accuracy)
best_rf_params <- select_best(rf_tuned_res, metric = "accuracy")

# Finalizar el workflow con los mejores hiperparámetros
final_rf_workflow <- finalize_workflow(
  rf_workflow_tuned,
  best_rf_params
)

# Entrenar el modelo final usando los mejores hiperparámetros
final_train_fit <- final_rf_workflow %>%
  fit(data = subset_train)

# Evaluar el modelo en el conjunto de prueba
final_test_fit <- final_rf_workflow %>%
  last_fit(split = subset_split)

# Recoger los resultados finales
final_results <- collect_metrics(final_test_fit)
print(rf_final_results)

# Obtener las predicciones en el conjunto de prueba
final_predictions <- final_test_fit %>%
  collect_predictions()  # Recoge las predicciones de la última evaluación

# Crear la matriz de confusión
final_conf_mat <- final_predictions %>%
  conf_mat(truth = group, estimate = .pred_class)

# Visualizar la matriz de confusión
autoplot(final_conf_mat, type = "heatmap")
```


Talking about metrics, we have accuracy, ROC AUC and brier class.

**Accuracy:** The accuracy of the initial random forest model was 79%, while the hyperparameter-tuned model achieved 79.17%. Although the difference is small, the hyperparameter-tuned model slightly improved the classification performance, indicating that fine-tuning the model parameters has led to a marginal enhancement in correctly classifying instances.

**ROC AUC:** The initial random forest model had a ROC AUC of 0.9481, which was slightly higher than the hyperparameter-tuned model's ROC AUC of 0.92. Despite this small difference, both values indicate strong classification performance, with the tuned model still demonstrating excellent discriminatory power. The ROC AUC values suggest that both models can effectively distinguish between classes, but the initial model performed slightly better in this regard.

**Brier Class:** The initial random forest model had a Brier Class score of 0.1987, which indicates moderate calibration of predicted probabilities. The hyperparameter-tuned model showed the same score of 0.1987.

In addition, as it is visualized in the confusion matrix, the Random Forest model shows better classification performance for group 4 compared to other models, with fewer misclassifications. While group 3 is sometimes confused with group 4 and group 2, the overall model performs well, especially for group 4, with more correct predictions than errors.

```{r}
# Visualizar la importancia de características para Random forest
modelo_rf_final <- rf_workflow %>%
  finalize_workflow(best_rf_params) %>%
  fit(subset_train) %>%
  extract_fit_parsnip()

vip(modelo_rf_final, geom = "col", aesthetics = list(fill = "lightblue3")) +
  labs(title = "Feature Importance - Random Forest Tuned",
       x = "Features",
       y = "Importance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # Reduce size and rotate labels
    axis.title = element_text(size = 8), # Adjust axis title size
    plot.title = element_text(size = 8), # Adjust plot title size
    plot.margin = margin(8, 8, 8, 8), # Adjust margins to save space
    axis.text.y = element_text(size = 6), # Reduce size of y-axis labels (importance values)
    axis.ticks = element_line(size = 0.5), # Reduce size of axis ticks
    legend.position = "none" # Hide legend if not necessary
  ) +
  coord_flip() # Flip coordinates if the variable names are still too large
```

The pathways involved in the random forest model with tuned hyperparameters is partially aligned with those of the initial random forest model. Here, pathways are associated with metabolism, intestinal absorption, gene expression regulation (silencing), transport regulation, apoptosis, liver metabolism and toxicity of acetaminophen, regulation of cell growth and immune system responses.
