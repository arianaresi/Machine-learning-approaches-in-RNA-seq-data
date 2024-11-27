
# **Part 3: Literature Review**

## **3.1 Primary Literature Analysis**

### **Multimetric feature selection for analyzing multicategory outcomes of colorectal cancer: random forest and multinomial logistic regression models (Feng, C. H., et al., 2022)**

The study aimed to classify the four different survival outcomes of colorectal cancer (CRC) patients using RNA-seq data, focusing on identifying the most effective set of features for classification. The categories included various survival outcomes, such as "alive with no progression" and "dead with progression."

The research utilized both Random Forest (RF) and Multinomial Logistic Regression (MLR) models, optimizing feature selection through multimetric approaches.

Due to validation, authors used cross-validation, including 5-fold and 10-fold approaches to assess accuracy and other performance metrics consistently.

The methodology presented by Feng, C. H., et al. (2022). demonstrates the effectiveness of combining Random Forest and Multinomial Logistic Regression for classifying multicategory outcomes using RNA-Seq data.

### **Using supervised learning methods for gene selection in RNA-Seq Case-Control studies (Wenric & Shemirani, 2018)**

This study investigates the use of supervised learning methods to prioritize genes in RNA-Seq data for case-control studies in the context of disease-related gene expression. The authors introduced two key innovations: 
1. The use of Random Forests to rank genes based on permutation importance
2. The Extreme Pseudo-Samples (EPS) method, which leverages Variational Autoencoders (VAE) to generate extreme pseudo-samples, enhancing gene selection.

To validate their approach, the authors applied survival analysis to cancer datasets, in addition, results were validated using separate training and testing datasets, with accuracy metrics evaluated for both.

Limitations include variability in gene rankings across iterations and computational trade-offs in Random Forests, as well as inconsistent performance of the EPS method across different datasets, and a small dataset.

The methodology presented by Wenric & Shemirani (2018) is directly relevant to this project, where machine learning is used to classify and prioritize genes based on RNA-Seq data from liver cancer in rats.

### **Machine learning model for predicting Major Depressive Disorder using RNA-Seq data: optimization of classification approach (Verma & Shakya, 2021)**

The primary objective was the classification of Major Depressive Disorder (MDD), including distinguishing between suicidal and non-suicidal MDD patients. 

The key innovation was the use of Random Forest (RF) and K-Nearest Neighbors (KNN) for gene classification and feature selection, coupled with PCA to reduce data dimensionality.

To validate their approach, the authors used the model's accuracy scores. As it is well-known, Random Forest has a significant limitation, the overfitting indicated by the high accuracy on the training data compared to the lower accuracy on the test data. 

The methodology presented by Verma & Shakya (2021) illustrates the effectiveness of combining PCA and Random Forest to handle high-dimensional RNA-Seq data, which is directly applicable to the present project.

## **3.2 Methods comparison**

In the literature, various machine learning techniques are utilized to analyze RNA-Seq data and classify biological outcomes. The studies discussed in this section employed different methods, each with its unique advantages and drawbacks.

**1. Random Forest:**

In all three studies examined, random forest played a key role in both feature selection and classification. This model excels at managing high-dimensional datasets, such as RNA-Seq, because it can accommodate a large number of features without overfitting.

Although random forest is generally resistant to overfitting, it can still experience this issue if not properly tuned, particularly with smaller datasets or when the number of trees is not optimized (as noted in the Verma & Shakya, 2021 study). Furthermore, random forest's lack of interpretability can be a disadvantage in biological contexts where understanding the relationships between features is essential.

**2. Multinomial Logistic Regression:**

In Feng, C. H. et al. (2022), authors integrates Random Forest with Multinomial Logistic Regression to classify various survival outcomes in colorectal cancer patients.
MLR is a straightforward and interpretable approach for multiclass classification, making it valuable when the aim is to uncover the relationships between features and multiple outcome categories. It offers clear insights into how each predictor influences the likelihood of class membership.

On the other hand, MLR may face challenges with high-dimensional data (like RNA-Seq data) and might necessitate dimensionality reduction to be effective.

In our study of liver cancer outcomes in rats using RNA-Seq data, we first utilized Multinomial Logistic Regression (MLR) to tackle the issue of multiclass classification. This approach was selected for its clarity and ability to differentiate between various results, in line with the objective of categorizing different cancer stages or treatment results. After that, we utilized Random Forest (RF) because of its ability to handle high-dimensional biological data effectively and its capability to detect non-linear correlations among features.
