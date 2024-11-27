## **Part 4. Implementation and results**

In this study, we aimed to classify four experimental groups of rats, two of which were administered a hepatotoxin to induce liver carcinoma and two were kept as healthy controls. The analysis was based on RNA-Seq data, focusing on pathway activity scores related to different feeding protocols and the administration of the hepatotoxin.

### **Decision Trees and Random Forest models**

We compared  the performance of Decision Trees (DT) and Random Forest (RF) for the classification task. The models were evaluated using several performance metrics: accuracy, ROC AUC, and Kappa.

Visual representations of the performance measurements show that Random Forest outperforms Decision Trees across all metrics, indicating it is more suitable for managing the intricacies of RNA-Seq data. This finding is consistent with prior research, which has demonstrated that Random Forest models are highly effective at managing high-dimensional biological data and identifying intricate non-linear connections among features. 

As a component of our investigation, we utilized *hyperparameter tuning* on the Random Forest model in order to enhance its performance even more. However, the metrics indicated that the tuned model performed similarly to the original. This could be due to the initial model being well-configured, or it might suggest that further hyperparameter optimization was not necessary, despite testing models with tree counts ranging from 100 to 800 and various combinations of variable and observation numbers, as the performance remained very similar. Future analysis could be useful with more data or computing power.

It is important to mention that Decision Trees were not tuned due to inferior performance compared to Random Forest's complexity and flexibility.

### **Biological interpretation**

Both models, Logistic Multinomial Regression and Random Forest, emphasize comparable biological processes, especially in metabolism, immune system responses, gene expression regulation, and adipogenesis. The Random Forest (tuned) model emphasizes gene regulation and metabolic control mechanisms, whereas the LMR model offers a wider perspective, encompassing immune system differentiation and disease response. Both models exhibit common pathways associated with cholesterol metabolism, immune activation, and gene silencing, highlighting the significance of these processes within the dataset's context.

So far, we have extensively discussed the pathways influencing the experimental groups of rats exposed to different diets and a hepatotoxic agent. However, certain genes stand out in our models and are prominently featured in these pathways. In the **Logistic Multinomial Regression model**, two key genes, Rb1 and TP53, emerge as noteworthy.

According to UniProt, **Rb1** is a tumor suppressor that plays a crucial role in regulating the G1/S transition of the cell cycle. Additionally, the Rat Genome Database indicates that this gene is expressed across various systems, including the hepatobiliary system and gastrointestinal system, aligning with the central focus of our current project.

**TP53**, according to UniProt, encodes a multifunctional transcription factor that regulates cell cycle arrest, DNA repair, or apoptosis by binding to its target DNA sequences. It functions as a tumor suppressor in various tumor types, inducing growth arrest or apoptosis depending on the physiological context and cell type. Additionally, it negatively regulates cell division by controlling the expression of genes involved in this process. Similar to **Rb1**, Tp53 is represented in multiple tissues and systems in the Rat Genome Database, including the hepatobiliary and gastrointestinal systems.

In the **tuned Random Forest model**, the **Pparg** gene is prominently featured within the pathways identified. According to UniProt, Pparg is described as a nuclear receptor that binds to peroxisome proliferators such as hypolipidemic drugs and fatty acids. Additionally, according to Rat Genome Database, it acts as a regulator of **BMAL1** transcription, a gene involved in cardiovascular circadian rhythms. Considering that dietary timing is a critical aspect of the experimental design and can influence molecular clock genes, the role of Pparg holds significant importance in the study.
Another gene that was highlighted as important in this model is Wnt1, which according to UniProt, is a transcription factor that plays a crucial role in cellular development and cell survival. Its expression is predominantly reported during embryonic development, as stated by RGD.

To conclude, the analysis highlights the interconnected biological processes identified across both Logistic Multinomial Regression and Random Forest models, with a particular focus on metabolism, immune response, and gene regulation. Key genes such as Rb1, TP53, Pparg and Wnt1 play central roles in these pathways.

Rb1 and TP53 underscore the importance of tumor suppression, cell cycle regulation, and tissue-specific expression in systems impacted by hepatotoxicity and dietary changes. Meanwhile, Pparg adds another layer of significance by linking gene regulation, metabolic control, and circadian rhythms, aligning closely with the experimental design's focus on dietary timing and its impact on molecular clock genes. The gene Wnt1, a crucial transcription factor involved in cellular development and survival, also emerges as significant, highlighting its role in tissue development and cellular functions.

These findings suggest a complex interaction of metabolic and regulatory mechanisms influenced by diet and hepatotoxic exposure, underscoring the pathways' relevance to understanding the experimental groups. Integrating these insights furthers our understanding of how specific genes and pathways contribute to the physiological responses observed in this study.

This study's clinical and research implications showcase the possibility of treating metabolic disorders, liver diseases, and cancers by focusing on genes such as Pparg, Rb1, Wnt1 and Tp53. Pparg's involvement in circadian rhythms presents possibilities for chrononutrition and chronotherapy, whereas Tp53 and Rb1, Wnt1 could offer insights for addressing hepatotoxicity or cellular stress. Future studies should prioritize investigations into gene interactions with circadian rhythms through time-series analysis, explore the effects of diet on Wnt1, Tp53 and Rb1 through mechanistic research, and explore into other genes and pathways associated with diet-induced physiological alterations. This will enhance our comprehension and could result in innovative treatment methods. 


## **References**

Ibm. (2024, August 15). Decision tree. IBM. https://www.ibm.com/topics/decision-trees

Ibm. (2024, October 25). Random Forest. IBM. https://www.ibm.com/topics/random-forest

Chidambaranathan-Reghupaty, S., Fisher, P. B., & Sarkar, D. (2020). Hepatocellular carcinoma (HCC): Epidemiology, etiology and molecular classification. Advances In Cancer Research, 1-61. https://doi.org/10.1016/bs.acr.2020.10.001

Wenric, S., & Shemirani, R. (2018). Using supervised learning methods for gene selection in RNA-Seq Case-Control studies. Frontiers in Genetics, 9. https://doi.org/10.3389/fgene.2018.00297

Verma, P., & Shakya, M. (2021). Machine learning model for predicting Major Depressive Disorder using RNA-Seq data: optimization of classification approach. Cognitive Neurodynamics, 16(2), 443–453. https://doi.org/10.1007/s11571-021-09724-8

Feng, C. H., et al. (2022). Multimetric feature selection for analyzing multicategory outcomes of colorectal cancer: Random forest and multinomial logistic regression models. Laboratory Investigation, 102(3), 236-244. https://doi.org/10.1038/s41374-021-00662-x

UniProt. (s. f.). https://www.uniprot.org/uniprotkb/P33568/entry

Rgd. (s. f.). RAT Genome Database. https://rgd.mcw.edu/
