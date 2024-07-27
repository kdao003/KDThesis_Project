#Evaluating the Model with GSE15061, GSE111085, TCGA-AML Data

#loading required libraries
library(xgboost)
library(tidyverse)
library(caret)
library(dplyr)
library(SHAPforxgboost)
library(ROSE)
library(pROC)
library(DiagrammeR) 


#Reading in the data
meta <- read.csv("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/May_to_July_ProjectWork/JuneAndJuly_Work/XG_Boost_Model/MutationData_with_QNE_ssGSEA_Scores.csv",
                 header = TRUE)

##Loading in Model 1 
#Model 1 (Had nrounds = 100)
model_1 <- readRDS("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/May_to_July_ProjectWork/JuneAndJuly_Work/XG_Boost_Model/70_Train_30_Test/70_30_WithRemovedFeatures/xgboost_model_70_30_RemoveFeatures.rds")


#Filtering data to only the 3 datasets of interest

GSE111085 <- meta %>%
  filter(startsWith(sample, "SRR")) %>%
  select(-c(dataset, sample))

GSE15061 <- meta %>%
  filter(startsWith(sample, "GSM")) %>%
  select(-c(dataset, sample))

TCGA_AML <- meta %>%
  filter(startsWith(sample, "TCGA")) %>%
  select(-c(dataset, sample))


##Data Preparation##
# Determining Categorical Columns
categorical_columns <- c("CBFB_MYH11_Present", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
                         "priorMDS")

# Convert the columns to factors
GSE111085[categorical_columns] <- lapply(GSE111085[categorical_columns], as.factor)
GSE15061[categorical_columns] <- lapply(GSE15061[categorical_columns], as.factor)
TCGA_AML[categorical_columns] <- lapply(TCGA_AML[categorical_columns], as.factor)

# Check structure of data
str(GSE111085)
str(GSE15061)
str(TCGA_AML)


#Dropping columns not of interest
#Will need to remove "prior_hematologic_disorder_diagnosis_indicator", "leukemia_french_american_british_morphology_code" since the models were not trained on this data
blacklist <- c("priorMalignancyNonMyeloid", "priorMalignancyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", 
               "isDenovo", "isTransformed", "dataset", "prior_hematologic_disorder_diagnosis_indicator", "leukemia_french_american_british_morphology_code",
               "priorMaligncyNonMyeloid", "priorMaligncyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", "X11q23_Present", "STAG2_Mutation",
               "SRSF2_Mutation", "U2AF1_Mutation", "SF3B1_Mutation", "RUNX1_RUNX1T1_Present","BCOR_Mutation", "EZH2_Mutation",
               "ZRSR2_Mutation", "GMP.like", "Promono.like")

sub_GSE111085 <- GSE111085[, !names(GSE111085) %in% blacklist, drop = FALSE]
sub_GSE15061 <- GSE15061[, !names(GSE15061) %in% blacklist, drop = FALSE]
sub_TCGA_AML <- TCGA_AML[, !names(TCGA_AML) %in% blacklist, drop = FALSE]

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(sub_GSE111085), c("sample", categorical_columns))
sub_GSE111085[numeric_columns] <- lapply(sub_GSE111085[numeric_columns], as.numeric)

numeric_columns <- setdiff(colnames(sub_GSE15061), c("sample", categorical_columns))
sub_GSE15061[numeric_columns] <- lapply(sub_GSE15061[numeric_columns], as.numeric)

numeric_columns <- setdiff(colnames(sub_TCGA_AML), c("sample", categorical_columns))
sub_TCGA_AML[numeric_columns] <- lapply(sub_TCGA_AML[numeric_columns], as.numeric)




# Convert the target variable to a factor with 'yes' as the positive class
sub_GSE111085$priorMDS <- factor(sub_GSE111085$priorMDS, levels = c('no', 'yes'))

sub_GSE15061$priorMDS <- factor(sub_GSE15061$priorMDS, levels = c('no', 'yes'))

sub_TCGA_AML$priorMDS <- factor(sub_TCGA_AML$priorMDS, levels = c('no', 'yes'))




###GSE11085###

#Separate features and target label
GSE111085_FeaturesTest <- sub_GSE111085[, !names(sub_GSE111085) %in% "priorMDS"]
GSE111085_LabelsTest <- sub_GSE111085$priorMDS


#Convert Features and Labels to numeric
GSE111085_Num_FeaturesTest <- data.frame(sapply(GSE111085_FeaturesTest, as.numeric))
GSE111085_Num_LabelsTest <- as.numeric(GSE111085_LabelsTest) - 1


#Adding filler column to Num_FeaturesTest. Columns are X11q23 and ZRSR2_Mutation. We do so since the models were trained with these columns
GSE111085_Num_FeaturesTest <- GSE111085_Num_FeaturesTest %>%
  mutate(X11q23_Present = NA) %>%
  mutate(ZRSR2_Mutation = NA) 

#Converting these columns to numeric
GSE111085_Num_FeaturesTest$X11q23_Present <- as.numeric(GSE111085_Num_FeaturesTest$X11q23_Present)
GSE111085_Num_FeaturesTest$ZRSR2_Mutation <- as.numeric(GSE111085_Num_FeaturesTest$ZRSR2_Mutation)

# Ensure that Num_FeaturesTest columns are in the same order as model_1$feature_names
GSE111085_Num_FeaturesTest <- GSE111085_Num_FeaturesTest %>%
  select(all_of(model_1$feature_names))

#Prepare data for xgboost
GSE111085_dtest <- xgb.DMatrix(data = as.matrix(GSE111085_Num_FeaturesTest), label = GSE111085_Num_LabelsTest)







###GSE15061###

# Separate features and target label
GSE15061_FeaturesTest <- sub_GSE15061[, !names(sub_GSE15061) %in% "priorMDS"]
GSE15061_LabelsTest <- sub_GSE15061$priorMDS


#Convert Features and Labels to numeric
GSE15061_Num_FeaturesTest <- data.frame(sapply(GSE15061_FeaturesTest, as.numeric))
GSE15061_Num_LabelsTest <- as.numeric(GSE15061_LabelsTest) - 1


#Adding filler column to Num_FeaturesTest. Columns are X11q23 and ZRSR2_Mutation. We do so since the models were trained with these columns
GSE15061_Num_FeaturesTest <- GSE15061_Num_FeaturesTest %>%
  mutate(X11q23_Present = NA) %>%
  mutate(ZRSR2_Mutation = NA) 

#Converting these columns to numeric
GSE15061_Num_FeaturesTest$X11q23_Present <- as.numeric(GSE15061_Num_FeaturesTest$X11q23_Present)
GSE15061_Num_FeaturesTest$ZRSR2_Mutation <- as.numeric(GSE15061_Num_FeaturesTest$ZRSR2_Mutation)

# Ensure that Num_FeaturesTest columns are in the same order as model_1$feature_names
GSE15061_Num_FeaturesTest <- GSE15061_Num_FeaturesTest %>%
  select(all_of(model_1$feature_names))

#Prepare data for xgboost
GSE15061_dtest <- xgb.DMatrix(data = as.matrix(GSE15061_Num_FeaturesTest), label = GSE15061_Num_LabelsTest)







###TCGA_AML###

# Separate features and target label
TCGA_AML_FeaturesTest <- sub_TCGA_AML[, !names(sub_TCGA_AML) %in% "priorMDS"]
TCGA_AML_LabelsTest <- sub_TCGA_AML$priorMDS


#Convert Features and Labels to numeric
TCGA_AML_Num_FeaturesTest <- data.frame(sapply(TCGA_AML_FeaturesTest, as.numeric))
TCGA_AML_Num_LabelsTest <- as.numeric(TCGA_AML_LabelsTest) - 1


#Adding filler column to Num_FeaturesTest. Columns are X11q23 and ZRSR2_Mutation. We do so since the models were trained with these columns
TCGA_AML_Num_FeaturesTest <- TCGA_AML_Num_FeaturesTest %>%
  mutate(X11q23_Present = NA) %>%
  mutate(ZRSR2_Mutation = NA) 

#Converting these columns to numeric
TCGA_AML_Num_FeaturesTest$X11q23_Present <- as.numeric(TCGA_AML_Num_FeaturesTest$X11q23_Present)
TCGA_AML_Num_FeaturesTest$ZRSR2_Mutation <- as.numeric(TCGA_AML_Num_FeaturesTest$ZRSR2_Mutation)

# Ensure that Num_FeaturesTest columns are in the same order as model_1$feature_names
TCGA_AML_Num_FeaturesTest <- TCGA_AML_Num_FeaturesTest %>%
  select(all_of(model_1$feature_names))

#Prepare data for xgboost
TCGA_AML_dtest <- xgb.DMatrix(data = as.matrix(TCGA_AML_Num_FeaturesTest), label = TCGA_AML_Num_LabelsTest)











###GSE111085 Model Prediction with model_1###

#model prediction
pred1 <- predict(model_1, GSE111085_dtest) 
pred1 <- ifelse(pred1 > 0.5,1,0)

# Convert predictions back to factor
predictions1 <- factor(pred1, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix1 <- confusionMatrix(predictions1, factor(GSE111085_Num_LabelsTest, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix1)

# Extract precision (Positive Predictive Value)
precision1 <- conf_matrix1$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall1 <- conf_matrix1$byClass["Sensitivity"]

# Calculate F1 score
f1_score1 <- 2 * (precision1 * recall1) / (precision1 + recall1)

# Print precision and F1 score
cat("Precision: ", precision1, "\n")
cat("F1 Score: ", f1_score1, "\n")

##Plots##
# Calculate SHAP values
shap_values <- shap.values(xgb_model = model_1, X_train = as.matrix(GSE111085_Num_FeaturesTest))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = model_1, X_train = as.matrix(GSE111085_Num_FeaturesTest))
shap.plot.summary(shap_long)

# Plot ROC curve
roc_curve <- roc(GSE111085_Num_LabelsTest, pred1)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)












###GSE15061 Model Prediction with model_1###

#model prediction
###Model Prediction with model_2###

#model prediction
pred2 <- predict(model_1, GSE15061_dtest) 
pred2 <- ifelse(pred2 > 0.5,1,0)

# Convert predictions back to factor
predictions2 <- factor(pred2, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix2 <- confusionMatrix(predictions2, factor(GSE15061_Num_LabelsTest, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix2)

# Extract precision (Positive Predictive Value)
precision2 <- conf_matrix2$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall2 <- conf_matrix2$byClass["Sensitivity"]

# Calculate F1 score
f1_score2 <- 2 * (precision2 * recall2) / (precision2 + recall2)

# Print precision and F1 score
cat("Precision: ", precision2, "\n")
cat("F1 Score: ", f1_score2, "\n")

##Plots##
# Calculate SHAP values
shap_values <- shap.values(xgb_model = model_1, X_train = as.matrix(GSE15061_Num_FeaturesTest))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = model_1, X_train = as.matrix(GSE15061_Num_FeaturesTest))
shap.plot.summary(shap_long)

# Plot ROC curve
roc_curve <- roc(GSE15061_Num_LabelsTest, pred2)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)













###TCGA_AML Model Prediction with model_1###

#model prediction
pred3 <- predict(model_1, TCGA_AML_dtest) 
pred3 <- ifelse(pred3 > 0.5,1,0)

# Convert predictions back to factor
predictions3 <- factor(pred3, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix3 <- confusionMatrix(predictions3, factor(TCGA_AML_Num_LabelsTest, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix3)

# Extract precision (Positive Predictive Value)
precision3 <- conf_matrix3$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall3 <- conf_matrix3$byClass["Sensitivity"]

# Calculate F1 score
f1_score3 <- 2 * (precision3 * recall3) / (precision3 + recall3)

# Print precision and F1 score
cat("Precision: ", precision3, "\n")
cat("F1 Score: ", f1_score3, "\n")

##Plots##
# Calculate SHAP values
shap_values <- shap.values(xgb_model = model_1, X_train = as.matrix(TCGA_AML_Num_FeaturesTest))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = model_1, X_train = as.matrix(TCGA_AML_Num_FeaturesTest))
shap.plot.summary(shap_long)


###Cannot Plot ROC or obtain AUC due to there only being 1 response variable (all of TCGA-AML are de novo AML))###
# Plot ROC curve
roc_curve <- roc(TCGA_AML_Num_LabelsTest, pred3)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)
